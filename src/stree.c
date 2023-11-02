/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp.h"

#define PROP_THRESHOLD 10

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))
#define SWAP_PMAT_INDEX(e,i) (((e)+(i))%((e)<<1))
#define SWAP_SCALER_INDEX(n,i) (((n)+((i)-1))%(2*(n)-2)) 
#define GET_HINDEX(t,p) (((node_is_mirror((p)) ? \
                          (p)->node_index : (p)->hybrid->node_index)) - \
                        ((t)->tip_count+(t)->inner_count))

/* species tree spr move related */
#define LINEAGE_A       16
#define LINEAGE_OTHER   32
#define NODE_SQUARE     64
#define NODE_MOVED      256

#define SNL_PUREA       16
#define SNL_MOVED       32

static long dbg_counter = 0;
static void debug_print_migmatrix(stree_t * stree);
static void debug_print_migrations(stree_t * stree);
static void debug_print_migrations_flip(stree_t * stree);
static void debug_print_bitmatrix(stree_t * stree);

/* the variables below are used in the propose_tau move. They are allocated only
   once in stree_alloc_internals() before MCMC starts, and deallocate after MCMC
   finishes at stree_fini().
   This is to avoid reallocation at every call of the propose_tau move. 

   '__gt_nodes' is an array to hold the gene tree nodes across all gene trees,
   which were affected by the tau change.

   '__aux' is to hold the old branches (prior changing them based on the new tau)
   of the nodes in '__gt_nodes' such that in case of rejection we can revert back

   '__mark_count[i]' holds how many gene tree nodes were affected for locus i'

   '__extra_count[i]' holds how many extra gene tree nodes were affected as a result
   that they are the children of nodes in __gt_nodes, but themselves did not appear
   in __gt_nodes. Those nodes will get different branch lengths, but their age will
   not change
*/
static gnode_t ** __gt_nodes = NULL;
static double * __aux = NULL;
static int * __mark_count = NULL;
static int * __extra_count = NULL;
static long * __gt_nodes_index = NULL;

static double * target_weight = NULL;
static snode_t ** target = NULL;

/* allocated and used only for species tree inference */
static gnode_t ** moved_space;
static unsigned int * moved_count;
static gnode_t ** gtarget_temp_space;
static gnode_t ** gtarget_space;

static snode_t ** snode_contrib_space;
static unsigned int * snode_contrib_count;
/* TODO: REMOVE */
static gnode_t * pruned_nodes[10000];
static gnode_t * gsources_list[10000];

static long * ts_indicator = NULL;

#if 1
double dbg_prop_a = -10;
double dbg_prop_b = -20;
long dbg_mig_idx = -1;
long dbg_mig_idx_prop = -1;
long mig_model_prop_count[256][256] = {0};
double mig_model_prop_acc[256][256] = {0};
#endif

static long debug_rjmcmc = 0;

#define SHRINK          1
#define EXPAND          2

/* hashtable for indexing species tree labels */
hashtable_t * species_hash(stree_t * tree)
{
  long i;

  /* using a hash table check whether there are duplicate nodes */
  hashtable_t * ht = hashtable_create(tree->tip_count);

  for (i = 0; i < tree->tip_count; ++i)
  {
    /* attempt to place the pair in the hash table and die with an
       error if a pair with the same label already exists */
    if (!hashtable_insert(ht,
                          (void *)(tree->nodes[i]),
                          hash_fnv(tree->nodes[i]->label),
                          hashtable_ptrcmp))
    {
       /* this should never happen because duplicate taxa were
          already checked for during tree parsing */
      fatal("Duplicate taxon or node label (%s) in species tree",
            tree->nodes[i]->label);
    }
  }

  return ht;
}

static int longint_len(long x)
{
  return x ? (int)floor(log10(abs(x))) + 1 : 1;
}

static double logPDFGamma(double x, double a, double b)
{
   /* gamma density: mean=a/b; var=a/b^2
   */
   if (x <= 0 || a <= 0 || b <= 0) {
      printf("x=%.6f a=%.6f b=%.6f", x, a, b);
      fatal("x a b outside range in logPDFGamma()");
   }
   if (a > 30000)
      fatal("large alpha in PDFGamma()");
   return a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;
}


void print_network_table(stree_t * stree, FILE * fp)
{
  long i;
  long hybrid_count = 0;
  long bidir_count = 0;

  fprintf(fp, "Species tree contains hybridization/introgression events.\n\n");
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->hybrid)
    {
      if  (node_is_hybridization(stree->nodes[i]))
      {
        hybrid_count++;
      }
      else if (node_is_bidirection(stree->nodes[i]))
      {
        if (stree->nodes[i]->prop_tau)
          bidir_count++;
      }
      else
      {
        fatal("Internal error when counting hybridization events");
      }
    }
  fprintf(fp, "Hybridization events: %ld\n", hybrid_count);
  fprintf(fp, "Bidirectional introgressions: %ld\n", bidir_count);
    
  fprintf(fp, "Label        Node  Child1  Child2  Parent\n");
  for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
  {
    char * label;
    if (stree->nodes[i]->label)
      label = xstrdup(stree->nodes[i]->label);
    else
      label = xstrdup("N/A");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    fprintf(fp, "%-12s", label);
    free(label);
    fprintf(fp, "  %3d", stree->nodes[i]->node_index);

    if (stree->nodes[i]->left)
      fprintf(fp, "  %6d", stree->nodes[i]->left->node_index);
    else
      fprintf(fp, "  %6s", "N/A");

    if (stree->nodes[i]->right)
      fprintf(fp, "  %6d", stree->nodes[i]->right->node_index);
    else
      fprintf(fp, "  %6s", "N/A");


    if (stree->nodes[i]->parent)
      fprintf(fp, "  %6d", stree->nodes[i]->parent->node_index);
    else
      fprintf(fp, "  %6s", "N/A");

    if (i >= stree->tip_count + stree->inner_count)
    {
      if (node_is_bidirection(stree->nodes[i]))
        fprintf(fp, "   Mirrored bidirection node");
      else
        fprintf(fp, "   Mirrored hybridization node");
      assert(stree->nodes[i]->hybrid);
      if (stree->nodes[i]->hybrid->label)
        fprintf(fp,
                " [Hybrid = %s (%d)]",
                stree->nodes[i]->hybrid->label,
                stree->nodes[i]->hybrid->node_index);
    }

    if (stree->nodes[i]->hybrid)
    {
      fprintf(fp,
              "   [tau = %ld, phi = %f, prop_tau = %d, has_phi = %ld]",
              stree->nodes[i]->htau,
              stree->nodes[i]->hphi,
              stree->nodes[i]->prop_tau,
              stree->nodes[i]->has_phi);

      #if 0

      /* old code before the introduction of has_phi */

      /* for bidirections print the mirror node phi, for hybridizations print
         both node phi. */
      if (i >= stree->tip_count+stree->inner_count)
      {
        if (node_is_bidirection(stree->nodes[i]))
          fprintf(fp,
                  "  phi_%s : %s -> %s",
                  stree->nodes[i]->label,
                  stree->nodes[i]->parent->label,
                  stree->nodes[i]->label);
        else
        {
          /* hybridization node */

          /* if main node htau==0 and mirror node htau==1 then that is the only
             case we use the phi from the main node */
          snode_t * tmpnode = stree->nodes[i];
          if (tmpnode->hybrid->htau == 0 && tmpnode->htau == 1)
          {
            fprintf(fp,
                    "  1-phi_%s : %s -> %s",
                    stree->nodes[i]->label,
                    stree->nodes[i]->parent->label,
                    stree->nodes[i]->label);
          }
          else
          {
            fprintf(fp,
                    "  phi_%s : %s -> %s",
                    stree->nodes[i]->label,
                    stree->nodes[i]->parent->label,
                    stree->nodes[i]->label);
          }
        }
      }
      else
      {
        if (!node_is_bidirection(stree->nodes[i]))
        {
          /* hybridization - main node */
          snode_t * tmpnode = stree->nodes[i];
          if (tmpnode->htau == 0 && tmpnode->hybrid->htau == 1)
          {
            fprintf(fp,
                    " phi_%s : %s -> %s",
                    stree->nodes[i]->label,
                    stree->nodes[i]->parent->label,
                    stree->nodes[i]->label);
          }
          else
          {
            fprintf(fp,
                    " 1-phi_%s : %s -> %s",
                    stree->nodes[i]->label,
                    stree->nodes[i]->parent->label,
                    stree->nodes[i]->label);
          }
        }
      }
      #else
      snode_t * tmpnode = stree->nodes[i];
      if (tmpnode->has_phi)
        fprintf(fp,
                "  phi_%s : %s -> %s",
                stree->nodes[i]->label,
                stree->nodes[i]->parent->label,
                stree->nodes[i]->label);

      #endif
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
}

void debug_print_network_node_attribs(stree_t * stree)
{
  long i,j;
  long msa_index;
  unsigned int total_nodes;
  char * s;

  int pad = 8;
  int minpad = 5;
  int dots = 3;
  int spc = 1;
  if (pad < minpad)
    fatal("Please set padding (pad) to at least %ld\n", minpad);
  if (dots >= pad)
    fatal("Please set number of dots (dots) to at most %ld\n", pad-1);
  
  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  /* print header */
  printf("The following table prints [seqin/events/leaves] for each locus at each species tree node\n");
  printf("%-*s", pad, "Locus");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("%-*s", spc, " ");
    
    s = NULL;
    xasprintf(&s, "%ld (%s)", i, stree->nodes[i]->label);
    if (strlen(s) > (size_t)pad)
    {
      for (j = 0; j < dots; ++j)
        s[strlen(s)-j-1] = '.';
      s[pad] = 0;
    }
      
    printf("%*s", pad, s);
  }
  printf("\n");

  for (msa_index = 0; msa_index < stree->locus_count; ++msa_index)
  {
    printf("%-*ld", pad, msa_index);

    for (i = 0; i < total_nodes; ++i)
    {
      printf("%-*s", spc, " ");
      s = NULL;
      xasprintf(&s,
                "%d/%d/%d",
                stree->nodes[i]->seqin_count[msa_index],
                stree->nodes[i]->event_count[msa_index],
                stree->nodes[i]->gene_leaves[msa_index]);
      if (strlen(s) > (size_t)pad)
      {
        for (j = 0; j < dots; ++j)
          s[strlen(s)-j-1] = '.';
        s[pad] = 0;
      }
      printf("%*s", pad, s);
      free(s);
    }
    printf("\n");
  }
}

void stree_show_pptable(stree_t * stree, int show_taus_and_thetas)
{
  long i, j;
  long nodes_count = stree->tip_count + stree->inner_count + stree->hybrid_count;
  int index_digits = longint_len(nodes_count);
  size_t maxlen = strlen("Species");

  for (i = 0; i < nodes_count; ++i)
    maxlen = MAX(maxlen, strlen(stree->nodes[i]->label));

  printf("\nMap of populations and ancestors (1 in map indicates ancestor):\n");

  /* print space before horizontal species enumeration */
  for (i = 0; i < index_digits; ++i)
    printf(" ");
  printf(" ");
  printf("Species");
  for (i = 0; i < (long)(maxlen - strlen("Species")); ++i)
    printf(" ");

  printf(" ");
  for (i = 0; i < nodes_count; ++i)
    printf("  %ld", i+1);
  printf("\n");

  for (i = 0; i < nodes_count; ++i)
  {
    printf("%*ld %-*s ", index_digits, i+1, (int)maxlen, stree->nodes[i]->label);

    for (j = 0; j < nodes_count; ++j)
      printf("  %*d", longint_len(j), stree->pptable[i][j]);

    if (show_taus_and_thetas)
    {
      printf("  tau = %f  theta = %f",
             stree->nodes[i]->tau, stree->nodes[i]->theta);
    }

    printf("\n");
  }

  printf("\n");
}

static void snode_clone(snode_t * snode, snode_t * clone, stree_t * clone_stree)
{
  unsigned int i;
  unsigned int msa_count = clone_stree->locus_count;

  if (clone->label)
    free(clone->label);

  clone->length = snode->length;
  clone->theta = snode->theta;
  clone->tau = snode->tau;
  clone->old_tau = snode->old_tau;
  clone->old_theta = snode->old_theta;
  clone->leaves = snode->leaves;
  clone->support = snode->support;
  clone->weight = snode->weight;
  clone->node_index = snode->node_index;
  clone->has_theta = snode->has_theta;
  clone->constraint = snode->constraint;
  clone->constraint_lineno = snode->constraint_lineno;
  clone->theta_step_index = snode->theta_step_index;

  if (!clone->mark)
    clone->mark = (int *)xmalloc((size_t)opt_threads*sizeof(int));
  memcpy(clone->mark, snode->mark, (size_t)opt_threads*sizeof(int));
  //clone->mark = snode->mark;

  /* points to relatives */
  if (snode->parent)
    clone->parent = clone_stree->nodes[snode->parent->node_index];
  else
    clone->parent = NULL;

  if (snode->left)
    clone->left = clone_stree->nodes[snode->left->node_index];
  else
    clone->left = NULL;

  if (snode->right)
    clone->right = clone_stree->nodes[snode->right->node_index];
  else
    snode->right = NULL;

  /* label */
  if (snode->label)
    clone->label = xstrdup(snode->label);
  else
    clone->label = NULL;

  /* gene leaves */
  if (!clone->gene_leaves)
    clone->gene_leaves = (unsigned int *)xmalloc(msa_count * sizeof(unsigned int));
  memcpy(clone->gene_leaves, snode->gene_leaves, msa_count * sizeof(unsigned int));

  /* data  - unused */
  clone->data = NULL;

  /* event doubly-linked lists */
  if (!clone->event)
  {
    clone->event = (dlist_t **)xcalloc(msa_count, sizeof(dlist_t *));
    for (i = 0; i < msa_count; ++i)
      clone->event[i] = dlist_create();
  }
  else
  {
    for (i = 0; i < msa_count; ++i)
      dlist_clear(clone->event[i], NULL);
  }

  /* event counts per locus */
  if (!clone->event_count)
    clone->event_count = (int *)xmalloc(msa_count * sizeof(int));
  memcpy(clone->event_count, snode->event_count, msa_count * sizeof(int));

  /* branch rate (per locus) */
  if (snode->brate)
  {
    if (!clone->brate)
      clone->brate = (double *)xmalloc(msa_count * sizeof(double));
    memcpy(clone->brate, snode->brate, msa_count * sizeof(double));
  }
    
  /* linked theta models */
  if (snode->linked_theta)
    clone->linked_theta = clone_stree->nodes[snode->linked_theta->node_index];
  else
    clone->linked_theta = NULL;

  /* seqin (incoming sequences to a population) counts per locus */
  if (!clone->seqin_count)
    clone->seqin_count = (int *)xmalloc(msa_count * sizeof(int));
  memcpy(clone->seqin_count, snode->seqin_count, msa_count * sizeof(int));

  /* per locus number of gene leaves at each clade */
  if (!clone->gene_leaves)
    clone->gene_leaves = (unsigned int*)xmalloc(msa_count*sizeof(unsigned int));
  memcpy(clone->gene_leaves,snode->gene_leaves,msa_count*sizeof(unsigned int));

  /* gene tree probability contributions for current population */
  if (!clone->logpr_contrib)
    clone->logpr_contrib = (double *)xmalloc(msa_count * sizeof(double));
  memcpy(clone->logpr_contrib, snode->logpr_contrib, msa_count*sizeof(double));

  /* old gene tree probability contributions for current population */
  if (!clone->old_logpr_contrib)
    clone->old_logpr_contrib = (double *)xmalloc(msa_count * sizeof(double));
  memcpy(clone->old_logpr_contrib,
         snode->old_logpr_contrib,
         msa_count * sizeof(double));

  if (!opt_est_theta)
  {
    clone->t2h_sum = snode->t2h_sum;
    clone->event_count_sum = snode->event_count_sum;
    clone->notheta_logpr_contrib = snode->notheta_logpr_contrib;
    clone->notheta_old_logpr_contrib = snode->notheta_old_logpr_contrib;

    if (!clone->t2h)
      clone->t2h = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
    memcpy(clone->t2h, snode->t2h, opt_locus_count * sizeof(double));

    if (!clone->old_t2h)
      clone->old_t2h = (double*)xmalloc((size_t)opt_locus_count*sizeof(double));
    memcpy(clone->old_t2h, snode->old_t2h, opt_locus_count * sizeof(double));
  }

  if (opt_migration)
  {
    for (i = 0; i < msa_count; ++i)
      clone->migevent_count[i] = snode->migevent_count[i];

    clone->mb_count = snode->mb_count;
    clone->mb_mrsum_isarray = snode->mb_mrsum_isarray;
    for (i = 0; i < clone_stree->inner_count; ++i)
    {
      clone->migbuffer[i].time = snode->migbuffer[i].time;
      clone->migbuffer[i].type = snode->migbuffer[i].type;
      clone->migbuffer[i].active_count= snode->migbuffer[i].active_count;
      memcpy(clone->migbuffer[i].mrsum,
             snode->migbuffer[i].mrsum,
             snode->migbuffer[i].active_count * sizeof(double));
    }
  }
}

static void gnode_clone(gnode_t * gnode,
                        gnode_t * clone,
                        gtree_t * clone_gtree,
                        stree_t * clone_stree)
{
  if (clone->label)
    free(clone->label);

  clone->length = gnode->length;
  clone->time = gnode->time;
  clone->old_time = gnode->old_time;
  clone->leaves = gnode->leaves;
  clone->node_index = gnode->node_index;
  clone->clv_valid = gnode->clv_valid;
  clone->clv_index = gnode->clv_index;
  clone->scaler_index = gnode->scaler_index;
  clone->pmatrix_index = gnode->pmatrix_index;
  clone->mark = gnode->mark;


  /* points to relatives */
  if (gnode->parent)
    clone->parent = clone_gtree->nodes[gnode->parent->node_index];
  else
    clone->parent = NULL;

  if (gnode->left)
    clone->left = clone_gtree->nodes[gnode->left->node_index];
  else
    clone->left = NULL;

  if (gnode->right)
    clone->right = clone_gtree->nodes[gnode->right->node_index];
  else
    clone->right = NULL;

  /* label */
  if (gnode->label)
    clone->label = xstrdup(gnode->label);
  else
    clone->label = NULL;

  /* data - unused */
  clone->data = NULL;

  /* populations */
  if (gnode->pop)
    clone->pop = clone_stree->nodes[gnode->pop->node_index];
  else
    clone->pop = NULL;

  if (gnode->old_pop)
    clone->old_pop = clone_stree->nodes[gnode->old_pop->node_index];
  else
    clone->old_pop = NULL;

  miginfo_clean(clone->mi,clone_gtree->msa_index); 
  if (gnode->mi)
    miginfo_clone(gnode->mi, &(clone->mi), clone_stree, clone_gtree->msa_index);
}

static void stree_clone(stree_t * stree, stree_t * clone)
{
  unsigned int i;
  unsigned nodes_count = stree->tip_count + stree->inner_count;

  assert(!opt_msci);

  /* clone node contents */
  for (i = 0; i < nodes_count; ++i)
    snode_clone(stree->nodes[i], clone->nodes[i], clone);

  /* clone pptable */
  for (i = 0; i < nodes_count; ++i)
    memcpy(clone->pptable[i], stree->pptable[i], nodes_count * sizeof(int));

  clone->root = clone->nodes[stree->root->node_index];

  /* notheta attributes */

  clone->notheta_logpr = stree->notheta_logpr;
  clone->notheta_old_logpr = stree->notheta_old_logpr;
  clone->notheta_hfactor = stree->notheta_hfactor;
  clone->notheta_sfactor = stree->notheta_sfactor;

  /* relaxed clock attributes */
  clone->locusrate_mubar = stree->locusrate_mubar;
  clone->locusrate_nubar = stree->locusrate_nubar;
  clone->nui_sum = stree->nui_sum;
  
  if (opt_migration)
  {
    for (i = 0; i < nodes_count; ++i)
      memcpy(clone->migcount_sum[i],
             stree->migcount_sum[i],
             nodes_count*sizeof(long));
  }
}

stree_t * stree_clone_init(stree_t * stree)
{
  unsigned int i;
  unsigned int nodes_count = stree->tip_count + stree->inner_count;
  long j;
  stree_t * clone;
  long ** migcount_sum = NULL;
  
  assert(!opt_msci);
  clone = (stree_t *)xcalloc(1, sizeof(stree_t));
  memcpy(clone, stree, sizeof(stree_t));

  /* create cloned species tree nodes */
  clone->nodes = (snode_t **)xmalloc(nodes_count * sizeof(snode_t *));
  clone->td = (snode_t **)xmalloc(nodes_count * sizeof(snode_t *));
  for (i = 0; i < nodes_count; ++i)
  {
    clone->nodes[i] = (snode_t *)xcalloc(1, sizeof(snode_t));
    snode_t * x = clone->nodes[i];

    if (opt_migration)
    {
      x->migevent_count = (long *)xcalloc((size_t)opt_locus_count,sizeof(long));
      x->migbuffer = (migbuffer_t *)xcalloc((size_t)(stree->inner_count),
                                            sizeof(migbuffer_t));
      size_t allocsize = stree->nodes[i]->mb_mrsum_isarray ? opt_locus_count : 1;
      for (j = 0; j < stree->inner_count; ++j)
        x->migbuffer[j].mrsum = (double *)xmalloc(allocsize*sizeof(double));

      x->mb_count = 0;

      x->mig_source = (dlist_t **)xcalloc((size_t)opt_locus_count,
                                          sizeof(dlist_t *));
      x->mig_target = (dlist_t **)xcalloc((size_t)opt_locus_count,
                                          sizeof(dlist_t *));
      for (j = 0; j < opt_locus_count; ++j)
      {
        x->mig_source[j] = dlist_create();
        x->mig_target[j] = dlist_create();
      }
    }
  }
  for (i = 0; i < nodes_count; ++i)
    snode_clone(stree->nodes[i], clone->nodes[i], clone);

  clone->pptable = (int **)xmalloc(nodes_count * sizeof(int *));
  for (i = 0; i < nodes_count; ++i)
  {
    clone->pptable[i] = (int *)xmalloc(nodes_count * sizeof(int));
    memcpy(clone->pptable[i], stree->pptable[i], nodes_count * sizeof(int));
  }
  clone->root = clone->nodes[stree->root->node_index];

  clone->mi_tbuffer = NULL;
  clone->migcount_sum = NULL;
  if (opt_migration)
  {
    void * mem = xmalloc((size_t)(nodes_count*nodes_count)*sizeof(long) +
                         (size_t)nodes_count*sizeof(long *));
    migcount_sum = (long **)mem;
    migcount_sum[0] = (long *)(migcount_sum+nodes_count);
    memcpy(migcount_sum[0],stree->migcount_sum[0],nodes_count*sizeof(long));
    for (i = 1; i < nodes_count; ++i)
    {
      migcount_sum[i] = (long *)(migcount_sum[i-1] + nodes_count);
      memcpy(migcount_sum[i],stree->migcount_sum[i],nodes_count*sizeof(long));
    }
    clone->migcount_sum = migcount_sum;

    clone->mi_tbuffer = (miginfo_t **)xcalloc((size_t)opt_threads,
                                              sizeof(miginfo_t *));
  }

  return clone;
}

static void gtree_clone(gtree_t * gtree,
                        gtree_t * clone_gtree,
                        stree_t * clone_stree)
{
  unsigned int i;
  unsigned nodes_count = gtree->tip_count + gtree->inner_count;

  for (i = 0; i < nodes_count; ++i)
    gnode_clone(gtree->nodes[i], clone_gtree->nodes[i], clone_gtree, clone_stree);

  clone_gtree->root = clone_gtree->nodes[gtree->root->node_index];

  clone_gtree->logl = gtree->logl;
  clone_gtree->logpr = gtree->logpr;
  clone_gtree->old_logl = gtree->old_logl;
  clone_gtree->old_logpr = gtree->old_logpr;

  /* locus rate */
  clone_gtree->rate_mui = gtree->rate_mui;
  clone_gtree->rate_nui = gtree->rate_nui;
  clone_gtree->lnprior_rates = gtree->lnprior_rates;

  if (opt_migration)
  {
    unsigned int m,n;
    for (m = 0; m < clone_stree->tip_count+clone_stree->inner_count; ++m)
      for (n = 0; n < clone_stree->tip_count+clone_stree->inner_count; ++n)
        clone_gtree->migcount[m][n] = gtree->migcount[m][n];
  }
}

gtree_t * gtree_clone_init(gtree_t * gtree, stree_t * clone_stree)
{
  unsigned int i;
  unsigned nodes_count = gtree->tip_count + gtree->inner_count;
  gtree_t * clone;
  long ** migcount = NULL;
  snode_t ** migpops = NULL;

  clone = (gtree_t *)xcalloc(1, sizeof(gtree_t));
  memcpy(clone, gtree, sizeof(gtree_t));

  /* create cloned gene tree nodes */
  clone->nodes = (gnode_t **)xmalloc(nodes_count * sizeof(gnode_t *));
  for (i = 0; i < nodes_count; ++i)
    clone->nodes[i] = (gnode_t *)xcalloc(1, sizeof(gnode_t));
  for (i = 0; i < nodes_count; ++i)
    gnode_clone(gtree->nodes[i], clone->nodes[i], clone, clone_stree);

  size_t alloc_size = MAX(4,nodes_count);
  /* TODO: Change to xmalloc for the first-touch numa policy */
  clone->travbuffer = (gnode_t **)xcalloc(alloc_size,sizeof(gnode_t*));
  clone->root = clone->nodes[gtree->root->node_index];

  if (opt_migration && !opt_simulate)
  {
    unsigned int total_nodes = clone_stree->tip_count+clone_stree->inner_count;

    void * mem = xmalloc((size_t)(total_nodes*total_nodes)*sizeof(long) +
                         (size_t)(total_nodes)*sizeof(long *));
    migcount = (long **)mem;
    migcount[0] = (long *)(migcount+total_nodes);
    memset(migcount[0],0,total_nodes*sizeof(long));
    for (i = 1; i < total_nodes; ++i)
    {
      migcount[i] = (long *)(migcount[i-1] + total_nodes);
      memset(migcount[i],0,total_nodes*sizeof(long));
    }
    migpops = (snode_t **)xcalloc((size_t)(total_nodes),sizeof(snode_t *));

    if (opt_exp_imrb)
      clone->rb_linked = (snode_t **)xmalloc((size_t)(nodes_count+1) *
                                             sizeof(snode_t *));
  }
  clone->migcount = migcount;
  clone->migpops = migpops;

  return clone;
}

static void events_clone(stree_t * stree,
                         stree_t * clone_stree,
                         gtree_t ** clone_gtree_list)
{
  unsigned int i, j;
  unsigned int msa_count = clone_stree->locus_count;
  unsigned stree_nodes_count = stree->tip_count + stree->inner_count;
  dlist_item_t * item;
  dlist_item_t * cloned;

  for (i = 0; i < stree_nodes_count; ++i)
  {
    for (j = 0; j < msa_count; ++j)
    {
      gtree_t * clone_gtree = clone_gtree_list[j];

      for (item = stree->nodes[i]->event[j]->head; item; item = item->next)
      {
        gnode_t * original_node = (gnode_t *)(item->data);
        unsigned int node_index = original_node->node_index;
        gnode_t * cloned_node = (gnode_t *)(clone_gtree->nodes[node_index]);

        cloned = dlist_append(clone_stree->nodes[i]->event[j], cloned_node);

        cloned_node->event = cloned;
      }
    }
  }
}

static void stree_label_recursive(snode_t * node)
{
  /* if node is a tip return */
  if (!node->left)
    return;

  if (node->left)
    stree_label_recursive(node->left);
  else
    fatal("Specified species tree is not binary");

  if (node->right)
    stree_label_recursive(node->right);
  else
    fatal("Specified species tree is not binary");

  if (opt_migration && node->label)
  {
    assert(!opt_est_stree);
    return;
  }

  if (node->label)
    free(node->label);

  node->label = (char *)xmalloc(strlen(node->left->label) +
                                strlen(node->right->label) + 1);

  /* concatenate species labels */
  node->label[0] = 0;
  strcat(node->label, node->left->label);
  strcat(node->label, node->right->label);
}

static void stree_label_network_recursive(snode_t * node)
{
  if (!node->left && !node->right) return;

  if (node->left)
    stree_label_network_recursive(node->left);

  if (node->right)
    stree_label_network_recursive(node->right);

  if (node->hybrid)
  {
    #if 0
    assert(node->hybrid->label || node->label)

    if (!node->hybrid->label)
      node->hybrid->label = xstrdup(node->label);
    else
      node->label = xstrdup(node->hybrid->label);
    #endif

    if (node_is_bidirection(node))
    {
      /* do nothing */
    }
    else
    {
      /* i guess do nothing here as well - everything will be sorted out by
         itself */
    }
  }

  /* keep node has a label, keep it */
  if (!node->label)
  {
    size_t len = 0;
    if (node->left)
      len += strlen(node->left->label);
    if (node->right)
      len += strlen(node->right->label);

    node->label = (char *)xmalloc((len+1)*sizeof(char));
    node->label[0] = 0;
    if (node->left)
      strcat(node->label, node->left->label);
    if (node->right)
      strcat(node->label, node->right->label);
  }
}

void stree_label(stree_t * stree)
{
  if  (opt_msci)
    stree_label_network_recursive(stree->root);
  else
    stree_label_recursive(stree->root);
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

/* return an array of per-locus number of sequences for each species (tip in
   species tree) */
static int ** populations_seqcount(stree_t * stree,
                                   msa_t ** msalist,
                                   list_t * maplist,
                                   int msa_count)
{
  unsigned int i, j, k;
  snode_t * node;
  hashtable_t * sht = NULL;
  hashtable_t * mht = NULL;

  assert(msa_count >= 0);

  /* create one hash table of species and one for sequence->species mappings */
  if (stree->tip_count > 1)
  {
    sht = species_hash(stree);
    mht = maplist_hash(maplist, sht);
  }

  /* create perloci sequence counters for tip and inner nodes */
  int ** seqcount = (int **)xmalloc(stree->tip_count * sizeof(int *));
  for (i = 0; i < stree->tip_count; ++i)
    seqcount[i] = (int *)xcalloc(msa_count, sizeof(int));

  /* one species case */
  if (stree->tip_count == 1)
  {
    for (j = 0; j < (unsigned int)(msa_count); ++j)
      seqcount[0][j] = msalist[j]->count;

    return seqcount;
  }

  /* go through the alignments and match each sequence with the corresponding
     species */
  for (i = 0; i < (unsigned int)(msa_count); ++i)
  {
    msa_t * msa = msalist[i];

    assert(msa->count >= 0);

    /* go through the sequences of the current locus and match each sequence
       with the corresponding population using its species tag */


    for (j = 0; j < (unsigned int)(msa->count); ++j)
    {
      /* first get the species tag */
      char * label = msa->label[j];
      label = strchr(label, '^');
      if (!label)
        fatal("Cannot find species tag on sequence %s of locus %d",
              msa->label[j], i);

      /* skip the '^' mark */
      label++;
      if (!(*label))
        fatal("Sequence %s of locus %d contains no label",
              msa->label[j], i);

      pair_t * pair;
      pair = hashtable_find(mht,
                            (void *)label,
                            hash_fnv(label),
                            cb_cmp_pairlabel);
      if (!pair)
        fatal("Cannot find a mapping to species for tag %s inside file %s",
              label, opt_mapfile);

      node = (snode_t *)(pair->data);

      #if 0
      if (opt_debug)
        printf("Matched %s -> %s\n", label, node->label);
      #endif

      /* increment sequence count for loci i in corresponding population */
      k = node->node_index;
      seqcount[k][i]++;
    }
  }

  hashtable_destroy(sht, NULL);
  hashtable_destroy(mht, cb_dealloc_pairlabel);

  return seqcount;
}

/* TODO: Terribly slow (quadratic to number of nodes) method for initializing
   tau on species tree nodes when we have hybridizations. The method basically
   iterates the species tree nodes array indefinitely, and at each iteration
   it sets the tau for the nodes whose parent(s) have been already processed.
   The iteration stops when no nodes need to be processed */
static void network_init_tau_iterative(stree_t * stree,
                                       double prop,
                                       long thread_index)
{
  long run = 1;
  long i;

  assert(opt_msci);
  assert(stree->root->tau && stree->root->tau != 1);

  while (run)
  {
    run = 0;
    for (i = 0; i < stree->inner_count; ++i)
    {
      snode_t * x = stree->nodes[stree->tip_count+i];
      if (!x->parent) continue;

      if (x->tau == 1) assert(x->parent->tau > 0);
      if (x->hybrid && x->tau)
        assert(x->hybrid->parent->tau > 0);

      if (x->hybrid && x->tau)
      {
        if (node_is_hybridization(x))
        {

          /* hybridization nodes */

          if (x->htau && x->parent->tau == 1)
          {
              run  = 1;
              continue;
          }
          if (x->hybrid->htau && x->hybrid->parent->tau == 1)
          {
            run = 1;
            continue;
          }

          if (x->htau == 0)
          {
            assert(x->parent->parent);
            assert(x->parent->parent->tau > 0);
            if (x->parent->parent->tau == 1)
            {
              run = 1;
              continue;
            }
          }

          if (x->hybrid->htau == 0)
          {
            assert(x->hybrid->parent->parent);
            assert(x->hybrid->parent->parent->tau > 0);
            if (x->hybrid->parent->parent->tau == 1)
            {
              run = 1;
              continue;
            }
          }
          double age1 = (x->htau) ?
                          x->parent->tau : x->parent->parent->tau;
          double age2 = (x->hybrid->htau) ?
                          x->hybrid->parent->tau : x->hybrid->parent->parent->tau;

          if (x->tau != 1)
            continue;

          x->tau = MIN(age1,age2) *
                   (prop + (1 - prop - 0.02)*legacy_rndu(thread_index));
          x->hybrid->tau = x->tau;
          if (x->htau == 0)
            x->parent->tau = x->tau;
          if (x->hybrid->htau == 0)
            x->hybrid->parent->tau = x->tau;
        }
        else
        {

          /* bidirectional introgression nodes */

          assert(node_is_bidirection(x));

          /* TODO: Account for parallel bidirectional introgressions among two
             lineages with no speciations in between. */

          assert(!node_is_mirror(x));
          if (x->parent->tau == 1 || x->right->hybrid->parent->tau == 1)
          {
            run = 1;
            continue;
          }

          assert(x->hybrid->parent);
          assert(x->hybrid->parent->hybrid);

          if (x->tau != 1) continue;

          assert(x->hybrid->tau == 1 &&
                 x->right->tau == 1 &&
                 x->right->hybrid->tau == 1);

          double age = MIN(x->parent->tau,x->right->hybrid->parent->tau)*
                       (prop + (1 - prop - 0.02)*legacy_rndu(thread_index));

          x->tau                = age;
          x->hybrid->tau        = age;
          x->right->tau         = age;
          x->right->hybrid->tau = age;

        }

      }
      else
      {
        if (x->parent->tau)
        {
          if (x->parent->tau == 1)
          {
            run = 1;
            continue;
          }
          else
          {
            if (x->tau > 0 && x->tau == 1)
            {
              if (x->prop_tau)
              {
                x->tau = x->parent->tau *
                         (prop + (1 - prop - 0.02)*legacy_rndu(thread_index));
              }
              else
              {
                run = 1;
                continue;
              }
            }
          }
        }
      }
    }
  }
}

static void stree_init_tau_recursive(snode_t * node,
                                     double prop,
                                     long thread_index)
{
  assert(!opt_msci);

  /* end recursion if node is a tip */
  if (!node->left && !node->right)
  {
    if (!node->parent->tau)
      node->theta = -1;

    return;
  }

  /* get species record associate with species tree node */
  double tau_parent = node->parent->tau;

  if (!node->parent->tau)
    node->theta = -1;

  if (node->parent->tau && node->tau > 0)
    node->tau = tau_parent * (prop + (1 - prop - 0.02)*legacy_rndu(thread_index));
  else
    node->tau = 0;

  stree_init_tau_recursive(node->left, prop, thread_index);
  stree_init_tau_recursive(node->right, prop, thread_index);
}

static void stree_init_tau(stree_t * stree, long thread_index)
{
  unsigned int i;
  unsigned int total_nodes;
  
  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  for (i = stree->tip_count; i < total_nodes; ++i)
    stree->nodes[i]->tau = 1;

   if (opt_method == METHOD_10)    /* method A10 */
   {
     if (opt_msci)
       fatal("Method A10 with hybridizations not implemented yet");

     double r = legacy_rndu(thread_index);
     int index = (int)(r * delimitation_getparam_count());
     delimitation_set(stree, index);

     char * s = delimitation_getparam_string();
     printf("Starting delimitation: %s\n", s);
   }
   else if (opt_method == METHOD_11)
   {
     if (opt_msci)
       fatal("Method A11 with hybridizations not implemented yet");

     double r = (long)(stree->tip_count*legacy_rndu(thread_index));
     if (r < stree->tip_count - 1)
       for (i = stree->tip_count; i < stree->tip_count * 2 - 1; ++i)
         stree->nodes[i]->tau = !stree->pptable[i][stree->tip_count + (long)r];
   }
   /* Initialize speciation times for each extinct species */

   double prop = (stree->root->leaves > PROP_THRESHOLD) ? 0.9 : 0.5;

   /* set the speciation time for root */
   if (stree->root->tau)
   {
     if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
       stree->root->tau = opt_tau_beta / (opt_tau_alpha - 1) *
                          (0.9 + 0.2*legacy_rndu(thread_index));
     else
     {
       #if 0
       stree->root->tau = opt_tau_alpha / opt_tau_beta;
       #else
       stree->root->tau = opt_tau_alpha / opt_tau_beta *
                          (0.9 + 0.2*legacy_rndu(thread_index));
       #endif
     }
   }

   /* recursively set the speciation time for the remaining inner nodes. For
      networks it is not necessary to check if root has both left and right */
   if (opt_msci)
   {
     network_init_tau_iterative(stree, prop, thread_index);
   }
   else
   {
     stree_init_tau_recursive(stree->root->left, prop, thread_index);
     stree_init_tau_recursive(stree->root->right, prop, thread_index);
   }

   /* check to see if everything is OK */
   if (opt_msci)
   {
     for (i = 0; i < stree->hybrid_count; ++i)
     {
       snode_t * h = stree->nodes[stree->tip_count+stree->inner_count+i];
       assert(h && h->hybrid);
       assert(node_is_mirror(h));
       if (node_is_hybridization(h))
         assert((h->parent->tau >= h->tau) && (h->hybrid->tau == h->tau));
       else
       {
         assert(node_is_bidirection(h));
         assert(h->hybrid->parent->tau > h->tau);
         assert(h->hybrid->tau == h->tau &&
                h->hybrid->right->tau == h->tau &&
                h->hybrid->right->hybrid->tau == h->tau);
         assert(h->hybrid->right->hybrid->parent->tau > h->tau);
       }
     }
   }
}

static int propose_phi(stree_t * stree,
                       gtree_t ** gtree,
                       snode_t * snode,
                       long thread_index)
{
  int accepted;
  int sequp_count;
  long i;
  double phinew;
  double phiold;
  double new_logpr;
  double old_logpr;
  double lnacceptance;
  double lnphiratio, lnphiratio1;
  double aratio, bratio;

  snode_t * pnode;   /* node that has the phi parameter */

  /* Note: snode is the main node */
  assert(!node_is_mirror(snode));

  /* Note: after our email exchange with Adam in June 2022 we introduced a flag
     (has_phi) on snode indicating the presence of the phi parameter. To
     minimize the changes to the phi proposal, we still pass the 'main' node to
     the function, even though it may not have the parameter. We keep all
     notations intact and we propose/accept/reject the new phi value on the
     node that has the parameter.
  */

  #if 0

  /* old code before the introduction of has_phi flag */
  phiold = snode->hphi;

  #else

  /* new correct code */
  pnode = snode->has_phi ? snode : snode->hybrid;
  assert(pnode->has_phi);
  phiold = pnode->hphi;


  #endif

  phinew = phiold + opt_finetune_phi*legacy_rnd_symmetrical(thread_index);
  phinew = reflect(phinew,0,1,thread_index);

  /* new correct code after the introduciton of has_phi flag
     determine the right ratios for each node */
  if (snode == pnode)
  {
    aratio = lnphiratio = log(phinew / phiold);
    bratio = lnphiratio1 = log((1-phinew) / (1-phiold));
  }
  else
  {
    aratio = lnphiratio1 = log(phinew / phiold);
    bratio = lnphiratio  = log((1-phinew) / (1-phiold));
  }

  if (opt_est_theta)
  {
    old_logpr = 0;
    new_logpr = 0;
    for (i = 0; i < stree->locus_count; ++i)
    {
      /* For bidirectional introgression we need to subtract the lineages
         coming from right. See issue #97 */
      sequp_count = snode->seqin_count[i];
      if (node_is_bidirection(snode))
        sequp_count -= snode->right->seqin_count[i];

      old_logpr += gtree[i]->logpr;
      new_logpr += gtree[i]->logpr +
                   sequp_count*lnphiratio +
                   snode->hybrid->seqin_count[i]*lnphiratio1;
    }
  }
  else
  {
    old_logpr = stree->notheta_logpr;
    new_logpr = stree->notheta_logpr;
    for (i = 0; i < stree->locus_count; ++i)
    {
      /* For bidirectional introgression we need to subtract the lineages
         coming from right. See issue #97 */
      sequp_count = snode->seqin_count[i];
      if (node_is_bidirection(snode))
        sequp_count -= snode->right->seqin_count[i];

      new_logpr += sequp_count*lnphiratio +
                   snode->hybrid->seqin_count[i]*lnphiratio1;
    }
  }

  lnacceptance = (opt_phi_alpha-1) * aratio +
                 (opt_phi_beta-1) * bratio +
                 new_logpr - old_logpr;

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */

    accepted = 1;

    #if 0
    /* old code before the introduction of has_phi flag */
    snode->hphi = phinew;
    snode->hybrid->hphi = 1-phinew;
    #else
    /* new correct code */
    pnode->hphi = phinew;
    pnode->hybrid->hphi = 1-phinew;
    #endif

    /* update logpr */
    if (opt_est_theta)
    {
      for (i = 0; i < stree->locus_count; ++i)
      {
        /* subtract from gene tree log-density the old MSCi contributions */
        gtree[i]->logpr -= snode->logpr_contrib[i] + 
                           snode->hybrid->logpr_contrib[i];

        /* For bidirectional introgression we need to subtract the lineages
           coming from right. See issue #97 */
        sequp_count = snode->seqin_count[i];
        if (node_is_bidirection(snode))
          sequp_count -= snode->right->seqin_count[i];

        /* update log-density contributions for the two populations */
        snode->logpr_contrib[i] += sequp_count*lnphiratio;
        snode->hybrid->logpr_contrib[i] += snode->hybrid->seqin_count[i] *
                                           lnphiratio1;

        /* add to the gene tree log-density the new phi contributions */
        gtree[i]->logpr += snode->logpr_contrib[i] +
                           snode->hybrid->logpr_contrib[i];
      }
    }
    else
    {
      /* subtract from total density (sum of densities of all gene trees) the
         old phi contributions from the two populations */
      stree->notheta_logpr -= snode->notheta_logpr_contrib;
      stree->notheta_logpr -= snode->hybrid->notheta_logpr_contrib;

      snode->notheta_logpr_contrib -= snode->hphi_sum;
      snode->hybrid->notheta_logpr_contrib -= snode->hybrid->hphi_sum;
      /* compute new contributions */
      for (i = 0; i < stree->locus_count; ++i)
      {
        /* For bidirectional introgression we need to subtract the lineages
           coming from right. See issue #97 */
        sequp_count = snode->seqin_count[i];
        if (node_is_bidirection(snode))
          sequp_count -= snode->right->seqin_count[i];

        snode->hphi_sum -= snode->notheta_phi_contrib[i];
        snode->hybrid->hphi_sum -= snode->hybrid->notheta_phi_contrib[i];

        snode->notheta_phi_contrib[i] += sequp_count*lnphiratio;
        snode->hybrid->notheta_phi_contrib[i] += snode->hybrid->seqin_count[i] *
                                                lnphiratio1;

        snode->hphi_sum += snode->notheta_phi_contrib[i];
        snode->hybrid->hphi_sum += snode->hybrid->notheta_phi_contrib[i];
      }
      snode->notheta_logpr_contrib += snode->hphi_sum;
      snode->hybrid->notheta_logpr_contrib += snode->hybrid->hphi_sum;

      /* add new contributions to total log-density */
      stree->notheta_logpr += snode->notheta_logpr_contrib;
      stree->notheta_logpr += snode->hybrid->notheta_logpr_contrib;
    }
  }
  else
  {
    /* rejected */
    accepted = 0;
  }
  return accepted;
}

double stree_propose_phi(stree_t * stree, gtree_t ** gtree)
{
  long i;
  long accepted = 0;

  long thread_index = 0;

  /* offset to start of hybridization nodes */
  long offset = stree->tip_count + stree->inner_count;

  /* go through hybridization nodes */
  for (i = 0; i < stree->hybrid_count; ++i)
    accepted += propose_phi(stree,
                            gtree,
                            stree->nodes[offset+i]->hybrid,
                            thread_index);

  return ((double)accepted / stree->hybrid_count);
}

static void stree_init_phi(stree_t * stree)
{
  long i;

  long offset = stree->tip_count + stree->inner_count;
  long thread_index = 0;

  if (opt_phi_alpha <= 0)
    fatal("Alpha value for 'phiprior' must be larger than 0");
  if (opt_phi_beta <= 0)
    fatal("Beta value for 'phiprior' must be larger than 0");

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * mnode = stree->nodes[offset+i];

    assert(node_is_bidirection(mnode) || node_is_mirror(mnode));

    /* set phi parameter to the mean */
    if (fabs(mnode->hphi + mnode->hybrid->hphi - 1) > 1e-10)
    {
      if (node_is_bidirection(mnode->hybrid))
      {
        /* bidirection (model D) */

        /* we set phi for the vertical branch to U|(0.7,0.9). The mnode is
           the mirror node, hence it is the horizontal branch */
        double a = 0.7; double b = 0.9;
        double r = (b-a)*legacy_rndu(thread_index) + a;
        
        mnode->hybrid->hphi = r;
        mnode->hphi = 1-r;
      }
      else
      {
        /* hybridization */

        /* for models A and C draw the value of phi from U(0,1). For model B
           set phi for the vertical branch to U(0.7,0.9) */
        if ((!mnode->htau && !mnode->hybrid->htau) ||
            (mnode->htau && mnode->hybrid->htau))
        {
          /* model A or C */
          mnode->hphi = legacy_rndu(thread_index);
          mnode->hybrid->hphi = 1 - mnode->hphi;
        }
        else
        {
          /* model B */
          double a = 0.7; double b = 0.9;
          double r = (b-a)*legacy_rndu(thread_index) + a;
          if (mnode->htau)
          {
            mnode->hphi = r;
            mnode->hybrid->hphi = 1 - r;
          }
          else
          {
            mnode->hybrid->hphi = r;
            mnode->hphi = 1 - r;
          }
        }
      }
    }
  }
}

static void msci_link_thetas_bfs(stree_t * stree)
{
  long i,j = 0;
  long head,tail;
  long total_nodes;
  snode_t * sibling;
  snode_t ** queue;
  snode_t ** hybrid;

  total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;
  queue = (snode_t **)xcalloc((size_t)(total_nodes+1),sizeof(snode_t *));
  hybrid = (snode_t **)xcalloc((size_t)(stree->hybrid_count),sizeof(snode_t *));
  
  /* order nodes using BFS */  
  head = tail = 0;
  queue[tail++] = stree->root;
  while (queue[head])
  {
    snode_t * x = queue[head++];

    if (x->left)
      queue[tail++] = x->left;

    if (x->right)
      queue[tail++] = x->right;

    if (x->hybrid && x->node_index < stree->tip_count+stree->inner_count)
      hybrid[j++] = x;

    printf(" %s", x->label);
  }
  printf("\n");
  free(queue);

  printf("Hybrid nodes: %ld\n", j);
  assert(j == stree->hybrid_count);

  /* now get list of hybrid nodes in bfs order */
  for (i = 0; i < j; ++i)
    printf(" %s", hybrid[i]->label);
  printf("\n");

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * snode = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;
    snode_t * mnode = snode->hybrid;

    if (!node_is_bidirection(snode))
    {
      /* hybridization */

      /* model A */
      if (snode->htau && mnode->htau) continue;

      /* model C */
      if (!snode->htau) 
      {
        /* sibling is linked to parent */
        sibling = (snode->parent->left == snode) ?
                    snode->parent->right : snode->parent->left;
        sibling->linked_theta = snode->parent->linked_theta ?
                                  snode->parent->linked_theta : snode->parent;
      }
      else
      {
        /* child is linked to hybrid */
        assert(snode->left && !snode->right);

        snode_t * child = snode->left;
        child->linked_theta = snode->linked_theta ? snode->linked_theta : snode;
      }

      if (!mnode->htau)
      {
        sibling = (mnode->parent->left == mnode) ?
                    mnode->parent->right : mnode->parent->left;
        sibling->linked_theta = mnode->parent->linked_theta ?
                            mnode->parent->linked_theta : mnode->parent;
      }
      else
      {
        assert(!mnode->left && !mnode->right && snode->left && !snode->right);
        snode_t * child = snode->left;
        child->linked_theta = mnode->linked_theta ? mnode->linked_theta : mnode;
      }
    }
    else
    {
      /* bidirection */
      snode->left->linked_theta = snode->linked_theta ?
                                    snode->linked_theta : snode;
    }
  }
}

static void init_theta_stepsize(stree_t * stree)
{
  long i,j;
  unsigned int total_nodes;
  long theta_params = 0;

  /* theta mode sets the number of steplengths:

     1: a single step length for all thetas
     2: two step lengths, one for thetas in tips and one for thetas in inner
        nodes
     3: one step length for each theta
  */

  total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  if (opt_finetune_theta_mode == 1 ||
      opt_linkedtheta == BPP_LINKEDTHETA_ALL ||
      stree->tip_count == 1)
  {
    /* single step length for all thetas */
    for (i = 0; i < total_nodes; ++i)
      stree->nodes[i]->theta_step_index = 0;
    opt_finetune_theta_count = 1;
  }
  else
  {
    double eps = opt_finetune_theta[0];

    free(opt_finetune_theta);

    if (opt_finetune_theta_mode == 2)
    {
      /* two step lengths, one for tips and one for inner nodes */
      opt_finetune_theta = (double *)xmalloc(2*sizeof(double));
      opt_finetune_theta[0] = eps;
      opt_finetune_theta[1] = eps;

      for (i = 0; i < total_nodes; ++i)
        if (stree->nodes[i]->linked_theta == NULL)
        {
          if (i < stree->tip_count)
            stree->nodes[i]->theta_step_index = 0;
          else
            stree->nodes[i]->theta_step_index = 1;
        }
      opt_finetune_theta_count = 2;
    }
    else
    {
      /* one step length for each theta */
      for (i = 0; i < total_nodes; ++i)
        if (stree->nodes[i]->linked_theta == NULL)
          theta_params++;
      
      opt_finetune_theta = (double *)xmalloc((size_t)(theta_params)*sizeof(double));
      opt_finetune_theta_count = theta_params;
      for (i = 0; i < theta_params; ++i)
        opt_finetune_theta[i] = eps;

      for (i = 0, j = 0; i < total_nodes; ++i)
        if (stree->nodes[i]->linked_theta == NULL)
          stree->nodes[i]->theta_step_index = j++;

    }
  }
}

static void init_theta_linkage(stree_t * stree)
{
  long i;
  unsigned int total_nodes;

  total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->linked_theta = NULL;

  if (opt_linkedtheta == BPP_LINKEDTHETA_NONE) return;

  if (opt_linkedtheta == BPP_LINKEDTHETA_ALL)
  {
    for (i = 0; i < total_nodes; ++i)
    {
      if (stree->nodes[i] != stree->root)
        stree->nodes[i]->linked_theta = stree->root;
    }
  }
  else if (opt_linkedtheta == BPP_LINKEDTHETA_INNER)
  {
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    {
      if (stree->nodes[i] != stree->root)
        stree->nodes[i]->linked_theta = stree->root;
    }
  }
  else if (opt_linkedtheta == BPP_LINKEDTHETA_MSCI)
  {
    msci_link_thetas_bfs(stree);
  }
}

static void stree_init_theta(stree_t * stree,
                             msa_t ** msalist,
                             list_t * maplist,
                             int msa_count,
                             FILE * fp_out,
                             long thread_index)
{
  long abort = 0;
  long warn = 0;
  unsigned int i, j;


  /* initialize population sizes for extinct populations and populations
     with more than one lineage at some locus */

     /* get an array of per-locus number of sequences for each species (tip in
        species tree) */
  int ** seqcount = populations_seqcount(stree, msalist, maplist, msa_count);

  /* Check number of sequences per locus with stated numbers from 'species&tree'
     tag in control file. Assume the following:

     X = specified max number of species in control file
     Y = maximum number in sequence file

     The behavious is the following:
     if (X == 1 && Y > 1) abort with error message
     if (X > 1 && Y <= 1) print warning on screen and in output file.Dont abort
     if (X > 1 && Y > 1) no error or warning as the numbers do not affect run

     See also issue #62 on GitHub
  */

  /* print table header on screen and in output file */
  fprintf(stdout, "\nPer-locus sequences in data and 'species&tree' tag:\n");
  fprintf(stdout, 
          "C.File | Data |                Status                | Population\n");
  fprintf(stdout,
          "-------+------+--------------------------------------+-----------\n");
  fprintf(fp_out, "\nPer-locus sequences in data and 'species&tree' tag:\n");
  fprintf(fp_out,
          "C.File | Data |                Status                | Population\n");
  fprintf(fp_out,
          "-------+------+--------------------------------------+-----------\n");

  for (i = 0; i < stree->tip_count; ++i)
  {
    int maxseqcount = 0;
    for (j = 0; j < opt_locus_count; ++j)
      if (seqcount[i][j] > maxseqcount)
        maxseqcount = seqcount[i][j];


    /* distinguish between cases and print corresponding status message */
    if (opt_sp_seqcount[i] == 1 && maxseqcount > 1)
    {
      abort = 1;
      fprintf(stdout,
              "%6ld | %4d | %-36s | %-10s\n",
              opt_sp_seqcount[i], maxseqcount,
              "[ERROR] Increase number in C.File",
              stree->nodes[i]->label);
      fprintf(fp_out,
              "%6ld | %4d | %-36s | %-10s\n",
              opt_sp_seqcount[i], maxseqcount,
              "[ERROR] Increase number in C.File",
              stree->nodes[i]->label);
    }
    else if (opt_sp_seqcount[i] > 1 && maxseqcount <= 1)
    {
      warn = 1;
      fprintf(stdout,
              "%6ld | %4d | %-36s | %-10s\n",
              opt_sp_seqcount[i], maxseqcount,
              "[WARNING] Parameter not identifiable",
              stree->nodes[i]->label);
      fprintf(fp_out,
              "%6ld | %4d | %-36s | %-10s\n",
              opt_sp_seqcount[i], maxseqcount,
              "[WARNING] Parameter not identifiable",
              stree->nodes[i]->label);
    }
    else
    {
      fprintf(stdout,
              "%6ld | %4d | %-36s | %-10s\n",
              opt_sp_seqcount[i], maxseqcount,
              "[OK]",
              stree->nodes[i]->label);
      fprintf(fp_out,
              "%6ld | %4d | %-36s | %-10s\n",
              opt_sp_seqcount[i], maxseqcount,
              "[OK]",
              stree->nodes[i]->label);
    }
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
  if (warn)
  {
    fprintf(stdout,
            "[Warning] Some parameters are not identifiable because there are "
            "0 or 1 sequences from some species, but the control file lists "
            "different numbers.\nThose entries are indicated in the table "
            "above. Posterior for theta parameters for those species will be "
            "given by the prior.\n\n");
    fprintf(fp_out,
            "[Warning] Some parameters are not identifiable because there are "
            "0 or 1 sequences from some species, but the control file lists "
            "different numbers.\nThose entries are indicated in the table "
            "above. Posterior for theta parameters for those species will be "
            "given by the prior.\n\n");
  }
  if (abort)
  {
    fprintf(fp_out,
            "[Error] Some populations consist of more than one sequence but "
            "control file states only one.\nPlease amend control file according"
            " to the table above.");
    fatal("[Error] Some populations consist of more than one sequence but "
          "control file states only one.\nPlease amend control file according "
          "to the table above.");
  }

  /* initialize 'has_theta' attribute */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (opt_est_theta)
      stree->nodes[i]->has_theta = 1;
    else
      stree->nodes[i]->has_theta = 0;

#if 0
   /* go through tip nodes and setup thetas only for those that have
      two sequences in some loci */
   for (i = 0; i < stree->tip_count; ++i)
   {
      snode_t * node = stree->nodes[i];

      for (j = 0; j < (unsigned int)msa_count; ++j)
         if (seqcount[i][j] >= 2)
            break;

      /* if no loci exists with two or more sequences of such species then move
         to the next tip node */
      if (j == (unsigned int)msa_count)
      {
         node->theta = -1;
         node->has_theta = 0;
         continue;
      }

      /* otherwise set theta around the mean of the inverse gamma prior */
      node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
         (0.9 + 0.2 * legacy_rndu());
   }
#else
   /* From Ziheng's email from 3.3.2018 (see also issue #62) this is changed
      to set estimation of thetas according to the 'species&tree' tag in the
      control file. */
   init_theta_linkage(stree);
   init_theta_stepsize(stree);
   for (i = 0; i < stree->tip_count; ++i)
   {
     snode_t * node = stree->nodes[i];

     for (j = 0; j < (unsigned int)msa_count; ++j)
       if (seqcount[i][j] >= 2)
         break;

     /*** Ziheng $$$ ***/
     /* If(opt_est_geneflow), all tip species have theta.
     *  If(opt_migration), under a fixed IM model, tips involved in migration have theta.
     *  Donor population must have theta so that coalescent rate is defined.
     *  Recipient population must have theta so that w = 4M/theta is defined.
     *  It is possible that either donor or recipient population has only 0 or 1 sequence.
     */
     long mig_i_donor_or_recipient = 0;
     if (!opt_est_geneflow && opt_migration)
     {
       for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
       {
         /*** Ziheng $$$ ***/
         mig_i_donor_or_recipient += opt_mig_bitmatrix[j][node->node_index];
         mig_i_donor_or_recipient += opt_mig_bitmatrix[node->node_index][j];
       }
     }
     /* if no loci exists with two or more sequences of such species then move
        to the next tip node */
     if (!opt_est_geneflow && (!opt_migration || !mig_i_donor_or_recipient))
     {
       if (opt_sp_seqcount[i] < 2)
       {
         node->theta = -1;
         node->has_theta = 0;
         continue;
       }
     }

     /* otherwise set theta around the mean of the prior */
     if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
       node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                     (0.9 + 0.2 * legacy_rndu(thread_index));
     else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
     {
       #if 0
       node->theta = opt_theta_alpha / opt_theta_beta;
       #else
       node->theta = opt_theta_alpha / opt_theta_beta *
                     (0.6+0.8*legacy_rndu(thread_index));
       #endif
     }
     else
     {
       assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
       double bmean = opt_theta_p / (opt_theta_p+opt_theta_q);
       double r = (0.8 + 0.2*legacy_rndu(thread_index))*bmean;
       node->theta = opt_theta_min + r*(opt_theta_max-opt_theta_min);
     }
   }
#endif

  /* go through inner nodes and setup thetas */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * node = stree->nodes[i];

    if (opt_msci && node->hybrid)
    {
      node->theta = node->hybrid->theta = -1;
      node->has_theta = node->hybrid->has_theta = 0;

      if (!node_is_bidirection(node))
      {
        /* node is a hybridization: we assign a theta to the nodes that
           compose it that have a 'tau-parent' (htau) annotation */

        if (node->htau)
        {
          if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
            node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                          (0.9 + 0.2 * legacy_rndu(thread_index));
          else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
            node->theta = opt_theta_alpha / opt_theta_beta *
                          (0.6+0.8*legacy_rndu(thread_index));
          else
          {
            assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
            double bmean = opt_theta_p / (opt_theta_p+opt_theta_q);
            double r = (0.8 + 0.2*legacy_rndu(thread_index))*bmean;
            node->theta = opt_theta_min + r*(opt_theta_max-opt_theta_min);
          }
          node->has_theta = 1;
        }
        else
        {
          node->theta = -1;
          node->has_theta = 0;
        }
        
        if (node->hybrid->htau)
        {
          if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
            node->hybrid->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                                  (0.9 + 0.2 * legacy_rndu(thread_index));
          else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
            node->hybrid->theta = opt_theta_alpha / opt_theta_beta *
                                  (0.6+0.8*legacy_rndu(thread_index));
          else
          {
            assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
            double bmean = opt_theta_p / (opt_theta_p+opt_theta_q);
            double r = (0.8 + 0.2*legacy_rndu(thread_index))*bmean;
            node->hybrid->theta = opt_theta_min+r*(opt_theta_max-opt_theta_min);
          }
          node->hybrid->has_theta = 1;
        }
        else
        {
          node->hybrid->theta = -1;
          node->hybrid->has_theta = 0;
        }
      }
      else
      {
        /* bidirectional introgression */

        if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
          node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                        (0.9 + 0.2 * legacy_rndu(thread_index));
        else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
          node->theta = opt_theta_alpha / opt_theta_beta *
                        (0.6+0.8*legacy_rndu(thread_index));
        else
        {
          assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
          double bmean = opt_theta_p / (opt_theta_p+opt_theta_q);
          double r = (0.8 + 0.2*legacy_rndu(thread_index))*bmean;
          node->theta = opt_theta_min + r*(opt_theta_max-opt_theta_min);
        }
        node->has_theta = 1;

        /* the mirrored nodes do not have a theta */
        node->hybrid->theta = -1;
        node->hybrid->has_theta = 0;
      }
    }
    else
    {
      /* 'orginary' inner nodes all have a theta. We discussed with Ziheng the
         case of inner nodes that have an incoming number of lineages equal to 1
         whether they should have a theta or not, and decided to keep it for
         code simplicity */
      if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
        node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                      (0.9 + 0.2 * legacy_rndu(thread_index));
      else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
      {
        #if 0
        node->theta = opt_theta_alpha / opt_theta_beta;
        #else
        node->theta = opt_theta_alpha / opt_theta_beta *
                      (0.6+0.8*legacy_rndu(thread_index));
        #endif
      }
      else
      {
        assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
        double bmean = opt_theta_p / (opt_theta_p+opt_theta_q);
        double r = (0.8 + 0.2*legacy_rndu(thread_index))*bmean;
        node->theta = opt_theta_min + r*(opt_theta_max-opt_theta_min);
      }
    }
  }

  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
  {
    if (stree->nodes[i]->linked_theta && stree->nodes[i]->has_theta)
      stree->nodes[i]->theta = stree->nodes[i]->linked_theta->theta;
  }

  /* deallocate seqcount */
  for (i = 0; i < stree->tip_count; ++i)
    free(seqcount[i]);
  free(seqcount);
}

/* bottom up filling of pptable */
static void stree_reset_pptable_tree(stree_t * stree)
{
  unsigned int i;
  snode_t * ancnode;
  snode_t * curnode;

  /* zero-out pptable */
  for (i=0; i < stree->tip_count + stree->inner_count; ++i)
    memset(stree->pptable[i],
           0,
           (stree->tip_count + stree->inner_count) * sizeof(int));

  for (i = 0; i < stree->tip_count; ++i)
  {
    for (curnode = stree->nodes[i]; curnode; curnode = curnode->parent)
      for (ancnode = curnode; ancnode; ancnode = ancnode->parent)
        stree->pptable[curnode->node_index][ancnode->node_index] = 1;
  }
}

int node_is_bidirection(const snode_t * node)
{
  assert(node && node->hybrid);

  const snode_t * hnode;    /* hybrid node */
  const snode_t * mnode;    /* mirror node */

  if (!node->left && !node->right)
    mnode = node;
  else
    mnode = node->hybrid;

  hnode = mnode->hybrid;

  if (hnode->left && hnode->right) return 1;

  assert(!hnode->right);

  return 0;
}


/* TODO: A simpler way would be just to check whether node->node_index >=
   stree->tip_count + stree->inner_count */
int node_is_mirror(const snode_t * node)
{
  assert(node);
  assert(node->hybrid);
  assert(node && node->hybrid);

  if (node->left || node->right)
    return 0;

  return 1;
}

int node_is_hybridization(const snode_t * node)
{
  assert(node && node->hybrid);

  if (!node) return 0;
  if (!node->hybrid) return 0;

  return !node_is_bidirection(node);

  #if 0
  if ((!node->left && !node->right) ||
      (!node->hybrid->left && !node->hybrid->right))
    return 1;
  #if 0
  /* check if it's bidirectional introgression */
  if (node->hybrid->parent->hybrid && node->hybrid->parent->parent == node->hybrid)
    return 0;
  #endif

  return 0;
  #endif
    
}

static void stree_reset_pptable_network_recursive(stree_t * stree,
                                                  snode_t * node,
                                                  snode_t * ancestor)
{
  stree->pptable[node->node_index][ancestor->node_index] = 1;

  if (node->left)
    stree_reset_pptable_network_recursive(stree,
                                          node->left,
                                          ancestor);
  if (node->right)
    stree_reset_pptable_network_recursive(stree,
                                          node->right,
                                          ancestor);

  if (node->hybrid)
  {
    /* check whether hybridization or bidirectional introgression */
    if (node_is_bidirection(node))
    {
      if (node_is_mirror(node))
      {
        /* move to the opposite child */
        stree_reset_pptable_network_recursive(stree,node->hybrid->left,ancestor);
      }
    }
    else
    {
      assert(!!node->left + !!node->right + !!node->hybrid->left + !!node->hybrid->right == 1);
      if (node->left)
        stree_reset_pptable_network_recursive(stree, node->left, ancestor);
      else if (node->hybrid->left)
        stree_reset_pptable_network_recursive(stree,node->hybrid->left,ancestor);
      else
        fatal("Internal error [hybrid] (stree_reset_pptable_network_recursive)");
    }
  }
}

/* top down filling of pptable */
static void stree_reset_pptable_network(stree_t * stree)
{
  unsigned int i;
  unsigned int nodes_count;
  
  nodes_count = stree->tip_count + stree->inner_count + stree->hybrid_count;

  /* zero-out pptable */
  for (i=0; i < nodes_count; ++i)
    memset(stree->pptable[i],
           0,
           nodes_count*sizeof(int));

  for (i=0; i < nodes_count; ++i)
    stree_reset_pptable_network_recursive(stree,stree->nodes[i],stree->nodes[i]);
}

void stree_reset_pptable(stree_t * stree)
{
  if (opt_msci)
    stree_reset_pptable_network(stree);
  else
    stree_reset_pptable_tree(stree);
}

void stree_alloc_internals(stree_t * stree, long * locus_seqcount, unsigned int gtree_inner_sum, long msa_count)
{
  long i;
   /* allocate traversal buffer to be the size of all nodes for all loci */
 //  unsigned int sum_count = 0;
   unsigned long sum_nodes = 2 * gtree_inner_sum + msa_count;
   //  for (i = 0; i < (unsigned int)msa_count; ++i)
   //    sum_count += (unsigned int)(msa[i]->count);

   //  sum_count = 2*sum_count - msa_count;
   //  __gt_nodes = (gnode_t **)xmalloc(sum_count * sizeof(gnode_t *));
   //  __aux = (double *)xmalloc((sum_count - msa_count)*sizeof(double *));
   __gt_nodes = (gnode_t **)xmalloc(sum_nodes * sizeof(gnode_t *));
   //__aux = (double *)xmalloc((sum_nodes - msa_count)*sizeof(double));
   __aux = (double *)xmalloc((sum_nodes) * sizeof(double));

   /* The following two arrays are used purely for the tau proposal.
      Entry i of marked_count indicates how many nodes from locus i are marked.
      Similarly, extra_count
      how many extra nodes where added in __gt_nodes whose branch lengths (and
      therefore) p-matrices need updating because their parent node's age was
      changed */
   __mark_count = (int *)xmalloc(msa_count * sizeof(int));
   __extra_count = (int *)xmalloc(msa_count * sizeof(int));

   __gt_nodes_index = (long *)xmalloc((size_t)msa_count * sizeof(long));
   __gt_nodes_index[0] = 0;
   for (i = 1; i < msa_count; ++i)
   {
     __gt_nodes_index[i] = __gt_nodes_index[i-1] + (2*locus_seqcount[i-1]-1);
   }
   assert(sum_nodes == (unsigned long)(__gt_nodes_index[msa_count-1] + 2*locus_seqcount[msa_count-1]-1));

   /* species tree inference */
   if (opt_est_stree)
   {
      unsigned int stree_nodes = stree->inner_count + stree->tip_count +
                                 stree->hybrid_count;

      target_weight = (double *)xmalloc((size_t)(stree_nodes) * sizeof(double));
      target = (snode_t **)xmalloc((size_t)(stree_nodes) * sizeof(snode_t *));

      /* TODO: memory is allocated for all loci to aid parallelization */
      moved_count = (unsigned int *)xcalloc(msa_count, sizeof(unsigned int));
      //    unsigned int sum_inner = 0;
      //    for (i = 0; i < (unsigned int)msa_count; ++i)
      //      sum_inner += msa[i]->count-1;
      //    moved_space = (gnode_t **)xmalloc(sum_inner * sizeof(gnode_t *));
      //    gtarget_space = (gnode_t **)xmalloc(sum_inner * sizeof(gnode_t *));
      moved_space = (gnode_t **)xmalloc(gtree_inner_sum * sizeof(gnode_t *));
      gtarget_space = (gnode_t **)xmalloc(gtree_inner_sum * sizeof(gnode_t *));

      //    unsigned int sum_nodes = 2*sum_inner + msa_count;
      gtarget_temp_space = (gnode_t **)xmalloc(sum_nodes * sizeof(gnode_t *));

      snode_contrib_space = (snode_t **)xmalloc((size_t)(msa_count*stree_nodes) *
         sizeof(snode_t *));
      snode_contrib_count = (unsigned int *)xmalloc((size_t)msa_count *
         sizeof(unsigned int));
   }
}

void stree_init_pptable(stree_t * stree)
{
  size_t pptable_size;
  unsigned int i;

#if 0
   snode_t * curnode;
   snode_t * ancnode;
#endif
  assert(opt_msci == !!stree->hybrid_count);

  pptable_size = stree->tip_count + stree->inner_count + stree->hybrid_count;

  /* initialize pptable, where pptable[i][j] indicates whether population j
     (that is, node with index j) is ancestral to population i */
  stree->pptable = (int **)xcalloc(pptable_size, sizeof(int *));
  for (i = 0; i < pptable_size; ++i)
    stree->pptable[i] = (int *)xcalloc(pptable_size, sizeof(int));

  stree_reset_pptable(stree);
#if 0
   for (i = 0; i < stree->tip_count; ++i)
   {
      for (curnode = stree->nodes[i]; curnode; curnode = curnode->parent)
         for (ancnode = curnode; ancnode; ancnode = ancnode->parent)
            stree->pptable[curnode->node_index][ancnode->node_index] = 1;
   }
#endif
}

static void network_init_hx(stree_t * stree)
{
  long i;

  #ifdef DEBUG_THREADS
  if (opt_threads == 1)
  {
    opt_threads = DEBUG_THREADS_COUNT;
    for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
    {
      snode_t * node = stree->nodes[i];
      node->hx = (long *)xcalloc((size_t)opt_threads,sizeof(long));
    }
    opt_threads = 1;
  }
  else
  {
    for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
    {
      snode_t * node = stree->nodes[i];
      node->hx = (long *)xcalloc((size_t)opt_threads,sizeof(long));
    }
  }
  #else
  for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
  {
    snode_t * node = stree->nodes[i];
    node->hx = (long *)xcalloc((size_t)opt_threads,sizeof(long));
  }
  #endif
}

static void stree_reset_leaves_network(stree_t * stree)
{
  /* NOTE: Assumes that pptable is computed correctly */
  unsigned int i,j;
  for (i = 0; i < stree->tip_count; ++i)
    stree->nodes[i]->leaves = 1;

  /* quadratic algorithm to count number of leaves for every inner node */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * snode = stree->nodes[i];

    for (j = 0; j < stree->tip_count; ++j)
      if (stree->pptable[j][snode->node_index])
        snode->leaves++;
  }

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * snode = stree->nodes[stree->tip_count+stree->inner_count+i];

    snode->leaves = snode->hybrid->leaves;
  }
}

static void stree_reset_leaves_tree_recursive(snode_t * snode)
{
  if (!snode->left)
  {
    snode->leaves = 1;
    return;
  }

  if (snode->left)
    stree_reset_leaves_tree_recursive(snode->left);
  if (snode->right)
    stree_reset_leaves_tree_recursive(snode->right);

  snode->leaves = 0;
  if (snode->left)
    snode->leaves = snode->left->leaves;
  if (snode->right)
    snode->leaves += snode->right->leaves;
}

static void stree_reset_leaves_tree(stree_t * stree)
{
  stree_reset_leaves_tree_recursive(stree->root);
}

void stree_reset_leaves(stree_t * stree)
{
  if (opt_msci)
    stree_reset_leaves_network(stree);
  else
    stree_reset_leaves_tree(stree);
}

static void msci_validate(stree_t * stree)
{
  unsigned int i;

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    /* hybridizations */

    snode_t * mnode = stree->nodes[stree->tip_count+stree->inner_count+i];

    assert(node_is_mirror(mnode));
    if (!node_is_hybridization(mnode)) continue;

    snode_t * hnode = mnode->hybrid;

    unsigned int nhindex = hnode->node_index;
    unsigned int phindex = hnode->parent->node_index;
    unsigned int pmindex = mnode->parent->node_index;

    if (stree->pptable[phindex][nhindex] || stree->pptable[pmindex][nhindex])
      fatal("[ERROR] "
            "Parental nodes of hybridization %s cannot be descendants "
            "(cannot hybridize from future past)", hnode->label);

    if (!hnode->htau && stree->pptable[pmindex][phindex])
      fatal("[ERROR] "
            "Parental nodes of hybridization %s have an ancestor-descendent "
            "relation, but the ancestor has no tau paremeter (tau-parent=no)",
            hnode->label);

    if (!mnode->htau && stree->pptable[phindex][pmindex])
      fatal("[ERROR] "
            "Parental nodes of hybridization %s have an ancestor-descendent "
            "relation, but the ancestor has no tau paremeter (tau-parent=no)",
            hnode->label);
  }


  for (i = 0; i < stree->hybrid_count; ++i)
  {
    /* bidirections */

    snode_t * mnode = stree->nodes[stree->tip_count+stree->inner_count+i];

    assert(node_is_mirror(mnode));
    if (node_is_hybridization(mnode)) continue;

    /* check for prop_tau, i.e. process each bidirection only once */
    snode_t * h1node = mnode->hybrid;
    if (!h1node->prop_tau) continue;

    snode_t * h2node = mnode->parent;
    assert(!h2node->prop_tau);

    unsigned int h1index = h1node->node_index;
    unsigned int h2index = h2node->node_index;

    if (stree->pptable[h1index][h2index] || stree->pptable[h2index][h1index])
      fatal("[ERROR] "
            "The two end-point nodes of bidirection %s <-> %s have an "
            "ancestor-descendent relation",
            h1node->label, h2node->label);
  }
}

void stree_init(stree_t * stree,
                msa_t ** msa,
                list_t * maplist,
                int msa_count,
                FILE * fp_out)
{
  unsigned int i, j;
  unsigned int nodes_count;
  long ** migcount_sum = NULL;

  long thread_index = 0;

  /* safety check */
  assert(msa_count > 0);
  assert(opt_msci == !!stree->hybrid_count);

  stree_init_pptable(stree);
  
  stree_reset_leaves(stree);

  if (opt_msci)
  {
    msci_validate(stree);
  }

  /* label each inner node of the species tree with the concatenated labels of
     its two children */
  stree_label(stree);

  /* Initialize population sizes */
  stree_init_theta(stree, msa, maplist, msa_count, fp_out, thread_index);

  if (stree->tip_count > 1 || opt_msci)
  {
    /* Initialize speciation times and create extinct species groups */
    stree_init_tau(stree, thread_index);
  }
  else
  {
    stree->nodes[0]->tau = 0;
  }
  nodes_count = stree->tip_count+stree->inner_count+stree->hybrid_count;

  stree->td = (snode_t **)xmalloc((size_t)nodes_count * sizeof(snode_t *));

  if (opt_msci)
    stree_init_phi(stree);

  stree->mi_tbuffer = NULL;
  stree->migcount_sum = NULL;
  if (opt_migration)
  {
    /* reset the number of migration events associated with each population */
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      stree->nodes[i]->migevent_count = (long *)xcalloc((size_t)opt_locus_count,
                                                        sizeof(long));
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode_t * x = stree->nodes[i];

      x->migbuffer  = (migbuffer_t *)xcalloc((size_t)(stree->inner_count),
                                             sizeof(migbuffer_t));
      size_t allocsize = x->mb_mrsum_isarray ? opt_locus_count : 1;
      for (j = 0; j < stree->inner_count; ++j)
        x->migbuffer[j].mrsum = (double *)xmalloc(allocsize*sizeof(double));
        
      x->mig_source = (dlist_t **)xmalloc((size_t)opt_locus_count *
                                          sizeof(dlist_t *));
      x->mig_target = (dlist_t **)xmalloc((size_t)opt_locus_count *
                                          sizeof(dlist_t *));
      for (j = 0; j < opt_locus_count; ++j)
      {
        x->mig_source[j] = dlist_create();
        x->mig_target[j] = dlist_create();
      }
    }
    stree->mi_tbuffer = (miginfo_t **)xcalloc((size_t)opt_threads,
                                              sizeof(miginfo_t *));


    void * mem = xmalloc((size_t)(nodes_count*nodes_count)*sizeof(long) +
                         (size_t)nodes_count*sizeof(long *));
    migcount_sum = (long **)mem;
    migcount_sum[0] = (long *)(migcount_sum+nodes_count);
    memset(migcount_sum[0],0,nodes_count*sizeof(long));
    for (i = 1; i < nodes_count; ++i)
    {
      migcount_sum[i] = (long *)(migcount_sum[i-1] + nodes_count);
      memset(migcount_sum[i],0,nodes_count*sizeof(long));
    }

    stree->migcount_sum = migcount_sum;
  }

  /* TODO: Perhaps move the hx allocations into wraptree. The problem is that
     species tree cloning functions do not call wraptree and that would require
     duplicate code */
  if (opt_msci)
    network_init_hx(stree);

  /* allocate space for keeping track of coalescent events at each species tree
     node for each locus */
  stree->locus_count = (unsigned int)msa_count;
  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
  {
    snode_t * snode = stree->nodes[i];

    snode->event = (dlist_t **)xcalloc(msa_count, sizeof(dlist_t *));
    snode->event_count = (int *)xcalloc(msa_count, sizeof(int));
    snode->seqin_count = (int *)xcalloc(msa_count, sizeof(int));
    snode->gene_leaves = (unsigned int *)xcalloc(msa_count,sizeof(unsigned int));
    /* TODO: The next two allocations might not be necessary when computing
       theta analytically */
    snode->logpr_contrib = (double*)xcalloc(msa_count, sizeof(double));
    snode->old_logpr_contrib = (double *)xcalloc(msa_count, sizeof(double));

    snode->t2h = NULL;
    snode->old_t2h = NULL;
    if (!opt_est_theta)
    {
      snode->t2h = (double*)xcalloc((size_t)msa_count, sizeof(double));
      snode->old_t2h = (double*)xcalloc((size_t)msa_count, sizeof(double));
      snode->t2h_sum = 0;
      snode->event_count_sum = 0;
    }

    for (j = 0; j < stree->locus_count; ++j)
      snode->event[j] = dlist_create();
  }

  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    for (i=0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
    {
      snode_t * snode = stree->nodes[i];
      snode->brate = NULL;
      
      if (opt_msci && snode->hybrid)
      {
        if (node_is_hybridization(snode) && !snode->htau) continue;
        if (node_is_bidirection(snode) && node_is_mirror(snode)) continue;
      }
      snode->brate = (double *)xmalloc((size_t)msa_count*sizeof(double));
    }
  }

  unsigned int sum_inner = 0;
  long * locus_seqcount = (long *)xmalloc((size_t)(msa_count) * sizeof(long));
  for (i = 0; i < (unsigned int)msa_count; ++i)
    locus_seqcount[i] = msa[i]->count;

  for (i = 0; i < (unsigned int)msa_count; ++i)
    sum_inner += msa[i]->count - 1;
  stree_alloc_internals(stree, locus_seqcount, sum_inner, msa_count);

  free(locus_seqcount);
}

void stree_fini()
{
  free(__gt_nodes);
  free(__aux);
  free(__mark_count);
  free(__extra_count);
  free(__gt_nodes_index);

  if (opt_est_stree)
  {
    free(target_weight);
    free(target);
    free(moved_count);
    free(moved_space);
    free(gtarget_temp_space);
    free(gtarget_space);
    free(snode_contrib_space);
    free(snode_contrib_count);
  }

  if (ts_indicator)
    free(ts_indicator);
}

static void all_partials_recursive(gnode_t * node,
                                   unsigned int * trav_size,
                                   gnode_t ** outbuffer)
{
  if (!node->left)
    return;

  all_partials_recursive(node->left,  trav_size, outbuffer);
  all_partials_recursive(node->right, trav_size, outbuffer);

  outbuffer[*trav_size] = node;
  *trav_size = *trav_size + 1;
}

static void gtree_all_partials(gnode_t * root,
                               gnode_t ** travbuffer,
                               unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return;

  all_partials_recursive(root, trav_size, travbuffer);
}

static int cb_cmp_double_asc(const void * a, const void * b)
{
  const double * x = a;
  const double * y = b;

  if (*x > *y) return 1;

  return -1;
}

static int propose_theta_gibbs_im(stree_t * stree,
                                 gtree_t ** gtree,
                                 locus_t ** locus,
                                 snode_t * snode,
                                 long thread_index)
{
  unsigned int i, j, k, n;
  unsigned int si;
  unsigned int ti = snode->node_index;
  long lcount = 0;
  long msa_index;
  long coal_sum = 0;
  long migcount_sum = 0;
  double heredity;
  double T2h_sum = 0;
  double a1, b1;
  size_t alloc_required;
  dlist_item_t * event;
  dlist_item_t * li;
  migbuffer_t * migbuffer;

  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  for (msa_index = 0; msa_index < opt_locus_count; ++msa_index)
  {
    /* make sure migbuffer is large enough */
    alloc_required = snode->migevent_count[msa_index] +
                     snode->event_count[msa_index] +
                     stree->inner_count+1;
    migbuffer_check_and_realloc(thread_index,alloc_required);
    migbuffer = global_migbuffer_r[thread_index];

    heredity = locus[msa_index]->heredity[0];

    coal_sum += snode->event_count[msa_index];
    for (j = 0; j < total_nodes; ++j)
    {
      if (stree->nodes[j] == snode) continue;

      si = stree->nodes[j]->node_index;
      migcount_sum += gtree[msa_index]->migcount[si][ti];
    }

    /* add taus and coalescence times in sortbuffer */
    migbuffer[0].time = snode->tau;
    migbuffer[0].type = EVENT_TAU;
    j = 1;
    for (event = snode->event[msa_index]->head; event; event = event->next)
    {
      gnode_t* gnode = (gnode_t*)(event->data);
      migbuffer[j].time   = gnode->time;
      migbuffer[j++].type = EVENT_COAL;
    }

    for (li = snode->mig_source[msa_index]->head; li; li = li->next)
    {
      migevent_t * me     = (migevent_t *)(li->data);
      migbuffer[j].time   = me->time;
      migbuffer[j++].type = EVENT_MIG_SOURCE;
    }
    for (li = snode->mig_target[msa_index]->head; li; li = li->next)
    {
      migevent_t * me     = (migevent_t *)(li->data);
      migbuffer[j].time   = me->time;
      migbuffer[j++].type = EVENT_MIG_TARGET;
    }

    /* add splitting of populations */
    for (k = 0; k < snode->mb_count; ++k)
    {
      migbuffer[j++] = snode->migbuffer[k];
    }

    if (snode->parent && !snode->mb_count)
    {
      printf("\nError when processing node %s\n", snode->label);
    }
    long epoch = 0;
    double mrsum = 0;
    assert(!snode->parent || snode->mb_count);
    if (snode->mb_count)
    {
      long idx = snode->migbuffer[epoch].active_count == 1 ? 0 : msa_index;
      mrsum = snode->migbuffer[epoch].mrsum[idx];
    }

    /* TODO: Probably split the following qsort case into two:
       in case snode->parent then sort j-2 elements, otherwise
       j-1 elements.
    */

    /* if there was at least one coalescent event, sort */
    if (j > 1)
      qsort(migbuffer + 1, j - 1, sizeof(migbuffer_t), cb_migbuf_asctime);

    for (k = 1, n = snode->seqin_count[msa_index]; k < j; ++k)
    {
      double t = migbuffer[k].time - migbuffer[k - 1].time;
      T2h_sum += n * (n - 1) * t / heredity;
      if (n>0 && snode->parent)  /* no need to count migration for stree root. */
        T2h_sum += 4 * n * mrsum * t / heredity;

      if (migbuffer[k].type == EVENT_COAL || migbuffer[k].type == EVENT_MIG_SOURCE)
        --n;
      else if (migbuffer[k].type == EVENT_MIG_TARGET)
        ++n;
      else if (migbuffer[k].type == EVENT_TAU && epoch < snode->mb_count-1)
      {
        ++epoch;
        long idx = snode->migbuffer[epoch].active_count == 1 ? 0 : msa_index;
        mrsum = snode->migbuffer[epoch].mrsum[idx];
      }
    }
  }

  a1 = opt_theta_alpha + coal_sum + migcount_sum;

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
    b1 = opt_theta_beta  + T2h_sum;
  else
  {
    assert(opt_theta_dist == BPP_THETA_PRIOR_GAMMA);
    b1 = opt_theta_beta * a1 * (a1+1) /
         (opt_theta_alpha*(opt_theta_alpha+1) + T2h_sum*opt_theta_beta);
  }
  
  /* save old theta in case of gamma prior */
  double oldtheta = snode->theta;

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
    snode->theta = 1/(legacy_rndgamma(thread_index,a1) / b1);
  else
  {
    assert(opt_theta_dist == BPP_THETA_PRIOR_GAMMA);
    snode->theta = legacy_rndgamma(thread_index,a1) / b1;
  }

  /* new theta */
  double newtheta = snode->theta;

  lcount = 1;
  stree->td[0] = snode;

  if (opt_linkedtheta)
  {
    for (i=0; i < total_nodes; ++i)
    {
      snode_t * x = stree->nodes[i];

      if (x->theta >= 0 && x->has_theta && x->linked_theta == snode)
      {
        x->theta = newtheta;
        stree->td[lcount++] = x;
      }
    }
  }

  double lnacceptance = 0;
  for (i = 0; i < opt_locus_count; ++i)
  {
    lnacceptance -= gtree[i]->logpr;
    for (j = 0; j < lcount; ++j)
    {
      snode_t * x = stree->td[j];

      gtree[i]->logpr -= x->logpr_contrib[i];
      gtree_update_logprob_contrib_mig(x,
                                       stree,
                                       gtree[i],
                                       locus[i]->heredity[0],
                                       i,
                                       thread_index);
      gtree[i]->logpr += x->logpr_contrib[i];
    }
    lnacceptance += gtree[i]->logpr;
  }

  /* gamma prior acceptance-rejetion step */
  if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
  {
    /* prior ratio */
    lnacceptance += (opt_theta_alpha-1) * log(newtheta / oldtheta) -
                    opt_theta_beta*(newtheta - oldtheta);

    /* proposal ratio */
    lnacceptance += (a1-1)*log(oldtheta / newtheta) - b1*(oldtheta - newtheta);

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
      return 1;         /* accept */

    /* reject */
    for (i = 0; i < lcount; ++i)
      stree->td[i]->theta = oldtheta;

    for (i = 0; i < opt_locus_count; ++i)
    {
      for (j = 0; j < lcount; ++j)
      {
        snode_t * x = stree->td[j];

        /* revert density contributions */
        gtree[i]->logpr -= x->logpr_contrib[i];
        gtree[i]->logpr += x->old_logpr_contrib[i];
        x->logpr_contrib[i] = x->old_logpr_contrib[i];
      }
    }
    return 0;
  }

  return 1;
}

static double rndCauchy(long thread_index)
{
   /* Standard Cauchy variate, generated using inverse CDF
   */
   return tan(BPP_PI*(legacy_rndu(thread_index) - 0.5));
}

static double rndt4(long thread_index)
{
   /* Student's t_4 variate, with d.f. = 4.
      This has variance 1, and is the standard t4 variate divided by sqrt(2).
      The standard t4 variate has variance 2.
   */
   double u, v, w, c2, r2, t4, sqrt2 = 0.7071067811865475244;

   for ( ; ; ) {
      u = 2 * legacy_rndu(thread_index) - 1;
      v = 2 * legacy_rndu(thread_index) - 1;
      w = u*u + v*v;
      if (w < 1) break;
   }
   c2 = u*u / w;
   r2 = 4 / sqrt(w) - 4;
   t4 = sqrt(r2*c2);
   if (legacy_rndu(thread_index) < 0.5) t4 = -t4;

   return t4 * sqrt2;
}

static double PDFt4(double x, double m, double s)
{
   /* This t4 PDF has mean m and variance s*s.  Note that the standard t4 has variance 2*s*s.
   */
   double z = (x - m) / s, pdf;

   pdf = 3 / (4 * 1.414213562*s)*pow(1 + z*z / 2, -2.5);

   return pdf;
}

static int propose_theta_gibbs(stree_t * stree,
                               gtree_t ** gtree,
                               locus_t ** locus,
                               snode_t * snode,
                               long thread_index)
{
  unsigned int j, k, n;
  long i;
  long lcount = 0;
  long msa_index;
  long coal_sum = 0;
  double T2h_sum = 0;
  double a1, b1;
  double heredity;
  double oldtheta,newtheta;
  double xm,xs,ym,ys;
  double znew,zold,ynew,yold;
  double * sortbuffer = global_sortbuffer_r[thread_index];
  dlist_item_t * event;

  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  /* TODO: For migration we'll need to account for incoming/outgoing lineages
     as in the gene tree density function */
  assert(!opt_migration);

  /* TODO: Speed this up by using coal_event_sum variable as in NoTheta */
  /* TODO: We can store T2h in the gtree_update_logprob_contrib function */
  for (msa_index = 0; msa_index < opt_locus_count; ++msa_index)
  {
    heredity = locus[msa_index]->heredity[0];

    coal_sum += snode->event_count[msa_index];

    sortbuffer[0] = snode->tau;
    j = 1;
    for (event = snode->event[msa_index]->head; event; event = event->next)
    {
      gnode_t* gnode = (gnode_t*)(event->data);
      sortbuffer[j++] = gnode->time;
    }
    if (snode->parent)
      sortbuffer[j++] = snode->parent->tau;

    /* if there was at least one coalescent event, sort */
    if (j > 1)
      qsort(sortbuffer + 1, j - 1, sizeof(double), cb_cmp_double_asc);

    /* skip the last step in case the last value of n was supposed to be 1 */
    if ((unsigned int)(snode->seqin_count[msa_index]) == j - 1) --j;
    for (k = 1, n = snode->seqin_count[msa_index]; k < j; ++k, --n)
    {
      T2h_sum += n * (n - 1) * (sortbuffer[k] - sortbuffer[k - 1]) / heredity;
    }
  }


  /* MSC(i) model */

  a1 = opt_theta_alpha + coal_sum;

  xm = xs = ym = ys = 0;
  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
  {
    b1 = opt_theta_beta  + T2h_sum;

    if (opt_theta_move == BPP_THETA_MG_INVG)
    {
      /* debugging step -- metropolized gibbs */
      a1 *= 1.01;
      b1 *= 0.99;
    }
    else if (opt_theta_move == BPP_THETA_MG_GAMMA)
    {
      a1 -= 2;
      double a = opt_theta_alpha-2;
      double b = (opt_theta_alpha-2)*(opt_theta_alpha-1) / opt_theta_beta;
      b1 = b*a1*(a1+1)/(a*(a+1)+T2h_sum*b);
    }
    else if (opt_theta_move == BPP_THETA_MG_CAUCHY ||
             opt_theta_move == BPP_THETA_MG_T4)
    {
      /* mean and variance for logx ~ InvG(a1,b1) */
      xm = b1/(a1-1);
      xs = b1/((a1-1)*sqrt(a1-2));
      ym = log(xm) - 0.5*xs*xs/(xm*xm);
      ys = xs/xm;
    }

  }
  else
  {
    assert(opt_theta_dist == BPP_THETA_PRIOR_GAMMA);

    if (opt_theta_move == BPP_THETA_MG_INVG)
    {
      /* invgamma */
      a1 += 2;
      b1  = opt_theta_alpha*(opt_theta_alpha+1) / opt_theta_beta + T2h_sum;
    }
    else
    {
      /* gamma (also cauchy and t4) */
      b1  = opt_theta_beta * a1 * (a1+1) /
            (opt_theta_alpha*(opt_theta_alpha+1) + T2h_sum*opt_theta_beta);
    }

    if (opt_theta_move == BPP_THETA_MG_CAUCHY || opt_theta_move == BPP_THETA_MG_T4)
    {
      /* mean and variance for logx ~ Gamma(a1,b1) */
      xm = a1/b1;
      xs = sqrt(a1)/b1;
      ym = log(xm) - 0.5*xs*xs/(xm*xm);
      ys = xs/xm;
    }
  }
  
  /* save old theta in case of metropolized gibbs */
  oldtheta = snode->theta;

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
  {
    if (opt_theta_move == BPP_THETA_MG_GAMMA)
      snode->theta = legacy_rndgamma(thread_index,a1) / b1; 
    else if (opt_theta_move == BPP_THETA_GIBBS ||
             opt_theta_move == BPP_THETA_MG_INVG)
      snode->theta = 1/(legacy_rndgamma(thread_index,a1) / b1);
  }
  else
  {
    assert(opt_theta_dist == BPP_THETA_PRIOR_GAMMA);

    if (opt_theta_move == BPP_THETA_MG_GAMMA)
      snode->theta = legacy_rndgamma(thread_index,a1) / b1; 
    else if (opt_theta_move == BPP_THETA_MG_INVG)
      snode->theta = 1/(legacy_rndgamma(thread_index,a1) / b1);
  }

  yold = ynew = znew = 0;
  if (opt_theta_move == BPP_THETA_MG_CAUCHY || opt_theta_move == BPP_THETA_MG_T4)
  {
    znew = opt_theta_move == BPP_THETA_MG_CAUCHY ?
             rndCauchy(thread_index) : rndt4(thread_index);
    ynew = znew*ys + ym;
    yold = log(oldtheta);
    snode->theta = exp(ynew);
  }

  /* new theta */
  newtheta = snode->theta;

  lcount = 1;
  stree->td[0] = snode;

  if (opt_linkedtheta)
  {
    for (i=0; i < total_nodes; ++i)
    {
      snode_t * x = stree->nodes[i];

      if (x->theta >= 0 && x->has_theta && x->linked_theta == snode)
      {
        x->theta = newtheta;
        stree->td[lcount++] = x;
      }
    }
  }
  double lnacceptance = 0;
  for (i = 0; i < opt_locus_count; ++i)
  {
    lnacceptance -= gtree[i]->logpr;
    for (j = 0; j < lcount; ++j)
    {
      snode_t * x = stree->td[j];

      gtree[i]->logpr -= x->logpr_contrib[i];
      gtree_update_logprob_contrib(x, locus[i]->heredity[0],i,thread_index);
      gtree[i]->logpr += x->logpr_contrib[i];
    }
    lnacceptance += gtree[i]->logpr;
  }

  /* gamma prior acceptance-rejetion step */
  if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA ||
      opt_theta_move != BPP_THETA_GIBBS)
  {
    /* prior ratio */
    if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
      lnacceptance += (opt_theta_alpha-1) * log(newtheta / oldtheta) -
                      opt_theta_beta*(newtheta - oldtheta);
    else
      lnacceptance += (-opt_theta_alpha-1) * log(newtheta / oldtheta) -
                      opt_theta_beta*(1/newtheta - 1/oldtheta);

    /* proposal ratio */
    if (opt_theta_move == BPP_THETA_MG_GAMMA)
      lnacceptance += (a1-1) * log(oldtheta / newtheta) -
                      b1*(oldtheta - newtheta);
    else if (opt_theta_move == BPP_THETA_MG_INVG)
      lnacceptance += (-a1 - 1) * log(oldtheta/newtheta) -
                      b1*(1/oldtheta - 1/newtheta);
    else if (opt_theta_move == BPP_THETA_MG_CAUCHY)
    {
      zold = (yold - ym) / ys;
      lnacceptance += log((1 + znew*znew) / (1+zold*zold)) + ynew - yold;
    }
    else if (opt_theta_move == BPP_THETA_MG_T4)
    {
      zold = (yold - ym) / ys;
      lnacceptance += ynew - yold;
      lnacceptance += log(PDFt4(zold,0,1)/PDFt4(znew,0,1));
    }

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
      return 1;         /* accept */

    /* reject */
    for (i = 0; i < lcount; ++i)
      stree->td[i]->theta = oldtheta;

    for (i = 0; i < opt_locus_count; ++i)
    {
      for (j = 0; j < lcount; ++j)
      {
        snode_t * x = stree->td[j];

        /* revert density contributions */
        gtree[i]->logpr -= x->logpr_contrib[i];
        gtree[i]->logpr += x->old_logpr_contrib[i];
        x->logpr_contrib[i] = x->old_logpr_contrib[i];
      }
    }
    return 0;
  }

  assert(opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA &&
         opt_theta_move == BPP_THETA_GIBBS);

  return 1;
}

static int propose_theta_slide(stree_t * stree,
                               gtree_t ** gtree,
                               locus_t ** locus,
                               snode_t * snode,
                               long thread_index)
{
  long i,j;
  long lcount = 0;
  double thetaold, logthetaold;
  double thetanew, logthetanew;
  double lnacceptance = 0;
  double minv = -99;
  double maxv =  99;
  double eps = opt_finetune_theta[snode->theta_step_index];

  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  thetaold = snode->theta;

  if (!opt_exp_theta)
  {
    /* original proposal for theta */
    thetanew = thetaold + eps*legacy_rnd_symmetrical(thread_index);
    if (opt_theta_dist == BPP_THETA_PRIOR_BETA)
      thetanew = reflect(thetanew, opt_theta_min, opt_theta_max, thread_index);
    else
    {
      if (thetanew < 0)
        thetanew = -thetanew;
    }
  }
  else
  {
    /* proposal on log-scale */

    if (opt_theta_dist == BPP_THETA_PRIOR_BETA)
    {
      minv = log(opt_theta_min);
      maxv = log(opt_theta_max);
    }

    logthetaold = log(thetaold);
    logthetanew = logthetaold + eps*legacy_rnd_symmetrical(thread_index);
    logthetanew = reflect(logthetanew, minv, maxv, thread_index);

    lnacceptance = logthetanew - logthetaold;
    thetanew = exp(logthetanew);
  }

  snode->theta = thetanew;

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
    lnacceptance += (-opt_theta_alpha - 1) * log(thetanew / thetaold) -
                    opt_theta_beta*(1 / thetanew - 1 / thetaold);
  else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
    lnacceptance += (opt_theta_alpha-1) * log(thetanew / thetaold) -
                    opt_theta_beta*(thetanew - thetaold);
  else
  {
    assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
    lnacceptance += (opt_theta_p - 1) *
                    log((thetanew-opt_theta_min) / (thetaold-opt_theta_min)) +
                    (opt_theta_q - 1) *
                    log((opt_theta_max-thetanew) / (opt_theta_max-thetaold));
  }

  lcount = 1;
  stree->td[0] = snode;

  if (opt_linkedtheta)
  {
    for (i=0; i < total_nodes; ++i)
    {
      snode_t * x = stree->nodes[i];

      if (x->theta >= 0 && x->has_theta && x->linked_theta == snode)
      {
        x->theta = thetanew;
        stree->td[lcount++] = x;
      }
    }
  }

  for (i = 0; i < opt_locus_count; ++i)
  {
    /* save a copy of old logpr */
    gtree[i]->old_logpr = gtree[i]->logpr;

    for (j = 0; j < lcount; ++j)
    {
      snode_t * x = stree->td[j];

      gtree[i]->logpr -= x->logpr_contrib[i];
      if (opt_migration)
        gtree_update_logprob_contrib_mig(x,
                                         stree,
                                         gtree[i],
                                         locus[i]->heredity[0],
                                         i,
                                         thread_index);
      else
        gtree_update_logprob_contrib(x, locus[i]->heredity[0],i,thread_index);
      gtree[i]->logpr += x->logpr_contrib[i];
    }

    lnacceptance += (gtree[i]->logpr - gtree[i]->old_logpr);
  }

  if (opt_debug_theta)
    printf("[Debug] (theta) lnacceptance = %f\n", lnacceptance);

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    return 1;
  }

  /* TODO: Since the rejections are higher than acceptances, it would
     make sense not to update gtree[i]->logpr in the first place, but
     only update it when proposal is accepted */

     /* reject */
  for (i = 0; i < opt_locus_count; ++i)
    gtree[i]->logpr = gtree[i]->old_logpr;

  for (i = 0; i < lcount; ++i)
    stree->td[i]->theta = thetaold;

  for (j = 0; j < lcount; ++j)
  {
    snode_t * x = stree->td[j];

    for (i = 0; i < opt_locus_count; ++i)
      x->logpr_contrib[i] = x->old_logpr_contrib[i];
  }

  return 0;
}

void stree_propose_theta(gtree_t ** gtree,
                         locus_t ** locus,
                         stree_t * stree,
                         double * acceptvec)
{
  unsigned int i;
  snode_t * snode;
  long tip_theta_count = 0;
  long inner_theta_count = 0;

  long thread_index = 0;

  /* reset acceptance vector */
  for (i = 0; i < opt_finetune_theta_count; ++i)
    acceptvec[i] = 0;

  /* propose theta for each node */
  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
  {
    snode = stree->nodes[i];
    if (snode->theta >= 0 && snode->has_theta && !snode->linked_theta)
    {
      long move = opt_theta_move;

      /* if mixed move, then decide between slide and gibbs */
      if (move == BPP_THETA_MIXED)
      {
        move = legacy_rndu(thread_index) < opt_theta_prop ?
               BPP_THETA_SLIDE : BPP_THETA_GIBBS;
        if (move == BPP_THETA_GIBBS && opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
          move = BPP_THETA_MG_GAMMA;
      }

      if (move == BPP_THETA_SLIDE)
      {
        acceptvec[snode->theta_step_index] += propose_theta_slide(stree,
                                                                  gtree,
                                                                  locus,
                                                                  snode,
                                                                  thread_index);
      }   
      else
      {
        if (opt_migration)
          acceptvec[snode->theta_step_index] += propose_theta_gibbs_im(stree,
                                                                       gtree,
                                                                       locus,
                                                                       snode,
                                                                       thread_index);
        else
          acceptvec[snode->theta_step_index] += propose_theta_gibbs(stree,
                                                                    gtree,
                                                                    locus,
                                                                    snode,
                                                                    thread_index);
      }
      if (i < stree->tip_count)
        tip_theta_count++;
      else
        inner_theta_count++;
    }
  }

  /* normalize acceptance vector */
  switch (opt_finetune_theta_mode)
  {
    case 1:
      /* one step size for all theta */
      acceptvec[0] /= tip_theta_count+inner_theta_count;
      break;

    case 2:
      /* one step size for tip theta and one for inner theta */
      acceptvec[0] /= tip_theta_count;
      acceptvec[1] /= inner_theta_count;
      break;

    case 3:
      /* one step size for each theta, already normalized */
      break;

    default:
      fatal("Internal error");
  }
}

static int cb_cmp_pspectime(const void * a, const void * b)
{
  
  snode_t * const * x = a;
  snode_t * const * y = b;

  if ((*x)->parent->tau - (*y)->parent->tau > 0) return 1;
  return -1;
}

/* Expensive but generic function that updates the sums of incoming migration
   rates for all time segments of each populations */
void stree_update_mig_subpops(stree_t * stree, long thread_index)
{
  long i,j,k,m;
  long total_nodes = stree->tip_count+stree->inner_count;
  snode_t ** epoch;
  double tstart,tend;

  /* add the times (taus) that split the population into different "migration
   * bands". Note that x->tau is not added to x->migbuffer
  */

  /* sort splits by their parental taus */
  epoch = (snode_t **)xmalloc((size_t)total_nodes * sizeof(snode_t *));
  for (k=0, i = 0; i < total_nodes; ++i)
    if (stree->nodes[i]->parent)
      epoch[k++] = stree->nodes[i];
  assert(k == total_nodes - 1);
  qsort(epoch, total_nodes-1, sizeof(snode_t *), cb_cmp_pspectime);

  for (i = 0; i < total_nodes; ++i)
  {
    snode_t * x = stree->nodes[i];
    x->mb_count = 0;
    if (!x->parent) continue;

    migbuffer_t * xm = x->migbuffer;

    /* first create the splits */
    for (j = 0; j < total_nodes-1; ++j)
    {
      snode_t * y = epoch[j];
      assert(y->parent);
      snode_t * z = epoch[j]->parent;

      unsigned int xi = x->node_index;
      unsigned int yi = y->node_index;
      unsigned int zi = z->node_index;

      if (!z->mark[thread_index] && ((z == x->parent) ||
          ((opt_mig_bitmatrix[yi][xi] || opt_mig_bitmatrix[zi][xi]) &&
          (z->tau > x->tau && z->tau < x->parent->tau))))
      {
        xm[x->mb_count].time  = z->tau;
        xm[x->mb_count].type  = EVENT_TAU;
        xm[x->mb_count].mrsum[0] = 0;
        xm[x->mb_count].active_count = 1;
        x->mb_count++;
        z->mark[thread_index] = 1;
      }

      /* if we reached x or its sister, then break as all following populations
         in epoch have parental tau older than x->parent->tau */
      if (z == x->parent) break;
    }
    for (j = 0; j < total_nodes; ++j)
      stree->nodes[j]->mark[thread_index] = 0;

    /* fill rates */
    tstart = x->tau;
    for (k = 0; k < x->mb_count; ++k)
    {
      tend = x->migbuffer[k].time;
      for (j = 0; j < total_nodes-1; ++j)
      {
        snode_t * y = epoch[j];

        unsigned int xi = x->node_index;
        unsigned int yi = y->node_index;

        if (opt_mig_bitmatrix[yi][xi] &&
            y->tau <= tstart && y->parent->tau >= tend)
        {
          long mindex = opt_migration_matrix[yi][xi];

          /* if variables rates then we have one mrsum per locus */
          if (opt_mig_specs[mindex].am)
          {
            /* if we reached the first migration with variables rates for the
               current segment, propagate mrsum[0] to the other loci */
            if (xm[k].active_count == 1)
              for (m = 1; m < opt_locus_count; ++m)
                xm[k].mrsum[m] = xm[k].mrsum[0];

            xm[k].active_count = opt_locus_count;
          }

          if (xm[k].active_count == 1)
          {
            /* if we didn't encounter a migration with variable rates yet */
            xm[k].mrsum[0] += opt_mig_specs[mindex].M;
          }
          else
          {
            /* we have already encountered a migration with variables rates, and
               we must decide if the current migration has variable rates */
            assert(xm[k].active_count == opt_locus_count);
            if (opt_mig_specs[mindex].am)
              for (m = 0; m < xm[k].active_count; ++m)
                xm[k].mrsum[m] += opt_mig_specs[mindex].Mi[m];
            else
              for (m = 0; m < xm[k].active_count; ++m)
                xm[k].mrsum[m] += opt_mig_specs[mindex].M;
          }
        }

      }
      tstart = tend;
    }
  }
  assert(!stree->root->mark[thread_index]);
  free(epoch);
}

#if 1
/* faster function that updates the sums of incoming migration rates for all
   time segments of a single population x when M_yx changes. It assumes taus
   haven't changed, and that only M_yx has changed */
void stree_update_mig_subpops_single(stree_t * stree,
                                     snode_t * x,
                                     snode_t * y,
                                     double oldM)
{
  long k,m;
  migbuffer_t * xm = x->migbuffer;
  unsigned int xi = x->node_index;
  unsigned int yi = y->node_index;
  double tstart,tend;

  assert(opt_mig_bitmatrix[yi][xi]);
  long mindex = opt_migration_matrix[yi][xi];
  assert(!opt_mig_specs[mindex].am);

  tstart = x->tau;
  for (k = 0; k < x->mb_count; ++k)
  {
    tend = x->migbuffer[k].time;

    if (y->tau <= tstart && y->parent->tau >= tend)
    {
      if (xm[k].active_count == 1)
      {
        xm[k].mrsum[0] -= oldM;
        xm[k].mrsum[0] += opt_mig_specs[mindex].M;
      }
      else
      {
        assert(xm[k].active_count == opt_locus_count);
        for (m = 0; m < xm[k].active_count; ++m)
        {
          xm[k].mrsum[m] -= oldM;
          xm[k].mrsum[m] += opt_mig_specs[mindex].M;
        }
      }
    }
    tstart = tend;
  }
}

/* faster function that updates the sums of incoming migration rates for all
   time segments of a single population x when the migration rate for a
   particular locus i M_{yx_i} has changed. It assumes taus haven't changed. */
void stree_update_mig_subpops_single_vrates(stree_t * stree,
                                            snode_t * x,
                                            snode_t * y,
                                            long msa_index,
                                            double oldMi)
{
  long k;
  migbuffer_t * xm = x->migbuffer;
  unsigned int xi = x->node_index;
  unsigned int yi = y->node_index;
  double tstart,tend;

  assert(opt_mig_bitmatrix[yi][xi]);
  long mindex = opt_migration_matrix[yi][xi];
  assert(opt_mig_specs[mindex].am);

  tstart = x->tau;
  for (k = 0; k < x->mb_count; ++k)
  {
    tend = x->migbuffer[k].time;

    if (y->tau <= tstart && y->parent->tau >= tend)
    {
      assert(xm[k].active_count == opt_locus_count);
      xm[k].mrsum[msa_index] -= oldMi;
      xm[k].mrsum[msa_index] += opt_mig_specs[mindex].Mi[msa_index];
    }
    tstart = tend;
  }
}
#endif

#if 1
void propose_tau_update_gtrees(locus_t ** loci,
                               gtree_t ** gtree,
                               stree_t * stree,
                               snode_t * snode,
                               double oldage,
                               double minage,
                               double maxage,
                               double minfactor,
                               double maxfactor,
                               long locus_start,
                               long locus_count,
                               snode_t ** affected,
                               unsigned int paffected_count,
                               unsigned int * ret_count_above,
                               unsigned int * ret_count_below,
                               double * ret_logl_diff,
                               double * ret_logpr_diff,
                               long thread_index)
{
  long i,m;
  unsigned int j,k;
  unsigned int locus_count_above;
  unsigned int locus_count_below;
  unsigned int count_above = 0;
  unsigned int count_below = 0;
  double logl_diff = 0;
  double logpr_diff = 0;
  double logpr = 0;

  if (opt_migration && !opt_est_theta)
    fatal("Integrating out thetas not yet implemented for IM model");

  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    k = 0;
    locus_count_above = locus_count_below = 0;

    if (opt_est_theta)
      logpr = gtree[i]->logpr;

    #ifdef OLD_CODE
    gnode_t ** gt_nodesptr = __gt_nodes + offset;
    double * oldbranches = __aux + offset;
    #else
    gnode_t ** gt_nodesptr = __gt_nodes + __gt_nodes_index[i];
    #endif

    /* traverse the gene tree nodes of the three populations, find the ones
       whose ages fall within the new age interval, update their age and mark
       them. Also, count how many of them are above the old species node age,
       and how many are below. Finally, update the gene tree probabilities */
    for (j = 0; j < paffected_count; ++j)
    {
      /* process events for current population */
      if (affected[j]->seqin_count[i] > 1)
      {
        dlist_item_t * event;
        for (event = affected[j]->event[i]->head; event; event = event->next)
        {
          gnode_t * node = (gnode_t *)(event->data);
          //if (node->time < minage) continue;
          if ((node->time < minage) || (node->time > maxage)) continue;

          gt_nodesptr[k++] = node;
          node->mark = FLAG_PARTIAL_UPDATE | FLAG_BRANCH_UPDATE;
          node->old_time = node->time;

          if (node->time >= oldage && snode != stree->root)
          {
            node->time = maxage + maxfactor*(node->time - maxage);
            locus_count_above++;
          }
          else
          {
            node->time = minage + minfactor*(node->time - minage);
            locus_count_below++;
          }
        }

        if (opt_est_theta)
          logpr -= affected[j]->logpr_contrib[i];

        double xtmp;
        if (opt_migration)
          xtmp = gtree_update_logprob_contrib_mig(affected[j],
                                                  stree,
                                                  gtree[i],
                                                  loci[i]->heredity[0],
                                                  i,
                                                  thread_index);
        else
          xtmp = gtree_update_logprob_contrib(affected[j],
                                              loci[i]->heredity[0],
                                              i,
                                              thread_index);

        if (opt_est_theta)
          logpr += xtmp;
      }
    }

    /* entry i of __mark_count holds the number of marked nodes for locus i */
    __mark_count[i] = k;
    #ifdef OLD_CODE
    offset += k;
    #endif

    if (opt_est_theta)
      logpr_diff += logpr - gtree[i]->logpr;

    if (opt_est_theta)
    {
      gtree[i]->old_logpr = gtree[i]->logpr;
      gtree[i]->logpr = logpr;
    }

    count_above += locus_count_above;
    count_below += locus_count_below;

    unsigned int branch_count = k;
    gnode_t ** branchptr = gt_nodesptr;

    /* go through the list of marked nodes, and append at the end of the list
       their children (only if they are unmarked to avoid duplicates). The final
       list represents the nodes for which branch lengths must be updated */
    int extra = 0;
    for (j = 0; j < k; ++j)
    {
      gnode_t * node = gt_nodesptr[j];

      /* if root is one of the marked nodes, we must not update its branch
         length as we will receive a segfaul. Therefore, move the root to the
         beginning of the list, incremenent the pointer to the next element and
         decrease branch count */
         //if (!node->parent && j > 0)
      if (!node->parent)
      {
        if (j)
        {
          SWAP(gt_nodesptr[0], gt_nodesptr[j]);
        }
        branchptr = &(gt_nodesptr[1]);
        --branch_count;
      }

      assert(node->left);
      assert(node->right);

      if (!node->left->mark)
      {
        branchptr[branch_count++] = node->left;
        node->left->mark = FLAG_BRANCH_UPDATE;
        extra++;
      }
      if (!node->right->mark)
      {
        branchptr[branch_count++] = node->right;
        node->right->mark = FLAG_BRANCH_UPDATE;
        extra++;
      }
    }

    __extra_count[i] = extra;
    #ifdef OLD_CODE
    offset += extra;
    #endif

    /* relaxed clock */
    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      snode_t * parents[2];
      long parents_count = 0;

      if (opt_msci && snode->hybrid)
      {
        if (node_is_hybridization(snode))
        {
          if (snode->htau && snode->hybrid->htau)
          {
            /* model H3
                 *
                / \
               *   *
              / \ / \
             /   +   \
            /    |    \
            
            */

            parents[0] = snode;
            parents[1] = snode->hybrid;
            parents_count = 2;
          }
          else if (!snode->htau && !snode->hybrid->htau)
          {
            /* model H1 
                *
               / \
              /   \
             *--+--*
            /   |   \

            */

            parents[0] = snode->parent;
            parents[1] = snode->hybrid->parent;
            parents_count = 2;
          }
          else
          {
            /* model H2
                   *                  *
                  / \                / \
                 /   *      or      *   \
                /   / \            / \   \
               *---+   \          /   +---*
              /    |    \        /    |    \

            */

            if (!snode->htau)
            {
              parents[0] = snode->parent;
              parents[1] = snode->hybrid;
            }
            else
            {
              assert(!snode->hybrid->htau);
              parents[0] = snode->hybrid->parent;
              parents[1] = snode;
              parents_count = 2;
            }
          }
        }
        else
        {
          assert(node_is_bidirection(snode));
          parents[0] = snode;
          parents[1] = snode->hybrid->parent;
          parents_count = 2;
          assert(snode->right->hybrid == snode->hybrid->parent);
        }
      }
      else
      {
        assert(paffected_count == 3);
        assert(affected[1]->parent == affected[2]->parent &&
               affected[1]->parent == affected[0]);
        parents[0] = affected[0];
        parents_count = 1;
      }

      /* Note: At this point all marked nodes have the FLAG_BRANCH_UPDATE set,
         and a subset of them have the FLAG_PARTIAL_UPDATE */
      for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
      {
        gnode_t * x = gtree[i]->nodes[j];

        if (x->mark || !x->parent) continue;
        for (m = 0; m < parents_count; ++m)
        {
          snode_t * p = parents[m];

          if (stree->pptable[x->pop->node_index][p->node_index] &&
              stree->pptable[p->node_index][x->parent->pop->node_index])
          {
            if (opt_msci)
            {
              long n;
              snode_t * mnode;
              snode_t * visited;
              for (n = 0; n < stree->hybrid_count; ++n)
              {
                mnode = stree->nodes[stree->tip_count+stree->inner_count+n];
                unsigned int hindex = GET_HINDEX(stree,mnode);

                if (x->hpath[hindex] == BPP_HPATH_NONE) continue;
                if (stree->pptable[p->node_index][mnode->node_index] &&
                    stree->pptable[p->node_index][mnode->hybrid->node_index]) continue;
                if (x->hpath[hindex] == BPP_HPATH_RIGHT)
                  visited = mnode;
                else
                  visited = mnode->hybrid;

                if (stree->pptable[visited->node_index][p->node_index]) continue;

                /* otherwise, it does not pass through p */
                break;
              }
              if (n == stree->hybrid_count)
              {
                assert(!x->parent->mark || x->parent->mark == FLAG_PARTIAL_UPDATE);
                x->mark = FLAG_BRANCH_UPDATE;
                x->parent->mark = FLAG_PARTIAL_UPDATE;
                branchptr[branch_count++] = x;
                ++extra;
                break;
              }
            }
            else
            {
              assert(!x->parent->mark || x->parent->mark == FLAG_PARTIAL_UPDATE);
              x->mark = FLAG_BRANCH_UPDATE;
              x->parent->mark = FLAG_PARTIAL_UPDATE;
              branchptr[branch_count++] = x;
              ++extra;
              break;
            }
          }
        }
      }
      __extra_count[i] = extra;
    }

    /* if at least one gene tree node age was changed, we need to recompute the
       log-likelihood */
    gtree[i]->old_logl = gtree[i]->logl;
    if (opt_debug_full)
    {
      gnode_t ** tmpbuf = (gnode_t **)xmalloc((size_t)(gtree[i]->tip_count+gtree[i]->inner_count)*
                                              sizeof(gnode_t *));
                                                
      k=0;
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      {
        gnode_t * tmp = gtree[i]->nodes[j];
        if (tmp->parent)
        {
          tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                               tmp->pmatrix_index);
          tmpbuf[k++] = tmp;
        }
      }
      locus_update_matrices(loci[i],gtree[i],tmpbuf,stree,i,k);

      gtree_all_partials(gtree[i]->root,tmpbuf,&k);
      for (j = 0; j < k; ++j)
      {
        tmpbuf[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                              tmpbuf[j]->clv_index);
        if (opt_scaling)
          tmpbuf[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                      tmpbuf[j]->scaler_index);
      }
    
      locus_update_partials(loci[i],tmpbuf,k);
    
      /* compute log-likelihood */
      assert(!gtree[i]->root->parent);
      double logl = locus_root_loglikelihood(loci[i],
                                             gtree[i]->root,
                                             loci[i]->param_indices,
                                             NULL);
      logl_diff += logl - gtree[i]->logl;
      gtree[i]->logl = logl;
      free(tmpbuf);
      
    }
    else if (k+extra)   /* Nasty bug!! When having relaxed clock +extra must be there */
    {
      for (j = 0; j < branch_count; ++j)
        branchptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                      branchptr[j]->pmatrix_index);
      locus_update_matrices(loci[i], gtree[i], branchptr, stree, i, branch_count);

      /* get list of nodes for which partials must be recomputed */
      unsigned int partials_count;
      gnode_t ** partials = gtree[i]->travbuffer;
      gtree_return_partials(gtree[i]->root,
                            gtree[i]->travbuffer,
                            &partials_count);

      for (j = 0; j < partials_count; ++j)
      {
        partials[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                partials[j]->clv_index);
        if (opt_scaling)
          partials[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                        partials[j]->scaler_index);
      }

      /* update partials */
      locus_update_partials(loci[i], partials, partials_count);

      /* evaluate log-likelihood */
      double logl = locus_root_loglikelihood(loci[i],
                                             gtree[i]->root,
                                             loci[i]->param_indices,
                                             NULL);

      logl_diff += logl - gtree[i]->logl;
      gtree[i]->logl = logl;
    }

    /* Test for checking whether all gene tree nodes can be marked. It seems to
       hold */
       //if (__mark_count[i] + __extra_count[i] >= 2*loci[i]->tips-1)
       //  assert(0);

    if (opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
    {
      double new_prior_rates = lnprior_rates(gtree[i],stree,i);
      logpr_diff += new_prior_rates - gtree[i]->lnprior_rates;
      gtree[i]->old_lnprior_rates = gtree[i]->lnprior_rates;
      gtree[i]->lnprior_rates = new_prior_rates;
      #if 0
      assert(new_prior_rates > gtree[i]->old_lnprior_rates - PLL_MISC_EPSILON &&
             new_prior_rates < gtree[i]->old_lnprior_rates + PLL_MISC_EPSILON);
      #endif
    }
  }
  *ret_logpr_diff = logpr_diff;
  *ret_logl_diff = logl_diff;
  *ret_count_above = count_above;
  *ret_count_below = count_below;
}
#endif

/* update migration times during tau proposal */
static long update_migs(locus_t * locus,
                        gtree_t * gtree,
                        stree_t * stree,
                        snode_t * snode,
                        snode_t ** affected,
                        unsigned int paffected_count,
                        double oldage,
                        double minage,
                        double maxage,
                        double minfactor,
                        double maxfactor,
                        unsigned int * lca,
                        unsigned int * lcb,
                        long thread_index)
{
  long i,j,k;
  long locus_count_below = 0,locus_count_above = 0;
  double minb,maxb;

  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    if (!x->mi || x->mi->count == 0) continue;

    /* go through migration events */
    for (j = 0; j < x->mi->count; ++j)
    {
      if (x->mi->me[j].time < minage || x->mi->me[j].time > maxage) continue;

      for (k = 0; k < paffected_count; ++k)
        if (x->mi->me[j].source == affected[k] || x->mi->me[j].target == affected[k])
          break;

      if (k == paffected_count) continue;

      x->mi->me[j].old_time = x->mi->me[j].time;

      /* mark populations to update logpr */
      assert(x->mi->me[j].source->mark[thread_index]);

      if (x->mi->me[j].time >= oldage && snode != stree->root)
      {
        x->mi->me[j].time = maxage + maxfactor*(x->mi->me[j].time - maxage);
        locus_count_above++;
      }
      else
      {
        x->mi->me[j].time = minage + minfactor*(x->mi->me[j].time - minage);
        locus_count_below++;
      }
    }

    /* TODO: we now check whether there is no age conflict. This can be done
       much more efficiently */
    if (!(locus_count_below+locus_count_above)) continue;

    *lca += locus_count_above;
    *lcb += locus_count_below;

    locus_count_above = locus_count_below = 0;

    for (j = 0; j < x->mi->count; ++j)
    {
      /* determine lower bound for migration event */
      minb = MAX(x->mi->me[j].source->tau, x->mi->me[j].target->tau);
      minb = MAX(minb, j ? x->mi->me[j-1].time : x->time);

      /* determine upper bound for migration event */
      maxb = MIN(x->mi->me[j].source->parent->tau, x->mi->me[j].target->parent->tau);
      if (j < x->mi->count-1)
        maxb = MIN(x->mi->me[j+1].time,maxb);
      else if (x->parent)
        maxb = MIN(x->parent->time,maxb);

      if (x->mi->me[j].time <= minb || x->mi->me[j].time >= maxb)
        return 0;
    }
  }
  return 1;
}


void propose_tau_update_gtrees_mig(locus_t ** loci,
                                   gtree_t ** gtree,
                                   stree_t * stree,
                                   snode_t * snode,
                                   double oldage,
                                   double minage,
                                   double maxage,
                                   double minfactor,
                                   double maxfactor,
                                   long locus_start,
                                   long locus_count,
                                   snode_t ** affected,
                                   unsigned int paffected_count,
                                   long * ret_mig_reject,
                                   unsigned int * ret_count_above,
                                   unsigned int * ret_count_below,
                                   double * ret_logl_diff,
                                   double * ret_logpr_diff,
                                   long thread_index)
{
  long i;
  unsigned int j,k;
  unsigned int locus_count_above;
  unsigned int locus_count_below;
  unsigned int count_above = 0;
  unsigned int count_below = 0;
  double logl_diff = 0;
  double logpr_diff = 0;
  double logpr = 0;
  assert(opt_migration);
  assert(opt_clock == BPP_CLOCK_GLOBAL);
  if (!opt_est_theta)
    fatal("Integrating out thetas not yet implemented for IM model");

  if (opt_exp_imrb)
  {
    /* this holds for both old and new version of improved rubberband 

       Note: We will use the per-locus affected lists instead of the
       paffected_count and affected variables */
    assert(paffected_count >= 3);
  }
  else
    assert(paffected_count == 3);

  *ret_mig_reject = 0;
  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    if (opt_exp_imrb)
    {
      /* new rubberband */
      affected = gtree[i]->rb_linked;
      paffected_count = gtree[i]->rb_lcount;
      assert(paffected_count >= 3);
    }

    /* reset marks for updating population msc density contribution */
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      stree->nodes[j]->mark[thread_index] = 0;

    /* mark affected populations with more than one sequences for updating their
       density contribution */
    for (j = 0; j < paffected_count; ++j)
    {
      /* TODO:
        I think the below for loop can be made to detect overlaps between populations
        and thus detect which ones really need to update */
      affected[j]->mark[thread_index] = 1;
      /* TODO: This appears to be correct but if I delete any of the two opt_mig_bitmatrix checks, it's wrong */
      for (k = 0; k < stree->tip_count+stree->inner_count; ++k)
        if (opt_mig_bitmatrix[k][affected[j]->node_index] || opt_mig_bitmatrix[affected[j]->node_index][k])
          stree->nodes[k]->mark[thread_index] = 1;
    }

    k = 0;
    locus_count_above = locus_count_below = 0;

    if (opt_est_theta)
      logpr = gtree[i]->logpr;

    gnode_t ** gt_nodesptr = __gt_nodes + __gt_nodes_index[i];

    /* traverse the gene tree nodes of the three populations, find the ones
       whose ages fall within the new age interval, update their age and mark
       them. Also, count how many of them are above the old species node age,
       and how many are below. Finally, update the gene tree probabilities */
    for (j = 0; j < paffected_count; ++j)
    {
      /* process events for current population */
        dlist_item_t * event;
        for (event = affected[j]->event[i]->head; event; event = event->next)
        {
          gnode_t * node = (gnode_t *)(event->data);
          if ((node->time < minage) || (node->time > maxage)) continue;

          gt_nodesptr[k++] = node;
          node->mark = FLAG_PARTIAL_UPDATE | FLAG_BRANCH_UPDATE;
          node->old_time = node->time;

          if (node->time >= oldage && snode != stree->root)
          {
            node->time = maxage + maxfactor*(node->time - maxage);
            locus_count_above++;
          }
          else
          {
            node->time = minage + minfactor*(node->time - minage);
            locus_count_below++;
          }
        }
    }
    if (!update_migs(loci[i],
                     gtree[i],
                     stree,
                     snode,
                     affected,
                     paffected_count,
                     oldage,
                     minage,
                     maxage,
                     minfactor,
                     maxfactor,
                     &locus_count_above,
                     &locus_count_below,
                     thread_index))
    {
      /* conflict was caused; reject! */
      assert(!opt_exp_imrb);
      *ret_mig_reject = 1;
      return;
    }
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    {
      if (stree->nodes[j]->mark[thread_index])
      {
        if (opt_est_theta)
          logpr -= stree->nodes[j]->logpr_contrib[i];

        double xtmp = gtree_update_logprob_contrib_mig(stree->nodes[j],
                                                       stree,
                                                       gtree[i],
                                                       loci[i]->heredity[0],
                                                       i,
                                                       thread_index);
        if (opt_est_theta)
          logpr += xtmp;

        /* reset mark */
        stree->nodes[j]->mark[thread_index] = 0;
      }
    }

    /* entry i of __mark_count holds the number of marked nodes for locus i */
    __mark_count[i] = k;

    if (opt_est_theta)
      logpr_diff += logpr - gtree[i]->logpr;

    if (opt_est_theta)
    {
      gtree[i]->old_logpr = gtree[i]->logpr;
      gtree[i]->logpr = logpr;
    }

    count_above += locus_count_above;
    count_below += locus_count_below;

    unsigned int branch_count = k;
    gnode_t ** branchptr = gt_nodesptr;

    /* go through the list of marked nodes, and append at the end of the list
       their children (only if they are unmarked to avoid duplicates). The final
       list represents the nodes for which branch lengths must be updated */
    int extra = 0;
    for (j = 0; j < k; ++j)
    {
      gnode_t * node = gt_nodesptr[j];

      /* if root is one of the marked nodes, we must not update its branch
         length as we will receive a segfaul. Therefore, move the root to the
         beginning of the list, incremenent the pointer to the next element and
         decrease branch count */
         //if (!node->parent && j > 0)
      if (!node->parent)
      {
        if (j)
        {
          SWAP(gt_nodesptr[0], gt_nodesptr[j]);
        }
        branchptr = &(gt_nodesptr[1]);
        --branch_count;
      }

      assert(node->left);
      assert(node->right);

      if (!node->left->mark)
      {
        branchptr[branch_count++] = node->left;
        node->left->mark = FLAG_BRANCH_UPDATE;
        extra++;
      }
      if (!node->right->mark)
      {
        branchptr[branch_count++] = node->right;
        node->right->mark = FLAG_BRANCH_UPDATE;
        extra++;
      }
    }

    __extra_count[i] = extra;

    /* if at least one gene tree node age was changed, we need to recompute the
       log-likelihood */
    gtree[i]->old_logl = gtree[i]->logl;
    if (opt_debug_full)
    {
      gnode_t ** tmpbuf = (gnode_t **)xmalloc((size_t)(gtree[i]->tip_count+gtree[i]->inner_count)*
                                              sizeof(gnode_t *));
                                                
      k=0;
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      {
        gnode_t * tmp = gtree[i]->nodes[j];
        if (tmp->parent)
        {
          tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                               tmp->pmatrix_index);
          tmpbuf[k++] = tmp;
        }
      }
      locus_update_matrices(loci[i],gtree[i],tmpbuf,stree,i,k);

      gtree_all_partials(gtree[i]->root,tmpbuf,&k);
      for (j = 0; j < k; ++j)
      {
        tmpbuf[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                              tmpbuf[j]->clv_index);
        if (opt_scaling)
          tmpbuf[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                      tmpbuf[j]->scaler_index);
      }
    
      locus_update_partials(loci[i],tmpbuf,k);
    
      /* compute log-likelihood */
      assert(!gtree[i]->root->parent);
      double logl = locus_root_loglikelihood(loci[i],
                                             gtree[i]->root,
                                             loci[i]->param_indices,
                                             NULL);
      logl_diff += logl - gtree[i]->logl;
      gtree[i]->logl = logl;
      free(tmpbuf);
      
    }
    else if (k+extra)   /* Nasty bug!! When having relaxed clock +extra must be there */
    {
      for (j = 0; j < branch_count; ++j)
        branchptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                      branchptr[j]->pmatrix_index);
      locus_update_matrices(loci[i], gtree[i], branchptr, stree, i, branch_count);

      /* get list of nodes for which partials must be recomputed */
      unsigned int partials_count;
      gnode_t ** partials = gtree[i]->travbuffer;
      gtree_return_partials(gtree[i]->root,
                            gtree[i]->travbuffer,
                            &partials_count);

      for (j = 0; j < partials_count; ++j)
      {
        partials[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                partials[j]->clv_index);
        if (opt_scaling)
          partials[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                        partials[j]->scaler_index);
      }

      /* update partials */
      locus_update_partials(loci[i], partials, partials_count);

      /* evaluate log-likelihood */
      double logl = locus_root_loglikelihood(loci[i],
                                             gtree[i]->root,
                                             loci[i]->param_indices,
                                             NULL);

      logl_diff += logl - gtree[i]->logl;
      gtree[i]->logl = logl;
    }

  }
  *ret_logpr_diff = logpr_diff;
  *ret_logl_diff = logl_diff;
  *ret_count_above = count_above;
  *ret_count_below = count_below;
}

static long propose_tau(locus_t ** loci,
                        snode_t * snode,
                        gtree_t ** gtree,
                        stree_t * stree,
                        unsigned int candidate_count,
                        long thread_index)
{
  unsigned int i, j, k;
  #ifdef OLD_CODE
  unsigned int offset = 0;
  #endif
  int theta_method = 2;   /* how we change theta */
  long accepted = 0;
  double oldage, newage;
  double minage = 0, maxage = 999;
  double lnacceptance = 0;
  double minfactor, maxfactor, thetafactor;
  double oldtheta = 0;
  double logpr = 0;
  double logpr_diff = 0;
  double logl_diff = 0;

  unsigned int count_above = 0;
  unsigned int count_below = 0;

  unsigned int paffected_count;        /* number of affected populations */

  oldage = snode->tau;

  /* TODO: Implement a generic version */
  if (opt_linkedtheta)
    theta_method = 0;

  /* compute minage and maxage bounds */
  if (opt_msci && snode->hybrid)
  {
    if (node_is_hybridization(snode))
    {
      /* hybridization */
      assert(!node_is_mirror(snode));
      assert(snode->left && !snode->right);

      minage = snode->left->tau;
      if (!snode->htau)
      {
        snode_t * sibling = (snode->parent->left == snode) ?
                              snode->parent->right : snode->parent->left;

        minage = MAX(minage,sibling->tau);
      }
      if (!snode->hybrid->htau)
      {
        snode_t * sibling = (snode->hybrid->parent->left == snode->hybrid) ?
                              snode->hybrid->parent->right : snode->hybrid->parent->left;

        minage = MAX(minage,sibling->tau);
      }
    }
    else
    {
      /* bidirectional introgression */
      assert(node_is_bidirection(snode));
      assert(!node_is_mirror(snode));
      assert(snode->left && snode->right);
      assert(snode->right->hybrid && snode->right->hybrid->left);
      assert(snode->htau && !snode->hybrid->htau &&
             snode->hybrid->parent->htau && !snode->right->htau);

      minage = MAX(snode->left->tau,snode->right->hybrid->left->tau);
    }
  }
  else
  {
    if (snode->left)
      minage = MAX(snode->left->tau, snode->right->tau);
  }

  if (opt_msci && snode->hybrid)
  {
    if (node_is_hybridization(snode))
    {
      assert(!node_is_mirror(snode));
      assert(snode->htau || snode->parent->parent);
      assert(snode->hybrid->htau || snode->hybrid->parent->parent);

      double p1age = snode->htau ?
                       snode->parent->tau : snode->parent->parent->tau;
      double p2age = snode->hybrid->htau ?
                       snode->hybrid->parent->tau : snode->hybrid->parent->parent->tau;

      maxage = MIN(p1age,p2age);
    }
    else
    {
      assert(node_is_bidirection(snode));
      assert(!node_is_mirror(snode));
      
      double p1age = snode->parent->tau;
      double p2age = snode->right->hybrid->parent->tau;

      maxage = MIN(p1age,p2age);
    }
  }
  else
  {
    if (snode->parent)
      maxage = snode->parent->tau;
  }

  /* propose new tau */
  newage = oldage + opt_finetune_tau * legacy_rnd_symmetrical(thread_index);
  newage = reflect(newage, minage, maxage, thread_index);
  snode->tau = newage;

  if (opt_msci && snode->hybrid)
  {
    if (node_is_hybridization(snode))
    {
      snode->hybrid->tau = snode->tau;

      if (!snode->htau)
        snode->parent->tau = snode->tau;
      if (!snode->hybrid->htau)
        snode->hybrid->parent->tau = snode->tau;
    }
    else
    {
      assert(node_is_bidirection(snode));
      assert(!node_is_mirror(snode));

      snode->hybrid->tau        = snode->tau;
      snode->right->tau         = snode->tau;
      snode->right->hybrid->tau = snode->tau;
    }
  }

  /* compute factors for multiplying associated gene tree nodes ages */
  minfactor = (newage - minage) / (oldage - minage);
  maxfactor = (newage - maxage) / (oldage - maxage);

  /* if we are dealing with the root population, add the following factor to
     the acceptance ratio */
  if (snode == stree->root)
  {
    if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
      lnacceptance = (-opt_tau_alpha - 1 - candidate_count + 1) *
                     log(newage / oldage) - opt_tau_beta*(1/newage - 1/oldage);
    else
      lnacceptance = (opt_tau_alpha-1 - candidate_count + 1) *
                     log(newage/oldage) - opt_tau_beta*(newage-oldage);
  }

  /* change theta as well */
  if (opt_est_theta && theta_method)
  {
    /* network code might have snode without a theta (i.e snode->theta == -1)
       when snode is a hybrid node and its two parents have no tau */
    if (snode->theta != -1)
    {
      oldtheta = snode->theta;
      if (theta_method == 1)
        thetafactor = newage / oldage;
      else if (theta_method == 2)
        thetafactor = (newage - minage) / (oldage - minage);
      else
        assert(0);

      /* check is needed for network code. When the two parents of snode have
         no tau, then snode has no theta */
      if (snode->theta != -1)
        snode->theta = oldtheta / thetafactor;
      if (opt_msci && snode->hybrid)
      {
        /* TODO : Need to update the hybrid as well */
        //snode->hybrid->theta = 
      }

      if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
        lnacceptance += -log(thetafactor) +
                        (-opt_theta_alpha - 1) * log(snode->theta / oldtheta) -
                        opt_theta_beta*(1 / snode->theta - 1 / oldtheta);
      else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
        lnacceptance += -log(thetafactor) +
                        (opt_theta_alpha-1) * log(snode->theta / oldtheta) -
                        opt_theta_beta*(snode->theta - oldtheta);
      else
      {
        assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
        lnacceptance += -log(thetafactor) + (opt_theta_p - 1) *
                       log((snode->theta-opt_theta_min) / (oldtheta-opt_theta_min)) +
                       (opt_theta_q - 1) *
                       log((opt_theta_max-snode->theta) / (opt_theta_max-oldtheta));
      }
    }
  }

  snode_t * affected[7];
  double old_logpr_contrib[7] = {0,0,0,0,0,0,0};

  if (opt_msci && snode->hybrid)
  {
    if (node_is_hybridization(snode))
    {
      if (snode->htau && snode->hybrid->htau)
      {
        /* model H3
             *
            / \
           *   *
          / \ / \
         /   +   \
        /    |    \
        
        */

        affected[0] = snode->left;
        affected[1] = snode;
        affected[2] = snode->hybrid;
        paffected_count = 3;
      }
      else if (!snode->htau && !snode->hybrid->htau)
      {
        /* model H1 
            *
           / \
          /   \
         *--+--*
        /   |   \

        */

        paffected_count = 0;
        for (i = 0; i < stree->locus_count; ++i)
        {
          assert(snode->event_count[i] == 0);
          assert(snode->hybrid->event_count[i] == 0);
        }

        affected[paffected_count++] = snode->parent;
        affected[paffected_count++] = snode->hybrid->parent;
        affected[paffected_count++] = snode->parent->left;
        affected[paffected_count++] = snode->parent->right;
        affected[paffected_count++] = snode->hybrid->parent->left;
        affected[paffected_count++] = snode->hybrid->parent->right;
        affected[paffected_count++] = snode->left;
      }
      else
      {
        /* model H2
               *                  *
              / \                / \
             /   *      or      *   \
            /   / \            / \   \
           *---+   \          /   +---*
          /    |    \        /    |    \

        */

        assert((!snode->htau && snode->hybrid->htau) ||
               (snode->htau && !snode->hybrid->htau));

        paffected_count = 0;

        /* assertions */
        if (!snode->htau)
          for (i = 0; i < stree->locus_count; ++i)
            assert(snode->event_count[i] == 0);
        if (!snode->hybrid->htau)
          for (i = 0; i < stree->locus_count; ++i)
            assert(snode->hybrid->event_count[i] == 0);

        if (!snode->htau)
        {
          affected[paffected_count++] = snode->parent;
          affected[paffected_count++] = snode->parent->left;
          affected[paffected_count++] = snode->parent->right;
          affected[paffected_count++] = snode->hybrid;
        }
        else
        {
          assert(!snode->hybrid->htau);

          affected[paffected_count++] = snode->hybrid->parent;
          affected[paffected_count++] = snode->hybrid->parent->left;
          affected[paffected_count++] = snode->hybrid->parent->right;
          affected[paffected_count++] = snode;
        }

        affected[paffected_count++] = snode->left;
      }
    }
    else
    {
      assert(node_is_bidirection(snode));
      for (i = 0; i < stree->locus_count; ++i)
      {
        assert(snode->hybrid->event_count[i] == 0);
        assert(snode->right->event_count[i] == 0);
      }

      paffected_count = 0;
      affected[paffected_count++] = snode;
      affected[paffected_count++] = snode->hybrid;
      affected[paffected_count++] = snode->right;
      affected[paffected_count++] = snode->right->hybrid;
      affected[paffected_count++] = snode->right->hybrid->left;
      affected[paffected_count++] = snode->left;
    }
  }
  else
  {
    paffected_count = 3;
    affected[0] = snode;
    affected[1] = snode->left;
    affected[2] = snode->right;
  }

  if (!opt_est_theta)
  {
    for (i = 0; i < paffected_count; ++i)
      old_logpr_contrib[i] = affected[i]->notheta_logpr_contrib;

    logpr = stree->notheta_logpr;
    for (j = 0; j < paffected_count; ++j)
      logpr -= affected[j]->notheta_logpr_contrib;
  }
  /* Note: we use the serial version when thetas are integrated out */
  if (opt_est_theta && opt_threads > 1)
  {
    thread_data_t tp;
    tp.locus = loci; tp.gtree = gtree; tp.stree = stree;
    tp.snode = snode;
    tp.oldage = oldage;
    tp.minage = minage;
    tp.maxage = maxage;
    tp.minfactor = minfactor;
    tp.maxfactor = maxfactor;
    tp.affected = affected;
    tp.paffected_count = paffected_count;
    tp.logl_diff = logl_diff;
    tp.logpr_diff = logpr_diff;
    threads_wakeup(THREAD_WORK_TAU,&tp);

    count_above = tp.count_above;
    count_below = tp.count_below;
    logl_diff = tp.logl_diff;
    logpr_diff = tp.logpr_diff;
  }
  else
    propose_tau_update_gtrees(loci,
                              gtree,
                              stree,
                              snode,
                              oldage,
                              minage,
                              maxage,
                              minfactor,
                              maxfactor,
                              0,
                              stree->locus_count,
                              affected,
                              paffected_count,
                              &count_above,
                              &count_below,
                              &logl_diff,
                              &logpr_diff,
                              0);

  if (!opt_est_theta)
  {
    for (j = 0; j < paffected_count; ++j)
      logpr += affected[j]->notheta_logpr_contrib;

    /* Note: At this point logpr_diff must be 0 unless
       we have a correlated clock with LN */
    #if 0
    if (opt_clock != BPP_CLOCK_CORR || opt_rate_prior != BPP_BRATE_PRIOR_LOGNORMAL)
      assert(logpr_diff == 0);
    #endif

    logpr_diff += logpr - stree->notheta_logpr;
    stree->notheta_old_logpr = stree->notheta_logpr;
    stree->notheta_logpr = logpr;
  }

  lnacceptance += logpr_diff + logl_diff + count_below*log(minfactor) +
                  count_above*log(maxfactor);

  if (opt_debug_tau)
    printf("[Debug] (tau) lnacceptance = %f\n", lnacceptance);

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */
    accepted++;

    for (i = 0; i < stree->locus_count; ++i)
    {
      k = __mark_count[i];

      /* The commented loop is probably sufficient */
      //gnode_t ** gt_nodesptr = __gt_nodes + __gt_nodes_index[i];
      //for (j = 0; j < k + __extra_count[i]; ++j)
      //  gt_nodesptr[j]->mark = 0;
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        gtree[i]->nodes[j]->mark = 0;
    }
  }
  else
  {
    /* rejected */
    snode->tau = oldage;
    if (opt_msci && snode->hybrid)
    {
      if (node_is_hybridization(snode))
      {
        snode->hybrid->tau = snode->tau;

        if (!snode->htau)
          snode->parent->tau = snode->tau;
        if (!snode->hybrid->htau)
          snode->hybrid->parent->tau = snode->tau;
      }
      else
      {
        assert(node_is_bidirection(snode));

        snode->hybrid->tau        = snode->tau;
        snode->right->tau         = snode->tau;
        snode->right->hybrid->tau = snode->tau;
      }
    }

    if (opt_est_theta && theta_method)
    {
      /* the -1 check is due to network code (when snode is hybridization with
         its two parents having no tau) */
      if (snode->theta != -1)
      {
        snode->theta = oldtheta;
        if (opt_msci && snode->hybrid)
        {
          /* TODO: Need to update hybrid theta */
        }
      }
    }

    for (i = 0; i < stree->locus_count; ++i)
    {
      k = __mark_count[i];
      gnode_t ** gt_nodesptr = __gt_nodes + __gt_nodes_index[i];

      /* restore gene tree node ages */
      for (j = 0; j < k; ++j)
        gt_nodesptr[j]->time = gt_nodesptr[j]->old_time;

      if (opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
      {
        gtree[i]->lnprior_rates = gtree[i]->old_lnprior_rates;
      }

      /* restore logpr contributions */
      if (opt_est_theta)
      {
        for (j = 0; j < paffected_count; ++j)
        {
          if (affected[j]->seqin_count[i] > 1)
          {
            if (opt_migration)
              gtree_update_logprob_contrib_mig(affected[j],
                                               stree,
                                               gtree[i],
                                               loci[i]->heredity[0],
                                               i,
                                               thread_index);
            else
              gtree_update_logprob_contrib(affected[j],
                                           loci[i]->heredity[0],
                                           i,
                                           thread_index);
          }
        }
      }
      else
      {
        for (j = 0; j < paffected_count; ++j)
          if (affected[j]->seqin_count[i] > 1)
            logprob_revert_notheta(affected[j], i);
      }

      /* get the list of nodes for which CLVs must be reverted, i.e. all marked
         nodes and all nodes whose left or right subtree has at least one marked
         node */
      unsigned int partials_count;
      gnode_t ** partials = gtree[i]->travbuffer;
      gtree_return_partials(gtree[i]->root,
                            gtree[i]->travbuffer,
                            &partials_count);

      /* revert CLV indices */
      if (opt_debug_full)
      {
        for (j = gtree[i]->tip_count; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
        {
          gnode_t * tmp = gtree[i]->nodes[j];
          tmp->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                          tmp->clv_index);
          if (opt_scaling)
            tmp->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                  tmp->scaler_index);
        }
      }
      else
      {
        for (j = 0; j < partials_count; ++j)
        {
          partials[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                  partials[j]->clv_index);
          if (opt_scaling)
            partials[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                      partials[j]->scaler_index);
        }
      }

      /* un-mark nodes */
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        gtree[i]->nodes[j]->mark = 0;

      /* restore branch lengths and pmatrices */
      int matrix_updates = __mark_count[i] + __extra_count[i];
      if (matrix_updates)
      {
        if (!gt_nodesptr[0]->parent)
        {
          --matrix_updates;
          gt_nodesptr++;
        }
        if (opt_debug_full)
        {
          for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
          {
            gnode_t * tmp = gtree[i]->nodes[j];
            if (tmp->parent)
              tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                   tmp->pmatrix_index);
          }
        }
        else
        {
          for (j = 0; j < matrix_updates; ++j)
            gt_nodesptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                            gt_nodesptr[j]->pmatrix_index);
        }
      }

      /* restore gene tree log-likelihood */
      gtree[i]->logl = gtree[i]->old_logl;

      /* restore gene tree log probability */
      if (opt_est_theta)
        gtree[i]->logpr = gtree[i]->old_logpr;
    }
    if (!opt_est_theta)
    {
      for (j = 0; j < paffected_count; ++j)
        affected[j]->notheta_logpr_contrib = old_logpr_contrib[j];

      stree->notheta_logpr = stree->notheta_old_logpr;
    }
  }
  return accepted;
}

#define RB_MARK 1
#define RB_INLIST 2
static long getlinkedpops(stree_t * stree,
                          gtree_t ** gtree_list,
                          snode_t * x,
                          double * bounds,
                          long round,
                          snode_t ** linked)
{
  long i,j,k,m;
  long lcount = 0;
  long mfound;
  long loc_lcount = 0;
  double tl = bounds[0];
  double tu = bounds[1];
  migevent_t * me;
  dlist_item_t * dli;
  snode_t ** loc_linked;
  const static long thread_index_zero = 0;

  snode_t * startx = x;

  /* calculate per locus lists of linked populations to x */
  for (m = 0; m < opt_locus_count; ++m)
  {
    x = startx;
    gtree_t * gtree = gtree_list[m];

    /* place one of the three affected nodes in the corresponding slot */
    gtree->rb_linked[round] = x;

    /* if it's the first round of calling this function zero-out the number of
       linked nodes for the current population. We set it to 3 since the first
       three slots are reserved for the focal node and its two children */
    if (round == 0)
      gtree->rb_lcount = 3;

    loc_linked = gtree->rb_linked + gtree->rb_lcount;
    loc_lcount = 0;

    /* part "auto-cleaning" the list from the previous MCMC step */
    loc_linked[0] = NULL;

    /* mark nodes already in the per locus list of linked pops */
    for (i = 3; i < gtree->rb_lcount; ++i)
      gtree->rb_linked[i]->mark[thread_index_zero] |= RB_MARK;

    k = 0;
    while (x)
    {
      /* get index of population for which we are looking for its linked pops */
      j = x->node_index;

      /* look for linked pops */
      for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      {
        /* skip populations already identified as linked */
        if (stree->nodes[i]->mark[thread_index_zero] & RB_MARK) continue;

        #if 0
        /* skip if no migration between the two populations */
        if (!(gtree->migcount[i][j] || gtree->migcount[j][i])) continue;
        #endif

        mfound = 0;
        if (gtree->migcount[i][j])
        {
          for (dli = x->mig_source[m]->head; dli; dli = dli->next)
          {
            me = (migevent_t *)(dli->data);
            if (me->target == stree->nodes[i] && me->time > tl && me->time < tu)
            {
              mfound = 1;
              break;
            }
          }
        }

        if (!mfound && gtree->migcount[j][i])
        {
          for (dli = x->mig_target[m]->head; dli; dli = dli->next)
          {
            me = (migevent_t *)(dli->data);
            if (me->source == stree->nodes[i] && me->time > tl && me->time < tu)
            {
              mfound = 1;
              break;
            }
          }
        }

        if (!mfound) continue;

        loc_linked[loc_lcount++] = stree->nodes[i];
        loc_linked[loc_lcount] = NULL;               /* part of auto-cleaning */

        stree->nodes[i]->mark[thread_index_zero] |= RB_MARK;
      }

      /* move to the next population in the list of linked pops and repeat */
      x = loc_linked[k++];
    }

    /* update number of linked populations for current locus */
    gtree->rb_lcount += loc_lcount;

    /* unmark the nodes already in the per locus list of linked pops */
    for (i = 3; i < gtree->rb_lcount; ++i)
      gtree->rb_linked[i]->mark[thread_index_zero] &= ~RB_MARK;

    /* update the list of linked populations */
    for (i = 0; i < loc_lcount; ++i)
    {
      if (loc_linked[i]->mark[thread_index_zero] & RB_INLIST) continue;

      loc_linked[i]->mark[thread_index_zero] |= RB_INLIST;
      linked[lcount++] = loc_linked[i];
    }
  }

  return lcount;
}

static long rb_bounds(stree_t * stree,
                      gtree_t ** gtree_list,
                      snode_t * x,
                      snode_t ** linked,
                      double * bounds)
{
  long i,j,m;
  long lcount = 0;
  double tl, tu;

  const static long thread_index_zero = 0;

  unsigned int total_nodes = stree->tip_count+stree->inner_count;

  /* sanity check -- remove */
  for (i = 0; i < total_nodes; ++i)
    assert(stree->nodes[i]->mark[thread_index_zero] == 0);

  /* 1. Set initial bounds */
  tu = (x->parent) ? x->parent->tau : 999;
  tl = MAX(x->left->tau, x->right->tau);
  
  bounds[0] = tl; bounds[1] = tu;

  /* 2. Find populations linked to X,V,W - both directly and indirectly linked */
  x->mark[thread_index_zero] = RB_MARK | RB_INLIST;
  x->left->mark[thread_index_zero] = RB_MARK | RB_INLIST;
  x->right->mark[thread_index_zero] = RB_MARK | RB_INLIST;
  linked[0] = x; linked[1] = x->left; linked[2] = x->right;
  lcount = 3;

  lcount += getlinkedpops(stree,gtree_list,x,       bounds,0,linked+lcount);
  lcount += getlinkedpops(stree,gtree_list,x->left, bounds,1,linked+lcount);
  lcount += getlinkedpops(stree,gtree_list,x->right,bounds,2,linked+lcount);

  /* cleans marks also on x, x->left and x->right */
  for (i = 0; i < lcount; ++i)
    linked[i]->mark[thread_index_zero] = 0;

  /* 3. and 4. merged -- get final upper and lower bounds. Skip first 3 nodes */
  for (i = 3; i < lcount; ++i)
  {
    if (linked[i]->tau < x->tau && linked[i]->parent->tau > x->tau)
    {
      if (linked[i]->parent->tau < tu)
        tu = linked[i]->parent->tau;
      if (linked[i]->tau > tl)
        tl = linked[i]->tau;
    }
    else if (linked[i]->parent->tau < x->tau)
    {
      /* t_a < t_X */
      if (linked[i]->parent->tau > tl)
        tl = linked[i]->parent->tau;
    }
    else
    {
      /* t_b > t_X */
      assert(linked[i]->tau > x->tau);
      if (linked[i]->tau < tu)
        tu = linked[i]->tau;
    }
  }
  bounds[0] = tl; bounds[1] = tu;

  /* 5. Delete unaffected linked populations. Skip first three nodes
        NOTE: j holds the last empty position */
  for (j=3, i=3; i < lcount; ++i)
  {
    if (linked[i]->tau > tu || linked[i]->parent->tau < tl)
    {
      linked[i] = NULL;
    }
    else
    {
      if (i != j)
      {
        linked[j] = linked[i];
        linked[i] = NULL;
      }
      ++j;
    }
  }
  lcount = j;

  /* update per locus lists */
  for (m = 0; m < opt_locus_count; ++m)
  {
    snode_t ** loc_linked = gtree_list[m]->rb_linked;
    long loc_lcount = gtree_list[m]->rb_lcount;

    for (j=3, i=3; i < loc_lcount; ++i)
    {
      if (loc_linked[i]->tau > tu || loc_linked[i]->parent->tau < tl)
      {
        loc_linked[i] = NULL;
      }
      else
      {
        if (i != j)
        {
          loc_linked[j] = loc_linked[i];
          loc_linked[i] = NULL;
        }
        ++j;
      }
    }
    gtree_list[m]->rb_lcount = j;
  }

  return lcount;
}

static long propose_tau_mig(locus_t ** loci,
                            snode_t * snode,
                            gtree_t ** gtree,
                            stree_t * stree,
                            unsigned int candidate_count,
                            long thread_index)
{
  unsigned int i, j, k;

  int theta_method = 2;   /* how we change theta */
  long accepted = 0;
  long mig_reject;        /* reject move due to conflict in gene tree ages */
  double oldage, newage;
  double minage = 0, maxage = 999;
  double lnacceptance = 0;
  double minfactor, maxfactor, thetafactor;
  double oldtheta = 0;
  double logpr = 0;
  double logpr_diff = 0;
  double logl_diff = 0;
  snode_t ** affected;
  snode_t * affected3[3];

  unsigned int count_above = 0;
  unsigned int count_below = 0;
  unsigned int paffected_count;        /* number of affected populations */

  /* TODO: Implement a generic version */
  if (opt_linkedtheta)
    theta_method = 0;

  oldage = snode->tau;

  /* compute minage and maxage bounds */
  if (snode->left)
    minage = MAX(snode->left->tau, snode->right->tau);

  if (snode->parent)
    maxage = snode->parent->tau;

  if (opt_exp_imrb)
  {
    double bounds[2];
    affected = (snode_t **)xcalloc((size_t)stree->tip_count+stree->inner_count+1, sizeof(snode_t *));

    paffected_count = rb_bounds(stree,gtree,snode,affected,bounds);

    minage = bounds[0];
    maxage = bounds[1];
  }
  else
  {
    affected = affected3;
    paffected_count = 3;
    affected[0] = snode;
    affected[1] = snode->left;
    affected[2] = snode->right;
  }

  /* propose new tau */
  newage = oldage + opt_finetune_tau * legacy_rnd_symmetrical(thread_index);
  newage = reflect(newage, minage, maxage, thread_index);
  snode->tau = newage;

  stree_update_mig_subpops(stree,thread_index);

  /* compute factors for multiplying associated gene tree nodes ages */
  minfactor = (newage - minage) / (oldage - minage);
  maxfactor = (newage - maxage) / (oldage - maxage);

  /* if we are dealing with the root population, add the following factor to
     the acceptance ratio */
  if (snode == stree->root)
  {
    if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
      lnacceptance = (-opt_tau_alpha - 1 - candidate_count + 1) *
                     log(newage / oldage) - opt_tau_beta*(1/newage - 1/oldage);
    else
      lnacceptance = (opt_tau_alpha-1 - candidate_count + 1) *
                     log(newage/oldage) - opt_tau_beta*(newage-oldage);
  }

  /* change theta as well */
  if (theta_method && opt_est_theta)
  {
    /* network code might have snode without a theta (i.e snode->theta == -1)
       when snode is a hybrid node and its two parents have no tau */
    if (snode->theta != -1)
    {
      oldtheta = snode->theta;
      if (theta_method == 1)
        thetafactor = newage / oldage;
      else if (theta_method == 2)
        thetafactor = (newage - minage) / (oldage - minage);
      else
        assert(0);

      /* check is needed for network code. When the two parents of snode have
         no tau, then snode has no theta */
      if (snode->theta != -1)
        snode->theta = oldtheta / thetafactor;

      if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
        lnacceptance += -log(thetafactor) +
                        (-opt_theta_alpha - 1) * log(snode->theta / oldtheta) -
                        opt_theta_beta*(1 / snode->theta - 1 / oldtheta);
      else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
        lnacceptance += -log(thetafactor) +
                        (opt_theta_alpha-1) * log(snode->theta / oldtheta) -
                        opt_theta_beta*(snode->theta - oldtheta);
      else
      {
        assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
        lnacceptance += -log(thetafactor) + (opt_theta_p - 1) *
                       log((snode->theta-opt_theta_min) / (oldtheta-opt_theta_min)) +
                       (opt_theta_q - 1) *
                       log((opt_theta_max-snode->theta) / (opt_theta_max-oldtheta));
      }
    }
  }

  if (!opt_est_theta)
  {
    logpr = stree->notheta_logpr;
    for (j = 0; j < paffected_count; ++j)
      logpr -= affected[j]->notheta_logpr_contrib;
  }
  #if 1
  /* change to if 0 to disable gtrees. Note: we use the serial version when
     thetas are integrated out */
  mig_reject = 0;
  if (opt_est_theta && opt_threads > 1)
  {
    thread_data_t tp;
    tp.locus = loci; tp.gtree = gtree; tp.stree = stree;
    tp.snode = snode;
    tp.oldage = oldage;
    tp.minage = minage;
    tp.maxage = maxage;
    tp.minfactor = minfactor;
    tp.maxfactor = maxfactor;
    tp.affected = affected;
    tp.paffected_count = paffected_count;
    tp.logl_diff = logl_diff;
    tp.logpr_diff = logpr_diff;
    threads_wakeup(THREAD_WORK_TAU_MIG,&tp);

    mig_reject  = tp.mig_reject;
    count_above = tp.count_above;
    count_below = tp.count_below;
    logl_diff = tp.logl_diff;
    logpr_diff = tp.logpr_diff;
  }
  else
    propose_tau_update_gtrees_mig(loci,
                                  gtree,
                                  stree,
                                  snode,
                                  oldage,
                                  minage,
                                  maxage,
                                  minfactor,
                                  maxfactor,
                                  0,
                                  stree->locus_count,
                                  affected,
                                  paffected_count,
                                  &mig_reject,
                                  &count_above,
                                  &count_below,
                                  &logl_diff,
                                  &logpr_diff,
                                  0);
  
  if (mig_reject)
    return -1;
  #endif

  if (!opt_est_theta)
  {
    for (j = 0; j < paffected_count; ++j)
      logpr += affected[j]->notheta_logpr_contrib;

    /* Note: At this point logpr_diff must be 0 unless
       we have a correlated clock with LN */
    #if 0
    if (opt_clock != BPP_CLOCK_CORR || opt_rate_prior != BPP_BRATE_PRIOR_LOGNORMAL)
      assert(logpr_diff == 0);
    #endif

    logpr_diff += logpr - stree->notheta_logpr;
    stree->notheta_old_logpr = stree->notheta_logpr;
    stree->notheta_logpr = logpr;
  }

  if (opt_exp_imrb)
    free(affected);

  #if 1
  /* Pseudo prior implementation */

  /* if no pseudoprior given then pseudoprior = prior and terms cancel out */
  if (opt_pseudop_exist)
  {
    snode_t * nodes[3];
    double atau[3], old_atau[3];
    double ptau[3], old_ptau[3];
    double amodel;
    double bmodel;
    migspec_t * spec;
      
    /* collect the tau boundaries for each of the two (if snode is root) or
       three  "migration bands" affected by the tau change. See diagram below */
    k = 0;
    if (snode->parent)
    {
      nodes[k] = snode;
      atau[k] = newage;             old_atau[k]   = oldage;
      ptau[k] = snode->parent->tau; old_ptau[k++] = snode->parent->tau;
    }

    nodes[k] = snode->left;
    atau[k] = snode->left->tau;     old_atau[k]   = snode->left->tau;
    ptau[k] = newage;               old_ptau[k++] = oldage;

    nodes[k] = snode->right;
    atau[k] = snode->right->tau;    old_atau[k]   = snode->right->tau;
    ptau[k] = newage;               old_ptau[k++] = oldage;

    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    {
      if (!stree->nodes[i]->parent) continue;

      unsigned int yi = stree->nodes[i]->node_index;

      for (j = 0; j < k; ++j)
      {

        unsigned int xi = nodes[j]->node_index;

        if (!opt_mig_bitmatrix[xi][yi] && !opt_mig_bitmatrix[yi][xi]) continue;

        /* check whether a flip (migration appearing or disappearing) occurs.

           The following figures show all possible tau configurations between two
           edges p-a and q-b and whether migration is possible for each case

           (a)       (b)         (c)           (d)           (e)           (f)
                                                  
           q o    p o              q o    p o      
             |      |                |      |       
           b o    a o       p o      |      |    q o    p o                  q o
                              | <--> |      | <--> |      |    q o    p o      |
        p o          q o      |    b o    a o      |      | <--> |      | <--> |
          |            |      |                    |      |    b o    a o      |
        a o          b o    a o                  b o    a o                  b o

           
          No migration in a and b (model H0) and migration in c,d,e,f (model H1)

          Now determine the model we were before the tau change (H0 or H1), and
          whether we flipped to the other model after the change.
        */

        double a = old_atau[j];
        double p = old_ptau[j];
        double b = stree->nodes[i]->tau;
        double q = stree->nodes[i]->parent->tau;

        /* get the model *before* tau change (0: no migration, 1: migration) */
        amodel = !(p < b || a > q);  /* check cases (a) and (b) */

        a = atau[j];
        p = ptau[j];

        /* now get the model *after* the change to tau */
        bmodel = !(p < b || a > q);

        /* flip = -1 if no flip, 1: migration appearing or 0: disappearing */
        long flip = (amodel == bmodel) ? -1 : bmodel;

        if (flip == 0)
        {
          /* M disappearing */

          if (opt_mig_bitmatrix[xi][yi])
          {
            spec = opt_mig_specs+opt_migration_matrix[xi][yi];
            if (spec->pseudo_a)
            {
              assert(spec->am == 0);  /* TODO: 5 params in mig specification */
              lnacceptance += logPDFGamma(spec->M, spec->pseudo_a, spec->pseudo_b);
              lnacceptance -= logPDFGamma(spec->M, spec->alpha, spec->beta);
            }
          }

          if (opt_mig_bitmatrix[yi][xi])
          {
            spec = opt_mig_specs+opt_migration_matrix[yi][xi];
            if (spec->pseudo_a)
            {
              assert(spec->am == 0);  /* TODO: 5 params in mig specification */
              lnacceptance += logPDFGamma(spec->M, spec->pseudo_a, spec->pseudo_b);
              lnacceptance -= logPDFGamma(spec->M, spec->alpha, spec->beta);
            }
          }
        }
        else if (flip == 1)
        {
          /* M appearing */

          if (opt_mig_bitmatrix[xi][yi])
          {
            spec = opt_mig_specs+opt_migration_matrix[xi][yi];
            if (spec->pseudo_a)
            {
              assert(spec->am == 0);  /* TODO: 5 params in mig specification */
              lnacceptance -= logPDFGamma(spec->M, spec->pseudo_a, opt_pseudo_beta);
              lnacceptance += logPDFGamma(spec->M, spec->alpha, spec->beta);
            }
          }

          if (opt_mig_bitmatrix[yi][xi])
          {
            spec = opt_mig_specs+opt_migration_matrix[yi][xi];
            if (spec->pseudo_a)
            {
              assert(spec->am == 0);  /* TODO: 5 params in mig specification */
              lnacceptance -= logPDFGamma(spec->M, spec->pseudo_a, spec->pseudo_b);
              lnacceptance += logPDFGamma(spec->M, spec->alpha, spec->beta);
            }
          }
        }
      }
    }
  }
  #endif

  lnacceptance += logpr_diff + logl_diff + count_below*log(minfactor) +
                  count_above*log(maxfactor);

  if (opt_debug_tau)
    printf("[Debug] (tau) lnacceptance = %f\n", lnacceptance);

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */
    accepted++;

    for (i = 0; i < stree->locus_count; ++i)
    {
      /* The commented loop is probably sufficient */
      //gnode_t ** gt_nodesptr = __gt_nodes + __gt_nodes_index[i];
      //for (j = 0; j < k + __extra_count[i]; ++j)
      //  gt_nodesptr[j]->mark = 0;
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        gtree[i]->nodes[j]->mark = 0;
    }
  }

  return accepted;
}

double stree_propose_tau(gtree_t ** gtree, stree_t * stree, locus_t ** loci)
{
  unsigned int i;
  unsigned int candidate_count = 0;
  long accepted = 0;

  long thread_index = 0;

  /* compute number of nodes with tau > 0 */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->tau > 0 && (!opt_msci || stree->nodes[i]->prop_tau))
      candidate_count++;


  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau > 0 && (!opt_msci || stree->nodes[i]->prop_tau))
      accepted += propose_tau(loci,
                              stree->nodes[i],
                              gtree,
                              stree,
                              candidate_count,
                              thread_index);
  }

  return ((double)accepted / candidate_count);
}

double stree_propose_tau_mig(stree_t ** streeptr,
                             gtree_t *** gtreeptr,
                             stree_t ** scloneptr,
                             gtree_t *** gcloneptr,
                             locus_t ** loci)
{
  unsigned int i,j, total_nodes;
  unsigned int candidate_count = 0;
  long accepted = 0;
  long rc;

  long thread_index = 0;

  int * candidate = NULL;
  assert(opt_migration && !opt_msci && opt_est_theta);

  stree_t * original_stree = *streeptr;
  gtree_t ** original_gtree = *gtreeptr;

  stree_t * stree = *scloneptr;
  gtree_t ** gtree = *gcloneptr;

  total_nodes = stree->tip_count + stree->inner_count;
  candidate = (int *)xcalloc((size_t)(total_nodes), sizeof(int));

  /* compute number of nodes with tau > 0 */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau > 0)
    {
      candidate_count++;
      candidate[i] = 1;
    }
  }

  for (i = 0; i < total_nodes; ++i)
  {
    if (candidate[i])
    {
      /* clone species tree and gene trees */
      stree_clone(original_stree,stree);
      for (j = 0; j < opt_locus_count; ++j)
        gtree_clone(original_gtree[j], gtree[j], stree);
      events_clone(original_stree, stree, gtree);
      
      rc = propose_tau_mig(loci,
                           stree->nodes[i],
                           gtree,
                           stree,
                           candidate_count,
                           thread_index);
      if (rc == -1) { /* gphocs rejection */ rc = 0; }
      if (rc)
      {
        accepted++;

        /* swap trees */
        SWAP(original_stree,stree);
        SWAP(original_gtree,gtree);
      }
    }
  }
  *streeptr = original_stree;
  *gtreeptr = original_gtree;

  *scloneptr = stree;
  *gcloneptr = gtree;

  free(candidate);

  return ((double)accepted / candidate_count);
}

void stree_rootdist(stree_t * stree,
                    list_t * maplist,
                    msa_t ** msalist,
                    unsigned int ** weights)
{
  unsigned int i, j, k, n;
  unsigned int msa_count = stree->locus_count;
  unsigned int lroot_index = stree->root->left->node_index;
  unsigned int jpop_index, kpop_index;
  unsigned int diff_count;
  unsigned int locus_used = 0;
  double diff_locus;
  double diff_pair;
  double md, vd;
  char * jseq;
  char * kseq;

  assert(msa_count);

  stree->root_age = (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA) ?
                      opt_tau_beta / (opt_tau_alpha - 1) * 4 :
                      QuantileGamma(0.9, opt_tau_alpha, opt_tau_beta);

  if (!opt_usedata) return;

  if (opt_diploid)
  {
    for (i = 0; i < stree->tip_count; ++i)
      if (stree->nodes[i]->diploid)
        return;
  }

  hashtable_t * sht = species_hash(stree);
  hashtable_t * mht = maplist_hash(maplist, sht);

  /* find max alignment length */
  size_t maxalloc = 0;
  for (i = 0; i < msa_count; ++i)
    if ((size_t)(msalist[i]->count) > maxalloc)
      maxalloc = msalist[i]->count;

  snode_t ** pop = (snode_t **)xmalloc(maxalloc * sizeof(snode_t *));

  md = vd = 0;

  for (i = 0; i < msa_count; ++i)
  {
    diff_count = 0;
    diff_locus = 0;

    /* find the population for each sequence */
    for (j = 0; j < (unsigned int)(msalist[i]->count); ++j)
    {
      char * label = msalist[i]->label[j];
      label = strchr(label, '^');
      if (!label)
        fatal("Cannot find species tag on sequence %s of locus %d",
              msalist[i]->label[j], i);

      /* skip the '^' mark */
      label++;
      if (!(*label))
        fatal("Sequence %s of locus %d contains no label",
              msalist[i]->label[j], i);

      pair_t * pair = hashtable_find(mht,
                                     (void *)label,
                                     hash_fnv(label),
                                     cb_cmp_pairlabel);

      if (!pair)
        fatal("!! Cannot find species mapping for sequence %s of locus %d",
              label, i);

      pop[j] = (snode_t *)(pair->data);
    }


    for (j = 0; j < (unsigned int)(msalist[i]->count); ++j)
    {
      for (k = j + 1; k < (unsigned int)(msalist[i]->count); ++k)
      {
        /* get population index for node j and k */
        jpop_index = pop[j]->node_index;
        kpop_index = pop[k]->node_index;

        /* check that the mRCA of the two nodes is the species tree root */
        if (stree->pptable[jpop_index][lroot_index] !=
           stree->pptable[kpop_index][lroot_index])
        {
          /* obtain the two sequences */
          jseq = msalist[i]->sequence[j];
          kseq = msalist[i]->sequence[k];
          diff_pair = 0;

          for (n = 0; n < (unsigned int)(msalist[i]->length); ++n)
          {
            if (jseq[n] != kseq[n])
              diff_pair += weights[i][n];
          }
          diff_pair /= msalist[i]->original_length;
          diff_locus += diff_pair;
          ++diff_count;
        }
      }
    }
    if (!diff_count) continue;

    locus_used++;

    diff_locus /= (2 * diff_count);
    vd += (diff_locus - md)*(diff_locus - md) * (locus_used - 1) / locus_used;
    md = (md * (locus_used - 1) + diff_locus) / locus_used;
  }

  vd /= msa_count;
  if (locus_used >= 2)
  {
    double theta = 2 * sqrt(vd);
    theta = sqrt(vd * 4 + 1) - 1;
    theta = (2 * sqrt(vd) + sqrt(vd * 4 + 1) - 1) / 2;
    if (md - theta / 2 > 0)
      stree->root_age = md - theta / 2;
    else
      stree->root_age = md;
  }
  else
    stree->root_age = md;

  //printf("root dist = %7.5f\n", stree->root_age);
  hashtable_destroy(sht, NULL);
  hashtable_destroy(mht, cb_dealloc_pairlabel);
  free(pop);
}

static void init_weights(stree_t * stree, int * feasible)
{
  unsigned int i;
  double sum = 0;

  for (i = stree->tip_count; i < stree->inner_count + stree->tip_count; ++i)
  {
    if (feasible && !feasible[i-stree->tip_count])
    {
      stree->nodes[i]->weight = 0;
      continue;
    }

    if (stree->nodes[i]->parent && stree->nodes[i]->tau > 0)
    {
      stree->nodes[i]->weight = 1.0 / sqrt(stree->nodes[i]->parent->tau -
                                stree->nodes[i]->tau);
      sum += stree->nodes[i]->weight;
    }
    else
      stree->nodes[i]->weight = 0;
  }
  for (i = stree->tip_count; i < stree->inner_count + stree->tip_count; ++i)
    if (stree->nodes[i]->weight)
      stree->nodes[i]->weight /= sum;
}

static int check_age_feasibility_recursive(snode_t * c,
                                           double tau,
                                           long constraint)
{
  int rc1 = 0;
  int rc2 = 0;

  if (c->constraint != constraint) return 0;

  if (c->parent->tau <= tau) return 0;

  if (c->tau < tau && c->parent->tau > tau) return 1;

  if (c->left)
    rc1 = check_age_feasibility_recursive(c->left,tau,constraint);
  if (c->right)
    rc2 = check_age_feasibility_recursive(c->right,tau,constraint);

  return (rc1 || rc2);
}

static int * fill_feasible_flags(stree_t * stree)
{
  unsigned int i;
  int * feasible = NULL;

  if (!opt_constraint_count) return NULL;
  
  feasible = (int *)xcalloc((size_t)(stree->inner_count), sizeof(int));

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * y = stree->nodes[i];

    if (!y->parent) continue;

    snode_t * sister = (y->parent->left == y) ?
                         y->parent->right : y->parent->left;


    if (!(y->constraint == y->left->constraint &&
          y->constraint == y->right->constraint &&
          y->constraint == sister->constraint))
      continue;

    /* now check if there exists a target branch that whose time interval
       overlaps with y->tau. Traverse the nodes on the path to the root and
       visit all compatible branching clades on that path */
    int flag = 0;
    snode_t * p = y;
    while (p->parent && !flag)
    {
      if (p->constraint != y->constraint) break;
      sister = (p->parent->left == p) ? p->parent->right : p->parent->left;
      assert(sister->constraint == y->constraint);

      flag = check_age_feasibility_recursive(sister,y->tau,y->constraint);
      p = p->parent; 
    }

    feasible[i-stree->tip_count] = flag;
  }

  return feasible;
}

/* Algorithm implemented according to Figure 1 in:
   Rannala, B., Yang, Z. Efficient Bayesian species tree inference under the 
   multispecies coalescent.  Systematic Biology, 2017, 66:823-842.
*/
long stree_propose_spr(stree_t ** streeptr,
                       gtree_t *** gtree_list_ptr,
                       stree_t ** scloneptr,
                       gtree_t *** gclonesptr,
                       locus_t ** loci)
{
  unsigned int i, j, k = 0;
  unsigned int branch_update_count;
  long target_count = 0;
  long source_count = 0;
  double r;
  double sum = 0;
  double lnacceptance = 0;

  long thread_index = 0;

  /* TODO: If relaxed clock, we need to detect gene tree branches that simply
     pass by (without coalescing) the pruned subtree
  */

  /* the following clones the species tree and gene trees, and then
     we work on a copy */
  stree_t * original_stree = *streeptr;
  gtree_t ** original_gtree_list = *gtree_list_ptr;

  stree_t * stree = *scloneptr;
  gtree_t ** gtree_list = *gclonesptr;

  stree_clone(original_stree, stree);

  for (i = 0; i < stree->locus_count; ++i)
    gtree_clone(original_gtree_list[i], gtree_list[i], stree);
  events_clone(original_stree, stree, gtree_list);

  double oldprior = lnprior_species_model(stree);

  int * feasible = fill_feasible_flags(stree);

  /* calculate the weight of each branch as the reciprocal of the square root
   * of its length */
  init_weights(stree,feasible);
  if (feasible)
    free(feasible);  /* no longer required */

  /* randomly select a branch according to weights */
  r = legacy_rndu(thread_index);
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count - 1; ++i)
  {
    sum += stree->nodes[i]->weight;
    if (r < sum) break;
  }

  /* selected node */
  snode_t * y = stree->nodes[i];

  assert(y != stree->root);

  if (opt_constraint_count && !stree->nodes[i]->weight) return 2;
  assert(stree->nodes[i]->weight);
  lnacceptance -= log(stree->nodes[i]->weight);

  /* parent of node */
  snode_t * x = y->parent;

#if 0
   /* sibling of y - not used */
   snode_t * c0 = (x->left == y) ? x->right : x->left;
#endif

  /* Randomly select children of y in randomly selected order */
  snode_t *a, *b;
  if ((int)(2 * legacy_rndu(thread_index)) == 0)
  {
    a = y->left;
    b = y->right;
  }
  else
  {
    a = y->right;
    b = y->left;
  }

  /* find all nodes that are candidates for becoming node C (Figure 1) */
  target_count = 0;
  for (i = 0, sum = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * c_cand;  /* candidate for node C */
    snode_t * z_cand;  /* candidate for node Z */
    snode_t * tmp;

    c_cand = stree->nodes[i];

    /* A C candidate node must fulfill the following FOUR properties:
       i) it is not a descendant of y,
       ii) is younger than y,
       iii) its parent is older than y
       iv) matching constraints between c and y */
    if (stree->pptable[i][y->node_index] ||
        c_cand->tau >= y->tau ||
        c_cand->parent->tau <= y->tau ||
        c_cand->constraint != y->constraint) continue;

    /* compute z_cand as the lowest common ancestor of c_cand and y */
    for (z_cand = c_cand->parent; z_cand; z_cand = z_cand->parent)
      if (stree->pptable[x->node_index][z_cand->node_index])
        break;

    /* compute the weight as the reciprocal of number of nodes on the shortest
       path between c_cand and y */
    target_weight[target_count] = 1; /* TODO: should this be 2? */
    for (tmp = y; tmp != z_cand; tmp = tmp->parent)
      target_weight[target_count]++;
    for (tmp = c_cand; tmp != z_cand; tmp = tmp->parent)
      target_weight[target_count]++;
    target_weight[target_count] = 1 / target_weight[target_count];
    sum += target_weight[target_count];


    target[target_count++] = c_cand;
  }

  /* normalize to weights to probabilities */
  for (i = 0; i < target_count; ++i)
    target_weight[i] /= sum;

  /* randomly select one node among the candidates to become node C */
  r = legacy_rndu(thread_index);
  for (i = 0, sum = 0; i < target_count - 1; ++i)
  {
    sum += target_weight[i];
    if (r < sum) break;
  }
  snode_t * c = target[i];

  /* constraint check for the quick-and-dirty method of applying constraints */
  if (a->constraint == c->constraint)
  {
    if (y->constraint != a->constraint)
    {
      b->constraint = y->constraint;
      y->constraint = a->constraint;
    }
  }
  else
  {
    if (c->left &&
        c->left->constraint == c->right->constraint &&
        c->left->constraint == a->constraint &&
        a->constraint == y->constraint)
    {
      /* if the move is accepted we need ot change the constraint id  of the new
       * Y* and of C. Since we work on a copy, we set it now */
      long tmp_const = c->constraint;
      c->constraint = y->constraint; y->constraint = tmp_const;
      /* SWAP(y->constraint,c->constraint); */
    }
    else
      assert(0);
  }

  lnacceptance -= log(target_weight[i]);

  /* now compute node Z, i.e. the LCA of C and Y */
  snode_t * z;
  for (z = c->parent; z; z = z->parent)
    if (stree->pptable[x->node_index][z->node_index])
      break;
  assert(z);

  /* now create two arrays that hold the nodes from A to Z and C to Z always
     excluding Z */
     /*
     snode_t * patha;
     snode_t * pathc;
     long patha_size = 0,pathc_size = 0;
     for (temp = y; temp != z; temp = temp->father)
       patha[.....

     we do not need them actually
     */

     /* perform SPR to modify gene tree topologies */
  gnode_t ** moved_nodes = moved_space;
  gnode_t ** gtarget_list = gtarget_temp_space;
  gnode_t ** gtarget_nodes = gtarget_space;
  gnode_t ** bl_list = __gt_nodes;
  snode_t ** snode_contrib = snode_contrib_space;
  for (i = 0; i < stree->locus_count; ++i)
  {
    snode_contrib_count[i] = 0;

    branch_update_count = 0;
    gtree_t * gtree = gtree_list[i];

    /* mark all nodes in gene tree paths starting from some tip and end at Z
       (excluding nodes in Z), but which also include at least one node in A */
    for (j = 0; j < gtree->tip_count; ++j)
    {
      gnode_t * tmp;
      unsigned int pop_index = gtree->nodes[j]->pop->node_index;
      if (!stree->pptable[pop_index][a->node_index]) continue;
      gtree->nodes[j]->mark = LINEAGE_A;
      for (tmp = gtree->nodes[j]->parent; tmp->mark == 0; tmp = tmp->parent)
      {
        if (tmp->pop == z || stree->pptable[z->node_index][tmp->pop->node_index])
          break;
        tmp->mark = LINEAGE_A;
        if (!tmp->parent) break;
      }
    }

    /* now mark all gene tree nodes that have ancestral populations in the path
       between A and Z (excluding both A and Z) */
    for (j = 0; j < gtree->tip_count; ++j)
    {
      gnode_t * node = gtree->nodes[j];
      snode_t * stmp;
      gnode_t * gtmp;
      unsigned int pop_index = node->pop->node_index;

      /* if A is an ancestor skip */
      if (stree->pptable[pop_index][a->node_index]) continue;

      /* if node has no ancestor between Y and Z (excluding) skip */
      for (stmp = y; stmp != z; stmp = stmp->parent)
        if (stree->pptable[pop_index][stmp->node_index])
          break;
      if (stmp == z) continue;

      node->mark |= LINEAGE_OTHER;

      for (gtmp = node->parent; !(gtmp->mark & LINEAGE_OTHER); gtmp = gtmp->parent)
      {
        if (gtmp->pop == z || stree->pptable[z->node_index][gtmp->pop->node_index])
          break;
        gtmp->mark |= LINEAGE_OTHER;
        if (!gtmp->parent) break;
      }
    }


    /* now identify Moved nodes */
    moved_count[i] = 0;
    for (j = gtree->tip_count; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      snode_t * pop_az;
      gnode_t * node = gtree->nodes[j];

      /* if node has no ancestor between Y and Z (excluding) skip */
      for (pop_az = y; pop_az != z; pop_az = pop_az->parent)
        //if (stree->pptable[pop_index][pop_az->node_index])
        if (node->pop == pop_az)
          break;
      if (pop_az == z) continue;

      /* Square nodes */
      if (node->pop == y &&
          (node->left->mark & LINEAGE_OTHER) &&
          (node->right->mark & LINEAGE_OTHER))
      {
        node->mark |= NODE_SQUARE;
        continue;
      }

      /* now we need to ensure that only one child has descendants in A only */
      int count = 0;
      gnode_t * pruned = NULL;
      gnode_t * intact = NULL;
      if (node->left->mark == LINEAGE_A)
      {
        count++;
        pruned = node->left;
        intact = node->right;
      }
      if (node->right->mark == LINEAGE_A)
      {
        count++;
        pruned = node->right;
        intact = node->left;
      }

      if (count != 1) continue;

      node->mark |= NODE_MOVED;
      moved_nodes[moved_count[i]] = node;
      pruned_nodes[moved_count[i]++] = pruned;

      node->mark |= FLAG_PARTIAL_UPDATE;
      if (node->parent)
        node->parent->mark |= FLAG_PARTIAL_UPDATE; /* required as the parent will change */

      snode_t * pop_cz = c;
      while (pop_cz->parent != z)
      {
        if (pop_cz->parent->tau >= node->time)
          break;
        pop_cz = pop_cz->parent;
      }

      /* now make a list of feasible target node */
      target_count = 0;
      for (k = 0; k < gtree->tip_count + gtree->inner_count; ++k)
      {
        gnode_t * tmp = gtree->nodes[k];

        if (tmp->time >= node->time || tmp->parent->time <= node->time)
          continue;

        if (stree->pptable[tmp->pop->node_index][pop_cz->node_index])
          gtarget_list[target_count++] = tmp;
      }

      if (!target_count)
        return 2;

      /* revolutionary methods */
      double twgt = 1;
      if (opt_revolutionary_spr_method)
      {
        unsigned int n;

        /* if more than one target nodes, select according to likelihood */
        if (target_count > 1)
        {
          /* allocate space for target nodes weights */
          double * tweight = (double *)xmalloc(target_count * sizeof(double));

          /* compute weights (logl) for each target node in gtarget_list and store them in tweight */
          revolutionary_spr_tselect_logl(pruned, gtarget_list, target_count, loci[i], tweight);

          /* normalize target weights */
          double maxw = tweight[0];
          for (n = 1; n < target_count; ++n)
            if (maxw < tweight[n])  maxw = tweight[n];
          for (sum = 0, n = 0; n < target_count; ++n)
          {
            tweight[n] = exp(tweight[n] - maxw);
            sum += tweight[n];
          }
 
          /* randomly select a target node according to weights */
          r = legacy_rndu(thread_index) * sum;
          double rsum = 0;
          for (n = 0; n < target_count - 1; ++n)  /* Z: no need for last comparison and rndu may be 1. */
          {
            rsum += tweight[n];
            if (r < rsum) break;
          }
          gtarget_nodes[moved_count[i] - 1] = gtarget_list[n];

          twgt = tweight[n]/sum;

          /* free weights */
          free(tweight);
        }
        else
          gtarget_nodes[moved_count[i] - 1] = gtarget_list[0];
      }
      else  /* randomly select a target from list */
        gtarget_nodes[moved_count[i] - 1] = gtarget_list[(int)(target_count*legacy_rndu(thread_index))];

      source_count = 1;
      gsources_list[0] = intact;
      for (k = 0; k < gtree->tip_count + gtree->inner_count; ++k)
      {
        gnode_t * tmp = gtree->nodes[k];

        if (tmp == intact || tmp == intact->parent)
          continue;  /* intact->father equals node */

        if (tmp->time >= node->time || tmp->parent->time <= node->time)
          continue;

        /* TODO: gsources_list is not required!!! */
        if (stree->pptable[tmp->pop->node_index][pop_az->node_index] && tmp->mark != LINEAGE_A)
          gsources_list[source_count++] = tmp;
      }

      if (opt_revolutionary_spr_method)
      {
        unsigned int n;
        double swgt = 1;

        /* if if more than one sources, select according to likelihood */
        if (source_count > 1)
        {
          double * tweight = (double *)xmalloc(source_count * sizeof(double));

          /* if a moved node appears in the sources, then replace it by
             descending until we find a non-moved node (not in LINEAGE_A) */
          for (n = 0; n < source_count; ++n)
          {
            while (gsources_list[n]->mark & NODE_MOVED)
            {
              gsources_list[n] = gsources_list[n]->left->mark & LINEAGE_A ?
                 gsources_list[n]->right : gsources_list[n]->left;
            }
          }

          /* compute weights and store in tweight */
          revolutionary_spr_tselect_logl(pruned, gsources_list, source_count, loci[i], tweight);
          /* normalize target node weights */
          double maxw = tweight[0];
          for (n = 1; n < source_count; ++n)
            if (maxw < tweight[n])  maxw = tweight[n];
          for (sum = 0, n = 0; n < source_count; ++n)
          {
            tweight[n] = exp(tweight[n] - maxw);
            sum += tweight[n];
          }
          /* we do not select randomly, as we already know the branch */
          gnode_t * srcnode = intact;
          while (srcnode->mark & NODE_MOVED)
            srcnode = srcnode->left->mark & LINEAGE_A ? srcnode->right : srcnode->left;

          for (n = 0; n < source_count - 1; ++n)
            if (gsources_list[n] == srcnode)
              break;

          swgt = tweight[n]/sum;

          /* free weights */
          free(tweight);
        }

        if (opt_revolutionary_spr_debug)
          printf("ntarget = %2ld weight = %9.6f  nsource = %2ld weight = %9.6f\n", target_count, twgt, source_count, swgt);
        
        assert(swgt > 0 && twgt > 0);

        lnacceptance += log(swgt / twgt);
      }
      else
        lnacceptance += log((double)target_count / source_count);
    }

    /* All moves nodes for current locus are now identified. Apply SPR to gene
       tree. */
    for (j = 0; j < moved_count[i]; ++j)
    {
      snode_t * pop_cz = c;
      while (pop_cz->parent != z)
      {
        if (pop_cz->parent->tau >= moved_nodes[j]->time)
          break;
        pop_cz = pop_cz->parent;
      }

      /* TODO: We probably don't need to keep the pruned nodes array above, but only check the 'mark' */
      gnode_t * node = pruned_nodes[j]->parent;
      assert(node == moved_nodes[j]);
      gnode_t * pruned = pruned_nodes[j];
      gnode_t * intact = (node->left == pruned) ? node->right : node->left;

#if 0
         if (!(intact->pop->mark & FLAG_POP_UPDATE))
         {
            intact->pop->mark |= FLAG_POP_UPDATE;
            snode_contrib[snode_contrib_count[i]++] = intact->pop;
         }

         if (!(node->pop->mark & FLAG_POP_UPDATE))
         {
            node->pop->mark |= FLAG_POP_UPDATE;
            snode_contrib[snode_contrib_count[i]++] = node->pop;
         }

         if (!(node->parent->pop->mark & FLAG_POP_UPDATE))
         {
            node->parent->pop->mark |= FLAG_POP_UPDATE;
            snode_contrib[snode_contrib_count[i]++] = node->parent->pop;
         }


         if (!(pruned->pop->mark & FLAG_POP_UPDATE))
         {
            pruned->pop->mark |= FLAG_POP_UPDATE;
            snode_contrib[snode_contrib_count[i]++] = pruned->pop;
         }
#endif

#if 0
         /* TODO: The above line should point to the same node as stored in moved_nodes. Let's check it */
         assert(node == moved_nodes[j]);
#endif

        /* link parent with child */
      intact->parent = node->parent;
      if (node->parent->left == node)
        node->parent->left = intact;
      else
        node->parent->right = intact;

      gnode_t * receiver;
      /* TODO: is the below loop necesary? target_nodes[i]->parent should ALWAYS have time larger than moved_nodes[i] */
      for (receiver = gtarget_nodes[j]; ; receiver = receiver->parent)
        if (receiver->parent->time > moved_nodes[j]->time)  /* TODO: moved_nodes[j] is node */
          break;

      /* regraft */
      if (receiver->parent->left == receiver)
        receiver->parent->left = node;
      else
        receiver->parent->right = node;

      node->parent = receiver->parent;
      if (node->left == pruned)
        node->right = receiver;
      else
        node->left = receiver;

#if 0
         /* TODO I THINK THE BELOW 1 CAN BE SET TO 0 */
#if 0
         if (!(node->parent->pop->mark & FLAG_POP_UPDATE))
         {
            node->parent->pop->mark |= FLAG_POP_UPDATE;
            snode_contrib[snode_contrib_count[i]++] = node->parent->pop;
         }
#endif

#if 0
         if (!(receiver->pop->mark & FLAG_POP_UPDATE))
         {
            receiver->pop->mark |= FLAG_POP_UPDATE;
            snode_contrib[snode_contrib_count[i]++] = receiver->pop;
         }
#endif
#endif

      receiver->parent = node;

      /* TODO : CHECK IF the following three IFs can fail, i.e. are they necessary? */

      /* three nodes need to have their branches updated */
      if (!(receiver->mark & FLAG_BRANCH_UPDATE))
      {
        bl_list[branch_update_count++] = receiver;
        receiver->mark |= FLAG_BRANCH_UPDATE;
      }

      if (!(node->mark & FLAG_BRANCH_UPDATE))
      {
        bl_list[branch_update_count++] = node;
        node->mark |= FLAG_BRANCH_UPDATE;
      }

      if (!(intact->mark & FLAG_BRANCH_UPDATE))
      {
        bl_list[branch_update_count++] = intact;
        intact->mark |= FLAG_BRANCH_UPDATE;
      }

      /* remove  gene node from list of coalescent events of its old population */
      unlink_event(node, i);

      node->pop->event_count[i]--;
      if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
      {
        node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
        snode_contrib[snode_contrib_count[i]++] = node->pop;
      }
      if (!opt_est_theta)
        node->pop->event_count_sum--;

      node->pop = pop_cz;
      if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
      {
        node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
        snode_contrib[snode_contrib_count[i]++] = node->pop;
      }

      dlist_item_append(node->pop->event[i], node->event);

      node->pop->event_count[i]++;
      if (!opt_est_theta)
        node->pop->event_count_sum++;

      /* update leaf counts */
      while (node)
      {
        node->leaves = node->left->leaves + node->right->leaves;
        node = node->parent;
      }

      for (node = intact->parent; node; node = node->parent)
        node->leaves = node->left->leaves + node->right->leaves;
    }

    /* Now process square nodes */
    for (j = gtree->tip_count; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      gnode_t * node = gtree->nodes[j];

      if (node->mark & NODE_SQUARE)
      {
        /* remove  gene node from list of coalescent events of its old population */
        unlink_event(node, i);

        node->pop->event_count[i]--;
        if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
        {
          node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }
        if (!opt_est_theta)
          node->pop->event_count_sum--;

        node->pop = b;
        if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
        {
          node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        dlist_item_append(node->pop->event[i], node->event);

        node->pop->event_count[i]++;
        if (!opt_est_theta)
          node->pop->event_count_sum++;
      }
      else if (node->pop == c && node->time > y->tau)
      {
        /* diamond nodes */

        if (opt_clock != BPP_CLOCK_GLOBAL)
        {
          /* flag diamond node branches to be updated */
          if (!(node->mark & FLAG_BRANCH_UPDATE) && node->parent)
          {
            bl_list[branch_update_count++] = node;
            node->mark |= FLAG_BRANCH_UPDATE;
            assert(node->parent);
            node->parent->mark |= FLAG_PARTIAL_UPDATE;

          }
        }
        /* remove  gene node from list of coalescent events of its old population */
        unlink_event(node, i);

        node->pop->event_count[i]--;
        if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
        {
          node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }
        if (!opt_est_theta)
          node->pop->event_count_sum--;

        node->pop = y;
        if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
        {
          node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        dlist_item_append(node->pop->event[i], node->event);

        node->pop->event_count[i]++;
        if (!opt_est_theta)
          node->pop->event_count_sum++;
      }
      //else if (node->mark == LINEAGE_A && node->time > y->tau && node->time < z->tau)
      else if ((node->mark & LINEAGE_A && !(node->mark & LINEAGE_OTHER)) && node->time > y->tau && node->time < z->tau)
      {
        /* circle and triangle nodes */

        /* find oldest population pop in C-Z (excluding Z) such that node->time > pop->tau */
        snode_t * pop;
        for (pop = c; pop->parent != z; pop = pop->parent)
          if (pop->parent->tau >= node->time)
            break;

        if (opt_clock != BPP_CLOCK_GLOBAL)
        {
          /* flag diamond node branches to be updated */
          if (!(node->mark & FLAG_BRANCH_UPDATE) && node->parent &&
              ((node->parent->time >= c->parent->tau) || (node->parent->time >= y->parent->tau)))
          {
            bl_list[branch_update_count++] = node;
            node->mark |= FLAG_BRANCH_UPDATE;
            assert(node->parent);
            node->parent->mark |= FLAG_PARTIAL_UPDATE;
          }
        }

        unlink_event(node, i);

        node->pop->event_count[i]--;
        if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
        {
          node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }
        if (!opt_est_theta)
          node->pop->event_count_sum--;

        if (pop == c)
          node->pop = y;
        else
          node->pop = pop;

        if (!(node->pop->mark[thread_index] & FLAG_POP_UPDATE))
        {
          node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        dlist_item_append(node->pop->event[i], node->event);

        node->pop->event_count[i]++;
        if (!opt_est_theta)
          node->pop->event_count_sum++;
      }
    }

    /* Flag populations Y,C and B for re-computing their contribution to the i-th gene tree probability only if:
         (a) they have not been already flagged in a previous step
         (b) there is more than one outgoing lineages (entering its parent population).
    */
    if (!(y->mark[thread_index] & FLAG_POP_UPDATE) && (y->seqin_count[i] - y->event_count[i] > 1))
      snode_contrib[snode_contrib_count[i]++] = y;
    if (!(c->mark[thread_index] & FLAG_POP_UPDATE) && (c->seqin_count[i] - c->event_count[i] > 1))
      snode_contrib[snode_contrib_count[i]++] = c;
    if (!(b->mark[thread_index] & FLAG_POP_UPDATE) && (b->seqin_count[i] - b->event_count[i] > 1))
      snode_contrib[snode_contrib_count[i]++] = b;

    moved_nodes += gtree->inner_count;
    gtarget_nodes += gtree->inner_count;
    gtarget_list += gtree->tip_count + gtree->inner_count;

    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      /* relaxed clock */
      for (j = 0; j < gtree_list[i]->tip_count+gtree_list[i]->inner_count; ++j)
      {
        gnode_t * node = gtree_list[i]->nodes[j];

        /* if root or parent already marked for branch update, skip */
        if (!node->parent || (node->mark & FLAG_BRANCH_UPDATE)) continue;

        if (stree->pptable[node->pop->node_index][b->node_index])
        {
          if (node->parent->time >= y->tau)
          {
            bl_list[branch_update_count++] = node;
            node->mark |= FLAG_BRANCH_UPDATE;
            node->parent->mark |= FLAG_PARTIAL_UPDATE;
          }
        }
        else if (stree->pptable[node->pop->node_index][c->node_index] &&
                 node->parent->time >= y->tau)
        {
          /* diamond node branches have already been flagged. Here we flag edges
             of nodes in C or descendants that coalesce above C */
          bl_list[branch_update_count++] = node;
          node->mark |= FLAG_BRANCH_UPDATE;
          node->parent->mark |= FLAG_PARTIAL_UPDATE;
        }
        else if (stree->pptable[node->pop->node_index][a->node_index] &&
                 ((node->parent->time >= c->parent->tau) || node->parent->time >= y->parent->tau))
        {
          bl_list[branch_update_count++] = node;
          node->mark |= FLAG_BRANCH_UPDATE;
          node->parent->mark |= FLAG_PARTIAL_UPDATE;
        }
        else if ((node->mark & NODE_MOVED) && 
                 ((node->time >= c->parent->tau) || (node->time >= y->parent->tau)))
        {
          bl_list[branch_update_count++] = node;
          node->mark |= FLAG_BRANCH_UPDATE;
          node->parent->mark |= FLAG_PARTIAL_UPDATE;
        }
        else if (node->pop == y &&
                 (!(node->mark & NODE_MOVED) && (node->mark & LINEAGE_OTHER)))
        {
            bl_list[branch_update_count++] = node;
            node->mark |= FLAG_BRANCH_UPDATE;
            node->parent->mark |= FLAG_PARTIAL_UPDATE;
        }
      }
    }

    __mark_count[i] = branch_update_count;
    bl_list += branch_update_count;
    snode_contrib += stree->tip_count + stree->inner_count;

    /* reset species tree marks */
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      stree->nodes[j]->mark[thread_index] = 0;
  } /* end of locus */

   /* update species tree */
#if 0
   printf("Y = %s\n", y->label);
   printf("A = %s\n", a->label);
   printf("B = %s\n", b->label);
   printf("C = %s\n", c->label);
   printf("Z = %s\n", z->label);
#endif

  /* make b child of y->parent */
  if (y->parent->left == y)
    y->parent->left = b;
  else
    y->parent->right = b;

  /* make y->parent parent of b */
  b->parent = y->parent;

  /* make y child of c->parent */
  if (c->parent->left == c)
    c->parent->left = y;
  else
    c->parent->right = y;

  /* make old c->parent parent of y */
  y->parent = c->parent;

  /* make y parent of c */
  c->parent = y;

  /* make c child of y */
  if (y->left == a)
    y->right = c;
  else
    y->left = c;

  assert(!opt_msci);
  for (i = 0; i < stree->locus_count; ++i)
    fill_seqin_counts(stree, NULL, i);

  /* TODO: Check whether reset_gene_leaves_count must operate on the whole gtree_list, or whether
     we can separate the 'reset_hybrid_gene_leaves_count() call inside the function */
  assert(!opt_msci);
  reset_gene_leaves_count(stree,gtree_list);
  stree_reset_pptable(stree);

  feasible = fill_feasible_flags(stree);
  init_weights(stree,feasible);
  if (feasible)
    free(feasible);

  /* probability of choosing focus branch in reverse move */
  lnacceptance += log(y->weight);

  /* probability of sampling target branches in reverse move */

  target_count = 0;
  for (i = 0, sum = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * c_cand;
    snode_t * z_cand;
    snode_t * tmp;

    c_cand = stree->nodes[i];

    if (stree->pptable[i][y->node_index] ||
        c_cand->tau >= y->tau ||
        c_cand->parent->tau <= y->tau)
      continue;

    if (c_cand == b)
      k = target_count;

    for (z_cand = c_cand->parent; z_cand; z_cand = z_cand->parent)
      if (stree->pptable[y->node_index][z_cand->node_index])
        break;  /* y is father of AC after move */

    target_weight[target_count] = 1;

    for (tmp = y; tmp != z_cand; tmp = tmp->parent)
      target_weight[target_count]++;

    for (tmp = c_cand; tmp != z_cand; tmp = tmp->parent)
      target_weight[target_count]++;

    target_weight[target_count] = 1 / target_weight[target_count];
    sum += target_weight[target_count++];
  }

  lnacceptance += log(target_weight[k] / sum);

  double newprior = lnprior_species_model(stree);

  lnacceptance += newprior - oldprior;

  bl_list = __gt_nodes;
  snode_contrib = snode_contrib_space;

  double logpr_notheta = stree->notheta_logpr;
  for (i = 0; i < stree->locus_count; ++i)
  {
    gtree_list[i]->old_logl = gtree_list[i]->logl;

    if (opt_debug_full)
    {
      k=0;
      for (j = 0; j < gtree_list[i]->tip_count + gtree_list[i]->inner_count; ++j)
      {
        gnode_t * tmp = gtree_list[i]->nodes[j];
        if (tmp->parent)
        {
          tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree_list[i]->edge_count,
                                               tmp->pmatrix_index);
          bl_list[k++] = tmp;
        }
      }
    
    
      locus_update_matrices(loci[i],gtree_list[i],bl_list,stree,i,k);
    
      gtree_all_partials(gtree_list[i]->root,bl_list,&k);
      for (j = 0; j < k; ++j)
      {
        bl_list[j]->clv_index = SWAP_CLV_INDEX(gtree_list[i]->tip_count,
                                               bl_list[j]->clv_index);
        if (opt_scaling)
          bl_list[j]->scaler_index = SWAP_SCALER_INDEX(gtree_list[i]->tip_count,
                                                       bl_list[j]->scaler_index);
      }
    
      locus_update_partials(loci[i],bl_list,k);
    
      /* compute log-likelihood */
      assert(!gtree_list[i]->root->parent);
      gtree_list[i]->logl = locus_root_loglikelihood(loci[i],
                                                     gtree_list[i]->root,
                                                     loci[i]->param_indices,
                                                     NULL);
    }
    else

    /* TODO: Perhaps change the below check to if (__mark_count[i]) ?? */
    //if (moved_count[i])
    if (__mark_count[i])
    {
      /* update branch lengths and transition probability matrices */
      for (j = 0; j < (unsigned int)__mark_count[i]; ++j)
      {
        bl_list[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree_list[i]->edge_count,
                                                    bl_list[j]->pmatrix_index);
      }

      locus_update_matrices(loci[i], gtree_list[i], bl_list, stree, i, __mark_count[i]);

      /* retrieve all nodes whose partials must be updates */
      unsigned int partials_count;
      gnode_t ** partials = gtree_list[i]->travbuffer;
      gtree_return_partials(gtree_list[i]->root,
                            gtree_list[i]->travbuffer,
                            &partials_count);

      /* point to the double-buffered partials space */
      for (j = 0; j < partials_count; ++j)
      {
        partials[j]->clv_index = SWAP_CLV_INDEX(gtree_list[i]->tip_count,
                                                partials[j]->clv_index);
        if (opt_scaling)
          partials[j]->scaler_index = SWAP_SCALER_INDEX(gtree_list[i]->tip_count,
                                                     partials[j]->scaler_index);
      }

      /* update conditional probabilities (partials) of affected nodes */
      locus_update_partials(loci[i], partials, partials_count);

      /* evaluate log-likelihood */
      gtree_list[i]->logl = locus_root_loglikelihood(loci[i],
                                                     gtree_list[i]->root,
                                                     loci[i]->param_indices,
                                                     NULL);
    }

    if (opt_est_theta)
      gtree_list[i]->old_logpr = gtree_list[i]->logpr;

    if (opt_debug_full)
    {
      /* This recomputes the gene tree probabilities from scratch. It can be used to verify that
         the code below, which only computes the gene tree probability for the changed components,
         is correct. */
    
      if (opt_est_theta)
        gtree_list[i]->old_logpr = gtree_list[i]->logpr;
    
      double logpr;
      
      if (opt_migration)
        logpr = gtree_logprob_mig(stree, gtree_list[i], loci[i]->heredity[0], i, thread_index);
      else
        logpr = gtree_logprob(stree, loci[i]->heredity[0], i, thread_index);
    
    
      if (opt_est_theta)
        gtree_list[i]->logpr = logpr;
      else
      {
        if (i == opt_locus_count - 1) 
          logpr += stree->notheta_hfactor+stree->notheta_sfactor;
        logpr_notheta = logpr;
      }
    }
    else
    {
      if (opt_migration && !opt_est_theta)
        fatal("Integrating out thetas not yet implemented for IM model");

       /* locate additional populations that need to be updated */

       /* find and mark those populations whose number of incoming lineages has
          changed due to the reset_gene_leaves_count() call, but were previously
          not marked for log-probability contribution update */
      for (j = 0; j < snode_contrib_count[i]; ++j)
        snode_contrib[j]->mark[thread_index] |= FLAG_POP_UPDATE;
      for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      {
        snode_t * snode = stree->nodes[j];
        if (!(snode->mark[thread_index] & FLAG_POP_UPDATE) &&
            (snode->seqin_count[i] != original_stree->nodes[j]->seqin_count[i]))
          snode_contrib[snode_contrib_count[i]++] = snode;
      }

      /* now update the log-probability contributions for the affected, marked
         populations */
      for (j = 0; j < snode_contrib_count[i]; ++j)
      {
        if (opt_est_theta)
          gtree_list[i]->logpr -= snode_contrib[j]->logpr_contrib[i];
        else
          logpr_notheta -= snode_contrib[j]->notheta_logpr_contrib;

        double xtmp;
        
        if (opt_migration)
          xtmp = gtree_update_logprob_contrib_mig(snode_contrib[j],
                                                  stree,
                                                  gtree_list[i],
                                                  loci[i]->heredity[0],
                                                  i,
                                                  thread_index);
        else
          xtmp = gtree_update_logprob_contrib(snode_contrib[j],
                                              loci[i]->heredity[0],
                                              i,
                                              thread_index);

        if (opt_est_theta)
          gtree_list[i]->logpr += snode_contrib[j]->logpr_contrib[i];
        else
          logpr_notheta += xtmp;
      }
    }

    /* 
       TODO: Several improvements can be made here
       1. Call this function only for the gene trees that are modified
       2. Only re-compute the prior for the affected gene tree nodes/edges
    */
    if (opt_clock == BPP_CLOCK_CORR)
    {
      double new_prior_rates = lnprior_rates(gtree_list[i],stree,i);
      lnacceptance += new_prior_rates - gtree_list[i]->lnprior_rates;
      gtree_list[i]->lnprior_rates = new_prior_rates;
    }

    /* reset markings on affected populations */
    for (j = 0; j < snode_contrib_count[i]; ++j)
      snode_contrib[j]->mark[thread_index] = 0;


    bl_list += __mark_count[i];

    for (j = 0; j < gtree_list[i]->tip_count + gtree_list[i]->inner_count; ++j)
      gtree_list[i]->nodes[j]->mark = 0;

    snode_contrib += stree->tip_count + stree->inner_count;

    if (opt_est_theta)
      lnacceptance += gtree_list[i]->logpr - gtree_list[i]->old_logpr +
                      gtree_list[i]->logl - gtree_list[i]->old_logl;
    else
      lnacceptance += gtree_list[i]->logl - gtree_list[i]->old_logl;
  }

  if (!opt_est_theta)
  {
    lnacceptance += logpr_notheta - stree->notheta_logpr;
    stree->notheta_logpr = logpr_notheta;
  }

  if (opt_debug_sspr)
    printf("[Debug] (SSPR) lnacceptance = %f\n", lnacceptance);

  /* in case of acceptance, cloned trees are substituted with the original ones,
     and species tree nodes are re-labeled, but all this is done in method_01.c
  */
  //return (lnacceptance >= 0 || legacy_rndu() < exp(lnacceptance));
  return (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance));
}

double lnprior_rates(gtree_t * gtree, stree_t * stree, long msa_index)
{
  /* Calculates the log of the prior of branch rates under the two rate-drift
     models:
     the independent rates (opt_clock=BPP_CLOCK_IND) and the geometric Brownian
     motion model (opt_clock=BPP_CLOCK_CORR). 

     BPP_CLOCK_IND:
     The algorithm cycles through the branches and adds up the log probabilities

     BPP_CLOCK_CORR:
     The root rate is mu or mean_rate. The algorithm cycles through the
     ancestral nodes and deals with the two daughter branches

  */

  long i;
  double z,r;
  double mui,nui;
  double alpha,beta;
  double logpr = 0;
  snode_t * snode;
  
  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  assert(opt_clock == BPP_CLOCK_IND || opt_clock == BPP_CLOCK_CORR);

  /* TODO: Implement autocorrelated clock */
  if (opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_GAMMA)
  {
    assert(!opt_msci);
    double v = gtree->rate_nui;

    for (i = stree->tip_count; i < total_nodes; ++i)
    {
      snode = stree->nodes[i];
      if (!snode->parent)
        assert(snode->brate[msa_index] == gtree->rate_mui);

      double m = snode->brate[msa_index];
      double alpha = m*m / v;
      double beta = alpha / m;
      double r1 = snode->left->brate[msa_index];
      double r2 = snode->right->brate[msa_index];
      logpr += -2 * lgamma(alpha) + 2 * alpha*log(beta) - beta*(r1 + r2) + (alpha - 1)*log(r1*r2);
    }
  }
  else if (opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
  {
    double Tinv[4], t1, t2, tA, detT;
    double rA, r1, r2, y1, y2, zz;

    assert(!opt_msci);

    mui = gtree->rate_mui;
    nui = gtree->rate_nui;

    for (i = stree->tip_count; i < total_nodes; ++i)
    {
      snode = stree->nodes[i];

      if (opt_msci && snode->hybrid)
      {
        assert(0);

        /* TODO: It is not clear how tA is calculated for hybridizations */
        if (node_is_hybridization(snode) && !snode->htau) continue;
        if (node_is_bidirection(snode) && node_is_mirror(snode)) continue;
      }

      tA = snode->parent ? (snode->parent->tau - snode->tau) / 2 : 0;
      t1 = (snode->tau - snode->left->tau) / 2;
      t2 = (snode->tau - snode->right->tau) / 2;

      detT = t1*t2 + tA*(t1+t2);
      Tinv[0] = (tA+t2) / detT;
      Tinv[1] = Tinv[2] = -tA / detT;
      Tinv[3] = (tA+t1) / detT;

      /* root node should have mui anyway */
      rA = snode->parent ? snode->brate[msa_index] : mui;
      r1 = snode->left->brate[msa_index];
      r2 = snode->right->brate[msa_index];
      y1 = log(r1/rA) + (tA+t1)*nui / 2;
      y2 = log(r2/rA) + (tA+t2)*nui / 2;
      zz = (y1*y1*Tinv[0] + 2*y1*y2*Tinv[1] + y2*y2*Tinv[3]);
      logpr -= zz / (2*nui) + log(detT*nui*nui)/2 + log(r1*r2);
    }
    assert(!opt_msci);
    logpr -= 0.5*log(2*BPP_PI) * stree->inner_count*2;
  }
  else if (opt_clock == BPP_CLOCK_IND && opt_rate_prior == BPP_BRATE_PRIOR_GAMMA)
  {
    /* clock == BPP_CLOCK_IND and gamma rate prior */

    /*
      We compute the following:
      
      Pr = \Prod_{i=1}^{k} \frac{\beta^\alpha}{\Gamma(\alpha)}
           r_{i}^{\alpha-1} e^{-\beta r_i}
      where k is the number of proposed rates (rates_count)
    */

    mui = gtree->rate_mui;
    nui = gtree->rate_nui;

    alpha = mui * mui / nui;
    beta  = mui / nui;
    long rates_count = 0;

    for (i = 0; i < total_nodes; ++i)
    {
      snode = stree->nodes[i];

      if (opt_msci && snode->hybrid)
      {
        if (node_is_hybridization(snode) && !snode->htau) continue;
        if (node_is_bidirection(snode) && node_is_mirror(snode)) continue;
      }

      r = snode->brate[msa_index];
      logpr += -beta*r + (alpha-1)*log(r);

      ++rates_count;
    }

    logpr += (alpha*log(beta) - lgamma(alpha)) * rates_count;
  }
  else if (opt_clock == BPP_CLOCK_IND && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
  {
    /* clock == BPP_CLOCK_IND and log-normal rate prior */

    /*
      We compute the following:

      Pr = \Prod_{i=1}^{k} \frac{1}{r_i * \sqrt{\sigma^2}}
           \exp{-\frac{1}{2\sigma^2}(\log\frac{r_i}{\mu} + \frac{sigma^2}{2})^2}

      where k is the number of proposed rates (rates_count)

    */

    mui = gtree->rate_mui;
    nui = gtree->rate_nui;
    long rates_count = 0;

    double logmui = log(mui);
    for (i = 0; i < total_nodes; ++i)
    {
      snode = stree->nodes[i];

      if (opt_msci && snode->hybrid)
      {
        if (node_is_hybridization(snode) && !snode->htau) continue;
        if (node_is_bidirection(snode) && node_is_mirror(snode)) continue;
      }

      double logr = log(snode->brate[msa_index]);
      z = logr - logmui + nui/2;
      logpr += -(z*z) / (2*nui) - logr;

      ++rates_count;
    }
    logpr -= 0.5*log(2*BPP_PI*nui) * rates_count;
  }

  return logpr;
}

double prop_locusrate_nui(gtree_t ** gtree,
                          stree_t * stree,
                          locus_t ** locus,
                          long thread_index)
{
  unsigned int j;
  unsigned int total_nodes;
  long i;
  long accepted = 0;
  double old_nu, old_lognu;
  double new_nu, new_lognu;
  double sum_old, sum_new;
  double alpha, beta;
  double logl = 0;
  double new_prior_rates = 0;
  double terma, termb;
  double lnacceptance = 0;
  gnode_t ** gnodeptr;

  assert(thread_index == 0);

  /* TODO: opt_finetune_locusrate */
  alpha = opt_vi_alpha;
  beta  = opt_vi_alpha / stree->locusrate_nubar;

  sum_old = sum_new = terma = termb = 0;
  /* TODO: Note Gamma-Dirichlet prior cannot be parallelized as we have to
     compute the sum of nu_i across loci */
  if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR ||
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR)
  {
    sum_old = stree->nui_sum;

    terma = opt_vi_alpha*opt_locus_count;
    termb = opt_vbar_beta/opt_locus_count;
  }

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (opt_clock == BPP_CLOCK_GLOBAL)
    {
      /* if molecular clock then swap pmatrices */

      gnodeptr = gtree[i]->nodes;
      total_nodes = gtree[i]->tip_count+gtree[i]->inner_count;

      for (j = 0; j < total_nodes; ++j)
        if (gnodeptr[j]->parent)
          gnodeptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                    gnodeptr[j]->pmatrix_index);
    }

    old_nu = gtree[i]->rate_nui;
    old_lognu = log(old_nu);

    double r = old_lognu + opt_finetune_nui * legacy_rnd_symmetrical(thread_index);
    new_lognu = reflect(r,-99,99,thread_index);
    gtree[i]->rate_nui = new_nu = exp(new_lognu);

    lnacceptance = new_lognu - old_lognu;

    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR ||
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR)
    {
      /* Gamma-Dirichlet prior */

      sum_new = sum_old + new_nu - old_nu;
      lnacceptance += (opt_vbar_alpha - terma)*log(sum_new / sum_old) -
                      termb*(sum_new - sum_old) +
                      (opt_vi_alpha - 1)*(new_lognu - old_lognu);
    }
    else if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    {
      /* Hierarchical conditional iid prior */

      /* gamma prior ratio */
      lnacceptance += (alpha-1)*(new_lognu-old_lognu) - beta*(new_nu-old_nu);
    }
    else
      assert(0);

    if (opt_clock == BPP_CLOCK_GLOBAL)
    {
      /* if molecular clock then recompute pmatrices, CLVs and log-L */
      locus_update_all_matrices(locus[i],gtree[i],stree,i);

      gnodeptr = gtree[i]->nodes;
      total_nodes = gtree[i]->tip_count+gtree[i]->inner_count;

      for (j = gtree[i]->tip_count; j < total_nodes; ++j)
      {
        gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                gnodeptr[j]->clv_index);
        if (opt_scaling)
          gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                     gnodeptr[j]->scaler_index);
      }

      locus_update_all_partials(locus[i],gtree[i]);

      logl = locus_root_loglikelihood(locus[i],
                                      gtree[i]->root,
                                      locus[i]->param_indices,
                                      NULL);

      lnacceptance += logl - gtree[i]->logl;
    }
    else
    {
      /* if relaxed clock then update rates prior */

      new_prior_rates = lnprior_rates(gtree[i],stree,i);
      lnacceptance += new_prior_rates - gtree[i]->lnprior_rates;
    }

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
    {
      /* accepted */
      accepted++;

      if (opt_clock == BPP_CLOCK_GLOBAL)
      {
        /* update log-likelihood */
        gtree[i]->logl = logl;
      }
      else
      {
        /* relaxed clock */

        /* update rates prior */
        gtree[i]->lnprior_rates = new_prior_rates;
      }

      if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR ||
          opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR)
        sum_old = sum_new;
    }
    else
    {
      /* rejected */

      gtree[i]->rate_nui = old_nu;

      if (opt_clock == BPP_CLOCK_GLOBAL)
      {
        /* reset selected locus */
        gnodeptr = gtree[i]->nodes;
        total_nodes = gtree[i]->tip_count+gtree[i]->inner_count;

        for (j = 0; j < total_nodes; ++j)
        {
          if (gnodeptr[j]->parent)
            gnodeptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                         gnodeptr[j]->pmatrix_index);
        }
        for (j = gtree[i]->tip_count; j < total_nodes; ++j)
        {
          gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                  gnodeptr[j]->clv_index);
          if (opt_scaling)
            gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                          gnodeptr[j]->scaler_index);
        }
      }
    }
  }
  if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR ||
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR)
    stree->nui_sum  = sum_old;
  return ((double)accepted / opt_locus_count);
}

double prop_locusrate_mui(gtree_t ** gtree,
                          stree_t * stree,
                          locus_t ** locus,
                          long thread_index)
{
  unsigned int j;
  unsigned int total_nodes;
  long i;
  long accepted = 0;
  double old_mui, old_logmui;
  double new_mui, new_logmui;
  double sum_old, sum_new;
  double alpha;
  double lnacceptance = 0;
  double logl = 0;
  double new_prior_rates = 0;
  double terma, termb;
  gnode_t ** gnodeptr;


  /* TODO: Separate this routine into different clock cases.
     
     If opt_clock == BPP_CLOCK_GLOBAL, then we propose the locus rate, which
     changes the likelihood but there is no prior for branch rates.

     If opt_clock != BPP_CLOCK_GLOBAL, then we propose mu and nu or loci,
     which changes the rate prior but not the likelihood

  */
  alpha = opt_mui_alpha;

  sum_old = sum_new = terma = termb = 0;
  /* TODO: Note Gamma-Dirichlet prior cannot be parallelized as we have to
     compute the sum of mu_i across loci */
  if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR)
  {
    sum_old = 0;
    for (i = 0; i < opt_locus_count; ++i)
      sum_old += gtree[i]->rate_mui;

    terma = opt_mui_alpha*opt_locus_count;
    termb = opt_mubar_beta/opt_locus_count;
  }

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (opt_clock == BPP_CLOCK_GLOBAL || opt_clock == BPP_CLOCK_CORR)
    {
      /* if molecular clock then swap pmatrices */

      gnodeptr = gtree[i]->nodes;
      total_nodes = gtree[i]->tip_count+gtree[i]->inner_count;

      for (j = 0; j < total_nodes; ++j)
        if (gnodeptr[j]->parent)
          gnodeptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                    gnodeptr[j]->pmatrix_index);
    }

    old_mui = gtree[i]->rate_mui;
    old_logmui = log(old_mui);

    double r = old_logmui + opt_finetune_mui * legacy_rnd_symmetrical(thread_index);
    new_logmui = reflect(r,-99,99,thread_index);
    gtree[i]->rate_mui = new_mui = exp(new_logmui);

    lnacceptance = new_logmui - old_logmui;

    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR)
    {
      /* Gamma-Dirichlet prior */

      sum_new = sum_old + new_mui - old_mui;
      lnacceptance += (opt_mubar_alpha - terma)*log(sum_new / sum_old) -
                      termb*(sum_new - sum_old) +
                      (opt_mui_alpha - 1)*(new_logmui - old_logmui);
    }
    else if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    {
      /* Hierarchical conditional iid prior */

      /* gamma prior ratio */
      lnacceptance += (alpha-1)*log(new_mui/old_mui) - 
                      (opt_mui_alpha/stree->locusrate_mubar)*(new_mui-old_mui);
    }
    else
    {
      /* Note: BPP_LOCRATE_PRIOR_DIR is not processed in this routine. Instead
         it is handled in function prop_locusrate_and_heredity(...) */
      assert(0);
    }

    if (opt_clock == BPP_CLOCK_GLOBAL || opt_clock == BPP_CLOCK_CORR)
    {
      if (opt_clock == BPP_CLOCK_CORR)
        stree->root->brate[i] = gtree[i]->rate_mui;

      /* if molecular clock then recompute pmatrices, CLVs and log-L */
      locus_update_all_matrices(locus[i],gtree[i],stree,i);

      gnodeptr = gtree[i]->nodes;
      total_nodes = gtree[i]->tip_count+gtree[i]->inner_count;

      for (j = gtree[i]->tip_count; j < total_nodes; ++j)
      {
        gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                gnodeptr[j]->clv_index);
        if (opt_scaling)
          gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                     gnodeptr[j]->scaler_index);
      }

      locus_update_all_partials(locus[i],gtree[i]);

      logl = locus_root_loglikelihood(locus[i],
                                      gtree[i]->root,
                                      locus[i]->param_indices,
                                      NULL);

      lnacceptance += logl - gtree[i]->logl;
    }

    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      /* if relaxed clock then update rates prior */

      new_prior_rates = lnprior_rates(gtree[i],stree,i);
      lnacceptance += new_prior_rates - gtree[i]->lnprior_rates;
    }

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
    {
      /* accepted */
      accepted++;

      if (opt_clock == BPP_CLOCK_GLOBAL || opt_clock == BPP_CLOCK_CORR)
      {
        /* update log-likelihood */
        gtree[i]->logl = logl;
      }

      if (opt_clock != BPP_CLOCK_GLOBAL)
      {
        /* relaxed clock */

        /* update rates prior */
        gtree[i]->lnprior_rates = new_prior_rates;
      }

      if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR)
        sum_old = sum_new;
    }
    else
    {
      /* rejected */
      gtree[i]->rate_mui = old_mui;
      if (opt_clock == BPP_CLOCK_CORR)
        stree->root->brate[i] = old_mui;

      if (opt_clock == BPP_CLOCK_GLOBAL || opt_clock == BPP_CLOCK_CORR)
      {
        /* reset selected locus */
        gnodeptr = gtree[i]->nodes;
        total_nodes = gtree[i]->tip_count+gtree[i]->inner_count;

        for (j = 0; j < total_nodes; ++j)
        {
          if (gnodeptr[j]->parent)
            gnodeptr[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                         gnodeptr[j]->pmatrix_index);
        }

        for (j = gtree[i]->tip_count; j < total_nodes; ++j)
        {
          gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                  gnodeptr[j]->clv_index);
          if (opt_scaling)
            gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                          gnodeptr[j]->scaler_index);
        }
      }
    }
  }
  return ((double)accepted / opt_locus_count);
}

long prop_locusrate_mubar(stree_t * stree, gtree_t ** gtree)
{
  long i;
  long accepted = 0;
  double old_rate, old_lograte;
  double new_rate, new_lograte;
  double lnacceptance = 0;

  const long thread_index = 0;

  old_rate = stree->locusrate_mubar;
  old_lograte = log(old_rate);

  /* now sample the universal mean and variance */
  double r = old_lograte + opt_finetune_mubar *
             legacy_rnd_symmetrical(thread_index);
  new_lograte = reflect(r,-99,99,thread_index);
  stree->locusrate_mubar = new_rate = exp(new_lograte);

  lnacceptance = new_lograte - old_lograte;

  lnacceptance += (opt_mubar_alpha-1)*log(new_rate/old_rate) -
                  opt_mubar_beta*(new_rate-old_rate);
  double a = opt_mui_alpha;
  double bnew = opt_mui_alpha / new_rate;
  double bold = opt_mui_alpha / old_rate;
  lnacceptance += opt_locus_count*a*log(bnew / bold);
  for (i = 0; i < opt_locus_count; ++i)
     lnacceptance -= (bnew - bold)*gtree[i]->rate_mui;

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */
    accepted++;
  }
  else
  {
    /* rejected */
    stree->locusrate_mubar = old_rate;
  }

  return accepted;
}

long prop_locusrate_nubar(stree_t * stree, gtree_t ** gtree)
{
  long i;
  long accepted = 0;
  double old_nu, old_lognu;
  double new_nu, new_lognu;
  double lnacceptance = 0;

  const long thread_index = 0;

  old_nu = stree->locusrate_nubar;
  old_lognu = log(old_nu);

  /* now sample the universal mean and variance */
  double r = old_lognu + opt_finetune_nubar *
             legacy_rnd_symmetrical(thread_index);
  new_lognu = reflect(r,-99,99,thread_index);
  stree->locusrate_nubar = new_nu = exp(new_lognu);

  lnacceptance = new_lognu - old_lognu;

  lnacceptance += (opt_vbar_alpha-1)*(new_lognu-old_lognu) -
                  opt_vbar_beta*(new_nu-old_nu);
  double a = opt_vi_alpha;
  double bnew = opt_vi_alpha / new_nu;
  double bold = opt_vi_alpha / old_nu;
  lnacceptance += opt_locus_count*a*log(bnew / bold);
  for (i = 0; i < opt_locus_count; ++i)
     lnacceptance -= (bnew - bold)*gtree[i]->rate_nui;

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */
    accepted++;
  }
  else
  {
    /* rejected */
    stree->locusrate_nubar = old_nu;
  }

  return accepted;
}

static double prior_logratio_rates_corr_ln(snode_t * node_changed,
                                           double old_rate,
                                           double new_rate,
                                           double variance,
                                           long msa_index)
{
  long i,j,cycles;
  double tA,t1,t2;
  double rA,r1,r2,y1,y2;
  double detT;
  double zz;
  double Tinv[4];
  double rates[2] = { old_rate, new_rate };
  double logterm[2] = {0,0};


  /* log-normal correlated clock */

  /* follows closely the notation in:
  
     Rannala B., Yang Z.
     Inferring Speciation Times under an Episodic Molecular Clock.
     Systematic Biology, 56(3):453-466, 2007.

  */ 

  /* for correlated log-normal model, rate at root branch is copied from mu_i, 
     so node_changed cannot be the root. */
  assert(opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL);
  assert(node_changed->parent);

  snode_t * snodes[2] = { node_changed->parent, node_changed };

  cycles = (node_changed->left) ? 2 : 1;

  /* cycle through parent and node or just parent if node is a tip */
  for (i = 0; i < cycles; ++i)
  {
    snode_t * node = snodes[i];   /* this is inode */

    tA = node->parent ? (node->parent->tau - node->tau) / 2 : 0;
    t1 = (node->tau - node->left->tau) / 2;
    t2 = (node->tau - node->right->tau) / 2;

    detT = t1*t2 + tA*(t1+t2);
    Tinv[0] = (tA + t2) / detT;
    Tinv[1] = Tinv[2] = -tA / detT;
    Tinv[3] = (tA+t1) / detT;

    /* now cycle through old and new rate */
    for (j = 0; j < 2; ++j)
    {
      node_changed->brate[msa_index] = rates[j];

      if (node->parent)
        rA = node->brate[msa_index];
      else
        rA = node->brate[msa_index];   /* we assume mu_i is in array, if node is root */
        
      assert(node->left && node->right);
      r1 = node->left->brate[msa_index];
      r2 = node->right->brate[msa_index];

      y1 = log(r1/rA) + (tA+t1)*variance/2;
      y2 = log(r2/rA) + (tA+t2)*variance/2;
      zz = y1*y1*Tinv[0] + 2*y1*y2*Tinv[1] + y2*y2*Tinv[3];
      zz = zz/(2*variance) + log(r1*r2);
      logterm[j] += zz;
    }
  }

  return (logterm[0] - logterm[1]);
}

static double prior_logratio_rates_corr(snode_t * node,
                                        gtree_t * gtree,
                                        double old_rate,
                                        double new_rate,
                                        long msa_index)
{
  /* Calculates the prior ratio when one branch rate is changed.
  If we are updating a terminal (tip) branch, we sum over 1 term (Eq. 7 in RY2007).
  If we are updating an inner branch, we sum over 2 terms.
  Rannala B., Yang. Z.: Speciation time and episodic molecular clock.
  Systematic Biology, 56(3): 453-466, 2007.
  */

  double logratio = 0;

  assert(opt_clock == BPP_CLOCK_CORR);
  assert(!opt_msci);
  assert(node->parent);

  if (opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
  {
    /* log-normal */
    logratio = prior_logratio_rates_corr_ln(node,
                                            old_rate,
                                            new_rate,
                                            gtree->rate_nui,
                                            msa_index);
  }
  else
  {
    /* gamma */

    /**** Ziheng-2020-04-06 ****/
    double m = node->parent->brate[msa_index];
    double v = gtree->rate_nui;
    double alpha = m*m / v;
    double beta = alpha / m;

    logratio = -beta*(new_rate - old_rate) + (alpha - 1)*log(new_rate / old_rate);
    if (node->left)
    {
      assert(node->right);

      double alpha = old_rate*old_rate / v;
      double beta = alpha / old_rate;
      double alphanew = new_rate*new_rate / v;
      double betanew = alphanew / new_rate;
      double r1 = node->left->brate[msa_index];
      double r2 = node->right->brate[msa_index];

      logratio += -2*lgamma(alphanew) + 2*alphanew*log(betanew) - betanew*(r1+r2) + (alphanew-1)*log(r1*r2);
      logratio -= -2*lgamma(alpha)    + 2*alpha*log(beta)       - beta*(r1+r2)    + (alpha-1)*log(r1*r2);

    }
  }

  return logratio;
}

static double prior_logratio_rates_iid(gtree_t * gtree,
                                       double old_rate,
                                       double new_rate)
{
  double a,b;
  double zold, znew;
  double ratio = 0;

  assert(opt_clock == BPP_CLOCK_IND);

  if (opt_clock == BPP_CLOCK_IND && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
  {
    /* iid rates and log-normal rate prior */

    zold = log(old_rate / gtree->rate_mui) + gtree->rate_nui / 2;
    znew = log(new_rate / gtree->rate_mui) + gtree->rate_nui / 2;

    ratio = -log(new_rate / old_rate) - 
            (znew*znew - zold*zold) / (2*gtree->rate_nui);
  }
  else if (opt_clock == BPP_CLOCK_IND && opt_rate_prior == BPP_BRATE_PRIOR_GAMMA)
  {
    /* iid rates and Gamma rate prior */

    a = gtree->rate_mui * gtree->rate_mui / gtree->rate_nui;
    b = gtree->rate_mui / gtree->rate_nui;

    ratio = -b*(new_rate - old_rate) + (a-1)*log(new_rate / old_rate);
  }
  else
    assert(0);

  return ratio;
}

static long fill_travbuffer_and_mark(gtree_t * gtree,
                                     stree_t * stree,
                                     gnode_t ** buffer,
                                     snode_t * target)
{
  /* 
     Fill buffer with nodes whose parent edges intersect with target
     population and return the number of matches. Also mark that the parents of
     those nodes need their CLVs updated.
  
  */

  long i;
  long count;
  snode_t * start;
  snode_t * end;

  count = 0;

  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    if (!x->parent) continue;

    start = x->pop;
    end = x->parent->pop;

    while (start != end)
    {
      assert(start && start->parent);

      if (start == target)
        break;

      start = start->parent;

      if (start->hybrid)
      {
        assert(!node_is_mirror(start));

        unsigned int hindex = GET_HINDEX(stree,start);
        assert(hindex >= 0 && hindex < stree->hybrid_count);

        /* find correct parent node according to hpath flag */
        assert(x->hpath[hindex] != BPP_HPATH_NONE);
        assert(start->left);
        if (x->hpath[hindex] == BPP_HPATH_RIGHT)
          start = start->hybrid;
      }
    }
    if (start == target)
    {
      /* add node x to the buffer to indicate that we want to update its
         branch length. Also mark that its parent CLV needs update */
      buffer[count++] = x;
      x->parent->mark |= FLAG_PARTIAL_UPDATE;
    }
  }
  return count;
}

static long prop_branch_rates(gtree_t * gtree,
                              stree_t * stree,
                              locus_t * locus,
                              unsigned int msa_index,
                              long * prop_count,
                              long thread_index)
{
  long j,k;
  long accepted = 0;
  long proposal_count = 0;
  long branch_count = 0;
  double diff;
  double old_rate, new_rate;
  double old_lograte, new_lograte;
  double lnacceptance = 0;
  snode_t * node;

  assert(opt_clock != BPP_CLOCK_GLOBAL);

  unsigned int partials_count;
  gnode_t ** partials = NULL;
  gnode_t ** updatelist = __gt_nodes + __gt_nodes_index[msa_index];

  for (j = 0; j < stree->tip_count+stree->inner_count+stree->hybrid_count; ++j)
  {
    node = stree->nodes[j];

    /* if root node and we use correlated clock, skip */
    if (!node->parent && opt_clock == BPP_CLOCK_CORR) continue;

    /* Mirror nodes in bidirectional introgression */
    if (opt_msci && node->hybrid)
    {
      if (node_is_hybridization(node) && !node->htau) continue;
      if (node_is_bidirection(node) && node_is_mirror(node)) continue;
    }
    proposal_count++;

    old_rate = node->brate[msa_index];
    old_lograte = log(old_rate);

    double r = old_lograte + opt_finetune_branchrate *
               legacy_rnd_symmetrical(thread_index);
    new_lograte = reflect(r,-99,99,thread_index);
    node->brate[msa_index] = new_rate = exp(new_lograte);

    lnacceptance = new_lograte - old_lograte;

    /* obtain a list of edges (represented by nodes) that intersect with
       species tree node and mark nodes whose CLV need update */
    if (opt_usedata)
      branch_count = fill_travbuffer_and_mark(gtree,stree,updatelist,node);

    /* if at least one gene tree branch needs updating, we must recompute the
       log-likelihood */
    gtree->old_logl = gtree->logl;
    if (branch_count)
    {
      /* swap pmatrices */
      for (k = 0; k < branch_count; ++k)
        updatelist[k]->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,
                                                       updatelist[k]->pmatrix_index);
      /* update necessary p-matrices */
      locus_update_matrices(locus,gtree,updatelist,stree,msa_index,branch_count);

      /* get the list of nodes for which CLVs must be reverted, i.e. all marked
         nodes and all nodes whose left or right subtree has at least one marked
         node */
      gtree_return_partials(gtree->root,
                            gtree->travbuffer,
                            &partials_count);
      partials = gtree->travbuffer;
      
      /* remove flags */
      for (k = 0 ; k < branch_count; ++k)
        updatelist[k]->parent->mark = 0;
        
      for (k = 0; k < partials_count; ++k)
      {
        partials[k]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                partials[k]->clv_index);
        if (opt_scaling)
          partials[k]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                        partials[k]->scaler_index);
      }

      /* update partials */
      locus_update_partials(locus, partials, partials_count);

      /* evaulate log-likelihood */
      double logl = locus_root_loglikelihood(locus,
                                             gtree->root,
                                             locus->param_indices,
                                             NULL);
      
      lnacceptance += logl - gtree->logl;

      gtree->logl = logl;
    }

    if (opt_clock == BPP_CLOCK_CORR)
      diff = prior_logratio_rates_corr(node,gtree,old_rate,new_rate,msa_index); 
    else
      diff = prior_logratio_rates_iid(gtree,old_rate,new_rate);

    lnacceptance += diff;

    if (opt_debug_br)
      printf("[Debug] (br) lnacceptance = %f\n", lnacceptance);

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
    {
      /* accepted */
      accepted++;
      gtree->lnprior_rates += diff;
    }
    else
    {
      /* rejected */

      /* restore old log-l */
      gtree->logl = gtree->old_logl;

      /* swap to old p-matrices */
      for (k = 0; k < branch_count; ++k)
        updatelist[k]->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,
                                                       updatelist[k]->pmatrix_index);

      /* need to reset clv indices to point to the old clv buffer */
      if (branch_count)
      {
        for (k = 0; k < partials_count; ++k)
        {
          partials[k]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                  partials[k]->clv_index);
          if (opt_scaling)
            partials[k]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                          partials[k]->scaler_index);
        }
      }

      /* now reset rate and pmatrices */
      node->brate[msa_index] = old_rate;
    }
  }  /* end species tree loop */
  *prop_count = proposal_count;
  return accepted;
}

double prop_branch_rates_serial(gtree_t ** gtree,
                                stree_t * stree,
                                locus_t ** locus)
{
  unsigned int i;
  long proposal_count = 0;
  long prop_count;
  long accepted = 0;

  #ifdef DEBUG_THREADS
  /* simulate parallel execution with DEBUG_THREADS_COUNT threads, by assigning
     each locus the index of the corresponding random number generator it would
     be assigned in a parallel execution */
  assert(DEBUG_THREADS_COUNT <= stree->locus_count);
  long * indices = (long *)xmalloc((size_t)(stree->locus_count) * sizeof(long));
  long loci_per_thread = opt_locus_count / DEBUG_THREADS_COUNT;
  long loci_remaining = opt_locus_count % DEBUG_THREADS_COUNT;

  long load = 0;
  long thread_index = -1;
  for (i = 0; i < stree->locus_count; ++i)
  {
    if  (load == 0)
    {
      load = loci_per_thread + (loci_remaining > 0 ? 1 : 0);
      if (loci_remaining)
        --loci_remaining;
      ++thread_index;
    }
    indices[i] = thread_index;
    --load;
  }
  #endif

  for (i = 0; i < stree->locus_count; ++i)
  {
    prop_count = 0;
    #ifdef DEBUG_THREADS
    accepted += prop_branch_rates(gtree[i],stree,locus[i],i,&prop_count,indices[i]);
    #else
    accepted += prop_branch_rates(gtree[i],stree,locus[i],i,&prop_count,0);
    #endif
    proposal_count += prop_count;
  }

  #ifdef DEBUG_THREADS
  free(indices);
  #endif

  if (!accepted)
    return 0;

  return ((double)accepted/proposal_count);
}

void prop_branch_rates_parallel(gtree_t ** gtree,
                                stree_t * stree,
                                locus_t ** locus,
                                long locus_start,
                                long locus_count,
                                long thread_index,
                                long * p_proposal_count,
                                long * p_accepted)
{
  unsigned int i;
  long proposal_count = 0;
  long prop_count;
  long accepted = 0;
  
  assert(locus_start >= 0);
  assert(locus_count > 0);

  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    prop_count = 0;
    accepted += prop_branch_rates(gtree[i],stree,locus[i],i,&prop_count,thread_index);
    proposal_count += prop_count;
  }

  *p_proposal_count = proposal_count;
  *p_accepted = accepted;
}

static double logpdf_power(double y, double b, double lambda)
{
  return log(lambda/b) + (lambda-1)*log(1-y/b);
}

long snl_scale_clade(gnode_t * node,
                     snode_t ** rway,
                     double ytaunew,
                     double taufactor,
                     long msa_index,
                     long thread_index)
{
  long i, scaled_count = 1;
  assert(node);
  node->mark |= FLAG_BRANCH_UPDATE;

  if (!node->left)
    return 0;

  assert(node->left && node->right);
  node->time *= taufactor;
  node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
  
  #if 1
  //node->parent->mark |= FLAG_BRANCH_UPDATE;
  //node->left->mark |= FLAG_BRANCH_UPDATE;
  //node->right->mark |= FLAG_BRANCH_UPDATE;
  node->mark |= FLAG_PARTIAL_UPDATE;
  #endif

  if (node->time > ytaunew)
  {
    for (i = 1; rway[i]; ++i)
      if (node->time < rway[i]->tau)
        break;
    if (node->pop != rway[i-1])
    {
      /* remove gene node from list of coal events of its old population */
      unlink_event(node,msa_index);
      node->pop->event_count[msa_index]--;
      if (!opt_est_theta)
        node->pop->event_count_sum--;

      node->pop = rway[i-1];
      node->pop->mark[thread_index] |= FLAG_POP_UPDATE;

      dlist_item_append(node->pop->event[msa_index], node->event);
      node->pop->event_count[msa_index]++;
      if (!opt_est_theta)
        node->pop->event_count_sum++;
    }
  }

  scaled_count += snl_scale_clade(node->left,
                                  rway,
                                  ytaunew,
                                  taufactor,
                                  msa_index,
                                  thread_index);
  scaled_count += snl_scale_clade(node->right,
                                  rway,
                                  ytaunew,
                                  taufactor,
                                  msa_index,
                                  thread_index);

  return scaled_count;
}

static void snl_paint_nodes_recursive(stree_t * stree,
                                      snode_t * a,
                                      gnode_t * node)
{
  unsigned int pop_index;

  assert(node);


  if (!node->left)
  {
    assert(!node->right);
    pop_index = node->pop->node_index;

    if (stree->pptable[pop_index][a->node_index])
      node->mark |= SNL_PUREA;

    return;
  }

  snl_paint_nodes_recursive(stree,a,node->left);
  snl_paint_nodes_recursive(stree,a,node->right);

  if ((node->left->mark & SNL_PUREA) && (node->right->mark & SNL_PUREA))
    node->mark |= SNL_PUREA;
  else if ((node->left->mark & SNL_PUREA) || (node->right->mark & SNL_PUREA))
  {
    node->mark |= SNL_MOVED;
    assert(node->time > a->parent->tau);
  }
}

long snl_expand_and_shrink(stree_t * stree,
                           stree_t * original_stree,
                           gtree_t ** gtree_list,
                           locus_t ** loci,
                           long movetype,
                           int downwards,
                           double ytaunew,
                           double taufactor,
                           double * lnacceptance,
                           snode_t * a,
                           snode_t * b,
                           snode_t * c,
                           snode_t * y,
                           snode_t ** rway)
{
  unsigned int i,j,k;
  long target_count = 0;
  long source_count = 0;
  long scaled_count = 0;
  double tau0, tau0new;
  long ndspecies;
  gnode_t * gtarget;
  long thread_index = 0;

  double oldprior = lnprior_species_model(stree);

  if (opt_migration && !opt_est_theta)
    fatal("Integrating out thetas not yet implemented for IM model");

  tau0 = stree->root->tau;
  ndspecies = 1;
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    if (stree->nodes[i]->tau > 0) ndspecies++;
  assert(ndspecies > 2);

  gnode_t ** moved_nodes = moved_space;
  gnode_t ** gtarget_list = gtarget_temp_space;
  gnode_t ** gtarget_nodes = gtarget_space;
  snode_t ** snode_contrib = snode_contrib_space;     /* TODO 5.8.2020 */

  for (i = 0; i < stree->locus_count; ++i)
  {
    snode_contrib_count[i] = 0;

    gtree_t * gtree = gtree_list[i];

    /* paint gene tree nodes */
    snl_paint_nodes_recursive(stree, a, gtree->root);

    /* now identify Moved nodes */
    moved_count[i] = 0;
    for (j = gtree->tip_count; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      gnode_t * node = gtree->nodes[j];

      if (!(node->mark & SNL_MOVED)) continue;
      assert(((node->left->mark & SNL_PUREA) && !(node->right->mark & SNL_PUREA)) ||
             (!(node->left->mark & SNL_PUREA) && (node->right->mark & SNL_PUREA)));


      /* now we need to ensure that only one child has descendants in A only */
      gnode_t * pruned = NULL;
      gnode_t * intact = NULL;
      if (node->left->mark & SNL_PUREA)
      {
        pruned = node->left;
        intact = node->right;
      }
      else
      {
        pruned = node->right;
        intact = node->left;
      }

      moved_nodes[moved_count[i]] = node;
      pruned_nodes[moved_count[i]++] = pruned;   /* TODO 4.8.2020 do we need it ? */

      node->mark |= FLAG_PARTIAL_UPDATE;
      if (node->parent)
        node->parent->mark |= FLAG_PARTIAL_UPDATE; /* required as the parent will change */


      /* find the new population for node */
      double tnew = node->time*taufactor;
      for (k = 1; rway[k]; ++k)
        if (tnew < rway[k]->tau) break;
      snode_t * newpop = rway[k-1];

      /* 
         List target branches on gene tree and sample one for reattachment.
         Feasible target branch must pass newpop, but there are two
         complications because newpop is for the new stree while pptable is for
         the old stree.
         if newpop = Y in new stree, branch must pass C on old stree (C->AC=Y)
         if newpop = B in new stree, branch must pass AB on old stree (AB->B)
      */

      if (newpop == y) newpop = c;
      else if (newpop == b) newpop = y;

      target_count = 0;
      for (k = 0; k < gtree->tip_count + gtree->inner_count; ++k)
      {
        gnode_t * tmp = gtree->nodes[k];

        if (tmp->time >= tnew || (tmp->parent && tmp->parent->time <= tnew))
          continue;

        if (!(tmp->mark & SNL_PUREA) &&
            stree->pptable[tmp->pop->node_index][newpop->node_index])
          gtarget_list[target_count++] = tmp;
      }

      /* if no C seqs exist at the locus, there won't be target branch on gene
         tree, and the move is disabled. This will never happen when tau_Y
         increases, and may happen when tau_Y decreases in the move.
      */
      if (!target_count)
      {
        if (taufactor > 1 && !downwards)
          fatal("Internal error - this should not happen for taufactor>1 - "
                "please contact us");
        return 2;
      }
      
      gtarget = gtarget_list[(int)(target_count*legacy_rndu(thread_index))];

      /* reset target, if it has flag SNL_MOVED, by tracing towards tips until
         flag = BLACK */
      while (gtarget->mark & SNL_MOVED)
        gtarget = (gtarget->left->mark & SNL_PUREA) ?
                    gtarget->right : gtarget->left;

      /* save gtarget branch */
      assert(gtarget && !(gtarget->mark & (SNL_MOVED | SNL_PUREA)));
      gtarget_nodes[moved_count[i] - 1] = gtarget;

      source_count = 1;
      /* gsources_list[0] = intact; */
      for (k = 0; k < gtree->tip_count + gtree->inner_count; ++k)
      {
        gnode_t * tmp = gtree->nodes[k];

        if (tmp == intact || tmp == intact->parent)
          continue;  /* intact->father equals node */

        if (tmp->time >= node->time || 
            (tmp->parent && tmp->parent->time <= node->time))
          continue;

        if (!(tmp->mark & SNL_PUREA) &&
            stree->pptable[tmp->pop->node_index][node->pop->node_index])
          source_count++;
      }
      *lnacceptance += log((double)target_count / source_count);
    }

    /* All moves nodes for current locus have been identified. Apply SPR on gene
       tree. For each moved node, prune off its pureA child and reattach it at
       gtarget
    */
    for (j = 0; j < moved_count[i]; ++j)
    {
      double tnew = moved_nodes[j]->time*taufactor;
      for (k = 1; rway[k]; ++k)
        if (tnew < rway[k]->tau) break;
      snode_t * newpop = rway[k-1];

      /* TODO: We probably don't need to keep the pruned nodes array above, but only check the 'mark' */
      /* TODO: 5.8.2020 REMOVE PRUNED NODES ARRAY */
      gnode_t * node = pruned_nodes[j]->parent;
      assert(node == moved_nodes[j]);
      node->time = tnew;
      gnode_t * pruned = pruned_nodes[j];
      gnode_t * intact = (node->left == pruned) ? node->right : node->left;

      /* possibly many insertions on same branch */
      for (gtarget = gtarget_nodes[j]; gtarget->parent; gtarget = gtarget->parent)
        if (gtarget->parent->time > tnew)
          break;


      /* check whether gene tree topology changes */
      if (gtarget != intact && gtarget != node)
      {
        /* topology changes */

        intact->parent = node->parent;
        if (node->parent)
        {
          if (node->parent->left == node)
            node->parent->left = intact;
          else
            node->parent->right = intact;
        }
        else
        {
          gtree->root = intact;

          /* TODO: This bug has led to serious headaches for many days.
          The root node points to an invalid pmatrix_index (= edge_count) which
          is not associated to a pmatrix. However, when we call SWAP_PMAT_INDEX
          it is changed to 0, resulting to two nodes having the same indices */
          SWAP(node->pmatrix_index,gtree->root->pmatrix_index);
        }

        if (gtarget->parent)
        {
          if (gtarget->parent->left == gtarget)
            gtarget->parent->left = node;
          else
            gtarget->parent->right = node;
        }
        else
        {
          gtree->root = node;

          SWAP(gtarget->pmatrix_index,gtree->root->pmatrix_index);
        }

        node->parent   = gtarget->parent;
        gtarget->parent = node;
        assert(gtarget->parent == node);

        if (node->left == intact)
          node->left = gtarget;
        else
          node->right = gtarget;
      }

      /* scaling ages of moved node PureA clades */
      scaled_count++;           /* this is for the moved node */
      scaled_count += snl_scale_clade(pruned,rway,ytaunew,taufactor,i,thread_index);
        
      gtarget->mark |= FLAG_BRANCH_UPDATE;
      node->mark    |= FLAG_BRANCH_UPDATE;
      intact->mark  |= FLAG_BRANCH_UPDATE;

      /* remove  gene node from list of coalescent events of its old population */
      unlink_event(node, i);
      node->pop->event_count[i]--;

      node->pop->mark[thread_index] |= FLAG_POP_UPDATE;
      if (!opt_est_theta)
        node->pop->event_count_sum--;

      node->pop = newpop;
      node->pop->mark[thread_index] |= FLAG_POP_UPDATE;

      dlist_item_append(node->pop->event[i], node->event);

      node->pop->event_count[i]++;
      if (!opt_est_theta)
        node->pop->event_count_sum++;

      /* update leaf counts for moved node */
      for (; node; node = node->parent)
        node->leaves = node->left->leaves + node->right->leaves;

      /* update leaf counts for intact node parent */
      for (node = intact->parent; node; node = node->parent)
        node->leaves = node->left->leaves + node->right->leaves;
    }

    /* rescale ages and reset populations if the whole tree is a pure-A clade */
    if (gtree->root->mark & SNL_PUREA)
      scaled_count += snl_scale_clade(gtree->root,rway,ytaunew,taufactor,i,thread_index);

    /* Now process square nodes: AB -> B */
    dlist_item_t * item = y->event[i]->head;
    while (item)
    {
      gnode_t * node = (gnode_t *)(item->data);
      item = item->next; 

      if (node->mark & (SNL_MOVED | SNL_PUREA)) continue;

      /* square node */
      unlink_event(node,i);
      node->pop->event_count[i]--;
      if (!opt_est_theta)
        node->pop->event_count_sum--;

      node->pop = b;

      dlist_item_append(node->pop->event[i], node->event);

      node->pop->event_count[i]++;
      if (!opt_est_theta)
        node->pop->event_count_sum++;

      y->mark[thread_index] |= FLAG_POP_UPDATE;
      b->mark[thread_index] |= FLAG_POP_UPDATE;
    }
    
    /* Now process diamond nodes: C -> AC */
    item = c->event[i]->head;
    while (item)
    {
      gnode_t * node = (gnode_t *)(item->data);
      item = item->next; 

      if (node->time <= ytaunew) continue;

      /* diamond node */
      unlink_event(node,i);
      node->pop->event_count[i]--;
      if (!opt_est_theta)
        node->pop->event_count_sum--;

      node->pop = y;

      dlist_item_append(node->pop->event[i], node->event);

      node->pop->event_count[i]++;
      if (!opt_est_theta)
        node->pop->event_count_sum++;

      y->mark[thread_index] |= FLAG_POP_UPDATE;
      c->mark[thread_index] |= FLAG_POP_UPDATE;
    }


    /* we must always recompute for Y */
    y->mark[thread_index] |= FLAG_POP_UPDATE;

    /* Flag populations Y,C and B for re-computing their contribution to the i-th gene tree probability only if:
         (a) they have not been already flagged in a previous step
         (b) there is more than one outgoing lineages (entering its parent population).
    */

    if (!(c->mark[thread_index] & FLAG_POP_UPDATE) && (c->seqin_count[i] - c->event_count[i] > 1))
      c->mark[thread_index] |= FLAG_POP_UPDATE;
      //snode_contrib[snode_contrib_count[i]++] = c;
    if (!(b->mark[thread_index] & FLAG_POP_UPDATE) && (b->seqin_count[i] - b->event_count[i] > 1))
      b->mark[thread_index] |= FLAG_POP_UPDATE;
      //snode_contrib[snode_contrib_count[i]++] = b;

    moved_nodes += gtree->inner_count;
    gtarget_nodes += gtree->inner_count;
    gtarget_list += gtree->tip_count + gtree->inner_count;

    /* TODO: 6.8.2020 THIS HAS NOT BEEN CHECKED YET */
    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      assert(opt_debug_full);
      /* relaxed clock */
      for (j = 0; j < gtree_list[i]->tip_count+gtree_list[i]->inner_count; ++j)
      {
        gnode_t * node = gtree_list[i]->nodes[j];

        /* if root or parent already marked for branch update, skip */
        if (!node->parent || (node->mark & FLAG_BRANCH_UPDATE)) continue;

        if (stree->pptable[node->pop->node_index][b->node_index])
        {
          if (node->time >= y->tau)
          {
            node->mark |= FLAG_BRANCH_UPDATE;
            node->parent->mark |= FLAG_PARTIAL_UPDATE;
          }
        }
        else if (stree->pptable[node->pop->node_index][y->node_index] &&
                 stree->pptable[y->node_index][node->parent->pop->node_index])
        {
          node->mark |= FLAG_BRANCH_UPDATE;
          node->parent->mark |= FLAG_PARTIAL_UPDATE;
        }
      }
    }

    /* fill snode_contrib array */
    for (j=0; j < stree->tip_count + stree->inner_count; ++j)
    {
      snode_t * stmp = stree->nodes[j];
      if (stree->pptable[stmp->node_index][a->node_index])
      {
        /* if a population of clade A was not flagged for recomputing MSC
          density, it means that it has no coalescent events (no nodes were
          rescaled). However, if more than 1 sequences pass it, we must
          recompute its density as its parent tau was changed */
        if (!(stmp->mark[thread_index] & FLAG_POP_UPDATE) &&
            stmp->seqin_count[i] > 1)
        {
          stmp->mark[thread_index] |= FLAG_POP_UPDATE;
        }
      }
      if (stmp->mark[thread_index] & FLAG_POP_UPDATE)
        snode_contrib[snode_contrib_count[i]++] = stmp;
      
      /* reset species tree marks */
      stmp->mark[thread_index] = 0;
    }    
    snode_contrib += stree->tip_count + stree->inner_count;

  } /* end of locus */

  /* update species tree */

  if (!y->parent)
  {
    assert(y == stree->root); assert(movetype == SHRINK);
    stree->root = b;

    /* TODO clock3 test brate swap strategies */
    if (opt_clock == BPP_CLOCK_CORR)
    {
      /* swap strategy */
      double * brate_tmp = y->brate;
      y->brate = b->brate;
      b->brate = brate_tmp;
    }
  }
  else
  {
    /* make b child of y->parent */
    if (y->parent->left == y)
      y->parent->left = b;
    else
      y->parent->right = b;
  }
  b->parent = y->parent;

  /* make y child of c->parent */
  if (!c->parent)
  {
    assert(c == stree->root); assert(movetype == EXPAND);
    stree->root = y;

    /* TODO clock3 test brate swap strategies */
    if (opt_clock == BPP_CLOCK_CORR)
    {
      /* swap stretegy */
      double * brate_tmp = y->brate;
      y->brate = c->brate;
      c->brate = brate_tmp;
    }
  }
  else
  {
    if (c->parent->left == c)
      c->parent->left = y;
    else
      c->parent->right = y;
  }
  y->parent = c->parent;
  c->parent = y;

  /* make c child of y */
  if (y->left == a)
    y->right = c;
  else
    y->left = c;

  /* set new tau */
  y->tau = ytaunew;

  /* scale nodes inside clade A on species tree by taufactor */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    /* note, pptable is invalid at the moment, but does not affect A node
     * ancestry */
    if (stree->pptable[stree->nodes[i]->node_index][a->node_index] &&
        stree->nodes[i]->tau)
    {
      stree->nodes[i]->tau *= taufactor;
      ++scaled_count;   /* counts a but not y */
    }
  }

  assert(!opt_msci);
  for (i = 0; i < stree->locus_count; ++i)
    fill_seqin_counts(stree,NULL,i);

  /* TODO: Check whether reset_gene_leaves_count must operate on the whole gtree_list, or whether
     we can separate the 'reset_hybrid_gene_leaves_count() call inside the function */
  assert(!opt_msci);
  reset_gene_leaves_count(stree,gtree_list);
  stree_reset_pptable(stree);
  stree_label(stree);

  /* TODO: 5.8.2020 At the moment, no consraints */
  assert(!opt_constraint_count);
  init_weights(stree,NULL);

  if (movetype == EXPAND)
  {
    if (downwards)
      *lnacceptance += log(c->parent->weight);
    else
      *lnacceptance += log(c->weight);
  }
  else
  {
    *lnacceptance += log(y->weight);
  }

  *lnacceptance += scaled_count*log(taufactor);

  double newprior = lnprior_species_model(stree);
  *lnacceptance += newprior - oldprior;

  /* prior on tau's (YR2010: Eq. 2) */
  tau0new = stree->root->tau;
  if (fabs(tau0new - tau0) > 1e-20)
  {
    if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
      *lnacceptance += (-opt_tau_alpha - 1 - (ndspecies - 2))*log(tau0new/tau0) -
                       opt_tau_beta*(1/tau0new - 1/tau0);
    else
      *lnacceptance += (opt_tau_alpha - 1 - (ndspecies - 2))*log(tau0new/tau0) -
                       opt_tau_beta*(tau0new - tau0);
  }

  snode_contrib = snode_contrib_space;
  double logpr_notheta = stree->notheta_logpr;
  for (i = 0; i < stree->locus_count; ++i)
  {
    gtree_list[i]->old_logl = gtree_list[i]->logl;

    gnode_t ** bl_list = __gt_nodes;

    if (opt_debug_full)
    {
        k=0;
        for (j = 0; j < gtree_list[i]->tip_count + gtree_list[i]->inner_count; ++j)
        {
          gnode_t * tmp = gtree_list[i]->nodes[j];
          if (tmp->parent)
          {
            tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree_list[i]->edge_count,
                                                 tmp->pmatrix_index);
            bl_list[k++] = tmp;
          }
        }
    
    
        locus_update_matrices(loci[i],gtree_list[i],bl_list,stree,i,k);
    
        gtree_all_partials(gtree_list[i]->root,bl_list,&k);
        for (j = 0; j < k; ++j)
        {
          bl_list[j]->clv_index = SWAP_CLV_INDEX(gtree_list[i]->tip_count,
                                                 bl_list[j]->clv_index);
          if (opt_scaling)
            bl_list[j]->scaler_index = SWAP_SCALER_INDEX(gtree_list[i]->tip_count,
                                                         bl_list[j]->scaler_index);
        }
    
        locus_update_partials(loci[i],bl_list,k);
    
        /* compute log-likelihood */
        assert(!gtree_list[i]->root->parent);
        gtree_list[i]->logl = locus_root_loglikelihood(loci[i],
                                                       gtree_list[i]->root,
                                                       loci[i]->param_indices,
                                                       NULL);
    }
    else
    {
        k = 0;
        for (j = 0; j < gtree_list[i]->tip_count + gtree_list[i]->inner_count; ++j)
        {
          gnode_t * tmp = gtree_list[i]->nodes[j];
    
          if (!tmp->parent) continue;
    
          if (tmp->mark & FLAG_BRANCH_UPDATE)
          {
            tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree_list[i]->edge_count,
                                                 tmp->pmatrix_index);
            bl_list[k++] = tmp;
          }
        }
    
        if (k)
        {
          locus_update_matrices(loci[i], gtree_list[i], bl_list, stree, i, k);
    
          /* retrieve all nodes whose partials must be updated */
          unsigned int partials_count;

          assert(!gtree_list[i]->root->parent);
          gnode_t ** partials = gtree_list[i]->travbuffer;
          gtree_return_partials(gtree_list[i]->root,
                                gtree_list[i]->travbuffer,
                                &partials_count);
    
          /* point to the double-buffered partials space */
          for (j = 0; j < partials_count; ++j)
          {
            partials[j]->clv_index = SWAP_CLV_INDEX(gtree_list[i]->tip_count,
                                                    partials[j]->clv_index);
            if (opt_scaling)
              partials[j]->scaler_index = SWAP_SCALER_INDEX(gtree_list[i]->tip_count,
                                                            partials[j]->scaler_index);
          }
    
          /* update conditional probabilities (partials) of affected nodes */
          locus_update_partials(loci[i], partials, partials_count);
    
          /* evaluate log-likelihood */
          assert(!gtree_list[i]->root->parent);
          gtree_list[i]->logl = locus_root_loglikelihood(loci[i],
                                                         gtree_list[i]->root,
                                                         loci[i]->param_indices,
                                                         NULL);
        }
    }

    if (opt_est_theta)
      gtree_list[i]->old_logpr = gtree_list[i]->logpr;

    if (opt_debug_full)
    {
        /* This recomputes the gene tree probabilities from scratch. It can be used to verify that
           the code below, which only computes the gene tree probability for the changed components,
           is correct. */
    
        if (opt_est_theta)
          gtree_list[i]->old_logpr = gtree_list[i]->logpr;
    
        
        double logpr;
        if (opt_migration)
          logpr = gtree_logprob_mig(stree, gtree_list[i], loci[i]->heredity[0], i, thread_index);
        else
          logpr = gtree_logprob(stree, loci[i]->heredity[0], i, thread_index);

    
    
        if (opt_est_theta)
          gtree_list[i]->logpr = logpr;
        else
        {
          if (i == opt_locus_count - 1) 
            logpr += stree->notheta_hfactor+stree->notheta_sfactor;
          logpr_notheta = logpr;
        }
    }
    else
    {
         /* locate additional populations that need to be updated */
    
         /* find and mark those populations whose number of incoming lineages has
            changed due to the reset_gene_leaves_count() call, but were previously
            not marked for log-probability contribution update */
        for (j = 0; j < snode_contrib_count[i]; ++j)
          snode_contrib[j]->mark[thread_index] |= FLAG_POP_UPDATE;
        for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
        {
          snode_t * snode = stree->nodes[j];
    
          if (!(snode->mark[thread_index] & FLAG_POP_UPDATE) &&
              (snode->seqin_count[i] != original_stree->nodes[j]->seqin_count[i]))
            snode_contrib[snode_contrib_count[i]++] = snode;
    
    
        }
    
        /* now update the log-probability contributions for the affected, marked
           populations */
        for (j = 0; j < snode_contrib_count[i]; ++j)
        {
          if (opt_est_theta)
            gtree_list[i]->logpr -= snode_contrib[j]->logpr_contrib[i];
          else
            logpr_notheta -= snode_contrib[j]->notheta_logpr_contrib;
    
          double xtmp;
          
          if (opt_migration)
            xtmp = gtree_update_logprob_contrib_mig(snode_contrib[j],
                                                    stree,
                                                    gtree_list[i],
                                                    loci[i]->heredity[0],
                                                    i,
                                                    thread_index);
          else
            xtmp = gtree_update_logprob_contrib(snode_contrib[j],
                                                loci[i]->heredity[0],
                                                i,
                                                thread_index);

    
          if (opt_est_theta)
            gtree_list[i]->logpr += snode_contrib[j]->logpr_contrib[i];
          else
            logpr_notheta += xtmp;
        }
    }

    /* 
       TODO: Several improvements can be made here
       1. Call this function only for the gene trees that are modified
       2. Only re-compute the prior for the affected gene tree nodes/edges
    */
    if (opt_clock == BPP_CLOCK_CORR)
    {
      double new_prior_rates = lnprior_rates(gtree_list[i],stree,i);
      *lnacceptance += new_prior_rates - gtree_list[i]->lnprior_rates;
      gtree_list[i]->lnprior_rates = new_prior_rates;
    }

    /* reset markings on affected populations */
    for (j = 0; j < snode_contrib_count[i]; ++j)
      snode_contrib[j]->mark[thread_index] = 0;

    for (j = 0; j < gtree_list[i]->tip_count + gtree_list[i]->inner_count; ++j)
      gtree_list[i]->nodes[j]->mark = 0;

    snode_contrib += stree->tip_count + stree->inner_count;

    if (opt_est_theta)
      *lnacceptance += gtree_list[i]->logpr - gtree_list[i]->old_logpr +
                       gtree_list[i]->logl - gtree_list[i]->old_logl;
    else
      *lnacceptance += gtree_list[i]->logl - gtree_list[i]->old_logl;
  }
  #if 0
  debug_consistency(stree, gtree_list);
  #endif

  double r = legacy_rndu(thread_index);

  if (opt_debug_snl)
    debug_snl_stage2(stree,gtree_list,logpr_notheta,r,*lnacceptance);

  if (!opt_est_theta)
  {
    *lnacceptance += logpr_notheta - stree->notheta_logpr;
    stree->notheta_logpr = logpr_notheta;
  }

  return (*lnacceptance >= -1e-10 || r < exp(*lnacceptance));
}

long stree_propose_stree_snl(stree_t ** streeptr,
                             gtree_t *** gtree_list_ptr,
                             stree_t ** scloneptr,
                             gtree_t *** gclonesptr,
                             locus_t ** loci)
{
  int downwards;
  unsigned int i;
  long movetype;
  long thread_index = 0;
  double r;
  double sum = 0;
  double lnacceptance = 0;
  double delta;
  double tau_new = 0;
  double tau_factor = 0;
  snode_t * y;
  snode_t * a;
  snode_t * b;
  snode_t * c = NULL;
  snode_t * x = NULL;
  snode_t * stmp;
  snode_t * target, * nextnode, * prevnode;
  snode_t ** rway;
  snode_t * lca = NULL;

  stree_t * original_stree = *streeptr;
  gtree_t ** original_gtree_list = *gtree_list_ptr;

  stree_t * stree = *scloneptr;
  gtree_t ** gtree_list = *gclonesptr;

  stree_clone(original_stree, stree);

  for (i = 0; i < stree->locus_count; ++i)
    gtree_clone(original_gtree_list[i], gtree_list[i], stree);
  events_clone(original_stree, stree, gtree_list);

  /* calculate the weight of each branch as the reciprocal of the square root
   * of its length */
  /* TODO: 31.7.2020 At the moment, no consraints */
  if (opt_constraint_count)
    fatal("\nERROR: Constraints are not yet implemented for the SNL move.\n"
          "Please edit the control file and replace the following line:\n\n"
          "  speciestree = 1\n\n"
          "to:\n\n"
          "  speciestree = 1 0\n\n"
          "to disable the SNL move and use constraints.");

  init_weights(stree, NULL);

  lnacceptance = 0;
  /* select move type */
  if (legacy_rndu(thread_index) < opt_prob_snl_shrink)
    movetype = SHRINK;
  else
    movetype = EXPAND;

  /* randomly select a branch according to weights */
  sum = 0;
  r = legacy_rndu(thread_index);
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count - 1; ++i)
  {
    sum += stree->nodes[i]->weight;
    if (r < sum) break;
  }
  assert(stree->nodes[i] != stree->root);
  assert(stree->nodes[i]->weight);
  lnacceptance -= log(stree->nodes[i]->weight);

  /* y-a for expand and y-c for shrink. Sample delta and set initial target for regrafting y-a */
  y = NULL;  a = NULL;  b = NULL;  c = NULL;  x = NULL;
  if (movetype == EXPAND)  /* (++ and +-) */
  {
    y = stree->nodes[i];
    x = y->parent;
    if ((int)(2 * legacy_rndu(thread_index)) == 0) { a = y->left;  b = y->right; }
    else { a = y->right; b = y->left; }
    delta = x->tau * (1 - pow(legacy_rndu(thread_index), 1. / opt_snl_lambda_expand));
    target = x;
  }
  else    /* (--) */
  {
    c = stree->nodes[i];
    y = c->parent;
    if (y->left == c) { a = y->right; b = y->left; }
    else { a = y->left;  b = y->right; }
    delta = c->tau * (1 - pow(legacy_rndu(thread_index), 1. / opt_snl_lambda_shrink));
    target = c;
  }
  /* only required for the first time an expand move takes a downward path, so
     that we don't go back the same path we came from */
  prevnode = y;
  downwards = 0;  /* if downwards once, then always downwards */
  if (movetype == SHRINK) downwards = 1;

  /* Find target branch for regrafting. Code works for both expand and shrink moves */
  while (1)
  {
    /* chose whether to take upward or downward path */
    if (!downwards && legacy_rndu(thread_index) < 0.5)  /* upwards */
    {
      if (target->parent)
      {
        double dist = target->parent->tau - target->tau;
        if (dist > delta)
        {
          tau_new = target->tau + delta;
          break;
        }
        else
        {
          prevnode = target;
          target = target->parent;
          delta -= dist;
        }
      }
      else
      {
        /* we reached the root node, hence we just add delta to root's tau */
        tau_new = target->tau + delta;
        break;
      }
    }
    else                  /* downwards */
    {
      if (!downwards)    /* (+-) */
      {
        /* first time a downwards move is selected. Ensure we do not visit same
           path again */
        nextnode = (target->left == prevnode) ? target->right : target->left;
        /* all future steps will be downwards */
        downwards = 1;
        /* store the common ancestor of x and target which we will need to
           compute the delta for the reverse move */
        lca = target;
      }
      else
      {
        /* randomly select which of the two */
        nextnode = legacy_rndu(thread_index) < 0.5 ? target->left : target->right;
      }
      double dist = target->tau - nextnode->tau;
      if (dist > delta)
      {
        tau_new = target->tau - delta;
        target = nextnode;
        break;
      }
      else
      {
        /* prevnode = target; */  /* only for consistency, but not needed anymore */
        target = nextnode;
        delta -= dist;
      }
    }
  }  /* while (1), loop for finding target branch */

  if ((movetype == EXPAND && !downwards) || movetype == SHRINK)
  {
    /* shrink or pure expand (no downward step) (-- or ++) */
    if (movetype == EXPAND)
    {
      /* pure expand (++): delta = tau_new - x->tau, delta* = target->tau - y->tau */
      assert(tau_new - x->tau < x->tau);
      if (target->tau - y->tau >= target->tau)
        return 2;

      lnacceptance += logpdf_power(target->tau - y->tau, target->tau, opt_snl_lambda_shrink);
      lnacceptance -= log(0.5);
      lnacceptance -= logpdf_power(tau_new - x->tau, x->tau, opt_snl_lambda_expand);
      lnacceptance += log(opt_prob_snl_shrink / (1 - opt_prob_snl_shrink));
    }
    else
    {
      /* shrink only (--): delta = c->tau - tau_new, delta* = y->tau - target->parent->tau */
      assert(c->tau - tau_new < c->tau);
      if (y->tau - target->parent->tau >= target->parent->tau)
        return 2;

      lnacceptance += logpdf_power(y->tau - target->parent->tau, target->parent->tau, opt_snl_lambda_expand);
      lnacceptance += log(0.5);
      lnacceptance -= logpdf_power(c->tau - tau_new, c->tau, opt_snl_lambda_shrink);
      lnacceptance += log((1 - opt_prob_snl_shrink) / opt_prob_snl_shrink);
    }
  }
  else
  {
    /* expand with both upward and downward steps (+-) */
    assert(movetype == EXPAND && downwards);
    assert(lca);
    /* compute distance */
    double dist = lca->tau - x->tau + lca->tau - tau_new;
    double dist_rev = lca->tau - y->tau + lca->tau - target->parent->tau;
    if (dist_rev >= target->parent->tau)
      return 2;

    lnacceptance += logpdf_power(dist_rev, target->parent->tau, opt_snl_lambda_expand);
    assert(dist < x->tau);  /*** Ziheng ??? ***/
    lnacceptance -= logpdf_power(dist, x->tau, opt_snl_lambda_expand);
  }

  tau_factor = tau_new / y->tau;
  assert(tau_factor);

  /* compute rway */
  rway = (snode_t **)xcalloc((size_t)(stree->inner_count+1),sizeof(snode_t *));
  rway[0] = y;  /* first in the list is Y*  */
  for (i = 1, stmp = target->parent; stmp; stmp = stmp->parent)
  {
    if (stmp == y) continue;
    rway[i++] = stmp;
  }
  rway[i] = NULL;

  /* statistics for shrink/expand */
  if (movetype == SHRINK)
    ++opt_debug_shrink_count;
  else
  {
    if (downwards)
      ++opt_debug_expshr_count;
    else
      ++opt_debug_expand_count;
  }

  if (opt_debug_snl)
    debug_snl_stage1(stree,
                     gtree_list,
                     y,
                     target,
                     a,
                     movetype,
                     downwards,
                     y->tau,
                     tau_new);

  long rc = snl_expand_and_shrink(stree,
                                  original_stree,
                                  gtree_list,
                                  loci,
                                  movetype,
                                  downwards,
                                  tau_new,
                                  tau_factor,
                                  &lnacceptance,
                                  a,
                                  b,
                                  target,
                                  y,
                                  rway);

  free(rway);
  return rc;
}

/* checks whether migration exists/valid (forward in time) */
long migration_valid(stree_t * stree, snode_t * from, snode_t * to)
{
  if (!opt_mig_bitmatrix[from->node_index][to->node_index])
    return 0;
  if (!from->parent || !to->parent)
    return 0;

  if (from->tau <= to->tau && from->parent->tau > to->tau)
    return 1;

  if (from->tau >= to->tau && from->tau <= to->parent->tau)
    return 1;

  return 0;
}

double migrate_gibbs(stree_t * stree,
                     gtree_t ** gtree,
                     locus_t ** locus,
                     unsigned int si,
                     unsigned int ti)
{
  long i,j,k,n;
  long asj = 0;
  long msa_index;
  double a1,b1;
  double bsj = 0;
  double heredity;
  double tstart,tend;
  size_t alloc_required;
  dlist_item_t * event;
  dlist_item_t * li;
  migbuffer_t * migbuffer;

  snode_t * src = stree->nodes[si];
  snode_t * tgt = stree->nodes[ti];

  const static long thread_index = 0;

  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  /* indicator: 1 if migration from si to ti is possible at time segment k */
  if (!ts_indicator)
    ts_indicator = (long *)xcalloc((size_t)total_nodes, sizeof(long));


  for (msa_index = 0; msa_index < opt_locus_count; ++msa_index)
  {
    alloc_required = tgt->migevent_count[msa_index] +
                     tgt->event_count[msa_index] +
                     stree->inner_count+1;
    migbuffer_check_and_realloc(thread_index,alloc_required);
    migbuffer = global_migbuffer_r[thread_index];

    heredity = locus[msa_index]->heredity[0];

    /* calculate Asj */
    asj += gtree[msa_index]->migcount[si][ti];

    /* add taus and coalescence times in sortbuffer */
    migbuffer[0].time = tgt->tau;
    migbuffer[0].type = EVENT_TAU;
    j = 1;
    for (event = tgt->event[msa_index]->head; event; event = event->next)
    {
      gnode_t* gnode = (gnode_t*)(event->data);
      migbuffer[j].time   = gnode->time;
      migbuffer[j++].type = EVENT_COAL;
    }

    for (li = tgt->mig_source[msa_index]->head; li; li = li->next)
    {
      migevent_t * me     = (migevent_t *)(li->data);
      migbuffer[j].time   = me->time;
      migbuffer[j++].type = EVENT_MIG_SOURCE;
    }
    for (li = tgt->mig_target[msa_index]->head; li; li = li->next)
    {
      migevent_t * me     = (migevent_t *)(li->data);
      migbuffer[j].time   = me->time;
      migbuffer[j++].type = EVENT_MIG_TARGET;
    }

    /* add splitting of populations */
    for (k = 0; k < tgt->mb_count; ++k)
    {
      migbuffer[j++] = tgt->migbuffer[k];
    }

    if (tgt->parent && !tgt->mb_count)
      fatal("\nError when processing node %s", tgt->label);
    /* TODO: Probably split the following qsort case into two:
       in case tgt->parent then sort j-2 elements, otherwise
       j-1 elements.
    */

    /* if there was at least one coalescent event, sort */
    if (j > 1)
      qsort(migbuffer + 1, j - 1, sizeof(migbuffer_t), cb_migbuf_asctime);

    
    /* fill indicator variable */
    for (k = 0, tstart = tgt->tau; k < tgt->mb_count; ++k)
    {
      tend = tgt->migbuffer[k].time;
      ts_indicator[k] = (src->tau <= tstart && src->parent->tau >= tend);
      tstart = tend;
    }

    long epoch = 0;
    assert(!tgt->parent || tgt->mb_count);
    long ind = ts_indicator[epoch];

    /* calculate Bsj */
    for (k = 1, n = tgt->seqin_count[msa_index]; k < j; ++k)
    {
      double t = migbuffer[k].time - migbuffer[k-1].time;
      if (n > 0 && tgt->parent)
        bsj += 4*n*ind*t;

      if (migbuffer[k].type == EVENT_COAL || migbuffer[k].type == EVENT_MIG_SOURCE)
        --n;
      else if (migbuffer[k].type == EVENT_MIG_TARGET)
        ++n;
      else if (migbuffer[k].type == EVENT_TAU && epoch < tgt->mb_count-1)
        ind = ts_indicator[++epoch];
    }
    bsj /= heredity;
  }
  bsj /= stree->nodes[ti]->theta;

  /* we have the new distribution */
  long mindex = opt_migration_matrix[si][ti];
  assert(mindex >= 0);
  migspec_t * spec = opt_mig_specs+mindex;

  a1 = spec->alpha + asj;
  b1 = spec->beta  + bsj;

  if (opt_mig_specs[mindex].am)
    assert(0);  /* TODO: Go through all loci and propose Mi? */
  else
    opt_mig_specs[mindex].M = legacy_rndgamma(thread_index,a1) / b1;

  stree_update_mig_subpops(stree,thread_index);

  /* update gene tree density */
  /* TODO: Update only necessary populations to improve computational speed */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree[i]->logpr = gtree_logprob_mig(stree,
                                        gtree[i],
                                        locus[i]->heredity[0],
                                        i,
                                        thread_index);
  }

  return 1;
}

static double prop_migrates_gibbs(stree_t * stree, gtree_t ** gtree, locus_t ** locus)
{
  long i,j;
  long accepted = 0;
  long total = 0;

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
    {
      if (!migration_valid(stree, stree->nodes[i], stree->nodes[j])) continue;

      ++total;

      accepted += migrate_gibbs(stree,gtree,locus,i,j);
    }
  }


  return accepted / (double)total;

}

static long prop_migrates_mbar_slide(migspec_t * spec)
{
  long i;
  long accepted = 0;
  double c,lnc;
  double mbar_old,mbar_new;
  double lnacceptance;
  
  assert(opt_est_theta);

  const static long thread_index_zero = 0;

  mbar_old = spec->M;

  lnc = opt_finetune_migrates * legacy_rnd_symmetrical(thread_index_zero);
  c = exp(lnc);

  mbar_new = mbar_old * c;
  lnacceptance = lnc + lnc*(spec->alpha-1) - (mbar_new-mbar_old)*spec->beta;

  double a = spec->am;
  double bnew = a / mbar_new;
  double bold = a / mbar_old;

  lnacceptance += opt_locus_count*a*log(bnew/bold);
  for (i = 0; i < opt_locus_count; ++i)
    lnacceptance -= (bnew - bold)*spec->Mi[i];

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
  {
    /* accepted */
    spec->M = mbar_new;
    accepted = 1;
  }

  return accepted;
}

static double prop_migrates_slide(stree_t * stree, gtree_t ** gtree, locus_t ** locus)
{
  long i,j,k;
  long accepted = 0;
  long total = 0;
  long mindex;
  double c,lnc;
  double rate_old,rate_new;
  double logpr_diff;
  double lnacceptance;
  migspec_t * spec;
  
  assert(opt_est_theta);

  const static long thread_index_zero = 0;

  /* TODO: Change the nested for loops with one loop over opt_mig_spec */
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
    {
      if (!migration_valid(stree, stree->nodes[i], stree->nodes[j])) continue;

      ++total;

      mindex = opt_migration_matrix[i][j];
      spec = opt_mig_specs+mindex;

      /* in case of Mbar */
      if (spec->am)
      {
        accepted += prop_migrates_mbar_slide(spec);
        continue;
      }
      
      /* in case of M */

      rate_old = spec->M;

      #if 1
      lnc = opt_finetune_migrates * legacy_rnd_symmetrical(thread_index_zero);
      c = exp(lnc);

      rate_new = rate_old * c;

      lnacceptance = lnc + lnc*(spec->alpha-1) - (rate_new-rate_old)*spec->beta;
      #else
      double lograte_old, lograte_new;
      double minv = -99, maxv = 99;

      lograte_old = log(rate_old);

      lograte_new = lograte_old + opt_finetune_migrates * legacy_rnd_symmetrical(thread_index_zero);
      lograte_new = reflect(lograte_new, minv, maxv, thread_index_zero);

      rate_new = exp(lograte_new);

      lnacceptance = lograte_new - lograte_old;
      lnacceptance += (spec->alpha-1)*log(rate_new/rate_old) -
                      spec->beta*(rate_new - rate_old);
      #endif

      spec->M = rate_new;

      /* TODO: improve the following */
      #if 0
      stree_update_mig_subpops(stree,thread_index_zero);
      #else
      stree_update_mig_subpops_single(stree, stree->nodes[j], stree->nodes[i], rate_old);
      #endif

      for (k = 0; k < opt_locus_count; ++k)
      {
        logpr_diff   = -stree->nodes[j]->logpr_contrib[k];
        logpr_diff  += gtree_update_logprob_contrib_mig(stree->nodes[j],
                                                        stree,
                                                        gtree[k],
                                                        locus[k]->heredity[0],
                                                        k,
                                                        thread_index_zero);
                                
        lnacceptance += logpr_diff;

        /* update gene tree density */
        gtree[k]->old_logpr = gtree[k]->logpr;
        gtree[k]->logpr += logpr_diff;
      }

      if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
      {
        /* accepted */
        accepted++;
      }
      else
      {
        spec->M = rate_old;

        /* TODO: improve the following */
        #if 0
        stree_update_mig_subpops(stree,thread_index_zero);
        #else
        stree_update_mig_subpops_single(stree, stree->nodes[j], stree->nodes[i], rate_new);
        #endif

        for (k = 0; k < opt_locus_count; ++k)
        {
          gtree[k]->logpr = gtree[k]->old_logpr;
          stree->nodes[j]->logpr_contrib[k] = stree->nodes[j]->old_logpr_contrib[k];
        }
      }
    }
  }
  if (!total) return 0;
  return accepted / (double)total;
}

#if 1
static double prop_mig_vrates_slide(stree_t * stree, gtree_t ** gtree, locus_t ** locus)
{
  long i,j,k;
  long accepted = 0;
  long total = 0;
  double alpha,beta;
  double c,lnc;
  double rate_old,rate_new;
  double logpr_diff;
  double lnacceptance;
  migspec_t * spec;
  
  assert(opt_est_theta);

  const static long thread_index_zero = 0;

  /* TODO: Change the nested for loops with one loop over opt_mig_spec */
  for (i = 0; i < opt_migration_count; ++i)
  {
    spec = opt_mig_specs+i;

    if (!spec->am || !migration_valid(stree, stree->nodes[spec->si], stree->nodes[spec->ti]))
      continue;

    total += opt_locus_count;;

    for (j = 0; j < opt_locus_count; ++j)
    {
      rate_old = spec->Mi[j];
      alpha = spec->am;
      beta = alpha / spec->M;

      #if 1
      lnc = opt_finetune_mig_Mi * legacy_rnd_symmetrical(thread_index_zero);
      c = exp(lnc);

      rate_new = rate_old * c;

      lnacceptance = lnc + lnc*(alpha-1) - (rate_new-rate_old)*beta;
      #else
      double lograte_old, lograte_new;
      double minv = -99, maxv = 99;

      lograte_old = log(rate_old);

      lograte_new = lograte_old + opt_finetune_migrates * legacy_rnd_symmetrical(thread_index_zero);
      lograte_new = reflect(lograte_new, minv, maxv, thread_index_zero);

      rate_new = exp(lograte_new);

      lnacceptance = lograte_new - lograte_old;
      lnacceptance += (alpha-1)*log(rate_new/rate_old)-beta*(rate_new-rate_old);
      #endif

      spec->Mi[j] = rate_new;

      #if 0
      /* TODO: improve the following. This is very expensive, we need to modify it to
         update only those node migbuffer elements affected by Mi[j] */
      stree_update_mig_subpops(stree,thread_index_zero);
      #else
      stree_update_mig_subpops_single_vrates(stree, stree->nodes[spec->ti], stree->nodes[spec->si], j, rate_old);
      #endif

      /* TODO: We probably change only one locus */
      for (k = 0; k < opt_locus_count; ++k)
      {
        logpr_diff   = -stree->nodes[spec->ti]->logpr_contrib[k];
        logpr_diff  += gtree_update_logprob_contrib_mig(stree->nodes[spec->ti],
                                                        stree,
                                                        gtree[k],
                                                        locus[k]->heredity[0],
                                                        k,
                                                        thread_index_zero);
                                
        lnacceptance += logpr_diff;

        /* update gene tree density */
        gtree[k]->old_logpr = gtree[k]->logpr;
        gtree[k]->logpr += logpr_diff;
      }

      if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
      {
        /* accepted */
        accepted++;
      }
      else
      {
        spec->Mi[j] = rate_old;

        #if 0
        /* TODO: improve the following, see previous TODO */
        stree_update_mig_subpops(stree,thread_index_zero);
        #else
        stree_update_mig_subpops_single_vrates(stree, stree->nodes[spec->ti], stree->nodes[spec->si], j, rate_new);
        #endif

        /* TODO: Same here, change only one locus */
        for (k = 0; k < opt_locus_count; ++k)
        {
          gtree[k]->logpr = gtree[k]->old_logpr;
          stree->nodes[spec->ti]->logpr_contrib[k] = stree->nodes[spec->ti]->old_logpr_contrib[k];
        }
      }
    }
  }
  if (!total) return 0;
  return accepted / (double)total;
}
double prop_mig_vrates(stree_t * stree, gtree_t ** gtree, locus_t ** locus)
{
  return prop_mig_vrates_slide(stree,gtree,locus);
}
#endif

double prop_migrates(stree_t * stree, gtree_t ** gtree, locus_t ** locus)
{
  if (opt_mrate_move == BPP_MRATE_GIBBS)
    return prop_migrates_gibbs(stree,gtree,locus);

  return prop_migrates_slide(stree,gtree,locus);
}

#if 0
double prop_migrates_mbar(stree_t * stree, gtree_t ** gtree)
{
  long i,j;
  long accepted = 0;
  long total = 0;
  double c,lnc;
  double mbar_old,mbar_new;
  double lnacceptance;
  migspec_t * spec;
  
  assert(opt_est_theta);

  const static long thread_index_zero = 0;

  for (i = 0; i < opt_migration; ++i)
  {
    spec = opt_mig_specs+i;

    if (!migration_valid(stree, stree->nodes[spec->si], stree->nodes[spec->ti]))
      continue;

    ++total;

    mbar_old = spec->M;

    lnc = opt_finetune_migrates * legacy_rnd_symmetrical(thread_index_zero);
    c = exp(lnc);

    mbar_new = mbar_old * c;
    lnacceptance = lnc + lnc*(spec->alpha-1) - (mbar_new-mbar_old)*spec->beta;

    double a = spec->am;
    double bnew = a / mbar_new;
    double bold = a / mbar_old;

    lnacceptance += opt_locus_count*a*log(bnew/bold);
    for (j = 0; j < opt_locus_count; ++j)
      lnacceptance -= (bnew - bold)*spec->Mi[j];

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
    {
      /* accepted */
      spec->M = mbar_new;
      accepted++;
    }
  }
  if (!total) return 0;
  return accepted / (double)total;
}
#endif

#define log_pdfgamma(x, a, b)  ( (a)*log(b) - lgamma(a) + ((a)-1)*log(x) - (b)*(x) )

#if 0
static long migspec_flip(stree_t * stree, snode_t * x, snode_t * y)
{
  /* flips M_{x->y} to M_{y->x} 
     
     XXX: Assumes that rate M_{y->x} does not exist */

  const long thread_index_zero = 0;
  migspec_t * spec;

  /* find index for rate M_{x->y} */
  assert(opt_migration_count);
  for (i = 0; i < opt_migration_count; ++i)
  {
    spec = opt_mig_specs+i;

    if (spec->si == x->node_index && spec->ti == y->node_index)
      break;
  }
  if (i == opt_migration_count)
    fatal("Cannot find migration");

  /* flip */
  SWAP(opt_migration_matrix[x->node_index][y->node_index],
       opt_migration_matrix[y->node_index][x->node_index]);
  SWAP(opt_mig_bitmatrix[x->node_index][y->node_index],
       opt_mig_bitmatrix[y->node_index][x->node_index]);

  stree_update_mig_subpops(stree,thread_index_zero);

  return i;
}
#endif

static void migspec_append(stree_t * stree, snode_t * x, snode_t * y)
{
  /* add M_{x->y} to the list */
  long thread_index = 0;
  const long thread_index_zero = 0;

  migspec_t * spec = opt_mig_specs+opt_migration_count;

  spec->Mi = NULL;
  spec->si = x->node_index;
  spec->ti = y->node_index;

  /* reset the rest */
  spec->alpha    = opt_mig_alpha;
  spec->beta     = opt_mig_beta;
  spec->pseudo_a = opt_pseudo_alpha;
  spec->pseudo_b = opt_pseudo_beta;

  /*** Ziheng $$$ ***/
#if 1
  spec->M        = legacy_rndgamma(thread_index,dbg_prop_a) / dbg_prop_b;
  #else
  spec->M        = spec->alpha / spec->beta;
  #endif

  assert(opt_mig_bitmatrix[spec->si][spec->ti] == 0);
  if (opt_migration_matrix[spec->si][spec->ti] != -1)
  {
    printf("[ERROR] adding migration M_%s->%s\n", stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label);
  }
  assert(opt_migration_matrix[spec->si][spec->ti] == -1);

  opt_mig_bitmatrix[spec->si][spec->ti] = 1;
  opt_migration_matrix[spec->si][spec->ti] = opt_migration_count++;

  //stree_update_mig_subpops_single(stree,y,x,0);
  stree_update_mig_subpops(stree,thread_index_zero);
}

static long migspec_remove(stree_t * stree, snode_t * x, snode_t * y)
{
  /* remove M_{x->y} from the list of migration events by:
     1. move it to the end of the list
     2. decrease opt_migration_count
     3. return it's original position in the list */
  long i;
  long pos = 0;
  const long thread_index_zero = 0;

  migspec_t * spec = NULL;

  assert(opt_migration_count);
  for (i = 0; i < opt_migration_count; ++i)
  {
    spec = opt_mig_specs+i;

    if (spec->si == x->node_index && spec->ti == y->node_index)
      break;
  }
  if (i == opt_migration_count)
    fatal("Cannot find migration");

  assert(i == opt_migration_matrix[spec->si][spec->ti]);
  /* update migbuffers  -- this is a bit of a hack */
  //double oldM = spec->M;
  //spec->M = 0;
  //stree_update_mig_subpops_single(stree,y,x,oldM);
  //spec->M = oldM;
  #if(DBG_TF)
  printf("     Setting M_%s->%s (%ld->%ld) (%ld) to -1\n", stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label, spec->si, spec->ti, opt_migration_matrix[spec->si][spec->ti]);
  #endif

  pos = i;
  if (i != opt_migration_count-1)
  {
    #if 0
    printf("spec->si: %d\n", spec->si);
    printf("spec->ti: %d\n", spec->ti);
    #endif
    /* swap */
    /* TODO: Save old record */
    SWAP(opt_mig_specs[i],opt_mig_specs[opt_migration_count-1]);

    /* reset pointer */
    spec = opt_mig_specs+opt_migration_count-1;

    SWAP(opt_migration_matrix[x->node_index][y->node_index],
         opt_migration_matrix[opt_mig_specs[i].si][opt_mig_specs[i].ti]);
    SWAP(opt_mig_bitmatrix[x->node_index][y->node_index],
         opt_mig_bitmatrix[opt_mig_specs[i].si][opt_mig_specs[i].ti]);
  }
  opt_migration_matrix[spec->si][spec->ti] = -1;

  --opt_migration_count;

  opt_mig_bitmatrix[spec->si][spec->ti] = 0;

  stree_update_mig_subpops(stree,thread_index_zero);

  #if 0
  printf("[DEBUG] Removed migration rate %s -> %s\n", x->label, y->label);
  #endif

  return pos;
}

#if 0
/* old version that i decided to abandon */
long dissolve_incoming_lineages(stree_t * stree,
                                gtree_t * gtree,
                                snode_t * x,
                                gnode_t ** output,
                                gnode_t ** deleted)
{
  long i,j,k;
  long simcount = 0;
  long delcount = 0;
  gnode_t * gnode;
  snode_t * curpop;


  #if 0
  /* TODO: check about pptable, is it computed and updated according to migration events? I assume here that yes */
  /* Note: no, the above does not hold. I check the pptable later and adjust it to account for migrations */
  if (!stree->pptable[gnode->pop->node_index][x->node_index])
    continue;
  #endif

  /* we have two cases when tracing the lineage backwards in time:
     1. Lineage enters population x from a descendant population y (in the path curpop<->x)
     2. Lineage migrates to other populations before migrating back to x

     For case (1) migrations may occur on the path curpop<->y before the
     lineage enters x.

  */

  /* deleted */
  deleted = (gnode_t **)xmalloc((size_t)gtree->inner_count * sizeof(gnode_t *));

  /* 1. mark lineages to dissolve */
  for (i = 0, k = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode = gtree->nodes[i];

    if (gnode->time > x->parent->tau) continue;

    curpop = gnode->pop;
    for (j = 0; j < gnode->mi->count && curpop != x; ++j)
    {
      if (gnode->mi->me[j].time > x->tau &&
          stree->pptable[curpop->node_index][x->node_index]) break;
      curpop = gnode->mi->me[j].target;
    }
    if (!stree->pptable[curpop->node_index][x->node_index]) continue;

    /* mark with red */
    gnode->parent->mark++;

    output[k++] = gnode;
  }

  /* 2. remove nodes with mark==2, update and repeat until no node is removed */
  while (1)
  {
    for (simcount = 0, i = 0; i < k; ++i)
    {
      if (output[i]->mark != 2)
      {
        /* keep */

        assert(output[i]->mark < 2);
        if (i != simcount)
          output[simcount] = output[i];
        simcount++;
      }
      else
      {
        /* delete */

        deleted[delcount++] = output[i]

        /* mark lineage to parent */
        output[i]->parent->mark++;
        assert(output[i]->parent->mark <= 2);
      }
    }
    if (simcount == k) break;
    k = simcount;
  }
  simcount = k;

  /* We have a list of nodes for deletion. To avoid corrupting the data
     structures, we must sort the nodes in such a way that any parent nodes to
     be deleted are processed before their respective child nodes. In other
     words, if a node x has a parent node y that is also scheduled for deletion,
     y must be deleted before x */
  for (i = 0; i < delcount; ++i)
    assert(output[i]->mark == 2);
  for (i = 0; i < delcount; ++i)
  {
    for (j = i; j < delcount; ++j)
    {
      gnode_t * p = deleted[j];

      if (!p->parent || p->parent->mark < 2)
      {
        /* swap positions */
        gnode_t * tmp = deleted[i];
        deleted[i] = p;
        deleted[j] = tmp;

        /* reset mark */
        p->mark = 0;

        break;
      }
    }
    /* if no swapping occured, break */
    if (j == delcount) break;
  }

  /* 3. clean marks */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    gtree->nodes[i]->mark = 0;


  /* disconnect nodes */
  for (j = 0; j < delcount; ++j)
  {
    gnode_t * p = deleted[j];
    gnode_t * father = p->parent;

    unlink_event(p, gtree->msa_index);
    p->pop->event_count[gtree->msa_index]--;
    if (!opt_est_theta)
      p->pop->event_count_sum--;

    /* decrease the number of incoming lineages to all populations in the path from
    p population (excluding) to the father population (including) */
    curpop = p->pop;
    if (p->mi && p->mi->count)
    {
      for (i = 0; i < p->mi->count; ++i)
      {
        p->mi->me[i].source->migevent_count[msa_index]--;
        p->mi->me[i].target->migevent_count[msa_index]--;

        /* decrease # of migration from s to t (forward in time) */
        long s = p->mi->me[i].target->node_index;
        long t = p->mi->me[i].source->node_index;
        gtree->migcount[s][t]--;

        if (p->mi->me[i].source != curpop)
        {
          for (pop = curpop->parent; pop != p->mi->me[i].source->parent; pop = pop->parent)
            pop->seqin_count[msa_index]--;
        }
        curpop = p->mi->me[i].target;
        migevent_unlink(p->mi->me+i,msa_index);
      }
      p->mi->count = 0;
    }
    for (pop = curpop->parent; pop != father->pop->parent; pop = pop->parent)
      pop->seqin_count[msa_index]--;
  }
  return simcount;
}
#endif

static gnode_t * subtree_prune_from_pop(stree_t * stree,
                                        gtree_t * gtree,
                                        snode_t * frompop,
                                        gnode_t * curnode,
                                        double tL)
{
  long i;
  long msa_index = gtree->msa_index;

  gnode_t * father;
  gnode_t * sibling = NULL;
  gnode_t * tmp;
  snode_t * pop;
  snode_t * curpop;
                      
  father  = curnode->parent;

  if (father)
    sibling = (father->left == curnode) ? father->right : father->left;

  /* this function assumes the lineages from curnode passes through frompop */

  /* correctness check */
  #if 1
  if (father)
  {
    if (opt_migration && !father->parent)
      assert(gtree->root == father);
  }
  #endif

  /* remove father from coalescent events of its population */
  if (father)
  {
    unlink_event(father,msa_index);
    father->pop->event_count[msa_index]--;
    if (!opt_est_theta)
      father->pop->event_count_sum--;
  }

  /* decrease the number of incoming lineages to all populations in the path from
  frompop population (excluding) to the father population (including) */
  curpop = frompop;
  snode_t * fatherpop_parent = father ? father->pop->parent : NULL;
  if (opt_migration && curnode->mi && curnode->mi->count)
  {
    /* skip migrations until we get into frompop */
    pop = curnode->pop;
    for (i = 0; i < curnode->mi->count; ++i)
    {
      /* TODO: This loop is identical to the loop in mig_dissolve_and_sim,
         with the exception of the checks for tU. I deemed them unnecessary
         here, since we know the lines from curnode definitely passes through
         the migration band as determined in mig_dissolve_and_sim. Perhaps
         check this again */
      if (curnode->mi->me[i].time > tL &&
          stree->pptable[pop->node_index][frompop->node_index]) break;
      pop = curnode->mi->me[i].target;
    }
    long remain_count = i;

    for (; i < curnode->mi->count; ++i)
    {
      curnode->mi->me[i].source->migevent_count[msa_index]--;
      curnode->mi->me[i].target->migevent_count[msa_index]--;

      /* decrease # of migration from s to t (forward in time) */
      long s = curnode->mi->me[i].target->node_index;
      long t = curnode->mi->me[i].source->node_index;
      gtree->migcount[s][t]--;

      if (curnode->mi->me[i].source != curpop)
      {
        for (pop = curpop->parent; pop != curnode->mi->me[i].source->parent; pop = pop->parent)
          pop->seqin_count[msa_index]--;
      }
      curpop = curnode->mi->me[i].target;
      migevent_unlink(curnode->mi->me+i,msa_index);
    }
    curnode->mi->count = remain_count;
  }
  for (pop = curpop->parent; pop != fatherpop_parent; pop = pop->parent)
    pop->seqin_count[msa_index]--;

  if (father)
  {
    if (!father->parent)
    {
      /* if father is the root, sibling becomes the new root */
      SWAP(gtree->root->pmatrix_index, sibling->pmatrix_index);
      gtree->root = sibling;
      sibling->parent = NULL;
    }
    else
    {
      if (father->parent->left == father)
        father->parent->left = sibling;
      else
        father->parent->right = sibling;

      sibling->parent = father->parent;

      /* update number of leaves all nodes from father's parent and up */
      for (tmp = father->parent; tmp; tmp = tmp->parent)
        tmp->leaves = tmp->left->leaves + tmp->right->leaves;
    }

    curnode->parent = NULL;
    father->left = father->right = NULL;
    father->parent = NULL;

    /* move migration events from father to sibling */
    if (opt_migration && father->mi && father->mi->count)
    {
      miginfo_check_and_extend(&(sibling->mi), father->mi->count);
      for (i = 0; i < father->mi->count; ++i)
      {
        /* TODO: miginfo_move is dangerous (see notes in function), and assumes
           caller maintains count */
        miginfo_move(father->mi->me+i,&(sibling->mi));
      }
      father->mi->count = 0;
    }
  }

  /* return deleted node (or NULL when pruning the root lineage) */
  return father;
}

void pruneoff(stree_t * stree,
              gtree_t * gtree,
              snode_t * frompop,
              gnode_t * curnode,
              double tL,
              gnode_t ** deleted,
              unsigned int * dcountptr)
{
  long i;
  long msa_index = gtree->msa_index;
  gnode_t * father;
  gnode_t * pruned = NULL;
  gnode_t * sibling = NULL;
  gnode_t * tmp;
  snode_t * pop;
  snode_t * curpop;

  unsigned int p = *dcountptr;
  unsigned int pinit = p;

  father = curnode->parent;
  if (father)
    sibling = (father->left == curnode) ? father->right : father->left;

  /* this function assumes the lineages from curnode passes through frompop */

  /* correctness check */
  #if 1
  if (father)
  {
    if (opt_migration && !father->parent)
      assert(gtree->root == father);
  }
  #endif

  /* remove ancestors from coalescent events of its population */
  pruned = curnode;
  while (pruned->parent)
  {
    if ((pruned->parent->left == curnode && pruned->parent->mark & FLAG_RED_LEFT) ||
        (pruned->parent->right == curnode && pruned->parent->mark & FLAG_RED_RIGHT))
      break;
    
    deleted[p++] = pruned->parent;

    /* remove ancestor from coalescent events of its population */
    unlink_event(pruned->parent,msa_index);
    pruned->parent->pop->event_count[msa_index]--;
    if (!opt_est_theta)
      pruned->parent->pop->event_count_sum--;

    pruned = pruned->parent;
  }
  if (pruned->parent)
  {
    unlink_event(pruned->parent,msa_index);
    pruned->parent->pop->event_count[msa_index]--;
    if (!opt_est_theta)
      pruned->parent->pop->event_count_sum--;

    sibling = (pruned->parent->left == pruned) ? pruned->parent->right : pruned->parent->left;
    deleted[p++] = pruned->parent;
    if (pruned->parent->parent)
    {
      if (sibling == pruned->parent->left && pruned->parent->mark & FLAG_RED_LEFT)
      {
        if (pruned->parent->parent->left == pruned->parent)
          pruned->parent->parent->mark |= FLAG_RED_LEFT;
        else
          pruned->parent->parent->mark |= FLAG_RED_RIGHT;
      }
      if (sibling == pruned->parent->right && pruned->parent->mark & FLAG_RED_RIGHT)
      {
        if (pruned->parent->parent->left == pruned->parent)
          pruned->parent->parent->mark |= FLAG_RED_LEFT;
        else
          pruned->parent->parent->mark |= FLAG_RED_RIGHT;
      }
    }
  }

  /* pruned is the last node on the string curnode---pruned---X where X is the
     father of pruned or NULL, and pruned can be curnode */

  /* decrease the number of incoming lineages to all populations in the path from
  frompop population (excluding) to the pruned node population (including) */
  curpop = frompop;
  snode_t * fatherpop_parent = pruned->parent ? pruned->parent->pop->parent : NULL;
  if (opt_migration && pruned->mi && pruned->mi->count)
  {
    /* skip migrations until we get into frompop */
    pop = pruned->pop;
    for (i = 0; i < pruned->mi->count; ++i)
    {
      /* TODO: This loop is identical to the loop in mig_dissolve_and_sim,
         with the exception of the checks for tU. I deemed them unnecessary
         here, since we know the lines from pruned definitely passes through
         the migration band as determined in mig_dissolve_and_sim. Perhaps
         check this again */
      if (pruned->mi->me[i].time > tL &&
          stree->pptable[pop->node_index][frompop->node_index]) break;
      pop = pruned->mi->me[i].target;
    }
    long remain_count = i;

    for (; i < pruned->mi->count; ++i)
    {
      pruned->mi->me[i].source->migevent_count[msa_index]--;
      pruned->mi->me[i].target->migevent_count[msa_index]--;

      /* decrease # of migration from s to t (forward in time) */
      long s = pruned->mi->me[i].target->node_index;
      long t = pruned->mi->me[i].source->node_index;
      gtree->migcount[s][t]--;

      if (pruned->mi->me[i].source != curpop)
      {
        for (pop = curpop->parent; pop != pruned->mi->me[i].source->parent; pop = pop->parent)
          pop->seqin_count[msa_index]--;
      }
      curpop = pruned->mi->me[i].target;
      migevent_unlink(pruned->mi->me+i,msa_index);
    }
    pruned->mi->count = remain_count;
  }
  for (pop = curpop->parent; pop != fatherpop_parent; pop = pop->parent)
    pop->seqin_count[msa_index]--;

  father = pruned->parent;

  /* concatenate string of nodes */
  if (curnode != pruned)
  {
    /* delete array of nodes */
    gnode_t * gnode = curnode->parent;
    while (gnode != father)
    {
      #if 0
      sibling = (gnode->parent->left == gnode) ?
                  gnode->parent->right : gnode->parent->left;
      #endif

      if (!gnode->parent)
      {
        #if 0
        /* TODO: This shouldn't happen */
        assert(0);
        /* if gnode is the root, sibling becomes the new root */
        SWAP(gtree->root->pmatrix_index, sibling->pmatrix_index);
        gtree->root = sibling;
        sibling->parent = NULL;
        #endif
      }
      else
      {
        #if 0
        if (gnode->parent->left == gnode)
          gnode->parent->left = sibling;
        else
          gnode->parent->right = sibling;

        sibling->parent = gnode->parent;

        /* update number of leaves all nodes from gnode's parent and up */
        for (tmp = gnode->parent; tmp; tmp = tmp->parent)
          tmp->leaves = tmp->left->leaves + tmp->right->leaves;
        #endif
      }

      curnode->parent = gnode->parent;
      assert(!(gnode->parent->mark & FLAG_RED_LEFT && gnode->parent->mark & FLAG_RED_RIGHT));
      if (gnode->parent->left == gnode)
      {
        assert(gnode->parent->mark & FLAG_RED_RIGHT);
        gnode->parent->left = curnode;
      }
      else
      {
        assert(gnode->parent->mark & FLAG_RED_LEFT);
        gnode->parent->right = curnode;
      }

      gnode_t * f = gnode->parent;

      gnode->parent = NULL;
      assert(gnode->mark & FLAG_RED_LEFT || gnode->mark & FLAG_RED_RIGHT);
      gnode->left  = NULL;
      gnode->right = NULL;

      gnode = f;
    }
  }
  
  if (father)
  {
    sibling = (father->left == pruned) ? father->right : father->left;
    /* TODO: Assert that the sibling of sibling is curnode */
    assert(sibling->parent->left == curnode || sibling->parent->right == curnode);
  }

  if (father)
  {
    if (!father->parent)
    {
      /* if father is the root, sibling becomes the new root */
      if (father == gtree->root)
      {
        SWAP(gtree->root->pmatrix_index, sibling->pmatrix_index);
        gtree->root = sibling;
      }
      sibling->parent = NULL;
    }
    else
    {
      if (father->parent->left == father)
        father->parent->left = sibling;
      else
        father->parent->right = sibling;

      sibling->parent = father->parent;

      /* update number of leaves all nodes from father's parent and up */
      for (tmp = father->parent; tmp; tmp = tmp->parent)
        tmp->leaves = tmp->left->leaves + tmp->right->leaves;
    }

    pruned->parent = NULL;
    father->left = father->right = NULL;
    father->parent = NULL;

    /* move migration events from father to sibling */
    if (opt_migration && father->mi && father->mi->count)
    {
      miginfo_check_and_extend(&(sibling->mi), father->mi->count);
      for (i = 0; i < father->mi->count; ++i)
      {
        /* TODO: miginfo_move is dangerous (see notes in function), and assumes
           caller maintains count */
        miginfo_move(father->mi->me+i,&(sibling->mi));
      }
      father->mi->count = 0;
    }
  }

  *dcountptr = p;

  for (i = pinit; i < p; ++i)
  {
    deleted[i]->old_pop = deleted[i]->pop;
    deleted[i]->pop = NULL;
  }


}

gnode_t * pruneoff2(stree_t * stree,
                    gtree_t * gtree,
                    snode_t * frompop,
                    gnode_t * curnode,
                    double tL)
{
  long i;
  long msa_index = gtree->msa_index;

  gnode_t * father;
  gnode_t * sibling = NULL;
  gnode_t * tmp;
  snode_t * pop;
  snode_t * curpop;
                      
  gnode_t * prevnode = curnode;                    
  father  = curnode->parent;

  while (father)
  {
    if ((father->left == prevnode  && father->mark & FLAG_RED_LEFT) ||
        (father->right == prevnode && father->mark & FLAG_RED_RIGHT))
      break;
    prevnode = father;
    father = father->parent;
  }
  if (father)
    sibling = (father->left == prevnode) ? father->right : father->left;

  /* this function assumes the lineages from curnode passes through frompop */

  /* correctness check */
  #if 1
  /* 12.6.2023 We do not swap the root pmatrix index, we fix it at the end */
  //if (father)
  //{
  //  if (opt_migration && !father->parent)
  //    assert(gtree->root == father);
  //}
  #endif

  /* remove father from coalescent events of its population */
  if (father)
  {
    unlink_event(father,msa_index);
    father->pop->event_count[msa_index]--;
    if (!opt_est_theta)
      father->pop->event_count_sum--;
  }

  /* decrease the number of incoming lineages to all populations in the path from
  frompop population (excluding) to the father population (including) */
  curpop = frompop;
  snode_t * fatherpop_parent = father ? father->pop->parent : NULL;
  if (opt_migration && prevnode->mi && prevnode->mi->count)
  {
    /* skip migrations until we get into frompop */
    pop = prevnode->pop;
    for (i = 0; i < prevnode->mi->count; ++i)
    {
      /* TODO: This loop is identical to the loop in mig_dissolve_and_sim,
         with the exception of the checks for tU. I deemed them unnecessary
         here, since we know the lines from prevnode definitely passes through
         the migration band as determined in mig_dissolve_and_sim. Perhaps
         check this again */
      if (prevnode->mi->me[i].time > tL &&
          stree->pptable[pop->node_index][frompop->node_index]) break;
      pop = prevnode->mi->me[i].target;
    }
    long remain_count = i;

    for (; i < prevnode->mi->count; ++i)
    {
      prevnode->mi->me[i].source->migevent_count[msa_index]--;
      prevnode->mi->me[i].target->migevent_count[msa_index]--;

      /* decrease # of migration from s to t (forward in time) */
      long s = prevnode->mi->me[i].target->node_index;
      long t = prevnode->mi->me[i].source->node_index;
      gtree->migcount[s][t]--;

      if (prevnode->mi->me[i].source != curpop)
      {
        for (pop = curpop->parent; pop != prevnode->mi->me[i].source->parent; pop = pop->parent)
          pop->seqin_count[msa_index]--;
      }
      curpop = prevnode->mi->me[i].target;
      migevent_unlink(prevnode->mi->me+i,msa_index);
    }
    prevnode->mi->count = remain_count;
  }
  for (pop = curpop->parent; pop != fatherpop_parent; pop = pop->parent)
    pop->seqin_count[msa_index]--;

  if (father)
  {
    if (!father->parent)
    {
      /* if father is the root, sibling becomes the new root */
      /* 12.6.2023 We do not swap the root pmatrix index, we fix it at the end */
      //SWAP(gtree->root->pmatrix_index, sibling->pmatrix_index);
      //gtree->root = sibling;
      sibling->parent = NULL;
    }
    else
    {
      if (father->parent->left == father)
        father->parent->left = sibling;
      else
        father->parent->right = sibling;

      sibling->parent = father->parent;

      /* update number of leaves all nodes from father's parent and up */
      for (tmp = father->parent; tmp; tmp = tmp->parent)
        tmp->leaves = tmp->left->leaves + tmp->right->leaves;
    }

    prevnode->parent = NULL;
    father->left = father->right = NULL;
    father->parent = NULL;

    /* move migration events from father to sibling */
    if (opt_migration && father->mi && father->mi->count)
    {
      miginfo_check_and_extend(&(sibling->mi), father->mi->count);
      for (i = 0; i < father->mi->count; ++i)
      {
        /* TODO: miginfo_move is dangerous (see notes in function), and assumes
           caller maintains count */
        miginfo_move(father->mi->me+i,&(sibling->mi));
      }
      father->mi->count = 0;
    }
  }

  /* return deleted node (or NULL when pruning the root lineage) */
  if (father)
  {
    sibling->mark |= FLAG_BRANCH_UPDATE;
  }
  return father;
}

#if 0
/* TODO: Probably delete */
void mig_simulate()
{
  snode_t * curpop = gnode->pop;
  for (j = 0; j < gnode->mi->count && curpop != x; ++j)
  {
    if (gnode->mi->me[j].time > x->tau &&
        stree->pptable[curpop->node_index][x->node_index]) break;
    curpop = gnode->mi->me[j].target;
  }
  assert(stree->pptable[curpop->node_index][x->node_index]);

  t = x->tau;
  if (gnode->mi && gnode->mi->count)
  {
    /* check whether we are entering population x via migration */
    long last_mig_index = gnode->mi->count-1;
    if (gnode->mi->me[last_mig_index].target == x)
      t = gnode->mi->me[last_mig_index].time;
  }

  
  
}
#endif

#if 1
/* determine whether the parental lineage associated with gnode passes (backwards in time) through the origin (youngest point) of snode.
   Note: a node (coalescent event) within snode means that the lineage does not originate in snode.

   This needs to be very fast as it is called many times for many nodes */
static long originates(gnode_t * gnode, snode_t * snode, gtree_t * gtree, stree_t * stree)
{
  long i;
  /* if it's a tip in snode return true, if tip in other node return false */
  if (gnode->node_index < gtree->tip_count)
    return gnode->pop == snode;

  if (gnode->time > snode->tau) return 0;

  if (gnode->pop == snode) { assert(0); } /* XXX: This should never happend because of the previous condition */

  /* if parent exists and is younger than snode origin return false */
  if (gnode->parent && gnode->parent->time < snode->tau)
    return 0;

  /* if our node has no migrations */
  if (!gnode->mi || !gnode->mi->count)
    return stree->pptable[gnode->pop->node_index][snode->node_index];

  /* if there are migrations */
  assert(gnode->mi && gnode->mi->count);

  snode_t * lopop = gnode->pop;
  for (i = 0; i < gnode->mi->count; ++i)
  {
    if (gnode->mi->me[i].time > snode->tau)
      break;

    lopop = gnode->mi->me[i].target;
  }

  return stree->pptable[lopop->node_index][snode->node_index];
}
#endif

#if 1
/* No change from the routine in gtree.c */
static migbuffer_t * wtimes_and_lineages(stree_t * stree,
                                         gtree_t * gtree,
                                         snode_t * snode,
                                         gnode_t * gnode,
                                         snode_t * affected1,
                                         snode_t * affected2,
                                         snode_t ** affected_nodemap,
                                         double t,
                                         long existing_mig_count,
                                         long * lineages_count,
                                         long * wtimes_count,
                                         long msa_index,
                                         long thread_index)
{
  long i,j,k;
  long lineages;
  long mig_count = 0;
  migbuffer_t * wtimes;

  /* make sure migbuffer is large enough */
  size_t alloc_required = snode->migevent_count[msa_index] +
                          snode->event_count[msa_index] +
                          stree->inner_count+1;
  migbuffer_check_and_realloc(thread_index,alloc_required);
  wtimes = global_migbuffer_r[thread_index];

  /* TODO: Note this was wrong. The seqin_count is not the number of lineages
  entering the population as we pruned those lineages */
  /* start with the lineages entering the current population */
  lineages = snode->seqin_count[msa_index];

  /* make a list of waiting times older than t and at the same time reduces
     number of lineages by the amount of coalescent events younger than t
     (including t) */
  k = 0;
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    
    #if 1
    //if (x == gnode && x->pop == snode && snode == affected && x->mark & FLAG_PRUNED && !(x->mark & FLAG_MIGRATE))
    //if (x->pop == affected && snode == affected && x->mark & FLAG_PRUNED && !(x->mark & FLAG_MIGRATE))
    #if 0
    if (snode == affected && x->mark & FLAG_PRUNED && !(x->mark & FLAG_MIGRATE))
    #else
    /* this is additional term for the flipping move */
    if (affected2)
    {
      assert(affected_nodemap);
      /* flipping move */
      if ((snode == affected1 || snode == affected2) && x->mark & FLAG_PRUNED && !(x->mark & FLAG_MIGRATE))
      {
        #if 0
        snode_t * lastpop = (x->mi && x->mi->count) ? x->mi->me[x->mi->count-1].target: x->pop;
        if (stree->pptable[lastpop->node_index][snode->node_index])
          --lineages;
        #endif

        #if 0
        snode_t * lastpop = (x->mi && x->mi->count) ? x->mi->me[x->mi->count-1].target: x->pop;
        if (!(x->mi && x->mi->count && x->mi->me[x->mi->count-1].time > MAX(affected1->tau,affected2->tau)) || lastpop == snode)
        {
          if (stree->pptable[lastpop->node_index][snode->node_index])
            --lineages;
        }
        #else
        assert(affected_nodemap[x->node_index]);
        if (affected_nodemap[x->node_index] == snode)
          --lineages;
        #endif
      }
    }
    else
    {
      assert(!affected_nodemap);
      /* add or remove case */
      if (snode == affected1 && x->mark & FLAG_PRUNED && !(x->mark & FLAG_MIGRATE))
        --lineages;
    }
    #endif
    {
      #if 0
      if (i >= gtree->tip_count && x->pop == snode)
       --lineages;
      continue;
      #endif
    }
    #else
    #if 1
    if ((x->mark & FLAG_PRUNED) &&
        x->pop == affected &&
        affected == snode)
          --lineages;
    #else
    if (x->mark & FLAG_PRUNED)
    {
      //if (x->pop == affected && affected == snode)
      if (!(x->mark & FLAG_MIGRATE))
      {
        if (affected == snode)
          --lineages;
      }
      else
      {
        /* it migrates, which means */
        if (snode != affected && enters(x,snode) && x==gnode)
          lineages--;
      }
    }
    #endif
    #endif


    /* skip pruned node */
    /* TODO: Note the last x != gtree->root is very important. This is because we may have
       x may be a tip and the sibling of the pruned subtree as well as a child of the root
       in which case it becomes the new root. However, it has x->left == NULL, x->right == NULL
       and x->parent == NULL */
    // if (!x->left && !x->right && !x->parent && x != gtree->root) continue;

    /* coalescent events */
    if (i >= gtree->tip_count)
    {
      if (x->pop == snode)
      {
        if (x->time > t)
        {
          wtimes[k].time   = x->time;
          wtimes[k].mrsum  = 0;
          wtimes[k++].type = EVENT_COAL;
        }
        else
          --lineages;
      }
    }

    #if 0
    /* if it's the node on which we are simulating the lineage, skip any migration events
       as they are those that were simulated so far */
    if (x == gnode)
    {
      assert(x->mark & FLAG_PRUNED);
      if (!(x->mark & FLAG_MIGRATE))
        continue;
    }
    #endif

    /* migrations */
    miginfo_t * mi = x->mi;
    if (!mi || !mi->count) continue;

    /* Note, if the node is flagged as pruned but there are migration events assigned on it,
       the last migration event is one entering population X (in the sense we have created
       or deleted a migration rate M_Y->X. In that case, we must skip the last migration event
       entering population X as it would otherwise increase the number of lineages in X */
    mig_count = mi->count;
    if (x == gnode)
    {
      /* account only for the migrations that existed prior to the start of simulation */
      #if 0
      assert((x->mark & FLAG_PRUNED) && (x->mark & FLAG_MIGRATE));
      #else
      assert(x->mark & FLAG_PRUNED);
      #endif
      mig_count = existing_mig_count;
    }

    for (j = 0; j < mig_count; ++j)
    {
      if (mi->me[j].source != snode && mi->me[j].target != snode) continue;

      /* skip migration entering to affected population */
      /* TODO: Check for flipping node if affected2 should be used */
      if (affected2)
      {
        if (x->mark & FLAG_MIGRATE && (snode == affected1 || snode == affected2) && j == mig_count-1) continue;
      }
      else
      {
        if (x->mark & FLAG_MIGRATE && snode == affected1 && j == mig_count-1) continue;
      }
      
      if (mi->me[j].time > t)
      {
        wtimes[k].time   = mi->me[j].time;
        wtimes[k].mrsum  = 0;
        wtimes[k++].type = (mi->me[j].source == snode) ?
                             EVENT_MIG_SOURCE : EVENT_MIG_TARGET;
      }
      else
      {
        if (mi->me[j].source == snode)
          --lineages;
        else
          ++lineages;
      }
    }
  }
  
  /* add taus dividing the population */
  for (i = 0; i < snode->mb_count; ++i)
  {
    if (snode->migbuffer[i].time > t)
    {
      wtimes[k++] = snode->migbuffer[i];
    }
  }

  qsort(wtimes,k,sizeof(migbuffer_t),cb_migbuf_asctime);

  *lineages_count = lineages;
  *wtimes_count = k;
  return wtimes;
}
#else
/* No change from the routine in gtree.c */
static migbuffer_t * wtimes_and_lineages(stree_t * stree,
                                         gtree_t * gtree,
                                         snode_t * snode,
                                         gnode_t * gnode,
                                         snode_t * affected,
                                         double t,
                                         long existing_mig_count,
                                         long * lineages_count,
                                         long * wtimes_count,
                                         long msa_index,
                                         long thread_index)
{
  long i,j,k;
  long lineages;
  long mig_count = 0;
  migbuffer_t * wtimes;

  /* make sure migbuffer is large enough */
  size_t alloc_required = snode->migevent_count[msa_index] +
                          snode->event_count[msa_index] +
                          stree->inner_count+1;
  migbuffer_check_and_realloc(thread_index,alloc_required);
  wtimes = global_migbuffer_r[thread_index];

  /* TODO: Note this was wrong. The seqin_count is not the number of lineages
  entering the population as we pruned those lineages */
  /* start with the lineages entering the current population */
  lineages = snode->seqin_count[msa_index];

  /* make a list of waiting times older than t and at the same time reduces
     number of lineages by the amount of coalescent events younger than t
     (including t) */
  k = 0;
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    #if 0
    if (x == gnode && x->pop == snode) --lineages;
    #else
    #if 1
    if ((x->mark & FLAG_PRUNED) &&
        x->pop == affected &&
        affected == snode)
          --lineages;
    #else
    if (x->mark & FLAG_PRUNED)
    {
      //if (x->pop == affected && affected == snode)
      if (!(x->mark & FLAG_MIGRATE))
      {
        if (affected == snode)
          --lineages;
      }
      else
      {
        /* it migrates, which means */
        if (snode != affected && enters(x,snode) && x==gnode)
          lineages--;
      }
    }
    #endif
    #endif
    /* skip pruned node */
    /* TODO: Note the last x != gtree->root is very important. This is because we may have
       x may be a tip and the sibling of the pruned subtree as well as a child of the root
       in which case it becomes the new root. However, it has x->left == NULL, x->right == NULL
       and x->parent == NULL */
    //if (!x->left && !x->right && !x->parent && x != gtree->root) continue;

    /* coalescent events */
    if (i >= gtree->tip_count)
    {
      if (x->pop == snode)
      {
        if (x->time > t)
        {
          wtimes[k].time   = x->time;
          wtimes[k].mrsum  = 0;
          wtimes[k++].type = EVENT_COAL;
        }
        else
          --lineages;
      }
    }

    /* if it's the node on which we are simulating the lineage, skip any migration events
       as they are those that were simulated so far */
    if (x == gnode)
    {
      continue;
    }

    /* migrations */
    miginfo_t * mi = x->mi;
    if (!mi || !mi->count) continue;

    /* Note, if the node is flagged as pruned but there are migration events assigned on it,
       the last migration event is one entering population X (in the sense we have created
       or deleted a migration rate M_Y->X. In that case, we must skip the last migration event
       entering population X as it would otherwise increase the number of lineages in X */
    mig_count = mi->count;
    if (snode == affected && x->mark & FLAG_PRUNED && x->mark & FLAG_MIGRATE)
    {
      assert(mi->count);
      assert(x->mi->me[mi->count-1].target  == affected);
      --mig_count;
    }
       

    for (j = 0; j < mig_count; ++j)
    {
      if (mi->me[j].source != snode && mi->me[j].target != snode) continue;
      
      if (mi->me[j].time > t)
      {
        wtimes[k].time   = mi->me[j].time;
        wtimes[k].mrsum  = 0;
        wtimes[k++].type = (mi->me[j].source == snode) ?
                             EVENT_MIG_SOURCE : EVENT_MIG_TARGET;
      }
      else
      {
        if (mi->me[j].source == snode)
          --lineages;
        else
          ++lineages;
      }
    }
  }
  
  /* add taus dividing the population */
  for (i = 0; i < snode->mb_count; ++i)
  {
    if (snode->migbuffer[i].time > t)
    {
      wtimes[k++] = snode->migbuffer[i];
    }
  }

  qsort(wtimes,k,sizeof(migbuffer_t),cb_migbuf_asctime);

  *lineages_count = lineages;
  *wtimes_count = k;
  return wtimes;
}
#endif

/* passing variable t and starting from population x */
static double simulate_coalescent_mig_from_point(stree_t * stree,
                                                 gtree_t * gtree,
                                                 gnode_t * gnode,
                                                 snode_t * affected1,
                                                 snode_t * affected2,
                                                 snode_t ** affected_nodemap,
                                                 double t,
                                                 long msa_index,
                                                 long thread_index)
{
  long i,j,k;
  long lineages;
  long epoch;
  long mpop_count;
  double crate,mrate,rate;
  double tnew;
  migbuffer_t * wtimes;
  snode_t * snode = affected1;

  snode_t ** migsource = gtree->migpops;

  assert(opt_migration);

  assert(gnode->mark & FLAG_PRUNED);
  #if 0
  long existing_mig_count = (gnode->mark & FLAG_MIGRATE) ? gnode->mi->count : 0;
  #else
  long existing_mig_count = gnode->mi ? gnode->mi->count : 0;
  #endif

  /* determine time to coalesce */
  while (1)
  {
    wtimes = wtimes_and_lineages(stree,gtree,snode,gnode,affected1,affected2,affected_nodemap,t,existing_mig_count,&lineages,&k,msa_index,thread_index);
    assert(lineages >= 0 && k >= 0);

    mrate = crate = rate = 0; epoch = 0;
    /* simulate until all waiting times are exhausted */
    for (i = 0; i < k; ++i)
    {
      /* TODO: The below 'if' line is wrong. Delete it.
         I think a simple if should do the trick, i.e.
           if (snode->migbuffer[j].time <= t) epoch++; */
      //if (i == 0 || wtimes[i].type == EVENT_TAU)
      {
        /* find current epoch and get total rate */
        assert(!snode->parent || snode->mb_count);
        for (j = 0; j < snode->mb_count; ++j)
          if (snode->migbuffer[j].time > t)
          {
            epoch = j;
            break;
          }
        /* find migrating populations */
        mpop_count = 0;
        for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        {
          snode_t * x = stree->nodes[j];
          if (x->parent && x != snode &&
              opt_mig_bitmatrix[x->node_index][snode->node_index] &&
              x->tau <= t && x->parent->tau > t)        /* Note: parent->tau > t is importnat */
            migsource[mpop_count++] = x;
        }
      }
      /* rate */
      long idx = snode->migbuffer[epoch].active_count == 1 ? 0 : msa_index;
      mrate = snode->parent ? 4*snode->migbuffer[epoch].mrsum[idx]/snode->theta : 0;

      crate = 2*lineages / snode->theta;
      rate = mrate + crate;

      if (!lineages && !mpop_count)
      {
        /* no lineages to coalesce and no pops to migrate */
        //assert(k==1 && rate == 0);
        t = wtimes[i].time;
        if (wtimes[i].type == EVENT_MIG_SOURCE)
          --lineages;
        else if (wtimes[i].type == EVENT_MIG_TARGET)
          ++lineages;
        else if (wtimes[i].type == EVENT_COAL)
          --lineages;

        continue;
      }

      /* generate time */
      tnew = legacy_rndexp(thread_index,1/rate);

      if (t+tnew <= wtimes[i].time)
      {
        t += tnew;
        break;
      }

      t = wtimes[i].time;

      if (wtimes[i].type == EVENT_MIG_SOURCE)
        --lineages;
      else if (wtimes[i].type == EVENT_MIG_TARGET)
        ++lineages;
      else if (wtimes[i].type == EVENT_COAL)
        --lineages;
    }
    
    if (i != k)
    {
      /* choose between coalescence and migration */
      double r = legacy_rndu(thread_index) * rate;

      /* if coalescence, stop the process */
      if (r >= mrate) break;

      /* choose migrating source (forward-time) population*/
      assert(mpop_count);
      double sum = 0;
      for (j = 0; j < mpop_count - 1; ++j)
      {
        assert(opt_mig_bitmatrix[migsource[j]->node_index][snode->node_index]);
        long mindex = opt_migration_matrix[migsource[j]->node_index][snode->node_index];
        if (opt_mig_specs[mindex].am)
          sum += 4*opt_mig_specs[mindex].Mi[msa_index] / snode->theta;
        else
          sum += 4*opt_mig_specs[mindex].M / snode->theta;
        if (r < sum) break;
      }

      /* store miginfo in current branch: snode is source in backward time */
      assert(migsource[j] != snode);
      miginfo_append(&(gnode->mi), snode, migsource[j], t, msa_index);

      snode = migsource[j];

      continue;
    }

    /* no migration, no coalescence */

    /* stop if we were in the root population */
    if (!snode->parent)
    {
      if (lineages == 0)
      {
      /* if this is the lineage going back to infinitity, i.e. first lineage
         we simulate that does not coalesce anywhere, return -1*/
        return -1;
      }
      assert(lineages == 1);

      /* generate time */
      rate = 2 / snode->theta;
      t += legacy_rndexp(thread_index, 1/rate);
      break;
    }
    snode = snode->parent;
  }

  return t;
}

static long target_branches(stree_t * stree,
                            gtree_t * gtree,
                            snode_t * snode,
                            gnode_t * curnode,
                            double t,
                            gnode_t ** outbuf)
{
  long i,j;
  long count = 0;
  snode_t * curpop = NULL;

  assert(!opt_msci);

  for (i = 0; i < gtree->tip_count; ++i)
  {
    gnode_t * gnode = gtree->nodes[i];

    /* skip lineages without migrations that do not pass through snode */
    if (!opt_migration)
      if (!stree->pptable[gnode->pop->node_index][snode->node_index])
        continue;


    while (!(gnode->mark & FLAG_MISC) && gnode->parent && gnode->parent->time < t)
    {
      gnode->mark |= FLAG_MISC;
      gnode = gnode->parent;
    }
    #if 0
    /* TODO: The following line is wrong, we need to check the flag instead */
    if ((gnode->mark & FLAG_MISC) || gnode == curnode) continue;
    #else
    //if ((gnode->mark & FLAG_MISC) || gnode->mark & FLAG_PRUNED) continue;
    if (gnode->mark & FLAG_MISC) continue;
    if (gnode->mark & FLAG_PRUNED && gnode->mark & FLAG_MIGRATE)
    {
      //if (!gnode->mi || !gnode->mi->count) continue;
      assert(gnode->mi && gnode->mi->count);

      long lastidx = gnode->mi->count-1;
      //assert(gnode->mi->me[idx].target == AFFECTED);
      double mt = gnode->mi->me[lastidx].time;

      if (mt < t) continue;
    }
    else if (gnode->mark & FLAG_PRUNED) continue;
    #endif

    
    curpop = gnode->pop;
    if (opt_migration && gnode->mi && gnode->mi->count)
    {
      assert(snode->tau < t && (!snode->parent || snode->parent->tau > t));
      assert(t > gnode->time && (!gnode->parent || gnode->parent->time > t));
      /* check if part of the lineage passes through snode at time t */
      for (j = 0; j < gnode->mi->count; ++j)
      {
        if (t < gnode->mi->me[j].time) break;
        curpop = gnode->mi->me[j].target;
      }
    }
    if (!stree->pptable[curpop->node_index][snode->node_index])
      continue;
    outbuf[count++] = gnode;
    gnode->mark |= FLAG_MISC;
  }

  /* clean marks */
  for (i = 0; i < gtree->tip_count; ++i)
  {
    gnode_t * gnode = gtree->nodes[i];

    while (gnode && (gnode->mark & FLAG_MISC))
    {
      gnode->mark &= ~FLAG_MISC;
      gnode = gnode->parent;
    }
  }

  return count;
}

static void subtree_regraft(stree_t * stree,
                            gtree_t * gtree,
                            gnode_t * curnode,
                            gnode_t * father,
                            gnode_t * target,
                            snode_t * newpop,
                            snode_t * affected,
                            long existing_migs_count,
                            double tnew)
{
  long i,k;
  long msa_index = gtree->msa_index;
  gnode_t * tmp;
  snode_t * pop;
  snode_t * curpop;

  father->parent = target->parent;
  target->parent = father;
  if (father->parent)
  {
    if (father->parent->left == target)
      father->parent->left = father;
    else
      father->parent->right = father;
  }

  father->left = target;
  father->right = curnode;

  father->leaves = father->left->leaves + father->right->leaves;
  curnode->parent = father;
  /* 12.6.2023 We do not swap the root pmatrix index, we fix it at the end */
  #if 0
  if (target == gtree->root)
  {
    /* father becomes new root */
    SWAP(gtree->root->pmatrix_index, father->pmatrix_index);
    gtree->root = father;
  }
  else
  #endif
  {
    for (tmp = father->parent; tmp; tmp = tmp->parent)
      tmp->leaves = tmp->left->leaves + tmp->right->leaves;
  }

  /* change the population for the current gene tree node */
  father->pop = newpop;
  father->time = tnew;

  /* move target migration events older than tnew onto the father branch */
  /* Assumes migration are sorted from younger one (which is the case)
     otherwise miginfo_move will have trouble  */
  if (opt_migration && target->mi && target->mi->count)
  {
    k = 0;
    for (i = 0; i < target->mi->count; ++i)
    {
      if (target->mi->me[i].time > tnew)
      {
        miginfo_move(target->mi->me+i, &(father->mi));

        /* keep track of how many migrations we copied */
        ++k;
      }
    }
    target->mi->count -= k;
  }

  /* now add the coalescent event to the new population, at the end */
  dlist_item_append(father->pop->event[msa_index],father->event);

  father->pop->event_count[msa_index]++;
  if (!opt_est_theta)
    father->pop->event_count_sum++;

  /* increase number of incoming lineages to all populations in the path from
   * affected population (excluding) to the father population (including) */
  curpop = affected;
  if (opt_migration && curnode->mi && curnode->mi->count)
  {
    for (i = existing_migs_count; i < curnode->mi->count; ++i)
    {
      curnode->mi->me[i].source->migevent_count[msa_index]++;
      curnode->mi->me[i].target->migevent_count[msa_index]++;

      /* increase # of migration from s to t (forward in time) */
      long s = curnode->mi->me[i].target->node_index;
      long t = curnode->mi->me[i].source->node_index;
      gtree->migcount[s][t]++;

      if (curnode->mi->me[i].source != curpop)
      {
        for (pop = curpop->parent; pop != curnode->mi->me[i].source->parent; pop = pop->parent)
          pop->seqin_count[msa_index]++;
      }
      curpop = curnode->mi->me[i].target;
    }
  }
  for (pop = curpop->parent; pop != father->pop->parent; pop = pop->parent)
    pop->seqin_count[msa_index]++;
}

static int cb_preorder_inner_only(gnode_t * gnode)
{
  if (!gnode->left)
    return 0;

  return 1; 
}

static int cb_postorder_inner_only(gnode_t * gnode)
{
  if (!gnode->left)
    return 0; 

  return 1;
}

void mig_dissolve_and_sim(stree_t * stree,
                          gtree_t * gtree,
                          locus_t * locus,
                          double tL,
                          double tU,
                          snode_t * x,
                          int msa_index,
                          long thread_index,
                          double * ret_logl_diff,
                          double * ret_logpr_diff)
{
  unsigned int i,j,k,count,p;
  unsigned int root_pmat_index;
  unsigned int total_nodes = gtree->tip_count+gtree->inner_count;
  double prevt;
  double tnew;
  double logl, logpr;
  gnode_t ** simnodes = NULL;
  gnode_t ** parents = NULL;
  long * migration = NULL;
  snode_t * prevpop = NULL;

  simnodes = (gnode_t **)xcalloc((size_t)total_nodes, sizeof(gnode_t *));
  migration = (long *)xcalloc((size_t)total_nodes, sizeof(long));

  /* store pmatrix index for root */
  root_pmat_index = gtree->root->pmatrix_index;

  /* we have two cases when tracing the lineage backwards in time:
     1. Lineage enters population x from a descendant population y (in the path curpop<->x)
     2. Lineage migrates to other populations before migrating back to x

     For case (1) migrations may occur on the path curpop<->y before the
     lineage enters x.
  */

  /* 1. find lineages to dissolve */
/*** Ziheng $$$ ***/
  #if(DEBUG_zy)
  printf("[DEBUG] -- mig_dissolve_and_sim -- START\n");
  debug_print_stree(stree);
  debug_print_gtree(gtree);
  printf("[DEBUG] -- mig_dissolve_and_sim -- END\n");
  #endif

  for (i = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->nodes[i];

    if (gnode->time > tU) continue;
    if (gnode->pop == x && gnode->time > tL) continue;
    if (gnode->parent && gnode->parent->time < tL) continue;

    prevt = gnode->time;
    prevpop = gnode->pop;
    if (opt_migration && gnode->mi && gnode->mi->count)
    {
      for (j = 0; j < gnode->mi->count; ++j)
      {
        if (gnode->mi->me[j].time > tL && prevt < tU &&
            stree->pptable[prevpop->node_index][x->node_index])
          break;

        prevt = gnode->mi->me[j].time;
        prevpop = gnode->mi->me[j].target;
      }
    }
    if (!stree->pptable[prevpop->node_index][x->node_index] ||
        prevt > tU) continue;

    if (gnode->parent)
    {
      gnode->parent->mark |= (gnode->parent->left == gnode) ? 
                               FLAG_RED_LEFT : FLAG_RED_RIGHT;
      
      /* unset FLAG_SIMULATE on parent if both its daughter lineages are red */
      if ((gnode->parent->mark & FLAG_SIMULATE) &&
          (gnode->parent->mark & FLAG_RED_LEFT) &&
          (gnode->parent->mark & FLAG_RED_RIGHT))
      {
        gnode->parent->mark &= ~FLAG_SIMULATE;
        if (gnode->parent->mark & FLAG_MIGRATE)
          gnode->parent->mark &= ~FLAG_MIGRATE;
      }
    }

    /* set SIMULATE flag only if at most one RED daughter lineage */
    if (!((gnode->mark & FLAG_RED_LEFT) && (gnode->mark & FLAG_RED_RIGHT)))
    {
      gnode->mark |= FLAG_SIMULATE;

      if (prevt > tL)
        gnode->mark |= FLAG_MIGRATE;
    }
  }

  /* color the lineages with both daughter lineages red */ 
  gtree_traverse(gtree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_postorder_inner_only,
                 gtree->travbuffer,
                 &total_nodes);
  for (i = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->travbuffer[i];

    if (gnode->left->mark & FLAG_RED_LEFT && gnode->left->mark & FLAG_RED_RIGHT)
      gnode->mark |= FLAG_RED_LEFT;

    if (gnode->right->mark & FLAG_RED_LEFT && gnode->right->mark & FLAG_RED_RIGHT)
      gnode->mark |= FLAG_RED_RIGHT;
  }

  /* 1b. Go through nodes with SIMULATE flag in preorder */
  gtree_traverse(gtree->root,
                 TREE_TRAVERSE_PREORDER,
                 cb_preorder_inner_only,
                 gtree->travbuffer,
                 &total_nodes);
  for (i = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->travbuffer[i];

    //assert(gnode->mark & FLAG_SIMULATE);

    if ((gnode->mark & FLAG_RED_LEFT) &&
        (gnode->mark & FLAG_RED_RIGHT))
    {
      /* unset FLAG_SIMULATE if both daughter lineages are red */
      if (gnode->mark & FLAG_SIMULATE)
        gnode->mark &= ~FLAG_SIMULATE;
      if (gnode->mark & FLAG_MIGRATE)
        gnode->mark &= ~FLAG_MIGRATE;
    }
    else if (gnode->mark & FLAG_RED_LEFT)
    {
      /* move FLAG_SIMULATE from gnode to right daughter lineages */
      if (gnode->mark & FLAG_SIMULATE)
      {
        gnode->mark &= ~FLAG_SIMULATE;
        gnode->right->mark |= FLAG_SIMULATE;
      }
      if (gnode->mark & FLAG_MIGRATE)
      {
        gnode->mark &= ~FLAG_MIGRATE;
        gnode->right->mark |= FLAG_MIGRATE;
      }
    }
    else if (gnode->mark & FLAG_RED_RIGHT)
    {
      /* move FLAG_SIMULATE from gnode to left daughter lineages */
      if (gnode->mark & FLAG_SIMULATE)
      {
        gnode->mark &= ~FLAG_SIMULATE;
        gnode->left->mark |= FLAG_SIMULATE;
      }
      if (gnode->mark & FLAG_MIGRATE)
      {
        gnode->mark &= ~FLAG_MIGRATE;
        gnode->left->mark |= FLAG_MIGRATE;
      }
    }
  }

  /* 1c. create the simnodes/migration arrays */
  total_nodes = gtree->tip_count + gtree->inner_count;
  for (i = 0, count = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->nodes[i];

    if (gnode->mark & FLAG_SIMULATE)
    {
      if (gnode->mark & FLAG_MIGRATE)
        migration[count] = 1;
      simnodes[count++] = gnode;
    }

    /* keep red marks */
    gnode->mark &= (FLAG_RED_LEFT | FLAG_RED_RIGHT);
  }

  parents  = (gnode_t **)xcalloc((size_t)count, sizeof(gnode_t *));

  /* 2. Prune the lineages, and store the deleted nodes */
  for (i = 0, p = 0; i < count; ++i)
  {
    gnode_t * deleted = pruneoff2(stree,gtree,x,simnodes[i],tL);
    if (deleted)
    {
      parents[p++] = deleted;
      deleted->old_pop = deleted->pop;
      deleted->pop = NULL;
    }

    simnodes[i]->mark |= FLAG_PRUNED;
    if  (migration[i])
      simnodes[i]->mark |= FLAG_MIGRATE;

    /* mark the current node for branch length update */
    simnodes[i]->mark |= FLAG_BRANCH_UPDATE;
    
  }

  /* unset color flags from nodes */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    gtree->nodes[i]->mark &= ~(FLAG_RED_LEFT | FLAG_RED_RIGHT);

  /* 3. re-simulate */
  for (p = 0, i = 0; i < count; ++i)
  {
    long existing_migs_count = simnodes[i]->mi ? simnodes[i]->mi->count : 0;
    tnew = simulate_coalescent_mig_from_point(stree,
                                              gtree,
                                              simnodes[i],
                                              x,
                                              NULL,
                                              NULL,
                                              simnodes[i]->mark & FLAG_MIGRATE ? simnodes[i]->mi->me[existing_migs_count-1].time : tL,               /* TODO: Check what we are passing here is correct */
                                              msa_index,
                                              thread_index);

    if (tnew < 0)
    {
      /* if this is the lineage going back to infinitity, i.e. first lineage
         we simulate that does not coalesce anywhere */

      /* delete flag */
      if (simnodes[i]->mark & FLAG_PRUNED)
        simnodes[i]->mark &= ~FLAG_PRUNED;
      if (simnodes[i]->mark & FLAG_MIGRATE)
        simnodes[i]->mark &= ~FLAG_MIGRATE;

      /* TODO: Update seqin counts */
      snode_t * pop;
      snode_t * curpop = x;
      if (opt_migration && simnodes[i]->mi && simnodes[i]->mi->count)
      {
        for (j = existing_migs_count; j < simnodes[i]->mi->count; ++j)
        {
          simnodes[i]->mi->me[j].source->migevent_count[msa_index]++;
          simnodes[i]->mi->me[j].target->migevent_count[msa_index]++;

          /* increase # of migration from s to t (forward in time) */
          long s = simnodes[i]->mi->me[j].target->node_index;
          long t = simnodes[i]->mi->me[j].source->node_index;
          gtree->migcount[s][t]++;

          if (simnodes[i]->mi->me[j].source != curpop)
          {
            for (pop = curpop->parent; pop != simnodes[i]->mi->me[j].source->parent; pop = pop->parent)
              pop->seqin_count[msa_index]++;
          }
          curpop = simnodes[i]->mi->me[j].target;
        }
      }
      for (pop = curpop->parent; pop != NULL; pop = pop->parent)
        pop->seqin_count[msa_index]++;
          continue;
    }
    
    snode_t * newpop = simnodes[i]->pop;
    if (opt_migration && simnodes[i]->mi && simnodes[i]->mi->count)
      newpop = simnodes[i]->mi->me[simnodes[i]->mi->count-1].target;

    /* find ancestral population at time tnew */
    while (newpop->parent && newpop->parent->tau < tnew)
      newpop = newpop->parent;
    
    /* get lineages entering that population */
    long lineages = target_branches(stree,gtree,newpop,simnodes[i],tnew,gtree->travbuffer);

    /* randomly pick one lineage */
    j = (long)(lineages*legacy_rndu(thread_index));
    assert(j < lineages);
    gnode_t * gtarget = gtree->travbuffer[j];

    /* regraft */
    subtree_regraft(stree, gtree, simnodes[i], parents[p], gtarget, newpop, x,existing_migs_count,tnew);

    /* check whether we regrafted on the parental edge of a node that is to be simulated.
       this can only happen for lineages that enter X via migration, i.e. lineages that
       were pruned only partially. If that's the case, replace the simulated node with its parent */
    if (gtarget->mark & FLAG_PRUNED && gtarget->mark & FLAG_MIGRATE)
    {
      gtarget->mark &= ~FLAG_PRUNED;
      gtarget->mark &= ~FLAG_MIGRATE;
      gtarget->parent->mark |= FLAG_PRUNED;
      gtarget->parent->mark |= FLAG_MIGRATE;

      for (j = i+1; j < count; ++j)
        if (simnodes[j] == gtarget)
        {
          simnodes[j] = gtarget->parent;
          break;
        }
    }

    /* mark the nodes where we need to update the parental branch lengths */
    parents[p]->left->mark |= FLAG_BRANCH_UPDATE;
    parents[p]->right->mark |= FLAG_BRANCH_UPDATE;
    if (parents[p]->parent)
      parents[p]->mark |= FLAG_BRANCH_UPDATE;
    ++p;

    /* delete flag */
    if (simnodes[i]->mark & FLAG_PRUNED)
      simnodes[i]->mark &= ~FLAG_PRUNED;
    if (simnodes[i]->mark & FLAG_MIGRATE)
      simnodes[i]->mark &= ~FLAG_MIGRATE;
  }

  /* replace correct root pmatrix */
  assert(gtree->root->pmatrix_index == root_pmat_index);
  if (gtree->root->parent)
  {
    /* find root */
    gnode_t * newroot = gtree->root;
    while (newroot->parent)
      newroot = newroot->parent;
    
    gtree->root->pmatrix_index = newroot->pmatrix_index;
    newroot->pmatrix_index = root_pmat_index;
    gtree->root = newroot;
  }

  /* unset potential branch length update flag for root */
  if (gtree->root->mark & FLAG_BRANCH_UPDATE)
  {
    assert(!gtree->root->parent);
    gtree->root->mark &= ~FLAG_BRANCH_UPDATE;
  }

  /* 1) place nodes with branches needing update to an array, and 2) mark all
     nodes that require CLV updates */
  for (k=0, i=0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * p = gtree->nodes[i];

    /* 1. place nodes for which we will update branch lengths into travbuffer */
    if (p->mark & FLAG_BRANCH_UPDATE)
    {
      assert(p->parent);
      gtree->travbuffer[k++] = p;
    }

    /* 2. mark nodes on the path p->parent <-> root for CLV update */
    /* TODO: I think it's sufficient to only set FLAG_PARTIAL_UPDATE only for
       q->parent as p->parent <-> root is visited in gtree_partials_return */
    gnode_t * q;
    for (q = p->parent; q && !(q->mark & FLAG_PARTIAL_UPDATE); q = q->parent)
      q->mark |= FLAG_PARTIAL_UPDATE;
  }
      
  #if 1
  /* swap pmatrix buffers */
  for (j = 0; j < k; ++j)
    gtree->travbuffer[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,
                                                          gtree->travbuffer[j]->pmatrix_index);
  locus_update_matrices(locus,gtree,gtree->travbuffer,stree,msa_index,k);

  /* find CLVs that must be updated */
  gnode_t ** partials = gtree->travbuffer;
  gtree_return_partials(gtree->root, partials, &k);
  for (i = 0; i < k; ++i)
  {
    partials[i]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                            partials[i]->clv_index);
    if (opt_scaling)
      partials[i]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                    partials[i]->scaler_index);
  }

  /* update partials */
  locus_update_partials(locus,gtree->travbuffer,k);

  /* compute log-l and log-pr */
  logl = locus_root_loglikelihood(locus,gtree->root,locus->param_indices,NULL);
  #else
  k=0;
  for (j = 0; j < gtree->tip_count+gtree->inner_count; ++j)
  {
    gnode_t * tmp = gtree->nodes[j];
    if (tmp->parent)
    {
      tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,tmp->pmatrix_index);
      gtree->travbuffer[k++] = tmp;
    }
  }

  locus_update_matrices(locus,gtree,gtree->travbuffer,stree,msa_index,k);

  gtree_all_partials(gtree->root,gtree->travbuffer,&k);
  for (j = 0; j < k; ++j)
  {
    gtree->travbuffer[j]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,gtree->travbuffer[j]->clv_index);
    if (opt_scaling)
      gtree->travbuffer[j]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,gtree->travbuffer[j]->scaler_index);
  }
  locus_update_partials(locus,gtree->travbuffer,k);

  logl = locus_root_loglikelihood(locus,gtree->root,locus->param_indices,NULL);
  #endif

  /* update logpr */
  logpr = (opt_est_theta) ? gtree->logpr : stree->notheta_logpr;

  if (opt_migration && !opt_est_theta)
    fatal("Integrating out thetas not yet implemented for IM model");

  assert(opt_migration);
  if (opt_migration)
  {
    /* TODO: Improve the following by updating only the necessary popultions */
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    {
      snode_t * pop = stree->nodes[j];
      if (opt_est_theta)
        logpr -= pop->logpr_contrib[msa_index];
      else
        logpr -= pop->notheta_logpr_contrib;

      logpr += gtree_update_logprob_contrib_mig(pop,
                                                stree,
                                                gtree,
                                                locus->heredity[0],
                                                msa_index,
                                                thread_index);
    }
  }
  #if 0
  else
  {
    for (pop = start; pop != end; pop = pop->parent)
    {
      if (opt_est_theta)
        logpr -= pop->logpr_contrib[msa_index];
      else
        logpr -= pop->notheta_logpr_contrib;

      logpr += gtree_update_logprob_contrib(pop,
                                            locus->heredity[0],
                                            msa_index,
                                            thread_index);
    }
  }
  #endif

  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    gtree->nodes[i]->mark = 0;

  *ret_logl_diff  = logl - gtree->logl;
  *ret_logpr_diff = logpr - gtree->logpr;

  gtree->logpr = logpr;
  gtree->logl  = logl;
  
  free(migration);
  free(simnodes);
  free(parents);
}

static long stree_migration_append(gtree_t *** gtreeptr,
                                   gtree_t *** gcloneptr,
                                   stree_t ** streeptr,
                                   stree_t ** scloneptr,
                                   locus_t ** locus,
                                   long si,
                                   long ti)
{
  long i;
  long rc = 0;
  double logl_diff = 0;
  double logpr_diff = 0;
  double lli_diff = 0;
  double lpi_diff = 0;
  double lnacceptance;
  snode_t * x;
  snode_t * y;

  gtree_t ** gtree = *gcloneptr;
  stree_t * stree = *scloneptr;

  stree_t * original_stree = *streeptr;
  gtree_t ** original_gtree = *gtreeptr;

  const long thread_index_zero = 0;

  x = original_stree->nodes[si];
  y = original_stree->nodes[ti];
  assert(x->parent && y->parent);

  /* identify lower and upper bounds tau_L and tau_U for the time interval
     during which migration is possible */
  double tau_L = MAX(x->tau, y->tau);
  double tau_U = MIN(x->parent->tau, y->parent->tau);

  /*** Ziheng $$$ ***/
  dbg_counter++;

  #if 0
  printf("x: %s    y: %s\n", x->label, y->label);
  printf("%s->parent: %s\n", x->label, x->parent->label);
  printf("%s->parent: %s\n", y->label, y->parent->label);
  printf("%s->tau: %f  %s->parent->tau: %f   %s->tau: %f  %s->parent->tau: %f\n", x->label, x->tau, x->label, x->parent->tau, y->label, y->tau, y->label, y->parent->tau);
  printf("tauL : %f  tauU: %f\n", tau_L, tau_U);
  #endif
  assert(tau_L < tau_U);

  stree_clone(original_stree,stree);
  for (i = 0; i < stree->locus_count; ++i)
    gtree_clone(original_gtree[i], gtree[i], stree);
  events_clone(original_stree, stree, gtree);

  x = stree->nodes[si];
  y = stree->nodes[ti];
  #if 0
  printf("\nAppend M_%s->%s \n",x->label, y->label);
  #endif
  /* create a record of the migration and increase opt_migration */
  migspec_append(stree,x,y);

  /* TODO: The variable y->mb_mrsum_isarray is unnecessary and wrong. We have an
    ->active_count in migbuffer which should be updated according to the relevant migration.
    As we may have multiple migrations from/to one node, the mb_mrsum_isarray is insufficient */
  /* no rate variation */
  y->mb_mrsum_isarray = 0;

  /* find all lineages entering population */
  #if 1
  for (i = 0; i < opt_locus_count; ++i)
  {
    mig_dissolve_and_sim(stree,
                         gtree[i],
                         locus[i],
                         tau_L,
                         tau_U,
                         y,
                         i,
                         thread_index_zero,
                         &lli_diff,
                         &lpi_diff);
    


    logl_diff  += lli_diff;
    logpr_diff += lpi_diff;
    assert(!gtree[i]->root->parent);
  }
  #else
  logl_diff = 0; logpr_diff = 0;
  #endif

  lnacceptance = logl_diff;
  /* lnacceptance = logl_diff + logpr_diff; */
  /*** Ziheng $$$ ***/

  #if 0
  lnacceptance += logprobM; /* log prior ratio */
  lnacceptance -= logprobM; /* log proposal ratio */
  #else
  migspec_t * spec = opt_mig_specs + opt_migration_count-1;
  lnacceptance += log_pdfgamma(spec->M, spec->alpha, spec->beta);   /* log prior ration */
  lnacceptance -= log_pdfgamma(spec->M, dbg_prop_a,  dbg_prop_b);   /* log proposal ratio */
  #endif

  /* TODO: This is constant for fixed phylogenies */
  long mrcount = stree->tip_count*(stree->tip_count+1)*(stree->tip_count-1)/3;

  /* proposal ratio: coin toss, possible migration rates */
  double term = (.5/.5) * (mrcount - opt_migration_count+1) / opt_migration_count;
  //printf("*APPEND* term: %f  %ld  %ld\n", term, mrcount, opt_migration_count);
  assert(term > 0);

  lnacceptance += log(term);

  if (debug_rjmcmc)
    printf("lnaccept(append): %f %f %f\n", logl_diff, logpr_diff, lnacceptance);

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
  {
    if (debug_rjmcmc)
      printf( ANSI_COLOR_RED " accept" ANSI_COLOR_RESET);
    /* accept */

    SWAP(original_stree,stree);
    SWAP(original_gtree,gtree);
    rc = 1;
    #if 0
    printf("\nAppend M_%s->%s " ANSI_COLOR_GREEN "SUCCESS\n" ANSI_COLOR_RESET, x->label, y->label);
    #endif
  }
  else
  {
    /* reject */
    #if 0
    printf("\nAppend M_%s->%s " ANSI_COLOR_RED "FAIL\n" ANSI_COLOR_RESET, x->label, y->label);
    #endif
    --opt_migration_count;

    migspec_t * spec = opt_mig_specs+opt_migration_count;
    opt_mig_bitmatrix[spec->si][spec->ti] = 0;
    opt_migration_matrix[spec->si][spec->ti] = -1;
  }

  *streeptr = original_stree;
  *gtreeptr = original_gtree;

  *scloneptr = stree;
  *gcloneptr = gtree;

  return rc;
}

static long stree_migration_remove(gtree_t *** gtreeptr,
                                   gtree_t *** gcloneptr,
                                   stree_t ** streeptr,
                                   stree_t ** scloneptr,
                                   locus_t ** locus,
                                   long si,
                                   long ti)
{
  long i;
  long pos = 0;
  long rc = 0;
  double logl_diff = 0;
  double logpr_diff = 0;
  double lli_diff = 0;
  double lpi_diff = 0;
  double lnacceptance;
  snode_t * x;
  snode_t * y;

  gtree_t ** gtree = *gcloneptr;
  stree_t * stree = *scloneptr;

  stree_t * original_stree = *streeptr;
  gtree_t ** original_gtree = *gtreeptr;

  const long thread_index_zero = 0;

  x = original_stree->nodes[si];
  y = original_stree->nodes[ti];
  assert(x->parent && y->parent);

  /* see comment below */
  assert(y->mb_mrsum_isarray == 0);

  /* identify lower and upper bounds tau_L and tau_U for the time interval
     during which migration is possible */
  double tau_L = MAX(x->tau, y->tau);
  double tau_U = MIN(x->parent->tau, y->parent->tau);

  dbg_counter++;
  assert(tau_L < tau_U);

  stree_clone(original_stree,stree);
  for (i = 0; i < stree->locus_count; ++i)
    gtree_clone(original_gtree[i], gtree[i], stree);
  events_clone(original_stree, stree, gtree);

  x = stree->nodes[si];
  y = stree->nodes[ti];

  /* create a record of the migration and decrease opt_migration */
  pos = migspec_remove(stree,x,y);

  #if (DBG_TF)
  debug_print_migmatrix(stree);
  #endif

  /* TODO: The variable y->mb_mrsum_isarray is unnecessary and wrong. We have an
    ->active_count in migbuffer which should be updated according to the relevant migration.
    As we may have multiple migrations from/to one node, the mb_mrsum_isarray is insufficient */
  /* no rate variation */
  y->mb_mrsum_isarray = 0;

  #if 1
  /* find all lineages entering population */
  for (i = 0; i < opt_locus_count; ++i)
  {
    mig_dissolve_and_sim(stree,
                         gtree[i],
                         locus[i],
                         tau_L,
                         tau_U,
                         y,
                         i,
                         thread_index_zero,
                         &lli_diff,
                         &lpi_diff);
    
    logl_diff  += lli_diff;
    logpr_diff += lpi_diff;
  }
  #else
  logl_diff = 0; logpr_diff = 0;
  #endif

  /* gene tree is simulated and thus proposal density and gene tree density are cancelled */
  lnacceptance = logl_diff;


  /* TODO: Prior ratio and proposal ratio
  lnacceptance += logprobM
  lnacceptance -= logprobM
  */
  #if 1
  migspec_t * spec = opt_mig_specs+opt_migration_count;

  lnacceptance += log_pdfgamma(spec->M, dbg_prop_a,  dbg_prop_b);  /* proposal ratio */
  lnacceptance -= log_pdfgamma(spec->M, spec->alpha, spec->beta);  /* prior ratio */
  #endif


  /* TODO: This is constant for fixed phylogenies */
  long mrcount = stree->tip_count*(stree->tip_count+1)*(stree->tip_count-1)/3;

  /* proposal ratio: coin toss, possible migration rates */
  double term = (.5/.5) * (opt_migration_count+1) / (mrcount-opt_migration_count);
  assert(term>0);
  //printf("term: %f\n", term);
  lnacceptance += log(term);

  /* TODO: We need to add the prior and proposal ratio as in _append() function */


  if (debug_rjmcmc)
    printf("lnaccept(remove): %f %f %f\n", logl_diff, logpr_diff, lnacceptance);

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
  {
    if (debug_rjmcmc)
      printf( ANSI_COLOR_RED " accept" ANSI_COLOR_RESET);
    /* accept */
    #if 0
    printf("\nRemove M_%s->%s " ANSI_COLOR_GREEN "SUCCESS\n" ANSI_COLOR_RESET, x->label, y->label);
    #endif

    SWAP(original_stree,stree);
    SWAP(original_gtree,gtree);
    rc = 1;

    #if (DEBUG_zy)
    printf("[DEBUG] stree_migration_remove() -- accepted\n");
    printf("[DEBUG] -- stree_migration_remove -- Listing species tree and gene trees\n");
    debug_print_stree(stree);
    debug_print_gtree(original_gtree[0]);
    //debug_print_gtree(original_gtree[1]);
    printf("[DEBUG] -- stree_migration_remove -- END\n");
    #endif
  }
  else
  {
    /* reject */

    #if 0
    printf("\nRemove M_%s->%s " ANSI_COLOR_RED "FAIL\n" ANSI_COLOR_RESET, x->label, y->label);
    #endif
    ++opt_migration_count;
    opt_migration_matrix[opt_mig_specs[opt_migration_count-1].si][opt_mig_specs[opt_migration_count-1].ti] = opt_migration_count-1;
    if (pos != opt_migration_count-1)
    {
      SWAP(opt_mig_specs[pos],opt_mig_specs[opt_migration_count-1]);

      unsigned int si1=opt_mig_specs[pos].si;
      unsigned int ti1=opt_mig_specs[pos].ti;
      unsigned int si2=opt_mig_specs[opt_migration_count-1].si;
      unsigned int ti2=opt_mig_specs[opt_migration_count-1].ti;
      SWAP(opt_migration_matrix[si1][ti1],opt_migration_matrix[si2][ti2]);
      #if 0
      SWAP(opt_mig_bitmatrix[si1][ti1],opt_mig_bitmatrix[si2][ti2]);
      #endif
    }
    migspec_t * spec = opt_mig_specs+pos;
    
    opt_mig_bitmatrix[spec->si][spec->ti] = 1;

    #if (DEBUG_zy)
    printf("[DEBUG] stree_migration_remove() -- rejected\n");
    printf("[DEBUG] -- stree_migration_remove -- Listing species tree and gene trees\n");
    debug_print_stree(stree);
    debug_print_gtree(original_gtree[0]);
    if(opt_msa_count>1) debug_print_gtree(original_gtree[1]);
    printf("[DEBUG] -- stree_migration_remove -- END\n");
    #endif
  }
  #if 0
  printf("---------------------------------------------\n");
  #endif

  *streeptr = original_stree;
  *gtreeptr = original_gtree;

  *scloneptr = stree;
  *gcloneptr = gtree;

  return rc;
}

void migdissolveandism2(stree_t * stree,
                        gtree_t * gtree,
                        locus_t * locus,
                        double tL,
                        double tU,
                        snode_t * x,
                        snode_t * y,
                        int msa_index,
                        long thread_index,
                        double * ret_logl_diff,
                        double * ret_logpr_diff)
{
  unsigned int i,j,k,count,p;
  unsigned int root_pmat_index;
  unsigned int total_nodes = gtree->tip_count+gtree->inner_count;
  double prevt;
  double tnew;
  double logl, logpr;
  gnode_t ** simnodes = NULL;
  gnode_t ** parents = NULL;
  long * migration = NULL;
  long * population = NULL;
  snode_t * prevpop = NULL;
  snode_t ** affected_nodemap = NULL;
  static long debug_prune_counter = 0;

  simnodes   = (gnode_t **)xcalloc((size_t)total_nodes, sizeof(gnode_t *));
  migration  = (long *)xcalloc((size_t)total_nodes, sizeof(long));
  population = (long *)xcalloc((size_t)total_nodes, sizeof(long));
  affected_nodemap = (snode_t **)xcalloc((size_t)total_nodes, sizeof(snode_t *));

  long pop_x = 1;
  long pop_y = 2;
  long pop_sim;

  /* store pmatrix index for root */
  root_pmat_index = gtree->root->pmatrix_index;

  /* 1. find lineages to dissolve */
  for (i = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->nodes[i];

    /* if any of the below conditions hold, the lineage is  not to be deleted */
    if (gnode->time > tU) continue;
    if ((gnode->pop == x || gnode->pop == y) && gnode->time > tL) continue;
    if (gnode->parent && gnode->parent->time < tL) continue;

    /* check if the lineage passes by x or y */
    prevt = gnode->time;
    prevpop = gnode->pop;
    if (opt_migration && gnode->mi && gnode->mi->count)
    {
      for (j = 0; j < gnode->mi->count; ++j)
      {
        if (gnode->mi->me[j].time > tL && prevt < tU &&
            (stree->pptable[prevpop->node_index][x->node_index] || stree->pptable[prevpop->node_index][y->node_index]))
          break;

        prevt = gnode->mi->me[j].time;
        prevpop = gnode->mi->me[j].target;
      }
    }
    if (!(stree->pptable[prevpop->node_index][x->node_index] || stree->pptable[prevpop->node_index][y->node_index]) ||
        prevt > tU) continue;

    #if 1
    /* assertion that lineage does not pass both x and y (that's impossible) */
    assert(!(stree->pptable[prevpop->node_index][x->node_index] &&
             stree->pptable[prevpop->node_index][y->node_index]));

    if (prevpop == x)
      pop_sim = pop_x;
    else if (prevpop == y)
      pop_sim = pop_y;
    else if (stree->pptable[prevpop->node_index][x->node_index])
      pop_sim = pop_x;
    else
      pop_sim = pop_y;
    #endif

    if (gnode->parent)
    {
      gnode->parent->mark |= (gnode->parent->left == gnode) ? 
                               FLAG_RED_LEFT : FLAG_RED_RIGHT;
      
      /* unset FLAG_SIMULATE on parent if both its daughter lineages are red */
      if ((gnode->parent->mark & FLAG_SIMULATE) &&
          (gnode->parent->mark & FLAG_RED_LEFT) &&
          (gnode->parent->mark & FLAG_RED_RIGHT))
      {
        gnode->parent->mark &= ~FLAG_SIMULATE;
        if (gnode->parent->mark & FLAG_MIGRATE)
          gnode->parent->mark &= ~FLAG_MIGRATE;

        population[gnode->parent->node_index] = 0;
      }
    }

    /* set SIMULATE flag only if at most one RED daughter lineage */
    if (!((gnode->mark & FLAG_RED_LEFT) && (gnode->mark & FLAG_RED_RIGHT)))
    {
      gnode->mark |= FLAG_SIMULATE;

      if (prevt > tL)
        gnode->mark |= FLAG_MIGRATE;

      population[gnode->node_index] = pop_sim;
    }
  }

  /* color the lineages with both daughter lineages red */ 
  gtree_traverse(gtree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_postorder_inner_only,
                 gtree->travbuffer,
                 &total_nodes);
  for (i = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->travbuffer[i];

    if (gnode->left->mark & FLAG_RED_LEFT && gnode->left->mark & FLAG_RED_RIGHT)
      gnode->mark |= FLAG_RED_LEFT;

    if (gnode->right->mark & FLAG_RED_LEFT && gnode->right->mark & FLAG_RED_RIGHT)
      gnode->mark |= FLAG_RED_RIGHT;
  }

  /* 1b. Go through nodes with SIMULATE flag in preorder */
  gtree_traverse(gtree->root,
                 TREE_TRAVERSE_PREORDER,
                 cb_preorder_inner_only,
                 gtree->travbuffer,
                 &total_nodes);
  for (i = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->travbuffer[i];

    //assert(gnode->mark & FLAG_SIMULATE);

    if ((gnode->mark & FLAG_RED_LEFT) &&
        (gnode->mark & FLAG_RED_RIGHT))
    {
      /* unset FLAG_SIMULATE if both daughter lineages are red */
      if (gnode->mark & FLAG_SIMULATE)
      {
        gnode->mark &= ~FLAG_SIMULATE;
        #if 1
        assert(population[gnode->node_index]);
        population[gnode->node_index] = 0;
        #endif
      }
      else
        assert(!population[gnode->node_index]);
      if (gnode->mark & FLAG_MIGRATE)
        gnode->mark &= ~FLAG_MIGRATE;

    }
    else if (gnode->mark & FLAG_RED_LEFT)
    {
      /* move FLAG_SIMULATE from gnode to right daughter lineages */
      if (gnode->mark & FLAG_SIMULATE)
      {
        #if 1
        assert(population[gnode->node_index]);
        if (gnode->right->mark & FLAG_SIMULATE)
          assert(population[gnode->node_index] == population[gnode->right->node_index]);
        population[gnode->right->node_index] = population[gnode->node_index];
        population[gnode->node_index] = 0;
        #endif
        gnode->mark &= ~FLAG_SIMULATE;
        gnode->right->mark |= FLAG_SIMULATE;
      }
      if (gnode->mark & FLAG_MIGRATE)
      {
        gnode->mark &= ~FLAG_MIGRATE;
        gnode->right->mark |= FLAG_MIGRATE;
      }
    }
    else if (gnode->mark & FLAG_RED_RIGHT)
    {
      /* move FLAG_SIMULATE from gnode to left daughter lineages */
      if (gnode->mark & FLAG_SIMULATE)
      {
        #if 1
        assert(population[gnode->node_index]);
        if (gnode->left->mark & FLAG_SIMULATE)
          assert(population[gnode->node_index] == population[gnode->left->node_index]);
        population[gnode->left->node_index] = population[gnode->node_index];
        population[gnode->node_index] = 0;
        #endif
        gnode->mark &= ~FLAG_SIMULATE;
        gnode->left->mark |= FLAG_SIMULATE;
      }
      if (gnode->mark & FLAG_MIGRATE)
      {
        gnode->mark &= ~FLAG_MIGRATE;
        gnode->left->mark |= FLAG_MIGRATE;
      }
    }
  }

  /* 1c. create the simnodes/migration arrays */
  total_nodes = gtree->tip_count + gtree->inner_count;
  for (i = 0, count = 0; i < total_nodes; ++i)
  {
    gnode_t * gnode = gtree->nodes[i];

    if (gnode->mark & FLAG_SIMULATE)
    {
      if (gnode->mark & FLAG_MIGRATE)
        migration[count] = 1;
      simnodes[count++] = gnode;
      #if 1
      assert(population[gnode->node_index]);
      #endif
    }

    /* keep red marks */
    gnode->mark &= (FLAG_RED_LEFT | FLAG_RED_RIGHT);
  }

  parents  = (gnode_t **)xcalloc((size_t)count, sizeof(gnode_t *));

  /* 2. Prune the lineages, and store the deleted nodes */
  for (i = 0, p = 0; i < count; ++i)
  {
    #if 0
    printf("\n\nBefore pruning:\n\n");
    debug_print_stree(stree);
    debug_print_gtree_detailed(gtree);
    #endif

    debug_prune_counter++;
    #if 0
    printf("[" ANSI_COLOR_MAGENTA "debug_prune_counter=%ld" ANSI_COLOR_RESET "] Pruning %d/%d (node: %d)\n", debug_prune_counter, i,count,simnodes[i]->node_index);
    if (population[simnodes[i]->node_index] != pop_x && population[simnodes[i]->node_index] != pop_y)
      printf("Assertion failing.. population[%d] = %ld\n", simnodes[i]->node_index, population[simnodes[i]->node_index]);
    #endif
    assert(population[simnodes[i]->node_index] == pop_x || population[simnodes[i]->node_index] == pop_y);
    gnode_t * deleted = pruneoff2(stree,gtree,population[simnodes[i]->node_index]==pop_x ? x : y,simnodes[i],tL);
    if (deleted)
    {
      parents[p++] = deleted;
      deleted->old_pop = deleted->pop;
      deleted->pop = NULL;
    }

    simnodes[i]->mark |= FLAG_PRUNED;
    if  (migration[i])
      simnodes[i]->mark |= FLAG_MIGRATE;

    /* mark the current node for branch length update */
    simnodes[i]->mark |= FLAG_BRANCH_UPDATE;
    
  }
  #if 0
  printf("\n\nPruned off tree:\n\n");
  debug_print_stree(stree);
  debug_print_gtree_detailed(gtree);
  #endif

  /* unset color flags from nodes */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    gtree->nodes[i]->mark &= ~(FLAG_RED_LEFT | FLAG_RED_RIGHT);

  /* create affected nodemap */
  for (i = 0; i < count; ++i)
  {
    assert(population[simnodes[i]->node_index] == pop_x || population[simnodes[i]->node_index] == pop_y);
    if (population[simnodes[i]->node_index] == pop_x)
      affected_nodemap[simnodes[i]->node_index] = x; 
    else
      affected_nodemap[simnodes[i]->node_index] = y; 
  }

  /* 3. re-simulate */
  for (p = 0, i = 0; i < count; ++i)
  {
    long existing_migs_count = simnodes[i]->mi ? simnodes[i]->mi->count : 0;
    #if 1
    assert(population[simnodes[i]->node_index] == pop_x || population[simnodes[i]->node_index] == pop_y);
    #endif
    tnew = simulate_coalescent_mig_from_point(stree,
                                              gtree,
                                              simnodes[i],
                                              population[simnodes[i]->node_index] == pop_x ? x : y,
                                              population[simnodes[i]->node_index] == pop_x ? y : x,
                                              affected_nodemap,
                                              simnodes[i]->mark & FLAG_MIGRATE ? simnodes[i]->mi->me[existing_migs_count-1].time : tL,               /* TODO: Check what we are passing here is correct */
                                              msa_index,
                                              thread_index);

    if (tnew < 0)
    {
      /* if this is the lineage going back to infinitity, i.e. first lineage
         we simulate that does not coalesce anywhere */

      /* delete flag */
      if (simnodes[i]->mark & FLAG_PRUNED)
        simnodes[i]->mark &= ~FLAG_PRUNED;
      if (simnodes[i]->mark & FLAG_MIGRATE)
        simnodes[i]->mark &= ~FLAG_MIGRATE;

      /* TODO: Update seqin counts */
      snode_t * pop;
      snode_t * curpop = population[simnodes[i]->node_index] == pop_x ? x : y;
      if (opt_migration && simnodes[i]->mi && simnodes[i]->mi->count)
      {
        for (j = existing_migs_count; j < simnodes[i]->mi->count; ++j)
        {
          simnodes[i]->mi->me[j].source->migevent_count[msa_index]++;
          simnodes[i]->mi->me[j].target->migevent_count[msa_index]++;

          /* increase # of migration from s to t (forward in time) */
          long s = simnodes[i]->mi->me[j].target->node_index;
          long t = simnodes[i]->mi->me[j].source->node_index;
          gtree->migcount[s][t]++;

          if (simnodes[i]->mi->me[j].source != curpop)
          {
            for (pop = curpop->parent; pop != simnodes[i]->mi->me[j].source->parent; pop = pop->parent)
              pop->seqin_count[msa_index]++;
          }
          curpop = simnodes[i]->mi->me[j].target;
        }
      }
      for (pop = curpop->parent; pop != NULL; pop = pop->parent)
        pop->seqin_count[msa_index]++;
          continue;
    }
    
    snode_t * newpop = simnodes[i]->pop;
    if (opt_migration && simnodes[i]->mi && simnodes[i]->mi->count)
      newpop = simnodes[i]->mi->me[simnodes[i]->mi->count-1].target;

    /* find ancestral population at time tnew */
    while (newpop->parent && newpop->parent->tau < tnew)
      newpop = newpop->parent;
    
    /* get lineages entering that population */
    long lineages = target_branches(stree,gtree,newpop,simnodes[i],tnew,gtree->travbuffer);

    /* randomly pick one lineage */
    j = (long)(lineages * legacy_rndu(thread_index));
    assert(j < lineages);
    gnode_t * gtarget = gtree->travbuffer[j];

    /* regraft */
    subtree_regraft(stree,
                    gtree,
                    simnodes[i],
                    parents[p],
                    gtarget,
                    newpop,
                    population[simnodes[i]->node_index] == pop_x ? x : y,
                    existing_migs_count,
                    tnew);

    /* check whether we regrafted on the parental edge of a node that is to be simulated.
       this can only happen for lineages that enter X via migration, i.e. lineages that
       were pruned only partially. If that's the case, replace the simulated node with its parent */
    if (gtarget->mark & FLAG_PRUNED && gtarget->mark & FLAG_MIGRATE)
    {
      gtarget->mark &= ~FLAG_PRUNED;
      gtarget->mark &= ~FLAG_MIGRATE;
      gtarget->parent->mark |= FLAG_PRUNED;
      gtarget->parent->mark |= FLAG_MIGRATE;

      for (j = i+1; j < count; ++j)
        if (simnodes[j] == gtarget)
        {
          population[gtarget->parent->node_index] = population[simnodes[j]->node_index];

          assert(population[simnodes[j]->node_index] == pop_x || population[simnodes[j]->node_index] == pop_y);
          assert(affected_nodemap[simnodes[j]->node_index] == x || affected_nodemap[simnodes[j]->node_index] == y);

          affected_nodemap[gtarget->parent->node_index] = affected_nodemap[simnodes[j]->node_index];
          simnodes[j] = gtarget->parent;

          break;
        }
    }

    /* mark the nodes where we need to update the parental branch lengths */
    parents[p]->left->mark |= FLAG_BRANCH_UPDATE;
    parents[p]->right->mark |= FLAG_BRANCH_UPDATE;
    if (parents[p]->parent)
      parents[p]->mark |= FLAG_BRANCH_UPDATE;
    ++p;

    /* delete flag */
    if (simnodes[i]->mark & FLAG_PRUNED)
      simnodes[i]->mark &= ~FLAG_PRUNED;
    if (simnodes[i]->mark & FLAG_MIGRATE)
      simnodes[i]->mark &= ~FLAG_MIGRATE;
  }

  /* replace correct root pmatrix */
  assert(gtree->root->pmatrix_index == root_pmat_index);
  if (gtree->root->parent)
  {
    /* find root */
    gnode_t * newroot = gtree->root;
    while (newroot->parent)
      newroot = newroot->parent;
    
    gtree->root->pmatrix_index = newroot->pmatrix_index;
    newroot->pmatrix_index = root_pmat_index;
    gtree->root = newroot;
  }

  /* unset potential branch length update flag for root */
  if (gtree->root->mark & FLAG_BRANCH_UPDATE)
  {
    assert(!gtree->root->parent);
    gtree->root->mark &= ~FLAG_BRANCH_UPDATE;
  }

  /* 1) place nodes with branches needing update to an array, and 2) mark all
     nodes that require CLV updates */
  for (k=0, i=0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * p = gtree->nodes[i];

    /* 1. place nodes for which we will update branch lengths into travbuffer */
    if (p->mark & FLAG_BRANCH_UPDATE)
    {
      assert(p->parent);
      gtree->travbuffer[k++] = p;
    }

    /* 2. mark nodes on the path p->parent <-> root for CLV update */
    /* TODO: I think it's sufficient to only set FLAG_PARTIAL_UPDATE only for
       q->parent as p->parent <-> root is visited in gtree_partials_return */
    gnode_t * q;
    for (q = p->parent; q && !(q->mark & FLAG_PARTIAL_UPDATE); q = q->parent)
      q->mark |= FLAG_PARTIAL_UPDATE;
  }
      
  #if 1
  /* swap pmatrix buffers */
  for (j = 0; j < k; ++j)
    gtree->travbuffer[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,
                                                          gtree->travbuffer[j]->pmatrix_index);
  locus_update_matrices(locus,gtree,gtree->travbuffer,stree,msa_index,k);

  /* find CLVs that must be updated */
  gnode_t ** partials = gtree->travbuffer;
  gtree_return_partials(gtree->root, partials, &k);
  for (i = 0; i < k; ++i)
  {
    partials[i]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                            partials[i]->clv_index);
    if (opt_scaling)
      partials[i]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                    partials[i]->scaler_index);
  }

  /* update partials */
  locus_update_partials(locus,gtree->travbuffer,k);

  /* compute log-l and log-pr */
  logl = locus_root_loglikelihood(locus,gtree->root,locus->param_indices,NULL);
  #else
  k=0;
  for (j = 0; j < gtree->tip_count+gtree->inner_count; ++j)
  {
    gnode_t * tmp = gtree->nodes[j];
    if (tmp->parent)
    {
      tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,tmp->pmatrix_index);
      gtree->travbuffer[k++] = tmp;
    }
  }

  locus_update_matrices(locus,gtree,gtree->travbuffer,stree,msa_index,k);

  gtree_all_partials(gtree->root,gtree->travbuffer,&k);
  for (j = 0; j < k; ++j)
  {
    gtree->travbuffer[j]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,gtree->travbuffer[j]->clv_index);
    if (opt_scaling)
      gtree->travbuffer[j]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,gtree->travbuffer[j]->scaler_index);
  }
  locus_update_partials(locus,gtree->travbuffer,k);

  logl = locus_root_loglikelihood(locus,gtree->root,locus->param_indices,NULL);
  #endif

  /* update logpr */
  logpr = (opt_est_theta) ? gtree->logpr : stree->notheta_logpr;

  if (opt_migration && !opt_est_theta)
    fatal("Integrating out thetas not yet implemented for IM model");

  assert(opt_migration);
  if (opt_migration)
  {
    /* TODO: Improve the following by updating only the necessary popultions */
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    {
      snode_t * pop = stree->nodes[j];
      if (opt_est_theta)
        logpr -= pop->logpr_contrib[msa_index];
      else
        logpr -= pop->notheta_logpr_contrib;

      logpr += gtree_update_logprob_contrib_mig(pop,
                                                stree,
                                                gtree,
                                                locus->heredity[0],
                                                msa_index,
                                                thread_index);
    }
  }
  #if 0
  else
  {
    for (pop = start; pop != end; pop = pop->parent)
    {
      if (opt_est_theta)
        logpr -= pop->logpr_contrib[msa_index];
      else
        logpr -= pop->notheta_logpr_contrib;

      logpr += gtree_update_logprob_contrib(pop,
                                            locus->heredity[0],
                                            msa_index,
                                            thread_index);
    }
  }
  #endif

  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    gtree->nodes[i]->mark = 0;

  *ret_logl_diff  = logl - gtree->logl;
  *ret_logpr_diff = logpr - gtree->logpr;

  gtree->logpr = logpr;
  gtree->logl  = logl;
  
  free(migration);
  free(simnodes);
  free(parents);
  free(population);
  free(affected_nodemap);


}

static long stree_migration_flip(gtree_t *** gtreeptr,
                                 gtree_t *** gcloneptr,
                                 stree_t ** streeptr,
                                 stree_t ** scloneptr,
                                 locus_t ** locus,
                                 long si,
                                 long ti)
{
  long i;
  long pos = 0;
  long rc = 0;
  double logl_diff = 0;
  double logpr_diff = 0;
  double lli_diff = 0;
  double lpi_diff = 0;
  double lnacceptance;
  snode_t * x;
  snode_t * y;

  gtree_t ** gtree = *gcloneptr;
  stree_t * stree = *scloneptr;

  stree_t * original_stree = *streeptr;
  gtree_t ** original_gtree = *gtreeptr;

  const long thread_index_zero = 0;

  x = original_stree->nodes[si];
  y = original_stree->nodes[ti];
  assert(x->parent && y->parent);

  /* see comment below */
  assert(y->mb_mrsum_isarray == 0);

  /* identify lower and upper bounds tau_L and tau_U for the time interval
     during which migration is possible */
  double tau_L = MAX(x->tau, y->tau);
  double tau_U = MIN(x->parent->tau, y->parent->tau);

  assert(tau_L < tau_U);

  stree_clone(original_stree,stree);
  for (i = 0; i < stree->locus_count; ++i)
    gtree_clone(original_gtree[i], gtree[i], stree);
  events_clone(original_stree, stree, gtree);

  x = stree->nodes[si];
  y = stree->nodes[ti];

  /* start of new code */
  pos = migspec_remove(stree,x,y);
  double old_mrate = opt_mig_specs[opt_migration_count].M;
  migspec_append(stree,y,x);
  opt_mig_specs[opt_migration_count-1].M = old_mrate;
  stree_update_mig_subpops(stree, thread_index_zero);

  /* TODO: The variable y->mb_mrsum_isarray is unnecessary and wrong. We have an
    ->active_count in migbuffer which should be updated according to the relevant migration.
    As we may have multiple migrations from/to one node, the mb_mrsum_isarray is insufficient */
  /* no rate variation */
  y->mb_mrsum_isarray = 0;

  /* find all lineages entering population */
  x->mb_mrsum_isarray = 0;
  for (i = 0; i < opt_locus_count; ++i)
  {
    migdissolveandism2(stree,
                       gtree[i],
                       locus[i],
                       tau_L,
                       tau_U,
                       y,
                       x,
                       i,
                       thread_index_zero,
                       &lli_diff,
                       &lpi_diff);
    logl_diff += lli_diff;
    logpr_diff += lpi_diff;
    assert(!gtree[i]->root->parent);
  }

  lnacceptance = logl_diff;
  /* lnacceptance = logl_diff + logpr_diff; */


  /* TODO: Prior ratio and proposal ratio
  lnacceptance += logprobM
  lnacceptance -= logprobM
  */

  /* end of new code */


  #if 1
  migspec_t * spec = opt_mig_specs+opt_migration_count-1;

  lnacceptance += log_pdfgamma(spec->M, dbg_prop_a,  dbg_prop_b);  /* proposal ratio */
  lnacceptance -= log_pdfgamma(spec->M, spec->alpha, spec->beta);  /* prior ratio */
  #endif

  if (debug_rjmcmc)
  {
    printf("lnacceptance(flip__): %f %f %f\n", logl_diff, logpr_diff, lnacceptance);
    if (dbg_counter > 10) assert(0);
  }
  //lnacceptance = -100;
  //if (lnacceptance >= -1e-10 || legacy_rndu(thread_index_zero) < exp(lnacceptance))
  if (lnacceptance >= -1e-10 || .5 < exp(lnacceptance))
  {
    if (debug_rjmcmc)
      printf( ANSI_COLOR_RED " accept" ANSI_COLOR_RESET);
    /* accept */
    #if 0
    printf("\nRemove M_%s->%s " ANSI_COLOR_GREEN "SUCCESS\n" ANSI_COLOR_RESET, x->label, y->label);
    #endif

    SWAP(original_stree,stree);
    SWAP(original_gtree,gtree);
    rc = 1;

    #if (DEBUG_zy)
    printf("[DEBUG] stree_migration_remove() -- accepted\n");
    printf("[DEBUG] -- stree_migration_remove -- Listing species tree and gene trees\n");
    debug_print_stree(stree);
    debug_print_gtree(original_gtree[0]);
    //debug_print_gtree(original_gtree[1]);
    printf("[DEBUG] -- stree_migration_remove -- END\n");
    #endif
  }
  else
  {
    /* reject */

    #if 0
    printf("\nRemove M_%s->%s " ANSI_COLOR_RED "FAIL\n" ANSI_COLOR_RESET, x->label, y->label);
    #endif
    migspec_t * spec = opt_mig_specs+opt_migration_count-1;
    opt_mig_bitmatrix[spec->si][spec->ti] = 0;
    //assert(opt_migration_matrix[spec->si][spec->ti] >= 0);
    //assert(opt_migration_matrix[spec->ti][spec->si] == -1);
    opt_migration_matrix[spec->ti][spec->si] = opt_migration_matrix[spec->si][spec->ti];
    opt_migration_matrix[spec->si][spec->ti] = -1;
    assert(opt_mig_bitmatrix[spec->ti][spec->si] == 0);
    SWAP(spec->si,spec->ti);

    if (pos != opt_migration_count-1)
    {
      SWAP(opt_mig_specs[pos],opt_mig_specs[opt_migration_count-1]);

      unsigned int si1=opt_mig_specs[pos].si;
      unsigned int ti1=opt_mig_specs[pos].ti;
      unsigned int si2=opt_mig_specs[opt_migration_count-1].si;
      unsigned int ti2=opt_mig_specs[opt_migration_count-1].ti;
      SWAP(opt_migration_matrix[si1][ti1],opt_migration_matrix[si2][ti2]);
      #if 0
      SWAP(opt_mig_bitmatrix[si1][ti1],opt_mig_bitmatrix[si2][ti2]);
      #endif
    }
    spec = opt_mig_specs+pos;
    
    opt_mig_bitmatrix[spec->si][spec->ti] = 1;

    #if (DEBUG_zy)
    printf("[DEBUG] stree_migration_remove() -- rejected\n");
    printf("[DEBUG] -- stree_migration_remove -- Listing species tree and gene trees\n");
    debug_print_stree(stree);
    debug_print_gtree(original_gtree[0]);
    if(opt_msa_count>1) debug_print_gtree(original_gtree[1]);
    printf("[DEBUG] -- stree_migration_remove -- END\n");
    #endif
  }
  #if 0
  printf("---------------------------------------------\n");
  #endif

  *streeptr = original_stree;
  *gtreeptr = original_gtree;

  *scloneptr = stree;
  *gcloneptr = gtree;

  return rc;
}

long stree_migration_flip_wrapper(gtree_t *** gtreeptr,
                                  gtree_t *** gcloneptr,
                                  stree_t ** streeptr,
                                  stree_t ** scloneptr,
                                  locus_t ** locus)
{
  long * indices;
  long i,j,k,r;
  const long thread_index_zero = 0;
  long rc = 0;
  migspec_t * spec;
  stree_t * stree = *streeptr;

  /* TODO: Duplicate code from _update(), delete */
  dbg_prop_a = 1.0 * opt_mig_alpha;
  dbg_prop_b = opt_mig_beta;

  long total = stree->tip_count+stree->inner_count;

  indices = (long *)xcalloc((size_t)(total*total), sizeof(long));

  /* go through the lower triangle of the matrix */
  k = 0;
  for (i = 0; i < total; ++i)
    for (j = 0; j < i; ++j)
    {
      if (opt_mig_bitmatrix[i][j] ^ opt_mig_bitmatrix[j][i])
        indices[k++] = opt_mig_bitmatrix[i][j] ?
          opt_migration_matrix[i][j] : opt_migration_matrix[j][i];
    }
  if (!k) goto l_unwind;

  r = (long)(k*legacy_rndu(thread_index_zero));
  spec = opt_mig_specs+indices[r];

  #if 0
  if (debug_rjmcmc)
  {
    printf("Flip__  ");
    debug_print_migrations_flip(stree);
    debug_print_migmatrix(stree);
    //debug_print_bitmatrix(stree);
  }
  printf(ANSI_COLOR_BLUE "Flipping M_%s -> %s\n" ANSI_COLOR_RESET, stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label);
  debug_print_stree(stree);
  debug_print_gtree(*gtreeptr[0]);
  debug_print_migmatrix(stree);
  #endif
  rc = stree_migration_flip(gtreeptr,
                            gcloneptr,
                            streeptr,
                            scloneptr,
                            locus,
                            spec->si,
                            spec->ti);
l_unwind:
  free(indices);
  
  #if 0
  if (debug_rjmcmc)
  {
    debug_print_migrations_flip(stree);
    debug_print_migmatrix(stree);
    //debug_print_bitmatrix(stree);
  }
  #endif
  return rc;
}


static long select_migpair(stree_t * stree, long * si, long * ti)
{
  long total_nodes = stree->tip_count + stree->inner_count;
  long total_pairs = 0;
  long i,j,k;
  long overlap;
  long rc = 0;
  long ** targets;
  long * tgt_count;

  const long thread_index_zero = 0;

  targets = (long **)xmalloc((size_t)total_nodes * sizeof(long *));
  for (i = 0; i < total_nodes; ++i)
    targets[i] = (long *)xmalloc((size_t)total_nodes * sizeof(long));

  tgt_count = (long *)xcalloc((size_t)total_nodes,sizeof(long));

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * s = stree->nodes[i];
    if (!s->parent) continue;

    double a = s->tau;
    double p = s->parent->tau;

    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
    {
      snode_t * t = stree->nodes[j];
      if (!t->parent || t == s) continue;

      /* check for time overlap */

      /* find populations with time overlap (as in propose_tau_mig) and no
         existing migration rate parameter */

      double b = t->tau;
      double q = t->parent->tau;

      overlap = !(p <= b || a >= q);

      if (overlap && !opt_mig_bitmatrix[s->node_index][t->node_index])
      {
        targets[s->node_index][tgt_count[s->node_index]++] = t->node_index;
        ++total_pairs;
      }
    }
  }

  /* if no pair found return with 0 */
  if (!total_pairs) goto l_unwind;

  /* randomly select a pair acrosss all possible pairs */
  long pair = (long)(total_pairs*legacy_rndu(thread_index_zero));

  for (i = 0, k = 0; i < total_nodes-1; ++i)
  {
    if (pair < k+tgt_count[i])
      break;
    k += tgt_count[i];
  }
  *si = i;
  *ti = targets[i][pair-k];
  rc = 1;

l_unwind:
  /* TODO: pre-allocate structure and deallocate them at the end of run */
  free(tgt_count);
  for (i = 0; i < total_nodes; ++i)
    free(targets[i]);
  free(targets);

  return rc;
}

static long lh_tree_counter(stree_t * stree)
{
  long i,j;
  long n,k;
  double prod = 1;
  long lv,rv;

  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * snode = stree->nodes[i];

    lv = snode->left->leaves == 1 ? 0 : snode->left->leaves-1;
    rv = snode->right->leaves == 1 ? 0 : snode->right->leaves-1;

    n = lv+rv;
    k = lv;

    if (n==k || k==0) continue;
    assert(n>k);

    
    double r = 1;
    for (j = n-k+1; j <= n; ++j)
      r *= j;
    for (j = 2; j <= k; ++j)
      r /= j;

    prod *= r;
  }

  return (long)prod;
}

static double lh_counter(long tip_count)
{
  long i;
  double lh = 1;

  for (i = 3; i <= tip_count; ++i)
    lh *= i*(i-1) / 2;

  return lh;
}

static long mig_models_count(long n)
{
  return n*(n+1)*(n-1)/3;
}


/*** Ziheng $$$ ***/
#define DEBUG_zy
int times = 0;

static void debug_print_migrations_flip(stree_t * stree)
{
  long i;
  printf("%ld: ", opt_migration_count);
  for (i = 0; i < opt_migration_count; ++i)
  {
    migspec_t * spec = opt_mig_specs+i;
    unsigned int si = spec->si;
    unsigned int ti = spec->ti;
    printf("  M_%s->%s=%f (%d,%d,%f,%f)", stree->nodes[si]->label, stree->nodes[ti]->label, spec->M, spec->si,spec->ti,spec->alpha,spec->beta);
  }
  printf("\n");
  //printf("  ");
}
static void debug_print_migrations(stree_t * stree)
{
  long i;
  printf("%ld: ", opt_migration_count);
  for (i = 0; i < opt_migration_count; ++i)
  {
    migspec_t * spec = opt_mig_specs+i;

    printf("  M_%s->%s=%f", stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label, spec->M);
  }
  printf("\n");
  //printf("  ");
}
static void debug_print_bitmatrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;

  printf("      ");
  for (i = 0; i < total_nodes; ++i)
    printf(" %2s", stree->nodes[i]->label);
  printf("\n");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("   ");
    printf(" %2s", stree->nodes[i]->label);
    for (j = 0; j < total_nodes; ++j)
    {
      if (opt_mig_bitmatrix[i][j] != 0)
        printf(ANSI_COLOR_RED " %2ld" ANSI_COLOR_RESET, opt_mig_bitmatrix[i][j]);
      else
        printf(" %2ld", opt_mig_bitmatrix[i][j]);
    }
    printf("\n");
  }
}
static void debug_print_migmatrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;

  printf("      ");
  for (i = 0; i < total_nodes; ++i)
    printf(" %2s", stree->nodes[i]->label);
  printf("\n");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("   ");
    printf(" %2s", stree->nodes[i]->label);
    for (j = 0; j < total_nodes; ++j)
    {
      if (opt_migration_matrix[i][j] != -1)
        printf(ANSI_COLOR_RED " %2ld" ANSI_COLOR_RESET, opt_migration_matrix[i][j]);
      else
        printf(" %2ld", opt_migration_matrix[i][j]);
    }
    printf("\n");
  }
}

static void debug_check_matrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;
  long count = 0;

  for (i = 0; i < total_nodes; ++i)
    for (j = 0; j < total_nodes; ++j)
      if (opt_migration_matrix[i][j] != -1)
        count++;

  assert(count == opt_migration_count);
}

long stree_migration_rj(gtree_t *** gtreeptr,
                            gtree_t *** gcloneptr,
                            stree_t ** streeptr,
                            stree_t ** scloneptr,
                            locus_t ** locus)
{
  long append = 0;
  long rc = 0;
  long si,ti;
  const long thread_index_zero = 0;
  migspec_t * spec;
  stree_t * stree = *streeptr;
  static long debug_counter = 0;

  dbg_prop_a = 1.0 * opt_mig_alpha;
  dbg_prop_b = opt_mig_beta;
  
  #if(DBG_TF)
  printf("HERE!\n");
  debug_print_migrations(stree);
  debug_print_migmatrix(stree);
  #endif
  debug_check_matrix(stree);
#if 1
  dbg_mig_idx = dbg_get_mig_idx(stree);
#endif

  append = legacy_rndu(thread_index_zero) < .5;
  if (!append && opt_migration_count == 0) return 0;

  ++debug_counter;
  if (append)
  {
    /* check if there is an available pair, and if yes, select it */
    rc = select_migpair(stree,&si,&ti);
    assert(rc || opt_migration_count);
    if (!rc)
      return 0;
  }

  if (debug_rjmcmc)
  {
    if (append)
      printf("Append  ");
    else
      printf("Remove  ");
    debug_print_migrations(stree);
    debug_print_migmatrix(stree);
  }
  #if 0
  /* TODO: Delete this, just for debugging */
  long i;
  for (i = 3; i < 10; ++i)
  {
    printf("LH for %ld species = %f\n", i, lh_counter(i));
  }
  for (i = 2; i < 10; ++i)
    printf("migmodels for %ld species = %ld\n", i, mig_models_count(i));
  printf("LH for tree: %ld\n", lh_tree_counter(stree));
  fatal("Finishing\n");
  #endif
  if (append)
  {
    #if 0
    printf("[DEBUG %ld] Append -- migration count: %ld\n", debug_counter, opt_migration_count);
    #endif
    rc = stree_migration_append(gtreeptr,
                                gcloneptr,
                                streeptr,
                                scloneptr,
                                locus,
                                si,
                                ti);
  }
  else
  {
    #if 0
    printf("[DEBUG %ld] Remove -- migration count: %ld\n", debug_counter, opt_migration_count);
    #endif
    /* randomly select a migration to delete */
    assert(opt_migration_count);
    long index = (long)(opt_migration_count*legacy_rndu(thread_index_zero));
    spec = opt_mig_specs+index;
    /* TODO: Pass migspec_index instead of si and ti */

    #if (DBG_TF)
    printf("     Preparing to remove M_%s->%s\n", stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label);
    #endif
    rc = stree_migration_remove(gtreeptr,
                                gcloneptr,
                                streeptr,
                                scloneptr,
                                locus,
                                spec->si,
                                spec->ti);
  }

  #if 1
  mig_model_prop_count[dbg_mig_idx][dbg_mig_idx_prop]++;
  if (rc)
    mig_model_prop_acc[dbg_mig_idx][dbg_mig_idx_prop]++;
  #endif
  #if (DBG_TF)
  if (append)
  {
    if (rc)
      printf("APPEND M_%s->%s - ACCEPT %ld\n", stree->nodes[si]->label, stree->nodes[ti]->label, dbg_counter);
    else
      printf("APPEND M_%s->%s - REJECT %ld\n", stree->nodes[si]->label, stree->nodes[ti]->label, dbg_counter);
  }
  else
  {
    if (dbg_counter == 10)
      debug_print_migmatrix(stree);
    if (rc)
      printf("REMOVE M_%s->%s - ACCEPT %ld\n", stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label, dbg_counter);
    else
      printf("REMOVE M_%s->%s - REJECT %ld\n", stree->nodes[spec->si]->label, stree->nodes[spec->ti]->label, dbg_counter);
  }
  #endif
  return rc;
}
