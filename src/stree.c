/*
    Copyright (C) 2016-2019 Tomas Flouri and Ziheng Yang

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

/* species tree spr move related */
#define LINEAGE_A       16
#define LINEAGE_OTHER   32
#define NODE_SQUARE     64
#define NODE_MOVED      256

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

void print_network_table(stree_t * stree)
{
  long i;

  printf("Species tree contains hybridization/introgression events.\n\n");
  printf("Label        Node-index  Child1-index  Child2-index  Parent-index\n");
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

    printf("%-12s", label);
    free(label);
    printf("  %9d", stree->nodes[i]->node_index);

    if (stree->nodes[i]->left)
      printf("  %12d", stree->nodes[i]->left->node_index);
    else
      printf("  %12s", "N/A");

    if (stree->nodes[i]->right)
      printf("  %12d", stree->nodes[i]->right->node_index);
    else
      printf("  %12s", "N/A");


    if (stree->nodes[i]->parent)
      printf("  %12d", stree->nodes[i]->parent->node_index);
    else
      printf("  %12s", "N/A");

    if (i >= stree->tip_count + stree->inner_count)
    {
      printf("   Mirrored hybridization node");
      assert(stree->nodes[i]->hybrid);

      if (stree->nodes[i]->hybrid->label)
        printf(" [Hybrid = %s (%d)]",
               stree->nodes[i]->hybrid->label,
               stree->nodes[i]->hybrid->node_index);
    }

    if (stree->nodes[i]->hybrid)
      printf("   [tau = %ld, phi = %f]", stree->nodes[i]->htau, stree->nodes[i]->hphi);

    printf("\n");
  }
  printf("\n");
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
    printf("  %ld", i);
  printf("\n");

  for (i = 0; i < nodes_count; ++i)
  {
    printf("%*ld %-*s ", index_digits, i, (int)maxlen, stree->nodes[i]->label);

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

   /* seqin (incoming sequences to a population) counts per locus */
   if (!clone->seqin_count)
      clone->seqin_count = (int *)xmalloc(msa_count * sizeof(int));
   memcpy(clone->seqin_count, snode->seqin_count, msa_count * sizeof(int));

   /* per locus number of gene leaves at each clade */
   if (!clone->gene_leaves)
      clone->gene_leaves = (unsigned int *)xmalloc(msa_count *
         sizeof(unsigned int));
   memcpy(clone->gene_leaves, snode->gene_leaves, msa_count * sizeof(unsigned int));

   /* gene tree probability contributions for current population */
   if (!clone->logpr_contrib)
      clone->logpr_contrib = (double *)xmalloc(msa_count * sizeof(double));
   memcpy(clone->logpr_contrib, snode->logpr_contrib, msa_count * sizeof(double));

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
         clone->old_t2h = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
      memcpy(clone->old_t2h, snode->old_t2h, opt_locus_count * sizeof(double));
   }
}

static void gnode_clone(gnode_t * gnode, gnode_t * clone, gtree_t * clone_gtree, stree_t * clone_stree)
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
}

static void stree_clone(stree_t * stree, stree_t * clone)
{
   unsigned int i;
   unsigned nodes_count = stree->tip_count + stree->inner_count;

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

}

stree_t * stree_clone_init(stree_t * stree)
{
   unsigned int i;
   unsigned nodes_count = stree->tip_count + stree->inner_count;
   stree_t * clone;

   clone = (stree_t *)xcalloc(1, sizeof(stree_t));
   memcpy(clone, stree, sizeof(stree_t));

   /* create cloned species tree nodes */
   clone->nodes = (snode_t **)xmalloc(nodes_count * sizeof(snode_t *));
   for (i = 0; i < nodes_count; ++i)
      clone->nodes[i] = (snode_t *)xcalloc(1, sizeof(snode_t));
   for (i = 0; i < nodes_count; ++i)
      snode_clone(stree->nodes[i], clone->nodes[i], clone);

   clone->pptable = (int **)xmalloc(nodes_count * sizeof(int *));
   for (i = 0; i < nodes_count; ++i)
   {
      clone->pptable[i] = (int *)xmalloc(nodes_count * sizeof(int));
      memcpy(clone->pptable[i], stree->pptable[i], nodes_count * sizeof(int));
   }
   clone->root = clone->nodes[stree->root->node_index];

   return clone;
}

static void gtree_clone(gtree_t * gtree, gtree_t * clone_gtree, stree_t * clone_stree)
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
}

gtree_t * gtree_clone_init(gtree_t * gtree, stree_t * clone_stree)
{
   unsigned int i;
   unsigned nodes_count = gtree->tip_count + gtree->inner_count;
   gtree_t * clone;

   clone = (gtree_t *)xcalloc(1, sizeof(gtree_t));
   memcpy(clone, gtree, sizeof(gtree_t));

   /* create cloned gene tree nodes */
   clone->nodes = (gnode_t **)xmalloc(nodes_count * sizeof(gnode_t *));
   for (i = 0; i < nodes_count; ++i)
      clone->nodes[i] = (gnode_t *)xcalloc(1, sizeof(gnode_t));
   for (i = 0; i < nodes_count; ++i)
      gnode_clone(gtree->nodes[i], clone->nodes[i], clone, clone_stree);

   clone->root = clone->nodes[gtree->root->node_index];

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

            dlist_item_t * cloned = dlist_append(clone_stree->nodes[i]->event[j],
               cloned_node);

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

  stree_label_recursive(node->left);
  stree_label_recursive(node->right);

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
            fatal("Internal error in populations_count");

         node = (snode_t *)(pair->data);

         if (opt_debug)
            printf("Matched %s -> %s\n", label, node->label);

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

          if (x->parent->htau && x->parent->tau == 1)
          {
              run  = 1;
              continue;
          }
          if (x->hybrid->parent->htau && x->hybrid->parent->tau == 1)
          {
            run = 1;
            continue;
          }

          if (x->parent->htau == 0)
          {
            assert(x->parent->parent);
            assert(x->parent->parent->tau > 0);
            if (x->parent->parent->tau == 1)
            {
              run = 1;
              continue;
            }
          }

          if (x->hybrid->parent->htau == 0)
          {
            assert(x->hybrid->parent->parent);
            assert(x->hybrid->parent->parent->tau > 0);
            if (x->hybrid->parent->parent->tau == 1)
            {
              run = 1;
              continue;
            }
          }
          double age1 = (x->parent->htau) ?
                          x->parent->tau : x->parent->parent->tau;
          double age2 = (x->hybrid->parent->htau) ?
                          x->hybrid->parent->tau : x->hybrid->parent->parent->tau;

          if (x->tau != 1)
            continue;

          x->tau = MIN(age1,age2) *
                   (prop + (1 - prop - 0.02)*legacy_rndu(thread_index));
          x->hybrid->tau = x->tau;
          if (x->parent->htau == 0)
            x->parent->tau = x->tau;
          if (x->hybrid->parent->htau == 0)
            x->hybrid->parent->tau = x->tau;
        }
        else
        {

          /* bidirectional introgression nodes */

          assert(node_is_bidirection(x));

          assert(x->parent->htau || 
                 (x->parent->hybrid && node_is_bidirection(x->parent)));

          /* TODO: Account for parallel bidirectional introgressions among two
             lineages with no speciations in between. The thing we need to check
             is that parent->htau may be 0 when visiting the bottom bidirection */

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

          double age = x->parent->tau*
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
              if (x->htau)
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

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
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
     stree->root->tau = opt_tau_beta / (opt_tau_alpha - 1) *
                        (0.9 + 0.2*legacy_rndu(thread_index));

   /* recursively set the speciation time for the remaining inner nodes. For
      networks it is not necessary to check if root has both left and right */
   if (opt_msci)
     network_init_tau_iterative(stree, prop, thread_index);
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
  long i;
  double phinew;
  double phiold;
  double new_logpr;
  double old_logpr;
  double lnacceptance;

  assert(!node_is_mirror(snode));

  phiold = snode->hphi;
  phinew = phiold + opt_finetune_phi*legacy_rnd_symmetrical(thread_index);
  phinew = reflect(phinew,0,1,thread_index);

  if (opt_est_theta)
  {
    old_logpr = 0;
    new_logpr = 0;
    for (i = 0; i < stree->locus_count; ++i)
    {
      old_logpr += gtree[i]->logpr;
      new_logpr += gtree[i]->logpr +
                   snode->seqin_count[i]*(log(phinew) - log(phiold)) +
                   snode->hybrid->seqin_count[i]*(log(1-phinew) - log(1-phiold));
    }
  }
  else
  {
    old_logpr = stree->notheta_logpr;
    new_logpr = stree->notheta_logpr;
    for (i = 0; i < stree->locus_count; ++i)
    {
      new_logpr += snode->seqin_count[i]*(log(phinew) - log(phiold)) +
                   snode->hybrid->seqin_count[i]*(log(1-phinew) - log(1-phiold));
    }
  }

  lnacceptance = (opt_phi_alpha-1) * (log(phinew) - log(phiold)) +
                 (opt_phi_beta-1) * (log(1-phinew) - log(1-phiold)) +
                 new_logpr - old_logpr;

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */

    accepted = 1;
    snode->hphi = phinew;
    snode->hybrid->hphi = 1-phinew;

    /* update logpr */
    if (opt_est_theta)
    {
      for (i = 0; i < stree->locus_count; ++i)
      {
        /* subtract from gene tree log-density the old phi contributions */
        gtree[i]->logpr = gtree[i]->logpr -
                          snode->logpr_contrib[i] -
                          snode->hybrid->logpr_contrib[i];

        /* update log-density contributions for the two populations */
        snode->logpr_contrib[i] += snode->seqin_count[i] *
                                   (log(phinew) - log(phiold));
        snode->hybrid->logpr_contrib[i] += snode->hybrid->seqin_count[i] *
                                           (log(1-phinew) - log(1-phiold));

        /* add to the gene tree log-density the new phi contributions */
        gtree[i]->logpr = gtree[i]->logpr +
                          snode->logpr_contrib[i] +
                          snode->hybrid->logpr_contrib[i];
      }
    }
    else
    {
      /* subtract from total density (sum of densities of all gene trees) the
         old phi contributions from the two populations */
      stree->notheta_logpr -= snode->notheta_logpr_contrib;
      stree->notheta_logpr -= snode->hybrid->notheta_logpr_contrib;

      /* compute new contributions */
      for (i = 0; i < stree->locus_count; ++i)
      {
        snode->notheta_logpr_contrib = snode->notheta_logpr_contrib +
                                       snode->seqin_count[i] * 
                                       (log(phinew) - log(phiold));
        snode->hybrid->notheta_logpr_contrib = snode->hybrid->notheta_logpr_contrib +
                                               snode->hybrid->seqin_count[i] * 
                                               (log(1-phinew)-log(1-phiold));
        
      }

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

  if (opt_phi_alpha <= 0)
    fatal("Alpha value for 'phiprior' must be larger than 0");
  if (opt_phi_beta <= 0)
    fatal("Beta value for 'phiprior' must be larger than 0");

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * mnode = stree->nodes[offset+i];

    assert(node_is_bidirection(mnode) || node_is_mirror(mnode));

    /* set phi parameter to the mean */
    mnode->hybrid->hphi = (opt_phi_alpha / (opt_phi_alpha + opt_phi_beta));
    mnode->hphi = 1 - mnode->hybrid->hphi;
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
   for (i = 0; i < stree->tip_count; ++i)
   {
     snode_t * node = stree->nodes[i];

     for (j = 0; j < (unsigned int)msa_count; ++j)
       if (seqcount[i][j] >= 2)
         break;

     /* if no loci exists with two or more sequences of such species then move
        to the next tip node */
     if (opt_sp_seqcount[i] < 2)
     {
       node->theta = -1;
       node->has_theta = 0;
       continue;
     }

     /* otherwise set theta around the mean of the inverse gamma prior */
     node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                   (0.9 + 0.2 * legacy_rndu(thread_index));
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

        if (node->parent->htau)
        {
          node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                        (0.9 + 0.2 * legacy_rndu(thread_index));
          node->has_theta = 1;
        }
        else
        {
          node->theta = -1;
          node->has_theta = 0;
        }
        
        if (node->hybrid->parent->htau)
        {
          node->hybrid->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                                (0.9 + 0.2 * legacy_rndu(thread_index));
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

        node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                      (0.9 + 0.2 * legacy_rndu(thread_index));
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
         node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                       (0.9 + 0.2 * legacy_rndu(thread_index));
    }
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

int node_is_bidirection(snode_t * node)
{
  assert(node && node->hybrid);

  snode_t * hnode;    /* hybrid node */
  snode_t * mnode;    /* mirror node */

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
int node_is_mirror(snode_t * node)
{
  assert(node);
  assert(node->hybrid);
  assert(node && node->hybrid);

  if (node->left || node->right)
    return 0;

  return 1;
}

int node_is_hybridization(snode_t * node)
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

void stree_init(stree_t * stree,
                msa_t ** msa,
                list_t * maplist,
                int msa_count,
                FILE * fp_out)
{
  unsigned int i, j;

  long thread_index = 0;

  /* safety check */
  assert(msa_count > 0);
  assert(opt_msci == !!stree->hybrid_count);

  stree_init_pptable(stree);

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

  if (opt_msci)
    stree_init_phi(stree);

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
    /* TODO: Implement branch rates for MSci */
    assert(!opt_msci);
    for (i=0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
      stree->nodes[i]->brate = (double *)xmalloc((size_t)msa_count*sizeof(double));
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
}

static int propose_theta(gtree_t ** gtree,
                         locus_t ** locus,
                         snode_t * snode,
                         long thread_index)
{
   long i;
   double thetaold;
   double thetanew;
   double lnacceptance;

   thetaold = snode->theta;

   thetanew = thetaold + opt_finetune_theta*legacy_rnd_symmetrical(thread_index);

   if (thetanew < 0)
      thetanew = -thetanew;

   snode->theta = thetanew;

   lnacceptance = (-opt_theta_alpha - 1) * log(thetanew / thetaold) -
      opt_theta_beta*(1 / thetanew - 1 / thetaold);

   for (i = 0; i < opt_locus_count; ++i)
   {
      /* save a copy of old logpr */
      gtree[i]->old_logpr = gtree[i]->logpr;

      gtree[i]->logpr -= snode->logpr_contrib[i];
      gtree_update_logprob_contrib(snode, locus[i]->heredity[0], i, thread_index);
      gtree[i]->logpr += snode->logpr_contrib[i];

      lnacceptance += (gtree[i]->logpr - gtree[i]->old_logpr);

   }

   if (opt_debug)
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

   snode->theta = thetaold;
   for (i = 0; i < opt_locus_count; ++i)
      snode->logpr_contrib[i] = snode->old_logpr_contrib[i];

   return 0;
}

double stree_propose_theta(gtree_t ** gtree, locus_t ** locus, stree_t * stree)
{
   unsigned int i;
   int theta_count = 0;
   long accepted = 0;
   snode_t * snode;

   long thread_index = 0;

   for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
   {
      snode = stree->nodes[i];
      if (snode->theta >= 0 && snode->has_theta)
      {
         accepted += propose_theta(gtree, locus, stree->nodes[i], thread_index);
         theta_count++;
      }
   }

   return ((double)accepted / theta_count);
}

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
  long i;
  unsigned int j,k;
  unsigned int locus_count_above;
  unsigned int locus_count_below;
  unsigned int count_above = 0;
  unsigned int count_below = 0;
  double logl_diff = 0;
  double logpr_diff = 0;
  double logpr = 0;

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
     double * oldbranches = __aux + __gt_nodes_index[i];
     #endif

     /* traverse the gene tree nodes of the three populations, find the ones
        whose ages fall within the new age interval, update their age and mark
        them. Also, count how many of them are above the old species node age,
        and how many are below. Finally, update the gene tree probabilities */
     for (j = 0; j < paffected_count; ++j)
     {
        /* process events for current population */
        if (affected[j]->event_count)
        {
           dlist_item_t * event;
           for (event = affected[j]->event[i]->head; event; event = event->next)
           {
              gnode_t * node = (gnode_t *)(event->data);
              //if (node->time < minage) continue;
              if ((node->time < minage) || (node->time > maxage)) continue;

              gt_nodesptr[k] = node;
              node->mark = FLAG_PARTIAL_UPDATE;
              oldbranches[k++] = node->time;

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

           double xtmp = gtree_update_logprob_contrib(affected[j],
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
              SWAP(oldbranches[0], oldbranches[j]);
           }
           branchptr = &(gt_nodesptr[1]);
           --branch_count;
        }

        assert(node->left);
        assert(node->right);

        if (!node->left->mark)
        {
           branchptr[branch_count++] = node->left;
           extra++;
        }
        if (!node->right->mark)
        {
           branchptr[branch_count++] = node->right;
           extra++;
        }
     }

     __extra_count[i] = extra;
     #ifdef OLD_CODE
     offset += extra;
     #endif

     /* if at least one gene tree node age was changed, we need to recompute the
        log-likelihood */
     gtree[i]->old_logl = gtree[i]->logl;
     if (k)
     {
        locus_update_matrices(loci[i], branchptr, branch_count);

        /* get list of nodes for which partials must be recomputed */
        unsigned int partials_count;
        gnode_t ** partials = gtree_return_partials(gtree[i]->root,
                                                    i,
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
        unsigned int param_indices[1] = { 0 };
        double logl = locus_root_loglikelihood(loci[i],
           gtree[i]->root,
           param_indices,
           NULL);

        logl_diff += logl - gtree[i]->logl;
        gtree[i]->logl = logl;
     }

     /* Test for checking whether all gene tree nodes can be marked. It seems to
        hold */
        //if (__mark_count[i] + __extra_count[i] >= 2*loci[i]->tips-1)
        //  assert(0);
  }
  *ret_logpr_diff = logpr_diff;
  *ret_logl_diff = logl_diff;
  *ret_count_above = count_above;
  *ret_count_below = count_below;
}
#endif

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

   /* compute minage and maxage bounds */
   if (opt_msci && snode->hybrid)
   {
     if (node_is_hybridization(snode))
     {
       /* hybridization */
       assert(!node_is_mirror(snode));
       assert(snode->left && !snode->right);

       minage = snode->left->tau;
       if (!snode->parent->htau)
       {
         snode_t * sibling = (snode->parent->left == snode) ?
                               snode->parent->right : snode->parent->left;

         minage = MAX(minage,sibling->tau);
       }
       if (!snode->hybrid->parent->htau)
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
              !snode->hybrid->parent->htau && !snode->right->htau);

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
       assert(snode->parent->htau || snode->parent->parent);
       assert(snode->hybrid->parent->htau || snode->hybrid->parent->parent);

       double p1age = snode->parent->htau ?
                        snode->parent->tau : snode->parent->parent->tau;
       double p2age = snode->hybrid->parent->htau ?
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

       if (!snode->parent->htau)
         snode->parent->tau = snode->tau;
       if (!snode->hybrid->parent->htau)
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
      lnacceptance = (-opt_tau_alpha - 1 - candidate_count + 1)*log(newage / oldage) -
      opt_tau_beta*(1 / newage - 1 / oldage);

   /* change theta as well */
   if (opt_est_theta)
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

        lnacceptance += -log(thetafactor) + (-opt_theta_alpha - 1) *
                        log(snode->theta / oldtheta) -
                        opt_theta_beta*(1 / snode->theta - 1 / oldtheta);
      }
   }

   snode_t * affected[7];
   double old_logpr_contrib[7] = {0,0,0,0,0,0,0};

   if (opt_msci && snode->hybrid)
   {
     if (node_is_hybridization(snode))
     {
       if (snode->parent->htau && snode->hybrid->parent->htau)
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
       else if (!snode->parent->htau && !snode->hybrid->parent->htau)
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

         assert((!snode->parent->htau && snode->hybrid->parent->htau) ||
                (snode->parent->htau && !snode->hybrid->parent->htau));

         paffected_count = 0;

         /* assertions */
         if (!snode->parent->htau)
           for (i = 0; i < stree->locus_count; ++i)
             assert(snode->event_count[i] == 0);
         if (!snode->hybrid->parent->htau)
           for (i = 0; i < stree->locus_count; ++i)
             assert(snode->hybrid->event_count[i] == 0);

         if (!snode->parent->htau)
         {
           affected[paffected_count++] = snode->parent;
           affected[paffected_count++] = snode->parent->left;
           affected[paffected_count++] = snode->parent->right;
           affected[paffected_count++] = snode->hybrid;
         }
         else
         {
           assert(!snode->hybrid->parent->htau);

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
   if (opt_threads > 1)
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

      logpr_diff = logpr - stree->notheta_logpr;
      stree->notheta_old_logpr = stree->notheta_logpr;
      stree->notheta_logpr = logpr;
   }

   lnacceptance += logpr_diff + logl_diff + count_below*log(minfactor) +
                   count_above*log(maxfactor);

   if (opt_debug)
      printf("[Debug] (tau) lnacceptance = %f\n", lnacceptance);

   if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
   {
      /* accepted */
      accepted++;

      for (i = 0; i < stree->locus_count; ++i)
      {
         k = __mark_count[i];
         gnode_t ** gt_nodesptr = __gt_nodes + __gt_nodes_index[i];

         for (j = 0; j < k; ++j)
            gt_nodesptr[j]->mark = 0;
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

          if (!snode->parent->htau)
            snode->parent->tau = snode->tau;
          if (!snode->hybrid->parent->htau)
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

      if (opt_est_theta)
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
         double * old_ageptr = __aux + __gt_nodes_index[i];

         /* restore gene tree node ages */
         for (j = 0; j < k; ++j)
            gt_nodesptr[j]->time = old_ageptr[j];

         /* restore logpr contributions */
         if (opt_est_theta)
           for (j = 0; j < paffected_count; ++j)
             gtree_update_logprob_contrib(affected[j],
                                          loci[i]->heredity[0],
                                          i,
                                          thread_index);
         else
         {
           for (j = 0; j < paffected_count; ++j)
             logprob_revert_notheta(affected[j], i);
         }

         /* get the list of nodes for which CLVs must be reverted, i.e. all marked
            nodes and all nodes whose left or right subtree has at least one marked
            node */
         unsigned int partials_count;
         gnode_t ** partials = gtree_return_partials(gtree[i]->root,
                                                     i,
                                                     &partials_count);

         /* revert CLV indices */
         for (j = 0; j < partials_count; ++j)
         {
            partials[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                    partials[j]->clv_index);
            if (opt_scaling)
              partials[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                    partials[j]->scaler_index);
         }

         /* un-mark nodes */
         for (j = 0; j < k; ++j)
            gt_nodesptr[j]->mark = 0;

         /* restore branch lengths and pmatrices */
         int matrix_updates = __mark_count[i] + __extra_count[i];
         if (matrix_updates)
         {
            if (!gt_nodesptr[0]->parent)
            {
               --matrix_updates;
               gt_nodesptr++;
            }
            if (matrix_updates)
               locus_update_matrices(loci[i], gt_nodesptr, matrix_updates);
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

double stree_propose_tau(gtree_t ** gtree, stree_t * stree, locus_t ** loci)
{
   unsigned int i;
   unsigned int candidate_count = 0;
   long accepted = 0;

   long thread_index = 0;

   /* compute number of nodes with tau > 0 */
   for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      if (stree->nodes[i]->tau > 0 && (!opt_msci || stree->nodes[i]->htau))
         candidate_count++;


   for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
   {
      if (stree->nodes[i]->tau > 0 && (!opt_msci || stree->nodes[i]->htau))
         accepted += propose_tau(loci,
                                 stree->nodes[i],
                                 gtree,
                                 stree,
                                 candidate_count,
                                 thread_index);
   }

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

   stree->root_age = opt_tau_beta / (opt_tau_alpha - 1) * 4;

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

static void init_weights(stree_t * stree)
{
   unsigned int i;
   double sum = 0;

   for (i = stree->tip_count; i < stree->inner_count + stree->tip_count; ++i)
   {
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

   /* calculate the weight of each branch as the reciprocal of the square root
    * of its length */
   init_weights(stree);

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

      /* A C candidate node must fulfill the following three properties:
         i) it is not a descendant of y,
         ii) is younger than y,
         iii) its parent is older than y */
      if (stree->pptable[i][y->node_index] ||
         c_cand->tau >= y->tau ||
         c_cand->parent->tau <= y->tau) continue;

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

            if (opt_revolutionary_spr_debug) {
               printf("ntarget = %2ld weight = %9.6f  nsource = %2ld weight = %9.6f\n", target_count, twgt, source_count, swgt);
            }
            
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

   init_weights(stree);

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
      if (moved_count[i])
      {
         /* update branch lengths and transition probability matrices */
         for (j = 0; j < (unsigned int)__mark_count[i]; ++j)
         {
            bl_list[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree_list[i]->edge_count,
               bl_list[j]->pmatrix_index);
         }

         locus_update_matrices(loci[i], bl_list, __mark_count[i]);

         /* retrieve all nodes whose partials must be updates */
         unsigned int partials_count;
         gnode_t ** partials = gtree_return_partials(gtree_list[i]->root,
            i,
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
         unsigned int param_indices[1] = { 0 };
         gtree_list[i]->logl = locus_root_loglikelihood(loci[i],
            gtree_list[i]->root,
            param_indices,
            NULL);
      }

      if (opt_est_theta)
         gtree_list[i]->old_logpr = gtree_list[i]->logpr;
#if 0
      /* This recomputes the gene tree probabilities from scratch. It can be used to verify that
         the code below, which only computes the gene tree probability for the changed components,
         is correct. */

      double logpr = gtree_logprob(stree, loci[i]->heredity[0], i);
#else

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

         double xtmp = gtree_update_logprob_contrib(snode_contrib[j], loci[i]->heredity[0], i, thread_index);

         if (opt_est_theta)
            gtree_list[i]->logpr += snode_contrib[j]->logpr_contrib[i];
         else
            logpr_notheta += xtmp;
      }

      /* reset markings on affected populations */
      for (j = 0; j < snode_contrib_count[i]; ++j)
         snode_contrib[j]->mark[thread_index] = 0;
#endif


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

   if (opt_debug)
      printf("[Debug] (SSPR) lnacceptance = %f\n", lnacceptance);

   /* in case of acceptance, cloned trees are substituted with the original ones,
      and species tree nodes are re-labeled, but all this is done in method_01.c
   */
   //return (lnacceptance >= 0 || legacy_rndu() < exp(lnacceptance));
   return (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance));
}
