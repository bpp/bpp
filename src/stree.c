/*
    Copyright (C) 2016-2017 Tomas Flouri and Ziheng Yang

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

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "bpp.h"

#define PROP_THRESHOLD 10

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))
#define SWAP_PMAT_INDEX(e,i) (((e)+(i))%((e)<<1))


/* species tree spr move related */
#define LINEAGE_A       16
#define LINEAGE_OTHER   32
#define NODE_SQUARE     64
#define NODE_MOVED      256

static gnode_t ** __gt_nodes = NULL;
static double * __aux = NULL;
static int * __mark_count = NULL;
static int * __extra_count = NULL;

static double * target_weight = NULL;
static snode_t ** target = NULL;

/* allocated and used only for species tree inference */
static gnode_t ** moved_space;
static unsigned int * moved_count;
static gnode_t ** gtarget_temp_space;
static gnode_t ** gtarget_space;
static unsigned int * target_count;

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
  return x ? (int)floor(log10(abs(x)))+1 : 1;
}

void stree_show_pptable(stree_t * stree)
{
  long i,j;
  long nodes_count = stree->tip_count + stree->inner_count;
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
  printf("  %ld",i);
  printf("\n");

  for (i = 0; i < nodes_count; ++i)
  {
    printf("%*ld %-*s ", index_digits, i, (int)maxlen, stree->nodes[i]->label);

    for (j = 0; j < nodes_count; ++j)
      printf("  %*d", longint_len(j), stree->pptable[i][j]);

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

  clone->length      = snode->length;
  clone->theta       = snode->theta;
  clone->tau         = snode->tau;
  clone->old_tau     = snode->old_tau;
  clone->old_theta   = snode->old_theta;
  clone->leaves      = snode->leaves;
  clone->mark        = snode->mark;
  clone->support     = snode->support;
  clone->weight      = snode->weight;
  clone->node_index  = snode->node_index;

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
    clone->gene_leaves = (unsigned int *)xmalloc(msa_count *
                                                 sizeof(unsigned int));
  memcpy(clone->gene_leaves,
         snode->gene_leaves,
         msa_count*sizeof(unsigned int)); 

  /* data  - unused */
  clone->data = NULL;

  /* event doubly-linked lists */
  if (!clone->event)
  {
    clone->event = (dlist_t **)xcalloc(msa_count,sizeof(dlist_t *));
    for (i = 0; i < msa_count; ++i)
      clone->event[i] = dlist_create();
  }
  else
  {
    for (i = 0; i < msa_count; ++i)
      dlist_clear(clone->event[i],NULL);
  }
    
  /* event counts per locus */
  if (!clone->event_count)
    clone->event_count = (int *)xmalloc(msa_count*sizeof(int));
  memcpy(clone->event_count,snode->event_count,msa_count*sizeof(int));

  /* seqin (incoming sequences to a population) counts per locus */
  if (!clone->seqin_count)
    clone->seqin_count = (int *)xmalloc(msa_count*sizeof(int));
  memcpy(clone->seqin_count,snode->seqin_count,msa_count*sizeof(int));

  /* per locus number of gene leaves at each clade */
  if (!clone->gene_leaves)
    clone->gene_leaves = (unsigned int *)xmalloc(msa_count *
                                                 sizeof(unsigned int));
  memcpy(clone->gene_leaves,snode->gene_leaves,msa_count*sizeof(unsigned int));

  /* gene tree probability contributions for current population */
  if (!clone->logpr_contrib)
    clone->logpr_contrib = (double *)xmalloc(msa_count*sizeof(double));
  memcpy(clone->logpr_contrib,snode->logpr_contrib,msa_count*sizeof(double));

  /* old gene tree probability contributions for current population */
  if (!clone->old_logpr_contrib)
    clone->old_logpr_contrib = (double *)xmalloc(msa_count*sizeof(double));
  memcpy(clone->old_logpr_contrib,
         snode->old_logpr_contrib,
         msa_count*sizeof(double));
}

static void gnode_clone(gnode_t * gnode, gnode_t * clone, gtree_t * clone_gtree, stree_t * clone_stree)
{
  if (clone->label)
    free(clone->label);

  clone->length = gnode->length;
  clone->time   = gnode->time;
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
    memcpy(clone->pptable[i], stree->pptable[i], nodes_count*sizeof(int));

  clone->root = clone->nodes[stree->root->node_index];
}

stree_t * stree_clone_init(stree_t * stree)
{
  unsigned int i;
  unsigned nodes_count = stree->tip_count + stree->inner_count;
  stree_t * clone;

  clone = (stree_t *)xcalloc(1, sizeof(stree_t));
  memcpy(clone, stree, sizeof(stree_t));

  /* create cloned species tree nodes */
  clone->nodes = (snode_t **)xmalloc(nodes_count*sizeof(snode_t *));
  for (i = 0; i < nodes_count; ++i)
    clone->nodes[i] = (snode_t *)xcalloc(1,sizeof(snode_t));
  for (i = 0; i < nodes_count; ++i)
    snode_clone(stree->nodes[i], clone->nodes[i], clone);

  clone->pptable = (int **)xmalloc(nodes_count*sizeof(int *));
  for (i = 0; i < nodes_count; ++i)
  {
    clone->pptable[i] = (int *)xmalloc(nodes_count*sizeof(int));
    memcpy(clone->pptable[i], stree->pptable[i], nodes_count*sizeof(int));
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
  clone->nodes = (gnode_t **)xmalloc(nodes_count*sizeof(gnode_t *));
  for (i = 0; i < nodes_count; ++i)
    clone->nodes[i] = (gnode_t *)xcalloc(1,sizeof(gnode_t));
  for (i = 0; i < nodes_count; ++i)
    gnode_clone(gtree->nodes[i], clone->nodes[i], clone, clone_stree);

  clone->root = clone->nodes[gtree->root->node_index];

  return clone;
}

static void events_clone(stree_t * stree,
                         stree_t * clone_stree,
                         gtree_t ** clone_gtree_list)
{
  unsigned int i,j;
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
  strcat(node->label,node->left->label);
  strcat(node->label,node->right->label);
}

void stree_label(stree_t * stree)
{
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
  unsigned int i,j,k;
  snode_t * node;

  assert(msa_count >= 0);

  /* create one hash table of species and one for sequence->species mappings */
  hashtable_t * sht = species_hash(stree);
  hashtable_t * mht = maplist_hash(maplist,sht);

  /* create perloci sequence counters for tip and inner nodes */
  int ** seqcount = (int **)xmalloc(stree->tip_count * sizeof(int *));
  for (i = 0; i < stree->tip_count; ++i)
    seqcount[i] = (int *)xcalloc(msa_count,sizeof(int));

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

  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);

  return seqcount;
}

static void stree_init_tau_recursive(snode_t * node,
                                     double prop)
{
  /* end recursion if node is a tip */
  if (!node->left) return;

  /* get species record associate with species tree node */
  double tau_parent = node->parent->tau;

  if (node->tau)
    node->tau = tau_parent * (prop + (1 - prop - 0.02)*legacy_rndu());
  else
    node->left->theta = node->right->theta = -1;

  stree_init_tau_recursive(node->left,prop);
  stree_init_tau_recursive(node->right,prop);
}

static void stree_init_tau(stree_t * stree)
{
  unsigned int i;

  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    stree->nodes[i]->tau = 1;
  
  if (opt_delimit && !opt_stree)    /* method A10 */
  {
    double r = legacy_rndu();
    int index = (int)(r * delimitation_getparam_count());
    delimitation_set(stree,index);

    char * s = delimitation_getparam_string();
    printf("Starting delimitation: %s\n", s);
  }
  /* Initialize speciation times for each extinct species */ 

  double prop = (stree->root->leaves > PROP_THRESHOLD) ? 0.9 : 0.5;

  /* set the speciation time for root */
  if (stree->root->tau)
    stree->root->tau = opt_tau_beta/(opt_tau_alpha-1) * (0.9 + 0.2*legacy_rndu());

  /* recursively set the speciation time for the remaining inner nodes */
  stree_init_tau_recursive(stree->root->left,prop);
  stree_init_tau_recursive(stree->root->right,prop);
}

static void stree_init_theta(stree_t * stree,
                             msa_t ** msalist,
                             list_t * maplist,
                             int msa_count)
{
  unsigned int i,j,k = 0;

  /* initialize population sizes for extinct populations and populations
     with more than one lineage at some locus */

  /* get an array of per-locus number of sequences for each species (tip in
     species tree) */
  int ** seqcount = populations_seqcount(stree,msalist,maplist,msa_count);

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
      continue;
    }

    /* otherwise set theta around the mean of the inverse gamma prior */
    node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                        (0.9 + 0.2 * legacy_rndu());
    ++k;
  }

  /* go through inner nodes and setup thetas */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * node = stree->nodes[i];

    node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                  (0.9 + 0.2 * legacy_rndu());
    ++k;
  }

  /* deallocate seqcount */
  for (i = 0; i < stree->tip_count; ++i)
    free(seqcount[i]);
  free(seqcount);
}

void stree_reset_pptable(stree_t * stree)
{
  unsigned int i;
  snode_t * ancnode;
  snode_t * curnode;

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    memset(stree->pptable[i],
           0,
           (stree->tip_count+stree->inner_count)*sizeof(int));

  for (i = 0; i < stree->tip_count; ++i)
  {
    for (curnode = stree->nodes[i]; curnode; curnode = curnode->parent)
      for (ancnode = curnode; ancnode; ancnode = ancnode->parent)
        stree->pptable[curnode->node_index][ancnode->node_index] = 1;
  }
}
void stree_alloc_internals(stree_t * stree, unsigned int gtree_inner_sum, long msa_count)
{
  /* allocate traversal buffer to be the size of all nodes for all loci */
//  unsigned int sum_count = 0;
  unsigned int sum_nodes = 2*gtree_inner_sum + msa_count;
//  for (i = 0; i < (unsigned int)msa_count; ++i)
//    sum_count += (unsigned int)(msa[i]->count);

//  sum_count = 2*sum_count - msa_count;
//  __gt_nodes = (gnode_t **)xmalloc(sum_count * sizeof(gnode_t *));
//  __aux = (double *)xmalloc((sum_count - msa_count)*sizeof(double *));
  __gt_nodes = (gnode_t **)xmalloc(sum_nodes * sizeof(gnode_t *));
  __aux = (double *)xmalloc((sum_nodes - msa_count)*sizeof(double *));

  /* The following two arrays are used purely for the tau proposal.
     Entry i of marked_count indicates how many nodes from locus i are marked.
     Similarly, extra_count 
     how many extra nodes where added in __gt_nodes whose branch lengths (and
     therefore) p-matrices need updating because their parent node's age was
     changed */
  __mark_count  = (int *)xmalloc(msa_count*sizeof(int));
  __extra_count = (int *)xmalloc(msa_count*sizeof(int));

  /* species tree inference */
  if (opt_stree)
  {
    target_weight = (double *)xmalloc((stree->tip_count + stree->inner_count) *
                                    sizeof(double));
    target = (snode_t **)xmalloc((stree->tip_count + stree->inner_count) *
                                 sizeof(snode_t *));

    /* TODO: memory is allocated for all loci to aid parallelization */
    moved_count = (unsigned int *)xcalloc(msa_count,sizeof(unsigned int));
    target_count = (unsigned int *)xcalloc(msa_count,sizeof(unsigned int));
//    unsigned int sum_inner = 0;
//    for (i = 0; i < (unsigned int)msa_count; ++i)
//      sum_inner += msa[i]->count-1;
//    moved_space = (gnode_t **)xmalloc(sum_inner * sizeof(gnode_t *));
//    gtarget_space = (gnode_t **)xmalloc(sum_inner * sizeof(gnode_t *));
    moved_space = (gnode_t **)xmalloc(gtree_inner_sum * sizeof(gnode_t *));
    gtarget_space = (gnode_t **)xmalloc(gtree_inner_sum * sizeof(gnode_t *));

//    unsigned int sum_nodes = 2*sum_inner + msa_count;
    gtarget_temp_space = (gnode_t **)xmalloc(sum_nodes * sizeof(gnode_t *));

    unsigned int stree_nodes = stree->inner_count + stree->tip_count;
    snode_contrib_space = (snode_t **)xmalloc((size_t)(msa_count*stree_nodes) *
                                              sizeof(snode_t *));
    snode_contrib_count = (unsigned int *)xmalloc((size_t)msa_count *
                                                  sizeof(unsigned int));
  }
}

void stree_init(stree_t * stree, msa_t ** msa, list_t * maplist, int msa_count)
{
  unsigned int i,j;
  #if 0
  snode_t * curnode;
  snode_t * ancnode;
  #endif

  assert(msa_count > 0);

  /* label each inner node of the species tree with the concatenated labels of
     its two children */
  stree_label(stree);

  /* Initialize population sizes */
  stree_init_theta(stree, msa, maplist, msa_count);

  /* Initialize speciation times and create extinct species groups */
  stree_init_tau(stree);

  /* allocate space for keeping track of coalescent events at each species tree
     node for each locus */
  stree->locus_count = (unsigned int)msa_count;
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    stree->nodes[i]->event = (dlist_t **)xcalloc(msa_count,sizeof(dlist_t *));
    stree->nodes[i]->event_count = (int *)xcalloc(msa_count,sizeof(int));
    stree->nodes[i]->seqin_count = (int *)xcalloc(msa_count,sizeof(int));
    stree->nodes[i]->gene_leaves = (unsigned int *)xcalloc(msa_count,
                                                           sizeof(unsigned int));
    stree->nodes[i]->logpr_contrib = (double *)xcalloc(msa_count,sizeof(double));
    stree->nodes[i]->old_logpr_contrib = (double *)xcalloc(msa_count,
                                                           sizeof(double));

    for (j = 0; j < stree->locus_count; ++j)
      stree->nodes[i]->event[j] = dlist_create();
  }

  /* initialize pptable, where pptable[i][j] indicates whether population j
     (that is, node with index j) is ancestral to population i */
  stree->pptable = (int **)xcalloc((stree->tip_count + stree->inner_count),
                                   sizeof(int *));
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    stree->pptable[i] = (int *)xcalloc((stree->tip_count + stree->inner_count),
                                       sizeof(int));

  
  stree_reset_pptable(stree);
  #if 0
  for (i = 0; i < stree->tip_count; ++i)
  {
    for (curnode = stree->nodes[i]; curnode; curnode = curnode->parent)
      for (ancnode = curnode; ancnode; ancnode = ancnode->parent)
        stree->pptable[curnode->node_index][ancnode->node_index] = 1;
  }
  #endif

  unsigned int sum_inner = 0;
  for (i = 0; i < (unsigned int)msa_count; ++i)
    sum_inner += msa[i]->count-1;
  stree_alloc_internals(stree,sum_inner,msa_count);
}

void stree_fini()
{
  free(__gt_nodes);
  free(__aux);
  free(__mark_count);
  free(__extra_count);

  if (opt_stree)
  {
    free(target_weight);
    free(target);
    free(moved_count);
    free(target_count);
    free(moved_space);
    free(gtarget_temp_space);
    free(gtarget_space);
    free(snode_contrib_space);
    free(snode_contrib_count);
  }
}

static int propose_theta(gtree_t ** gtree, locus_t ** locus, snode_t * snode)
{
  long i;
  double thetaold;
  double thetanew;
  double acceptance;

  thetaold = snode->theta;
  
  thetanew = thetaold + opt_finetune_theta * legacy_rnd_symmetrical();
 
  if (thetanew < 0)
    thetanew = -thetanew;

  snode->theta = thetanew;

  acceptance = (-opt_theta_alpha-1) * log(thetanew/thetaold) -
               opt_theta_beta*(1/thetanew - 1/thetaold);

  for (i = 0; i < opt_locus_count; ++i)
  {
    /* save a copy of old logpr */
    gtree[i]->old_logpr = gtree[i]->logpr;

    gtree[i]->logpr -= snode->logpr_contrib[i];
    gtree_update_logprob_contrib(snode,locus[i]->heredity[0],i);
    gtree[i]->logpr += snode->logpr_contrib[i];
    
    acceptance += (gtree[i]->logpr - gtree[i]->old_logpr);

  }

  if (opt_debug)
    printf("[Debug] (theta) lnacceptance = %f\n", acceptance);

  if (acceptance >= 0 || legacy_rndu() < exp(acceptance))
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

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    if (snode->theta >= 0)
    {
      accepted += propose_theta(gtree, locus, stree->nodes[i]);
      theta_count++;
    }
  }

  return ((double)accepted/theta_count);
}

static long propose_tau(locus_t ** loci,
                        snode_t * snode,
                        gtree_t ** gtree,
                        stree_t * stree,
                        unsigned int candidate_count)
{
  unsigned int i,j,k;
  unsigned int offset = 0;
  int changetheta = 1;
  int theta_method = 2;
  long accepted = 0;
  double oldage, newage;
  double minage = 0, maxage = 999;
  double acceptance = 0;
  double minfactor,maxfactor,thetafactor;
  double oldtheta;
  double logpr = 0;
  double logpr_diff = 0;
  double logl_diff = 0;

  unsigned int count_above = 0;
  unsigned int count_below = 0;
  unsigned int locus_count_above;
  unsigned int locus_count_below;

  oldage = snode->tau;

  /* compute minage and maxage bounds */
  if (snode->left)
    minage = MAX(snode->left->tau,snode->right->tau);

  if (snode->parent)
    maxage = snode->parent->tau;

  /* propose new tau */
  newage = oldage + opt_finetune_tau * legacy_rnd_symmetrical();
  newage = reflect(newage,minage,maxage);
  snode->tau = newage;

  /* compute factors for multiplying associated gene tree nodes ages */
  minfactor = (newage - minage) / (oldage - minage);
  maxfactor = (newage - maxage) / (oldage - maxage);

  /* if we are dealing with the root population, add the following factor to
     the acceptance ratio */
  if (snode == stree->root)
    acceptance = (-opt_tau_alpha-1 - candidate_count+1)*log(newage/oldage) -
                 opt_tau_beta*(1/newage - 1/oldage);

  /* change theta as well */
  if (changetheta)
  {
    oldtheta = snode->theta;
    if (theta_method == 1)
      thetafactor = newage / oldage;
    else if (theta_method == 2)
      thetafactor = (newage - minage) / (oldage - minage);
    else
      assert(0);

    snode->theta = oldtheta / thetafactor;

    acceptance += -log(thetafactor) + (-opt_theta_alpha-1) * 
                  log(snode->theta/oldtheta) -
                  opt_theta_beta*(1/snode->theta - 1/oldtheta);
  }

  snode_t * affected[3] = {snode, snode->left, snode->right};
  for (i = 0; i < stree->locus_count; ++i)
  {
    k = 0;
    locus_count_above = locus_count_below = 0;

    logpr = gtree[i]->logpr;

    gnode_t ** gt_nodesptr = __gt_nodes + offset;
    double * oldbranches = __aux + offset;

    /* traverse the gene tree nodes of the three populations, find the ones
       whose ages fall within the new age interval, update their age and mark
       them. Also, count how many of them are above the old species node age,
       and how many are below. Finally, update the gene tree probabilities */
    for (j = 0; j < 3; ++j)
    {
      /* process events for current population */
      if (affected[j]->event_count)
      {
        dlist_item_t * event;
        for (event = affected[j]->event[i]->head; event; event = event->next)
        {
          gnode_t * node = (gnode_t *)(event->data);
          if (node->time < minage) continue;

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
        logpr -= affected[j]->logpr_contrib[i];
        logpr += gtree_update_logprob_contrib(affected[j],loci[i]->heredity[0],i);
      }
    }

    /* entry i of __mark_count holds the number of marked nodes for locus i */
    __mark_count[i] = k;
    offset += k;

    logpr_diff += logpr - gtree[i]->logpr;

    gtree[i]->old_logpr = gtree[i]->logpr;
    gtree[i]->logpr = logpr;

    count_above += locus_count_above;
    count_below += locus_count_below;

    /* TODO: Do the integrated-out theta option */

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
          SWAP(gt_nodesptr[0],gt_nodesptr[j]);
          SWAP(oldbranches[0],oldbranches[j]);
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
    offset += extra;

    /* if at least one gene tree node age was changed, we need to recompute the
       log-likelihood */
    gtree[i]->old_logl = gtree[i]->logl;
    if (k)
    {
      locus_update_matrices_jc69(loci[i],branchptr,branch_count);

      /* get list of nodes for which partials must be recomputed */
      unsigned int partials_count;
      gnode_t ** partials = gtree_return_partials(gtree[i]->root,
                                                  i,
                                                  &partials_count);

      for (j = 0; j < partials_count; ++j)
        partials[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                partials[j]->clv_index);

      /* update partials */
      locus_update_partials(loci[i],partials,partials_count);

      /* evaluate log-likelihood */
      unsigned int param_indices[1] = {0};
      double logl = locus_root_loglikelihood(loci[i],
                                             gtree[i]->root,
                                             param_indices,
                                             NULL);

      logl_diff += logl - gtree[i]->logl;
      gtree[i]->logl = logl;
    }
  }

  acceptance += logpr_diff + logl_diff + count_below*log(minfactor) +
                count_above*log(maxfactor);

  if (opt_debug)
    printf("[Debug] (tau) lnacceptance = %f\n", acceptance);

  if (acceptance >= 0 || legacy_rndu() < exp(acceptance))
  {
    /* accepted */
    accepted++;

    offset = 0;
    for (i = 0; i < stree->locus_count; ++i)
    {
      k = __mark_count[i];
      gnode_t ** gt_nodesptr = __gt_nodes + offset;

      for (j = 0; j < k; ++j)
        gt_nodesptr[j]->mark = 0;

      offset += __mark_count[i] + __extra_count[i];
    }

  }
  else
  {
    /* rejected */
    offset = 0;
    snode->tau = oldage;
    snode->theta = oldtheta;
    for (i = 0; i < stree->locus_count; ++i)
    {
      k = __mark_count[i];
      gnode_t ** gt_nodesptr = __gt_nodes + offset;
      double * old_ageptr = __aux + offset;

      /* restore gene tree node ages */
      for (j = 0; j < k; ++j)
        gt_nodesptr[j]->time = old_ageptr[j];

      /* restore logpr contributions */
      for (j = 0; j < 3; ++j)
        gtree_update_logprob_contrib(affected[j],loci[i]->heredity[0],i);


      /* get the list of nodes for which CLVs must be reverted, i.e. all marked
         nodes and all nodes whose left or right subtree has at least one marked
         node */
      unsigned int partials_count;
      gnode_t ** partials = gtree_return_partials(gtree[i]->root,
                                                  i,
                                                  &partials_count);

      /* revert CLV indices */
      for (j = 0; j < partials_count; ++j)
        partials[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                partials[j]->clv_index);

      /* un-mark nodes */
      for (j = 0; j < k; ++j)
        gt_nodesptr[j]->mark = 0;

      /* restore branch lengths and pmatrices */
      int matrix_updates = __mark_count[i]+__extra_count[i];
      if (matrix_updates)
      {
        if (!gt_nodesptr[0]->parent)
        {
          --matrix_updates;
          gt_nodesptr++;
        }
        if (matrix_updates)
          locus_update_matrices_jc69(loci[i],gt_nodesptr,matrix_updates);

      }

      /* restore gene tree log-likelihood */
      gtree[i]->logl = gtree[i]->old_logl;

      /* restore gene tree log probability */
      gtree[i]->logpr = gtree[i]->old_logpr;

      /* move offset to the batch of nodes for the next locus */
      offset += __mark_count[i] + __extra_count[i];
    }
  }
  return accepted;
}

double stree_propose_tau(gtree_t ** gtree, stree_t * stree, locus_t ** loci)
{
  unsigned int i;
  unsigned int candidate_count = 0;
  long accepted = 0;

  /* compute number of nodes with tau > 0 */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->tau > 0)
      candidate_count++;


  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau > 0)
      accepted += propose_tau(loci,stree->nodes[i],gtree,stree,candidate_count);
  }

  return ((double)accepted/candidate_count);
}

void stree_rootdist(stree_t * stree,
                    list_t * maplist,
                    msa_t ** msalist,
                    unsigned int ** weights)
{
  unsigned int i,j,k,n;
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

  stree->root_age = opt_tau_beta / (opt_tau_alpha-1)*4;
  printf("Root age beginning = %f\n", stree->root_age);

  if (!opt_usedata) return;

  if (opt_diploid)
  {
    for (i = 0; i < stree->tip_count; ++i)
      if (stree->nodes[i]->diploid)
        return;
  }

  hashtable_t * sht = species_hash(stree);
  hashtable_t * mht = maplist_hash(maplist,sht);

  /* find max alignment length */
  size_t maxalloc = 0;
  for (i = 0; i < msa_count; ++i)
    if ((size_t)(msalist[i]->count) > maxalloc)
      maxalloc = msalist[i]->count;

  snode_t ** pop = (snode_t **)xmalloc(maxalloc*sizeof(snode_t *));

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
      for (k = j+1; k < (unsigned int)(msalist[i]->count); ++k)
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
          diff_pair  /= msalist[i]->original_length;
          diff_locus += diff_pair;
          ++diff_count;
        }
      }
    }
    if (!diff_count) continue;

    locus_used++;

    diff_locus /= (2*diff_count);
    vd += (diff_locus-md)*(diff_locus-md) * (locus_used-1)/locus_used;
    md = (md * (locus_used-1) + diff_locus)/locus_used;
  }
  
  vd /= msa_count;
  if (locus_used >= 2)
  {
    double theta = 2*sqrt(vd);
    theta = sqrt(vd*4+1) - 1;
    theta = (2*sqrt(vd) + sqrt(vd*4+1) - 1)/2;
    if (md - theta/2 > 0)
      stree->root_age = md-theta/2;
    else
      stree->root_age = md;
  }
  else
    stree->root_age = md;

  printf("root dist = %7.5f\n", stree->root_age);
  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);
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
   Rannala, B., Yang, Z. Efficient Bayesian Species Tree Inference under the
   Multispecies Coalescent.
   Systematic Biology, 2017
*/
long stree_propose_spr(stree_t ** streeptr,
                       gtree_t *** gtree_list_ptr,
                       stree_t ** scloneptr,
                       gtree_t *** gclonesptr,
                       locus_t ** loci)
{
  unsigned int i,j,k=0;
  unsigned int branch_update_count;
  long target_count = 0;
  long source_count = 0;
  double r;
  double sum = 0;
  double lnacceptance = 0;

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
  r = legacy_rndu();
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
  if ((int)(2*legacy_rndu()) == 0)
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
  r = legacy_rndu();
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

      for (gtmp=node->parent; !(gtmp->mark & LINEAGE_OTHER); gtmp=gtmp->parent)
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

      /* experimental methods */
      double twgt = 0;
      if (opt_experimental_method)
      {
        unsigned int n;

        /* if more than one target nodes select according to likelihood */
        if (target_count > 1)
        {
          /* allocate space for target nodes weights */
          double * tweight = (double *)xmalloc(target_count * sizeof(double));

          /* compute weights for each target node in gtarget_list and store them
             in tweight */
          experimental_tselect_logl(pruned,
                                    gtarget_list,
                                    target_count,
                                    loci[i],
                                    tweight);

          /* normalize target node weights */
          for (sum=0,n=0; n < target_count; ++n)
            sum += tweight[n];

          if (sum == 0)
          {
            for (n = 0; n < target_count; ++n)
              tweight[n] = 1.0/n;
          }
          else
          {
            for (n = 0; n < target_count; ++n)
              tweight[n] /= sum;
          }

          /* randomly select a target node according to weights */
          r = legacy_rndu();
          for (sum=0,n=0; n < target_count; ++n)
          {
            sum += tweight[n];
            if (r < sum) break;
          }

          if (n == target_count)
          {
            for (n = 0; n < target_count; ++n)
              printf("%d: %f\n", n, tweight[n]);
          }

          assert(n != target_count);

          gtarget_nodes[moved_count[i]-1] = gtarget_list[n];

          twgt = tweight[n];

          /* free weights */
          free(tweight);     
        }
        else
        {
          twgt = 1;
          gtarget_nodes[moved_count[i]-1] = gtarget_list[0];
        }

      }
      else
      {
        /* randomly select a target from list */
        gtarget_nodes[moved_count[i]-1] = gtarget_list[(int)(target_count*legacy_rndu())];
      }

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
        if (stree->pptable[tmp->pop->node_index][pop_az->node_index] &&
            tmp->mark != LINEAGE_A)
        {
          gsources_list[source_count++] = tmp;
        }
      }

      if (opt_experimental_method)
      {
        unsigned int n;
        double swgt;

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
          experimental_tselect_logl(pruned,
                                    gsources_list,
                                    source_count,
                                    loci[i],
                                    tweight);

          /* normalize target node weights */
          for (sum=0,n=0; n < source_count; ++n)
            sum += tweight[n];

          if (sum == 0)
          {
            for (n = 0; n < source_count; ++n)
              tweight[n] = 1.0/n;
          }
          else
          {
            for (n = 0; n < source_count; ++n)
              tweight[n] /= sum;
          }

          /* we do not select randomly, as we already know the branch */
          gnode_t * srcnode = intact;
          while (srcnode->mark & NODE_MOVED)
            srcnode = srcnode->left->mark & LINEAGE_A ?
                        srcnode->right : srcnode->left;

          for (n = 0; n < source_count; ++n)
            if (gsources_list[n] == srcnode)
              break;
          assert(n != source_count);

          swgt = tweight[n];

          /* free weights */
          free(tweight);
        }
        else
        {
          swgt = 1;
        }

        lnacceptance += log(swgt/twgt);
      }
      else
      {
        lnacceptance += log((double)target_count / source_count);
      }
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
      assert(node==moved_nodes[j]);
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
      unlink_event(node,i);

      node->pop->event_count[i]--;
      if (!(node->pop->mark & FLAG_POP_UPDATE))
      {
        node->pop->mark |= FLAG_POP_UPDATE;
        snode_contrib[snode_contrib_count[i]++] = node->pop;
      }

      node->pop = pop_cz;
      if (!(node->pop->mark & FLAG_POP_UPDATE))
      {
        node->pop->mark |= FLAG_POP_UPDATE;
        snode_contrib[snode_contrib_count[i]++] = node->pop;
      }

      dlist_item_append(node->pop->event[i],node->event);

      node->pop->event_count[i]++;

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
        unlink_event(node,i);

        node->pop->event_count[i]--;
        if (!(node->pop->mark & FLAG_POP_UPDATE))
        {
          node->pop->mark |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        node->pop = b;
        if (!(node->pop->mark & FLAG_POP_UPDATE))
        {
          node->pop->mark |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        dlist_item_append(node->pop->event[i],node->event);

        node->pop->event_count[i]++;
      }
      else if (node->pop == c && node->time > y->tau)
      {
        /* diamond nodes */
        
        /* remove  gene node from list of coalescent events of its old population */
        unlink_event(node,i);

        node->pop->event_count[i]--;
        if (!(node->pop->mark & FLAG_POP_UPDATE))
        {
          node->pop->mark |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        node->pop = y;
        if (!(node->pop->mark & FLAG_POP_UPDATE))
        {
          node->pop->mark |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        dlist_item_append(node->pop->event[i],node->event);

        node->pop->event_count[i]++;
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

        unlink_event(node,i);

        node->pop->event_count[i]--;
        if (!(node->pop->mark & FLAG_POP_UPDATE))
        {
          node->pop->mark |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }
        
        if (pop == c)
          node->pop = y;
        else
          node->pop = pop;

        if (!(node->pop->mark & FLAG_POP_UPDATE))
        {
          node->pop->mark |= FLAG_POP_UPDATE;
          snode_contrib[snode_contrib_count[i]++] = node->pop;
        }

        dlist_item_append(node->pop->event[i],node->event);

        node->pop->event_count[i]++;
      }
    }

    /* Flag populations Y,C and B for re-computing their contribution to the i-th gene tree probability only if:
         (a) they have not been already flagged in a previous step
         (b) there is more than one outgoing lineages (entering its parent population).
    */
    if (!(y->mark & FLAG_POP_UPDATE) && (y->seqin_count[i]-y->event_count[i]>1))
      snode_contrib[snode_contrib_count[i]++] = y;
    if (!(c->mark & FLAG_POP_UPDATE) && (c->seqin_count[i]-c->event_count[i]>1))
      snode_contrib[snode_contrib_count[i]++] = c;
    if (!(b->mark & FLAG_POP_UPDATE) && (b->seqin_count[i]-b->event_count[i]>1))
      snode_contrib[snode_contrib_count[i]++] = b;

    moved_nodes += gtree->inner_count;
    gtarget_nodes += gtree->inner_count;
    gtarget_list += gtree->tip_count + gtree->inner_count;

    __mark_count[i] = branch_update_count;
    bl_list += branch_update_count;
    snode_contrib += stree->tip_count + stree->inner_count;

    /* reset species tree marks */
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      stree->nodes[j]->mark = 0;
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

  for (i = 0; i < stree->locus_count; ++i)
    fill_seqin_counts(stree,i);

  reset_gene_leaves_count(stree);
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
  for (i=0; i < stree->locus_count; ++i)
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
      
      locus_update_matrices_jc69(loci[i],bl_list,__mark_count[i]);
      
      /* retrieve all nodes whose partials must be updates */
      unsigned int partials_count;
      gnode_t ** partials = gtree_return_partials(gtree_list[i]->root,
                                                  i,
                                                  &partials_count);
      
      /* point to the double-buffered partials space */
      for (j = 0; j < partials_count; ++j)
        partials[j]->clv_index = SWAP_CLV_INDEX(gtree_list[i]->tip_count,
                                                partials[j]->clv_index);

      /* update conditional probabilities (partials) of affected nodes */
      locus_update_partials(loci[i],partials,partials_count);

      /* evaluate log-likelihood */
      unsigned int param_indices[1] = {0};
      gtree_list[i]->logl = locus_root_loglikelihood(loci[i],
                                                     gtree_list[i]->root,
                                                     param_indices,
                                                     NULL);
    }

    gtree_list[i]->old_logpr = gtree_list[i]->logpr;
    #if 0
    /* This recomputes the gene tree probabilities from scratch. It can be used to verify that
       the code below, which only computes the gene tree probability for the changed components,
       is correct. */
    
      double logpr = gtree_logprob(stree,loci[i]->heredity[0],i);
    #else

    /* locate additional populations that need to be updated */

    /* find and mark those populations whose number of incoming lineages has
       changed due to the reset_gene_leaves_count() call, but were previously
       not marked for log-probability contribution update */
    for (j = 0; j < snode_contrib_count[i]; ++j)
      snode_contrib[j]->mark |= FLAG_POP_UPDATE;
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
    {
      snode_t * snode = stree->nodes[j];
      if (!(snode->mark & FLAG_POP_UPDATE) &&
          (snode->seqin_count[i] != original_stree->nodes[j]->seqin_count[i]))
        snode_contrib[snode_contrib_count[i]++] = snode;
    }

    /* now update the log-probability contributions for the affected, marked
       populations */
    for (j = 0; j < snode_contrib_count[i]; ++j)
    {
      gtree_list[i]->logpr -= snode_contrib[j]->logpr_contrib[i];
      gtree_update_logprob_contrib(snode_contrib[j],loci[i]->heredity[0],i);
      gtree_list[i]->logpr += snode_contrib[j]->logpr_contrib[i];
    }

    /* reset markings on affected populations */
    for (j = 0; j < snode_contrib_count[i]; ++j)
      snode_contrib[j]->mark = 0;
    #endif


    bl_list += __mark_count[i];

    for (j = 0; j < gtree_list[i]->tip_count + gtree_list[i]->inner_count; ++j)
      gtree_list[i]->nodes[j]->mark = 0;

    snode_contrib += stree->tip_count + stree->inner_count;
    
    lnacceptance += gtree_list[i]->logpr - gtree_list[i]->old_logpr +
                    gtree_list[i]->logl  - gtree_list[i]->old_logl;
  }

  if (opt_debug)
    printf("[Debug] (SSPR) lnacceptance = %f\n", lnacceptance);

  /* in case of acceptance, cloned trees are substituted with the original ones,
     and species tree nodes are re-labeled, but all this is done in method_01.c
  */
  return (lnacceptance >= 0 || legacy_rndu() < exp(lnacceptance));
}
