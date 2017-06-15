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

static int cb_cmp_pairlabel(void * a, void * b)
{
  pair_t * pair = (pair_t *)a;
  char * label  = (char *)b;

  return (!strcmp(pair->label,label));
}

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
      fatal("Duplicate taxon or node label (%s) in file %s",
            tree->nodes[i]->label,
            opt_streefile);
    }
  }

  return ht;
}

static void stree_label(stree_t * stree, int msa_count)
{
  unsigned int i;

  /* labels at inner nodes are concatenated labels of their children */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    /* get the species structures of the two children */
    snode_t * node  = stree->nodes[i];
    snode_t * lnode = stree->nodes[i]->left;
    snode_t * rnode = stree->nodes[i]->right;

    /* allocate necessary memory for label */
    if (node->label)
      free(node->label);
    node->label = (char *)xmalloc(strlen(lnode->label)+
                                  strlen(rnode->label)+1);

    /* concatenate species labels */
    node->label[0] = 0;
    strcat(node->label,lnode->label);
    strcat(node->label,rnode->label);
  }
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

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
        fatal("Cannot find");

      node = (snode_t *)(pair->data);

#ifdef DEBUG_STREE_INIT
      printf("Matched %s -> %s\n", label, node->label);
#endif

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

  node->tau = tau_parent * (prop + (1 - prop - 0.02)*legacy_rndu());

  stree_init_tau_recursive(node->left,prop);
  stree_init_tau_recursive(node->right,prop);
}

static void stree_init_tau(stree_t * stree)
{
  /* Initialize speciation times for each extinct species */ 

  double prop = (stree->root->leaves > PROP_THRESHOLD) ? 0.9 : 0.5;

  /* set the speciation time for root */
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
  /* TODO: Note, that we compute the theta for the root first to be in line
     with the original bpp code */
  stree->root->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                       (0.9 + 0.2 * legacy_rndu());
  ++k;
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count-1; ++i)
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

void stree_init(stree_t * stree, msa_t ** msa, list_t * maplist, int msa_count)
{
  unsigned int i,j;
  snode_t * curnode;
  snode_t * ancnode;

  assert(msa_count > 0);

  /* label each inner node of the species tree with the concatenated labels of
     its two children */
  stree_label(stree, msa_count);

  /* allocate space for keeping track of coalescent events at each species tree
     node for each locus */
  stree->locus_count = (unsigned int)msa_count;
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    stree->nodes[i]->event = (dlist_t **)xcalloc(msa_count,sizeof(dlist_t *));
    stree->nodes[i]->event_count = (int *)xcalloc(msa_count,sizeof(int));
    stree->nodes[i]->seqin_count = (int *)xcalloc(msa_count,sizeof(int));
    stree->nodes[i]->gene_leaves = (unsigned int *)xcalloc(msa_count,sizeof(unsigned int));
    stree->nodes[i]->logpr_contrib = (double *)xcalloc(msa_count,sizeof(double));
    stree->nodes[i]->old_logpr_contrib = (double *)xcalloc(msa_count,sizeof(double));

    for (j = 0; j < stree->locus_count; ++j)
      stree->nodes[i]->event[j] = dlist_create();
  }

  /* Initialize population sizes */
  stree_init_theta(stree, msa, maplist, msa_count);

  /* Initialize speciation times and create extinct species groups */
  stree_init_tau(stree);

  /* initialize pptable, where pptable[i][j] indicates whether population j
     (that is, node with index j) is ancestral to population i */
  stree->pptable = (int **)xcalloc((stree->tip_count + stree->inner_count),
                                   sizeof(int *));
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    stree->pptable[i] = (int *)xcalloc((stree->tip_count + stree->inner_count),
                                       sizeof(int));

  
  for (i = 0; i < stree->tip_count; ++i)
  {
    //curnode = stree->nodes[i];
    //ancnode = stree->nodes[i]->parent;

    for (curnode = stree->nodes[i]; curnode; curnode = curnode->parent)
      for (ancnode = curnode; ancnode; ancnode = ancnode->parent)
        stree->pptable[curnode->node_index][ancnode->node_index] = 1;
  }
}

static int propose_theta(gtree_t ** gtree, int locus_count, snode_t * snode)
{
  int i;
  double thetaold;
  double thetanew;
  const double finetune = 0.001;
  double acceptance;

  thetaold = snode->theta;
  
  thetanew = thetaold + finetune * legacy_rnd_symmetrical();
 
  if (thetanew < 0)
    thetanew = -thetanew;

  snode->theta = thetanew;

  acceptance = (-opt_theta_alpha-1) * log(thetanew/thetaold) -
               opt_theta_beta*(1/thetanew - 1/thetaold);

  for (i = 0; i < locus_count; ++i)
  {
    /* save a copy of old logpr */
    gtree[i]->old_logpr = gtree[i]->logpr;

    gtree[i]->logpr -= snode->logpr_contrib[i];
    gtree_update_logprob_contrib(snode,i);
    gtree[i]->logpr += snode->logpr_contrib[i];
    
    acceptance += (gtree[i]->logpr - gtree[i]->old_logpr);

  }

//#ifdef DEBUG_PROPOSAL_THETA
#if 1
  printf("THETAlnacceptance = %f\n", acceptance);
#endif

  if (acceptance >= 0 || legacy_rndu() < exp(acceptance))
  {
    return 1;
  }

  /* TODO: Since the rejections are higher than acceptances, it would
     make sense not to update gtree[i]->logpr in the first place, but
     only update it when proposal is accepted */

  /* reject */
  for (i = 0; i < locus_count; ++i)
    gtree[i]->logpr = gtree[i]->old_logpr;

  snode->theta = thetaold;
  for (i = 0; i < locus_count; ++i)
    snode->logpr_contrib[i] = snode->old_logpr_contrib[i];

  return 0;
}

void stree_propose_theta(gtree_t ** gtree, stree_t * stree)
{
  unsigned int i;
  int accepted = 0;
  int theta_count = 0;
  snode_t * snode;

  /* TODO: this kind of loop is for backwards compatibility with the old bpp
     since the root node is the first inner node and not the last one */
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode = stree->nodes[i];
    if (snode->theta)
    {
      accepted += propose_theta(gtree, stree->locus_count, stree->nodes[i]);
      theta_count++;
    }
  }

  if (stree->root->theta)
  {
    accepted += propose_theta(gtree, stree->locus_count, stree->root);
    theta_count++;
  }
    
  for (i = 0; i < stree->inner_count-1; ++i)
  {
    snode = stree->nodes[stree->tip_count+i];
    if (snode->theta)
    {
      accepted += propose_theta(gtree, stree->locus_count, snode);
      theta_count++;
    }
  }
}
