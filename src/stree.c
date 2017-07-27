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

static gnode_t ** __gt_nodes = NULL;
static double * __aux = NULL;
static int * __mark_count = NULL;
static int * __extra_count = NULL;

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

static void stree_label(stree_t * stree)
{
  stree_label_recursive(stree->root);
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

void stree_init(stree_t * stree, msa_t ** msa, list_t * maplist, int msa_count)
{
  unsigned int i,j;
  snode_t * curnode;
  snode_t * ancnode;

  assert(msa_count > 0);

  /* label each inner node of the species tree with the concatenated labels of
     its two children */
  stree_label(stree);

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
    for (curnode = stree->nodes[i]; curnode; curnode = curnode->parent)
      for (ancnode = curnode; ancnode; ancnode = ancnode->parent)
        stree->pptable[curnode->node_index][ancnode->node_index] = 1;
  }
  

  /* allocate traversal buffer to be the size of all nodes for all loci */
  unsigned int sum_count = 0;
  for (i = 0; i < (unsigned int)msa_count; ++i)
    sum_count += (unsigned int)(msa[i]->count);

  sum_count = 2*sum_count - msa_count;
  __gt_nodes = (gnode_t **)xmalloc(sum_count * sizeof(gnode_t *));
  __aux = (double *)xmalloc((sum_count - msa_count)*sizeof(double *));

  /* The following two arrays are used purely for the tau proposal.
     Entry i of marked_count indicates how many nodes from locus i are marked.
     Similarly, extra_count 
     how many extra nodes where added in __gt_nodes whose branch lengths (and
     therefore) p-matrices need updating because their parent node's age was
     changed */
  __mark_count  = (int *)xmalloc(msa_count*sizeof(int));
  __extra_count = (int *)xmalloc(msa_count*sizeof(int));
}

void stree_fini()
{
  free(__gt_nodes);
  free(__aux);
  free(__mark_count);
  free(__extra_count);
}

static int propose_theta(gtree_t ** gtree, int locus_count, snode_t * snode)
{
  int i;
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

  for (i = 0; i < locus_count; ++i)
  {
    /* save a copy of old logpr */
    gtree[i]->old_logpr = gtree[i]->logpr;

    gtree[i]->logpr -= snode->logpr_contrib[i];
    gtree_update_logprob_contrib(snode,i);
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
  for (i = 0; i < locus_count; ++i)
    gtree[i]->logpr = gtree[i]->old_logpr;

  snode->theta = thetaold;
  for (i = 0; i < locus_count; ++i)
    snode->logpr_contrib[i] = snode->old_logpr_contrib[i];

  return 0;
}

double stree_propose_theta(gtree_t ** gtree, stree_t * stree)
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
      accepted += propose_theta(gtree, stree->locus_count, stree->nodes[i]);
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
          node->mark = 1;
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
        logpr += gtree_update_logprob_contrib(affected[j],i);
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
        gtree_update_logprob_contrib(affected[j],i);


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
