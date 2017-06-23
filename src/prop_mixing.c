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

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))

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

long proposal_mixing(gtree_t ** gtree, stree_t * stree, locus_t ** locus)
{
  unsigned i,j,k;
  unsigned int theta_count, tau_count;
  double lnc,c;
  double finetune = 0.3;
  double lnacceptance;
  long accepted = 0;


  /* TODO: Account for integrated out theta, and also for rj-MCMC */
  theta_count = stree->tip_count + stree->inner_count;

  /* TODO: This will not be the same for rj-MCMC where collapsed species tree
   * nodes will have a tau=0 */
  tau_count   = stree->inner_count;

  lnc = finetune * legacy_rnd_symmetrical();
  c = exp(lnc);

  /* sum of inner nodes for all loci */
  for (i=0,k=0; i < stree->locus_count; ++i)
    k += gtree[i]->inner_count; 

  lnacceptance = (theta_count + tau_count + k)*lnc;

  /* TODO: skip this for integrated-out theta */
  /* TODO: This loop separation is for having the same traversal as old bpp */
  snode_t ** snodes = (snode_t **)xmalloc((stree->tip_count+stree->inner_count)*
                                          sizeof(snode_t *));
  for (i = 0; i < stree->tip_count; ++i) snodes[i] = stree->nodes[i];
  snodes[i++] = stree->root;
  for (j = stree->tip_count; j < stree->tip_count + stree->inner_count-1; ++j)
    snodes[i++] = stree->nodes[j];

  /* change the thetas */
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (snodes[i]->theta <= 0) continue;

    snodes[i]->old_theta = snodes[i]->theta;
    snodes[i]->theta *= c;
    lnacceptance += (-opt_theta_alpha-1)*lnc -
                    opt_theta_beta*(1/snodes[i]->theta - 1/snodes[i]->old_theta);
  }

  /* change the taus */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (snodes[i]->tau == 0) continue;
    snodes[i]->old_tau = snodes[i]->tau;
    snodes[i]->tau *= c;
    if (!snodes[i]->parent)
    {
      lnacceptance += (-opt_tau_alpha-1 - tau_count+1)*lnc -
                      opt_tau_beta*(1/snodes[i]->tau - 1/snodes[i]->old_tau);
    }
  }
  
  for (i = 0; i < stree->locus_count; ++i)
  {
    /* go through all gene nodes */
    for (j = gtree[i]->tip_count; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      gtree[i]->nodes[j]->old_time = gtree[i]->nodes[j]->time;
      gtree[i]->nodes[j]->time *= c;
    }

    /* update branch lengths */
    for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      if (gtree[i]->nodes[j]->parent)
        gtree[i]->nodes[j]->length = gtree[i]->nodes[j]->parent->time - gtree[i]->nodes[j]->time;
    }

    /* update pmatrices */
    /* TODO: Remove this allocation */
    gnode_t ** gt_nodes = (gnode_t **)xmalloc((gtree[i]->tip_count + gtree[i]->inner_count) * sizeof(gnode_t *));
    k=0;
    for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      if (gtree[i]->nodes[j]->parent)
        gt_nodes[k++] = gtree[i]->nodes[j];
    locus_update_matrices_jc69(locus[i],gt_nodes,k);

    gtree_all_partials(gtree[i]->root,gt_nodes,&k);
    for (j = 0; j < k; ++j)
      gt_nodes[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,gt_nodes[j]->clv_index);

    locus_update_partials(locus[i],gt_nodes,k);

    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    double logl = locus_root_loglikelihood(locus[i],gtree[i]->root,param_indices,NULL);


    double logpr = gtree_logprob(stree,i);
    if (opt_debug > 1)
    {
      printf("old logl: %f\n", gtree[i]->logl);
      printf("new logl: %f\n", logl);
      printf("old logpr: %f\n", gtree[i]->logpr);
      printf("new logpr: %f\n", logpr);
    }

    lnacceptance += logl - gtree[i]->logl + logpr - gtree[i]->logpr;

    gtree[i]->old_logpr = gtree[i]->logpr;
    gtree[i]->old_logl = gtree[i]->logl;
    gtree[i]->logpr = logpr;
    gtree[i]->logl = logl;

    free(gt_nodes);
    
  }

  if (opt_debug)
    printf("[Debug] (mixing) lnacceptance = %f\n", lnacceptance);


  if (lnacceptance >= 0 || legacy_rndu() < exp(lnacceptance))
  {
    /* accept */
    accepted = 1;
  }
  else
  {

    /* revert thetas and logpr contributions */
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      for (j = 0; j < stree->locus_count; ++j)
        snodes[i]->logpr_contrib[j] = snodes[i]->old_logpr_contrib[j];

      if (snodes[i]->theta <= 0) continue;
      snodes[i]->theta = snodes[i]->old_theta;
    }

    /* revert taus */
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    {
      if (snodes[i]->tau == 0) continue;
      snodes[i]->tau = snodes[i]->old_tau;
    }

    /* go through all loci */
    for (i = 0; i < stree->locus_count; ++i)
    {
      /* restore logl and logpr */
      gtree[i]->logl  = gtree[i]->old_logl;
      gtree[i]->logpr = gtree[i]->old_logpr;

      gnode_t ** gnodeptr = gtree[i]->nodes;
      /* revert CLV indices and coalescent event ages */
      for (j = gtree[i]->tip_count; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
      {
        gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                gnodeptr[j]->clv_index);
        gnodeptr[j]->time = gnodeptr[j]->old_time;
      }

      
      gnode_t ** gt_nodes = (gnode_t **)xmalloc((gtree[i]->tip_count +
                                                gtree[i]->inner_count) * 
                                                sizeof(gnode_t *));

      /* revert branch lengths and p-matrices */
      k=0;
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        if (gtree[i]->nodes[j]->parent)
          gt_nodes[k++] = gtree[i]->nodes[j];
      locus_update_matrices_jc69(locus[i],gt_nodes,k);

      free(gt_nodes);
      
    }
  }

  return accepted;

}
