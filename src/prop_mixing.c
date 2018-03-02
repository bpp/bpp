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
  unsigned int theta_count=0;
  unsigned int tau_count=0;
  double lnc,c;
  double logpr = 0;
  double lnacceptance;
  long accepted = 0;

  double * notheta_old_logpr = NULL;

  size_t nodes_count = stree->tip_count+stree->inner_count; 


  if (!opt_est_theta)
  {
    notheta_old_logpr = (double *)xmalloc(nodes_count * sizeof(double));
    for (i = 0; i < nodes_count; ++i)
      notheta_old_logpr[i] = stree->nodes[i]->notheta_logpr_contrib;
  }
  
  /* TODO: Account for method 11 / rj-MCMC */
  if (opt_est_theta)
  {
    /* TODO: Precompute how many theta parameters we have for A00 */
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      if (stree->nodes[i]->theta > 0)
        theta_count++;
  }
  else
    theta_count = 0;

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->tau > 0)
      tau_count++;

  lnc = opt_finetune_mix * legacy_rnd_symmetrical();
  c = exp(lnc);

  /* sum of inner nodes for all loci */
  for (i=0,k=0; i < stree->locus_count; ++i)
    k += gtree[i]->inner_count; 

  lnacceptance = (theta_count + tau_count + k)*lnc;

  /* TODO: skip this for integrated-out theta */
  /* TODO: This loop separation is for having the same traversal as old bpp */

  /* TODO: Why is this allocation here? Perhaps no longer needed? */
  snode_t ** snodes = (snode_t **)xmalloc((stree->tip_count+stree->inner_count)*
                                          sizeof(snode_t *));
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    snodes[i] = stree->nodes[i];

  /* change the thetas */
  if (opt_est_theta)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      if (snodes[i]->theta <= 0) continue;

      snodes[i]->old_theta = snodes[i]->theta;
      snodes[i]->theta *= c;
      lnacceptance += (-opt_theta_alpha-1)*lnc -
                     opt_theta_beta*(1/snodes[i]->theta-1/snodes[i]->old_theta);
    }
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
  
  if (!opt_est_theta)
    logpr = stree->notheta_logpr;

  for (i = 0; i < stree->locus_count; ++i)
  {
    gtree_t * gt = gtree[i];

    /* go through all gene nodes */
    for (j = gt->tip_count; j < gt->tip_count + gt->inner_count; ++j)
    {
      gt->nodes[j]->old_time = gt->nodes[j]->time;
      gt->nodes[j]->time *= c;
    }

    /* update branch lengths */
    for (j = 0; j < gt->tip_count+gt->inner_count; ++j)
    {
      if (gt->nodes[j]->parent)
        gt->nodes[j]->length = gt->nodes[j]->parent->time - gt->nodes[j]->time;
    }

    /* update pmatrices */
    /* TODO: Remove this allocation */
    gnode_t ** gt_nodes = (gnode_t **)xmalloc((gt->tip_count + gt->inner_count)*
                                              sizeof(gnode_t *));
    k=0;
    for (j = 0; j < gt->tip_count + gt->inner_count; ++j)
      if (gt->nodes[j]->parent)
        gt_nodes[k++] = gt->nodes[j];
    locus_update_matrices_jc69(locus[i],gt_nodes,k);

    gtree_all_partials(gt->root,gt_nodes,&k);
    for (j = 0; j < k; ++j)
      gt_nodes[j]->clv_index = SWAP_CLV_INDEX(gt->tip_count,gt_nodes[j]->clv_index);

    locus_update_partials(locus[i],gt_nodes,k);

    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    double logl = locus_root_loglikelihood(locus[i],gt->root,param_indices,NULL);


    if (opt_est_theta)
      logpr = gtree_logprob(stree,locus[i]->heredity[0],i);
    else
    {
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      {
        logpr -= stree->nodes[j]->notheta_logpr_contrib;
        logpr += gtree_update_logprob_contrib(stree->nodes[j],locus[i]->heredity[0],i);
      }
    }
        

    if (opt_est_theta)
      lnacceptance += logl - gt->logl + logpr - gt->logpr;
    else
      lnacceptance += logl - gt->logl;

    if (opt_est_theta)
    {
      gt->old_logpr = gt->logpr;
      gt->logpr = logpr;
    }

    gt->old_logl = gt->logl;
    gt->logl = logl;

    free(gt_nodes);
    
  }

  if (!opt_est_theta)
    lnacceptance += logpr - stree->notheta_logpr;

  if (opt_debug)
    printf("[Debug] (mixing) lnacceptance = %f\n", lnacceptance);


  if (lnacceptance >= -1e-10 || legacy_rndu() < exp(lnacceptance))
  {
    /* accept */
    accepted = 1;

    if (!opt_est_theta)
      stree->notheta_logpr = logpr;
  }
  else
  {
    /* revert thetas and logpr contributions */
    if (opt_est_theta)
    {
      for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      {
        for (j = 0; j < stree->locus_count; ++j)
          snodes[i]->logpr_contrib[j] = snodes[i]->old_logpr_contrib[j];

        if (snodes[i]->theta <= 0) continue;

        /* TODO: Note that, it is both faster and more precise to restore the old
           value from memory, than re-computing it with a division. For now, we
           use the division here to be compatible with the old bpp */
#if 0
        snodes[i]->theta = snodes[i]->old_theta;
#else
        snodes[i]->theta /= c;
#endif
      }
    }
    else
    {
      for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      {
        for (j = 0; j < opt_locus_count; ++j)
          logprob_revert_notheta(stree->nodes[i],j);
        stree->nodes[i]->notheta_logpr_contrib = notheta_old_logpr[i];
      }
    }

    /* revert taus */
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    {
      if (snodes[i]->tau == 0) continue;

      /* TODO: Same here */
#if 0
      snodes[i]->tau = snodes[i]->old_tau;
#else
      snodes[i]->tau /= c;
#endif
    }

    /* go through all loci */
    for (i = 0; i < stree->locus_count; ++i)
    {
      /* restore logl and logpr */
      gtree[i]->logl  = gtree[i]->old_logl;
      if (opt_est_theta)
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
  free(snodes);

  if (!opt_est_theta)
    free(notheta_old_logpr);

  return accepted;

}
