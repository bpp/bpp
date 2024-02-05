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

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))
#define SWAP_SCALER_INDEX(n,i) (((n)+((i)-1))%(2*(n)-2)) 
#define SWAP_PMAT_INDEX(e,i) (((e)+(i))%((e)<<1))

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

void prop_mixing_update_gtrees(locus_t ** locus,
                               gtree_t ** gtree,
                               stree_t * stree,
                               long locus_start,
                               long locus_count,
                               double c,
                               long thread_index,
                               double * ret_lnacceptance,
                               double * ret_logpr)
{
  long i;
  unsigned int j,k;
  double lnacceptance = 0;
  double logpr = 0;
  size_t nodes_count = stree->tip_count+stree->inner_count+stree->hybrid_count; 

  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    gtree_t * gt = gtree[i];

    /* go through all gene nodes */
    for (j = gt->tip_count; j < gt->tip_count + gt->inner_count; ++j)
    {
      gt->nodes[j]->old_time = gt->nodes[j]->time;
      gt->nodes[j]->time *= c;
    }

    /* migration */
    if (opt_migration)
    {
      for (j = 0; j < gt->tip_count + gt->inner_count; ++j)
      {
        if (!gt->nodes[j]->mi || !gt->nodes[j]->mi->count) continue;
        miginfo_t * mi = gt->nodes[j]->mi;

        for (k = 0; k < mi->count; ++k)
        {
          mi->me[k].old_time = mi->me[k].time;
          mi->me[k].time *= c;
        }
      }
    }

    /* update branch lengths */
    for (j = 0; j < gt->tip_count+gt->inner_count; ++j)
    {
      if (gt->nodes[j]->parent)
      {
        /* TODO: 28.1.2019 */
        gt->nodes[j]->pmatrix_index = SWAP_PMAT_INDEX(gt->edge_count,gt->nodes[j]->pmatrix_index);
        /* gt->nodes[j]->length = gt->nodes[j]->parent->time - gt->nodes[j]->time; */
      }
    }

    /* update pmatrices */
    /* TODO: Remove this allocation */
    gnode_t ** gt_nodes = (gnode_t **)xmalloc((gt->tip_count + gt->inner_count)*
                                              sizeof(gnode_t *));
    k=0;
    for (j = 0; j < gt->tip_count + gt->inner_count; ++j)
      if (gt->nodes[j]->parent)
        gt_nodes[k++] = gt->nodes[j];
    locus_update_matrices(locus[i],gtree[i],gt_nodes,stree,i,k);

    gtree_all_partials(gt->root,gt_nodes,&k);
    for (j = 0; j < k; ++j)
    {
      gt_nodes[j]->clv_index = SWAP_CLV_INDEX(gt->tip_count,gt_nodes[j]->clv_index);
      if (opt_scaling)
        gt_nodes[j]->scaler_index = SWAP_SCALER_INDEX(gt->tip_count,
                                                      gt_nodes[j]->scaler_index);
    }

    locus_update_partials(locus[i],gt_nodes,k);

    /* compute log-likelihood */
    double logl = locus_root_loglikelihood(locus[i],gt->root,locus[i]->param_indices,NULL);


    if (opt_est_theta)
    {
      if (opt_migration)
      {
        logpr = gtree_logprob_mig(stree,
                                  gtree[i],
                                  locus[i]->heredity[0],
                                  i,
                                  thread_index);
      }
      else
        logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index);
    }
    else
    {
      if (opt_migration)
        fatal("Integrating out thetas for IM model not implemented yet");

      for (j = 0; j < nodes_count; ++j)
      {
        logpr -= stree->nodes[j]->notheta_logpr_contrib;
        logpr += gtree_update_logprob_contrib(stree->nodes[j],
                                              locus[i]->heredity[0],
                                              i,
                                              thread_index);
      }
    }

    if (opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
    {
      double new_prior_rates = lnprior_rates(gtree[i], stree, i);
      lnacceptance += new_prior_rates - gtree[i]->lnprior_rates;
      gtree[i]->old_lnprior_rates = gtree[i]->lnprior_rates;
      gtree[i]->lnprior_rates = new_prior_rates;
      #if 0
      assert(new_prior_rates > gtree[i]->old_lnprior_rates - PLL_MISC_EPSILON &&
             new_prior_rates < gtree[i]->old_lnprior_rates + PLL_MISC_EPSILON);
      #endif
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
  /* return values */
  *ret_lnacceptance = lnacceptance;
  if (!opt_est_theta)
    *ret_logpr = logpr;
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

  const long thread_index = 0;

  double * notheta_old_logpr = NULL;

  size_t nodes_count = stree->tip_count+stree->inner_count+stree->hybrid_count; 


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
    for (i = 0; i < nodes_count; ++i)
      if (stree->nodes[i]->theta > 0 && !stree->nodes[i]->linked_theta)
        theta_count++;
  }
  else
    theta_count = 0;

  if (opt_msci)
  {
    for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    {
      if (stree->nodes[i]->tau > 0 && stree->nodes[i]->prop_tau)
        tau_count++;
    }
  }
  else
  {
    for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
      if (stree->nodes[i]->tau > 0)
        tau_count++;
  }

  lnc = opt_finetune_mix * legacy_rnd_symmetrical(thread_index);
  c = exp(lnc);

  /* sum of inner nodes for all loci */
  for (i=0,k=0; i < stree->locus_count; ++i)
    k += gtree[i]->inner_count; 

  lnacceptance = (theta_count + tau_count + k)*lnc;

  /* account for migration events */
  if (opt_migration)
  {
    for (i=0,k=0; i < stree->locus_count; ++i)
      for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
        if (gtree[i]->nodes[j]->mi)
          k += gtree[i]->nodes[j]->mi->count;
    lnacceptance += k*lnc;
  }

  /* TODO: skip this for integrated-out theta */
  /* TODO: This loop separation is for having the same traversal as old bpp */

  /* TODO: Why is this allocation here? Perhaps no longer needed? */
  snode_t ** snodes = (snode_t **)xmalloc((nodes_count)*
                                          sizeof(snode_t *));
  for (i = 0; i < nodes_count; ++i)
    snodes[i] = stree->nodes[i];

  /* change the thetas */
  if (opt_est_theta)
  {
    for (i = 0; i < nodes_count; ++i)
    {
      if (snodes[i]->theta <= 0) continue;

      snodes[i]->old_theta = snodes[i]->theta;
      snodes[i]->theta *= c;

      if (snodes[i]->linked_theta) continue;
      
      if (opt_theta_prior == BPP_THETA_PRIOR_INVGAMMA)
        lnacceptance += (-opt_theta_alpha-1)*lnc -
                     opt_theta_beta*(1/snodes[i]->theta-1/snodes[i]->old_theta);
      else if (opt_theta_prior == BPP_THETA_PRIOR_GAMMA)
        lnacceptance += (opt_theta_alpha-1)*lnc -
                       opt_theta_beta*(snodes[i]->theta - snodes[i]->old_theta);
    }
  }

  /* change the taus */
  if (opt_msci)
  {
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    {
      if (snodes[i]->tau == 0 || snodes[i]->prop_tau == 0) continue;
      snodes[i]->old_tau = snodes[i]->tau;
      snodes[i]->tau *= c;
      if (snodes[i]->hybrid)
      {
        if (node_is_hybridization(snodes[i]))
        {
          assert(!node_is_mirror(snodes[i]));
          snodes[i]->hybrid->tau = snodes[i]->tau;
          if (snodes[i]->htau == 0)
            snodes[i]->parent->tau = snodes[i]->tau;
          if (snodes[i]->hybrid->htau == 0)
            snodes[i]->hybrid->parent->tau = snodes[i]->tau;
        }
        else
        {
          assert(node_is_bidirection(snodes[i]));
          assert(!node_is_mirror(snodes[i]));

          assert(snodes[i]->htau && !snodes[i]->hybrid->htau &&
                 snodes[i]->hybrid->parent->htau && !snodes[i]->right->htau);
          assert(snodes[i]->prop_tau && !snodes[i]->hybrid->parent->prop_tau);
          snodes[i]->hybrid->tau        = snodes[i]->tau;
          snodes[i]->right->tau         = snodes[i]->tau;
          snodes[i]->right->hybrid->tau = snodes[i]->tau;
          
        }
      }
      if (!snodes[i]->parent)
      {
        if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
        {
          /* inv-gamma prior */
          lnacceptance += (-opt_tau_alpha-1 - tau_count+1)*lnc -
                          opt_tau_beta*(1/snodes[i]->tau - 1/snodes[i]->old_tau);
        }
        else
        {
          /* gamma prior */
          lnacceptance += (opt_tau_alpha-1 - tau_count+1)*lnc -
                          opt_tau_beta*(snodes[i]->tau - snodes[i]->old_tau);
        }
      }
    }
  }
  else
  {
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    {
      if (snodes[i]->tau == 0) continue;
      snodes[i]->old_tau = snodes[i]->tau;
      snodes[i]->tau *= c;
      if (!snodes[i]->parent)
      {
        if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
        {
          /* inv-gamma prior */
          lnacceptance += (-opt_tau_alpha-1 - tau_count+1)*lnc -
                          opt_tau_beta*(1/snodes[i]->tau - 1/snodes[i]->old_tau);
        }
        else
        {
          /* gamma prior */
          lnacceptance += (opt_tau_alpha-1 - tau_count+1)*lnc -
                          opt_tau_beta*(snodes[i]->tau - snodes[i]->old_tau);
        }
      }
    }
  }
  
  if (!opt_est_theta)
    logpr = stree->notheta_logpr;

  if (opt_migration)
  {
    #if 1
    stree_update_mig_subpops(stree, thread_index);
    #else
    /* TODO: Note, the following loop should be a correct and faster alternative
       than calling the stree_update_mig_subpops function. However, those
       multiplications cause the migbuffer[j].time to differ numerically with
       the taus in the species tree. We need to include an snode pointer in the
       structure and copy the values instead of multiplying */
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    {
      snode_t * x = stree->nodes[i];
      for (j = 0; j < x->mb_count; ++j)
        x->migbuffer[j].time *= c;
    }
    #endif
  }

  /* update gene trees with either parallel or serial code. Note, for thetas
     integrated out we use the serial version */
  if (opt_est_theta && opt_threads > 1)
  {
    /* TODO: It seems logpr_change is not used with multiple threads and 
       integrated out thetas. This must be fixed. */
    assert(opt_est_theta);

    thread_data_t td;
    td.locus = locus; td.gtree = gtree; td.stree = stree;
    td.c = c;
    threads_wakeup(THREAD_WORK_MIXING,&td);
    lnacceptance += td.lnacceptance;
  }
  else
  {
    double lnacc_contrib = 0;
    double logpr_change = 0;
    prop_mixing_update_gtrees(locus,
                              gtree,
                              stree,
                              0,
                              stree->locus_count,
                              c,
                              0,
                              &lnacc_contrib,
                              &logpr_change);
    lnacceptance += lnacc_contrib;
    if (!opt_est_theta)
      logpr += logpr_change;
                              
  }

  #if 0
  if (opt_est_theta && !(opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL))
    assert(logpr == 0);
  #endif

  if (!opt_est_theta)
    lnacceptance += logpr - stree->notheta_logpr;

  if (opt_debug_mix)
    printf("[Debug] (mixing) lnacceptance = %f\n", lnacceptance);


  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
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
      for (i = 0; i < nodes_count; ++i)
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
      for (i = 0; i < nodes_count; ++i)
      {
        for (j = 0; j < opt_locus_count; ++j)
          logprob_revert_notheta(stree->nodes[i],j);
        stree->nodes[i]->notheta_logpr_contrib = notheta_old_logpr[i];
      }
    }

    /* revert taus */
    if (opt_msci)
    {
      for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
      {
        if (snodes[i]->tau == 0 || snodes[i]->prop_tau == 0) continue;

        snodes[i]->tau /= c;
        if (snodes[i]->hybrid)
        {
          
          if (node_is_hybridization(snodes[i]))
          {
            assert(!node_is_mirror(snodes[i]));
            snodes[i]->hybrid->tau = snodes[i]->tau;
            if (snodes[i]->htau == 0)
              snodes[i]->parent->tau = snodes[i]->tau;
            if (snodes[i]->hybrid->htau == 0)
              snodes[i]->hybrid->parent->tau = snodes[i]->tau;
          }
          else
          {
            assert(node_is_bidirection(snodes[i]));
            assert(!node_is_mirror(snodes[i]));
            assert(snodes[i]->prop_tau && !snodes[i]->hybrid->parent->prop_tau);
            snodes[i]->hybrid->tau        = snodes[i]->tau;
            snodes[i]->right->tau         = snodes[i]->tau;
            snodes[i]->right->hybrid->tau = snodes[i]->tau;
          }
        }
      }
    }
    else
    {
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
        if (opt_scaling)
          gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                gnodeptr[j]->scaler_index);
        gnodeptr[j]->time = gnodeptr[j]->old_time;
      }
      if (opt_migration)
      {
        for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
        {
          if (!gtree[i]->nodes[j]->mi || !gtree[i]->nodes[j]->mi->count)
            continue;

          miginfo_t * mi = gtree[i]->nodes[j]->mi;
          for (k = 0; k < mi->count; ++k)
            mi->me[k].time = mi->me[k].old_time;
        }
      }

      /* revert trans prob matrices */
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        if (gtree[i]->nodes[j]->parent)
          gtree[i]->nodes[j]->pmatrix_index = SWAP_PMAT_INDEX(gtree[i]->edge_count,
                                                              gtree[i]->nodes[j]->pmatrix_index);
      
      if (opt_clock == BPP_CLOCK_CORR && opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
        gtree[i]->lnprior_rates = gtree[i]->old_lnprior_rates;
    }
    if (opt_migration)
    {
      #if 1
      stree_update_mig_subpops(stree, thread_index);
      #else
    /* TODO: Note, the following loop should be a correct and faster alternative
       than calling the stree_update_mig_subpops function. However, those
       divisions cause the migbuffer[j].time to differ numerically with
       the taus in the species tree. We need to include an snode pointer in the
       structure and copy the values instead of dividing */
      for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      {
        snode_t * x = stree->nodes[i];
        for (j = 0; j < x->mb_count; ++j)
          x->migbuffer[j].time /= c;
      }
      #endif
    }
  }
  free(snodes);

  if (!opt_est_theta)
    free(notheta_old_logpr);

  return accepted;

}
