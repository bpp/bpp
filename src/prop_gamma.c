/*
    Copyright (C) 2016-2019 Tomas Flouri, Bruce Rannala and Ziheng Yang

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
#define SWAP_PMAT_INDEX(e,i) (i) = (((e)+(i))%((e)<<1))
#define SWAP_SCALER_INDEX(n,i) (((n)+((i)-1))%(2*(n)-2))

/* TODO: The below two functions are duplicated several times. Move them as
   public functions in gtree.c */
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

static long propose_alpha(stree_t * stree,
                          locus_t * locus,
                          gtree_t * gtree,
                          long msa_index,
                          long thread_index)
{
  unsigned int m,n;
  long accepted = 0;
  double lnacceptance;
  double logl;
  double alpha_old,alpha_new;
  double loga_old, loga_new;
  double * old_rates;
  gnode_t ** gt_nodes;

  double minv = -99;
  double maxv =  99;

  /* allocate temporary space for gene tree traversal */
  gt_nodes = (gnode_t **)xmalloc((gtree->tip_count+gtree->inner_count) *
                                 sizeof(gnode_t *));

  alpha_old = locus->rates_alpha;
  loga_old  = log(alpha_old);

  loga_new  = loga_old + opt_finetune_alpha * legacy_rnd_symmetrical(thread_index);
  loga_new  = reflect(loga_new, minv, maxv, thread_index);

  alpha_new = exp(loga_new);

  lnacceptance = loga_new - loga_old;

  locus->rates_alpha = alpha_new;

  /* store old rates */
  old_rates = (double *)xmalloc((size_t)(locus->rate_cats) * sizeof(double));
  memcpy(old_rates,locus->rates,(size_t)(locus->rate_cats) * sizeof(double));

  /* update locus->rates with new rates */
  pll_compute_gamma_cats(locus->rates_alpha,
                         locus->rate_cats,
                         locus->rates,
                         PLL_GAMMA_RATES_MEAN);

  /* swap pmatrix indices to new buffers, and update pmatrices. No need to
     recompute eigen decomposition */
  for (m=0,n=0; m < gtree->tip_count+gtree->inner_count; ++m)
  {
    gnode_t * p = gtree->nodes[m];
    if (p->parent)
    {
      SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
      gt_nodes[n++] = p;
    }
  }
  locus_update_matrices(locus,gt_nodes,stree,msa_index,n);


  /* get postorder traversal of inner nodes, swap CLV indidces to point to new
     buffer, and update partials */
  gtree_all_partials(gtree->root,gt_nodes,&n);
  for (m = 0; m < n; ++m)
  {
    gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                            gt_nodes[m]->clv_index);
    if (opt_scaling)
      gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                    gt_nodes[m]->scaler_index);
  }
  locus_update_partials(locus,gt_nodes,n);

  /* compute log-likelihood */
  logl = locus_root_loglikelihood(locus,gtree->root,locus->param_indices,NULL);

  lnacceptance += (logl - gtree->logl);

  /* prior rate */
  //lnacceptance += (opt_alpha_alpha-1) * (loga_new - loga_old) - 
  lnacceptance += (opt_alpha_alpha-1) * log(alpha_new/alpha_old) - 
                  (opt_alpha_beta) * (alpha_new - alpha_old);

  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* accepted */
    accepted = 1;
    gtree->logl = logl;
  }
  else
  {
    /* rejected */

    locus->rates_alpha = alpha_old;

    /* revert pmatrices */
    for (m = 0; m < gtree->tip_count + gtree->inner_count; ++m)
    {
      gnode_t * p = gtree->nodes[m];
      if (p->parent)
        SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
    }

    /* revert CLV */
    for (m = 0; m < n; ++m)
    {
      gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                              gt_nodes[m]->clv_index);
      if (opt_scaling)
        gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                      gt_nodes[m]->scaler_index);
    }

    /* revert old rates */
    pll_set_category_rates(locus,old_rates);
  }
  free(gt_nodes);
  free(old_rates);

  return accepted;
}

double locus_propose_alpha_serial(stree_t * stree, locus_t ** locus, gtree_t ** gtree)
{
  long i;
  long accepted = 0;
  long thread_index = 0;
  long candidates = 0;

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (locus[i]->dtype == BPP_DATA_DNA && locus[i]->rate_cats > 1)
    {
      ++candidates;
      accepted += propose_alpha(stree,locus[i],gtree[i],i,thread_index);
    }
  }

  if (!accepted)
    return 0;

  return ((double)accepted/candidates);
}

void locus_propose_alpha_parallel(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree,
                                  long locus_start,
                                  long locus_count,
                                  long thread_index,
                                  long * p_proposal_count,
                                  long * p_accepted)
{
  unsigned int i;
  long accepted = 0;
  long candidates = 0;

  assert(locus_start >= 0);
  assert(locus_count > 0);

  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    if (locus[i]->dtype == BPP_DATA_DNA && locus[i]->rate_cats > 1)
    {
      ++candidates;
      accepted += propose_alpha(stree,locus[i],gtree[i],i,thread_index);
    }
  }

  *p_proposal_count = candidates;
  *p_accepted = accepted;
}
