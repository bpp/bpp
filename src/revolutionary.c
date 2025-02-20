/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

void revolutionary_spr_tselect_logl(gnode_t * mnode, gnode_t ** target_list, long target_count, locus_t * locus, double * weights)
{
  /*
                      *              if t1 is the selected target node, then
                     / \             anode subtree is placed on branch t1<->x
                    /   \
            moved  *     \
                  / \     \
                 /   \     * x
        anode   *         / \
               / \       /   \
              /   \     /     \
               red      t1   t2
            lineages    targets
  */
   long i, j, n;
   double logl;
   double tsum;
   gnode_t * tnode;
   const double * mclv, *tclv;  /* CLV of moved and target node */
   double * nptr;
   const double * mmat = locus->pmatrix[mnode->pmatrix_index];       /* transition probability matrix for moved node */
   /* space for new p-matrix for target: */
   double * tmat = pll_aligned_alloc(locus->states * locus->states_padded * locus->rate_cats * sizeof(double), locus->alignment);
   /* allocate storage space for CLV subtree root, which is to be computed: */
  double * clv = pll_aligned_alloc(locus->sites * locus->states_padded * locus->rate_cats * sizeof(double), locus->alignment);
  /* temp space for normalized CLV of moved node:  */
  //double * nmclv = pll_aligned_alloc(locus->sites * locus->states_padded * locus->rate_cats * sizeof(double), locus->alignment);

  /* space for normalized CLV of target node.  Ziheng: nclv is temp space, for scaled tclv for target. */
  double * ntclv = pll_aligned_alloc(locus->sites * locus->states_padded * locus->rate_cats * sizeof(double), locus->alignment);
  double tlength[1] = { mnode->parent->time * 0.314 }; /* if(opt_revolutionary_spr_method == 2) */
  unsigned int matrix_indices[1] = { 0 };
  double * matrices[1] = { tmat };
  assert(locus->states == 4);

  /* go through the list of target nodes (edge from target node to its parent is the target branch) */
  /* Ziheng: get CLV of moved node.  This is scaled outside the loop for target branches.  */
  mclv = locus->clv[mnode->clv_index];
  
  for (n = 0; n < target_count; ++n)
  {
    if (opt_revolutionary_spr_debug > 2)
    {
      printf("Target %ld\n", n);
      printf("pat       Conditional Probabilities\n");
      printf("          A        C        G        T\n\n");
    }

    /* get target node */
    tnode = target_list[n];
    if (opt_revolutionary_spr_method == 1)
       tlength[0] = mnode->parent->time - tnode->time;

    /* get CLV of target node */
    tclv = locus->clv[tnode->clv_index];

    /* now construct the normalized CLV */
    /*** Ziheng 2018.04.15 ***
    This assumes JC with one rate class.  Think later about what to do if we use GTR+G4. 
    ***/
    nptr = ntclv;
    for (i = 0; i < locus->sites; ++i)
    {
      for (tsum=0, j = 0; j < locus->states; ++j)
        tsum += tclv[j];
      for (j = 0; j < locus->states; ++j)
        nptr[j] = tclv[j]/tsum;

      tclv += locus->states;
      nptr += locus->states;
    }

    /* get transition probability matrix for moved node */
    mmat = locus->pmatrix[mnode->pmatrix_index];  /* Ziheng: mmat is already calculated outside this routine */
    /* construct transition probability matrix tmat for target node */
    pll_core_update_pmatrix_4x4_jc69(matrices, locus->states, locus->rate_cats, NULL, tlength, matrix_indices, locus->param_indices, 1, locus->attributes);
                                     
    /* TODO: Account for scalers - currently disabled */
    /* update conditional probabilities vector clv using vectors mclv and ntclv and matrices mmat and (new) tmat */
    pll_core_update_partial_ii(locus->states, locus->sites, locus->rate_cats, clv, NULL, mclv, ntclv, mmat, tmat, NULL, NULL, locus->attributes);

    /* compute log-likelihood of tree having root with conditional probabilities vector clv */
    logl = pll_core_root_loglikelihood(locus->states, locus->sites, locus->rate_cats, clv, NULL, locus->frequencies, locus->rate_weights, locus->pattern_weights, locus->param_indices, NULL, locus->attributes);
    weights[n] = logl;

    /* if debugging information enabled */
    if (opt_revolutionary_spr_debug > 2)
    {
      long j;
      for (j = 0; j < locus->sites; ++j)
      {
        long k;
        printf("      m: ");
        for (k = 0; k < locus->states; ++k)
          printf(" %lf", mclv[j*4+k]);
        printf("\n");
        printf("%3d   t: ", locus->pattern_weights[j]);
        for (k = 0; k < locus->states; ++k)
          printf(" %lf", ntclv[j*4+k]);
        printf("\n");
        printf("     y*: ");
        for (k=0; k < locus->states; ++k)
          printf(" %lf", clv[j*4+k]);
        printf("\n\n");
      }
      printf("  logl: %12.6e\n\n", weights[n]);
    }
  }

  /* deallocate */
  pll_aligned_free(clv);
  pll_aligned_free(ntclv);
  pll_aligned_free(tmat);
}

void rev_spr_tselect(gnode_t * mnode,
                     double t,
                     gnode_t ** targets,
                     unsigned int target_count,
                     locus_t * locus,
                     double * weights)
{
  unsigned int i,j,k;
  size_t matsize;
  size_t clvsize;
  double * mmat;
  double * tmat;
  unsigned int matrix_indices[2] = {0,1};
  double length[1];
  double * matrices[1];

  assert(locus->states == 4);

  matsize = locus->states*locus->states_padded*locus->rate_cats*sizeof(double);
  clvsize = locus->sites*locus->states_padded*locus->rate_cats*sizeof(double);

  /* allocate space for normalized CLV vectors and for subtree root */
  double * ntclv = pll_aligned_alloc(clvsize,locus->alignment);
  double * clv = pll_aligned_alloc(clvsize,locus->alignment);

  /* allocate memory for the two p-matrices (moved node and target) */
  mmat = pll_aligned_alloc(matsize,locus->alignment);
  tmat = pll_aligned_alloc(matsize,locus->alignment);

  /* TODO: we assume scaling for numerical underflow is disabled */
  assert(opt_scaling == 0);

  /* compute p-matrix for moved node for new branch length */
  matrices[0] = mmat;  length[0] = t - mnode->time;
  pll_core_update_pmatrix_4x4_jc69(matrices,
                                   locus->states,
                                   locus->rate_cats,
                                   NULL,
                                   length,
                                   matrix_indices,
                                   locus->param_indices,
                                   1,
                                   locus->attributes);

  for (i = 0; i < target_count; ++i)
  {
    gnode_t * tnode = targets[i];

    /* compute p-matrices for target node from new branch length */
    matrices[0] = tmat;  length[0] = t - tnode->time;
    pll_core_update_pmatrix_4x4_jc69(matrices,
                                     locus->states,
                                     locus->rate_cats,
                                     NULL,
                                     length,
                                     matrix_indices,
                                     locus->param_indices,
                                     1,
                                     locus->attributes);

    /* get CLV vectors for moved (mclv) and target (tclv) node */
    const double * mclv = locus->clv[mnode->clv_index];
    const double * tclv = locus->clv[tnode->clv_index];

    /* construct normalized CLVs */
    double * tptr = ntclv;
    for (j = 0; j < locus->sites; ++j)
    {
      double tsum = 0;
      for (k=0; k < locus->states; ++k)
        tsum += tclv[k];
      for (k = 0; k < locus->states; ++k)
        tptr[k] = tclv[k]/tsum;

      tclv += locus->states;
      tptr += locus->states;
    }

    /* update conditional probabilities vector clv using vectors mclv and ntclv
     * and matrices mmat and tmat. Numerical scaling is assumed disabled */
    pll_core_update_partial_ii(locus->states,
                               locus->sites,
                               locus->rate_cats,
                               clv,
                               NULL,
                               mclv,
                               ntclv,
                               mmat,
                               tmat,
                               NULL,
                               NULL,
                               locus->attributes);
    
    /* compute log-likelihood of tree with clv as root node clv */
    weights[i] = pll_core_root_loglikelihood(locus->states,
                                             locus->sites,
                                             locus->rate_cats,
                                             clv,
                                             NULL,
                                             locus->frequencies,
                                             locus->rate_weights,
                                             locus->pattern_weights,
                                             locus->param_indices,
                                             NULL,
                                             locus->attributes);
  }

  /* deallocate */
  pll_aligned_free(clv);
  pll_aligned_free(ntclv);
  pll_aligned_free(mmat);
  pll_aligned_free(tmat);
}
