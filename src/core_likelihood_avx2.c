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

double pll_core_root_loglikelihood_avx2(unsigned int states,
                                        unsigned int sites,
                                        unsigned int rate_cats,
                                        const double * clv,
                                        const unsigned int * scaler,
                                        double * const * frequencies,
                                        const double * rate_weights,
                                        const unsigned int * pattern_weights,
                                        const unsigned int * freqs_indices,
                                        double * persite_lnl)
{
  unsigned int i,j,k;
  double logl = 0;

  const double * freqs = NULL;

  double term, term_r;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  __m256d xmm0, xmm1, xmm3;

  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      xmm3 = _mm256_setzero_pd();

      for (k = 0; k < states_padded; k += 4)
      {
        /* load frequencies for current rate matrix */
        xmm0 = _mm256_load_pd(freqs);

        /* load clv */
        xmm1 = _mm256_load_pd(clv);

        /* multiply with frequencies */
        xmm3 = _mm256_fmadd_pd(xmm0, xmm1, xmm3);

        freqs += 4;
        clv += 4;
      }

      /* add up the elements of xmm2 */
      xmm1 = _mm256_hadd_pd(xmm3,xmm3);

      term_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[2];

      term += term_r * rate_weights[j];
    }

    /* compute site log-likelihood and scale if necessary */
    term = log(term);
    if (scaler && scaler[i])
      term += scaler[i] * log(PLL_SCALE_THRESHOLD);

    term *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = term;

    logl += term;
  }
  return logl;
}

void pll_core_root_likelihood_vec_avx2(unsigned int states,
                                       unsigned int sites,
                                       unsigned int rate_cats,
                                       const double * clv,
                                       const unsigned int * scaler,
                                       double * const * frequencies,
                                       const double * rate_weights,
                                       const unsigned int * pattern_weights,
                                       const unsigned int * freqs_indices,
                                       double * persite_lh)
{
  unsigned int i,j,k;

  const double * freqs = NULL;

  double term, term_r;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  __m256d xmm0, xmm1, xmm3;

  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      xmm3 = _mm256_setzero_pd();

      for (k = 0; k < states_padded; k += 4)
      {
        /* load frequencies for current rate matrix */
        xmm0 = _mm256_load_pd(freqs);

        /* load clv */
        xmm1 = _mm256_load_pd(clv);

        /* multiply with frequencies */
        xmm3 = _mm256_fmadd_pd(xmm0, xmm1, xmm3);

        freqs += 4;
        clv += 4;
      }

      /* add up the elements of xmm2 */
      xmm1 = _mm256_hadd_pd(xmm3,xmm3);

      term_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[2];

      term += term_r * rate_weights[j];
    }

    persite_lh[i] = term;
    #if 0
    /* compute site log-likelihood and scale if necessary */
    term = log(term);
    if (scaler && scaler[i])
      term += scaler[i] * log(PLL_SCALE_THRESHOLD);

    term *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = term;

    logl += term;
    #endif
  }
}
