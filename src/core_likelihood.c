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

double pll_core_root_loglikelihood(unsigned int states,
                                   unsigned int sites,
                                   unsigned int rate_cats,
                                   const double * clv,
                                   const unsigned int * scaler,
                                   double * const * frequencies,
                                   const double * rate_weights,
                                   const unsigned int * pattern_weights,
                                   const unsigned int * freqs_indices,
                                   double * persite_lnl,
                                   unsigned int attrib)
{
  unsigned int i,j,k;
  double logl = 0;
  const double * freqs = NULL;

  double term, term_r;
  double site_lk;

  unsigned int states_padded = states;

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    if (states == 4)
    {
      return pll_core_root_loglikelihood_4x4_sse(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 freqs_indices,
                                                 persite_lnl);
    }
    else
    {
      return pll_core_root_loglikelihood_sse(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             freqs_indices,
                                             persite_lnl);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+1) & 0xFFFFFFFE;
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    if (states == 4)
    {
      return pll_core_root_loglikelihood_4x4_avx(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 freqs_indices,
                                                 persite_lnl);
    }
    else
    {
      return pll_core_root_loglikelihood_avx(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             freqs_indices,
                                             persite_lnl);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
  {
    if (states == 4)
    {
      return pll_core_root_loglikelihood_4x4_avx(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 freqs_indices,
                                                 persite_lnl);
    }
    else
    {
      return pll_core_root_loglikelihood_avx2(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             freqs_indices,
                                             persite_lnl);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif

  /* iterate through sites */
  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      term_r = 0;
      for (k = 0; k < states; ++k)
      {
        term_r += clv[k] * freqs[k];
      }

      term += term_r * rate_weights[j];

      clv += states_padded;
    }

    site_lk = term;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(site_lk);
    if (scaler && scaler[i])
      site_lk += scaler[i] * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = site_lk;

    logl += site_lk;
  }
  return logl;
}

void pll_core_root_likelihood_vector(unsigned int states,
                                     unsigned int sites,
                                     unsigned int rate_cats,
                                     const double * clv,
                                     const unsigned int * scaler,
                                     double * const * frequencies,
                                     const double * rate_weights,
                                     const unsigned int * pattern_weights,
                                     const unsigned int * freqs_indices,
                                     double * persite_lh,
                                     unsigned int attrib)
{
  unsigned int i,j,k;
  //double logl = 0;
  const double * freqs = NULL;

  double term, term_r;
  //double site_lk;

  unsigned int states_padded = states;

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    if (states == 4)
    {
      return pll_core_root_likelihood_vec_4x4_sse(sites,
                                                  rate_cats,
                                                  clv,
                                                  scaler,
                                                  frequencies,
                                                  rate_weights,
                                                  pattern_weights,
                                                  freqs_indices,
                                                  persite_lh);
    }
    else
    {
      return pll_core_root_likelihood_vec_sse(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             freqs_indices,
                                             persite_lh);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+1) & 0xFFFFFFFE;
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    if (states == 4)
    {
      return pll_core_root_likelihood_vec_4x4_avx(sites,
                                                  rate_cats,
                                                  clv,
                                                  scaler,
                                                  frequencies,
                                                  rate_weights,
                                                  pattern_weights,
                                                  freqs_indices,
                                                  persite_lh);
    }
    else
    {
      return pll_core_root_likelihood_vec_avx(states,
                                              sites,
                                              rate_cats,
                                              clv,
                                              scaler,
                                              frequencies,
                                              rate_weights,
                                              pattern_weights,
                                              freqs_indices,
                                              persite_lh);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
  {
    if (states == 4)
    {
      return pll_core_root_likelihood_vec_4x4_avx(sites,
                                                  rate_cats,
                                                  clv,
                                                  scaler,
                                                  frequencies,
                                                  rate_weights,
                                                  pattern_weights,
                                                  freqs_indices,
                                                  persite_lh);
    }
    else
    {
      return pll_core_root_likelihood_vec_avx2(states,
                                               sites,
                                               rate_cats,
                                               clv,
                                               scaler,
                                               frequencies,
                                               rate_weights,
                                               pattern_weights,
                                               freqs_indices,
                                               persite_lh);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif

  /* iterate through sites */
  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      term_r = 0;
      for (k = 0; k < states; ++k)
      {
        term_r += clv[k] * freqs[k];
      }

      term += term_r * rate_weights[j];

      clv += states_padded;
    }

    persite_lh[i] = term;
    #if 0
    site_lk = term;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(site_lk);
    if (scaler && scaler[i])
      site_lk += scaler[i] * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = site_lk;

    logl += site_lk;
    #endif
  }
}

