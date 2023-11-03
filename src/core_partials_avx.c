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

#if defined(__x86_64__) || defined(_M_AMD64)

static void fill_parent_scaler(unsigned int scaler_size,
                               unsigned int * parent_scaler,
                               const unsigned int * left_scaler,
                               const unsigned int * right_scaler)
{
  unsigned int i;

  if (!left_scaler && !right_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);
  else if (left_scaler && right_scaler)
  {
    memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    for (i = 0; i < scaler_size; ++i)
      parent_scaler[i] += right_scaler[i];
  }
  else
  {
    if (left_scaler)
      memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    else
      memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * scaler_size);
  }
}

void pll_core_create_lookup_avx(unsigned int states,
                                unsigned int rate_cats,
                                double * ttlookup,
                                const double * left_matrix,
                                const double * right_matrix,
                                const unsigned int * tipmap,
                                unsigned int tipmap_size)
{
  if (states == 4)
  {
    pll_core_create_lookup_4x4_avx(rate_cats,
                                   ttlookup,
                                   left_matrix,
                                   right_matrix);
    return;
  }

  if (states == 20)
  {
    pll_core_create_lookup_20x20_avx(rate_cats,
                                   ttlookup,
                                   left_matrix,
                                   right_matrix,
                                   tipmap,
                                   tipmap_size);
    return;
  }

  unsigned int i,j,k,n,m;
  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int maxstates = tipmap_size;
  unsigned int index = 0;

  unsigned int log2_maxstates = (unsigned int)ceil(log2(maxstates));
  unsigned int span_padded = states_padded*rate_cats;

  /* precompute first the entries that contain only one 1 */
  double termj = 0;
  double termk = 0;

  const double * jmat;
  const double * kmat;
  double * lookup;

  /* go through all pairs j,k of states for the two tips; i is the inner
     node state */
  for (j = 0; j < maxstates; ++j)
  {
    for (k = 0; k < maxstates; ++k)
    {
      jmat = left_matrix;
      kmat = right_matrix;

      /* find offset of state-pair in the precomputation table */
      lookup = ttlookup;
      lookup += ((j << log2_maxstates) + k)*span_padded;

      /* precompute the likelihood for each state and each rate */
      for (n = 0; n < rate_cats; ++n)
      {
        index = 0;
        for (i = 0; i < states; ++i)
        {
          termj = 0;
          termk = 0;

          unsigned int jstate = tipmap[j];
          unsigned int kstate = tipmap[k];

          /* decompose basecall into the encoded residues and set the appropriate
             positions in the tip vector */
          for (m = 0; m < states; ++m)
          {
            if (jstate & 1)
              termj += jmat[m];

            if (kstate & 1)
              termk += kmat[m];

            jstate >>= 1;
            kstate >>= 1;
          }

          jmat += states_padded;
          kmat += states_padded;
          lookup[index++] = termj*termk;
        }
        /* this is to avoid valgrind warnings on accessing uninitialized memory
           when using AVX and states are not a multiple of 4 */
        if (states_padded-states)
          memset(lookup+index, 0, (states_padded-states)*sizeof(double));

        lookup += states_padded;
      }
    }
  }
}

void pll_core_create_lookup_20x20_avx(unsigned int rate_cats,
                                      double * ttlookup,
                                      const double * left_matrix,
                                      const double * right_matrix,
                                      const unsigned int * tipmap,
                                      unsigned int tipmap_size)
{
  unsigned int i,j,k,n,m;
  unsigned int states = 20;
  unsigned int states_padded = 20;
  unsigned int maxstates = tipmap_size;

  unsigned int log2_maxstates = (unsigned int)ceil(log2(maxstates));
  unsigned int span_padded = states_padded*rate_cats;

  /* precompute first the entries that contain only one 1 */
  double terml = 0;
  double termr = 0;

  double * lookupl = pll_aligned_alloc(span_padded*maxstates*sizeof(double),
                                       PLL_ALIGNMENT_AVX);

  double * lookupr = pll_aligned_alloc(span_padded*maxstates*sizeof(double),
                                       PLL_ALIGNMENT_AVX);

  const double * lmat;
  const double * rmat;
  double * lookup;
  double * ll = lookupl;
  double * lr = lookupr;

  for (j = 0; j < maxstates; ++j)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    unsigned int state = tipmap[j];

    int ss = PLL_POPCOUNT(state) == 1 ? PLL_CTZ(state) : -1;

    for (n = 0; n < rate_cats; ++n)
    {
      for (i = 0; i < states; ++i)
      {
        if (ss != -1)
        {
          /* special case for non-ambiguous states */
          terml = lmat[ss];
          termr = rmat[ss];
        }
        else
        {
          terml = 0;
          termr = 0;
          for (m = 0; m < states; ++m)
          {
            if ((state>>m) & 1)
            {
              terml += lmat[m];
              termr += rmat[m];
            }
          }
        }

        lmat += states_padded;
        rmat += states_padded;

        ll[i] = terml;
        lr[i] = termr;
      }

      ll += states_padded;
      lr += states_padded;
    }
  }

  assert(ll == lookupl + maxstates * span_padded);
  assert(lr == lookupr + maxstates * span_padded);

  /* go through all pairs j,k of states for the two tips; i is the inner
     node state */
  for (j = 0; j < maxstates; ++j)
  {
    for (k = 0; k < maxstates; ++k)
    {
      /* find offset of state-pair in the precomputation table */
      lookup = ttlookup;
      lookup += ((j << log2_maxstates) + k)*span_padded;

      ll = lookupl + j * span_padded;
      lr = lookupr + k * span_padded;

      /* precompute the likelihood for each state and each rate */
      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states; i += 4)
        {
          __m256d v_terml = _mm256_load_pd(ll+i);
          __m256d v_termr = _mm256_load_pd(lr+i);

          __m256d v_prod = _mm256_mul_pd(v_terml,v_termr);

          _mm256_store_pd(lookup+i, v_prod);
        }

        lookup += states_padded;
        ll += states_padded;
        lr += states_padded;
      }
    }
  }

  pll_aligned_free(lookupl);
  pll_aligned_free(lookupr);
}

void pll_core_create_lookup_4x4_avx(unsigned int rate_cats,
                                    double * lookup,
                                    const double * left_matrix,
                                    const double * right_matrix)
{
  unsigned int j,k,n;
  unsigned int maxstates = 16;
  unsigned int states = 4;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
  __m256i jmask,kmask;

  const double * jmat;
  const double * kmat;

  lookup += maxstates*states*rate_cats;

  for (j = 1; j < maxstates; ++j)
  {
    jmask = _mm256_set_epi64x(
               ((j >> 3) & 1) ? ~0 : 0,
               ((j >> 2) & 1) ? ~0 : 0,
               ((j >> 1) & 1) ? ~0 : 0,
               (j & 1) ? ~0 : 0);
    
    lookup += states*rate_cats;

    for (k = 1; k < maxstates; ++k)
    {
      jmat = left_matrix;
      kmat = right_matrix;

      kmask = _mm256_set_epi64x(
                 ((k >> 3) & 1) ? ~0 : 0,
                 ((k >> 2) & 1) ? ~0 : 0,
                 ((k >> 1) & 1) ? ~0 : 0,
                 (k & 1) ? ~0 : 0);

      for (n = 0; n < rate_cats; ++n)
      {
        xmm0 = _mm256_maskload_pd(jmat,jmask);
        ymm0 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        xmm1 = _mm256_maskload_pd(jmat,jmask);
        ymm1 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        xmm2 = _mm256_maskload_pd(jmat,jmask);
        ymm2 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        xmm3 = _mm256_maskload_pd(jmat,jmask);
        ymm3 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        /* compute x */
        xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
        xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

        xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
        xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

        xmm0 = _mm256_add_pd(xmm4,xmm5);
        xmm1 = _mm256_add_pd(xmm6,xmm7);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
        xmm4 = _mm256_add_pd(xmm2,xmm3);

        /* compute y */
        ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
        ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

        ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
        ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

        ymm0 = _mm256_add_pd(ymm4,ymm5);
        ymm1 = _mm256_add_pd(ymm6,ymm7);

        ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
        ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
        ymm4 = _mm256_add_pd(ymm2,ymm3);

        /* compute x*y */
        xmm0 = _mm256_mul_pd(xmm4,ymm4);

        _mm256_store_pd(lookup, xmm0);

        lookup += states;
      }
    }
  }
}

void pll_core_update_partial_ii_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const double * left_clv,
                                        const double * right_clv,
                                        const double * left_matrix,
                                        const double * right_matrix,
                                        const unsigned int * left_scaler,
                                        const unsigned int * right_scaler,
                                        unsigned int attrib)
{
  unsigned int states = 4;
  unsigned int n,k,i;

  const double * lmat;
  const double * rmat;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

  unsigned int span = states * rate_cats;


  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }

  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      /* compute vector of x */
      xmm4 = _mm256_load_pd(lmat);
      xmm5 = _mm256_load_pd(left_clv);
      xmm0 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm5 = _mm256_load_pd(right_clv);
      ymm0 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm4 = _mm256_load_pd(lmat);
      xmm1 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm1 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm4 = _mm256_load_pd(lmat);
      xmm2 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm2 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm4 = _mm256_load_pd(lmat);
      xmm3 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm3 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      /* compute x */
      xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
      xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

      xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
      xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

      xmm0 = _mm256_add_pd(xmm4,xmm5);
      xmm1 = _mm256_add_pd(xmm6,xmm7);

      xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
      xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
      xmm4 = _mm256_add_pd(xmm2,xmm3);

      /* compute y */
      ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
      ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

      ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
      ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

      ymm0 = _mm256_add_pd(ymm4,ymm5);
      ymm1 = _mm256_add_pd(ymm6,ymm7);

      ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
      ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
      ymm4 = _mm256_add_pd(ymm2,ymm3);

      /* compute x*y */
      xmm0 = _mm256_mul_pd(xmm4,ymm4);

      /* check if scaling is needed for the current rate category */
      __m256d v_cmp = _mm256_cmp_pd(xmm0, v_scale_threshold, _CMP_LT_OS);
      const unsigned int rate_mask = _mm256_movemask_pd(v_cmp);

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          xmm0 = _mm256_mul_pd(xmm0,v_scale_factor);
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      _mm256_store_pd(parent_clv, xmm0);

      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span;
      for (i = 0; i < span; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

void pll_core_update_partial_tt_avx(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const unsigned char * right_tipchars,
                                    const double * lookup,
                                    unsigned int tipstates_count,
                                    unsigned int attrib)
{
  unsigned int j,k,n;
  unsigned int log2_maxstates = (unsigned int)ceil(log2(tipstates_count));
  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;
  const double * offset;

  if (states == 4)
  {
    pll_core_update_partial_tt_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_tipchars,
                                       lookup,
                                       attrib);
    return;
  }

  size_t scaler_size = (attrib & PLL_ATTRIB_RATE_SCALERS) ?
                                                        sites*rate_cats : sites;

  if (parent_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);

  for (n = 0; n < sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    offset = lookup;
    offset += ((j << log2_maxstates) + k)*span_padded;

    memcpy(parent_clv, offset, span_padded*sizeof(double));

    parent_clv += span_padded;
  }
}

void pll_core_update_partial_tt_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const unsigned char * left_tipchars,
                                        const unsigned char * right_tipchars,
                                        const double * lookup,
                                        unsigned int attrib)
{
  unsigned int j,k,n;
  unsigned int states = 4;
  unsigned int span = states*rate_cats;
  const double * offset;

  size_t scaler_size = (attrib & PLL_ATTRIB_RATE_SCALERS) ?
                                                        sites*rate_cats : sites;

  if (parent_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);

  for (n = 0; n < sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    offset = lookup;
    offset += ((j << 4) + k)*span;

    //    memcpy(parent_clv, offset, span*sizeof(double));
    for (k = 0; k < span; k += 4)
    {
      __m256d v_clv = _mm256_load_pd(offset + k);
      _mm256_store_pd(parent_clv + k, v_clv);
    }

    parent_clv += span;
  }
}

void pll_core_update_partial_ti_avx(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * right_scaler,
                                    const unsigned int * tipmap,
                                    unsigned int tipmap_size,
                                    unsigned int attrib)
{
  unsigned int i,j,k,n;

  const double * lmat;
  const double * rmat;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;

  unsigned int lstate;

  /* dedicated functions for 4x4 matrices (DNA) */
  if (states == 4)
  {
    pll_core_update_partial_ti_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       right_scaler,
                                       attrib);
    return;
  }

  /* dedicated functions for 20x20 matrices (AA) */
  if (states == 20)
  {
    pll_core_update_partial_ti_20x20_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       right_scaler,
                                       tipmap,
                                       tipmap_size,
                                       attrib);
    return;
  }

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  __m256i mask;

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    scale_mask = init_mask;

    lstate = tipmap[left_tipchars[n]];

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_rclv;

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* set position of least significant bit in character state */
        register int lsb = 0;

        /* iterate over quadruples of columns */
        for (j = 0; j < states_padded; j += 4)
        {

          /* set mask */
          mask = _mm256_set_epi64x(
                    ((lstate >> (lsb+3)) & 1) ? ~0 : 0,
                    ((lstate >> (lsb+2)) & 1) ? ~0 : 0,
                    ((lstate >> (lsb+1)) & 1) ? ~0 : 0,
                    ((lstate >> (lsb+0)) & 1) ? ~0 : 0);

          if ((lstate >> lsb) & 0b1111)
          {
            v_mat    = _mm256_maskload_pd(lm0,mask);
            v_terma0 = _mm256_add_pd(v_terma0,v_mat);

            v_mat    = _mm256_maskload_pd(lm1,mask);
            v_terma1 = _mm256_add_pd(v_terma1,v_mat);

            v_mat    = _mm256_maskload_pd(lm2,mask);
            v_terma2 = _mm256_add_pd(v_terma2,v_mat);

            v_mat    = _mm256_maskload_pd(lm3,mask);
            v_terma3 = _mm256_add_pd(v_terma3,v_mat);
          }

          lsb += 4;

          lm0 += 4;
          lm1 += 4;
          lm2 += 4;
          lm3 += 4;

          v_rclv    = _mm256_load_pd(right_clv+j);

          /* row 0 */
          v_mat    = _mm256_load_pd(rm0);
          v_termb0 = _mm256_add_pd(v_termb0,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm0 += 4;

          /* row 1 */
          v_mat    = _mm256_load_pd(rm1);
          v_termb1 = _mm256_add_pd(v_termb1,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm1 += 4;

          /* row 2 */
          v_mat    = _mm256_load_pd(rm2);
          v_termb2 = _mm256_add_pd(v_termb2,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm2 += 4;

          /* row 3 */
          v_mat    = _mm256_load_pd(rm3);
          v_termb3 = _mm256_add_pd(v_termb3,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm3 += 4;
        }

        /* point pmatrix to the next four rows */ 
        lmat = lm3;
        rmat = rm3;

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);

      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      right_clv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

void pll_core_update_partial_ti_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const unsigned char * left_tipchar,
                                        const double * right_clv,
                                        const double * left_matrix,
                                        const double * right_matrix,
                                        const unsigned int * right_scaler,
                                        unsigned int attrib)
{
  unsigned int states = 4;
  unsigned int i,k,n;

  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;

  const double * lmat;
  const double * rmat;

  unsigned int span = states * rate_cats;
  unsigned int lstate;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
  __m256i mask;

  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(64*rate_cats*sizeof(double),
                                      PLL_ALIGNMENT_AVX);
  if (!lookup)
    fatal("Cannot allocate space for precomputation.");

  /* skip first entry of lookup table as it is never used */
  double * ptr = lookup + 4*rate_cats;

  for (i = 1; i < 16; ++i)
  {
    lmat = left_matrix;

    /* mask the entries of pmatrix row to be loaded */
    mask = _mm256_set_epi64x(
              ((i >> 3) & 1) ? ~0 : 0,
              ((i >> 2) & 1) ? ~0 : 0,
              ((i >> 1) & 1) ? ~0 : 0,
              ((i >> 0) & 1) ? ~0 : 0);

    for (k = 0; k < rate_cats; ++k)
    {
      /* masked load of the four matrix rows */
      xmm0 = _mm256_maskload_pd(lmat+0, mask);
      xmm1 = _mm256_maskload_pd(lmat+4, mask);
      xmm2 = _mm256_maskload_pd(lmat+8, mask);
      xmm3 = _mm256_maskload_pd(lmat+12,mask);

      /* create a vector containing the sums of xmm0, xmm1, xmm2, xmm3 */
      xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
      xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

      xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
      xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

      xmm0 = _mm256_add_pd(xmm4,xmm5);
      xmm1 = _mm256_add_pd(xmm6,xmm7);

      xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
      xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
      xmm4 = _mm256_add_pd(xmm2,xmm3);

      /* store the result (four sums) to the lookup table */
      _mm256_store_pd(ptr,xmm4);

      /* move pointers */
      ptr  += 4;
      lmat += 16;
    }
  }

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;

    /* update the parent scaler with the scaler of the right child */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }

  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  /* iterate over sites and compute CLV entries */
  for (n = 0; n < sites; ++n)
  {
    rmat = right_matrix;

    scale_mask = init_mask;

    lstate = left_tipchar[n];

    unsigned int loffset = rate_cats*lstate*4;

    for (k = 0; k < rate_cats; ++k)
    {
      ymm4 = _mm256_load_pd(rmat);
      ymm5 = _mm256_load_pd(right_clv);
      ymm0 = _mm256_mul_pd(ymm4,ymm5);
      rmat += states;

      ymm4 = _mm256_load_pd(rmat);
      ymm1 = _mm256_mul_pd(ymm4,ymm5);
      rmat += states;

      ymm4 = _mm256_load_pd(rmat);
      ymm2 = _mm256_mul_pd(ymm4,ymm5);
      rmat += states;

      ymm4 = _mm256_load_pd(rmat);
      ymm3 = _mm256_mul_pd(ymm4,ymm5);
      rmat += states;

      /* load x from precomputed lookup table into xmm4 */
      xmm4 = _mm256_load_pd(lookup+loffset);

      /* compute y */
      ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
      ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

      ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
      ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

      ymm0 = _mm256_add_pd(ymm4,ymm5);
      ymm1 = _mm256_add_pd(ymm6,ymm7);

      ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
      ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
      ymm4 = _mm256_add_pd(ymm2,ymm3);

      /* compute x*y */
      xmm0 = _mm256_mul_pd(xmm4,ymm4);

      /* check if scaling is needed for the current rate category */
      __m256d v_cmp = _mm256_cmp_pd(xmm0, v_scale_threshold, _CMP_LT_OS);
      const unsigned int rate_mask = _mm256_movemask_pd(v_cmp);

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          xmm0 = _mm256_mul_pd(xmm0, v_scale_factor);
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      _mm256_store_pd(parent_clv, xmm0);

      parent_clv += states;
      right_clv  += states;
      loffset    += 4;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span;
      for (i = 0; i < span; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
  pll_aligned_free(lookup);
}

void pll_core_update_partial_ti_20x20_avx(unsigned int sites,
                                          unsigned int rate_cats,
                                          double * parent_clv,
                                          unsigned int * parent_scaler,
                                          const unsigned char * left_tipchar,
                                          const double * right_clv,
                                          const double * left_matrix,
                                          const double * right_matrix,
                                          const unsigned int * right_scaler,
                                          const unsigned int * tipmap,
                                          unsigned int tipmap_size,
                                          unsigned int attrib)
{
  unsigned int states = 20;
  unsigned int states_padded = states;
  unsigned int maxstates = tipmap_size;
  unsigned int i,j,k,n,m;

  const double * lmat;
  const double * rmat;

  unsigned int span_padded = states_padded * rate_cats;
  unsigned int lstate;

  __m256d xmm0,xmm1,xmm2,xmm3;

  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(maxstates*span_padded*sizeof(double),
                                      PLL_ALIGNMENT_AVX);
  if (!lookup)
    fatal("Cannot allocate space for precomputation.");

  double * ptr = lookup;

  /* precompute left-side values and store them in lookup table */
  for (j = 0; j < maxstates; ++j)
  {
    lmat = left_matrix;

    unsigned int state = tipmap[j];

    int ss = PLL_POPCOUNT(state) == 1 ? PLL_CTZ(state) : -1;

    for (n = 0; n < rate_cats; ++n)
    {
      for (i = 0; i < states; ++i)
      {
        double terml;
        if (ss != -1)
        {
          /* special case for non-ambiguous states */
          terml = lmat[ss];
        }
        else
        {
          terml = 0;
          for (m = 0; m < states; ++m)
          {
            if ((state>>m) & 1)
            {
              terml += lmat[m];
            }
          }
        }

        lmat += states;

        ptr[i] = terml;
      }

      ptr += states;
    }
  }

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* iterate over sites and compute CLV entries */
  for (n = 0; n < sites; ++n)
  {
    rmat = right_matrix;

    scale_mask = init_mask;

    lstate = (unsigned int) left_tipchar[n];

    unsigned int loffset = lstate*span_padded;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_rclv;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
        for (j = 0; j < states_padded; j += 4)
        {
          v_rclv    = _mm256_load_pd(right_clv+j);

          /* row 0 */
          v_mat    = _mm256_load_pd(rm0);
          v_termb0 = _mm256_add_pd(v_termb0,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm0 += 4;

          /* row 1 */
          v_mat    = _mm256_load_pd(rm1);
          v_termb1 = _mm256_add_pd(v_termb1,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm1 += 4;

          /* row 2 */
          v_mat    = _mm256_load_pd(rm2);
          v_termb2 = _mm256_add_pd(v_termb2,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm2 += 4;

          /* row 3 */
          v_mat    = _mm256_load_pd(rm3);
          v_termb3 = _mm256_add_pd(v_termb3,
                                   _mm256_mul_pd(v_mat,v_rclv));
          rm3 += 4;
        }

        /* point pmatrix to the next four rows */
        rmat = rm3;

        /* load x from precomputed lookup table */
        __m256d v_terma_sum = _mm256_load_pd(lookup+loffset);
        loffset += 4;

        /* compute termb */
        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      rmat -= displacement;

      parent_clv += states_padded;
      right_clv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
  pll_aligned_free(lookup);
}

void pll_core_update_partial_ii_avx(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const double * left_clv,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * left_scaler,
                                    const unsigned int * right_scaler,
                                    unsigned int attrib)
{
  unsigned int i,j,k,n;

  const double * lmat;
  const double * rmat;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;

  /* dedicated functions for 4x4 matrices */
  if (states == 4)
  {
    pll_core_update_partial_ii_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_clv,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       left_scaler,
                                       right_scaler,
                                       attrib);
    return;
  }

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_lclv;
        __m256d v_rclv;

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
        for (j = 0; j < states_padded; j += 4)
        {
          v_lclv    = _mm256_load_pd(left_clv+j);
          v_rclv    = _mm256_load_pd(right_clv+j);

          /* row 0 */
          v_mat    = _mm256_load_pd(lm0);
          v_terma0 = _mm256_add_pd(v_terma0,
                                   _mm256_mul_pd(v_mat,v_lclv));
          v_mat    = _mm256_load_pd(rm0);
          v_termb0 = _mm256_add_pd(v_termb0,
                                   _mm256_mul_pd(v_mat,v_rclv));
          lm0 += 4;
          rm0 += 4;

          /* row 1 */
          v_mat    = _mm256_load_pd(lm1);
          v_terma1 = _mm256_add_pd(v_terma1,
                                   _mm256_mul_pd(v_mat,v_lclv));
          v_mat    = _mm256_load_pd(rm1);
          v_termb1 = _mm256_add_pd(v_termb1,
                                   _mm256_mul_pd(v_mat,v_rclv));
          lm1 += 4;
          rm1 += 4;

          /* row 2 */
          v_mat    = _mm256_load_pd(lm2);
          v_terma2 = _mm256_add_pd(v_terma2,
                                   _mm256_mul_pd(v_mat,v_lclv));
          v_mat    = _mm256_load_pd(rm2);
          v_termb2 = _mm256_add_pd(v_termb2,
                                   _mm256_mul_pd(v_mat,v_rclv));
          lm2 += 4;
          rm2 += 4;

          /* row 3 */
          v_mat    = _mm256_load_pd(lm3);
          v_terma3 = _mm256_add_pd(v_terma3,
                                   _mm256_mul_pd(v_mat,v_lclv));
          v_mat    = _mm256_load_pd(rm3);
          v_termb3 = _mm256_add_pd(v_termb3,
                                   _mm256_mul_pd(v_mat,v_rclv));
          lm3 += 4;
          rm3 += 4;
        }

        /* point pmatrix to the next four rows */ 
        lmat = lm3;
        rmat = rm3;

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      left_clv   += states_padded;
      right_clv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

#endif
