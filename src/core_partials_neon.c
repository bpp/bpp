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

#ifdef __aarch64__

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

void pll_core_create_lookup_neon(unsigned int states,
                                 unsigned int rate_cats,
                                 double * ttlookup,
                                 const double * left_matrix,
                                 const double * right_matrix,
                                 const unsigned int * tipmap,
                                 unsigned int tipmap_size)
{
  if (states == 4)
  {
    pll_core_create_lookup_4x4_neon(rate_cats,
                                    ttlookup,
                                    left_matrix,
                                    right_matrix);
    return;
  }

  unsigned int i,j,k,n,m;
  unsigned int states_padded = (states+1) & 0xFFFFFFFE;
  unsigned int maxstates = tipmap_size;

  float64x2_t xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6;

  unsigned int log2_maxstates = (unsigned int)ceil(log2(maxstates));
  unsigned int span_padded = states_padded*rate_cats;

  size_t displacement = (states_padded - states) * (states_padded);

  /* precompute first the entries that contain only one 1 */
  const double * jmat;
  const double * kmat;
  double * lookup;

  /* go through all pairs j,k of states for the two tips; i is the inner
     node state */
  xmm0 = vdupq_n_f64(0);

  for (j = 0; j < maxstates; ++j)
  {
    unsigned int jstate = tipmap[j];
    for (k = 0; k < maxstates; ++k)
    {
      unsigned int kstate = tipmap[k];
      jmat = left_matrix;
      kmat = right_matrix;

      /* find offset of state-pair in the precomputation table */
      lookup = ttlookup;
      lookup += ((j << log2_maxstates) + k)*span_padded;

      /* precompute the likelihood for each state and each rate */
      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states_padded; i += 2)
        {
          float64x2_t v_termj0 = vdupq_n_f64(0);
          float64x2_t v_termj1 = vdupq_n_f64(0);
          float64x2_t v_termk0 = vdupq_n_f64(0);
          float64x2_t v_termk1 = vdupq_n_f64(0);

          const double * jm0 = jmat;
          const double * jm1 = jm0 + states_padded;

          const double * km0 = kmat;
          const double * km1 = km0 + states_padded;

          /* set position of least significant bit in character state */
          register int lsb = 0;

          /* decompose basecall into the encoded residues and set the appropriate
             positions in the tip vector */
          for (m = 0; m < states_padded; m += 2)
          {
            /* set mask for left matrix */
            PLL_ALIGN_HEADER(16) double mask[2] PLL_ALIGN_FOOTER(16) =
             {
              ((jstate >> (lsb+0)) & 1) ? 1 : 0,
              ((jstate >> (lsb+1)) & 1) ? 1 : 0
             };
            xmm1 = vld1q_f64(mask);
            xmm2 = vreinterpretq_f64_u64(vcgtq_f64(xmm1,xmm0));

            /* set mask for right matrix */
            mask[0] = ((kstate >> (lsb+0)) & 1) ? 1 : 0;
            mask[1] = ((kstate >> (lsb+1)) & 1) ? 1 : 0;
            xmm1 = vld1q_f64(mask);
            xmm3 = vreinterpretq_f64_u64(vcgtq_f64(xmm1,xmm0));

            lsb += 2;

            /* load row 0 */
            xmm4 = vld1q_f64(jm0);
            xmm5 = vreinterpretq_f64_s64(
                     vandq_s64(vreinterpretq_s64_f64(xmm4),
                               vreinterpretq_s64_f64(xmm2)));
            v_termj0 = vaddq_f64(v_termj0,xmm5);

            xmm4 = vld1q_f64(km0);
            xmm5 = vreinterpretq_f64_s64(
                     vandq_s64(vreinterpretq_s64_f64(xmm4),
                               vreinterpretq_s64_f64(xmm3)));
            v_termk0 = vaddq_f64(v_termk0,xmm5);

            jm0 += 2;
            km0 += 2;

            /* load row 1 */
            xmm4 = vld1q_f64(jm1);
            xmm5 = vreinterpretq_f64_s64(
                     vandq_s64(vreinterpretq_s64_f64(xmm4),
                               vreinterpretq_s64_f64(xmm2)));
            v_termj1 = vaddq_f64(v_termj1,xmm5);

            xmm4 = vld1q_f64(km1);
            xmm5 = vreinterpretq_f64_s64(
                     vandq_s64(vreinterpretq_s64_f64(xmm4),
                               vreinterpretq_s64_f64(xmm3)));
            v_termk1 = vaddq_f64(v_termk1,xmm5);

            jm1 += 2;
            km1 += 2;
          }

          jmat = jm1;
          kmat = km1;

          xmm4 = vpaddq_f64(v_termj0,v_termj1);
          xmm5 = vpaddq_f64(v_termk0,v_termk1);
          xmm6 = vmulq_f64(xmm4,xmm5);
          vst1q_f64(lookup+i,xmm6);
        }

        /* reset pointers to the start of the next p-matrix, as the vectorization
           assumes a square states_padded * states_padded matrix, even though the
           real matrix is states * states_padded */
        jmat -= displacement;
        kmat -= displacement;
        ///* this is to avoid valgrind warnings on accessing uninitialized memory
        //   when using SSE and states are not a multiple of 2 */
        //if (states_padded-states)
        //  memset(lookup+index, 0, (states_padded-states)*sizeof(double));

        lookup += states_padded;
      }
    }
  }
}

#if 0
static void pprint_sse(float64x2_t x)
{
  double * p = (double *) & x;

  printf("%f ", *p++);
  printf("%f ", *p++);
}

static void pshow_sse(char * name, float64x2_t x)
{
  printf("%s: ", name);
  pprint_sse(x);
  printf("\n");
}
#endif

void pll_core_create_lookup_4x4_neon(unsigned int rate_cats,
                                     double * lookup,
                                     const double * left_matrix,
                                     const double * right_matrix)
{
  unsigned int j,k,n;
  unsigned int maxstates = 16;

  float64x2_t ymm0,ymm1,ymm2,ymm3,ymm4;
  float64x2_t xmm4,xmm5,xmm6,xmm7;
  float64x2_t xmm0,xmm1,xmm2,xmm3;

  const double * jmat;
  const double * kmat;

  ymm4 = vdupq_n_f64(0);

  /* skip entries for j = 0 */
  lookup += maxstates*4*rate_cats;

  for (j = 1; j < maxstates; ++j)
  {
    /* masks for state j */
    PLL_ALIGN_HEADER(16) double mask[2] PLL_ALIGN_FOOTER(16) =
     {
       ((j >> 0) & 1) ? 1 : 0,
       ((j >> 1) & 1) ? 1 : 0
     };
    xmm1 = vld1q_f64(mask);
    xmm0 = vreinterpretq_f64_u64(vcgtq_f64(xmm1,ymm4));

    mask[0] = ((j >> 2) & 1) ? 1 : 0;
    mask[1] = ((j >> 3) & 1) ? 1 : 0;
    xmm2 = vld1q_f64(mask);
    xmm1 = vreinterpretq_f64_u64(vcgtq_f64(xmm2,ymm4));

    /* skip entry for k = 0 */
    lookup += 4*rate_cats;

    for (k = 1; k < maxstates; ++k)
    {
      /* masks for state k */
      mask[0] = ((k >> 0) & 1) ? 1 : 0;
      mask[1] = ((k >> 1) & 1) ? 1 : 0;
      xmm3 = vld1q_f64(mask);
      xmm2 = vreinterpretq_f64_u64(vcgtq_f64(xmm3,ymm4));

      mask[0] = ((k >> 2) & 1) ? 1 : 0;
      mask[1] = ((k >> 3) & 1) ? 1 : 0;
      xmm4 = vld1q_f64(mask);
      xmm3 = vreinterpretq_f64_u64(vcgtq_f64(xmm4,ymm4));

      jmat = left_matrix;
      kmat = right_matrix;

      for (n = 0; n < rate_cats; ++n)
      {
        /* load row0 from left matrix  */
        ymm0 = vld1q_f64(jmat);
        ymm1 = vld1q_f64(jmat+2);

        /* mask row */
        xmm4 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm0)));
        xmm5 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm1)));

        /* point to row1 */
        jmat += 4;

        /* load row1 from left matrix */
        ymm0 = vld1q_f64(jmat);
        ymm1 = vld1q_f64(jmat+2);

        /* mask row */
        xmm6 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm0)));
        xmm7 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm1)));

        /* point to row2 */
        jmat += 4;

        /* horizontally add the two rows */
        ymm0 = vpaddq_f64(xmm4,xmm5);
        ymm1 = vpaddq_f64(xmm6,xmm7);
        
        /* create vector containing sums of row0 and row1 (left matrix) */
        ymm2 = vpaddq_f64(ymm0,ymm1);

        /* load row0 from right matrix  */
        ymm0 = vld1q_f64(kmat);
        ymm1 = vld1q_f64(kmat+2);

        /* mask row */
        xmm4 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm2)));
        xmm5 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm3)));

        /* point to row1 */
        kmat += 4;

        /* load row1 from left matrix */
        ymm0 = vld1q_f64(kmat);
        ymm1 = vld1q_f64(kmat+2);

        /* mask row */
        xmm6 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm2)));
        xmm7 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm3)));

        /* point to row2 */
        kmat += 4;

        /* horizontally add the two rows */
        ymm0 = vpaddq_f64(xmm4,xmm5);
        ymm1 = vpaddq_f64(xmm6,xmm7);
        
        /* create vector containing sums of row0 and row1 (right matrix) */
        ymm3 = vpaddq_f64(ymm0,ymm1);

        /* multiply the sums from left and right matrix */
        ymm0 = vmulq_f64(ymm2,ymm3);
        vst1q_f64(lookup,ymm0);

        /* point to the next two entries to fill */
        lookup += 2;

        /* load row2 from left matrix  */
        ymm0 = vld1q_f64(jmat);
        ymm1 = vld1q_f64(jmat+2);

        /* mask row */
        xmm4 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm0)));
        xmm5 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm1)));

        /* point to row1 */
        jmat += 4;

        /* load row3 from left matrix */
        ymm0 = vld1q_f64(jmat);
        ymm1 = vld1q_f64(jmat+2);

        /* mask row */
        xmm6 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm0)));
        xmm7 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm1)));

        /* point to row0 of next p-matrix */
        jmat += 4;

        /* horizontally add the two rows */
        ymm0 = vpaddq_f64(xmm4,xmm5);
        ymm1 = vpaddq_f64(xmm6,xmm7);
        
        /* create vector containing sums of row2 and row3 (left matrix) */
        ymm2 = vpaddq_f64(ymm0,ymm1);

        /* load row2 from right matrix  */
        ymm0 = vld1q_f64(kmat);
        ymm1 = vld1q_f64(kmat+2);

        /* mask row */
        xmm4 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm2)));
        xmm5 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm3)));

        /* point to row3 */
        kmat += 4;

        /* load row3 from left matrix */
        ymm0 = vld1q_f64(kmat);
        ymm1 = vld1q_f64(kmat+2);

        /* mask row */
        xmm6 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm0),
                           vreinterpretq_s64_f64(xmm2)));
        xmm7 = vreinterpretq_f64_s64(
                 vandq_s64(vreinterpretq_s64_f64(ymm1),
                           vreinterpretq_s64_f64(xmm3)));

        /* point to row0 of next p-matrix */
        kmat += 4;

        /* horizontally add the two rows */
        ymm0 = vpaddq_f64(xmm4,xmm5);
        ymm1 = vpaddq_f64(xmm6,xmm7);
        
        /* create vector containing sums of row0 and row1 (right matrix) */
        ymm3 = vpaddq_f64(ymm0,ymm1);

        /* multiply the sums from left and right matrix */
        ymm0 = vmulq_f64(ymm2,ymm3);
        vst1q_f64(lookup,ymm0);

        /* point to the beginning of next state pair */
        lookup += 2;
      }
    }
  }
}

void pll_core_update_partial_tt_neon(unsigned int states,
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
  unsigned int states_padded = (states+1) & 0xFFFFFFFE;
  unsigned int span_padded = states_padded * rate_cats;
  const double * offset;

  if (states == 4)
  {
    pll_core_update_partial_tt_4x4_neon(sites,
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

void pll_core_update_partial_tt_4x4_neon(unsigned int sites,
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

    memcpy(parent_clv, offset, span*sizeof(double));

    parent_clv += span;
  }
}

void pll_core_update_partial_ii_4x4_neon(unsigned int sites,
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
  unsigned int span = states * rate_cats;
  unsigned int n,k,i;

  const double * lmat;
  const double * rmat;

  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;

  float64x2_t v_scale_threshold = vdupq_n_f64(PLL_SCALE_THRESHOLD);
  float64x2_t v_scale_factor = vdupq_n_f64(PLL_SCALE_FACTOR);

  float64x2_t xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8,xmm9,xmm10;

  if (parent_scaler)
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0x3 : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }
  else
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }

  /* 
     perform the following matrix multiplications:

  +-----+-----+-----+-----+     +-----+-----+-----+-----+
  | a1  | a2  | a3  | a4  |     | b1  | b2  | b3  | b4  |
  +-----+-----+-----+-----+     +-----+-----+-----+-----+
  | a5  | a6  | a7  | a8  |     | b5  | b6  | b7  | b8  |
  +-----+-----+-----+-----+     +-----+-----+-----+-----+
  | a9  | a10 | a11 | a12 |     | b9  | b10 | b11 | b12 |
  +-----+-----+-----+-----+     +-----+-----+-----+-----+
  | a13 | a14 | a15 | a16 |     | b13 | b14 | b15 | b16 |
  +-----+-----+-----+-----+     +-----+-----+-----+-----+
    
              x                             x

    +----+----+----+----+         +----+----+----+----+
    | c1 | c2 | c3 | c4 |         | d1 | d2 | d3 | d4 |
    +----+----+----+----+         +----+----+----+----+

  */

  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      /* do the computation on the first two rows of the two matrices */

      /* compute left */
      xmm0 = vld1q_f64(left_clv);          /* needed */
      xmm2 = vld1q_f64(lmat);
      xmm3 = vmulq_f64(xmm0,xmm2);

      xmm1 = vld1q_f64(left_clv+2);        /* needed */
      xmm2 = vld1q_f64(lmat+2);
      xmm4 = vmulq_f64(xmm1,xmm2);

      /* calculate (a1*c1 + a3*c3 | a2*c2 + a4*c4) */
      xmm2 = vaddq_f64(xmm3,xmm4);         /* needed (1) */

      /* compute right */
      xmm3 = vld1q_f64(right_clv);         /* needed */
      xmm4 = vld1q_f64(rmat);
      xmm6 = vmulq_f64(xmm3,xmm4);

      xmm4 = vld1q_f64(right_clv+2);       /* needed */
      xmm7 = vld1q_f64(rmat+2);
      xmm8 = vmulq_f64(xmm4,xmm7);

      /* calculate (b1*d1 + b3*d3 | b2*d2 + b4*d4) */
      xmm5 = vaddq_f64(xmm6,xmm8);         /* needed (2) */

      rmat += states;
      lmat += states;

      /* compute left */
      xmm6 = vld1q_f64(lmat);
      xmm7 = vmulq_f64(xmm0,xmm6);

      xmm6 = vld1q_f64(lmat+2);
      xmm8 = vmulq_f64(xmm1,xmm6);

      /* calculate (a5*c1 + a7*c3 | a6*c2 + a8*c4) */
      xmm6 = vaddq_f64(xmm7,xmm8);         /* needed (3) */

      /* compute right */
      xmm7 = vld1q_f64(rmat);
      xmm8 = vmulq_f64(xmm3,xmm7);

      xmm7 = vld1q_f64(rmat+2);
      xmm9 = vmulq_f64(xmm4,xmm7);
      
      /* calculate (b5*d1 + b7*d3 | b6*d2 + b8*d4) */
      xmm7 = vaddq_f64(xmm8,xmm9);         /* needed (4) */
      
      xmm8 = vpaddq_f64(xmm2,xmm6);
      xmm2 = vpaddq_f64(xmm5,xmm7);
      xmm9 = vmulq_f64(xmm8,xmm2);

      rmat += states;
      lmat += states;

      /* do the computation on the last two rows of the two matrices */

      /* compute left */
      xmm5 = vld1q_f64(lmat);
      xmm6 = vmulq_f64(xmm0,xmm5);

      xmm5 = vld1q_f64(lmat+2);
      xmm7 = vmulq_f64(xmm1,xmm5);

      /* calculate (a9*c1 + a11*c3 | a10*c2 + a12*c4) */
      xmm2 = vaddq_f64(xmm6,xmm7);         /* needed (1) */

      /* compute right */
      xmm6 = vld1q_f64(rmat);
      xmm7 = vmulq_f64(xmm3,xmm6);

      xmm6 = vld1q_f64(rmat+2);
      xmm8 = vmulq_f64(xmm4,xmm6);
      
      /* calculate (b9*d1 + b11*d3 | b10*d2 + b12*d4) */
      xmm5 = vaddq_f64(xmm7,xmm8);         /* needed (2) */
      
      rmat += states;
      lmat += states;

      /* compute left */
      xmm6 = vld1q_f64(lmat);
      xmm7 = vmulq_f64(xmm0,xmm6);

      xmm6 = vld1q_f64(lmat+2);
      xmm8 = vmulq_f64(xmm1,xmm6);

      /* calculate (a13*c1 + a15*c3 | a14*c2 + a16*c4) */
      xmm6 = vaddq_f64(xmm7,xmm8);         /* needed (3) */

      /* compute right */
      xmm7 = vld1q_f64(rmat);
      xmm8 = vmulq_f64(xmm3,xmm7);

      xmm7 = vld1q_f64(rmat+2);
      xmm10 = vmulq_f64(xmm4,xmm7);

      /* calculate (b13*d1 + b15*d3 | b14*d2 + b16*d4) */
      xmm7 = vaddq_f64(xmm8,xmm10);         /* needed (4) */

      rmat += states;
      lmat += states;
      
      xmm8 = vpaddq_f64(xmm2,xmm6);
      xmm2 = vpaddq_f64(xmm5,xmm7);
      xmm10 = vmulq_f64(xmm8,xmm2);

      /* check if scaling is needed for the current rate category */
      float64x2_t v_cmp = vreinterpretq_f64_u64(vcltq_f64(xmm9, v_scale_threshold));
      #if 0
      unsigned int rate_mask = _mm_movemask_pd(v_cmp);
      #else
      uint64x2_t high_bits = vshrq_n_u64(vreinterpretq_u64_f64(v_cmp),63);
      unsigned int rate_mask = vgetq_lane_u64(high_bits,0) |
                               (vgetq_lane_u64(high_bits,1) << 1);
      #endif
      v_cmp = vreinterpretq_f64_u64(vcltq_f64(xmm10, v_scale_threshold));
      #if 0
      rate_mask = rate_mask & _mm_movemask_pd(v_cmp);
      #else
      high_bits = vshrq_n_u64(vreinterpretq_u64_f64(v_cmp),63);
      rate_mask &= vgetq_lane_u64(high_bits,0) |
                   (vgetq_lane_u64(high_bits,1) << 1);
      #endif

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0x3)
        {
          xmm9 = vmulq_f64(xmm9,v_scale_factor);
          xmm10 = vmulq_f64(xmm10,v_scale_factor);
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      vst1q_f64(parent_clv, xmm9);
      vst1q_f64(parent_clv+2, xmm10);

      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0x3)
    {
      parent_clv -= span;
      for (i = 0; i < span; i += 2)
      {
        float64x2_t v_prod = vld1q_f64(parent_clv + i);
        v_prod = vmulq_f64(v_prod,v_scale_factor);
        vst1q_f64(parent_clv + i, v_prod);
      }
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

void pll_core_update_partial_ii_neon(unsigned int states,
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

  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  /* dedicated functions for 4x4 matrices */
  if (states == 4)
  {
    pll_core_update_partial_ii_4x4_neon(sites,
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

  unsigned int span_padded = states_padded * rate_cats;
  size_t displacement = (states_padded - states) * (states_padded);

  float64x2_t xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6;

  /* scaling stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  float64x2_t v_scale_threshold = vdupq_n_f64(PLL_SCALE_THRESHOLD);
  float64x2_t v_scale_factor = vdupq_n_f64(PLL_SCALE_FACTOR);

  if (parent_scaler)
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0x3 : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }
  else
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0x3;

      for (i = 0; i < states_padded; i += 2)
      {
        float64x2_t v_terma0 = vdupq_n_f64(0);
        float64x2_t v_terma1 = vdupq_n_f64(0);
        float64x2_t v_termb0 = vdupq_n_f64(0);
        float64x2_t v_termb1 = vdupq_n_f64(0);

        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;

        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;

        for (j = 0; j < states_padded; j += 2)
        {
          /* load left and right clvs */
          xmm0 = vld1q_f64(left_clv+j);
          xmm1 = vld1q_f64(right_clv+j);

          /* row 0 */
          xmm2 = vld1q_f64(lm0);
          xmm3 = vmulq_f64(xmm2,xmm0);
          v_terma0 = vaddq_f64(v_terma0, xmm3);

          xmm2 = vld1q_f64(rm0);
          xmm3 = vmulq_f64(xmm2,xmm1);
          v_termb0 = vaddq_f64(v_termb0,xmm3);

          lm0 += 2;
          rm0 += 2;

          /* row 1 */
          xmm2 = vld1q_f64(lm1);
          xmm3 = vmulq_f64(xmm2,xmm0);
          v_terma1 = vaddq_f64(v_terma1,xmm3);

          xmm2 = vld1q_f64(rm1);
          xmm3 = vmulq_f64(xmm2,xmm1);
          v_termb1 = vaddq_f64(v_termb1,xmm3);

          lm1 += 2;
          rm1 += 2;
        }
        
        lmat = lm1;
        rmat = rm1;

        xmm4 = vpaddq_f64(v_terma0,v_terma1);
        xmm5 = vpaddq_f64(v_termb0,v_termb1);
        xmm6 = vmulq_f64(xmm4,xmm5);

        /* check if scaling is needed for the current rate category */
        float64x2_t v_cmp = vreinterpretq_f64_u64(vcltq_f64(xmm6, v_scale_threshold));
        #if 0
        rate_mask = rate_mask & _mm_movemask_pd(v_cmp);
        #else
        uint64x2_t high_bits = vshrq_n_u64(vreinterpretq_u64_f64(v_cmp),63);
        rate_mask &= vgetq_lane_u64(high_bits,0) |
                     (vgetq_lane_u64(high_bits,1) << 1);
        #endif

        vst1q_f64(parent_clv+i,xmm6);
      }

      /* reset pointers to the start of the next p-matrix, as the vectorization
         assumes a square states_padded * states_padded matrix, even though the
         real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0x3)
        {
          for (i = 0; i < states_padded; i += 2)
          {
            float64x2_t v_prod = vld1q_f64(parent_clv + i);
            v_prod = vmulq_f64(v_prod, v_scale_factor);
            vst1q_f64(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      parent_clv += states_padded;
      left_clv   += states_padded;
      right_clv  += states_padded;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0x3)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 2)
      {
        float64x2_t v_prod = vld1q_f64(parent_clv + i);
        v_prod = vmulq_f64(v_prod,v_scale_factor);
        vst1q_f64(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

void pll_core_update_partial_ti_4x4_neon(unsigned int sites,
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
  unsigned int span = states * rate_cats;
  unsigned int i,k,n;
  unsigned int lstate;

  const double * lmat;
  const double * rmat;

  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;

  float64x2_t v_scale_threshold = vdupq_n_f64(PLL_SCALE_THRESHOLD);
  float64x2_t v_scale_factor = vdupq_n_f64(PLL_SCALE_FACTOR);

  float64x2_t xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
  float64x2_t ymm0,ymm1,ymm2,ymm3;

  PLL_ALIGN_HEADER(16) double mask[2] PLL_ALIGN_FOOTER(16) = {0,0};
  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(64*rate_cats*sizeof(double),
                                      PLL_ALIGNMENT_NEON);
  if (!lookup)
    fatal("Cannot allocate space for precomputation.");

  /* skip first entry of lookup table as it is never used */
  double * ptr = lookup + span;

  ymm3 = vdupq_n_f64(0);

  /* iterate all ambiguities skipping 0 */
  for (i = 1; i < 16; ++i)
  {
    lmat = left_matrix;

    /* mask the entries of pmatrix row to be loaded */
    mask[0] = ((i >> 0) & 1) ? 1 : 0;
    mask[1] = ((i >> 1) & 1) ? 1 : 0;
    xmm1 = vld1q_f64(mask);
    xmm0 = vreinterpretq_f64_u64(vcgtq_f64(xmm1,ymm3));

    mask[0] = ((i >> 2) & 1) ? 1 : 0;
    mask[1] = ((i >> 3) & 1) ? 1 : 0;
    xmm2 = vld1q_f64(mask);
    xmm1 = vreinterpretq_f64_u64(vcgtq_f64(xmm2,ymm3));

    for (k = 0; k < rate_cats; ++k)
    {
      /* load row0 from matrix */
      ymm0 = vld1q_f64(lmat);
      ymm1 = vld1q_f64(lmat+2);

      /* mask row0 from matrix */
      xmm4 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm0),
                         vreinterpretq_s64_f64(xmm0)));
      xmm5 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm1),
                         vreinterpretq_s64_f64(xmm1)));
      xmm6 = vaddq_f64(xmm4,xmm5);     /* a1+a3 | a2+a4 */

      lmat += 4;

      /* load row1 from left matrix */
      ymm0 = vld1q_f64(lmat);
      ymm1 = vld1q_f64(lmat+2);

      /* mask row */
      xmm4 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm0),
                         vreinterpretq_s64_f64(xmm0)));
      xmm5 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm1),
                         vreinterpretq_s64_f64(xmm1)));
      xmm7 = vaddq_f64(xmm4,xmm5);     /* a5+a7 | a3+a8 */

      ymm2 = vpaddq_f64(xmm6,xmm7);
      vst1q_f64(ptr,ymm2);

      /* point to row2 */
      lmat += 4;

      /* load row2 from left matrix */
      ymm0 = vld1q_f64(lmat);
      ymm1 = vld1q_f64(lmat+2);

      /* mask row */
      xmm4 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm0),
                         vreinterpretq_s64_f64(xmm0)));
      xmm5 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm1),
                         vreinterpretq_s64_f64(xmm1)));
      xmm6 = vaddq_f64(xmm4,xmm5);     /* a9+a11 | a10+a12 */

      lmat += 4;

      /* load row3 from left matrix */
      ymm0 = vld1q_f64(lmat);
      ymm1 = vld1q_f64(lmat+2);

      /* mask row */
      xmm4 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm0),
                         vreinterpretq_s64_f64(xmm0)));
      xmm5 = vreinterpretq_f64_s64(
               vandq_s64(vreinterpretq_s64_f64(ymm1),
                         vreinterpretq_s64_f64(xmm1)));
      xmm7 = vaddq_f64(xmm4,xmm5);     /* a13+a15 | a14+a16 */

      ymm2 = vpaddq_f64(xmm6,xmm7);
      vst1q_f64(ptr+2,ymm2);

      /* move pointers */
      ptr  += 4;
      lmat += 4;
    }
  }

  if (parent_scaler)
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0x3 : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* update the parent scaler with the scaler of the right child */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }
  else
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }

  /* iterate over sites and compute CLV entries */
  for (n = 0; n < sites; ++n)
  {
    rmat = right_matrix;

    scale_mask = init_mask;

    lstate = left_tipchar[n];

    unsigned int loffset = lstate*span;

    for (k = 0; k < rate_cats; ++k)
    {
      /* load right child CLV */
      xmm0 = vld1q_f64(right_clv);
      xmm1 = vld1q_f64(right_clv+2);

      /* load right pmatrix row0 */
      xmm2 = vld1q_f64(rmat);
      xmm3 = vld1q_f64(rmat+2);

      xmm4 = vmulq_f64(xmm0,xmm2);
      xmm5 = vmulq_f64(xmm1,xmm3);
      xmm6 = vaddq_f64(xmm4,xmm5);    /* a1*c1 + a3*c3 | a2*c2 + a4*c4 */

      rmat += states;

      /* load right pmatrix row1 */
      xmm2 = vld1q_f64(rmat);
      xmm3 = vld1q_f64(rmat+2);

      xmm4 = vmulq_f64(xmm0,xmm2);
      xmm5 = vmulq_f64(xmm1,xmm3);
      xmm7 = vaddq_f64(xmm4,xmm5);    /* a5*c1 + a7*c3 | a6*c2 + a8*c4 */

      rmat += states;

      /* create a1*c2+a2*c2+a3*c3+a4*c4 | a5*c1+a6*c2+a7*c3+a8*c4 */
      xmm4 = vpaddq_f64(xmm6,xmm7);

      /* load precomputed lookup table into xmm2 */
      xmm2 = vld1q_f64(lookup+loffset);
      xmm3 = vmulq_f64(xmm4,xmm2);

      /* load right pmatrix row2 */
      xmm2 = vld1q_f64(rmat);
      xmm8 = vld1q_f64(rmat+2);

      xmm4 = vmulq_f64(xmm0,xmm2);
      xmm5 = vmulq_f64(xmm1,xmm8);
      xmm6 = vaddq_f64(xmm4,xmm5);    /* a1*c1 + a3*c3 | a2*c2 + a4*c4 */

      rmat += states;

      /* load right pmatrix row3 */
      xmm2 = vld1q_f64(rmat);
      xmm8 = vld1q_f64(rmat+2);

      xmm4 = vmulq_f64(xmm0,xmm2);
      xmm5 = vmulq_f64(xmm1,xmm8);
      xmm7 = vaddq_f64(xmm4,xmm5);    /* a5*c1 + a7*c3 | a6*c2 + a8*c4 */

      rmat += states;

      /* create a1*c2+a2*c2+a3*c3+a4*c4 | a5*c1+a6*c2+a7*c3+a8*c4 */
      xmm4 = vpaddq_f64(xmm6,xmm7);

      /* load precomputed lookup table into xmm2 */
      xmm2 = vld1q_f64(lookup+loffset+2);
      xmm8 = vmulq_f64(xmm4,xmm2);

//      for (i = 0; i < states; ++i)
//        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);

      /* check if scaling is needed for the current rate category */
      float64x2_t v_cmp = vreinterpretq_f64_u64(vcltq_f64(xmm3, v_scale_threshold));
      #if 0
      unsigned int rate_mask = _mm_movemask_pd(v_cmp);
      #else
      uint64x2_t high_bits = vshrq_n_u64(vreinterpretq_u64_f64(v_cmp),63);
      unsigned int rate_mask = vgetq_lane_u64(high_bits,0) |
                               (vgetq_lane_u64(high_bits,1) << 1);
      #endif
      v_cmp = vreinterpretq_f64_u64(vcltq_f64(xmm8, v_scale_threshold));
      #if 0
      rate_mask = rate_mask & _mm_movemask_pd(v_cmp);
      #else
      high_bits = vshrq_n_u64(vreinterpretq_u64_f64(v_cmp),63);
      rate_mask &= vgetq_lane_u64(high_bits,0) |
                   (vgetq_lane_u64(high_bits,1) << 1);
      #endif

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0x3)
        {
          xmm3 = vmulq_f64(xmm3,v_scale_factor);
          xmm8 = vmulq_f64(xmm8,v_scale_factor);
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      vst1q_f64(parent_clv,xmm3);
      vst1q_f64(parent_clv+2,xmm8);

      parent_clv += states;
      right_clv  += states;
      loffset    += 4;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0x3)
    {
      parent_clv -= span;
      for (i = 0; i < span; i += 2)
      {
        float64x2_t v_prod = vld1q_f64(parent_clv + i);
        v_prod = vmulq_f64(v_prod,v_scale_factor);
        vst1q_f64(parent_clv + i, v_prod);
      }
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
  pll_aligned_free(lookup);
}

void pll_core_update_partial_ti_neon(unsigned int states,
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
  unsigned int states_padded = (states+1) & 0xFFFFFFFE;
  unsigned int span_padded = states_padded * rate_cats;

  const double * lmat;
  const double * rmat;

  unsigned int lstate;

  /* dedicated functions for 4x4 matrices */
  if (states == 4)
  {
    pll_core_update_partial_ti_4x4_neon(sites,
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

  size_t displacement = (states_padded - states) * (states_padded);
  float64x2_t xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

  xmm7 = vdupq_n_f64(0);

  /* scaling stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  float64x2_t v_scale_threshold = vdupq_n_f64(PLL_SCALE_THRESHOLD);
  float64x2_t v_scale_factor = vdupq_n_f64(PLL_SCALE_FACTOR);

  if (parent_scaler)
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0x3 : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }
  else
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    lstate = tipmap[left_tipchars[n]];

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0x3;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 2)
      {
        float64x2_t v_terma0 = vdupq_n_f64(0);
        float64x2_t v_terma1 = vdupq_n_f64(0);
        float64x2_t v_termb0 = vdupq_n_f64(0);
        float64x2_t v_termb1 = vdupq_n_f64(0);

        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;

        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;

        /* set position of least significant bit in character state */
        register int lsb = 0;

        for (j = 0; j < states_padded; j += 2)
        {
          /* set mask */
          PLL_ALIGN_HEADER(16) double mask[2] PLL_ALIGN_FOOTER(16) =
           {
             ((lstate >> (lsb+0)) & 1) ? 1 : 0,
             ((lstate >> (lsb+1)) & 1) ? 1 : 0,
           };
          xmm1 = vld1q_f64(mask);
          xmm0 = vreinterpretq_f64_u64(vcgtq_f64(xmm1,xmm7));

          lsb += 2;

          /* load clv */
          xmm2 = vld1q_f64(right_clv+j);

          /* row 0 */
          xmm3 = vld1q_f64(lm0);
          xmm4 = vreinterpretq_f64_s64(
                   vandq_s64(vreinterpretq_s64_f64(xmm3),
                             vreinterpretq_s64_f64(xmm0)));
          v_terma0 = vaddq_f64(v_terma0,xmm4);

          xmm3 = vld1q_f64(rm0);
          xmm4 = vmulq_f64(xmm3,xmm2);
          v_termb0 = vaddq_f64(v_termb0,xmm4);

          lm0 += 2;
          rm0 += 2;

          /* row 1 */
          xmm3 = vld1q_f64(lm1);
          xmm4 = vreinterpretq_f64_s64(
                   vandq_s64(vreinterpretq_s64_f64(xmm3),
                             vreinterpretq_s64_f64(xmm0)));
          v_terma1 = vaddq_f64(v_terma1,xmm4);

          xmm3 = vld1q_f64(rm1);
          xmm4 = vmulq_f64(xmm3,xmm2);
          v_termb1 = vaddq_f64(v_termb1,xmm4);

          lm1 += 2;
          rm1 += 2;
        }

        lmat = lm1;
        rmat = rm1;

        xmm4 = vpaddq_f64(v_terma0,v_terma1);
        xmm5 = vpaddq_f64(v_termb0,v_termb1);
        xmm6 = vmulq_f64(xmm4,xmm5);

        /* check if scaling is needed for the current rate category */
        float64x2_t v_cmp = vreinterpretq_f64_u64(vcltq_f64(xmm6, v_scale_threshold));
        #if 0
        rate_mask = rate_mask & _mm_movemask_pd(v_cmp);
        #else
        uint64x2_t high_bits = vshrq_n_u64(vreinterpretq_u64_f64(v_cmp),63);
        rate_mask &= vgetq_lane_u64(high_bits,0) |
                     (vgetq_lane_u64(high_bits,1) << 1);
        #endif

        vst1q_f64(parent_clv+i,xmm6);
      }

      /* reset pointers to the start of the next p-matrix, as the vectorization
         assumes a square states_padded * states_padded matrix, even though the
         real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0x3)
        {
          for (i = 0; i < states_padded; i += 2)
          {
            float64x2_t v_prod = vld1q_f64(parent_clv + i);
            v_prod = vmulq_f64(v_prod, v_scale_factor);
            vst1q_f64(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      parent_clv += states_padded;
      right_clv  += states_padded;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0x3)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 2)
      {
        float64x2_t v_prod = vld1q_f64(parent_clv + i);
        v_prod = vmulq_f64(v_prod,v_scale_factor);
        vst1q_f64(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

#endif
