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

void pll_core_update_partial_ii(unsigned int states,
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

  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int site_scale;
  unsigned int init_mask;

  const double * lmat;
  const double * rmat;

  unsigned int span = states * rate_cats;

  /* init scaling-related stuff */
  if (parent_scaler)
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 1 : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;

    /* add up the scale vectors of the two children if available */
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
    site_scale = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_scale = 1;
      for (i = 0; i < states; ++i)
      {
        double terma = 0;
        double termb = 0;
        for (j = 0; j < states; ++j)
        {
          terma += lmat[j] * left_clv[j];
          termb += rmat[j] * right_clv[j];
        }
        parent_clv[i] = terma*termb;

        rate_scale &= (parent_clv[i] < PLL_SCALE_THRESHOLD);

        lmat += states;
        rmat += states;
      }

      /* check if scaling is needed for the current rate category */
      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_scale)
        {
          for (i = 0; i < states; ++i)
            parent_clv[i] *= PLL_SCALE_FACTOR;
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        site_scale = site_scale && rate_scale;

      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }
    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (site_scale)
    {
      parent_clv -= span;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}
