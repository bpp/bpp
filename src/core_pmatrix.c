/*
    Copyright (C) 2016-2018 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

int pll_core_update_pmatrix_4x4_jc69(double ** pmatrix,
                                     unsigned int states,
                                     unsigned int rate_cats,
                                     const double * rates,
                                     const double * branch_lengths,
                                     const unsigned int * matrix_indices,
                                     const unsigned int * params_indices,
                                     unsigned int count,
                                     unsigned int attrib)
{
  unsigned int i,n;
  unsigned int states_padded = states;

  double * pmat;


  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;

      double t = branch_lengths[i];

      if (t < 1e-100)
      {
        pmat[0]  = 1;
        pmat[1]  = 0;
        pmat[2]  = 0;
        pmat[3]  = 0;

        pmat[4]  = 0;
        pmat[5]  = 1;
        pmat[6]  = 0;
        pmat[7]  = 0;

        pmat[8]  = 0;
        pmat[9]  = 0;
        pmat[10] = 1;
        pmat[11] = 0;

        pmat[12] = 0;
        pmat[13] = 0;
        pmat[14] = 0;
        pmat[15] = 1;
      }
      else
      {
        double a =  (1 + 3*exp(-4*t/3) ) / 4;
        double b = (1 - a) / 3;

        pmat[0]  = a;
        pmat[1]  = b;
        pmat[2]  = b;
        pmat[3]  = b;

        pmat[4]  = b;
        pmat[5]  = a;
        pmat[6]  = b;
        pmat[7]  = b;

        pmat[8]  = b;
        pmat[9]  = b;
        pmat[10] = a;
        pmat[11] = b;

        pmat[12] = b;
        pmat[13] = b;
        pmat[14] = b;
        pmat[15] = a;
      }

    }
  }

  return BPP_SUCCESS;
}

int pll_core_update_pmatrix(double ** pmatrix,
                            unsigned int states,
                            unsigned int rate_cats,
                            const double * rates,
                            const double * branch_lengths,
                            const unsigned int * matrix_indices,
                            const unsigned int * params_indices,
                            double * const * eigenvals,
                            double * const * eigenvecs,
                            double * const * inv_eigenvecs,
                            unsigned int count,
                            unsigned int attrib)
{
  unsigned int i,n,j,k,m;
  unsigned int states_padded = states;
  double * expd;
  double * temp;

  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;


  expd = (double *)xmalloc(states * sizeof(double));
  temp = (double *)xmalloc(states*states*sizeof(double));

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;

      evecs = eigenvecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals = eigenvals[params_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (!branch_lengths[i])
      {
        for (j = 0; j < states; ++j)
          for (k = 0; k < states; ++k)
            pmat[j*states_padded + k] = (j == k) ? 1 : 0;
      }
      else
      {
        /* NOTE: in order to deal with numerical issues in cases when Qt -> 0, we
         * use a trick suggested by Ben Redelings and explained here:
         * https://github.com/xflouris/libpll/issues/129#issuecomment-304004005
         * In short, we use expm1() to compute (exp(Qt) - I), and then correct
         * for this by adding an identity matrix I in the very end */

        /* exponentiate eigenvalues */
        for (j = 0; j < states; ++j)
          expd[j] = expm1(evals[j] * rates[n] * branch_lengths[i]);

        for (j = 0; j < states; ++j)
          for (k = 0; k < states; ++k)
            temp[j*states+k] = inv_evecs[j*states_padded+k] * expd[k];

        for (j = 0; j < states; ++j)
        {
          for (k = 0; k < states; ++k)
          {
            pmat[j*states_padded+k] = (j==k) ? 1.0 : 0;
            for (m = 0; m < states; ++m)
            {
              pmat[j*states_padded+k] +=
                  temp[j*states+m] * evecs[m*states_padded+k];
            }
          }
        }
      }
      #ifdef DEBUG
      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k)
          assert(pmat[j*states_padded+k] >= 0);
      #endif
    }
  }

  free(expd);
  free(temp);
  return BPP_SUCCESS;
}
