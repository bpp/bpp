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

static int mytqli(double *d, double *e, const unsigned int n, double **z)
{
  unsigned int     m, l, iter, i, k;
  double  s, r, p, g, f, dd, c, b;

  for (i = 2; i <= n; i++)
    e[i - 2] = e[i - 1];

  e[n - 1] = 0.0;

  for (l = 1; l <= n; l++)
    {
      iter = 0;
      do
        {
          for (m = l; m <= n - 1; m++)
            {
              dd = fabs(d[m - 1]) + fabs(d[m]);
              if (fabs(e[m - 1]) + dd == dd)
                break;
            }
          if (m != l)
           {
             assert(iter < 30);

             g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);
             r = sqrt((g * g) + 1.0);
             g = d[m - 1] - d[l - 1] + e[l - 1] / (g + ((g < 0)?-fabs(r):fabs(r)));/*MYSIGN(r, g));*/
             s = c = 1.0;
             p = 0.0;

             for (i = m - 1; i >= l; i--)
               {
                 f = s * e[i - 1];
                 b = c * e[i - 1];
                 if (fabs(f) >= fabs(g))
                   {
                     c = g / f;
                     r = sqrt((c * c) + 1.0);
                     e[i] = f * r;
                     c *= (s = 1.0 / r);
                   }
                 else
                   {
                     s = f / g;
                     r = sqrt((s * s) + 1.0);
                     e[i] = g * r;
                     s *= (c = 1.0 / r);
                   }
                 g = d[i] - p;
                 r = (d[i - 1] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++)
                   {
                     f = z[i][k-1];
                     z[i][k-1] = s * z[i - 1][k - 1] + c * f;
                     z[i - 1][k - 1] = c * z[i - 1][k - 1] - s * f;
                   }
               }

             d[l - 1] = d[l - 1] - p;
             e[l - 1] = g;
             e[m - 1] = 0.0;
           }
        }
      while (m != l);
    }



    return (1);
 }


static void mytred2(double **a, const unsigned int n, double *d, double *e)
{
  unsigned int     l, k, j, i;
  double  scale, hh, h, g, f;

  for (i = n; i > 1; i--)
    {
      l = i - 1;
      h = 0.0;
      scale = 0.0;

      if (l > 1)
        {
          for (k = 1; k <= l; k++)
            scale += fabs(a[k - 1][i - 1]);
          if (scale == 0.0)
            e[i - 1] = a[l - 1][i - 1];
          else
            {
              for (k = 1; k <= l; k++)
                {
                  a[k - 1][i - 1] /= scale;
                  h += a[k - 1][i - 1] * a[k - 1][i - 1];
                }
              f = a[l - 1][i - 1];
              g = ((f > 0) ? -sqrt(h) : sqrt(h)); /* diff */
              e[i - 1] = scale * g;
              h -= f * g;
              a[l - 1][i - 1] = f - g;
              f = 0.0;
              for (j = 1; j <= l; j++)
                {
                  a[i - 1][j - 1] = a[j - 1][i - 1] / h;
                  g = 0.0;
                  for (k = 1; k <= j; k++)
                    g += a[k - 1][j - 1] * a[k - 1][i - 1];
                  for (k = j + 1; k <= l; k++)
                    g += a[j - 1][k - 1] * a[k - 1][i - 1];
                  e[j - 1] = g / h;
                  f += e[j - 1] * a[j - 1][i - 1];
                }
              hh = f / (h + h);
              for (j = 1; j <= l; j++)
                {
                  f = a[j - 1][i - 1];
                  g = e[j - 1] - hh * f;
                  e[j - 1] = g;
                  for (k = 1; k <= j; k++)
                    a[k - 1][j - 1] -= (f * e[k - 1] + g * a[k - 1][i - 1]);
                }
            }
        }
      else
        e[i - 1] = a[l - 1][i - 1];
      d[i - 1] = h;
    }
  d[0] = 0.0;
  e[0] = 0.0;

  for (i = 1; i <= n; i++)
    {
      l = i - 1;
      if (d[i - 1] != 0.0)
        {
          for (j = 1; j <= l; j++)
            {
                g = 0.0;
                for (k = 1; k <= l; k++)
                  g += a[k - 1][i - 1] * a[j - 1][k - 1];
                for(k = 1; k <= l; k++)
                  a[j - 1][k - 1] -= g * a[i - 1][k - 1];
            }
        }
      d[i - 1] = a[i - 1][i - 1];
      a[i - 1][i - 1] = 1.0;
      for (j = 1; j <= l; j++)
        a[i - 1][j - 1] = a[j - 1][i - 1] = 0.0;
    }
}

/* TODO: Add code for SSE/AVX. Perhaps allocate qmatrix in one chunk to avoid the
complex checking when to dealloc */
static double ** create_ratematrix(double * params,
                                   double * frequencies,
                                   long states)
{
  long i,j,k;
  double ** qmatrix;

  /* normalize substitution parameters */
  long params_count = (states*states - states) / 2;
  double * params_normalized = (double *)xmalloc(sizeof(double) * params_count);

  memcpy(params_normalized,params,params_count*sizeof(double));

  if (params_normalized[params_count - 1] > 0.0)
    for (i = 0; i < params_count; ++i)
      params_normalized[i] /= params_normalized[params_count - 1];

  /* allocate qmatrix */
  qmatrix = (double **)xmalloc(states*sizeof(double *));
  for (i = 0; i < states; ++i)
    qmatrix[i] = (double *)xmalloc(states*sizeof(double));

  /* construct a matrix equal to sqrt(pi) * Q sqrt(pi)^-1 in order to ensure
     it is symmetric */

  for (i = 0; i < states; ++i)
    qmatrix[i][i] = 0;

  k = 0;
  for (i = 0; i < states; ++i)
  {
    for (j = i+1; j < states; ++j)
    {
      double factor = params_normalized[k++];
      qmatrix[i][j] = qmatrix[j][i] = factor * sqrt(frequencies[i] * frequencies[j]);
      qmatrix[i][i] -= factor * frequencies[j];
      qmatrix[j][j] -= factor * frequencies[i];
    }
  }


  double mean = 0;
  for (i = 0; i < states; ++i)
    mean += frequencies[i] * (-qmatrix[i][i]);
  for (i = 0; i < states; ++i)
    for (j = 0; j < states; ++j)
      qmatrix[i][j] /= mean;

  free(params_normalized);

  return qmatrix;
}

void pll_update_eigen(double * eigenvecs,
                      double * inv_eigenvecs,
                      double * eigenvals,
                      double * freqs,
                      double * subst_params,
                      long states,
                      long states_padded)
{
  long i,j,k;
  double *e, *d;
  double ** a;

  a = create_ratematrix(subst_params,
                        freqs,
                        states);
  if (!a)
    fatal("Unable to allocate enough memory.");

  d = (double *)xmalloc(states*sizeof(double));
  e = (double *)xmalloc(states*sizeof(double));

  mytred2(a, states, d, e);
  mytqli(d, e, states, a);

  /* store eigen vectors */
  for (i = 0; i < states; ++i)
    memcpy(eigenvecs + i*states_padded, a[i], states*sizeof(double));

  /* store eigen values */
  memcpy(eigenvals, d, states*sizeof(double));

  /* store inverse eigen vectors */
  for (k = 0, i = 0; i < states; ++i)
  {
    for (j = i; j < states_padded*states; j += states_padded)
      inv_eigenvecs[k++] = eigenvecs[j];

    /* account for padding */
    k += states_padded - states;
  }

  assert(k == states_padded*states);

  /* multiply the inverse eigen vectors from the left with sqrt(pi)^-1 */
  for (i = 0; i < states; ++i)
    for (j = 0; j < states; ++j)
      inv_eigenvecs[i*states_padded+j] /= sqrt(freqs[i]);

  /* multiply the eigen vectors from the right with sqrt(pi) */
  for (i = 0; i < states; ++i)
    for (j = 0; j < states; ++j)
      eigenvecs[i*states_padded+j] *= sqrt(freqs[j]);

  free(d);
  free(e);
  for (i = 0; i < states; ++i)
    free(a[i]);
  free(a);
}

int pll_core_update_pmatrix_4x4_f81(double ** pmatrix,
                                    unsigned int states,
                                    unsigned int rate_cats,
                                    const double * rates,
                                    const double * branch_lengths,
                                    const double ** frequencies,
                                    const unsigned int * matrix_indices,
                                    const unsigned int * param_indices,
                                    unsigned int count,
                                    unsigned int attrib)
{
  unsigned int i,j,k,n,m;
  unsigned int states_padded = states;
  double t;

  double * pmat;
  const double * freqs;
  double beta;

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;
      freqs = frequencies[param_indices[n]];
      t = branch_lengths[i]*rates[n];

      /* compute beta */
      for (beta=1,j = 0; j < 4; ++j)
        beta -= freqs[j]*freqs[j];
      beta = 1./beta;
      
      double e = exp(-beta*t);
      double em1 = expm1(-beta*t);

      /* fill pmatrix */
      for (m=0,j = 0; j < 4; ++j)
        for (k = 0; k < 4; ++k)
          if (j==k)
            pmat[m++]  = e - freqs[k]*em1;
          else
            pmat[m++]  = -freqs[k]*em1;
    }
  }
  return BPP_SUCCESS;
}

int pll_core_update_pmatrix_4x4_k80(double ** pmatrix,
                                    unsigned int states,
                                    unsigned int rate_cats,
                                    double kappa,
                                    const double * rates,
                                    const double * branch_lengths,
                                    const unsigned int * matrix_indices,
                                    const unsigned int * param_indices,
                                    unsigned int count,
                                    unsigned int attrib)
{
  unsigned int i,j,k,n,m;
  unsigned int states_padded = states;
  double t;
  double e1,e2;

  double * pmat;

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;
      t = branch_lengths[i]*rates[n];
      e1 = expm1(-4*t / (kappa+2));

      if (fabs(kappa-1) < 1e-20)
      {
        for (m=0, j = 0; j < 4; ++j)
          for (k = 0; k < 4; ++k)
            if (j == k)
              pmat[m++] = 1. + 3/4.*e1;
            else
              pmat[m++] = -e1/4;
      }
      else
      {
        e2 = expm1(-2 * t*(kappa+1)/(kappa+2));

        pmat[0]  = 1 + (e1 + 2*e2)/4;       /* AA */
        pmat[1]  = -e1/4;                   /* AC */
        pmat[2]  = (e1 - 2*e2)/4;           /* AG */
        pmat[3]  = -e1/4;                   /* AT */

        pmat[4]  = -e1/4;                   /* CA */ 
        pmat[5]  = 1 + (e1 + 2*e2)/4;       /* CC */ 
        pmat[6]  = -e1/4;                   /* CG */ 
        pmat[7]  = (e1 - 2*e2)/4;           /* CT */ 
        
        pmat[8]  = (e1 - 2*e2)/4;           /* GA */
        pmat[9]  = -e1/4;                   /* GC */
        pmat[10] = 1 + (e1 + 2*e2)/4;       /* GG */
        pmat[11] = -e1/4;                   /* GT */

        pmat[12] = -e1/4;                   /* TA */
        pmat[13] = (e1 - 2*e2)/4;           /* TC */
        pmat[14] = -e1/4;                   /* TG */
        pmat[15] = 1 + (e1 + 2*e2)/4;       /* TT */
      }
    }
  }
  return BPP_SUCCESS;
}

int pll_core_update_pmatrix_4x4_tn93(double ** pmatrix,
                                     unsigned int states,
                                     unsigned int rate_cats,
                                     double a1t,
                                     double a2t,
                                     const double * rates,
                                     const double * branch_lengths,
                                     const double ** frequencies,
                                     const unsigned int * matrix_indices,
                                     const unsigned int * param_indices,
                                     unsigned int count,
                                     unsigned int attrib)
{
  unsigned int i,n;
  unsigned int states_padded = states;
  double beta;
  double bt;
  double e1,e2,e3;

  double * pmat;
  const double * freqs;


  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;
      freqs = frequencies[param_indices[n]];

      double A = freqs[0];
      double C = freqs[1];
      double G = freqs[2];
      double T = freqs[3];
      double Y = T + C;
      double R = A + G;

      beta = 0.5 / (Y*R + a1t*T*C + a2t*A*G);
      bt = branch_lengths[i]*rates[n]*beta;

      e1 = expm1(-bt);
      e2 = expm1(-(R*a2t + Y*bt));
      e3 = expm1(-(Y*a1t + R*bt));

      pmat[0]  = 1 + Y*A / R*e1 + G / R*e2;
      pmat[1]  = -C*e1;
      pmat[2]  = Y*G / R*e1 - G / R*e2;
      pmat[3]  = -T*e1;

      pmat[4]  = -A*e1;
      pmat[5]  = 1 + (R*C*e1 + T*e3) / Y;
      pmat[6]  = -G*e1;
      pmat[7]  = (R*e1 - e3)*T / Y;

      pmat[8]  = Y*A / R*e1 - A / R*e2;
      pmat[9]  = -C*e1;
      pmat[10] = 1 + Y*G / R*e1 + A / R*e2;
      pmat[11] = -T*e1;

      pmat[12] = -A*e1;
      pmat[13] = (R*e1 - e3)*C / Y;
      pmat[14] = -G*e1;
      pmat[15] = 1 + (R*T*e1 + C*e3) / Y;
    }
  }

  return BPP_SUCCESS;
}
int pll_core_update_pmatrix_4x4_hky(double ** pmatrix,
                                     unsigned int states,
                                     unsigned int rate_cats,
                                     const double * rates,
                                     const double * kappa,
                                     const double * branch_lengths,
                                     const double ** frequencies,
                                     const unsigned int * matrix_indices,
                                     const unsigned int * param_indices,
                                     unsigned int count,
                                     unsigned int attrib)
{
  unsigned int i,j,k,n,m;
  unsigned int states_padded = states;
  double u,w,x,y,z;
  double beta;
  double pijk[4];

  double * pmat;
  const double * freqs;


  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;
      freqs = frequencies[param_indices[n]];


      double t = branch_lengths[i];
      if (rate_cats > 1)
        t *= rates[n];

      double pi_ag_sum  = freqs[0] + freqs[2];
      double pi_ct_sum  = freqs[1] + freqs[3];
      double pi_ag_prod = freqs[0] * freqs[2];
      double pi_ct_prod = freqs[1] * freqs[3];

      pijk[0] = freqs[0] + freqs[2];
      pijk[1] = freqs[1] + freqs[3];
      pijk[2] = freqs[0] + freqs[2];
      pijk[3] = freqs[1] + freqs[3];

      beta = 0.5 / (pi_ag_sum*pi_ct_sum + kappa[n]*(pi_ag_prod + pi_ct_prod));

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
        m=0;
        for (j = 0; j < 4; ++j)
        {
          for (k = 0; k < 4; ++j)
          {
            u = 1.0/pijk[k] - 1;
            w = -beta * (1+pijk[k]*(kappa[n]-1));
            x = exp(-beta*t);
            y = exp(w*t);
            z = (pijk[k] - freqs[k]) / pijk[k];

            if (j == k)
              pmat[m++] = freqs[k] + freqs[k]*u*x + z*y;
            else if ((j == 0 && k == 2) || (j == 2 && k == 0) ||
                     (j == 1 && k == 3) || (j == 3 && k == 1))
              pmat[m++] = freqs[k] + freqs[k]*u*x - (freqs[k]/pijk[k])*y;
            else
              pmat[m++] = freqs[k]*(1-x);
          }
        }
      }
    }
  }

  return BPP_SUCCESS;
}

int pll_core_update_pmatrix_4x4_jc69(double ** pmatrix,
                                     unsigned int states,
                                     unsigned int rate_cats,
                                     const double * rates,
                                     const double * branch_lengths,
                                     const unsigned int * matrix_indices,
                                     const unsigned int * param_indices,
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
        #if 0
        double a =  (1 + 3*exp(-4*t/3) ) / 4;
        double b = (1 - a) / 3;
        #endif

        double exptm1 = expm1(-4*t/3);
        double a = 1 + 3/4. * exptm1;
        double b = -exptm1/4;

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

void bpp_core_update_pmatrix(locus_t * locus,
                             gnode_t ** traversal,
                             unsigned int count)
{
  unsigned int i,n,j,k,m;
  unsigned int states = locus->states;
  unsigned int states_padded = states;
  unsigned int rate_cats = locus->rate_cats;
  unsigned int attrib = locus->attributes;
  double t;
  const double * rates = locus->rates;
  double * const * eigenvals = locus->eigenvals;
  double * const * eigenvecs = locus->eigenvecs;
  double * const * inv_eigenvecs = locus->inv_eigenvecs;
  double * expd;
  double * temp;

  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;

  gnode_t * node;


  expd = (double *)xmalloc(states * sizeof(double));
  temp = (double *)xmalloc(states*states*sizeof(double));

  unsigned int * param_indices = locus->param_indices;


  for (i = 0; i < count; ++i)
  {
    node = traversal[i];

    t = node->length = (node->parent->time - node->time)*locus->mut_rates[0];

    assert(t >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = locus->pmatrix[node->pmatrix_index] + n*states*states_padded;

      evecs = eigenvecs[param_indices[n]];
      inv_evecs = inv_eigenvecs[param_indices[n]];
      evals = eigenvals[param_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (t < 1e-100)
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
          expd[j] = expm1(evals[j] * rates[n] * t);

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
}

int pll_core_update_pmatrix(double ** pmatrix,
                            unsigned int states,
                            unsigned int rate_cats,
                            const double * rates,
                            const double * branch_lengths,
                            const unsigned int * matrix_indices,
                            const unsigned int * param_indices,
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

      evecs = eigenvecs[param_indices[n]];
      inv_evecs = inv_eigenvecs[param_indices[n]];
      evals = eigenvals[param_indices[n]];

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
