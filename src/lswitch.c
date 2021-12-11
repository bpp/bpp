/*
    Copyright (C) 2016-2021 Tomas Flouri and Ziheng Yang

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

/* number of parameters (phi_X,phi_Y,theta_X,theta_Y) in a single BDI */
#define TPARAMS         4

#define IDX_PHI1        0
#define IDX_PHI2        1
#define IDX_THETA1      2
#define IDX_THETA2      3

#define ALG_COG0        0
#define ALG_COGN        1
#define ALG_BG          2
#define ALG_DEFAULT     2

static long algorithm;

static double mean[TPARAMS];
static double var[TPARAMS];
static double hyper[2*TPARAMS];
static double sum_lnphi[2], sum_ln1mphi[2];
static double sum_theta[2], sum_lntheta[2];
static long hparams;

static double lnlike_msci(double x[], int np)
{
   /* loglikelihood for fitting beta(p,q) & gamma(a, b) to phi_x, phi_y, theta_x, theta_y, using
      phi_sln[4] = \sum log phi, for phi_x, phi_y, phi_z, phi_w;
      phi_sln1[4] = sum log (1 - phi);
      theta_s[4] = \sum theta, for theta_x, theta_y, theta_z, theta_w;
      theta_sln[4] = \sum log theta.
      x[] are the 8 or 16 hyperparameters.
   */
   double a, b, p, q, lnp = 0, sumx, sumlnx, sumln1x;

   p = x[0]; q = x[1];  sumlnx = sum_lnphi[0]; sumln1x = sum_ln1mphi[0];  /* phi_x */
   lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
   p = x[2]; q = x[3];  sumlnx = sum_lnphi[1]; sumln1x = sum_ln1mphi[1];  /* phi_y */
   lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
   
   if (opt_est_theta)
   {
     a = x[4]; b = x[5];  sumx = sum_theta[0]; sumlnx = sum_lntheta[0];    /* theta_x */
     lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
     a = x[6]; b = x[7];  sumx = sum_theta[1]; sumlnx = sum_lntheta[1];    /* theta_y */
     lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
   }
   #if 0
   if (_D0DD1) {
      p = x[ 8]; q = x[ 9];  sumlnx = sum_lnphi[2]; sumln1x = sum_ln1mphi[2];  /* phi_z */
      lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
      p = x[10]; q = x[11];  sumlnx = sum_lnphi[3]; sumln1x = sum_ln1mphi[3];  /* phi_w */
      lnp += lgamma(p + q) - lgamma(p) - lgamma(q) + (p - 1) * sumlnx + (q - 1) * sumln1x;
      a = x[12]; b = x[13];  sumx = sum_theta[2]; sumlnx = sum_lntheta[2];    /* theta_z */
      lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
      a = x[14]; b = x[15];  sumx = sum_theta[3]; sumlnx = sum_lntheta[3];    /* theta_w */
      lnp += a * log(b) - lgamma(a) + (a - 1) * sumlnx - b * sumx;
   }
   #endif
   return (-lnp);
}

static void switch_tower(double * sample)
{
  sample[0] = 1 - sample[0];
  sample[1] = 1 - sample[1];

  SWAP(sample[2], sample[3]);
}

static void fit_beta_moments(double pq[2], double m, double v)
{
  double z = m*(1-m)/v - 1;
  if (z <= 0)
  {
    printf("v > m*(1-m) in fit_beta_moments\n");
    z = 0.01;
  }

  pq[0] = m*z; pq[1] = (1-m)*z;
}

static void fit_gamma_moments(double ab[2], double m, double v)
{
  ab[0] = m*m / v;
  ab[1] = m / v;
}

static double ln_beta_ratio(double xnew, double x, double p, double q)
{
  double lnpd = 0;
  if (fabs(xnew - x) > 1e-200)
    lnpd = (p-1)*log(xnew/x) + (q-1)*log((1-xnew)/(1-x));
  return lnpd;
}

static double ln_gamma_ratio(double xnew, double x, double a, double b)
{
  double lnpd = 0;
  if (fabs(xnew - x) > 1e-200)
    lnpd = (a-1)*log(xnew/x) - b*(xnew-x);
  return lnpd;
}

#if 0
static void print_towers(long * tower)
{
  long i;

  for (i = 0; i < opt_samples; ++i)
  {
    printf("%ld", tower[i]);
    if ((i + 1) % 100 == 0) printf(" [%5ld]\n", i+1);
  }
}
#endif

static void update_summary(double ** matrix, long * pindex, long * tower)
{
  long i,j;
  double mold;
  double p[TPARAMS];

  for (i = 0; i < TPARAMS; ++i)
    mean[i] = var[i] = 0;
  for (i = 0; i < TPARAMS/2; ++i)
    sum_lnphi[i] = sum_ln1mphi[i] = sum_theta[i] = sum_lntheta[i] = 0;

  for (i = 0; i < opt_samples; ++i)
  {
    p[0] = matrix[pindex[0]][i];   /* phi x */
    p[1] = matrix[pindex[1]][i];   /* phi y */
    p[2] = (pindex[2] >= 0) ? matrix[pindex[2]][i] : -1;   /* theta x */
    p[3] = (pindex[3] >= 0) ? matrix[pindex[3]][i] : -1;   /* theta y */

    if (tower[i]) switch_tower(p);
    
    /* progressively calculate mean and variance */
    for (j = 0; j < TPARAMS; ++j)
    {
      if (p[j] == -1) continue;

      mold = mean[j];
      mean[j] += (p[j] - mold) / (i + 1.);
      var[j]  += (p[j] - mold) * (p[j] - mean[j]);
    }

    sum_lnphi[0] += log(p[0]);    sum_ln1mphi[0] += log(1 - p[0]);
    sum_lnphi[1] += log(p[1]);    sum_ln1mphi[1] += log(1 - p[1]);

    if (p[2] >= 0)
    {
      sum_theta[0] += p[2];
      sum_lntheta[0] += log(p[2]);
    }
    if (p[3] >= 0)
    {
      sum_theta[1] += p[3];
      sum_lntheta[1] += log(p[3]);
    }
  }

  for (i = 0; i < TPARAMS; ++i)  var[i] /= opt_samples-1;
  for (i = 0; i < TPARAMS/2; ++i)
  {
    sum_lnphi[i] /= opt_samples;   sum_ln1mphi[i] /= opt_samples;

    if (p[2+i] == -1) continue;

    sum_theta[i] /= opt_samples;   sum_lntheta[i] /= opt_samples;
  }
}

static long compare_towers(long * tower, double ** matrix, long * pindex)
{
  long i,j;
  long changes = 0;
  double p[TPARAMS];
  double pnew[TPARAMS];
  double score_diff = 0;

  for (i = 0; i < opt_samples; ++i)
  {
    pnew[0] = p[0] = matrix[pindex[0]][i];   /* phi x */
    pnew[1] = p[1] = matrix[pindex[1]][i];   /* phi y */
    pnew[2] = p[2] = (pindex[2] >= 0) ? matrix[pindex[2]][i] : -1;   /* theta x */
    pnew[3] = p[3] = (pindex[3] >= 0) ? matrix[pindex[3]][i] : -1;   /* theta y */

    switch_tower(p);
    if (!tower[i])
      for (j = 0; j < TPARAMS; ++j)
        SWAP(p[j],pnew[j]);

    /* diff of eucledian distances between mean and each set of parameters */
    if (algorithm == ALG_COG0)
    {
      for (j = 0, score_diff = 0; j < TPARAMS; ++j)
        if (p[j] >= 0)
          score_diff += (p[j] - pnew[j])*(p[j] + pnew[j] - 2*mean[j]);
    }
    else if (algorithm == ALG_COGN)
    {
      for (j = 0, score_diff = 0; j < TPARAMS; ++j)
        if (p[j] >= 0)
          score_diff += (p[j] - pnew[j])*(p[j] + pnew[j] - 2*mean[j]) / (2*var[j]);
    }
    else
    {
      assert(algorithm == ALG_BG);
      score_diff  =  ln_beta_ratio(pnew[0], p[0], hyper[0], hyper[1]);    /* phi_x   */
      score_diff +=  ln_beta_ratio(pnew[1], p[1], hyper[2], hyper[3]);    /* phi_y   */
      if (opt_est_theta)
      {
        score_diff += ln_gamma_ratio(pnew[2], p[2], hyper[4], hyper[5]);    /* theta_x */
        score_diff += ln_gamma_ratio(pnew[3], p[3], hyper[6], hyper[7]);    /* theta_y */
      }
    }

    if (score_diff > 0)
    {
      tower[i] = !tower[i];
      ++changes;
    }
  }
  return changes;
}

static void init_tower(const char * phi1_label,
                       const char * phi2_label,
                       double ** matrix,
                       long * pindex,
                       long * tower)
{
  long i;

  for (i = 0; i < opt_samples; ++i)
  {
    if (matrix[pindex[0]][i] <= 0 || matrix[pindex[0]][i] >= 1)
      fatal("phi_%s out of range", phi1_label);
    if (matrix[pindex[1]][i] <= 0 || matrix[pindex[1]][i] >= 1)
      fatal("phi_%s out of range", phi2_label);

    tower[i] = (matrix[pindex[0]][i] < .5 || matrix[pindex[1]][i] < .5) ? 0 : 1;
  }

  update_summary(matrix, pindex, tower);
  memset(hyper,0,hparams*sizeof(double));

  if (algorithm == ALG_BG)
  {
    fit_beta_moments(hyper+0, mean[0], var[0]);         /* p q for phi_x   */
    fit_beta_moments(hyper+2, mean[1], var[1]);         /* p q for phi_y   */
    if (opt_est_theta)
    {
      fit_gamma_moments(hyper+4, mean[2], var[2]);         /* a b for theta_x */
      fit_gamma_moments(hyper+6, mean[3], var[3]);         /* a b for theta_y */
    }
  }
}

static void update_matrix(double ** matrix, long * pindex, long * tower)
{
  long i,j;
  double p[TPARAMS];

  for (i = 0; i < opt_samples; ++i)
  {
    if (tower[i])
    {
      for (j = 0; j < TPARAMS; ++j)
        if (pindex[j] >= 0)
          p[j] = matrix[pindex[j]][i];
        else
          p[j] = -1;

      switch_tower(p);

      for (j = 0; j < TPARAMS; ++j)
        if (p[j] >= 0)
          matrix[pindex[j]][i] = p[j];
    }
  }
}

static void write_output(FILE * fp_out,
                         const char * filename,
                         const char * header,
                         double ** matrix,
                         long col_count)
{
  long i,j;

  printf("Printing processed sample into %s\n\n", filename);

  fprintf(fp_out, "%s\n", header);
  for (i = 0; i < opt_samples; ++i)
  {
    /* write sample */
    fprintf(fp_out, "%ld", (i+1)*opt_samplefreq);

    for (j = 0; j < col_count; ++j)
    {
      if (j == col_count-1 && opt_usedata)
        fprintf(fp_out, "\t%.3f", matrix[j][i]);
      else
        fprintf(fp_out, "\t%.6f", matrix[j][i]);

    }
    fprintf(fp_out, "\n");
  }
}

void lswitch(stree_t * stree, const char * header, double ** matrix, long col_count)
{
  long i,j;
  long mpoint_count;
  long theta_count = 0;
  long tau_count = 0;
  long * pindex;
  long * tower;
  long rounds = 100;
  FILE * fp_out;

  double bounds[16][2];
  double e = 1e-7;
  double space[10000];
  double lnL = 0;

  char * outfile = NULL;
  xasprintf(&outfile, "%s.processed", opt_mcmcfile);
  fp_out = xopen(outfile, "w");

  assert(opt_method == METHOD_00);
  assert(opt_msci);

  /* set default algorithm */
  algorithm = ALG_DEFAULT;

  /* allocate space for storing the parameter indices (columns) */
  pindex = (long *)xmalloc(TPARAMS*sizeof(long));
  
  /* allocate space for the indicator variable */
  tower = (long *)xcalloc(opt_samples,sizeof(long));

  hparams = opt_est_theta ? 2*TPARAMS : TPARAMS;

  /* find positions */

  /* find theta count */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau)
      ++tau_count;
    if (opt_est_theta && stree->nodes[i]->theta >= 0)
      ++theta_count;
  }

  /* go through all BDI events */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (!(stree->nodes[i]->hybrid &&
          node_is_bidirection(stree->nodes[i]) &&
          stree->nodes[i]->prop_tau))
      continue;

    printf("Resolving potential unidentifiability for BDI %s <-> %s\n",
           stree->nodes[i]->label, stree->nodes[i]->hybrid->parent->label);

    /* phi1 and phi2 indices */
    pindex[0] = theta_count + tau_count +
                stree->nodes[i]->hybrid->node_index -
                stree->tip_count - stree->inner_count;
    pindex[1] = theta_count + tau_count +
                stree->nodes[i]->hybrid->parent->hybrid->node_index -
                stree->tip_count - stree->inner_count;

    /* find the first theta column index in MCMC sample */
    if (opt_est_theta && stree->nodes[i]->theta >= 0)
    {
      long index = 0;
      for (j = 0; j < stree->nodes[i]->node_index; ++j)
        if (stree->nodes[j]->theta >= 0) ++index;

      pindex[2] = index;
    }
    else
    {
      pindex[2] = -1;
    }

    /* find the second theta column index in MCMC sample */
    if (opt_est_theta && stree->nodes[i]->hybrid->parent->theta >= 0)
    {
      long index = 0;
      for (j = 0; j < stree->nodes[i]->hybrid->parent->node_index; ++j)
        if (stree->nodes[j]->theta >= 0) ++index;

      pindex[3] = index;
    }
    else
    {
      pindex[3] = -1;
    }

    #if 0
    printf("Indices: Phi1: %ld Phi2: %ld Theta1: %ld Theta2: %ld\n",
           pindex[0], pindex[1], pindex[2], pindex[3]);
    printf("%f %f", matrix[pindex[0]][0], matrix[pindex[1]][0]);
    if (pindex[2] >= 0)
      printf(" %f", matrix[pindex[2]][0]);
    if (pindex[3] >= 0)
      printf(" %f", matrix[pindex[3]][0]);
    printf("\n");
    #endif

    /* initialize tower array */
    init_tower(stree->nodes[i]->label,
               stree->nodes[i]->hybrid->parent->label,
               matrix,
               pindex,
               tower);

    #if 0
    printf("Means: %f %f %f %f\n", mean[0], mean[1], mean[2], mean[3]);
    printf(" Vars: %f %f %f %f\n", var[0], var[1], var[2], var[3]);
    #endif
    if (algorithm == ALG_BG)
    {
      lnL = lnlike_msci(hyper,hparams);
      printf("lnL0: %f\n", lnL);
    }

    for (j = 0; j < hparams; ++j) { bounds[j][0] = 0.5; bounds[j][1] = 99999; }
    for (j = 0; j < rounds; ++j)
    {
      mpoint_count = compare_towers(tower, matrix, pindex);
      printf("Round %2ld, %2ld points moved...\n", j, mpoint_count);
      update_summary(matrix, pindex, tower);
      if (algorithm == ALG_BG) { double lnL = lnlike_msci(hyper, hparams); printf("  lnL = %f\n", lnL); }
      if (!mpoint_count) break;

      if (algorithm == ALG_BG)
      {

        int k = ming2(NULL, &lnL, lnlike_msci, NULL, hyper, bounds, space, e, hparams);
        for (k = 0; k < 8; ++k)
          printf("  %f", hyper[k]);
        printf("\n");
      }
    }

    mpoint_count = 0;
    for (j = 0; j < opt_samples; ++j)
      mpoint_count += (tower[j] > 0);
    printf("\n%4ld points reflected\n", mpoint_count);
    
    update_matrix(matrix, pindex, tower);
  }

  /* output */
  write_output(fp_out, outfile, header, matrix, col_count);

  free(pindex);
  free(tower);
  free(outfile);
  fclose(fp_out);
}
