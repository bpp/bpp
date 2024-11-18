/*
    Copyright (C) 2016-2024 Tomas Flouri, Xiyun Jiao, Bruce Rannala and Ziheng Yang

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

#define DIST_GAMMA      0
#define DIST_INVGAMMA   1
#define DIST_BETA       2

typedef struct Distribution_s
{
  double y;
  double pdf;
} Distribution;

int compare(const void * pa, const void * pb)
{
  const Distribution * a = (const Distribution *)pa;
  const Distribution * b = (const Distribution *)pb;

  if (a->y < b->y) return -1;
  if (a->y > b->y) return 1;
  
  return 0;
}

static void posterior_interval(int nbins,
                               int kmode,
                               double tail,
                               double * x,
                               double * cdf,
                               double * HPD,
                               double * EqualTail)
{
  /* this is the function to calculate the HPD and equal-tail intervals using estimated cdf */
  int lindex, uindex, l0, u0;
  double pinterval;
  /* calculate HPD interval, current interval is (l0, u0) */
  lindex = 0; uindex = nbins - 1;
  for (l0 = 0, u0 = 1; l0 < kmode + 1 && l0 < nbins - 1; l0++) {
    /* u0 = l0 + 1; */  /* use u0 from previous round */
    for (; u0 < nbins; u0++) {
      if ((pinterval = cdf[u0] - cdf[l0]) > 1 - tail)
        break;
      if (l0 == 0 || u0 - l0 < uindex - lindex) {
        lindex = l0;  uindex = u0;
      }
    }
    if (u0 == nbins) break;
  }
  HPD[0] = x[lindex];
  HPD[1] = x[uindex];

  /* calculate equal-tail intercal */
  for (l0 = 0; l0 < nbins - 1; l0++)
    if (cdf[l0] > tail / 2) break;
  for (u0 = nbins - 1; u0 > 0; u0--)
    if (cdf[u0] < 1 - tail / 2) break;
  EqualTail[0] = x[l0];
  EqualTail[1] = x[u0];
}

void conditional_to_marginal_M(double * ai1,
                               double * bi1,
                               double * ai2,
                               double * bi2,
                               int nsamples,
                               int nbins,
                               int cond_dist1,
                               double tail,
                               double * out_wmean,
                               double * out_wstdev,
                               double * out_mean,
                               double * out_sd,
                               double * out_et025,
                               double * out_et975,
                               double * out_hpd025,
                               double * out_hpd975,
                               double * out_c,
                               double * out_effu,
                               double * out_effy)
{
  /* this uses the conditional gamma or inverse-gamma distributions with parameters
  (ai_1, bi_1) and (ai_2, bi_2) for tweo variables theta and w
  to recover the marginal distribution of the product M = theta * w / 4,
  and calculates the mean, variance, CI and efficiency */
  int i, k, l, j, jmode, npoints;
  double* ui1, * ui2, * ui, * y1, * y2, * lny1, * lny2, * y, * cdf;
  double lbin1, lbin2, ubin1, ubin2, bin_width_1, bin_width_2, a1, b1, a2, b2;
  double mean1, mean2, meany, meany2;
  double variance1, variance2, variancey;
  double sd1, sd2, sdy;
  double m1, m2, v1, v2, mv1, mv2, vm1, vm2, my, my2, vmy;
  double rho1_u;
  double E_u, E_y;
  double T_u, c, lnconst1, lnconst2;
  double HPD[2], EqualTail[2];
  #if 0
  char timestr[32];
  #endif

  #ifdef A1B1_DEBUG
  printf("\n*** conditional_to_marginal_M ***\n");
  #endif
  #if 0
  starttimer();
  #endif
  npoints = nbins * nbins;
  ui1 = (double*)xmalloc((3 * nsamples + 4 * nbins + 2 * npoints)*sizeof(double));
  //if (ui1 == NULL) fatal("oom");
  ui2 = ui1 + nsamples;
  ui = ui2 + nsamples;
  y1 = ui + nsamples;
  y2 = y1 + nbins;
  lny1 = y2 + nbins;
  lny2 = lny1 + nbins;
  y = lny2 + nbins;
  cdf = y + npoints;

  /* calculate mean and variance of theta, w, and M */
  meany2 = 0; mv1 = mv2 = 0;
  for (i = 0; i < nsamples; i++)
  {
    if (ai1[i] <= 0. || ai2[i] <= 0. || bi1[i] <= 0. || bi2[i] <= 0.)
      fatal("record %4d: ai_1 = %9.6f bi_1 = %9.6f ai_2 = %9.6f bi_2 = %9.6f!",
            i + 1, ai1[i], bi1[i], ai2[i], bi2[i]);
    if (cond_dist1 == DIST_GAMMA)
    {
      /* gamma */
      m1 = ai1[i] / bi1[i];
      v1 = ai1[i] / (bi1[i] * bi1[i]);
    }
    else if (cond_dist1 == DIST_INVGAMMA)
    {
      /* inverse gamma */
      if (ai1[i] <= 2.)
        fatal("ai1<=2!");
      m1 = bi1[i] / (ai1[i] - 1.);
      v1 = (bi1[i] * bi1[i]) / ((ai1[i] - 1.) * (ai1[i] - 1.) * (ai1[i] - 2.));
    }
    else
      assert(0);

    m2 = ai2[i] / bi2[i];
    v2 = ai2[i] / (bi2[i] * bi2[i]);

    my = m1 * m2;
    my2 = (v1 + m1 * m1) * (v2 + m2 * m2);

    ui1[i] = m1; ui2[i] = m2;
    ui[i] = my;

    mv1 += (v1 - mv1) / (i + 1.);                   /* mean of the conditional variance for theta */
    mv2 += (v2 - mv2) / (i + 1.);                   /* mean of the conditional variance for w */

    meany2 += (my2 - meany2) / (i + 1.);  /* mean of the conditional mean of M^2 */
  }

  /* calculate means and stdevs */
  mean1 = mean2 = meany = 0;
  for (i = 0; i < nsamples; ++i)
  {
    mean1 += ui1[i];
    mean2 += ui2[i];
    meany += ui[i];
  }
  mean1 /= nsamples;
  mean2 /= nsamples;
  meany /= nsamples;

  sd1 = sd2 = sdy = 0;
  for (i = 0; i < nsamples; ++i)
  {
    sd1 += (ui1[i]-mean1)*(ui1[i]-mean1);
    sd2 += (ui2[i]-mean2)*(ui2[i]-mean2);
    sdy += (ui[i]-meany)*(ui[i]-meany);
  }
  sd1 = sqrt(sd1/nsamples);
  sd2 = sqrt(sd2/nsamples);
  sdy = sqrt(sdy/nsamples);

  #if 0
  double rho1_u_1, rho1_u_2;
  double E_u_1, E_u_2;
  E_u_1 = 1/eff_ict(ui1, nsamples, mean1, sd1, &rho1_u_1);
  E_u_2 = 1/eff_ict(ui2, nsamples, mean2, sd2, &rho1_u_2);
  #endif
  T_u   = eff_ict(ui,  nsamples, meany, sdy, &rho1_u);
  vm1 = sd1*sd1;   vm2 = sd2*sd2;   vmy = sdy*sdy;

  #ifdef A1B1_DEBUG
  printf("\ntheta: mean&sd = %9.6f\t%9.6f\n", mean1, sqrt(vm1+mv1));
  printf("w: mean&sd = %9.6f\t%9.6f\n", mean2, sqrt(vm2+mv2));
  #endif
  *out_wmean = mean2;
  *out_wstdev = sqrt(vm2+mv2);

  mean1 /= 2.; mean2 /= 2.;
  variance1 = (vm1 + mv1) / 4.; variance2 = (vm2 + mv2) / 4.;
  sd1 = sqrt(variance1); sd2 = sqrt(variance2);

  variancey = meany2 - meany * meany;
  meany /= 4.; vmy /= 16.; variancey /= 16.;
  sdy = sqrt(variancey);

  /* calculate the efficiency of M using the efficiency of the corresponding conditional mean sample */
  c = variancey / vmy; 
  E_u = c / T_u;
  E_y = 1. / (1. + (T_u - 1.) / c);

  #ifdef A1B1_DEBUG
  printf("\nM: mean&sd = %9.6f\t%9.6f\n", meany, sdy);
  printf("M: c = %9.6f\t efficiency of u = %9.6f\t efficiency of y = %9.6f\n", c, E_u, E_y);
  #endif
  #if 0
  printf("Processed samples, time used: %s\n", printtime(timestr));
  #endif

  /* estimate marginal PDF or histogram of M = theta * w / 4;
     the points are defined on a rectangle area,
     and the pdf for each point is the average of the pdf for theta/2 times the pdf for w/2 */
  Distribution* pdfy = xmalloc(npoints*sizeof(Distribution));
  //if (pdfy == NULL) fatal("oom");

  lbin1 = mean1 - 4. * sd1; ubin1 = mean1 + 4. * sd1;
  lbin1 = (lbin1 > 0.) ? lbin1 : 0.;
  bin_width_1 = (ubin1 - lbin1) / nbins;
  lbin2 = mean2 - 4. * sd2; ubin2 = mean2 + 4. * sd2;
  lbin2 = (lbin2 > 0.) ? lbin2 : 0.;
  bin_width_2 = (ubin2 - lbin2) / nbins;
  memset(y1, 0, nbins * sizeof(double));
  memset(y2, 0, nbins * sizeof(double));
  for (k = 0; k < nbins; k++) {
    y1[k] = lbin1 + (k + 0.5) * bin_width_1;
    y2[k] = lbin2 + (k + 0.5) * bin_width_2;
    lny1[k] = log(y1[k]);
    lny2[k] = log(y2[k]);
  }

  for (j = 0; j < npoints; j++)
    pdfy[j].pdf = 0;

  for (k = 0; k < nbins; k++)
    for (l = 0; l < nbins; l++)
      pdfy[k * nbins + l].y = y1[k] * y2[l];

  for (i = 0; i < nsamples; i++) {
    a1 = ai1[i]; b1 = bi1[i];
    a2 = ai2[i]; b2 = bi2[i];
    if (cond_dist1 == DIST_GAMMA)
      lnconst1 = a1 * (log(2.) + log(b1)) - lgamma(a1);
    else
      lnconst1 = a1 * (log(b1) - log(2.)) - lgamma(a1);

    lnconst2 = a2 * (log(2.) + log(b2)) - lgamma(a2);

    if (cond_dist1 == DIST_GAMMA)
      for (k = 0; k < nbins; k++)
        for (l = 0; l < nbins; l++)
          pdfy[k * nbins + l].pdf += exp(lnconst1 + (a1 - 1.) * lny1[k] - 2. * b1 * y1[k] + lnconst2 + (a2 - 1.) * lny2[l] - 2. * b2 * y2[l]);
    else
      for (k = 0; k < nbins; k++)
        for (l = 0; l < nbins; l++)
          pdfy[k * nbins + l].pdf += exp(lnconst1 + (-a1 - 1.) * lny1[k] - b1 / (2. * y1[k]) + lnconst2 + (a2 - 1.) * lny2[l] - 2. * b2 * y2[l]);
  }

  qsort(pdfy, npoints, sizeof(Distribution), compare);

  for (j = 0; j < npoints; j++)
  {
    pdfy[j].pdf *= (bin_width_1 * bin_width_2) / nsamples;
    y[j] = pdfy[j].y;
  }

  for (j = 1, jmode = 0; j < npoints; j++)
    if (pdfy[j].pdf > pdfy[jmode].pdf) jmode = j;  /* mode */
  for (j = 1, cdf[0] = pdfy[0].pdf; j < npoints; j++)
    cdf[j] = pdfy[j].pdf + cdf[j - 1];           /* estimated cdf for M */

  #ifdef A1B1_DEBUG
  printf("\nM: sum pdf = %9.6f\n", cdf[npoints - 1]);
  #endif

  /* calculate the HPD and equal-tail intervals for M = theta * w / 4 */
  posterior_interval(npoints, jmode, tail, y, cdf, HPD, EqualTail);

  #ifdef A1B1_DEBUG
  printf("M: HPD: (%8.6f, %8.6f) equal-tail: (%8.6f, %8.6f)\n", HPD[0], HPD[1], EqualTail[0], EqualTail[1]);
  #endif
  # if 0
  printf("Calculated pdf over grid, time used: %s\n", printtime(timestr));
  #endif

  
  *out_c = c;
  *out_mean = meany;
  *out_sd = sdy;
  *out_et025 = EqualTail[0];
  *out_et975 = EqualTail[1];
  *out_hpd025 = HPD[0];
  *out_hpd975 = HPD[1];
  *out_effu = E_u;
  *out_effy = E_y;

  free(ui1);
  free(pdfy);
}

void conditional_to_marginal(double * ai,
                             double * bi,
                             long nsamples,
                             long nbins,
                             long cond_dist,
                             double tail,
                             double * out_mean,
                             double * out_sd,
                             double * out_et025,
                             double * out_et975,
                             double * out_hpd025,
                             double * out_hpd975,
                             double * out_c,
                             double * out_effu,
                             double * out_effy)
{
  /* this uses the conditional gamma or inverse-gamma distributions with parameters
  (ai, bi) to recover the marginal distribution, and calculates the mean, variance, CI and efficiency */
  long i, k, kmode;
  double* ui, * pdf, * cdf, * y, * lny;
  double lbin, ubin, bin_width, a, b;
  double EqualTail[2], HPD[2];
  double lnconst, m, v, mv, vm, rho1_u, E_u, E_y, T_u, c;
  double mean,sd,variance;
  #if 0
  char timestr[32];
  #endif

  #ifdef A1B1_DEBUG
  printf("\n*** conditional_to_marginal ***\n");
  #endif
  #if 0
  starttimer();
  #endif
  ui = (double*)xmalloc((nsamples + 4 * nbins) * sizeof(double));
  pdf = ui + nsamples;
  cdf = pdf + nbins;
  y = cdf + nbins;
  lny = y + nbins;

  assert(cond_dist == DIST_GAMMA ||
         cond_dist == DIST_INVGAMMA ||
         cond_dist == DIST_BETA);

  /* calculate the conditional mean and variance using ai and bi and then estimate the mean of variance */
  mv = 0;
  for (i = 0; i < nsamples; i++)
  {
    if (ai[i] <= 0. || bi[i] <= 0.)
      fatal("Error: a1 and b1 must be > 0\n"
            "  line %ld: ai = %9.6f bi = %9.6f", i+1, ai[i], bi[i]);

    if (cond_dist == DIST_GAMMA)
    {
      m = ai[i] / bi[i];
      v = ai[i] / (bi[i] * bi[i]);
    }
    else if (cond_dist == DIST_INVGAMMA)
    {
      if (ai[i] <= 2)
        fatal("Error: a1 must be > 2\n"
              "  line %ld: ai = %9.6f", i+1, ai[i]);
      m = bi[i] / (ai[i] - 1);
      v = bi[i] * bi[i] / ((ai[i] - 1) * (ai[i] - 1) * (ai[i] - 2));
    }
    else
    {
      /* beta */
      m = ai[i] / (ai[i] + bi[i]);
      v = ai[i]*bi[i] / ((ai[i]+bi[i]) * (ai[i]+bi[i]) * (ai[i]+bi[i]+1.));
    }
    ui[i] = m;

    mv += (v - mv) / (i + 1);                   /* mean of the conditional variance */
  }

  /* calculate the efficiency using the efficiency of the conditional mean sample */
  //E_u = Eff_IntegratedCorrelationTime(ui, nsamples, &mean, &vm, &rho1_u);

  /* calculate ui mean and stdev */
  mean = 0; sd = 0;
  for (i = 0; i < nsamples; ++i)
    mean += ui[i];
  mean /= nsamples;

  for (i = 0; i < nsamples; ++i)
    sd += (ui[i]-mean)*(ui[i]-mean);
  sd = sqrt(sd/nsamples);

  T_u = eff_ict(ui, nsamples, mean, sd, &rho1_u);
  vm = sd*sd;
  variance = mv + vm; sd = sqrt(variance);
  c = variance / vm; 
  E_u = c / T_u;
  E_y = 1. / (1. + (T_u - 1.) / c);

  /* estimate marginal PDF or histogram */
  lbin = mean - 4*sd; ubin = mean + 4*sd;
  lbin = (lbin > 0) ? lbin : 0;
  if (cond_dist == DIST_BETA)
    ubin = (ubin < 1) ? ubin : 1;
  bin_width = (ubin - lbin) / nbins;
  memset(y,   0, nbins*sizeof(double));
  memset(lny, 0, nbins*sizeof(double));
  for (k = 0; k < nbins; k++)
  {
    y[k] = lbin + (k + 0.5) * bin_width;
    lny[k] = log(y[k]);
  }

  memset(pdf, 0, nbins * sizeof(double));
  for (i = 0; i < nsamples; i++) {
    a = ai[i]; b = bi[i];
    /* constant for all the bins */
    if (cond_dist == DIST_BETA)
      lnconst = lgamma(a + b) - lgamma(a) - lgamma(b);
    else
      lnconst = a * log(b) - lgamma(a);

    if (cond_dist == DIST_GAMMA)
      for (k = 0; k < nbins; k++)
        pdf[k] += exp(lnconst + (a - 1)*lny[k] - b*y[k]);
    else if (cond_dist == DIST_INVGAMMA)
      for (k = 0; k < nbins; k++)
        pdf[k] += exp(lnconst + (-a - 1)*lny[k] - b/y[k]);
    else
      for (k = 0; k < nbins; k++)
        pdf[k] += exp(lnconst + (a - 1.) * lny[k] + (b - 1.) * log(1. - y[k]));
  }

  for (k = 0; k < nbins; k++)
    pdf[k] *= bin_width / nsamples;

  for (k = 1, kmode = 0; k < nbins; k++)
    if (pdf[k] > pdf[kmode]) kmode = k;  /* mode */
  for (k = 1, cdf[0] = pdf[0]; k < nbins; k++)
    cdf[k] = pdf[k] + cdf[k - 1];

  #if 0
  printf("\nsum pdf = %9.6f, pdf/histogram, time used: %s\n", cdf[nbins - 1], printtime(timestr));
  #endif
  
  #ifdef A1B1_DEBUG
  printf("\nsum pdf = %9.6f, pdf/histogram\n", cdf[nbins - 1]);
  #endif

  /* calculate HPD and equal-tail intervals */
  posterior_interval(nbins, kmode, tail, y, cdf, HPD, EqualTail);

  #ifdef A1B1_DEBUG
  printf("mean&sd = %9.6f\t%9.6f\t HPD: (%8.6f, %8.6f) equal-tail: (%8.6f, %8.6f)\n", mean, sd, HPD[0], HPD[1], EqualTail[0], EqualTail[1]);
  printf("c = %9.6f\t efficiency of u = %9.6f\t efficiency of y = %9.6f\n", c, E_u, E_y);
  #endif
  #if 0
  printf("Time used: %s\n", printtime(timestr));
  #endif

  *out_mean = mean;
  *out_sd = sd;
  *out_hpd025 = HPD[0];
  *out_hpd975 = HPD[1];
  *out_et025 = EqualTail[0];
  *out_et975 = EqualTail[1];
  *out_effu = E_u;
  *out_effy = E_y;
  *out_c = c;

  free(ui);
}
