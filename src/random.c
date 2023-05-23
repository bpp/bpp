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

#define mBactrian  0.90
#define sBactrian  sqrt(1-mBactrian*mBactrian)
#define aBox 0.5
#define bBox (sqrt(12 - 3*aBox*aBox) - aBox) / 2
#define aAirplane 1.0
#define aStrawhat 1.0

/* legacy random number generators */
static unsigned int * z_rndu = NULL;

void legacy_init(void)
{
   int seed = (int)opt_seed;
   long i;

   /* z_rndu = (unsigned int)opt_seed; */
   if (sizeof(int) != 4)
      fatal("oh-oh, we are in trouble.  int not 32-bit?  rndu() assumes 32-bit int.");

   if (seed <= 0) {
      FILE *frand = fopen("/dev/urandom", "r");
      if (frand) {
         if (fread(&seed, sizeof(int), 1, frand) != 1)
            fatal("failure to read white noise...");
         fclose(frand);
         seed = abs(seed * 2 - 1);
      }
      else {
         seed = abs(1234 * (int)time(NULL) + 1);
      }

      FILE *fseed;
      fseed = fopen("SeedUsed", "w");
      if (fseed == NULL) fatal("can't open file SeedUsed.");
      fprintf(fseed, "%d\n", seed);
      fclose(fseed);
   }

   assert(opt_threads >= 1);

   if (z_rndu) free(z_rndu);
   z_rndu = (unsigned int *)xmalloc((size_t)opt_threads * sizeof(unsigned int));
   for (i = 0; i < opt_threads; ++i)
     z_rndu[i] = (unsigned int)seed;
}

void legacy_fini(void)
{
  free(z_rndu);
}

unsigned int get_legacy_rndu_status(long index)
{
  return z_rndu[index];
}

unsigned int * get_legacy_rndu_array(void)
{
  return z_rndu;
}

void set_legacy_rndu_status(long index, unsigned int x)
{
  z_rndu[index] = x;
}

void set_legacy_rndu_array(unsigned int * x)
{
  if (z_rndu)
    free(z_rndu);
  z_rndu = x;
}

double legacy_rndu(long index)
{
/* 32-bit integer assumed.
   From Ripley (1987) p. 46 or table 2.4 line 2. 
   This may return 0 or 1, which can be a problem.
*/

   /* the below random number generator is the one used until v4.0.6.
      Change if 0 to if 1 to use it */
   #if 0
   z_rndu[index] = z_rndu[index]*69069 + 1;
   if(z_rndu[index] == 0 || z_rndu[index] == 4294967295)  z_rndu[index] = 13;
   return z_rndu[index]/4294967295.0;
   #else
   z_rndu[index] = z_rndu[index] * 69069 + 1;
   if (z_rndu[index] == 0)  z_rndu[index] = 12345671;
   return ldexp((double)(z_rndu[index]), -32);
   #endif
}

static double rndTriangle(long index)
{
	double u, z;
/* Standard Triangle variate, generated using inverse CDF  */
	u = legacy_rndu(index);
	if(u > 0.5)
		z =  sqrt(6.0) - 2.0*sqrt(3.0*(1.0 - u));
   else
		z = -sqrt(6.0) + 2.0*sqrt(3.0*u);
	return z;
}

static double getRoot(double(*f)(double), double(*df)(double), double initVal)
{
   double x, newx = initVal;
   int nIter = 0;
   do {
      x = newx;
      newx = x - (*f)(x) / (*df)(x);
      nIter++;
   } while ((fabs(x - newx) > 1e-10) && nIter < 100);

   if (fabs(x - newx) > 1e-10) {
      fatal("root finder didn't converge");
   }
   return(newx);
}

static double BStrawhat(double b) {
   return 5 * b*b*b - 15 * b + 10 * aStrawhat - 2 * aStrawhat*aStrawhat*aStrawhat;
}

static double dBStrawhat(double b) {
   return 15 * b*b - 15;
}

static double rndStrawhat(long index)
{
   static int firsttime = 1;
   static double bStrawhat;
   double z;

   if (firsttime) {
      bStrawhat = getRoot(&BStrawhat, &dBStrawhat, 2.0);
      firsttime = 0;
   }
   if (legacy_rndu(index) < aStrawhat / ((3 * bStrawhat - 2 * aStrawhat))) {
      /* sample from Strawhat part */
      z = aStrawhat * pow(legacy_rndu(index), 1.0 / 3.0);
   }
   else {
      /* sample from the box part */
      z = legacy_rndu(index) * (bStrawhat - aStrawhat) + aStrawhat;
   }
   return (legacy_rndu(index) < 0.5 ? -z : z);
}

static double rndBactrianTriangle(long index)
{
/* This returns a variate from the 1:1 mixture of two Triangle Tri(-m, 1-m^2) and Tri(m, 1-m^2),
   which has mean 0 and variance 1. 
*/
   double z = mBactrian + rndTriangle(index)*sBactrian;
   if (legacy_rndu(index) < 0.5) z = -z;
   return (z);
}

static double rndLaplace(long index)
{
   /* Standard Laplace variate, generated using inverse CDF  */
   double u, r;
   u = legacy_rndu(index) - 0.5;
   r = log(1 - 2 * fabs(u)) * 0.70710678118654752440;
   return (u >= 0 ? -r : r);
}

static double rndBactrianLaplace(long index)
{
   /* This returns a variate from the 1:1 mixture of two Laplace Lap(-m, 1-m^2) and Lap(m, 1-m^2),
      which has mean 0 and variance 1.
   */
   double z = mBactrian + rndLaplace(index)*sBactrian;
   if (legacy_rndu(index) < 0.5) z = -z;
   return (z);
}

double rndNormal(long index)
{
/* Standard normal variate, using the Box-Muller method (1958), improved by 
   Marsaglia and Bray (1964).  The method generates a pair of N(0,1) variates, 
   but only one is used.
   Johnson et al. (1994), Continuous univariate distributions, vol 1. p.153.
*/
   double u, v, s;

   for (; ;) {
      u = 2*legacy_rndu(index) - 1;
      v = 2*legacy_rndu(index) - 1;
      s = u*u + v*v;
      if (s>0 && s<1) break;
   }
   s = sqrt(-2*log(s)/s);
   return (u*s);  /* (v*s) is the other N(0,1) variate, wasted. */
}


double legacy_rnd_symmetrical(long index)
{
  #if 0
  return rndBactrianTriangle(index);
  #else
  return rndBactrianLaplace(index);
  /* return rndStrawhat(index); */
  #endif
}

double legacy_rndgamma (long index, double a)
{
/* This returns a random variable from gamma(a, 1).
   Marsaglia and Tsang (2000) A Simple Method for generating gamma variables", 
   ACM Transactions on Mathematical Software, 26 (3): 363-372.
   This is not entirely safe and is noted to produce zero when a is small (0.001).
 */
   double a0 = a, c, d, u, v, x, smallv=1E-300;

   if (a < 1) a++;

   d = a - 1.0 / 3.0;
   c = (1.0 / 3.0) / sqrt(d);

   for (; ; ) {
      do {
         x = rndNormal(index);
         v = 1.0 + c * x;
      } while (v <= 0);

      v *= v * v;
      u = legacy_rndu(index);

      if (u < 1 - 0.0331 * x * x * x * x)
         break;
      if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
         break;
   }
   v *= d;

   if (a0 < 1)    /* this may cause underflow if a is small, like 0.01 */
      v *= pow(legacy_rndu(index), 1 / a0);
   if (v == 0)   /* underflow */
      v = smallv;
   return v;
}

double legacy_rndbeta (long index, double p, double q)
{
/* this generates a random beta(p,q) variate
*/
   double gamma1, gamma2;
   gamma1 = legacy_rndgamma(index,p);
   gamma2 = legacy_rndgamma(index,q);
   return gamma1/(gamma1+gamma2);
}

void legacy_rnddirichlet(long index, double * output, double * alpha, long k)
{
  /* generates a random variate from the dirichlet distribution with K cats */

  long i;
  double s = 0;
  
  for (i = 0; i < k; ++i)
    s += output[i] = legacy_rndgamma(index,alpha[i]);

  for (i = 0; i < k; ++i)
    output[i] /= s;
}

long legacy_rndpoisson(long index, double m)
{
   /* m is the rate parameter of the poisson
      Numerical Recipes in C, 2nd ed. pp. 293-295
   */
   static double sq, alm, g, oldm = -1;
   double em, t, y;

   /* search from the origin
      if (m<5) {
         if (m!=oldm) { oldm=m; g=exp(-m); }
         y=rndu();  sq=alm=g;
         for (em=0; ; ) {
            if (y<sq) break;
            sq+= (alm*=m/ ++em);
         }
      }
   */
   if (m < 12) {
      if (m != oldm) { oldm = m; g = exp(-m); }
      em = -1; t = 1;
      for (; ;) {
         em++; t *= legacy_rndu(index);
         if (t <= g) break;
      }
   }
   else {
      if (m != oldm) {
         oldm = m;  sq = sqrt(2 * m);  alm = log(m);
         g = m*alm - lgamma(m + 1);
      }
      do {
         do {
            y = tan(3.141592654*legacy_rndu(index));
            em = sq*y + m;
         } while (em < 0);
         em = floor(em);
         t = 0.9*(1 + y*y)*exp(em*alm - lgamma(em + 1) - g);
      } while (legacy_rndu(index) > t);
   }
   return ((long)em);
}
