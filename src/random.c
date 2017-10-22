/*
    Copyright (C) 2016 Tomas Flouri and Ziheng Yang

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

#define mBactrian  0.95
#define sBactrian  sqrt(1-mBactrian*mBactrian)

/* legacy random number generators */
static unsigned int z_rndu = 666;

void legacy_init()
{
  z_rndu = (unsigned int)opt_seed;
}

double legacy_rndu()
{
/* 32-bit integer assumed.
   From Ripley (1987) p. 46 or table 2.4 line 2. 
   This may return 0 or 1, which can be a problem.
*/
   z_rndu = z_rndu*69069 + 1;
   if(z_rndu==0 || z_rndu==4294967295)  z_rndu = 13;
   return z_rndu/4294967295.0;
}

static double rndTriangle()
{
	double u, z;
/* Standard Triangle variate, generated using inverse CDF  */
	u = legacy_rndu();
	if(u > 0.5)
		z =  sqrt(6.0) - 2.0*sqrt(3.0*(1.0 - u));
   else
		z = -sqrt(6.0) + 2.0*sqrt(3.0*u);
	return z;
}

static double rndBactrianTriangle()
{
/* This returns a variate from the 1:1 mixture of two Triangle Tri(-m, 1-m^2) and Tri(m, 1-m^2),
   which has mean 0 and variance 1. 
*/
   double z = mBactrian + rndTriangle()*sBactrian;
   if (legacy_rndu() < 0.5) z = -z;
   return (z);
}

static double rndNormal (void)
{
/* Standard normal variate, using the Box-Muller method (1958), improved by 
   Marsaglia and Bray (1964).  The method generates a pair of N(0,1) variates, 
   but only one is used.
   Johnson et al. (1994), Continuous univariate distributions, vol 1. p.153.
*/
   double u, v, s;

   for (; ;) {
      u = 2*legacy_rndu() - 1;
      v = 2*legacy_rndu() - 1;
      s = u*u + v*v;
      if (s>0 && s<1) break;
   }
   s = sqrt(-2*log(s)/s);
   return (u*s);  /* (v*s) is the other N(0,1) variate, wasted. */
}


double legacy_rnd_symmetrical()
{
  return rndBactrianTriangle();
}

double legacy_rndgamma (double a)
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
         x = rndNormal( );
         v = 1.0 + c * x;
      } while (v <= 0);

      v *= v * v;
      u = legacy_rndu( );

      if (u < 1 - 0.0331 * x * x * x * x)
         break;
      if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
         break;
   }
   v *= d;

   if (a0 < 1)    /* this may cause underflow if a is small, like 0.01 */
      v *= pow(legacy_rndu( ), 1 / a0);
   if (v == 0)   /* underflow */
      v = smallv;
   return v;
}

double legacy_rndbeta (double p, double q)
{
/* this generates a random beta(p,q) variate
*/
   double gamma1, gamma2;
   gamma1 = legacy_rndgamma(p);
   gamma2 = legacy_rndgamma(q);
   return gamma1/(gamma1+gamma2);
}
