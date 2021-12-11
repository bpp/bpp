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

#define square(a) ((a)*(a))
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define swap2(a,b,y) { y=a; a=b; b=y; }

static int noisy = 0;
static int Iround = 0;
static int NFunCall = 0;
static int AlwaysCenter = 0;
static double SIZEp = 0;
static double Small_Diff = 1e-6;  /* reasonable values 1e-5, 1e-7 */

static int zero(double x[], int n)
{
   for (int i = 0; i < n; i++) x[i] = 0; 
   return (0);
}

static int xtoy(double x[], double y[], int n)
{
   for (int i = 0; i < n; i++) y[i] = x[i];
   return(0);
}

static int identity(double x[], int n)
{
   for (int i = 0; i < n; i++) { 
     for (int j = 0; j < n; j++)
       x[i * n + j] = 0;
     x[i * n + i] = 1;
   }  
   return (0);
}

static double distance(double x[], double y[], int n)
{
   double t = 0;
   for (int i = 0; i < n; i++) t += square(x[i] - y[i]);
   return(sqrt(t));
}

static double innerp(double x[], double y[], int n)
{
   double t = 0;  for (int i = 0; i < n; i++)  t += x[i] * y[i];  return(t);
}

static double norm(double x[], int n)
{
   double t = 0;  for (int i = 0; i < n; i++)  t += x[i] * x[i];  return sqrt(t);
}

static int H_end(double x0[], double x1[], double f0, double f1, double e1, double e2, int n)
{   
/*  Himmelblau termination rule.   return 1 for stop, 0 otherwise.
*/
   double r;
   if ((r = norm(x0, n)) < e2)  r = 1;
   r *= e1;
   if (distance(x1, x0, n) >= r)  return (0);
   r = fabs(f0);  if (r < e2) r = 1;
   r *= e1;
   if (fabs(f1 - f0) >= r) return (0);
   return (1);
}

static double fun_LineSearch(double t, double(*fun)(double x[], int n), double x0[], double p[], double x[], int n)
{
   int i;
   for (i = 0; i < n; i++) x[i] = x0[i] + t * p[i];
   return((*fun)(x, n));
}

static double LineSearch2(double(*fun)(double x[], int n),
                          double *f, double x0[], double p[], double step,
                          double limit, double e, double space[], int n)
{
   /* linear search using quadratic interpolation
      from x0[] in the direction of p[],
                   x = x0 + a*p        a ~(0,limit)
      returns (a).    *f: f(x0) for input and f(x) for output

      x0[n] x[n] p[n] space[n]

      adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
      optimization: An introduction.  Van Nostrand Reinhold Company, New York.
      pp. 62-73.
      step is used to find the bracket and is increased or reduced as necessary,
      and is not terribly important.
   */
   int ii = 0, maxround = 10, status, i, nsymb = 0;
   double *x = space, factor = 4, smallv = 1e-10, smallgapa = 0.2;
   double a0, a1, a2, a3, a4 = -1, a5, a6, f0, f1, f2, f3, f4 = -1, f5, f6;

   /* look for bracket (a1, a2, a3) with function values (f1, f2, f3)
      step length step given, and only in the direction a>=0
   */

   if (noisy > 2)
      printf("\n%3d h-m-p %7.4f %6.4f %8.4f ", Iround + 1, step, limit, norm(p, n));

   if (step <= 0 || limit < smallv || step >= limit) {
      if (noisy > 2)
         printf("\nh-m-p:%20.8e%20.8e%20.8e %12.6f\n", step, limit, norm(p, n), *f);
      return (0);
   }
   a0 = a1 = 0; f1 = f0 = *f;
   a2 = a0 + step; f2 = fun_LineSearch(a2, fun, x0, p, x, n);
   if (f2 > f1) {  /* reduce step length so the algorithm is decreasing */
      for (; ;) {
         step /= factor;
         if (step < smallv) return (0);
         a3 = a2;    f3 = f2;
         a2 = a0 + step;  f2 = fun_LineSearch(a2, fun, x0, p, x, n);
         if (f2 <= f1) break;
         if (noisy > 2) { printf("-"); nsymb++; }
      }
   }
   else {       /* step length is too small? */
      for (; ;) {
         step *= factor;
         if (step > limit) step = limit;
         a3 = a0 + step;  f3 = fun_LineSearch(a3, fun, x0, p, x, n);
         if (f3 >= f2) break;

         if (noisy > 2) { printf("+"); nsymb++; }
         a1 = a2; f1 = f2;    a2 = a3; f2 = f3;
         if (step >= limit) {
            if (noisy > 2) for (; nsymb < 5; nsymb++) printf(" ");
            if (noisy > 2) printf(" %12.6f%3c %6.4f %5d", *f = f3, 'm', a3, NFunCall);
            *f = f3; return(a3);
         }
      }
   }

   /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
   for (ii = 0; ii < maxround; ii++) {
      /* a4 is the minimum from the parabola over (a1,a2,a3)  */
      a4 = (a2 - a3)*f1 + (a3 - a1)*f2 + (a1 - a2)*f3;
      if (fabs(a4) > 1e-100)
         a4 = ((a2*a2 - a3*a3)*f1 + (a3*a3 - a1*a1)*f2 + (a1*a1 - a2*a2)*f3) / (2 * a4);
      if (a4 > a3 || a4 < a1) {   /* out of range */
         a4 = (a1 + a2) / 2;
         status = 'N';
      }
      else {
         if ((a4 <= a2 && a2 - a4 > smallgapa*(a2 - a1)) || (a4 > a2 && a4 - a2 > smallgapa*(a3 - a2)))
            status = 'Y';
         else
            status = 'C';
      }
      f4 = fun_LineSearch(a4, fun, x0, p, x, n);
      if (noisy > 2) putchar(status);
      if (fabs(f2 - f4) < e*(1 + fabs(f2))) {
         if (noisy > 2)
            for (nsymb += ii + 1; nsymb < 5; nsymb++) printf(" ");
         break;
      }

      /* possible multiple local optima during line search */
      if (noisy > 2 && ((a4<a2&&f4>f1) || (a4 > a2&&f4 > f3))) {
         printf("\n\na %12.6f %12.6f %12.6f %12.6f", a1, a2, a3, a4);
         printf("\nf %12.6f %12.6f %12.6f %12.6f\n", f1, f2, f3, f4);

         for (a5 = a1; a5 <= a3; a5 += (a3 - a1) / 20) {
            printf("\t%.6e ", a5);
            if (n < 5) for (i = 0; i < n; i++)  printf("\t%.6f", x0[i] + a5*p[i]);
            printf("\t%.6f\n", fun_LineSearch(a5, fun, x0, p, x, n));
         }
         puts("Linesearch2 a4: multiple optima?");
      }
      if (a4 <= a2) {    /* fig 2.2.10 */
         if (a2 - a4 > smallgapa*(a2 - a1)) {
            if (f4 <= f2) { a3 = a2; a2 = a4;  f3 = f2; f2 = f4; }
            else { a1 = a4; f1 = f4; }
         }
         else {
            if (f4 > f2) {
               a5 = (a2 + a3) / 2; f5 = fun_LineSearch(a5, fun, x0, p, x, n);
               if (f5 > f2) { a1 = a4; a3 = a5;  f1 = f4; f3 = f5; }
               else { a1 = a2; a2 = a5;  f1 = f2; f2 = f5; }
            }
            else {
               a5 = (a1 + a4) / 2; f5 = fun_LineSearch(a5, fun, x0, p, x, n);
               if (f5 >= f4)
               {
                  a3 = a2; a2 = a4; a1 = a5;  f3 = f2; f2 = f4; f1 = f5;
               }
               else {
                  a6 = (a1 + a5) / 2; f6 = fun_LineSearch(a6, fun, x0, p, x, n);
                  if (f6 > f5)
                  {
                     a1 = a6; a2 = a5; a3 = a4;  f1 = f6; f2 = f5; f3 = f4;
                  }
                  else { a2 = a6; a3 = a5; f2 = f6; f3 = f5; }
               }
            }
         }
      }
      else {                     /* fig 2.2.9 */
         if (a4 - a2 > smallgapa*(a3 - a2)) {
            if (f2 >= f4) { a1 = a2; a2 = a4;  f1 = f2; f2 = f4; }
            else { a3 = a4; f3 = f4; }
         }
         else {
            if (f4 > f2) {
               a5 = (a1 + a2) / 2; f5 = fun_LineSearch(a5, fun, x0, p, x, n);
               if (f5 > f2) { a1 = a5; a3 = a4;  f1 = f5; f3 = f4; }
               else { a3 = a2; a2 = a5;  f3 = f2; f2 = f5; }
            }
            else {
               a5 = (a3 + a4) / 2; f5 = fun_LineSearch(a5, fun, x0, p, x, n);
               if (f5 >= f4)
               {
                  a1 = a2; a2 = a4; a3 = a5;  f1 = f2; f2 = f4; f3 = f5;
               }
               else {
                  a6 = (a3 + a5) / 2; f6 = fun_LineSearch(a6, fun, x0, p, x, n);
                  if (f6 > f5)
                  {
                     a1 = a4; a2 = a5; a3 = a6;  f1 = f4; f2 = f5; f3 = f6;
                  }
                  else { a1 = a5; a2 = a6;  f1 = f5; f2 = f6; }
               }
            }
         }
      }
   }

   if (f2 > f0 && f4 > f0)  a4 = 0;
   if (f2 <= f4) { *f = f2; a4 = a2; }
   else         *f = f4;
   if (noisy > 2) printf(" %12.6f%3d %6.4f %5d", *f, ii, a4, NFunCall);

   return (a4);
}





static int gradientB(int n, double x[], double f0, double g[],
                     double(*fun)(double x[], int n),
                     double space[], int xmark[])
{
   /* f0=fun(x) is always provided.
   xmark=0: central; 1: upper; -1: down
   */
   int i, j;
   double *x0 = space, *x1 = space + n, eh0 = Small_Diff, eh;  /* eh0=1e-6 || 1e-7 */

   for (i = 0; i < n; i++) {
      eh = eh0*(fabs(x[i]) + 1);
      if (xmark[i] == 0 && (AlwaysCenter || SIZEp < 1)) {   /* central */
         for (j = 0; j < n; j++)  x0[j] = x1[j] = x[j];
         eh = pow(eh, .67);   x0[i] -= eh;  x1[i] += eh;
         g[i] = ((*fun)(x1, n) - (*fun)(x0, n)) / (eh*2.0);
      }
      else {                                              /* forward or backward */
         for (j = 0; j < n; j++)  x1[j] = x[j];
         if (xmark[i]) eh *= -xmark[i];
         x1[i] += eh;
         g[i] = ((*fun)(x1, n) - f0) / eh;
      }
   }
   return(0);
}

#define BFGS
/*
#define SR1
#define DFP
*/

int ming2(FILE *fout, double *f, double(*fun)(double x[], int n),
          int(*dfun)(double x[], double *f, double dx[], int n),
          double x[], double xb[][2], double space[], double e, int n)
{
   /* n-variate minimization with bounds using the BFGS algorithm
        g0[n] g[n] p[n] x0[n] y[n] s[n] z[n] H[n*n] C[n*n] tv[2*n]
        xmark[n],ix[n]
      Size of space should be (check carefully?)
         #define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(double))
      nfree: # free variables
      xmark[i]=0 for inside space; -1 for lower boundary; 1 for upper boundary.
      x[] has initial values at input and returns the estimates in return.
      ix[i] specifies the i-th free parameter

   */
   int i, j, i1, i2, it, maxround = 10000, fail = 0, *xmark, *ix, nfree;
   int Ngoodtimes = 2, goodtimes = 0;
   double smallv = 1.e-30, sizep0 = 0;     /* small value for checking |w|=0 */
   double f0, *g0, *g, *p, *x0, *y, *s, *z, *H, *C, *tv;
   double w, v, alpha, am, h, maxstep = 8;

   if (n == 0) return(0);
   g0 = space;   g = g0 + n;  p = g + n;   x0 = p + n;
   y = x0 + n;     s = y + n;   z = s + n;   H = z + n;  C = H + n*n, tv = C + n*n;
   xmark = (int*)(tv + 2 * n);  ix = xmark + n;

   for (i = 0; i < n; i++) { xmark[i] = 0; ix[i] = i; }
   for (i = 0, nfree = 0; i < n; i++) {
      if (x[i] <= xb[i][0]) { x[i] = xb[i][0]; xmark[i] = -1; continue; }
      if (x[i] >= xb[i][1]) { x[i] = xb[i][1]; xmark[i] = 1; continue; }
      ix[nfree++] = i;
   }
   if (noisy > 2 && nfree < n && n < 50) {
      printf("\n"); for (j = 0; j < n; j++) printf(" %9.6f", x[j]);  printf("\n");
      for (j = 0; j < n; j++) printf(" %9.5f", xb[j][0]);
      printf("\n");
      for (j = 0; j < n; j++) printf(" %9.5f", xb[j][1]);
      printf("\n");
      if (nfree < n && noisy >= 3) printf("warning: ming2, %d paras at boundary.", n - nfree);
   }

   f0 = *f = (*fun)(x, n);
   xtoy(x, x0, n);
   SIZEp = 99;
   if (noisy > 2) {
      printf("\nIterating by ming2\nInitial: fx= %12.6f\nx=", f0);
      for (i = 0; i < n; i++) printf(" %8.5f", x[i]);
      printf("\n");
   }

   if (dfun)  (*dfun) (x0, &f0, g0, n);
   else       gradientB(n, x0, f0, g0, fun, tv, xmark);

   identity(H, nfree);
   for (Iround = 0; Iround < maxround; Iround++) {
      if (fout) {
         fprintf(fout, "\n%3d %7.4f %13.6f  x: ", Iround, sizep0, f0);
         for (i = 0; i < n; i++) fprintf(fout, "%8.5f  ", x0[i]);
         fflush(fout);
      }

      for (i = 0, zero(p, n); i < nfree; i++)  for(j=0; j<nfree; j++)
         p[ix[i]] -= H[i*nfree + j] * g0[ix[j]];
      sizep0 = SIZEp;
      SIZEp = norm(p, n);      /* check this */

      for (i = 0, am = maxstep; i < n; i++) {  /* max step length */
         if (p[i] > 0 && (xb[i][1] - x0[i]) / p[i] < am) am = (xb[i][1] - x0[i]) / p[i];
         else if (p[i] < 0 && (xb[i][0] - x0[i]) / p[i] < am) am = (xb[i][0] - x0[i]) / p[i];
      }

      if (Iround == 0) {
         h = fabs(2 * f0*.01 / innerp(g0, p, n));  /* check this?? */
         h = min2(h, am / 2000);

      }
      else {
         h = norm(s, nfree) / SIZEp;
         h = max2(h, am / 500);
      }
      h = max2(h, 1e-5);   h = min2(h, am / 5);
      *f = f0;
      alpha = LineSearch2(fun, f, x0, p, h, am, min2(1e-3, e), tv, n); /* n or nfree? */

      if (alpha <= 0) {
         if (fail) {
            if (AlwaysCenter) { Iround = maxround;  break; }
            else { AlwaysCenter = 1; identity(H, n); fail = 1; }
         }
         else
         {
            if (noisy > 2) printf(".. ");
            identity(H, nfree);
            fail = 1;
         }
      }
      else {
         fail = 0;
         for (i = 0; i < n; i++)  x[i] = x0[i] + alpha*p[i];
         w = min2(2, e * 1000); if (e<1e-4 && e>1e-6) w = 0.01;

         if (Iround == 0 || SIZEp < sizep0 || (SIZEp < .001 && sizep0 < .001)) goodtimes++;
         else  goodtimes = 0;
         if ((n == 1 || goodtimes >= Ngoodtimes) && SIZEp < (e > 1e-5 ? 1 : .001)
            && H_end(x0, x, f0, *f, e, e, n))
            break;
      }
      if (dfun)
         (*dfun) (x, f, g, n);
      else
         gradientB(n, x, *f, g, fun, tv, xmark);

      /* modify the working set */
      for (i = 0; i < n; i++) {         /* add constraints, reduce H */
         if (xmark[i]) continue;
         if (fabs(x[i] - xb[i][0]) < 1e-6 && -g[i] < 0)  xmark[i] = -1;
         else if (fabs(x[i] - xb[i][1]) < 1e-6 && -g[i] > 0)  xmark[i] = 1;
         if (xmark[i] == 0) continue;
         xtoy(H, C, nfree*nfree);
         for (it = 0; it < nfree; it++) if (ix[it] == i) break;
         for (i1 = it; i1 < nfree - 1; i1++) ix[i1] = ix[i1 + 1];
         for (i1 = 0, nfree--; i1 < nfree; i1++) for (i2 = 0; i2 < nfree; i2++)
            H[i1*nfree + i2] = C[(i1 + (i1 >= it))*(nfree + 1) + i2 + (i2 >= it)];
      }
      for (i = 0, it = 0, w = 0; i < n; i++) {  /* delete a constraint, enlarge H */
         if (xmark[i] == -1 && -g[i] > w) { it = i; w = -g[i]; }
         else if (xmark[i] == 1 && -g[i] < -w) { it = i; w = g[i]; }
      }
      if (w > 10 * SIZEp / nfree) {          /* *** */
         xtoy(H, C, nfree*nfree);
         for (i1 = 0; i1 < nfree; i1++) for (i2 = 0; i2 < nfree; i2++) 
           H[i1 * (nfree + 1) + i2] = C[i1 * nfree + i2];
         for (i1 = 0; i1 < nfree; i1++) H[i1*(nfree + 1) + nfree] = H[nfree*(nfree + 1) + i1] = 0;
         H[(nfree + 1)*(nfree + 1) - 1] = 1;
         xmark[it] = 0;   ix[nfree++] = it;
      }

      if (noisy > 2) {
         printf(" | %d/%d", n - nfree, n);
         /* for (i = 0; i < n; i++)  if (xmark[i]) printf ("%4d", i+1); */
      }
      for (i = 0, f0 = *f; i < nfree; i++)
      {
         y[i] = g[ix[i]] - g0[ix[i]];  s[i] = x[ix[i]] - x0[ix[i]];
      }
      for (i = 0; i < n; i++) { g0[i] = g[i]; x0[i] = x[i]; }


      /* renewal of H varies with different algorithms   */
#if (defined SR1)
      /*   Symmetrical Rank One (Broyden, C. G., 1967) */
      for (i = 0, w = .0; i < nfree; i++) {
         for (j = 0, v = .0; j < nfree; j++) v += H[i*nfree + j] * y[j];
         z[i] = s[i] - v;
         w += y[i] * z[i];
      }
      if (fabs(w) < smallv) { identity(H, nfree); fail = 1; continue; }
      for (i = 0; i < nfree; i++)  for (j = 0; j < n; j++)  H[i*nfree + j] += z[i] * z[j] / w;
#elif (defined DFP)
      /* Davidon (1959), Fletcher and Powell (1963). */
      for (i = 0, w = v = 0.; i < nfree; i++) {
         for (j = 0, z[i] = 0; j < nfree; j++) z[i] += H[i*nfree + j] * y[j];
         w += y[i] * z[i];  v += y[i] * s[i];
      }
      if (fabs(w) < smallv || fabs(v) < smallv) { identity(H, nfree); fail = 1; continue; }
      for (i = 0; i < nfree; i++) for (j = 0; j < nfree; j++)
         H[i*nfree + j] += s[i] * s[j] / v - z[i] * z[j] / w;
#else /* BFGS */
      for (i = 0, w = v = 0.; i < nfree; i++) {
         for (j = 0, z[i] = 0.; j < nfree; j++) z[i] += H[i*nfree + j] * y[j];
         w += y[i] * z[i];    v += y[i] * s[i];
      }
      if (fabs(v) < smallv) { identity(H, nfree); fail = 1; continue; }
      for (i = 0; i < nfree; i++) for (j = 0; j < nfree; j++)
         H[i*nfree + j] += ((1 + w / v)*s[i] * s[j] - z[i] * s[j] - s[i] * z[j]) / v;
#endif
   }    /* for (Iround,maxround)  */

   /* try to remove this after updating LineSearch2() */
   *f = (*fun)(x, n);
   if (noisy > 2) printf("\n");

   if (Iround == maxround) {
      if (fout) fprintf(fout, "\ncheck convergence!\n");
      return(-1);
   }
   if (nfree == n) {
      xtoy(H, space, n*n);  /* H has variance matrix, or inverse of Hessian */
      return(1);
   }
   return(0);
}
