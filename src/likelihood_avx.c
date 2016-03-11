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

static double lmatrix[4*4] __attribute__ ((aligned (32)));
static double rmatrix[4*4] __attribute__ ((aligned (32)));

void pmat_jc69(double * pmatrix, double t)
{
  if (t < -0.0001)
    printf ("\nt = %.5f in pijJC69", t);
  
  if (t < 1e-100)
  {
    pmatrix[0]  = 1;
    pmatrix[1]  = 0;
    pmatrix[2]  = 0;
    pmatrix[3]  = 0;

    pmatrix[4]  = 0;
    pmatrix[5]  = 1;
    pmatrix[6]  = 0;
    pmatrix[7]  = 0;

    pmatrix[8]  = 0;
    pmatrix[9]  = 0;
    pmatrix[10] = 1;
    pmatrix[11] = 0;

    pmatrix[12] = 0;
    pmatrix[13] = 0;
    pmatrix[14] = 0;
    pmatrix[15] = 1;
  }
  else
  {
    double a =  (1 + 3*exp(-4*t/3) ) / 4;
    double b = (1 - a) / 3;

    pmatrix[0]  = a;
    pmatrix[1]  = b;
    pmatrix[2]  = b;
    pmatrix[3]  = b;

    pmatrix[4]  = b;
    pmatrix[5]  = a;
    pmatrix[6]  = b;
    pmatrix[7]  = b;

    pmatrix[8]  = b;
    pmatrix[9]  = b;
    pmatrix[10] = a;
    pmatrix[11] = b;

    pmatrix[12] = b;
    pmatrix[13] = b;
    pmatrix[14] = b;
    pmatrix[15] = a;
  }
}

static void tip_tip(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int n,i;

  __m256d ymm0,ymm1,ymm2;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  for (n = 0; n < sites; ++n)
  {
    double * lmat = lmatrix + com.z[lchild][n] * states;
    double * rmat = rmatrix + com.z[rchild][n] * states;

    //for (i = 0; i < states; ++i)
    //  clv[i] = lmat[i] * rmat[i];
    
    ymm0 = _mm256_load_pd(lmat);
    ymm1 = _mm256_load_pd(rmat);
    ymm2 = _mm256_mul_pd(ymm0,ymm1);

    _mm256_store_pd(clv,ymm2);


    clv += states;
  }
}
#if 0
static void inner_inner(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h,j,k;
  double x,y;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  double * lclv = nodes[lchild].conP;
  double * rclv = nodes[rchild].conP;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix;
    double * rmat = rmatrix;
    for (j = 0; j < states; ++j)
    {
      x = 0;
      y = 0;
      for (k = 0; k < states; ++k)
      {
        x += lmat[k] * lclv[k];
        y += rmat[k] * rclv[k];
      }

      lmat += states;
      rmat += states;

      clv[j] = x*y;
    }
    clv += states;
    lclv += states;
    rclv += states;
  }
}
#else
static void inner_inner(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h,j,k;
  double x,y;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  double * lclv = nodes[lchild].conP;
  double * rclv = nodes[rchild].conP;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix;
    double * rmat = rmatrix;

    /* compute vector of x */
    xmm4 = _mm256_load_pd(lmat);
    xmm5 = _mm256_load_pd(lclv);
    xmm0 = _mm256_mul_pd(xmm4,xmm5);

    ymm4 = _mm256_load_pd(rmat);
    ymm5 = _mm256_load_pd(rclv);
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

    _mm256_store_pd(clv, xmm0);

    clv += states;
    lclv += states;
    rclv += states;
  }
}
#endif
static void tip_inner(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h,j,k;
  double y;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  double * rclv = nodes[rchild].conP;
  
  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm4;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix + com.z[lchild][h] * states;
    double * rmat = rmatrix;

    xmm4 = _mm256_load_pd(lmat);

    ymm4 = _mm256_load_pd(rmat);
    ymm5 = _mm256_load_pd(rclv);
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

    _mm256_store_pd(clv, xmm0);

    clv += states;
    rclv += states;
  }
}

int ConditonalPNode (int inode)
{
  int states = 4;
  int sites = com.npatt;
  int i;
  int child, lchild, rchild;

  /* recursive call the ConditionalPNode on all children of the current node */
  for (i = 0; i < nodes[inode].nson; ++i)
  {
    /* get id of i-th child of current node */
    child = nodes[inode].sons[i];

    /* if child is an inner node then recursively compute its conditional
       probabilities vector */
    if (nodes[child].nson > 0 && (!mcmc.saveconP || !com.oldconP[child]))
      ConditonalPNode(nodes[inode].sons[i]);
  }

  /* initialize CLV entries of current node to 1 */
  int n = sites * states;
  for (i = 0; i < n; ++i)
    nodes[inode].conP[i] = 1;

  if (nodes[inode].nson == 0) return (0);

  lchild = nodes[inode].sons[0];
  rchild = nodes[inode].sons[1];

  int ltip = (nodes[lchild].nson == 0);
  int rtip = (nodes[rchild].nson == 0);

  if (ltip && rtip)
    tip_tip(lchild,rchild,nodes[inode].conP);
  else if (ltip)
    tip_inner(lchild,rchild,nodes[inode].conP);
  else if (rtip)
    tip_inner(rchild,lchild,nodes[inode].conP);
  else
    inner_inner(lchild,rchild,nodes[inode].conP);

  return (0);
}

#if 0
double lnpD_locus (int locus)
{
/* This calculates ln{Di|Gi, Bi} using nodes[].age and tree.root.
   Note that nodes[].branch may not be kept up to date in the program.
*/
   int  h, i, haserr=0;
   double lnL=0, fh = 0;
   int states = 4;

   if(!mcmc.usedata) return(0);
   for(i=0; i<tree.nnode; i++) {
      if(i==tree.root) continue;   
      nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
      if(data.est_locusrate)
         nodes[i].branch *= data.locusrate[locus];
      if(nodes[i].branch < 0) {
         printf("branch length = %.6f < 0\n", nodes[i].branch);
         exit(-1);
      }
   }
   ConditonalPNode(tree.root);
   /*
   for(h=0; h<com.npatt; h++) {
      for(i=0,fh=0; i<com.ncode; i++) {
         fh += com.pi[i]*nodes[tree.root].conP[h*com.ncode+i];
      }
      if(fh<1E-300)  lnL += com.fpatt[h]*(-500);
      else           lnL += com.fpatt[h]*log(fh);
   }
   */
   
   int vec = com.npatt & 0xFFFFFFFC
   double * clv = nodes[tree.root].conP;
   double * fpatt = com.fpatt;

   __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8,xmm9;

   xmm8 = _mm256_load_pd(com.pi);
   xmm9 = _mm256_setzero_pd();
   fh = 0;
   for (h = 0; h < vec; h += 4)
   {
     xmm4 = _mm256_load_pd(clv+0);
     xmm0 = _mm256_mul_pd(xmm4,xmm8);

     xmm4 = _mm256_load_pd(clv+4);
     xmm1 = _mm256_mul_pd(xmm4,xmm8);
     
     xmm4 = _mm256_load_pd(clv+8);
     xmm2 = _mm256_mul_pd(xmm4,xmm8);

     xmm4 = _mm256_load_pd(clv+12);
     xmm3 = _mm256_mul_pd(xmm4,xmm8);

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

     xmm0 = _mm256_log_pd(xmm4);
     xmm1 = _mm256_load_pd(fpatt);
     xmm2 = _mm256_mul_pd(xmm0,xmm1);

     xmm9 = _mm256_add_pd(xmm9,xmm2);

     clv += 4*states;
     fpatt += 4;


   }
   double* res = (double*)&xmm9;
   fh += res[0] + res[1] + res[2] + res[3];


   for(h=0; h<com.npatt-vec; h++) {
      for(i=0,i<com.ncode; i++) {
         fh += com.pi[i]*clv[h*com.ncode+i];
      }
      if(fh<1E-300)  lnL += com.fpatt[h]*(-500);
      else           lnL += com.fpatt[h]*log(fh);
   }
   return (lnL);
}
#endif
