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

static double lmatrix[4*4] __attribute__ ((aligned (16)));
static double rmatrix[4*4] __attribute__ ((aligned (16)));

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

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  __m128d xmm0,xmm1,xmm2;

  for (n = 0; n < sites; ++n)
  {
    double * lmat = lmatrix + (com.z[lchild][n] << 2); //*states;
    double * rmat = rmatrix + (com.z[rchild][n] << 2); // *states;

    //for (i = 0; i < states; ++i)
    //  clv[i] = lmat[i] * rmat[i];

    xmm0 = _mm_load_pd(lmat);
    xmm1 = _mm_load_pd(rmat);
    xmm2 = _mm_mul_pd(xmm0,xmm1);

    _mm_store_pd(clv,xmm2);

    xmm0 = _mm_load_pd(lmat+2);
    xmm1 = _mm_load_pd(rmat+2);
    xmm2 = _mm_mul_pd(xmm0,xmm1);

    _mm_store_pd(clv+2,xmm2);

    clv += states;
  }
}

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

  __m128d xmm0,xmm1,xmm2,xmm3,xmm4;
  __m128d xmm5,xmm6,xmm7,xmm8,xmm9;

  /* 

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
    | c1 | c2 | c3 | c4 |         | c1 | c2 | c3 | c4 |
    +----+----+----+----+         +----+----+----+----+

  */

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix;
    double * rmat = rmatrix;

    /* compute left */
    xmm0 = _mm_load_pd(lclv);              /* we need it */
    xmm2 = _mm_load_pd(lmat);
    xmm3 = _mm_mul_pd(xmm0,xmm2);

    xmm1 = _mm_load_pd(lclv+2);            /* we need it */
    xmm2 = _mm_load_pd(lmat+2);
    xmm4 = _mm_mul_pd(xmm1,xmm2);

    /* calculate (a1*c1 + a2*c2 | a3*c3 + a4*c4) */
    xmm2 = _mm_hadd_pd(xmm3,xmm4);         /* we need it (1) */

    /* compute right */
    xmm3 = _mm_load_pd(rclv);              /* we need it */
    xmm4 = _mm_load_pd(rmat);
    xmm6 = _mm_mul_pd(xmm3,xmm4);

    xmm4 = _mm_load_pd(rclv+2);            /* we need it */
    xmm7 = _mm_load_pd(rmat+2);
    xmm8 = _mm_mul_pd(xmm4,xmm7);

    /* calculate (b1*d1 + b2*d2 | b3*d3 + b4*d4) */
    xmm5 = _mm_hadd_pd(xmm6,xmm8);         /* we need it (2) */

    rmat += states;
    lmat += states;

    /* compute left */
    xmm6 = _mm_load_pd(lmat);
    xmm7 = _mm_mul_pd(xmm0,xmm6);

    xmm6 = _mm_load_pd(lmat+2);
    xmm8 = _mm_mul_pd(xmm1,xmm6);

    /* calculate (a5*c1 + a6*c2 | a7*c3 + a8*c4) */
    xmm6 = _mm_hadd_pd(xmm7,xmm8);         /* we need it (3) */

    /* compute right */
    xmm7 = _mm_load_pd(rmat);
    xmm8 = _mm_mul_pd(xmm3,xmm7);

    xmm7 = _mm_load_pd(rmat+2);
    xmm9 = _mm_mul_pd(xmm4,xmm7);
    
    /* calculate (b5*d1 + b6*d2 | b7*d3 + b8*d4) */
    xmm7 = _mm_hadd_pd(xmm8,xmm9);         /* we need it (4) */
    
    xmm8 = _mm_hadd_pd(xmm2,xmm6);
    xmm2 = _mm_hadd_pd(xmm5,xmm7);
    xmm9 = _mm_mul_pd(xmm8,xmm2);
    _mm_store_pd(clv, xmm9);

    rmat += states;
    lmat += states;

    /* compute left */
    xmm5 = _mm_load_pd(lmat);
    xmm6 = _mm_mul_pd(xmm0,xmm5);

    xmm5 = _mm_load_pd(lmat+2);
    xmm7 = _mm_mul_pd(xmm1,xmm5);

    /* calculate (a9*c1 + a10*c2 | a11*c3 + a12*c4) */
    xmm2 = _mm_hadd_pd(xmm6,xmm7);         /* we need it (1) */

    /* compute right */
    xmm6 = _mm_load_pd(rmat);
    xmm7 = _mm_mul_pd(xmm3,xmm6);

    xmm6 = _mm_load_pd(rmat+2);
    xmm8 = _mm_mul_pd(xmm4,xmm6);
    
    /* calculate (b9*d1 + b10*d2 | b11*d3 + b12*d4) */
    xmm5 = _mm_hadd_pd(xmm7,xmm8);         /* we need it (2) */
    
    rmat += states;
    lmat += states;

    /* compute left */
    xmm6 = _mm_load_pd(lmat);
    xmm7 = _mm_mul_pd(xmm0,xmm6);

    xmm6 = _mm_load_pd(lmat+2);
    xmm8 = _mm_mul_pd(xmm1,xmm6);

    /* calculate (a13*c1 + a14*c2 | a15*c3 + a16*c4) */
    xmm6 = _mm_hadd_pd(xmm7,xmm8);         /* we need it (3) */

    /* compute right */
    xmm7 = _mm_load_pd(rmat);
    xmm8 = _mm_mul_pd(xmm3,xmm7);

    xmm7 = _mm_load_pd(rmat+2);
    xmm9 = _mm_mul_pd(xmm4,xmm7);
    
    /* calculate (b13*d1 + b14*d2 | b15*d3 + b16*d4) */
    xmm7 = _mm_hadd_pd(xmm8,xmm9);         /* we need it (4) */
    
    xmm8 = _mm_hadd_pd(xmm2,xmm6);
    xmm2 = _mm_hadd_pd(xmm5,xmm7);
    xmm9 = _mm_mul_pd(xmm8,xmm2);
    _mm_store_pd(clv+2, xmm9);

    clv  += states;
    lclv += states;
    rclv += states;
  }
}

static void tip_inner(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h,j,k;
  double y;

  __m128d xmm0,xmm1,xmm2,xmm3,xmm4;
  __m128d xmm5,xmm6,xmm7,xmm8,xmm9;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  double * rclv = nodes[rchild].conP;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix + (com.z[lchild][h] << 2); // *states;
    double * rmat = rmatrix;

    /*
    for (j = 0; j < states; ++j)
    {
      y = 0;
      for (k = 0; k < states; ++k)
        y += rmat[k] * rclv[k];

      rmat += states;

      clv[j] = y*lmat[j];
    }
    */

    xmm0 = _mm_load_pd(rclv);           /* we need it */
    xmm2 = _mm_load_pd(rmat);
    xmm3 = _mm_mul_pd(xmm0,xmm2);

    xmm1 = _mm_load_pd(rclv+2);         /* we need it */
    xmm2 = _mm_load_pd(rmat+2);
    xmm4 = _mm_mul_pd(xmm1,xmm2);

    /* calculate (b1*d1 + b2*d2 | b3*d3 + b4*d4) */
    xmm2 = _mm_hadd_pd(xmm3,xmm4);      /* we need it (1) */

    rmat += states;

    xmm3 = _mm_load_pd(rmat);
    xmm4 = _mm_mul_pd(xmm0,xmm3);

    xmm5 = _mm_load_pd(rmat+2);
    xmm6 = _mm_mul_pd(xmm1,xmm5);

    /* calculate (b5*d1 + b6*d2 | b7*d3 + b8*d4) */
    xmm3 = _mm_hadd_pd(xmm4,xmm6);     /* we need it (2) */

    rmat += states;

    xmm4 = _mm_load_pd(rmat);
    xmm5 = _mm_mul_pd(xmm0,xmm4);

    xmm6 = _mm_load_pd(rmat+2);
    xmm7 = _mm_mul_pd(xmm1,xmm6);

    /* calculate (b9*d1 + b10*d2 | b11*d3 + b12*d4) */
    xmm4 = _mm_hadd_pd(xmm5,xmm7);     /* we need it (3) */

    rmat += states;

    xmm5 = _mm_load_pd(rmat);
    xmm6 = _mm_mul_pd(xmm0,xmm5);

    xmm7 = _mm_load_pd(rmat+2);
    xmm8 = _mm_mul_pd(xmm1,xmm7);

    /* calculate (b13*d1 + b14*d2 | b15*d3 + b16*d4) */
    xmm5 = _mm_hadd_pd(xmm6,xmm8);     /* we need it (4) */
    
    /* calculate (b1*d1 + b2*d2 + b3*d3 + b4*d3 |
                  b5*d1 + b6*d2 + b7*d3 + b8*d4 ) */
    xmm0 = _mm_hadd_pd(xmm2,xmm3);

    /* calculate (b9*d1 + b10*d2 + b11*d3 + b12*d3 |
                  b13*d1 + b14*d2 + b15*d3 + b16*d4 ) */
    xmm1 = _mm_hadd_pd(xmm4,xmm5);

    xmm2 = _mm_load_pd(lmat);
    xmm3 = _mm_load_pd(lmat+2);
    xmm4 = _mm_mul_pd(xmm0,xmm2);
    xmm5 = _mm_mul_pd(xmm1,xmm3);

    _mm_store_pd(clv,xmm4);
    _mm_store_pd(clv+2,xmm5);

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
