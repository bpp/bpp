/*
    Copyright (C) 2016-2024 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#define SHRINK 1
#define EXPAND 2

void sort2(double * t)
{
  double x;

  if (t[0] > t[1])
  {
    x = t[0];
    t[0] = t[1];
    t[1] = x;
  }
  assert(t[0] <= t[1]);
}

void sort3(double * t)
{
  double x;

  if (t[0] > t[1])
  {
    x = t[0];
    t[0] = t[1];
    t[1] = x;
  }

  if (t[1] > t[2])
  {
    x = t[1];
    t[1] = t[2];
    t[2] = x;
  }

  if (t[0] > t[1])
  {
    x = t[0];
    t[0] = t[1];
    t[1] = x;
  }

  assert(t[0] <= t[1]);
  assert(t[1] <= t[2]);
}

void debug_linked_notheta3(stree_t * stree, gtree_t ** gtree, double bpp_logPG, const char * move, int only_pg, long lmodel)
{
  double diff = 1e-5;
  double con0 = 0;
  double con1 = 0;
  double con2 = 0;
  double con3 = 0;
  double con4 = 0;
  double con5 = 0;
  double con6 = 0;
  double con7 = 0;

  double hphi_sum_0 = 0;
  double hphi_sum_1 = 0;
  double hphi_sum_2 = 0;
  double hphi_sum_3 = 0;
  double hphi_sum_4 = 0;
  double hphi_sum_5 = 0;
  double hphi_sum_6 = 0;
  double hphi_sum_7 = 0;

  double tauR = stree->root->tau;
  double tauX = stree->nodes[4]->tau;
  double tauY = stree->nodes[5]->tau;
  double tauH = stree->nodes[6]->tau;

  double phi = stree->nodes[6]->hphi;
  double a = opt_theta_alpha;
  double b = opt_theta_beta;

  double logPG = 0;

  /* contributions A and B are always 0 */
  con0 = 0;
  con1 = 0;
  
  /* locally used variables */
  double C20 = 0;
  double C21 = 0;
  double C30 = 0;
  double C31 = 0;
  double C40 = 0;
  double C41 = 0;
  double C50 = 0;
  double C51 = 0;
  double C60 = 0;
  double C61 = 0;
  double C70 = 0;
  double C71 = 0;
  double t0[3];
  double t1[3];
  long cc0 = 0; /* coalescent count for locus 0 */
  long cc1 = 0;
  gnode_t * x;

  /* contribution C (2) */
  if (stree->nodes[2]->coal_count[0])
  {
    x = (gnode_t *)(stree->nodes[2]->coalevent[0]->head->data);
    C20 = 2*x->time;
  }
  else
    C20 = 2*tauH;
  if (stree->nodes[2]->coal_count[1])
  {
    x = (gnode_t *)(stree->nodes[2]->coalevent[1]->head->data);
    C21 = 2*x->time;
  }
  else
    C21 = 2*tauH;
  cc0 = stree->nodes[2]->coal_count[0];
  cc1 = stree->nodes[2]->coal_count[1];
  con2 = a*log(b) - lgamma(a) + lgamma(a+cc0+cc1) - (a+cc0+cc1)*log(b+C20+C21);
  
  cc0 = cc1 = 0;

  /* contribution R (3) */
  assert(stree->nodes[3]->coal_count[0]);
  assert(stree->nodes[3]->coal_count[1]);
  assert(stree->nodes[3]->coal_count[0] <= 3);
  assert(stree->nodes[3]->coal_count[1] <= 3);
  
  /* locus 0 */
  x = (gnode_t *)(stree->nodes[3]->coalevent[0]->head->data);
  t0[0] = x->time;
  if (stree->nodes[3]->coal_count[0] > 1)
  {
    x = (gnode_t *)(stree->nodes[3]->coalevent[0]->head->next->data);
    t0[1] = x->time;
  }
  if (stree->nodes[3]->coal_count[0] > 2)
  {
  x = (gnode_t *)(stree->nodes[3]->coalevent[0]->head->next->next->data);
    t0[2] = x->time;
  }

  if (stree->nodes[3]->coal_count[0] == 1)
  {
    C30 = 2*(t0[0]-tauR);
  }
  else if (stree->nodes[3]->coal_count[0] == 2)
  {
    sort2(t0);
    C30 = 2*(3*(t0[0]-tauR)+t0[1]-t0[0]);
  }
  else if (stree->nodes[3]->coal_count[0] == 3)
  {
    sort3(t0);
    C30 = 2*(6*(t0[0]-tauR)+3*(t0[1]-t0[0])+(t0[2]-t0[1]));
  }
  else
    assert(0);

  /* locus 1 */
  x = (gnode_t *)(stree->nodes[3]->coalevent[1]->head->data);
  t1[0] = x->time;
  if (stree->nodes[3]->coal_count[1] > 1)
  {
    x = (gnode_t *)(stree->nodes[3]->coalevent[1]->head->next->data);
    t1[1] = x->time;
  }
  if (stree->nodes[3]->coal_count[1] > 2)
  {
    x = (gnode_t *)(stree->nodes[3]->coalevent[1]->head->next->next->data);
    t1[2] = x->time;
  }

  if (stree->nodes[3]->coal_count[1] == 1)
  {
    C31 = 2*(t1[0]-tauR);
  }
  if (stree->nodes[3]->coal_count[1] == 2)
  {
    sort2(t1);
    C31 = 2*(3*(t1[0]-tauR)+t1[1]-t1[0]);
  }
  else if (stree->nodes[3]->coal_count[1] == 3)
  {
    sort3(t1);
    C31 = 2*(6*(t1[0]-tauR)+3*(t1[1]-t1[0])+(t1[2]-t1[1]));
  }
  cc0 = stree->nodes[3]->coal_count[0];
  cc1 = stree->nodes[3]->coal_count[1];
  con3 = a*log(b) - lgamma(a) + lgamma(a+cc0+cc1) - (a+cc0+cc1)*log(b+C30+C31);

  /* contribution X (4) */
  assert(stree->nodes[4]->seqin_count[0] >= 1);
  assert(stree->nodes[4]->seqin_count[1] >= 1);
  assert(stree->nodes[4]->seqin_count[0] <= 3);
  assert(stree->nodes[4]->seqin_count[1] <= 3);

  if (stree->nodes[4]->coal_count[0] == 0)
  {
    /* no coalescence */
    cc0 = 0;
    if (stree->nodes[4]->seqin_count[0] == 1)
    {
      C40 = 0;
    }
    else if (stree->nodes[4]->seqin_count[0] == 2)
    {
      C40 = 2*(tauR-tauX);
    }
    else if (stree->nodes[4]->seqin_count[0] == 3)
    {
      C40 = 2*3*(tauR-tauX);
    }
  }
  else 
  {
    /* at least one coalescence */

    x = (gnode_t *)(stree->nodes[4]->coalevent[0]->head->data);
    t0[0] = x->time;
    if (stree->nodes[4]->coal_count[0] > 1)
    {
      /* two coal events */
      assert(stree->nodes[4]->seqin_count[0] == 3);
      x = (gnode_t *)(stree->nodes[4]->coalevent[0]->head->next->data);
      t0[1] = x->time;
      sort2(t0);
      cc0=2;
      C40 = 2*(3*(t0[0]-tauX)+t0[1]-t0[0]);
    }
    else
    {
      /* one coal event */
      cc0 = 1;
      assert(stree->nodes[4]->coal_count[0] == 1);
      if (stree->nodes[4]->seqin_count[0] == 2)
      {
        C40 = 2*(t0[0]-tauX);
      }
      else
      {
        assert(stree->nodes[4]->seqin_count[0] == 3);
        C40 = 2*(3*(t0[0]-tauX)+tauR-t0[0]);
      }
    }
  }

  if (stree->nodes[4]->coal_count[1] == 0)
  {
    /* no coalescence */
    cc1 = 0;
    if (stree->nodes[4]->seqin_count[1] == 1)
    {
      C41 = 0;
    }
    else if (stree->nodes[4]->seqin_count[1] == 2)
    {
      C41 = 2*(tauR-tauX);
    }
    else if (stree->nodes[4]->seqin_count[1] == 3)
    {
      C41 = 2*3*(tauR-tauX);
    }
  }
  else 
  {
    /* at least one coalescence */

    x = (gnode_t *)(stree->nodes[4]->coalevent[1]->head->data);
    t1[0] = x->time;
    if (stree->nodes[4]->coal_count[1] > 1)
    {
      /* two coal events */
      assert(stree->nodes[4]->seqin_count[1] == 3);
      x = (gnode_t *)(stree->nodes[4]->coalevent[1]->head->next->data);
      t1[1] = x->time;
      sort2(t1);
      cc1=2;
      C41 = 2*(3*(t1[0]-tauX)+t1[1]-t1[0]);
    }
    else
    {
      /* one coal event */
      cc1 = 1;
      assert(stree->nodes[4]->coal_count[1] == 1);
      if (stree->nodes[4]->seqin_count[1] == 2)
      {
        C41 = 2*(t1[0]-tauX);
      }
      else
      {
        assert(stree->nodes[4]->seqin_count[1] == 3);
        C41 = 2*(3*(t1[0]-tauX)+tauR-t1[0]);
      }
    }
  }

  con4 = a*log(b) - lgamma(a) + lgamma(a+cc0+cc1) - (a+cc0+cc1)*log(b+C40+C41);

  /* contributions to Y (5) */
  assert(stree->nodes[5]->seqin_count[0] >= 1);
  assert(stree->nodes[5]->seqin_count[1] >= 1);
  assert(stree->nodes[5]->seqin_count[0] <= 3);
  assert(stree->nodes[5]->seqin_count[1] <= 3);

  if (stree->nodes[5]->coal_count[0] == 0)
  {
    /* no coalescence */
    cc0 = 0;
    if (stree->nodes[5]->seqin_count[0] == 1)
    {
      C50 = 0;
    }
    else if (stree->nodes[5]->seqin_count[0] == 2)
    {
      C50 = 2*(tauR-tauY);
    }
    else if (stree->nodes[5]->seqin_count[0] == 3)
    {
      C50 = 2*3*(tauR-tauY);
    }
  }
  else 
  {
    /* at least one coalescence */

    x = (gnode_t *)(stree->nodes[5]->coalevent[0]->head->data);
    t0[0] = x->time;
    if (stree->nodes[5]->coal_count[0] > 1)
    {
      /* two coal events */
      assert(stree->nodes[5]->seqin_count[0] == 3);
      x = (gnode_t *)(stree->nodes[5]->coalevent[0]->head->next->data);
      t0[1] = x->time;
      sort2(t0);
      cc0=2;
      C50 = 2*(3*(t0[0]-tauY)+t0[1]-t0[0]);
    }
    else
    {
      /* one coal event */
      cc0 = 1;
      assert(stree->nodes[5]->coal_count[0] == 1);
      if (stree->nodes[5]->seqin_count[0] == 2)
      {
        C50 = 2*(t0[0]-tauY);
      }
      else
      {
        assert(stree->nodes[5]->seqin_count[0] == 3);
        C50 = 2*(3*(t0[0]-tauY)+tauR-t0[0]);
      }
    }
  }

  if (stree->nodes[5]->coal_count[1] == 0)
  {
    /* no coalescence */
    cc1 = 0;
    if (stree->nodes[5]->seqin_count[1] == 1)
    {
      C51 = 0;
    }
    else if (stree->nodes[5]->seqin_count[1] == 2)
    {
      C51 = 2*(tauR-tauY);
    }
    else if (stree->nodes[5]->seqin_count[1] == 3)
    {
      C51 = 2*3*(tauR-tauY);
    }
  }
  else 
  {
    /* at least one coalescence */

    x = (gnode_t *)(stree->nodes[5]->coalevent[1]->head->data);
    t1[0] = x->time;
    if (stree->nodes[5]->coal_count[1] > 1)
    {
      /* two coal events */
      assert(stree->nodes[5]->seqin_count[1] == 3);
      x = (gnode_t *)(stree->nodes[5]->coalevent[1]->head->next->data);
      t1[1] = x->time;
      sort2(t1);
      cc1=2;
      C51 = 2*(3*(t1[0]-tauY)+t1[1]-t1[0]);
    }
    else
    {
      /* one coal event */
      cc1 = 1;
      assert(stree->nodes[5]->coal_count[1] == 1);
      if (stree->nodes[5]->seqin_count[1] == 2)
      {
        C51 = 2*(t1[0]-tauY);
      }
      else
      {
        assert(stree->nodes[5]->seqin_count[1] == 3);
        C51 = 2*(3*(t1[0]-tauY)+tauR-t1[0]);
      }
    }
  }

  con5 = a*log(b) - lgamma(a) + lgamma(a+cc0+cc1) - (a+cc0+cc1)*log(b+C50+C51);

  /* contribution H (6) */
  cc0 = C60 = 0;
  if (stree->nodes[2]->coal_count[0] == 1)
  {
    assert(stree->nodes[6]->seqin_count[0] == 1 || stree->nodes[6]->seqin_count[0] == 0);
    x = (gnode_t *)(stree->nodes[2]->coalevent[0]->head->data);
    assert(x->hpath[0] == BPP_HPATH_LEFT || x->hpath[0] == BPP_HPATH_RIGHT);

    if (x->hpath[0] == BPP_HPATH_LEFT)
    {
      assert(stree->nodes[6]->seqin_count[0] == 1);
      hphi_sum_6 = log(phi);
      con6 = log(phi);
    }
  }
  else
  {
    if (stree->nodes[6]->seqin_count[0] == 0)
      hphi_sum_6 = 0;
    else if (stree->nodes[6]->seqin_count[0] == 1)
    {
      hphi_sum_6 = log(phi);
      con6 = log(phi);
    }
    else
    {
      assert(stree->nodes[6]->seqin_count[0] == 2);
      hphi_sum_6 = 2*log(phi);
      if (stree->nodes[6]->coal_count[0] == 1)
      {
        x = (gnode_t *)(stree->nodes[6]->coalevent[0]->head->data);
        t0[0] = x->time;
        cc0 = 1;
        C60 = 2*(t0[0]-tauH);
      }
      else
      {
        assert(stree->nodes[6]->coal_count[0] == 0);
        C60 = 2*(tauY-tauH);
      }
    }
  }

  cc1 = C61 = 0;
  if (stree->nodes[2]->coal_count[1] == 1)
  {
    assert(stree->nodes[6]->seqin_count[1] == 1 || stree->nodes[6]->seqin_count[1] == 0);
    x = (gnode_t *)(stree->nodes[2]->coalevent[1]->head->data);
    assert(x->hpath[0] == BPP_HPATH_LEFT || x->hpath[0] == BPP_HPATH_RIGHT);

    if (x->hpath[0] == BPP_HPATH_LEFT)
    {
      assert(stree->nodes[6]->seqin_count[1] == 1);
      hphi_sum_6 += log(phi);
      con6 = log(phi);
    }
  }
  else
  {
    if (stree->nodes[6]->seqin_count[1] == 0)
    {
      ;
    }
    else if (stree->nodes[6]->seqin_count[1] == 1)
    {
      hphi_sum_6 += log(phi);
    }
    else
    {
      hphi_sum_6 += 2*log(phi);
      assert(stree->nodes[6]->seqin_count[1] == 2);
      if (stree->nodes[6]->coal_count[1] == 1)
      {
        x = (gnode_t *)(stree->nodes[6]->coalevent[1]->head->data);
        t1[0] = x->time;
        cc1 = 1;
        C61 = 2*(t1[0]-tauH);
      }
      else
      {
        assert(stree->nodes[6]->coal_count[1] == 0);
        C61 = 2*(tauY-tauH);
      }
    }
  }
  con6 = a*log(b) - lgamma(a) + lgamma(a+cc0+cc1) - (a+cc0+cc1)*log(b+C60+C61);
  con6 += hphi_sum_6;

  /* contribution H (7) */
  cc0 = C70 = 0;
  if (stree->nodes[2]->coal_count[0] == 1)
  {
    assert(stree->nodes[7]->seqin_count[0] == 1 || stree->nodes[7]->seqin_count[0] == 0);
    x = (gnode_t *)(stree->nodes[2]->coalevent[0]->head->data);
    assert(x->hpath[0] == BPP_HPATH_LEFT || x->hpath[0] == BPP_HPATH_RIGHT);

    if (x->hpath[0] == BPP_HPATH_RIGHT)
    {
      assert(stree->nodes[7]->seqin_count[0] == 1);
      hphi_sum_7 = log(1-phi);
      con7 = log(1-phi);
    }
  }
  else
  {
    if (stree->nodes[7]->seqin_count[0] == 0)
      hphi_sum_7 = 0;
    else if (stree->nodes[7]->seqin_count[0] == 1)
    {
      hphi_sum_7 = log(1-phi);
      con7 = log(1-phi);
    }
    else
    {
      assert(stree->nodes[7]->seqin_count[0] == 2);
      hphi_sum_7 = 2*log(1-phi);
      if (stree->nodes[7]->coal_count[0] == 1)
      {
        x = (gnode_t *)(stree->nodes[7]->coalevent[0]->head->data);
        t0[0] = x->time;
        cc0 = 1;
        C70 = 2*(t0[0]-tauH);
      }
      else
      {
        assert(stree->nodes[7]->coal_count[0] == 0);
        C70 = 2*(tauX-tauH);
      }
    }
  }

  cc1 = C71 = 0;
  if (stree->nodes[2]->coal_count[1] == 1)
  {
    assert(stree->nodes[7]->seqin_count[1] == 1 || stree->nodes[7]->seqin_count[1] == 0);
    x = (gnode_t *)(stree->nodes[2]->coalevent[1]->head->data);
    assert(x->hpath[0] == BPP_HPATH_LEFT || x->hpath[0] == BPP_HPATH_RIGHT);

    if (x->hpath[0] == BPP_HPATH_RIGHT)
    {
      assert(stree->nodes[7]->seqin_count[1] == 1);
      hphi_sum_7 += log(1-phi);
      con7 = log(1-phi);
    }
  }
  else
  {
    if (stree->nodes[7]->seqin_count[1] == 0)
    {
      ;
    }
    else if (stree->nodes[7]->seqin_count[1] == 1)
    {
      hphi_sum_7 += log(1-phi);
    }
    else
    {
      hphi_sum_7 += 2*log(1-phi);
      assert(stree->nodes[7]->seqin_count[1] == 2);
      if (stree->nodes[7]->coal_count[1] == 1)
      {
        x = (gnode_t *)(stree->nodes[7]->coalevent[1]->head->data);
        t1[0] = x->time;
        cc1 = 1;
        C71 = 2*(t1[0]-tauH);
      }
      else
      {
        assert(stree->nodes[7]->coal_count[1] == 0);
        C71 = 2*(tauX-tauH);
      }
    }
  }
  con7 = a*log(b) - lgamma(a) + lgamma(a+cc0+cc1) - (a+cc0+cc1)*log(b+C70+C71);
  con7 += hphi_sum_7;

  logPG = con0+con1+con2+con3+con4+con5+con6+con7+6*log(2);



  if (!only_pg)
  {
    if (fabs(con0-stree->nodes[0]->notheta_logpr_contrib) > diff)
    {
      debug_print_gtree(gtree[0]);
      debug_print_gtree(gtree[1]);
      fatal("opt_debug_counter: %ld  move: %s con0: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[0]->notheta_logpr_contrib, con0);
    }
    if (fabs(con1-stree->nodes[1]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con1: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[1]->notheta_logpr_contrib, con1);
    if (fabs(con2-stree->nodes[2]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con2: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[2]->notheta_logpr_contrib, con2);
    if (fabs(con3-stree->nodes[3]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con3: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[3]->notheta_logpr_contrib, con3);
    if (fabs(con4-stree->nodes[4]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con4: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[4]->notheta_logpr_contrib, con4);
    if (fabs(con5-stree->nodes[5]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con5: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[5]->notheta_logpr_contrib, con5);
    if (fabs(con6-stree->nodes[6]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con6: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[6]->notheta_logpr_contrib, con6);
    if (fabs(con7-stree->nodes[7]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s con7: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[7]->notheta_logpr_contrib, con7);
  }

  if (fabs(logPG-bpp_logPG) > diff)
      fatal("opt_debug_counter: %ld  move: %s logPG: %f correct: %f\n",
            opt_debug_counter, move, bpp_logPG, logPG);
}

void debug_linked_notheta(stree_t * stree, gtree_t * gtree, double bpp_logPG, const char * move, int only_pg, long lmodel)
{
  double conA  = 0, hphi_sum_A  = 0;   /* snode 0 */
  double conB  = 0, hphi_sum_B  = 0;   /* snode 1 */
  double conR  = 0, hphi_sum_R  = 0;   /* snode 2 */
  double conX  = 0, hphi_sum_X  = 0;   /* snode 3 */
  double conY  = 0, hphi_sum_Y  = 0;   /* snode 4 */
  double conYM = 0, hphi_sum_YM = 0;   /* snode 5 */
  double phi = stree->nodes[5]->hphi;  /* horizontal branch phi */

  double tauR = stree->root->tau;
  double tauX = stree->nodes[3]->tau;
  double t = gtree->nodes[2]->time;

  double hpath_flag = gtree->nodes[1]->hpath[0];

  double a = opt_theta_alpha;
  double b = opt_theta_beta;

  double logPG = 0;
  double C2;

  double diff = 1e-5;


  if (t > tauR)
  {
    /* case 1 */
    if (hpath_flag == 1)
    {
      /* case 1: no migration of sequence b to A */
      C2 = 2*(t-tauR);

      if (lmodel == BPP_LINKEDTHETA_MSCI)
      {
        conA  = 0;
        conB  = log(1-phi);
        conR  = a*log(b) - lgamma(a) + lgamma(a+1) - (a+1)*log(b+C2);
        conX  = 0;
        conY  = 0;
        conYM = 0;

        hphi_sum_A  = 0;
        hphi_sum_B  = 0;
        hphi_sum_R  = 0;
        hphi_sum_X  = 0;
        hphi_sum_Y  = log(1-phi);
        hphi_sum_YM = 0;

        logPG = conB + conR + log(2);
      }
      else if (lmodel == BPP_LINKEDTHETA_ALL)
      {
        conA  = 0;
        conB  = 0;
        conR  = a*log(b) - lgamma(a) + lgamma(a+1) - (a+1)*log(b+C2) + log(1-phi);
        conX  = 0;
        conY  = 0;
        conYM = 0;

        hphi_sum_A  = 0;
        hphi_sum_B  = 0;
        hphi_sum_R  = 0;
        hphi_sum_X  = 0;
        hphi_sum_Y  = log(1-phi);
        hphi_sum_YM = 0;
        
        logPG = conR + log(2);
      }
      else
        assert(0);
    }
    else
    {
      assert(hpath_flag == 2);
      /* case 2: migration of sequence b to A */
      if (lmodel == BPP_LINKEDTHETA_MSCI)
      {
        C2 = 2*(t-tauR);
        conA  = a*log(b) - a*log(b+2*(tauR-tauX));
        conB  = 0;
        conR  = a*log(b) - lgamma(a) + lgamma(a+1) - (a+1)*log(b+C2);
        conX  = 0;
        conY  = 0;
        conYM = log(phi);

        hphi_sum_A  = 0;
        hphi_sum_B  = 0;
        hphi_sum_R  = 0;
        hphi_sum_X  = 0;
        hphi_sum_Y  = 0;
        hphi_sum_YM = log(phi);

        logPG = conA + conR + log(2) + conYM;
      }
      else if (lmodel == BPP_LINKEDTHETA_ALL)
      {
        C2 = 2*(t-tauR) + 2*(tauR-tauX);

        conA  = 0;
        conB  = 0;
        conR  = a*log(b) - lgamma(a) + lgamma(a+1) - (a+1)*log(b+C2) + log(phi);
        conX  = 0;
        conY  = 0;
        conYM = 0;

        hphi_sum_A  = 0;
        hphi_sum_B  = 0;
        hphi_sum_R  = 0;
        hphi_sum_X  = 0;
        hphi_sum_Y  = 0;
        hphi_sum_YM = log(phi);

        logPG = conR + log(2);
      }
      else
        assert(0);
    }
  }
  else
  {
    /* case 3: coalescence in population X */
    assert(hpath_flag == 2);

    C2 = 2*(t-tauX);

    if (lmodel == BPP_LINKEDTHETA_MSCI)
    {
      conA  = a*log(b) - lgamma(a) + lgamma(a+1) - (a+1)*log(b+C2);
      conB  = 0;
      conR  = 0;
      conX  = 0;
      conY  = 0;
      conYM = log(phi);

      hphi_sum_A  = 0;
      hphi_sum_B  = 0;
      hphi_sum_R  = 0;
      hphi_sum_X  = 0;
      hphi_sum_Y  = 0;
      hphi_sum_YM = log(phi);

      logPG = conA + log(2) + conYM;
    }
    else if (lmodel == BPP_LINKEDTHETA_ALL)
    {
      conA  = 0;
      conB  = 0;
      conR  = a*log(b) - lgamma(a) + lgamma(a+1) - (a+1)*log(b+C2) + log(phi);
      conX  = 0;
      conY  = 0;
      conYM = 0;

      hphi_sum_A  = 0;
      hphi_sum_B  = 0;
      hphi_sum_R  = 0;
      hphi_sum_X  = 0;
      hphi_sum_Y  = 0;
      hphi_sum_YM = log(phi);

      logPG = conR + log(2);
    }
    else 
      assert(0);
  }

  if (!only_pg)
  {
    if (fabs(conA-stree->nodes[0]->notheta_logpr_contrib) > diff)
    {
      debug_print_gtree(gtree);
      fatal("opt_debug_counter: %ld  move: %s conA: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[0]->notheta_logpr_contrib, conA);
    }
    if (fabs(conB-stree->nodes[1]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s conB: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[1]->notheta_logpr_contrib, conB);
    if (fabs(conR-stree->nodes[2]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s conR: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[2]->notheta_logpr_contrib, conR);
    if (fabs(conX-stree->nodes[3]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s conX: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[3]->notheta_logpr_contrib, conX);
    if (fabs(conY-stree->nodes[4]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s conY: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[4]->notheta_logpr_contrib, conY);
    if (fabs(conYM-stree->nodes[5]->notheta_logpr_contrib) > diff)
      fatal("opt_debug_counter: %ld  move: %s conYM: %f correct: %f\n",
            opt_debug_counter, move, stree->nodes[5]->notheta_logpr_contrib, conYM);
  }

  if (fabs(logPG-bpp_logPG) > diff)
      fatal("opt_debug_counter: %ld  move: %s logPG: %f correct: %f\n",
            opt_debug_counter, move, bpp_logPG, logPG);
}

void debug_print_migrations(stree_t * stree)
{
  long i;
  for (i = 0; i < opt_migration_count; ++i)
  {
    migspec_t * spec = opt_mig_specs+i;
    unsigned int si = spec->si;
    unsigned int ti = spec->ti;
    printf("  M_%s->%s=%f (%d,%d,%f,%f)",
           stree->nodes[si]->label,
           stree->nodes[ti]->label,
           spec->M,
           spec->si,
           spec->ti,
           spec->alpha,
           spec->beta);
  }
  printf("\n");
}

void debug_print_migmatrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;

  printf("      ");
  for (i = 0; i < total_nodes; ++i)
    printf(" %2s", stree->nodes[i]->label);
  printf("\n");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("   ");
    printf(" %2s", stree->nodes[i]->label);
    for (j = 0; j < total_nodes; ++j)
    {
      if (opt_migration_matrix[i][j] != -1)
        printf(ANSI_COLOR_RED " %2ld" ANSI_COLOR_RESET, opt_migration_matrix[i][j]);
      else
        printf(" %2ld", opt_migration_matrix[i][j]);
    }
    printf("\n");
  }
}


void debug_print_bitmatrix(stree_t * stree)
{
  long i,j;
  long total_nodes = stree->tip_count + stree->inner_count;

  printf("      ");
  for (i = 0; i < total_nodes; ++i)
    printf(" %2s", stree->nodes[i]->label);
  printf("\n");
  for (i = 0; i < total_nodes; ++i)
  {
    printf("   ");
    printf(" %2s", stree->nodes[i]->label);
    for (j = 0; j < total_nodes; ++j)
    {
      if (opt_mig_bitmatrix[i][j] != 0)
        printf(ANSI_COLOR_RED " %2ld" ANSI_COLOR_RESET, opt_mig_bitmatrix[i][j]);
      else
        printf(" %2ld", opt_mig_bitmatrix[i][j]);
    }
    printf("\n");
  }
}

void debug_im_check_sum(stree_t * stree, gtree_t ** gtree)
{
  long i,j,k;
  long migsum;

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      if (opt_mig_bitmatrix[i][j])
      {
        migsum = 0;
        for (k = 0; k < opt_locus_count; ++k)
          migsum += gtree[k]->migcount[i][j];

        if (migsum != stree->migcount_sum[i][j])
        {
          fprintf(stderr, "[DEBUG] Error in function %s\n", __FUNCTION__);
          fatal("migsum_count[%ld][%ld]: %ld   correct: %ld   locus: %ld\n",
                stree->migcount_sum[i][j], migsum, k);
          for (k = 0; k < opt_locus_count; ++k)
            debug_print_gtree(gtree[k]);
        }
        assert(migsum == stree->migcount_sum[i][j]);
      }
}

void debug_print_stree(stree_t * stree)
{
  long i;
  long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  printf("Total nodes: %ld\n", total_nodes);

  printf("Label        Node  Child1  Child2  Parent  Linked     Theta       Tau  Leaves  nin[0]  "
         "prop_tau      hphi  htau  has_phi\n");
  for (i = 0; i <total_nodes; ++i)
  {
    char * label;
    if (stree->nodes[i]->label)
      label = xstrdup(stree->nodes[i]->label);
    else
      label = xstrdup("N/A");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    printf("%-12s", label);
    free(label);
    printf(" %4d", stree->nodes[i]->node_index);

    if (stree->nodes[i]->left)
      printf("  %6d", stree->nodes[i]->left->node_index);
    else
      printf("  %6s", "N/A");
    if (stree->nodes[i]->right)
      printf("  %6d", stree->nodes[i]->right->node_index);
    else
      printf("  %6s", "N/A");
    if (stree->nodes[i]->parent)
      printf("  %6d", stree->nodes[i]->parent->node_index);
    else
      printf("  %6s", "N/A");

    if (stree->nodes[i]->linked_theta)
      printf("  %1s:%2d", "Yes", stree->nodes[i]->linked_theta->node_index);
    else
      printf("  %6s", "No");

    double theta = 0;
    if (stree->nodes[i]->linked_theta)
      theta = stree->nodes[i]->linked_theta->theta;
    else
      theta = stree->nodes[i]->theta;

    if (theta < 0)
      printf("  %8s", "-");
    else
      printf("  %.6f", theta);

    if (stree->nodes[i]->tau)
      printf("  %.6f", stree->nodes[i]->tau);
    else
      printf("  %8s", "-");

    printf("  %6d", stree->nodes[i]->leaves);
    printf("  %6d", stree->nodes[i]->seqin_count[0]);
    printf("  %8d", stree->nodes[i]->prop_tau);
    
    if (opt_msci && stree->nodes[i]->hybrid)
      printf("  %.6f", stree->nodes[i]->hphi);
    else
      printf("  %8s", "-");

    printf("  %4ld", stree->nodes[i]->htau);

    if (opt_msci)
      printf("  %7ld\n", stree->nodes[i]->has_phi);
    else
      printf("\n");
  }
}

void debug_print_gtree(gtree_t * gtree)
{
  long i,j;

  printf("Label        Node-index  Child1-index  Child2-index  Parent-index    Mi  "
         "Pmat-index    Pop     time\n");
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    char * label;
    if (gtree->nodes[i]->label)
      label = xstrdup(gtree->nodes[i]->label);
    else
      label = xstrdup("N/A");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    printf("%-12s", label);
    free(label);
    printf("  %9d", gtree->nodes[i]->node_index);

    if (gtree->nodes[i]->left)
      printf("  %12d", gtree->nodes[i]->left->node_index);
    else
      printf("  %12s", "N/A");

    if (gtree->nodes[i]->right)
      printf("  %12d", gtree->nodes[i]->right->node_index);
    else
      printf("  %12s", "N/A");


    if (gtree->nodes[i]->parent)
      printf("  %12d", gtree->nodes[i]->parent->node_index);
    else
      printf("  %12s", "N/A");

    if (gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)
      printf("  %4ld", gtree->nodes[i]->mi->count);
    else
      printf("  %4s", "-");

    printf("  %10d", gtree->nodes[i]->pmatrix_index);

    printf(" %5s", gtree->nodes[i]->pop->label);
    printf(" %.12f", gtree->nodes[i]->time);

    if (gtree->root == gtree->nodes[i])
      printf(" (root)");
    printf("\n");
  }
  if (opt_migration)
  {
    long migfound = 0;
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;
      migfound++;
    }

    if (migfound)
    {
      printf("\nMigration events:\n");
      printf("Label  Node  #Migs  Migrations (backwards in time)\n");
    }

    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;

      gnode_t * x = gtree->nodes[i];

      char * label;
      if (x->label)
        label = xstrdup(x->label);
      else
        label = xstrdup("N/A");

      /* shorten label if over 12 characters */
      if (strlen(label) > 6)
      {
        label[4] = '.'; label[5] = '.'; label[6] = 0;
      }

      printf("%-6s", label);
      free(label);
      printf(" %4d", gtree->nodes[i]->node_index);
      printf("  %5ld", x->mi->count);

      printf("  %s -> %s   %.15f\n",
             x->mi->me[0].source->label,
             x->mi->me[0].target->label,
             x->mi->me[0].time);
      for (j = 1; j < x->mi->count; ++j)
        printf("                    %s -> %s   %f\n",
               x->mi->me[j].source->label,
               x->mi->me[j].target->label,
               x->mi->me[j].time);



    }
  }
}

void debug_print_gtree_detailed(gtree_t * gtree)
{
  long i,j;

  printf("Label        Node-index  Child1-index  Child2-index  Parent-index    Mi  "
         "Pmat-index    Pop            time   Flags\n");
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    char * label;
    if (gtree->nodes[i]->label)
      label = xstrdup(gtree->nodes[i]->label);
    else
      label = xstrdup("  -");

    /* shorten label if over 12 characters */
    if (strlen(label) > 12)
    {
      label[9] = '.'; label[10] = '.'; label[11] = '.'; label[12] = 0;
    }

    if (gtree->nodes[i]->mark & FLAG_PRUNED)
      printf(ANSI_COLOR_RED);
    printf("%-12s", label);
    free(label);
    char * str_nodeindex = NULL;
    xasprintf(&str_nodeindex,
              gtree->nodes[i] == gtree->root ? "* %d" : "%d",
              gtree->nodes[i]->node_index);
    printf("  %9s", str_nodeindex);
    free(str_nodeindex);

    if (gtree->nodes[i]->left)
      printf("  %12d", gtree->nodes[i]->left->node_index);
    else
      printf("  %12s", "  -");

    if (gtree->nodes[i]->right)
      printf("  %12d", gtree->nodes[i]->right->node_index);
    else
      printf("  %12s", "  -");


    if (gtree->nodes[i]->parent)
      printf("  %12d", gtree->nodes[i]->parent->node_index);
    else
      printf("  %12s", "  -");

    if (gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)
      printf("  %4ld", gtree->nodes[i]->mi->count);
    else
      printf("  %4s", "-");

    printf("  %10d", gtree->nodes[i]->pmatrix_index);

    if (gtree->nodes[i]->pop)
      printf("  %5s", gtree->nodes[i]->pop->label);
    else
      printf("  %5s", "-");

    printf("  %.12f", gtree->nodes[i]->time);
    if (gtree->nodes[i]->mark)
      printf("   ");

    int comma = 0;
    if (gtree->nodes[i]->mark & FLAG_AGE_UPDATE)
    {
      if (comma)
        printf(",");
      printf("AGE_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_POP_UPDATE)
    {
      if (comma)
        printf(",");
      printf("POP_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_BRANCH_UPDATE)
    {
      if (comma)
        printf(",");
      printf("BRANCH_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_MISC)
    {
      if (comma)
        printf(",");
      printf("MISC");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_PRUNED)
    {
      if (comma)
        printf(",");
      printf("PRUNED");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_MIGRATE)
    {
      if (comma)
        printf(",");
      printf("MIGRATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_PARTIAL_UPDATE)
    {
      if (comma)
        printf(",");
      printf("PARTIAL_UPDATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_RED_LEFT)
    {
      if (comma)
        printf(",");
      printf("RED_LEFT");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_RED_RIGHT)
    {
      if (comma)
        printf(",");
      printf("RED_RIGHT");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_SIMULATE)
    {
      if (comma)
        printf(",");
      printf("SIMULATE");
      comma = 1;
    }
    if (gtree->nodes[i]->mark & FLAG_PRUNED)
      printf(ANSI_COLOR_RESET);

    printf("\n");
  }
  if (opt_migration)
  {
    long migfound = 0;
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;
      migfound++;
    }

    if (migfound)
    {
      printf("\nMigration events:\n");
      printf("Label  Node  #Migs  Migrations (backwards in time)\n");
    }

    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      if (!(gtree->nodes[i]->mi && gtree->nodes[i]->mi->count)) continue;

      gnode_t * x = gtree->nodes[i];

      char * label;
      if (x->label)
        label = xstrdup(x->label);
      else
        label = xstrdup("N/A");

      /* shorten label if over 12 characters */
      if (strlen(label) > 6)
      {
        label[4] = '.'; label[5] = '.'; label[6] = 0;
      }

      printf("%-6s", label);
      free(label);
      printf(" %4d", gtree->nodes[i]->node_index);
      printf("  %5ld", x->mi->count);

      printf("  %s -> %s   %f\n", x->mi->me[0].source->label, x->mi->me[0].target->label, x->mi->me[0].time);
      for (j = 1; j < x->mi->count; ++j)
        printf("                    %s -> %s   %f\n", x->mi->me[j].source->label, x->mi->me[j].target->label, x->mi->me[j].time);



    }
  }
}

void debug_check_relations(stree_t * stree, gtree_t * gtree, long msa_index)
{
  long i,j;

  assert(gtree->root && !gtree->root->parent);

  /* marks */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    if (!gtree->nodes[i]->parent)
      assert(gtree->nodes[i] == gtree->root);
    assert(!gtree->nodes[i]->mark);
  }

  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    while (x && !x->mark)
    {
      x->mark = 1;
      
      if (x->parent)
      {
        assert(x->parent->left && x->parent->right);
        assert(x->parent->left == x || x->parent->right == x);
      }
      x = x->parent;
    }
  }

  /* marks */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    assert(gtree->nodes[i]->mark);
    gtree->nodes[i]->mark = 0;
  }

  /* check leaves */
  for (i = 0; i < gtree->tip_count; ++i)
    assert(gtree->nodes[i]->leaves == 1);
  for (i = gtree->tip_count; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];
    assert(x->leaves == x->left->leaves + x->right->leaves);
  }

  /* check seqin count, coalescences etc */
  long sum = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    sum += stree->nodes[i]->seqin_count[msa_index];
  }
  assert(sum == gtree->tip_count);
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * snode = stree->nodes[i];

    long lout = snode->left->seqin_count[msa_index] -
                snode->left->coal_count[msa_index];
    long rout = snode->right->seqin_count[msa_index] -
                snode->right->coal_count[msa_index];

    if (opt_migration)
    {
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      {
        lout += gtree->migcount[snode->left->node_index][j];
        rout += gtree->migcount[snode->right->node_index][j];

        lout -= gtree->migcount[j][snode->left->node_index];
        rout -= gtree->migcount[j][snode->right->node_index];
      }
    }

    assert(snode->seqin_count[msa_index] == lout+rout);
  }
  snode_t * sroot = stree->root;
  assert(sroot->seqin_count[msa_index] - sroot->coal_count[msa_index] == 1);

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode_t * x = stree->nodes[i];

    for (j = 1; j < x->mb_count; ++j)
      assert(x->migbuffer[j].time != x->migbuffer[j-1].time);
  }

  /* check ages */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    if (!x->parent)
    {
      assert(x->left->time < x->time && x->right->time < x->time);
    }
    else
    {
      assert(x->parent->time > x->time);

      if (x->left)
      {
        assert(x->left->time < x->time && x->right->time < x->time);
      }
    }

    if (x->pop->parent)
      assert(x->time < x->pop->parent->tau);
    if (x->left)
      assert(x->time > x->pop->tau);
  }

  /* check for duplicate pmatrix index */
  /* TODO: Also check for the corresponding index, i.e. pmatindex+edge_count */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
  {
    for (j = i+1; j < gtree->tip_count+gtree->inner_count; ++j)
    {
      if (gtree->nodes[i]->parent && gtree->nodes[j]->parent)
        assert(gtree->nodes[i]->pmatrix_index != gtree->nodes[j]->pmatrix_index);
    }
  }
}

void debug_check_leaves(gtree_t ** gtree)
{
  long i,j;

  for (i = 0; i < opt_locus_count; ++i)
  {
    for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      gnode_t * node = gtree[i]->nodes[j];

      if (node->left)
      {
        /* inner */
        assert(node->leaves == node->left->leaves + node->right->leaves);
      }
      else
      {
        /* tip */
        assert(node->leaves == 1);
      }
    }
  }
}
static void debug_validate_logpg_theta(stree_t * stree,
                                       gtree_t ** gtree,
                                       locus_t ** locus,
                                       const char * move)
{
  long i,j;
  double logpr, old_logpr;
  const long thread_index = 0;

  logpr = 0;

  for (i = 0; i < opt_locus_count; ++i)
  {
    old_logpr = gtree[i]->logpr;
    if (opt_migration)
      logpr = gtree_logprob_mig(stree,gtree[i],locus[i]->heredity[0],i,thread_index);
    else
      logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index);

    for (j = 0; j < stree->tip_count+stree->inner_count+stree->hybrid_count; ++j)
    {
      snode_t * x = stree->nodes[j];

      if (fabs(x->logpr_contrib[i] - x->old_logpr_contrib[i]) > 1e-5)
        fatal("[FATAL-%ld] move: %s locus: %ld snode: %d (%s)"
              "logpr_contrib: %f old_contrib: %f\n",
              opt_debug_counter, move, i, j, stree->nodes[j]->label,
              x->logpr_contrib[i], x->old_logpr_contrib[i]);
      if (fabs(x->C2ji[i] - x->old_C2ji[i]) > 1e-9)
        fatal("[FATAL-%ld] move: %s locus: %ld snode: %d "
              "C2ji: %f old_C2ji: %f\n",
              opt_debug_counter, move, i, j,
              x->C2ji[i], x->old_C2ji[i]);
    }

    if (fabs(logpr - old_logpr) > 1e-5)
      fatal("FATAL-%ld] move: %s locus %ld logpr: %.15f old_logpr: %.15f\n",
            opt_debug_counter, move, i, logpr, old_logpr);
  }
}

static void debug_validate_logpg_notheta(stree_t * stree,
                                         locus_t ** locus,
                                         const char * move)
{
  size_t hnodes_total;
  long i,j;
  double logpr, old_logpr;
  const long thread_index = 0;

  hnodes_total = (size_t)(2*stree->hybrid_count);
  old_logpr = stree->notheta_logpr;
  logpr = 0;

  /* temporary storage */
  double * phi_contrib = (double *)xmalloc(hnodes_total *
                                           (size_t)opt_locus_count *
                                           sizeof(double));
  double * old_phi_contrib = (double *)xmalloc(hnodes_total *
                                               (size_t)opt_locus_count *
                                               sizeof(double));
  double * phi_sum = (double *)xmalloc(hnodes_total * sizeof(double));

  double * phi_tmp = phi_contrib;
  double * old_phi_tmp = old_phi_contrib;
  double * phi_sum_tmp = phi_sum;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+i];

    /* copy old values from mirror hybrid node */
    memcpy(phi_tmp, x->notheta_phi_contrib, opt_locus_count * sizeof(double));
    memcpy(old_phi_tmp,
           x->notheta_old_phi_contrib,
           opt_locus_count * sizeof(double));
    *phi_sum_tmp = x->hphi_sum;

    /* zero-out the arrays */
    memset(x->notheta_phi_contrib,0,opt_locus_count*sizeof(double));
    memset(x->notheta_old_phi_contrib,0,opt_locus_count*sizeof(double));
    x->hphi_sum = 0;

    /* move storage pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;

    x = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;

    /* copy old values from hybrid node */
    memcpy(phi_tmp, x->notheta_phi_contrib, opt_locus_count * sizeof(double));
    memcpy(old_phi_tmp,
           x->notheta_old_phi_contrib,
           opt_locus_count * sizeof(double));
    *phi_sum_tmp = x->hphi_sum;

    /* zero-out the arrays */
    memset(x->notheta_phi_contrib,0,opt_locus_count*sizeof(double));
    memset(x->notheta_old_phi_contrib,0,opt_locus_count*sizeof(double));
    x->hphi_sum = 0;

    /* move storage pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;
  }


  /* recompute population contributions */
  for (i = 0; i < opt_locus_count; ++i)
    logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index);

  /* add constant part */
  logpr += stree->notheta_hfactor+stree->notheta_sfactor;

  phi_sum_tmp = phi_sum;
  phi_tmp = phi_contrib;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    /* mirror hybrid node */
    snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+i];

    /* validate phi contributions */
    for (j = 0; j < opt_locus_count; ++j)
      if (fabs(x->notheta_phi_contrib[j] - phi_tmp[j]) > 1e-9)
        fatal("[FATAL-%ld] move: %s hybrid: %d mirror: YES    "
              "phi_contrib[%ld]: %f  old_contrib[%ld]: %f    "
              "seqin_count: %d phi: %f\n",
              opt_debug_counter,
              move,
              i,
              j,
              x->notheta_phi_contrib[j],
              j,
              phi_tmp[j],
              x->seqin_count[j],
              x->hphi);


    if (fabs(x->hphi_sum - *phi_sum_tmp) > 1e-9)
    {
      fatal("[FATAL-%ld] move: %s hybrid: %d mirror: YES    "
            "hphi_sum (correct): %f  old hphi_sum (wrong): %f\n",
            opt_debug_counter,
            move,
            x->hphi_sum,
            *phi_sum_tmp);
    }
    ++phi_sum_tmp;
    phi_tmp += opt_locus_count;

    /* normal hybrid node */
    x = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;

    /* validate phi contributions */
    for (j = 0; j < opt_locus_count; ++j)
      if (fabs(x->notheta_phi_contrib[j] - phi_tmp[j]) > 1e-9)
        fatal("[FATAL-%ld] move: %s hybrid: %d mirror: NO    "
              "phi_contrib[%ld]: %f  old_contrib[%ld]: %f    "
              "seqin_count: %d phi: %f\n",
              opt_debug_counter,
              move,
              i,
              j,
              x->notheta_phi_contrib[j],
              j,
              phi_tmp[j],
              x->seqin_count[j],
              x->hphi);

    if (fabs(x->hphi_sum - *phi_sum_tmp) > 1e-9)
    {
      fatal("[FATAL-%ld] move: %s hybrid: %d mirror: YES    "
            "hphi_sum (correct): %f  old hphi_sum (wrong): %f\n",
            opt_debug_counter,
            move,
            x->hphi_sum,
            *phi_sum_tmp);
    }
    ++phi_sum_tmp;
    phi_tmp += opt_locus_count;
  }

  if (fabs(logpr - old_logpr) > 1e-5)
    fatal("[FATAL-%ld] move: %s logpr (correct): %f  old logpr (wrong): %f\n",
          opt_debug_counter, move, logpr, old_logpr);

  /* restore */
  phi_tmp = phi_contrib;
  old_phi_tmp = old_phi_contrib;
  phi_sum_tmp = phi_sum;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    /* normal hybrid node */
    snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+i];

    /* restore elements */
    memcpy(x->notheta_phi_contrib, phi_tmp, opt_locus_count * sizeof(double));
    memcpy(x->notheta_old_phi_contrib,
           old_phi_tmp,
           opt_locus_count * sizeof(double));
    x->hphi_sum = *phi_sum_tmp;

    /* move pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;

    /* mirror hybrid node */
    x = stree->nodes[stree->tip_count+stree->inner_count+i]->hybrid;

    /* restore elements */
    memcpy(x->notheta_phi_contrib, phi_tmp, opt_locus_count * sizeof(double));
    memcpy(x->notheta_old_phi_contrib,
           old_phi_tmp,
           opt_locus_count * sizeof(double));
    x->hphi_sum = *phi_sum_tmp;

    /* move pointers */
    phi_tmp += opt_locus_count;
    old_phi_tmp += opt_locus_count;
    phi_sum_tmp++;
  }

  free(phi_contrib);
  free(old_phi_contrib);
  free(phi_sum);
}

void debug_validate_logpg(stree_t * stree,
                          gtree_t ** gtree,
                          locus_t ** locus,
                          const char * move)
{
  if (opt_est_theta)
    debug_validate_logpg_theta(stree,gtree,locus,move);
  else
  {
    #if 0
    debug_validate_logpg_notheta(stree,locus,move);
    #else
    long i,j;
    double logpr_sum = 0;
    double * old_hphi_sum;
    double * logpr_contrib;
    double ** phi_contrib;
    double ** C2j;
    long total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

    double notheta_old_logpr = stree->notheta_logpr;
    old_hphi_sum = (double *)xmalloc((size_t)total_nodes * sizeof(double));
    logpr_contrib = (double *)xmalloc((size_t)total_nodes * sizeof(double));
    C2j = (double **)xmalloc((size_t)total_nodes * sizeof(double *));
    phi_contrib = (double **)xmalloc((size_t)total_nodes * sizeof(double *));

    stree->notheta_logpr = stree->notheta_sfactor+stree->notheta_hfactor;

    for (i = 0; i < total_nodes; ++i)
    {
      C2j[i] = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
      phi_contrib[i] = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
      snode_t * x = stree->nodes[i];
      /* save C2j */
      memcpy(C2j[i],
             x->C2ji,
             opt_locus_count*sizeof(double));

      /* save old phi contrib */
      if (opt_msci)
        memcpy(phi_contrib[i],
               x->notheta_phi_contrib,
               opt_locus_count*sizeof(double));

      /* save contribs */
      logpr_contrib[i] = x->notheta_logpr_contrib;

      old_hphi_sum[i] = x->hphi_sum;
    }

    for (j = 0; j < total_nodes; ++j)
      for (i = 0; i < opt_locus_count; ++i)
        gtree_update_C2j(stree->nodes[j], locus[i]->heredity[0], i,0);

    for (j = 0; j < total_nodes; ++j)
    {
      if (!stree->nodes[j]->linked_theta || stree->nodes[j]->hybrid)
        logpr_sum += update_logpg_contrib(stree,stree->nodes[j]);
    }
    stree->notheta_logpr += logpr_sum;
    stree->notheta_old_logpr = 0;

    /* validate */

    /* logpr */
    if (fabs(stree->notheta_logpr - notheta_old_logpr) > 1e-9)
    {
      fatal("[FATAL-%ld] move: %s "
            "notheta_logpr (correct): %f notheta_old_logpr (wrong): %f\n",
            opt_debug_counter,
            move,
            stree->notheta_logpr,
            notheta_old_logpr);
    }
    for (i = 0; i < total_nodes; ++i)
    {
      snode_t * x = stree->nodes[i];

      /* hphi_sum */
      if (fabs(x->hphi_sum - old_hphi_sum[i]) > 1e-9)
      {
        fatal("[FATAL-%ld] move: %s "
              "hphi_sum (correct): %f  old hphi_sum (wrong): %f\n",
              opt_debug_counter,
              move,
              x->hphi_sum,
              old_hphi_sum[i]);
      }

      if (fabs(x->notheta_logpr_contrib - logpr_contrib[i]) > 1e-9)
      {
        fatal("[FATAL-%ld] move: %s "
              "notheta_logpr_contrib (correct): %f  old contrib (wrong): %f\n",
              opt_debug_counter,
              move,
              x->notheta_logpr_contrib,
              logpr_contrib[i]);
      }

      for (j = 0; j < opt_locus_count; ++j)
      {
        /* C2ji */
        if (fabs(x->C2ji[j] - C2j[i][j]) > 1e-9)
        {
          fatal("[FATAL-%ld] move: %s "
                "C2ji[%ld] (correct): %f  old_C2ji[%ld] (wrong): %f\n",
                opt_debug_counter,
                move,
                j,
                x->C2ji[j],
                j,
                C2j[j]);
        }

        /* phi contrib */
        if (opt_msci)
          if (fabs(x->notheta_phi_contrib[j] - phi_contrib[i][j]) > 1e-9)
          {
            fatal("[FATAL-%ld] move: %s "
                  "notheta_phi_contrib[%ld] (correct): %f  old_phi_contrib[%ld] (wrong): %f\n",
                  opt_debug_counter,
                  move,
                  j,
                  x->notheta_phi_contrib[j],
                  j,
                  x->notheta_old_phi_contrib[j]);
          }
      }
    }



    free(old_hphi_sum);
    for (i = 0; i < total_nodes; ++i)
    {
      free(phi_contrib[i]);
      free(C2j[i]);
    }
    free(phi_contrib);
    free(C2j);
    
    #endif
  }
}

void debug_snl_stage1(stree_t * stree,
                      gtree_t ** gtree_list,
                      snode_t * y,
                      snode_t * c,
                      snode_t * a,
                      long movetype,
                      int downwards,
                      double tauy,
                      double newtauy)
{
  long i;
  double taufactor;
  char * gtree_newick = NULL;
  char * stree_newick = NULL;
  char * (*cb_stree_branch)(const snode_t *) = NULL;
  char * (*cb_gtree_branch)(const gnode_t *) = NULL;
  char * stree_desc = NULL;
  char * gtree_desc = NULL;

  if (opt_debug_counter < opt_debug_start || opt_debug_counter > opt_debug_end)
    return;

  if (opt_debug_snl > 1)
  {
    fprintf(stdout,
            "\n[DBG-%ld] ------ Proposal: SNL Number: %ld ------\n",
            opt_debug_snl,
            opt_debug_counter);

    switch (opt_debug_snl)
    {
      case 2:
        cb_stree_branch = cb_stree_print_none;
        cb_gtree_branch = cb_gtree_print_none;
        stree_desc = (char *)xmalloc(sizeof(char)); stree_desc[0] = '\0';
        gtree_desc = (char *)xmalloc(sizeof(char)); gtree_desc[0] = '\0';
        break;
      case 3:
        cb_stree_branch = cb_stree_print_tau;
        cb_gtree_branch = cb_gtree_print_age;
        stree_desc = xstrdup(" (with taus)");
        gtree_desc = xstrdup(" (with ages)");
        break;
      case 4:
      default:
        cb_stree_branch = cb_stree_print_node_tau;
        cb_gtree_branch = cb_gtree_print_node_age;
        stree_desc = xstrdup(" (with node indices and taus)");
        gtree_desc = xstrdup(" (with node indices and ages)");
        break;
    }

    stree_newick = stree_export_newick(stree->root,cb_stree_branch);
    fprintf(stdout,
            "[DBG-%ld] Species tree before move%s: %s\n",
            opt_debug_snl, stree_desc, stree_newick);
    free(stree_newick);
    free(stree_desc);

    taufactor = newtauy / tauy;
    fprintf(stdout, "[DBG-%ld] Move: %s  |  nodes Y: %s C: %s A: %s  |  "
            "tau_Y: %f tau*_Y: %f tau_Y/tau*_Y: %f\n",
            opt_debug_counter,
            movetype == EXPAND && !downwards ? 
              "EXPAND" : movetype == SHRINK ? "SHRINK" : "EXPAND+SHRINK",
            y->label,
            c->label,
            a->label,
            tauy,
            newtauy,
            taufactor);

    fprintf(stdout,
            "[DBG-%ld] Gene trees before proposal%s:\n",
            opt_debug_snl, gtree_desc);
    free(gtree_desc);
    for (i = 0; i < opt_locus_count; ++i)
    {
      gtree_newick = gtree_export_newick(gtree_list[i]->root,
                                               cb_gtree_branch);
      printf("[DBG-%ld]  %ld  ->  %s\n", opt_debug_snl, i+1, gtree_newick);
      free(gtree_newick);
    }
  }
}

void debug_gage_stage2(stree_t * stree,
                       gtree_t ** gtree,
                       double logpr_notheta,
                       double r,
                       double lnacceptance)
{
  if (opt_debug_gage == 1)
  {
  }
}

void debug_snl_stage2(stree_t * stree,
                      gtree_t ** gtree,
                      double logpr_notheta,
                      double r,
                      double lnacceptance)
{
  long i;
  char * gtree_newick = NULL;
  char * stree_newick = NULL;
  char * (*cb_stree_branch)(const snode_t *) = NULL;
  char * (*cb_gtree_branch)(const gnode_t *) = NULL;
  char * stree_desc = NULL;
  char * gtree_desc = NULL;

  if (opt_debug_counter < opt_debug_start || opt_debug_counter > opt_debug_end)
    return;

  if (opt_debug_snl == 1)
  {
    long rc = (lnacceptance >= -1e-10 || r < exp(lnacceptance));
    fprintf(stdout,
            "[DBG-%ld] Proposal: SNL Number: %ld  lnacceptance: %f  |  %s\n",
            opt_debug_snl, opt_debug_counter, lnacceptance,
            rc ? "ACCEPTED" : "REJECTED");
    return;
  }

  if (opt_debug_snl > 1)
  {
    switch (opt_debug_snl)
    {
      case 2:
        cb_stree_branch = cb_stree_print_none;
        cb_gtree_branch = cb_gtree_print_none;
        stree_desc = (char *)xmalloc(sizeof(char)); stree_desc[0] = '\0';
        gtree_desc = (char *)xmalloc(sizeof(char)); gtree_desc[0] = '\0';
        break;
      case 3:
        cb_stree_branch = cb_stree_print_tau;
        cb_gtree_branch = cb_gtree_print_age;
        stree_desc = xstrdup(" (with taus)");
        gtree_desc = xstrdup(" (with ages)");
        break;
      case 4:
      default:
        cb_stree_branch = cb_stree_print_node_tau;
        cb_gtree_branch = cb_gtree_print_node_age;
        stree_desc = xstrdup(" (with node indices and taus)");
        gtree_desc = xstrdup(" (with node indices and ages)");
        break;
    }

    stree_newick = stree_export_newick(stree->root,cb_stree_branch);
    fprintf(stdout,
            "[DBG-%ld] Proposed species tree%s: %s\n",
            opt_debug_snl, stree_desc, stree_newick);
    free(stree_newick);
    free(stree_desc);

    fprintf(stdout,
            "[DBG-%ld] Proposed gene trees%s:\n",
            opt_debug_snl, gtree_desc);
    free(gtree_desc);
    for (i = 0; i < opt_locus_count; ++i)
    {
      gtree_newick = gtree_export_newick(gtree[i]->root,cb_gtree_branch);
      fprintf(stdout,
              "[DBG-%ld]  %ld  ->  %s\n",
              opt_debug_snl, i+1, gtree_newick);
      free(gtree_newick);
    }

    fprintf(stdout,
            "[DBG-%ld] log-L and gene tree probabilities:\n",
            opt_debug_snl);
    for (i = 0; i < opt_locus_count; ++i)
    {
      if (opt_est_theta)
        fprintf(stdout,
                "[DBG-%ld]   %ld  ->  log-L: %f  old log-L: %f  log-L diff: %f   "
                "|   log-gP: %f  old log-gP: %f  log-gP diff: %f\n",
                opt_debug_snl, i+1,
                gtree[i]->logl, gtree[i]->old_logl,
                gtree[i]->logl - gtree[i]->old_logl,
                gtree[i]->logpr, gtree[i]->old_logpr,
                gtree[i]->logpr - gtree[i]->old_logpr);
      else
        fprintf(stdout,
                "[DBG-%ld]   %ld  ->  log-L: %f  old log-L: %f  log-L diff: %f\n",
                opt_debug_snl, i+1,
                gtree[i]->logl, gtree[i]->old_logl,
                gtree[i]->logl - gtree[i]->old_logl);
    }
    if (!opt_est_theta)
    {
      fprintf(stdout,
              "[DBG-%ld] log-gP: %f  old log-gP: %f  log-gP diff: %f\n",
              opt_debug_snl, logpr_notheta, stree->notheta_logpr,
              logpr_notheta - stree->notheta_logpr);
    }

    fprintf(stdout,
            "[DBG-%ld] SNL log-acceptance ratio (lnacceptance): %f\n",
            opt_debug_snl, lnacceptance);
    long rc = (lnacceptance >= -1e-10 || r < exp(lnacceptance));
    if (rc)
      fprintf(stdout, "[DBG-%ld] Proposal accepted", opt_debug_snl);
    else
      fprintf(stdout, "[DBG-%ld] Proposal rejected", opt_debug_snl);
    fprintf(stdout,
           "\n[DBG-%ld] ---------------------------------------------\n",
           opt_debug_snl);
  }
}

/* callback function to print the population label for each inner gtree node */
char * cb_gtree_print_pop(const gnode_t * node)
{
  char * s = NULL;
  if (node->left)      /* has a child <=> inner node */
  {
    if (node->pop->label)
      xasprintf(&s, "[%s]", node->pop->label);
    else
      xasprintf(&s, "[%d]", node->pop->node_index);
  }
  else
  {
    /* tip node */
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

/* alternative callback function to print the population label AND age for each
   inner gtree node */
char * cb_gtree_print_pop_and_age(const gnode_t * node)
{
  char * s = NULL;
  if (node->left)      /* has a child <=> inner node */
  {
    if (node->pop->label)
      xasprintf(&s, "[%s]:%f", node->pop->label, node->time);
    else
      xasprintf(&s, "[%d]:%f", node->pop->node_index, node->time);
  }
  else
  {
    /* tip node */
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

void debug_bruce(stree_t * stree,
                 gtree_t ** gtree,
                 const char * move,
                 long iter,
                 FILE * fp_out)
{
  long i;
  char * newick;

  fprintf(fp_out, "%ld - %s - Gene trees:\n", iter, move);
  for (i = 0; i < opt_locus_count; ++i)
  {
    newick = gtree_export_newick(gtree[i]->root,cb_gtree_print_pop);
    fprintf(fp_out, "  %s\n", newick);
    free(newick);
  }
}

void debug_migration_internals(stree_t * stree, gtree_t * gtree, long msa_index)
{
  long i,j;

  long * migevent_count;
  int * seqin_count;
  long ** migcount;
  snode_t * curpop;
  snode_t * pop;

  /* return; */

  if (!opt_migration) return;

  assert(!gtree->root->parent);

  long total_nodes = stree->tip_count+stree->inner_count;

  migevent_count = (long *)xcalloc((size_t)total_nodes, sizeof(long));
  seqin_count = (int *)xcalloc((size_t)total_nodes, sizeof(int));

  
  migcount = (long **)xcalloc((size_t)total_nodes, sizeof(long *));
  for (i = 0; i < total_nodes; ++i)
    migcount[i] = (long *)xcalloc((size_t)total_nodes, sizeof(long));


  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gnode_t * x = gtree->nodes[i];

    if (i < gtree->tip_count)
      seqin_count[x->pop->node_index]++;

    curpop = x->pop;
    if (x->mi)
    {
      for (j = 0; j < x->mi->count; ++j)
      {
        long s = x->mi->me[j].source->node_index;
        long t = x->mi->me[j].target->node_index;

        assert(x->mi->me[j].target->parent);
        assert(x->mi->me[j].target->tau < x->mi->me[j].time);
        assert(x->mi->me[j].target->parent->tau > x->mi->me[j].time);

        migevent_count[s]++;
        migevent_count[t]++;
        migcount[t][s]++;

        if (curpop != x->mi->me[j].source)
        {
          for (pop = curpop->parent; pop != x->mi->me[j].source->parent; pop = pop->parent)
            seqin_count[pop->node_index]++;
        }
        curpop = x->mi->me[j].target;
      }
    }
    snode_t * end = x->parent ? x->parent->pop->parent : NULL;
    for (pop = curpop->parent; pop != end; pop = pop->parent)
      seqin_count[pop->node_index]++;
  }

  /* check */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * snode = stree->nodes[i];
    unsigned int node_index = snode->node_index;

   
    if (snode->seqin_count[msa_index] != seqin_count[node_index])
    {
      printf("\n");
      printf("[DEBUG]: Locus %ld snode %s\n", msa_index, snode->label);
      printf("seqin_count: %d    wrong: %d\n",
             seqin_count[node_index], snode->seqin_count[msa_index]);
      debug_print_gtree(gtree);
      fatal("Exiting");
    }
    //assert(snode->seqin_count[msa_index] == seqin_count[node_index]);
    if (snode->migevent_count[msa_index] != migevent_count[node_index])
    {
      printf("\n");
      printf("[DEBUG]: Locus %ld snode %s\n", msa_index, snode->label);
      printf("migevent_count: %ld    wrong: %ld\n",
             migevent_count[node_index], snode->migevent_count[msa_index]);
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        printf("  node %ld migevent_count %ld  wrong %ld  \n",
               j,migevent_count[j],stree->nodes[j]->migevent_count[msa_index]);
      debug_print_gtree(gtree);
      fatal("Exiting");
    }
    //assert(snode->migevent_count[msa_index] == migevent_count[node_index]);
  }
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      assert(gtree->migcount[i][j] == migcount[i][j]);

  /* deallocate */
  free(seqin_count);
  free(migevent_count);
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    free(migcount[i]);
  free(migcount);
}

void debug_consistency(stree_t * stree, gtree_t ** gtree_list, const char * msg)
{
  long i,j;
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gtree = gtree_list[i];
  
    for (j = 0; j < gtree->tip_count; ++j)
    {
      gnode_t * x = gtree->nodes[j];
      if (x->left || x->right)
        fatal("[ERROR] Locus %ld [%d] gnode %ld should be tip but has children",
              i, gtree->tip_count, x->node_index);
    }
    for (j = gtree->tip_count; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      gnode_t * x = gtree->nodes[j];

      if (x->time <= x->pop->tau)
        fatal("[FATAL-%ld] %s: Locus %ld [%d] gnode %ld has age %f but its population (%s) starts at %f",
              opt_debug_counter, msg, i, gtree->tip_count, x->node_index, x->time, x->pop->label, x->pop->tau);

      if (x->pop->parent && x->pop->parent->tau <= x->time)
        fatal("[FATAL-%ld] %s: Locus %ld [%d] gnode %ld [%s] has age %f but its population's parent (%s) starts at %f",
              opt_debug_counter, msg, i, gtree->tip_count, x->node_index, x->pop->label, x->time, x->pop->parent->label, x->pop->parent->tau);
      
      if (x->parent && x->parent->time <= x->time)
        fatal("[FATAL-%ld] %s: Locus %ld [%d] gnode %ld has age %f but its parent has older age %f",
              opt_debug_counter, msg, i, gtree->tip_count, x->node_index, x->time, x->parent->time);

      if (x->left->time >= x->time)
        fatal("[FATAL-%ld] %s: Locus %ld [%d] gnode %ld has age %f but its left daughter has older age %f",
              opt_debug_counter, msg, i, gtree->tip_count, x->node_index, x->time, x->left->time);
      if (x->right->time >= x->time)
        fatal("[FATAL-%ld] %s:  Locus %ld [%d] gnode %ld has age %f but its right daughter has older age %f",
              opt_debug_counter, msg, i, gtree->tip_count, x->node_index, x->time, x->right->time);
    }

    /* chekc if each node is in the correct population */
    for (j = 0; j < gtree->tip_count+gtree->inner_count; ++j)
    {
      gnode_t * x = gtree->nodes[j];
      snode_t * pop = gtree->nodes[j]->pop;

      if (x->time < pop->tau)
      {
        printf("node_index: %d time: %f pop: %s pop->tau: %f\n", x->node_index, x->time, x->pop->label, x->pop->tau);
        debug_print_gtree(gtree);
      }
      assert(x->time >= pop->tau);
      if (pop->parent && pop->parent->tau)
        assert(x->time < pop->parent->tau);
    }
    unsigned int event_count = 0;
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    {
      event_count += stree->nodes[j]->coal_count[i];
      dlist_item_t * item = stree->nodes[j]->coalevent[i]->head;
      int pop_event_count = 0;
      while (item)
      {
        gnode_t * x = (gnode_t *)(item->data);
        assert(x->pop == stree->nodes[j]);
        ++pop_event_count;

        item = item->next;
      }
      assert(pop_event_count == stree->nodes[j]->coal_count[i]);
    }
    assert(event_count == gtree->inner_count);
  }
}

void debug_write_migs_header(FILE * fp_debug, stree_t * stree)
{
  long j,k;

  fprintf(fp_debug, "Order:");
  for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
    for (k = 0; k < stree->tip_count+stree->inner_count; ++k)
      if (opt_mig_bitmatrix[j][k])
        fprintf(fp_debug, " MC_%s->%s", stree->nodes[j]->label, stree->nodes[k]->label);
  fprintf(fp_debug,"\n");
}
void debug_write_migs(FILE * fp_debug, stree_t * stree, gtree_t ** gtree, long sample)
{
  int print_perlocus = 0;
  int print_sample_number = 0;
  long i,j,k;
  
  long * perlocus  = (long *)xcalloc((size_t)opt_locus_count, sizeof(long));
  long * perthread = (long *)xcalloc((size_t)opt_locus_count, sizeof(long));

  if (print_sample_number)
    fprintf(fp_debug, "%ld\n", sample);
  for (i = 0; i < opt_locus_count; ++i)
  {
    for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      for (k = 0; k < stree->tip_count+stree->inner_count; ++k)
        if (opt_mig_bitmatrix[j][k])
        {
          if (print_perlocus)
            fprintf(fp_debug," %5ld", gtree[i]->migcount[j][k]);
          perlocus[i] += gtree[i]->migcount[j][k];
        }
    if (print_perlocus)
      fprintf(fp_debug, "  |  %ld\n", perlocus[i]);
  }

  #if 0
  fprintf(fp_debug, "Per-thread migrations:\n");
  #endif
  thread_info_t * ti = threads_ti();
  double mean = 0;
  double stdev = 0;
  for (i = 0; i < opt_threads; ++i)
  {
    thread_info_t * tip = ti+i;
    for (j = tip->locus_first; j < tip->locus_first+tip->locus_count; ++j)
      perthread[i] += perlocus[j];
    if (!i)
      fprintf(fp_debug, "%5ld", perthread[i]);
    else
      fprintf(fp_debug, " %5ld", perthread[i]);
    mean += perthread[i];
  }
  mean /= opt_threads;
  for (i = 0; i < opt_threads; ++i)
    stdev += (perthread[i] - mean)*(perthread[i] - mean);
  stdev = sqrt(stdev / opt_threads);
  
  #if 0
  fprintf(fp_debug, "  | %.2f %.2f\n", mean, stdev);
  #endif

  free(perlocus);
  free(perthread);
  #if 0
  fprintf(fp_debug,"\n\n");
  #else
  fprintf(fp_debug,"\n");
  #endif
}

void debug_print_wsji(stree_t * stree, int prec, const char * prefix_msg, FILE * fp)
{
  long nodes_count = stree->tip_count+stree->inner_count;
  long i,j,k;
  assert(opt_migration);

  long cnt = 0;

  printf("%s Wsji matrix:\n", prefix_msg);
  printf("%s    ", prefix_msg);
  for (i = 0; i < nodes_count; ++i)
    printf("%4ld", i);
  printf("\n");
  for (i = 0; i < nodes_count; ++i)
  {
    printf("%s %2ld [", prefix_msg, i);
    /* row exists only if migration rate W_{i->j} exists */
    if (!opt_est_geneflow)
    {
      for (j = 0; j < nodes_count; ++j)
        if (opt_mig_bitmatrix[i][j])
          break;
      if (j == nodes_count)
      {
        for (j = 0; j < nodes_count; ++j)
          printf(" N/A");
        printf(" ]\n");
        continue;
      }
    }
    
    if (i == stree->root->node_index)
    {
      for (j = 0; j < nodes_count; ++j)
        printf(" N/A");
      printf(" ]\n");
      continue;
    }
    for (j = 0; j < nodes_count; ++j)
    {
      if (!opt_est_geneflow && !opt_mig_bitmatrix[i][j])
      {
        printf(" N/A");
        continue;
      }
      
      if (j == i || j == stree->root->node_index)
      {
        printf(" N/A");
        continue;
      }

      assert(stree->Wsji[i][j]);
      printf(ANSI_COLOR_RED " %3ld" ANSI_COLOR_RESET, ++cnt);
      //for (k = 0; k < opt_locus_count; ++k)
      //  printf(" %.*f", prec, stree->Wsji[i][j][k]);
    }
    printf(" ]\n");
  }

  cnt = 0;
  for (i = 0; i < nodes_count; ++i)
  {
    if (!stree->Wsji[i]) continue;
    for (j = 0; j < nodes_count; ++j)
      if (stree->Wsji[i][j])
      {
        printf("%s " ANSI_COLOR_RED "%3ld" ANSI_COLOR_RESET ":", prefix_msg, ++cnt);
        for (k = 0; k < opt_locus_count; ++k)
          printf(" %.*f", prec, stree->Wsji[i][j][k]);
        printf("\n");
      }
  }
}
