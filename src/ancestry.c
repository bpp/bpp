/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#define RMEAN_UPDATE(m,r,c) (m) = ((m)*((r)-1)+(c)) / (double)(r);

static double *** ancestry_mean;
static long *** ancestry_seq_count;
static double *** ancestry_seq_mean;

static FILE * ancestry_openfile(const char * suffix, const char * mode)
{
  char * filename;
  FILE * fp;

  xasprintf(&filename, "%s.%s", opt_jobname, suffix);
  fp = xopen(filename, mode);
  free(filename);

  return fp;
}

static void ancestry_destroy(stree_t * stree)
{
  unsigned int i,j;

  if (stree->ancestry_tips)
    free(stree->ancestry_tips);

  if (!stree->ancestry_count) return;
  for (i = 0; i < stree->tip_count; ++i)
  {
    for (j = 0; j < opt_locus_count; ++j)
      free(stree->ancestry_count[i][j]);
    free(stree->ancestry_count[i]);
  }
  free(stree->ancestry_count);
}

void ancestry_init(stree_t * stree, gtree_t ** gtree_list)
{
  unsigned int i,j,k;
  unsigned int tips = stree->tip_count;
  unsigned int hybs = stree->hybrid_count;
  snode_t * tip;
  snode_t * hnode;
  snode_t * mnode;
  snode_t * pnode;
  long ** tip_counts = NULL;
  double ** tip_means = NULL;

  stree->ancestry_count = (long ***)xmalloc((size_t)tips*sizeof(long **));
  ancestry_mean = (double ***)xmalloc((size_t)tips*sizeof(double **));
  for (i = 0; i < stree->tip_count; ++i)
  {
    tip = stree->nodes[i];
    stree->ancestry_count[i] = (long **)xcalloc((size_t)hybs,sizeof(long *));
    ancestry_mean[i] = (double **)xcalloc((size_t)hybs,sizeof(double *));
    tip_counts = stree->ancestry_count[i];
    tip_means = ancestry_mean[i];

    for (j = 0; j < hybs; ++j)
    {
      mnode = stree->nodes[stree->tip_count+stree->inner_count+j];
      hnode = mnode->hybrid;
      pnode = hnode->has_phi ? hnode : mnode;

      tip_counts[j] = stree->pptable[tip->node_index][pnode->node_index] ?
                  (long *)xcalloc((size_t)opt_locus_count,sizeof(long)) : NULL;
      tip_means[j] = stree->pptable[tip->node_index][pnode->node_index] ?
                  (double *)xcalloc((size_t)opt_locus_count,sizeof(double)) : NULL;
    }
  }

  stree->ancestry_tips = (snode_t **)xcalloc((size_t)(tips+1),sizeof(snode_t *));

  /* collect tips that have root-paths passing hybridizations */
  for (i = 0, k = 0; i < tips; ++i)
  {
    tip = stree->nodes[i];

    for (j = 0; j < hybs; ++j)
    {
      hnode = stree->nodes[stree->tip_count+stree->inner_count+j]->hybrid;

      if (stree->pptable[tip->node_index][hnode->node_index])
      {
        stree->ancestry_tips[k++] = tip;
        break;
      }
    }
  }

  /* allocate per-sequence per-locus tables */
  ancestry_seq_count = (long ***)xmalloc((size_t)opt_locus_count*sizeof(long **));
  ancestry_seq_mean = (double ***)xmalloc((size_t)opt_locus_count*sizeof(double **));
  for (i = 0; i < (unsigned int)opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    ancestry_seq_count[i] = (long **)xcalloc((size_t)gt->tip_count,sizeof(long *));
    ancestry_seq_mean[i] = (double **)xcalloc((size_t)gt->tip_count,sizeof(double *));
    for (j = 0; j < gt->tip_count; ++j)
    {
      gnode_t * gtip = gt->nodes[j];
      if (!gtip->pop) continue;

      /* check if this gene tree tip belongs to a species tip in ancestry_tips */
      int found = 0;
      unsigned int m = 0;
      while (stree->ancestry_tips[m])
      {
        if (gtip->pop == stree->ancestry_tips[m])
        {
          found = 1;
          break;
        }
        m++;
      }
      if (found)
      {
        ancestry_seq_count[i][j] = (long *)xcalloc((size_t)hybs,sizeof(long));
        ancestry_seq_mean[i][j] = (double *)xcalloc((size_t)hybs,sizeof(double));
      }
    }
  }
}

void ancestry_reset(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m = 0;
  snode_t * stip;

  while ((stip = stree->ancestry_tips[m++]))
  {
    i = stip->node_index;
    for (j = 0; j < stree->hybrid_count; ++j)
    {
      if (stree->ancestry_count[i][j] == NULL) continue;

      for (k = 0; k < opt_locus_count; ++k)
        stree->ancestry_count[i][j][k] = 0;
    }
  }

  /* zero per-sequence tables */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (ancestry_seq_count[i][j] == NULL) continue;

      for (k = 0; k < stree->hybrid_count; ++k)
        ancestry_seq_count[i][j][k] = 0;
    }
  }
}

void ancestry_reset_mean(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m = 0;
  snode_t * stip;

  while ((stip = stree->ancestry_tips[m++]))
  {
    i = stip->node_index;
    for (j = 0; j < stree->hybrid_count; ++j)
    {
      if (ancestry_mean[i][j] == NULL) continue;

      for (k = 0; k < opt_locus_count; ++k)
        ancestry_mean[i][j][k] = 0;
    }
  }

  /* zero per-sequence mean tables */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (ancestry_seq_mean[i][j] == NULL) continue;

      for (k = 0; k < stree->hybrid_count; ++k)
        ancestry_seq_mean[i][j][k] = 0;
    }
  }
}
void ancestry_print(stree_t * stree)
{
  long i,j,k,m = 0;
  snode_t * stip;
  snode_t * mnode;
  snode_t * hnode;
  snode_t * pnode;

  /* go through tip nodes with root-paths passing hybridizations */
  while ((stip = stree->ancestry_tips[m++]))
  {
    i = stip->node_index;
    /* go through all hybridizations */
    for (j = 0; j < stree->hybrid_count; ++j)
    {
      if (ancestry_mean[i][j] == NULL) continue;

      mnode = stree->nodes[stree->tip_count+stree->inner_count+j];
      snode_t * hnode = mnode->hybrid;
      pnode = hnode->has_phi ? hnode : mnode;
      printf("  phi_%s->%s: ", pnode->parent->label, pnode->label);

      /* go through all loci */
      double mean = 0;
      long seqin_count = 0;
      for (k = 0; k < opt_locus_count; ++k)
      {
        //printf(" %f", ancestry_mean[i][j][k]);
        mean += ancestry_mean[i][j][k];
        seqin_count += stip->seqin_count[k];
      }
      printf(": %f", mean / seqin_count);
    }
  }
}

void ancestry_write(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m;
  snode_t * stip;
  snode_t * mnode;
  snode_t * hnode;
  snode_t * pnode;
  FILE * fp_tip = NULL;
  FILE * fp_locus = NULL;
  FILE * fp_seq = NULL;

  if (!opt_msci || !ancestry_mean || !stree || !stree->ancestry_tips) return;

  fp_tip = ancestry_openfile("ancestry.per-tip.txt", "w");
  fp_locus = ancestry_openfile("ancestry.per-tip-per-locus.txt", "w");

  /* per-tip header */
  fprintf(fp_tip, "tip");
  for (j = 0; j < stree->hybrid_count; ++j)
  {
    mnode = stree->nodes[stree->tip_count+stree->inner_count+j];
    hnode = mnode->hybrid;
    pnode = hnode->has_phi ? hnode : mnode;
    fprintf(fp_tip, ",%s->%s",
            pnode->parent ? pnode->parent->label : "",
            pnode->label ? pnode->label : "");
  }
  fprintf(fp_tip, "\n");

  /* per-tip rows */
  m = 0;
  while ((stip = stree->ancestry_tips[m++]))
  {
    i = stip->node_index;

    double sum_lineages = 0.0;
    for (k = 0; k < opt_locus_count; ++k)
      sum_lineages += stree->nodes[i]->seqin_count[k];

    fprintf(fp_tip, "%s", stip->label ? stip->label : "");
    for (j = 0; j < stree->hybrid_count; ++j)
    {
      double sum_counts = 0.0;
      double frac = 0.0;

      if (ancestry_mean[i][j])
      {
        for (k = 0; k < opt_locus_count; ++k)
          sum_counts += ancestry_mean[i][j][k];
      }
      if (sum_lineages > 0.0)
        frac = sum_counts / sum_lineages;

      fprintf(fp_tip, ",%f", frac);
    }
    fprintf(fp_tip, "\n");
  }

  /* per-locus header */
  fprintf(fp_locus, "Locus_number");
  m = 0;
  while ((stip = stree->ancestry_tips[m++]))
  {
    fprintf(fp_locus, ",%s", stip->label ? stip->label : "");
    for (j = 0; j < stree->hybrid_count; ++j)
    {
      mnode = stree->nodes[stree->tip_count+stree->inner_count+j];
      hnode = mnode->hybrid;
      pnode = hnode->has_phi ? hnode : mnode;
      fprintf(fp_locus, ",%s->%s",
              pnode->parent ? pnode->parent->label : "",
              pnode->label ? pnode->label : "");
    }
  }
  fprintf(fp_locus, "\n");

  /* per-locus rows */
  for (k = 0; k < opt_locus_count; ++k)
  {
    fprintf(fp_locus, "%ld", k + 1);
    m = 0;
    while ((stip = stree->ancestry_tips[m++]))
    {
      i = stip->node_index;
      fprintf(fp_locus, ",%s", stip->label ? stip->label : "");

      long denom = stree->nodes[i]->seqin_count[k];
      for (j = 0; j < stree->hybrid_count; ++j)
      {
        double frac = 0.0;
        if (denom > 0 && ancestry_mean[i][j])
          frac = ancestry_mean[i][j][k] / (double)denom;
        fprintf(fp_locus, ",%f", frac);
      }
    }
    fprintf(fp_locus, "\n");
  }

  /* per-sequence per-locus file */
  fp_seq = ancestry_openfile("ancestry.per-sequence-per-locus.txt", "w");

  /* header */
  fprintf(fp_seq, "Locus,Sequence");
  for (j = 0; j < stree->hybrid_count; ++j)
  {
    mnode = stree->nodes[stree->tip_count+stree->inner_count+j];
    hnode = mnode->hybrid;
    pnode = hnode->has_phi ? hnode : mnode;
    fprintf(fp_seq, ",%s->%s",
            pnode->parent ? pnode->parent->label : "",
            pnode->label ? pnode->label : "");
  }
  fprintf(fp_seq, "\n");

  /* rows: one per (locus, sequence) pair */
  for (k = 0; k < opt_locus_count; ++k)
  {
    gtree_t * gt = gtree_list[k];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (ancestry_seq_mean[k][j] == NULL) continue;

      fprintf(fp_seq, "%ld,%s", k + 1,
              gt->nodes[j]->label ? gt->nodes[j]->label : "");
      for (i = 0; i < stree->hybrid_count; ++i)
        fprintf(fp_seq, ",%f", ancestry_seq_mean[k][j][i]);
      fprintf(fp_seq, "\n");
    }
  }
  fclose(fp_seq);

  fclose(fp_tip);
  fclose(fp_locus);
}

void ancestry_update(stree_t * stree, gtree_t ** gtree_list)
{
  int altpath,phipath;
  long i,j,k,m = 0;
  snode_t * mnode;
  snode_t * hnode;
  snode_t * stip;
  gtree_t * gtree;

  
  /* go through tip nodes with root-paths passing hybridizations */
  while ((stip = stree->ancestry_tips[m++]))
  {
    /* go through all loci */
    for (i = 0; i < opt_locus_count; ++i)
    {
      gtree = gtree_list[i];
      /* go through sequences of current stip */
      for (j = 0; j < gtree->tip_count; ++j)
      {
        gnode_t * x = gtree->nodes[j];
        if (x->pop != stip) continue;

        while (x)
        {
          for (k = 0; k < stree->hybrid_count; ++k)
          {
            if (x->hpath[k] == BPP_HPATH_NONE) continue;

            mnode = stree->nodes[stree->tip_count+stree->inner_count+k];
            hnode = mnode->hybrid;
            /* pnode = hnode->has_phi ? hnode : mnode; */
            phipath = hnode->has_phi ? BPP_HPATH_LEFT : BPP_HPATH_RIGHT;
            altpath = mnode->has_phi ? BPP_HPATH_LEFT : BPP_HPATH_RIGHT;

            assert(phipath != altpath);
            assert(hnode->has_phi ^ mnode->has_phi);

            if (x->hpath[k] == phipath)
            {
              stree->ancestry_count[stip->node_index][k][i]++;
              if (ancestry_seq_count[i][j])
                ancestry_seq_count[i][j][k] = 1;
            }
          }
          x = x->parent;
        }
      }
    }
  }
}

void ancestry_update_mean(stree_t * stree, gtree_t ** gtree_list, long round)
{
  long i,j,k,m = 0;
  snode_t * stip;

  /* go through tip nodes with root-paths passing hybridizations */
  while ((stip = stree->ancestry_tips[m++]))
  {
    i = stip->node_index;
    /* go through all hybridizations */
    for (j = 0; j < stree->hybrid_count; ++j)
    {
      if (ancestry_mean[i][j] == NULL) continue;

      /* go through all loci */
      for (k = 0; k < opt_locus_count; ++k)
        RMEAN_UPDATE(ancestry_mean[i][j][k],round,stree->ancestry_count[i][j][k]);
    }
  }

  /* update per-sequence per-locus running means */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (ancestry_seq_mean[i][j] == NULL) continue;

      for (k = 0; k < stree->hybrid_count; ++k)
        RMEAN_UPDATE(ancestry_seq_mean[i][j][k],round,ancestry_seq_count[i][j][k]);
    }
  }
}

