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

static double *** migflow_mean;
static long *** migflow_seq_count;
static double *** migflow_seq_mean;

static FILE * migflow_openfile(const char * suffix, const char * mode)
{
  char * filename;
  FILE * fp;

  xasprintf(&filename, "%s.%s", opt_jobname, suffix);
  fp = xopen(filename, mode);
  free(filename);

  return fp;
}

static void migflow_destroy(stree_t * stree)
{
  unsigned int i,j;

  if (stree->migflow_tips)
    free(stree->migflow_tips);

  if (!stree->migflow_count) return;
  for (i = 0; i < stree->tip_count; ++i)
  {
    for (j = 0; j < (unsigned int)opt_migration_count; ++j)
      free(stree->migflow_count[i][j]);
    free(stree->migflow_count[i]);
  }
  free(stree->migflow_count);
}

void migflow_init(stree_t * stree, gtree_t ** gtree_list)
{
  unsigned int i,j;
  unsigned int tips = stree->tip_count;
  long bands = opt_migration_count;

  if (opt_est_geneflow)
    fatal("Migration flow tracking (ancestry = 1) is not supported "
          "with gene flow estimation (geneflow = 1)");

  stree->migflow_count = (long ***)xmalloc((size_t)tips*sizeof(long **));
  migflow_mean = (double ***)xmalloc((size_t)tips*sizeof(double **));
  for (i = 0; i < tips; ++i)
  {
    stree->migflow_count[i] = (long **)xmalloc((size_t)bands*sizeof(long *));
    migflow_mean[i] = (double **)xmalloc((size_t)bands*sizeof(double *));

    for (j = 0; j < (unsigned int)bands; ++j)
    {
      stree->migflow_count[i][j] =
        (long *)xcalloc((size_t)opt_locus_count,sizeof(long));
      migflow_mean[i][j] =
        (double *)xcalloc((size_t)opt_locus_count,sizeof(double));
    }
  }

  /* all species tips are relevant in MSC-M */
  stree->migflow_tips = (snode_t **)xcalloc((size_t)(tips+1),sizeof(snode_t *));
  for (i = 0; i < tips; ++i)
    stree->migflow_tips[i] = stree->nodes[i];

  /* allocate per-sequence per-locus tables */
  migflow_seq_count = (long ***)xmalloc((size_t)opt_locus_count*sizeof(long **));
  migflow_seq_mean = (double ***)xmalloc((size_t)opt_locus_count*sizeof(double **));
  for (i = 0; i < (unsigned int)opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    migflow_seq_count[i] = (long **)xcalloc((size_t)gt->tip_count,sizeof(long *));
    migflow_seq_mean[i] = (double **)xcalloc((size_t)gt->tip_count,sizeof(double *));
    for (j = 0; j < gt->tip_count; ++j)
    {
      gnode_t * gtip = gt->nodes[j];
      if (!gtip->pop) continue;

      migflow_seq_count[i][j] =
        (long *)xcalloc((size_t)bands,sizeof(long));
      migflow_seq_mean[i][j] =
        (double *)xcalloc((size_t)bands,sizeof(double));
    }
  }
}

void migflow_reset(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m = 0;
  snode_t * stip;

  while ((stip = stree->migflow_tips[m++]))
  {
    i = stip->node_index;
    for (j = 0; j < opt_migration_count; ++j)
    {
      for (k = 0; k < opt_locus_count; ++k)
        stree->migflow_count[i][j][k] = 0;
    }
  }

  /* zero per-sequence tables */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (migflow_seq_count[i][j] == NULL) continue;

      for (k = 0; k < opt_migration_count; ++k)
        migflow_seq_count[i][j][k] = 0;
    }
  }
}

void migflow_reset_mean(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m = 0;
  snode_t * stip;

  while ((stip = stree->migflow_tips[m++]))
  {
    i = stip->node_index;
    for (j = 0; j < opt_migration_count; ++j)
    {
      for (k = 0; k < opt_locus_count; ++k)
        migflow_mean[i][j][k] = 0;
    }
  }

  /* zero per-sequence mean tables */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (migflow_seq_mean[i][j] == NULL) continue;

      for (k = 0; k < opt_migration_count; ++k)
        migflow_seq_mean[i][j][k] = 0;
    }
  }
}

void migflow_print(stree_t * stree)
{
  long i,j,k,m = 0;
  snode_t * stip;

  /* go through tip nodes */
  while ((stip = stree->migflow_tips[m++]))
  {
    i = stip->node_index;
    /* go through all migration bands */
    for (j = 0; j < opt_migration_count; ++j)
    {
      migspec_t * spec = opt_mig_specs+j;
      printf("  mig_%s->%s: ",
             stree->nodes[spec->si]->label,
             stree->nodes[spec->ti]->label);

      /* go through all loci */
      double mean = 0;
      long seqin_count = 0;
      for (k = 0; k < opt_locus_count; ++k)
      {
        mean += migflow_mean[i][j][k];
        seqin_count += stip->seqin_count[k];
      }
      printf(": %f", mean / seqin_count);
    }
  }
}

void migflow_write(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m,r;
  snode_t * stip;
  FILE * fp_tip = NULL;
  FILE * fp_locus = NULL;
  FILE * fp_seq = NULL;

  if (!opt_migration || !migflow_mean || !stree || !stree->migflow_tips) return;

  /* build reverse mapping: perm[original_index] = sorted_index */
  long * perm = (long *)xmalloc((size_t)opt_locus_count * sizeof(long));
  for (k = 0; k < opt_locus_count; ++k)
    perm[gtree_list[k]->original_index] = k;

  fp_tip = migflow_openfile("migflow.per-tip.txt", "w");
  fp_locus = migflow_openfile("migflow.per-tip-per-locus.txt", "w");

  /* per-tip header */
  fprintf(fp_tip, "tip");
  for (j = 0; j < opt_migration_count; ++j)
  {
    migspec_t * spec = opt_mig_specs+j;
    fprintf(fp_tip, ",%s->%s",
            stree->nodes[spec->si]->label,
            stree->nodes[spec->ti]->label);
  }
  fprintf(fp_tip, "\n");

  /* per-tip rows */
  m = 0;
  while ((stip = stree->migflow_tips[m++]))
  {
    i = stip->node_index;

    double sum_lineages = 0.0;
    for (k = 0; k < opt_locus_count; ++k)
      sum_lineages += stree->nodes[i]->seqin_count[k];

    fprintf(fp_tip, "%s", stip->label ? stip->label : "");
    for (j = 0; j < opt_migration_count; ++j)
    {
      double sum_counts = 0.0;
      double frac = 0.0;

      for (k = 0; k < opt_locus_count; ++k)
        sum_counts += migflow_mean[i][j][k];

      if (sum_lineages > 0.0)
        frac = sum_counts / sum_lineages;

      fprintf(fp_tip, ",%f", frac);
    }
    fprintf(fp_tip, "\n");
  }

  /* per-locus header */
  fprintf(fp_locus, "Locus_number");
  m = 0;
  while ((stip = stree->migflow_tips[m++]))
  {
    fprintf(fp_locus, ",%s", stip->label ? stip->label : "");
    for (j = 0; j < opt_migration_count; ++j)
    {
      migspec_t * spec = opt_mig_specs+j;
      fprintf(fp_locus, ",%s->%s",
              stree->nodes[spec->si]->label,
              stree->nodes[spec->ti]->label);
    }
  }
  fprintf(fp_locus, "\n");

  /* per-locus rows (iterate in original locus order) */
  for (r = 0; r < opt_locus_count; ++r)
  {
    k = perm[r];
    fprintf(fp_locus, "%ld", r + 1);
    m = 0;
    while ((stip = stree->migflow_tips[m++]))
    {
      i = stip->node_index;
      fprintf(fp_locus, ",%s", stip->label ? stip->label : "");

      long denom = stree->nodes[i]->seqin_count[k];
      for (j = 0; j < opt_migration_count; ++j)
      {
        double frac = 0.0;
        if (denom > 0)
          frac = migflow_mean[i][j][k] / (double)denom;
        fprintf(fp_locus, ",%f", frac);
      }
    }
    fprintf(fp_locus, "\n");
  }

  /* per-sequence per-locus file */
  fp_seq = migflow_openfile("migflow.per-sequence-per-locus.txt", "w");

  /* header */
  fprintf(fp_seq, "Locus,Sequence");
  for (j = 0; j < opt_migration_count; ++j)
  {
    migspec_t * spec = opt_mig_specs+j;
    fprintf(fp_seq, ",%s->%s",
            stree->nodes[spec->si]->label,
            stree->nodes[spec->ti]->label);
  }
  fprintf(fp_seq, "\n");

  /* rows: one per (locus, sequence) pair (iterate in original locus order) */
  for (r = 0; r < opt_locus_count; ++r)
  {
    k = perm[r];
    gtree_t * gt = gtree_list[k];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (migflow_seq_mean[k][j] == NULL) continue;

      fprintf(fp_seq, "%ld,%s", r + 1,
              gt->nodes[j]->label ? gt->nodes[j]->label : "");
      for (i = 0; i < opt_migration_count; ++i)
        fprintf(fp_seq, ",%f", migflow_seq_mean[k][j][i]);
      fprintf(fp_seq, "\n");
    }
  }
  fclose(fp_seq);

  fclose(fp_tip);
  fclose(fp_locus);
  free(perm);
}

void migflow_update(stree_t * stree, gtree_t ** gtree_list)
{
  long i,j,k,m = 0;
  long band;
  long ft_src, ft_tgt;
  snode_t * stip;
  gtree_t * gtree;

  /* go through tip nodes */
  while ((stip = stree->migflow_tips[m++]))
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
          if (x->mi && x->mi->count > 0)
          {
            for (k = 0; k < x->mi->count; ++k)
            {
              /* map migevent to band index (backward -> forward time swap) */
              ft_src = x->mi->me[k].target->node_index;
              ft_tgt = x->mi->me[k].source->node_index;
              band = opt_migration_matrix[ft_src][ft_tgt];
              assert(band >= 0);

              /* binary indicator: count each sequence at most once per band */
              if (migflow_seq_count[i][j] &&
                  migflow_seq_count[i][j][band] == 0)
              {
                migflow_seq_count[i][j][band] = 1;
                stree->migflow_count[stip->node_index][band][i]++;
              }
            }
          }
          x = x->parent;
        }
      }
    }
  }
}

void migflow_update_mean(stree_t * stree, gtree_t ** gtree_list, long round)
{
  long i,j,k,m = 0;
  snode_t * stip;

  /* go through tip nodes */
  while ((stip = stree->migflow_tips[m++]))
  {
    i = stip->node_index;
    /* go through all migration bands */
    for (j = 0; j < opt_migration_count; ++j)
    {
      /* go through all loci */
      for (k = 0; k < opt_locus_count; ++k)
        RMEAN_UPDATE(migflow_mean[i][j][k],round,stree->migflow_count[i][j][k]);
    }
  }

  /* update per-sequence per-locus running means */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_t * gt = gtree_list[i];
    for (j = 0; j < gt->tip_count; ++j)
    {
      if (migflow_seq_mean[i][j] == NULL) continue;

      for (k = 0; k < opt_migration_count; ++k)
        RMEAN_UPDATE(migflow_seq_mean[i][j][k],round,migflow_seq_count[i][j][k]);
    }
  }
}
