/*
    Copyright (C) 2016-2019 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static void print_pretty_phylip(FILE * fp,
                                msa_t * msa,
                                int pad,
                                int every,
                                unsigned int * weights)
{
  long i,j;
  if (!msa) return;

  fprintf(fp, "%d %d P\n", msa->count, msa->length);
  for (i = 0; i < msa->count; ++i)
  {
    fprintf(fp, "%-*s", pad, msa->label[i]);
    for (j = 0; j < (long)(msa->length); ++j)
    {
      /* note: prints an extra space before the sequence */
      if (j % every == 0)
        fprintf(fp, " ");
      fprintf(fp, "%c", msa->sequence[i][j]);
    }
    fprintf(fp, "\n");
  }
  assert(weights);
  fprintf(fp, "%d", weights[0]);
  for (i = 1; i < msa->length; ++i)
    fprintf(fp, " %d", weights[i]);
  fprintf(fp, "\n");
}

void msa_print_phylip(FILE * fp,
                      msa_t ** msa,
                      long count,
                      unsigned int ** weights)
{
  int every = 10;  /* separate data in sequence every 10 bases */
  int pad = 4;
  int maxlen = 0;
  long i,j,k;

  /* find length of longest sequence label */
  for (i = 0; i < count; ++i)
  {
    for (j = 0; j < msa[i]->count; ++j)
    {
      k = strlen(msa[i]->label[j]);
      if (k > maxlen)
        maxlen = k;
    }
  }

  for (i = 0; i < count; ++i)
  {
    print_pretty_phylip(fp, msa[i], maxlen+pad, every, weights[i]);
    fprintf(fp,"\n");
  }
}

void msa_count_ambiguous_sites(msa_t * msa, const unsigned int * map)
{
  int i,j;
  int amb;

  msa->amb_sites_count = 0;

  if (msa->dtype == BPP_DATA_AA) return;
  assert(msa->dtype == BPP_DATA_DNA);

  for (i = 0; i < msa->length; ++i)
  {
    amb = 0;
    for (j = 0; j < msa->count; ++j)
      amb |= map[(int)(msa->sequence[j][i])];

    if (amb)
      msa->amb_sites_count++;
  }
}
static int * mark_ambiguous_sites(msa_t * msa, const unsigned int * map)
{
  int i,j;
  int amb;
  int * ambvector = (int *)xcalloc(msa->length, sizeof(int));

  msa->amb_sites_count = 0;

  for (i = 0; i < msa->length; ++i)
  {
    amb = 0;
    for (j = 0; j < msa->count; ++j)
      amb |= map[(int)(msa->sequence[j][i])];

    if (amb)
    {
      ambvector[i] = 1;
      msa->amb_sites_count++;
    }
  }

  return ambvector;
}

static int remove_ambiguous(msa_t * msa, int * ambiguous)
{
  int i,j,k;
  int amb_count = 0;

  for (i = 0; i < msa->length; ++i)
    amb_count += ambiguous[i];

  /* if all sites contain ambigous characters exit with error */
  if (amb_count == msa->length)
    return 0;

  /* we will move all ambiguous sites to the right end of the alignment */

  i = 0;  j = msa->length-1;
  for (i = 0, j = msa->length-1; i < msa->length && j >= 0; ++i, --j)
  {
    /* find next ambiguous site from left to right */
    while (i < msa->length && !ambiguous[i])
      ++i;

    /* find next non-ambigous site from right to left */
    while (j >= 0 && ambiguous[j])
      --j;

    /* if we moved all sites containing ambiguities to the right, break */
    if (j < i) break;

    /* swap sites */
    for (k = 0; k < msa->count; ++k)
      SWAP(msa->sequence[k][i],msa->sequence[k][j]);

    /* swap mark */
    SWAP(ambiguous[j],ambiguous[i]);
  }

  /* all ambiguous sites should now be at the right end. We will place
     terminating zeroes to the new sequence lengths and update the length
     variable */
  msa->length -= amb_count;

  for (k = 0; k < msa->count; ++k)
    msa->sequence[k][msa->length] = 0;

  return 1;
}

int msa_remove_ambiguous(msa_t * msa)
{
  int * ambiguous;
  int rc;

  /* get a vector indicating which sites are ambiguous */
  ambiguous = mark_ambiguous_sites(msa,pll_map_amb);

  /* remove ambiugous sites from alignment */
  rc = remove_ambiguous(msa,ambiguous);

  free(ambiguous);

  return rc;
}

void msa_destroy(msa_t * msa)
{
  int i;

  if (msa->label)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->label[i])
        free(msa->label[i]);
    free(msa->label);
  }

  if (msa->sequence)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->sequence[i])
        free(msa->sequence[i]);
    free(msa->sequence);
  }

  free(msa);
}

#define COLUMNS 5

static int longint_len(long x)
{
  return x ? (int)floor(log10(abs(x)))+1 : 1;
}

void msa_summary(msa_t ** msa_list, int msa_count)
{
  long i,j,k;

  char * labels[COLUMNS] = {"Locus",
                            "Sequences",
                            "Length",
                            "Ambiguous sites",
                            "JC69 Compressed"};

  int col_len[COLUMNS];

  /* initialize length of each column according to labels */
  for (i = 0; i < COLUMNS; ++i)
    col_len[i] = strlen(labels[i]) + 2;

  /* 'locus count' column length */
  col_len[0] = MAX(col_len[0], longint_len(msa_count)+2);

  /* 'sequences' column length */
  k = 0;
  for (i = 0; i < msa_count; ++i)
    k = MAX(k,msa_list[i]->count);
  col_len[1] = MAX(col_len[1], longint_len(k)+2);

  /* 'length' column length */
  k = 0;
  for (i = 0; i < msa_count; ++i)
    k = MAX(k,msa_list[i]->length);
  col_len[2] = MAX(col_len[2], longint_len(k)+2);

  /* 'Ambiguous sites' column length */
  k = 0;
  for (i = 0; i < msa_count; ++i)
    k = MAX(k,msa_list[i]->amb_sites_count);
  col_len[3] = MAX(col_len[3], longint_len(k)+2);

  /* 'JC69 compressed length' column length */
  k = 0;
  for (i = 0; i < msa_count; ++i)
    k = MAX(k,msa_list[i]->length);
  col_len[4] = MAX(col_len[4], longint_len(k)+2);

  printf("\n");

  /* print header column with labels centered */
  for (i = 0; i < COLUMNS; ++i)
  {
    long blanks = (col_len[i] - strlen(labels[i]))/2;
    for (j = 0; j < blanks; ++j) printf(" ");

    printf("%s", labels[i]);

    blanks = col_len[i] - blanks - strlen(labels[i]);
    for (j = 0; j < blanks; ++j) printf(" ");

    if (i != COLUMNS-1)
      printf("|");
  }
  printf("\n");

  /* print separator */
  for (i = 0; i < COLUMNS; ++i)
  {
    for (j = 0; j < col_len[i]; ++j)
      printf("-");
    if (i != COLUMNS-1)
      printf("+");
  }
  printf("\n");

  /* print table rows */
  for (i = 0; i < msa_count; ++i)
  {
    printf("%*ld |", col_len[0]-1, i+1);
    printf("%*d |", col_len[1]-1, msa_list[i]->count);
    printf("%*d |", col_len[2]-1, msa_list[i]->original_length);
    printf("%*d |", col_len[3]-1, msa_list[i]->amb_sites_count);
    printf("%*d ", col_len[4]-1, msa_list[i]->length);
    printf("\n");
  }
  printf("\n");
}
