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

void msa_print(msa_t * msa)
{
  if (!msa) return;

  int i;

  printf("%d %d\n", msa->count, msa->length);
  for (i = 0; i < msa->count; ++i)
    printf("%s %s\n", msa->label[i], msa->sequence[i]);
}

static int * mark_ambiguous_sites(msa_t * msa, const unsigned int * map)
{
  int i,j;
  int amb;
  int * ambvector = (int *)xcalloc(msa->length, sizeof(int));

  for (i = 0; i < msa->length; ++i)
  {
    amb = 0;
    for (j = 0; j < msa->count; ++j)
      amb |= map[(int)(msa->sequence[j][i])];

    if (amb)
      ambvector[i] = 1;
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
    while (!ambiguous[i])
      ++i;

    /* find next non-ambigous site from right to left */
    while (ambiguous[j])
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
