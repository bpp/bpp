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

#define SWAP(x,y) do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)

static void vecswap(int i, int j, int n, char ** x)
{
  while (n-- > 0)
  {
    SWAP(x[i], x[j]);
    ++ i; ++ j;
  }
}

static void ssort1(char ** x, int n, int depth)
{
  int           a, b, c, d, r, v;

  if (n <= 1) return;

  a = rand() % n;

  SWAP(x[0], x[a]);

  v = x[0][depth];

  a = b = 1;
  c = d = n-1;

  for (;;)
  {
    while (b <= c && (r = x[b][depth]-v) <= 0)
    {
      if (r == 0)
      {
        SWAP(x[a], x[b]);
        ++a;
      }
      ++b;
    }
    while (b <= c && (r = x[c][depth]-v) >= 0)
    {
      if (r == 0)
      {
        SWAP(x[c], x[d]);
        --d;
      }
      --c;
    }
    if (b > c) break;
    SWAP (x[b], x[c]);
    ++b; --c;
  }
  r = MIN(a,b-a); vecswap(0,b-r,r,x);
  r = MIN(d-c,n-d-1); vecswap(b,n-r,r,x);
  r = b-a; ssort1(x,r,depth);

  if (x[r][depth] != 0)
    ssort1 (x + r, a + n - d - 1, depth + 1);

  r = d - c; ssort1(x+n-r,r,depth);
}

static void pllssort1main(char ** x, int n)
{
  ssort1(x, n, 0);
}


void compress(alignment_t * locus)
{
  int i,j;
  char ** sites;
  unsigned int * weight = locus->weight;
  int dups;

  /* allocate space for transposed alignment */
  sites = (char **)xmalloc(locus->seq_len * sizeof(char *));

  /* transpose alignment */
  for (i = 0; i < locus->seq_len; ++i)
  {
    for (j = 0; j < locus->seq_count; ++j)
      sites[i][j] = locus->seq[j][i];
    sites[i][j] = 0;
  }

  /* sort the sites */
  pllssort1main(locus->seq, locus->seq_len);


  for (i = 1; i < locus->seq_len; ++i) 
  {
    if (!strcmp(sites[i], sites[i-1]))
    {
      ++dups;
    }
  }
}

void alignment_print(alignment_t * alignment)
{
  if (!alignment) return;

  int i;

  printf("%d %d\n", alignment->seq_count, alignment->seq_len);
  for (i = 0; i < alignment->seq_count; ++i)
    printf("%s %s\n", alignment->label[i], alignment->seq[i]);
}

void alignment_destroy(alignment_t * alignment)
{
  int i;

  for (i = 0; i < alignment->seq_count; ++i)
  {
    free(alignment->seq[i]);
    free(alignment->label[i]);
  }
  free(alignment->seq);
  free(alignment->label);

  if (alignment->weight)
    free(alignment->weight);

  free(alignment);
}
