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

static void vecswap(int i, int j, int n, char ** x, int * oi)
{
  while (n--)
  {
    SWAP(x[i],x[j]);
    if (oi)
      SWAP(oi[i],oi[j]);
    ++i; ++j;
  }
}

static void ssort1(char ** x, int n, int depth, int * oi)
{
  int a,b,c,d,r,v;

  if (n <= 1) return;

  a = rand() % n;

  SWAP(x[0], x[a]);
  if (oi)
    SWAP(oi[0], oi[a]);

  v = x[0][depth];

  a = b = 1;
  c = d = n-1;

  while (1)
  {
    while (b <= c && (r = x[b][depth]-v) <= 0)
    {
      if (r == 0)
      {
        SWAP(x[a], x[b]);
        if (oi)
          SWAP(oi[a], oi[b]);
        ++a;
      }
      ++b;
    }
    while (b <= c && (r = x[c][depth]-v) >= 0)
    {
      if (r == 0)
      {
        SWAP(x[c], x[d]);
        if (oi)
          SWAP(oi[c], oi[d]);
        --d;
      }
      --c;
    }
    if (b > c) break;
    SWAP(x[b], x[c]);
    if (oi)
      SWAP(oi[b], oi[c]);
    ++b; --c;
  }

  r = MIN(a,b-a); vecswap(0,b-r,r,x,oi);
  r = MIN(d-c,n-d-1); vecswap(b,n-r,r,x,oi);
  r = b-a; ssort1(x,r,depth,oi);

  if (x[r][depth] != 0)
  {
    if (oi)
      ssort1 (x + r, a + n - d - 1, depth + 1, oi + r);
    else
      ssort1 (x + r, a + n - d - 1, depth + 1, NULL);

  }

    r = d - c; 
    if (oi)
      ssort1(x+n-r,r,depth,oi+n-r);
    else
      ssort1(x+n-r,r,depth,NULL);
}

static void remap_range(const unsigned int * map,
                        unsigned char * charmap)
{
  unsigned int oldmap[ASCII_SIZE];
  unsigned int i,j;
  unsigned char k = 1;

  memcpy(oldmap, map, ASCII_SIZE * sizeof(unsigned int));
  memset(charmap, 0, ASCII_SIZE * sizeof(unsigned char));

  for (i = 0; i < ASCII_SIZE; ++i)
    if (oldmap[i])
    {
      charmap[i] = k;

      for (j = i+1; j < ASCII_SIZE; ++j)
        if (oldmap[i] == oldmap[j])
        {
          charmap[j] = k;
          oldmap[j] = 0;
        }

      ++k;
    }
}

static unsigned int findmax(const unsigned int * map)
{
  int i;
  unsigned int max = 0;

  for (i = 0; i < ASCII_SIZE; ++i)
    if (map[i] > max)
      max = map[i];

  return max;
}

static void encode(char ** sequence,
                   const unsigned char * map,
                   int count,
                   int len)
{
  int i,j;
  char * p;

  for (i = 0; i < count; ++i)
  {
    p = sequence[i];
    j = len;
    while (j--)
    {
      *p = map[(int)(*p)];
      ++p;
    }
  }
}

static unsigned char ** encode_jc69(char ** column,
                                    int count,
                                    int len)
{
  int i,j;
  char * p;
  unsigned char ** jc69_invmaps;
  unsigned char sitemap[16];
  int convert;

  /* allocate memory for inverse map in order to later decode sequence data */
  jc69_invmaps = (unsigned char **)xcalloc(count,sizeof(unsigned char *)); 

  /* go through the sites */
  for (i = 0; i < count; ++i)
  {
    convert = 1;
    p = column[i];
    j = len;

    /* we re-encode only sites that have no ambiguities with the exception of
       gaps */
    while (j--)
    {
      convert &= (*p > 15) ? 0 : pll_map_validjc69[(unsigned int)*p];
      ++p;
    }


    /* if no ambiguities apart gaps were found, encode the column and create
       an inverse map to decode back later */
    if (convert)
    {
      p = column[i];
      jc69_invmaps[i] = (unsigned char *)xcalloc(16,sizeof(unsigned char));
      memset(sitemap,0,16*sizeof(unsigned char));
      sitemap[15] = jc69_invmaps[i][15] = 15;  /* gaps do not get converted */
      unsigned char code = 1;
      for ( j = 0; j < len; ++j)
      {
        unsigned int c = (unsigned int)(p[j]);
        if (!sitemap[c])
        {
          sitemap[c] = code;
          jc69_invmaps[i][code] = p[j];
          code++;
        }

        p[j] = sitemap[c];
      }
    }
  }

  /* return the inverse maps */
  return jc69_invmaps;
}

unsigned int * compress_site_patterns(char ** sequence,
                                      const unsigned int * map,
                                      int count,
                                      int * length,
                                      int attrib)
{
  int i,j;
  char * memptr;
  char ** column;
  unsigned int * weight;
  unsigned char ** jc69_invmaps = NULL;

  unsigned char charmap[ASCII_SIZE];
  unsigned char inv_charmap[ASCII_SIZE];

  /* check that at least one sequence is given */
  if (!count) return NULL;

  /* a map must be given */
  if (!map) return NULL;

  /* a zero can never be used as a state */
  if (map[0]) return NULL;

  /* if map states are out of the BYTE range, remap */
  if (findmax(map) >= ASCII_SIZE)
  {
    remap_range(map,charmap);

    /* for now only DNA with pll_map_nt */
    assert(0);
  }
  else
  {
    for (i = 0; i < ASCII_SIZE; ++i)
      charmap[i] = (unsigned char)(map[i]);
  }

  /* create inverse charmap to decode states back to characters when
     compression is finished */
  for (i = 0; i < ASCII_SIZE; ++i)
    if (map[i])
      inv_charmap[charmap[i]] = (unsigned char)i;

  /* encode sequences using charmap */
  encode(sequence,charmap,count,*length);

  /* allocate memory for columns */
  column = (char **)xmalloc((size_t)(*length)*sizeof(char *));

  /* allocate memory for the alignment */
  memptr = column[0] = (char *)xmalloc((size_t)(*length) *
                                       (size_t)(count+1) *
                                       sizeof(char));

  /* map memory to each column */
  for (i = 1; i < *length; ++i)
    column[i] = column[i-1] + (count+1);

  /* allocate space for weight vector */
  weight = (unsigned int *)xmalloc((size_t)(*length)*sizeof(unsigned int));

  /* split alignment into columns instead of rows */
  for (i = 0; i < (*length); ++i)
  {
    for (j = 0; j < count; ++j)
      column[i][j] = sequence[j][i];
    column[i][j] = 0;
  }

  /* allocate space for storing original indices (before sorting sites) */
  int * oi = (int *)xmalloc(*length * sizeof(int));
  for (i = 0; i < *length; ++i)
    oi[i] = i;

    /* do the jc69 now */
  if (attrib == COMPRESS_JC69)
    jc69_invmaps = encode_jc69(column,*length,count);

  /* sort the columns and keep original indices */
  ssort1(column, *length, 0, oi);

  /* we have at least one unique site with weight 1 (the first site) */
  int compressed_length = 1;
  size_t ref = 0;
  weight[ref] = 1;

  /* find all unique columns and set their weights */
  int * compressed_oi = (int *)xmalloc(*length * sizeof(int));

  compressed_oi[0] = oi[0];
  for (i = 1; i < *length; ++i)
  {
    if (strcmp(column[i],column[i-1]))
    {
      column[ref+1] = column[i];
      compressed_oi[ref+1] = oi[i];
      ++ref;
      ++compressed_length;
      weight[ref] = 1;
    }
    else
      weight[ref]++;
  }

  /* decode the jc69 encoding */
  if (attrib == COMPRESS_JC69)
  {
    for (i=0; i < compressed_length; ++i)
    {
      unsigned char * sitemap = jc69_invmaps[compressed_oi[i]];
      if (sitemap)
        for (j=0; j < count; ++j)
          column[i][j] = sitemap[(unsigned int)column[i][j]];
    }
    for (i = 0; i < *length; ++i)
      if (jc69_invmaps[i])
        free(jc69_invmaps[i]);
    free(jc69_invmaps);
  }

  /* copy the unique columns over the original sequences */
  for (i = 0; i < compressed_length; ++i)
    for (j = 0; j < count; ++j)
      sequence[j][i] = column[i][j];

  /* add terminating zero */
  for (j = 0; j < count; ++j)
    sequence[j][compressed_length] = 0;

  /* deallocate memory */
  free(memptr);
  free(column);

  /* adjust weight vector size to compressed length */
  unsigned int * mem = (unsigned int *)xmalloc((size_t)compressed_length *
                                               sizeof(unsigned int));
  if (mem)
  {
    /* copy weights */
    for (i = 0; i < compressed_length; ++i)
      mem[i] = weight[i];

    /* free and re-point */
    free(weight);
    weight = mem;
  }

  /* update length */
  *length = compressed_length;

  /* decode sequences using inv_charmap */
  encode(sequence,inv_charmap,count,compressed_length);

  free(oi);
  free(compressed_oi);

  return weight;
}

unsigned long * compress_site_patterns_diploid(char ** sequence,
                                               const unsigned int * map,
                                               int count,
                                               int * length,
                                               unsigned int ** wptr,
                                               int attrib)
{
  int i,j;
  char * memptr;
  char ** column;
  unsigned int * weight;
  unsigned long * mapping;
  unsigned char ** jc69_invmaps = NULL;

  unsigned char charmap[ASCII_SIZE];
  unsigned char inv_charmap[ASCII_SIZE];

  /* check that at least one sequence is given */
  if (!count) return NULL;

  /* a map must be given */
  if (!map) return NULL;

  /* a zero can never be used as a state */
  if (map[0]) return NULL;

  /* if map states are out of the BYTE range, remap */
  if (findmax(map) >= ASCII_SIZE)
  {
    remap_range(map,charmap);

    /* for now only DNA with pll_map_nt */
    assert(0);
  }
  else
  {
    for (i = 0; i < ASCII_SIZE; ++i)
      charmap[i] = (unsigned char)(map[i]);
  }

  /* create inverse charmap to decode states back to characters when
     compression is finished */
  for (i = 0; i < ASCII_SIZE; ++i)
    if (map[i])
      inv_charmap[charmap[i]] = (unsigned char)i;

  /* encode sequences using charmap */
  encode(sequence,charmap,count,*length);

  /* allocate memory for columns */
  column = (char **)xmalloc((size_t)(*length)*sizeof(char *));

  /* allocate memory for the alignment */
  memptr = column[0] = (char *)xmalloc((size_t)(*length) *
                                       (size_t)(count+1) *
                                       sizeof(char));

  /* map memory to each column */
  for (i = 1; i < *length; ++i)
    column[i] = column[i-1] + (count+1);

  /* allocate space for weight vector */
  weight = (unsigned int *)xmalloc((size_t)(*length)*sizeof(unsigned int));

  /* allocate space for mapping vector */
  mapping = (unsigned long *)xmalloc((size_t)(*length)*sizeof(unsigned long));

  /* split alignment into columns instead of rows */
  for (i = 0; i < (*length); ++i)
  {
    for (j = 0; j < count; ++j)
      column[i][j] = sequence[j][i];
    column[i][j] = 0;
  }

  /* allocate space for storing original indices (before sorting sites) */
  int * oi = (int *)xmalloc(*length * sizeof(int));
  for (i = 0; i < *length; ++i)
    oi[i] = i;

    /* do the jc69 now */
  if (attrib == COMPRESS_JC69)
    jc69_invmaps = encode_jc69(column,*length,count);

  /* sort the columns and keep original indices */
  ssort1(column, *length, 0, oi);

  /*first site in uncompressed alignment maps to first site in compressed */
  int compressed_length = 1;
  size_t ref = 0;
  mapping[0] = 0;
  weight[ref] = 1;

  /* find all unique columns and set their mappings A2->A3 */
  int * compressed_oi = (int *)xmalloc(*length * sizeof(int));

  compressed_oi[0] = oi[0];
  for (i = 1; i < *length; ++i)
  {
    if (strcmp(column[i],column[i-1]))
    {
      column[ref+1] = column[i];
      compressed_oi[ref+1] = oi[i];
      ++ref;
      ++compressed_length;
      weight[ref] = 1;
    }
    else
      weight[ref]++;

    /* map original index i in uncompressed alignment to index in compressed */
    mapping[oi[i]] = ref;
  }

  /* decode the jc69 encoding */
  if (attrib == COMPRESS_JC69)
  {
    for (i=0; i < compressed_length; ++i)
    {
      unsigned char * sitemap = jc69_invmaps[compressed_oi[i]];
      if (sitemap)
        for (j=0; j < count; ++j)
          column[i][j] = sitemap[(unsigned int)column[i][j]];
    }
    for (i = 0; i < *length; ++i)
      if (jc69_invmaps[i])
        free(jc69_invmaps[i]);
    free(jc69_invmaps);
  }

  /* copy the unique columns over the original sequences */
  for (i = 0; i < compressed_length; ++i)
    for (j = 0; j < count; ++j)
      sequence[j][i] = column[i][j];

  /* add terminating zero */
  for (j = 0; j < count; ++j)
    sequence[j][compressed_length] = 0;

  /* deallocate memory */
  free(memptr);
  free(column);


  /* adjust weight vector size to compressed length */
  unsigned int * mem = (unsigned int *)xmalloc((size_t)compressed_length *
                                               sizeof(unsigned int));
  if (mem)
  {
    /* copy weights */
    for (i = 0; i < compressed_length; ++i)
      mem[i] = weight[i];

    /* free and re-point */
    free(weight);
    weight = mem;
  }
  *wptr = weight;

  /* update length */
  *length = compressed_length;

  /* decode sequences using inv_charmap */
  encode(sequence,inv_charmap,count,compressed_length);

  free(oi);
  free(compressed_oi);

  return mapping;
}
