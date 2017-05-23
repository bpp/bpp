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
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "bpp.h"

static int mark_ambiguous_sites(msa_t * msa,
                                const unsigned int * map)
{
  int i,j,k;
  unsigned int c;
  int clean_count = 0;

  int len = msa->length;
  int count = msa->count;

  char ** seq = msa->sequence;

  /* iterate sites */
  for (i = 0; i < len; ++i)
  {
    /* iterate sequences */
    for (j = 0; j < count; ++j)
    {
      /* convert character to state code */
      c = (map[(unsigned int)(seq[j][i])]);

      /* get number of bits set */
      k = __builtin_popcount(c);

      /* if no bits are set then invalid state */
      if (!k)
        fatal("Illegal state code \%c\" in tip %s",
              seq[j][i], msa->label[j]);

      /* if ambiguity, set the whole site to an ambiguous state and increase
         the number of sites to be removed */
      if (k > 1)
      {
        c = seq[j][i];
        for (k = 0; k < count; ++k)
          seq[k][i] = c;
        clean_count++;
        break;
      }
    }
  }

  return clean_count;
}

static void set_tipstates(locus_t * locus,
                          msa_t * msa,
                          const unsigned int * map)
{
  int i,j,k;
  unsigned int c;
  int clean_count = 0;
  
  int len = msa->length;
  int count = msa->count;

  /* if --cleandata then get number of sites to be removed */
  if (opt_cleandata)
    clean_count = mark_ambiguous_sites(msa,map);

  /* allocate necessary space for storing tip sequences */
  locus->tipstates = (unsigned char **)xmalloc(count*sizeof(unsigned char *));
  for (i = 0; i < count; ++i)
    locus->tipstates[i] = (unsigned char *)xmalloc((len-clean_count) *
                                                   sizeof(unsigned char));
    
  /* iterate sequences */
  for (i = 0; i < count; ++i)
  {
    char * seq = msa->sequence[i];

    /* go through the current sequence */
    for (j=0, k=0; j < len; ++j)
    {
      /* convert sequence character to state code */
      c = map[(unsigned int)(seq[j])];

      /* seq[j] is an illegal character */
      if (!c)
        fatal("Illegal state code \%c\" in tip %s",
              seq[j], msa->label[i]);

      /* if cleandata enabled and ambiguous state then skip */
      if (opt_cleandata && __builtin_popcount(c) > 1)

      /* store state */
      locus->tipstates[i][k++] = c;
    }

  }
}

void init_locus(locus_t * locus, msa_t * msa)
{
  /* convert alignments from ASCII character representation to state nucleotide
     state representation and store them in the locus structure. If cleandata
     is specified, remove all sites that contain an ambiguity at some sequence
  */

  set_tipstates(locus, msa, pll_map_nt);

}
