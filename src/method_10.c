/*
    Copyright (C) 2016-2017 Tomas Flouri and Ziheng Yang

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

static stree_t * load_tree(void)
{
  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  stree_t * stree = stree_parse_newick(opt_streefile);
  if (!stree)
    fatal("Error while reading tree file %s\n.", opt_streefile);

  return stree;
}

void cmd_a10()
{
  int i;
  int msa_count;
//  double logl,logpr;
//  double logl_sum = 0;
  msa_t ** msa_list;
  FILE * fp_mcmc;

  if (opt_samples < 1)
    fatal("--samples must be a positive integer greater than zero");

  if (opt_burnin < 0)
    fatal("--burnin must be a positive integer");

  /* load species tree */
  stree_t * stree = load_tree();

  /* parse the phylip file */
  printf("Parsed tree\n");

  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  assert(fd);

  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);

  phylip_close(fd);

  /* print the alignments */
  for (i = 0; i < msa_count; ++i)
    msa_print(msa_list[i]);

  /* compress it */
  unsigned int ** weights = (unsigned int **)xmalloc(msa_count *
                                                     sizeof(unsigned int *));
  for (i = 0; i < msa_count; ++i)
  {
    int ol = msa_list[i]->length;
    weights[i] = compress_site_patterns(msa_list[i]->sequence,
                                        pll_map_nt,
                                        msa_list[i]->count,
                                        &(msa_list[i]->length),
                                        COMPRESS_JC69);
    printf("Locus %d: original length %d, after compression %d\n", i, ol, msa_list[i]->length);
  }

  /* parse map file */
  printf("Parsing map file...\n");
  list_t * map_list = yy_parse_map(opt_mapfile);
  maplist_print(map_list);

  if (!(fp_mcmc = fopen(opt_mcmcfile, "w")))
    fatal("Cannot open file %s for writing...");

  printf("Total species delimitations: %ld\n", delimitations_count(stree));
  delimitations_enumerate(stree);

//  histories();
  printf("METHOD 10\n");
}
