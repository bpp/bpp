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

#if 0
static int cb_full_traversal(gnode_t * node)
{
  if (!node->left)
    return 0;

  return 1;
}
#endif

void cmd_a00()
{
  int i,j;
  int msa_count;
  double logl,logpr;
  double logl_sum = 0;
  msa_t ** msa_list;

  if (opt_mcmc_steps <= 0)
    fatal("--mcmc_steps must be a positive integer greater than zero");

  if (opt_mcmc_burnin < 1 || opt_mcmc_burnin > opt_mcmc_steps)
    fatal("--mcmc_burnin must be a positive integer smaller or equal to --mcmc_steps");

  /* load species tree */
  stree_t * stree = load_tree();

  /* init random number generator */
  srand48(opt_seed);

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

  /* parse map file */
  printf("Parsing map file...\n");
  list_t * map_list = yy_parse_map(opt_mapfile);
  maplist_print(map_list);


  /* call MCMC */

  stree_init(stree,msa_list,map_list,msa_count);

  gtree_t ** gtree = gtree_init(stree,msa_list,map_list,msa_count);

  locus_t ** locus = (locus_t **)xcalloc(msa_count, sizeof(locus_t *));

  gtree_update_branch_lengths(gtree, msa_count);
  for (i = 0; i < msa_count; ++i)
  {
    msa_t * msa = msa_list[i];
    double frequencies[4] = {0.25, 0.25, 0.25, 0.25};

    /* create the locus structure */
    locus[i] = locus_create(gtree[i]->tip_count,        /* # tip sequence */
                            2*gtree[i]->inner_count,    /* # CLV vectors */
                            4,                          /* # states */
                            msa->length,                /* sequence length */
                            1,                          /* # subst matrices */
                            gtree[i]->edge_count,       /* # prob matrices */
                            1,                          /* # rate categories */
                            0,                          /* # scale buffers */
                            PLL_ATTRIB_ARCH_SSE);       /* attributes */

    /* set frequencies for model with index 0 */
    pll_set_frequencies(locus[i],0,frequencies);

    /* set tip sequences */
    for (j = 0; j < (int)(gtree[i]->tip_count); ++j)
      pll_set_tip_states(locus[i], j, pll_map_nt, msa_list[i]->sequence[j]);

    /* compute the conditional probabilities for each inner node */
    locus_update_matrices_jc69(locus[i],gtree[i]->nodes,gtree[i]->edge_count);
    locus_update_partials(locus[i],
                          gtree[i]->nodes+gtree[i]->tip_count,
                          gtree[i]->inner_count);

    /* optionally, show root CLV 

    pll_show_clv(locus[i], gtree[i]->root->clv_index, PLL_SCALE_BUFFER_NONE, 9);

    */

    /* now that we computed the CLVs, calculate the log-likelihood for the
       current gene tree */
    unsigned int param_indices[1] = {0};
    logl = locus_root_loglikelihood(locus[i],
                                           gtree[i]->root,
                                           param_indices,
                                           NULL);
    logl_sum += logl;

    /* store current log-likelihood in each gene tree structure */
    gtree[i]->logl = logl;
    logpr = gtree_logprob(stree,i);
    gtree[i]->logpr = logpr;
  }

  for (i = 0; i < msa_count; ++i)
    gtree_propose_ages(locus[i], gtree[i], stree, i);


  for (i = 0; i < msa_count; ++i)
    locus_destroy(locus[i]);
  free(locus);

  /* deallocate gene trees */
  for (i = 0; i < msa_count; ++i)
    gtree_destroy(gtree[i],NULL);
  free(gtree);
  gtree_fini(msa_count);

  /* deallocate tree */
  stree_destroy(stree,NULL);

  /* deallocate alignments */
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  /* deallocate maplist */
  list_clear(map_list,map_dealloc);
  free(map_list);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}
