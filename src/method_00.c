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

static int cb_full_traversal(gnode_t * node)
{
  if (!node->left)
    return 0;

  return 1;
}

void cmd_a00()
{
  int i,j;
  int msa_count;
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

  /* debugging */
  stree_init(stree,msa_list,map_list,msa_count);

  gtree_t ** gtree = gtree_init(stree,msa_list,map_list,msa_count);

  locus_t ** locus = (locus_t **)xcalloc(msa_count, sizeof(locus_t *));

  for (i = 0; i < msa_count; ++i)
  {
    msa_t * msa = msa_list[i];
    double frequencies[4] = {0.25, 0.25, 0.25, 0.25};

    /* create the locus structure */
    locus[i] = locus_create(gtree[i]->tip_count,        /* # tip sequence */
                            gtree[i]->inner_count,      /* # CLV vectors */
                            4,                          /* # states */
                            msa->length,                /* sequence length */
                            1,                          /* # subst matrices */
                            2,                          /* # prob matrices */
                            1,                          /* # rate categories */
                            0,                          /* # scale buffers */
                            PLL_ATTRIB_ARCH_CPU);       /* attributes */

    /* set frequencies for model with index 0 */
    pll_set_frequencies(locus[i],0,frequencies);

    /* set tip sequences */
    for (j = 0; j < (int)(gtree[i]->tip_count); ++j)
      pll_set_tip_states(locus[i], j, pll_map_nt, msa_list[i]->sequence[j]);

    /* get a postorder traversal of inner nodes for gene tree i */
    unsigned int trav_size;
    gnode_t ** trav = (gnode_t **)xmalloc(gtree[i]->inner_count * sizeof(gnode_t *));
    gtree_traverse(gtree[i]->root,
                   TREE_TRAVERSE_POSTORDER,
                   cb_full_traversal,
                   trav,
                   &trav_size);


    /* compute the conditional probabilities for each inner node */
    for (j = 0; j < (int)trav_size; ++j)
    {
      double branch_length[2];
      unsigned int matrix_indices[2] = {0,1};
      unsigned int param_indices[1] = {0};
      double rates[1] = {1};
      
      branch_length[0] = trav[j]->left->length;
      branch_length[1] = trav[j]->right->length;

      /* update two transition probability matrices with indices 0 and 1,
         and use the substitution matrix and frequencies with index 0 
         (param_indices) */
      pll_core_update_pmatrix_4x4_jc69(locus[i]->pmatrix,
                                       locus[i]->states,
                                       locus[i]->rate_cats,
                                       rates,
                                       branch_length,
                                       matrix_indices,
                                       param_indices,
                                       2,
                                       PLL_ATTRIB_ARCH_CPU);

      /* optionally, show pmatrices 

      pll_show_pmatrix(locus[i], 0, 5);
      pll_show_pmatrix(locus[i], 1, 5);

      */

      /* now compute the conditional probabilities for the current node */
      pll_core_update_partial_ii(locus[i]->states,
                                 locus[i]->sites,
                                 locus[i]->rate_cats,
                                 locus[i]->clv[trav[j]->clv_index],
                                 NULL,
                                 locus[i]->clv[trav[j]->left->clv_index],
                                 locus[i]->clv[trav[j]->right->clv_index],
                                 locus[i]->pmatrix[0],
                                 locus[i]->pmatrix[1],
                                 NULL,
                                 NULL,
                                 PLL_ATTRIB_ARCH_CPU);
    
    }

    assert(gtree[i]->root == trav[gtree[i]->inner_count-1]);


    /* optionally, show root CLV 

    pll_show_clv(locus[0], gtree[0]->root->clv_index, PLL_SCALE_BUFFER_NONE, 9);

    */

    /* now that we computed the CLVs, calculate the log-likelihood for the
       current gene tree */
    unsigned int param_indices[1] = {0};
    double logl = pll_core_root_loglikelihood(locus[i]->states,
                                locus[i]->sites,
                                locus[i]->rate_cats,
                                locus[i]->clv[gtree[i]->root->clv_index],
                                NULL,
                                locus[i]->frequencies,
                                locus[i]->rate_weights,
                                locus[i]->pattern_weights,
                                param_indices,
                                NULL,
                                0);

    printf("logL gene tree %d : %f\n", i,logl);
    free(trav);
  }


  for (i = 0; i < msa_count; ++i)
    locus_destroy(locus[i]);
  free(locus);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");

  /* deallocate gene trees */
  for (i = 0; i < msa_count; ++i)
    gtree_destroy(gtree[i],NULL);
  free(gtree);

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
