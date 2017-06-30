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

#define PI  3.1415926535897932384626433832795
#define PROP_COUNT 5

static double pj_optimum = 0.3;

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

static void reset_finetune_onestep(double * pjump, double * param)
{
  double maxstep = 99;

  if (*pjump < 0.001)
    *param /= 100;
  else if (*pjump > 0.999)
    *param = MIN(maxstep, *param * 100);
  else
  {
    *param *= tan(PI/2*(*pjump)) / tan(PI/2*pj_optimum);
    *param = MIN(maxstep, *param);
  }
  
}

static void reset_finetune(double * pjump)
{
  int i;

  fprintf(stdout, "Current Pjump:    ");
  for (i = 0; i < PROP_COUNT; ++i)
    fprintf(stdout, " %8.5f", pjump[i]);
  fprintf(stdout, "\n");

  fprintf(stdout, "Current finetune: ");
  fprintf(stdout, " %8.5f", opt_finetune_gtage);
  fprintf(stdout, " %8.5f", opt_finetune_gtspr);
  fprintf(stdout, " %8.5f", opt_finetune_theta);
  fprintf(stdout, " %8.5f", opt_finetune_tau);
  fprintf(stdout, " %8.5f\n", opt_finetune_mix);

  reset_finetune_onestep(pjump+0,&opt_finetune_gtage);
  reset_finetune_onestep(pjump+1,&opt_finetune_gtspr);
  reset_finetune_onestep(pjump+2,&opt_finetune_theta);
  reset_finetune_onestep(pjump+3,&opt_finetune_tau);
  reset_finetune_onestep(pjump+4,&opt_finetune_mix);

  fprintf(stdout, "New finetune:     ");
  fprintf(stdout, " %8.5f", opt_finetune_gtage);
  fprintf(stdout, " %8.5f", opt_finetune_gtspr);
  fprintf(stdout, " %8.5f", opt_finetune_theta);
  fprintf(stdout, " %8.5f", opt_finetune_tau);
  fprintf(stdout, " %8.5f\n", opt_finetune_mix);
}

static void mcmc_printheader(FILE * fp, stree_t * stree)
{
  int print_labels = 1;
  unsigned int i,j=0;
  fprintf(fp, "Gen");

  /* TODO: Account for integrated out theta */
  /* TODO: If number of species > 10 do not print labels */

  if (stree->tip_count > 10)
    print_labels = 0;

  /* 1. Print thetas */

  /* first print thetas for tips */
  for (i = 0; i < stree->tip_count; ++i)
    if (stree->nodes[i]->theta >= 0)
    {
      if (print_labels)
        fprintf(fp, "\ttheta_%d%s", ++j, stree->nodes[i]->label);
      else
        fprintf(fp, "\ttheta_%d", ++j);
    }

  /* then for root */
  if (stree->root->theta >= 0)
  {
    if (print_labels)
      fprintf(fp, "\ttheta_%d%s", ++j, stree->root->label);
    else
      fprintf(fp, "\ttheta_%d", ++j);
  }

  /* and finally for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (stree->nodes[i] == stree->root) continue;

    if (stree->nodes[i]->theta >= 0)
    {
      if (print_labels)
        fprintf(fp, "\ttheta_%d%s", ++j, stree->nodes[i]->label);
      else
        fprintf(fp, "\ttheta_%d", ++j);
    }
  }

  /* 2. Print taus */
  j = stree->tip_count;

  /* first print tau for root */
  if (stree->root->tau)
  {
    if (print_labels)
      fprintf(fp, "\ttau_%d%s", ++j, stree->root->label);
    else
      fprintf(fp, "\ttau_%d", ++j);
  }

  /* and finally for the remaining inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (stree->nodes[i] == stree->root) continue;

    if (stree->nodes[i]->tau)
    {
      if (print_labels)
        fprintf(fp, "\ttau_%d%s", ++j, stree->nodes[i]->label);
      else
        fprintf(fp, "\ttau_%d", ++j);
    }
  }

  /* 3. Print log likelihood */
  fprintf(fp, "\tlnL\n"); 
}

static void mcmc_logsample(FILE * fp,
                           int step,
                           stree_t * stree,
                           gtree_t ** gtree)
{
  unsigned int i;
  double logl = 0;

  fprintf(fp, "%d", step);

  /* 1. Print thetas */

  /* first print thetas for tips */
  for (i = 0; i < stree->tip_count; ++i)
    if (stree->nodes[i]->theta >= 0)
      fprintf(fp, "\t%.5g", stree->nodes[i]->theta);

  /* then for root */
  if (stree->root->theta >= 0)
      fprintf(fp, "\t%.5g", stree->root->theta);

  /* and finally for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (stree->nodes[i] == stree->root) continue;

    if (stree->nodes[i]->theta >= 0)
      fprintf(fp, "\t%.5g", stree->nodes[i]->theta);
  }

  /* 2. Print taus */

  /* first print tau for root */
  if (stree->root->tau)
    fprintf(fp, "\t%.5g", stree->root->tau);

  /* and finally for the remaining inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (stree->nodes[i] == stree->root) continue;

    if (stree->nodes[i]->tau)
      fprintf(fp, "\t%.5g", stree->nodes[i]->tau);
  }

  /* print log-likelihood */
  for (i = 0; i < stree->locus_count; ++i)
    logl += gtree[i]->logl;

  fprintf(fp, "\t%.3f\n", logl);
}

void cmd_a00()
{
  int i,j;
  int msa_count;
  double logl,logpr;
  double logl_sum = 0;
  msa_t ** msa_list;
  FILE * fp_mcmc;

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

  if (!(fp_mcmc = fopen(opt_mcmcfile, "w")))
    fatal("Cannot open file %s for writing...");

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
                            PLL_ATTRIB_ARCH_CPU);       /* attributes */

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

  /* start of MCMC loop */

  long steps = opt_burnin + opt_samples*opt_samplefreq;

  double * pjump = (double *)xcalloc(PROP_COUNT, sizeof(double));
  long ft_round = 0;

  /* print header in mcmc file */
  mcmc_printheader(fp_mcmc,stree);

  /* start of MCMC loop */
  for (i = 0; i < steps; ++i)
  {
    /* TODO: STOPPING CRITERION FOR DEBUGGING */
    //if (i == 10001) assert(0);

    /* reset finetune parameters */
    if (i == opt_burnin || (opt_finetune_reset && opt_burnin >= 200 &&
        i < opt_burnin && ft_round >= 100 && (i%(opt_burnin/4)==0)))
    {
      reset_finetune(pjump);
      for (j = 0; j < PROP_COUNT; ++j)
        pjump[j] = 0;

      /* reset pjump and number of steps since last finetune reset to zero */
      ft_round = 0;
      memset(pjump,0,PROP_COUNT*sizeof(double));
    }

    ++ft_round;

    /* perform proposals sequentially */   
    double ratio;

    /* proposal on gene tree ages */
    ratio = gtree_propose_ages(locus, gtree, stree);
    pjump[0] = (pjump[0]*(ft_round-1) + ratio) / (double)ft_round;

    ratio = gtree_propose_spr(locus,gtree,stree);
    pjump[1] = (pjump[1]*(ft_round-1) + ratio) / (double)ft_round;

    ratio = stree_propose_theta(gtree,stree);
    pjump[2] = (pjump[2]*(ft_round-1) + ratio) / (double)ft_round;

    ratio = stree_propose_tau(gtree,stree,locus);
    pjump[3] = (pjump[3]*(ft_round-1) + ratio) / (double)ft_round;

    ratio = proposal_mixing(gtree,stree,locus);
    pjump[4] = (pjump[4]*(ft_round-1) + ratio) / (double)ft_round;

    /* log into file */
    if (opt_log_samples && i >= opt_burnin && 
        (i-opt_burnin+1) % opt_samplefreq == 0)
    {
      /* TODO: print parameters in output mcmc file */
      mcmc_logsample(fp_mcmc,i-opt_burnin+1,stree,gtree);
    }

    /* TODO: print on screen */
  }

  free(pjump);

  fclose(fp_mcmc);

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
  stree_fini();

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
