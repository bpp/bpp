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

const static int rate_matrices = 1;

static double pj_optimum = 0.3;

static stree_t * load_tree(void)
{
  stree_t * stree;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing species tree...");

  assert(opt_streenewick);

  stree = stree_parse_newick_string(opt_streenewick);

  if (!stree)
    fatal("Error while reading species tree");

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

  fprintf(stdout, "\nCurrent Pjump:    ");
  for (i = 0; i < PROP_COUNT + (opt_est_locusrate == MUTRATE_ESTIMATE); ++i)
    fprintf(stdout, " %8.5f", pjump[i]);
  fprintf(stdout, "\n");

  fprintf(stdout, "Current finetune: ");
  fprintf(stdout, " %8.5f", opt_finetune_gtage);
  fprintf(stdout, " %8.5f", opt_finetune_gtspr);
  fprintf(stdout, " %8.5f", opt_finetune_theta);
  fprintf(stdout, " %8.5f", opt_finetune_tau);
  fprintf(stdout, " %8.5f", opt_finetune_mix);
  if (opt_est_locusrate == MUTRATE_ESTIMATE)
    fprintf(stdout, " %8.5f\n", opt_finetune_locusrate);
  else
    fprintf(stdout, "\n");

  reset_finetune_onestep(pjump+0,&opt_finetune_gtage);
  reset_finetune_onestep(pjump+1,&opt_finetune_gtspr);
  reset_finetune_onestep(pjump+2,&opt_finetune_theta);
  reset_finetune_onestep(pjump+3,&opt_finetune_tau);
  reset_finetune_onestep(pjump+4,&opt_finetune_mix);

  if (opt_est_locusrate == MUTRATE_ESTIMATE)
    reset_finetune_onestep(pjump+5,&opt_finetune_locusrate);

  fprintf(stdout, "New finetune:     ");
  fprintf(stdout, " %8.5f", opt_finetune_gtage);
  fprintf(stdout, " %8.5f", opt_finetune_gtspr);
  fprintf(stdout, " %8.5f", opt_finetune_theta);
  fprintf(stdout, " %8.5f", opt_finetune_tau);
  fprintf(stdout, " %8.5f", opt_finetune_mix);
  if (opt_est_locusrate == MUTRATE_ESTIMATE)
    fprintf(stdout, " %8.5f\n", opt_finetune_locusrate);
  else
    fprintf(stdout, "\n");
}

static char * cb_serialize_branch(const snode_t * node)
{
  char * s = NULL;

  /* inner node */
  if (node->left)
  {
    if (node->parent)
    {
      if (node->theta > 0)
        xasprintf(&s, " #%f: %f", node->theta, node->parent->tau - node->tau);
      else
        xasprintf(&s, ": %f", node->parent->tau - node->tau);
    }
    else
    {
      if (node->theta > 0)
        xasprintf(&s, " #%f", node->theta);
    }
      
  }
  else
  {
    if (node->theta > 0)
      xasprintf(&s, "%s #%f: %f",
               node->label, node->theta, node->parent->tau - node->tau);
    else
      xasprintf(&s, "%s: %f", node->label, node->parent->tau - node->tau);
  }
    

  return s;
}

static void mcmc_printheader(FILE * fp, stree_t * stree)
{
  int print_labels = 1;
  unsigned int i;

  if (opt_method == METHOD_10)          /* species delimitation */
    fprintf(fp, "Gen\tnp\ttree");
  else
    fprintf(fp, "Gen");

  /* TODO: Account for integrated out theta */
  /* TODO: If number of species > 10 do not print labels */

  if (stree->tip_count > 10)
    print_labels = 0;

  /* 1. Print thetas */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->theta >= 0)
    {
      if (print_labels)
        fprintf(fp, "\ttheta_%d%s", i+1, stree->nodes[i]->label);
      else
        fprintf(fp, "\ttheta_%d", i+1);
    }

  /* 2. Print taus for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau)
    {
      if (print_labels)
        fprintf(fp, "\ttau_%d%s", i+1, stree->nodes[i]->label);
      else
        fprintf(fp, "\ttau_%d", i+1);
    }
  }

  /* 3. Print log likelihood */
  fprintf(fp, "\tlnL\n"); 
}

static void mcmc_printinitial(FILE * fp, stree_t * stree)
{
  char * newick = stree_export_newick(stree->root, cb_serialize_branch);
  fprintf(fp, "%s\n", newick);
  free(newick);
}

static void mcmc_logsample(FILE * fp,
                           int step,
                           stree_t * stree,
                           gtree_t ** gtree,
                           long dparam_count)
{
  unsigned int i;
  double logl = 0;

  if (opt_method == METHOD_01)          /* species tree inference */
  {
    char * newick = stree_export_newick(stree->root, cb_serialize_branch);
    fprintf(fp, "%s\n", newick);
    free(newick);
    return;
  }

  fprintf(fp, "%d", step);

  if  (opt_method == METHOD_10)         /* species delimitation */
  {
    fprintf(fp, "\t%ld", dparam_count);
    fprintf(fp, "\t%s", delimitation_getparam_string());
  }

  /* 1. Print thetas */

  /* first print thetas for tips */
  for (i = 0; i < stree->tip_count; ++i)
    if (stree->nodes[i]->theta >= 0)
      fprintf(fp, "\t%.5g", stree->nodes[i]->theta);

  /* then for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    if (stree->nodes[i]->theta >= 0)
      fprintf(fp, "\t%.5g", stree->nodes[i]->theta);

  /* 2. Print taus for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    if (stree->nodes[i]->tau)
      fprintf(fp, "\t%.5g", stree->nodes[i]->tau);

  /* print log-likelihood */
  for (i = 0; i < stree->locus_count; ++i)
    logl += gtree[i]->logl;

  fprintf(fp, "\t%.3f\n", logl);
}

static FILE * resume(stree_t ** ptr_stree,
                     gtree_t *** ptr_gtree,
                     locus_t *** ptr_locus,
                     double ** ptr_pjump,
                     unsigned long * ptr_curstep,
                     long * ptr_ft_round,

                     long * ptr_dparam_count,
                     long * ptr_ft_round_rj,
                     double * ptr_pjump_rj,

                     long * ptr_ft_round_spr,
                     long * ptr_pjump_slider,
                     double * ptr_mean_logl,
                     double * ptr_mean_root_age,
                     double * ptr_mean_root_theta,
                     stree_t ** ptr_sclone, 
                     gtree_t *** ptr_gclones)
{
  long i,j;
  FILE * fp_mcmc;
  long mcmc_offset;

  if (sizeof(BYTE) != 1)
    fatal("Checkpoint does not work on systems with sizeof(char) <> 1");

  /* load data from checkpoint file */
  checkpoint_load(ptr_gtree,
                  ptr_locus,
                  ptr_stree,
                  ptr_pjump,
                  ptr_curstep,
                  ptr_ft_round,
                  &mcmc_offset,
                  ptr_dparam_count,
                  ptr_ft_round_rj,
                  ptr_pjump_rj,
                  ptr_ft_round_spr,
                  ptr_pjump_slider,
                  ptr_mean_logl,
                  ptr_mean_root_age,
                  ptr_mean_root_theta);

  /* truncate MCMC file to specific offset */
  checkpoint_truncate(mcmc_offset);

  gtree_t ** gtree = *ptr_gtree;
  stree_t  * stree = *ptr_stree;

  gtree_alloc_internals(gtree,opt_locus_count);
  reset_gene_leaves_count(stree);
  stree_reset_pptable(stree);


  /* compute MSC density */
  for (i = 0; i < opt_locus_count; ++i)
    gtree[i]->logpr = gtree_logprob(stree,i);

  for (i = 0; i < opt_locus_count; ++i)
    printf("Gene tree %ld - logl: %f   logp: %f\n", i, gtree[i]->logl, gtree[i]->logpr);

  /* set old_pop to NULL */
  for (i = 0; i < opt_locus_count; ++i)
    for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      gtree[i]->nodes[j]->old_pop = NULL;

  for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
    stree->nodes[j]->mark = 0;


  /* set method */
  if (!opt_stree && !opt_delimit)
    opt_method = METHOD_00;
  else if (!opt_stree)
    opt_method = METHOD_10;
  else if (!opt_delimit)
    opt_method = METHOD_01;
  else
    fatal("Method 11 not yet implemented");

  /* open truncated MCMC file for appending */
  if (!(fp_mcmc = fopen(opt_mcmcfile, "a")))
    fatal("Cannot open file %s for appending...");

  if (opt_method == METHOD_01)
  {
    *ptr_sclone = stree_clone_init(stree);
    *ptr_gclones = (gtree_t **)xmalloc((size_t)opt_locus_count*sizeof(gtree_t *));
    for (i = 0; i < opt_locus_count; ++i)
      (*ptr_gclones)[i] = gtree_clone_init(gtree[i], *ptr_sclone);
  }
  if (opt_method == METHOD_10)          /* species delimitation */
  {
    /* quite ugly hack to resume species delimitation from a checkpoint.

     The idea is to:
     1) Store the taus in old_tau (0s in tau indicate a collapsed species)
     2) Call delimitations_init which computes all delimitation models,
        does all allocations, but has the side-effect of resetting the
        taus to a random delimitation
     3) Restore taus from old_tau
     4) Find delimittion string and return index using delimit_getindex()
     5) Set delimitation index using delimit_setindex()
    */
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      stree->nodes[i]->old_tau = stree->nodes[i]->tau;

    long dmodels_count = delimitations_init(stree);
    printf("Number of delimitation models: %ld\n", dmodels_count);
    rj_init(gtree,stree,opt_locus_count);

    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      stree->nodes[i]->tau = stree->nodes[i]->old_tau;
    
    long dindex = delimit_getindex(stree);
    delimit_setindex(dindex);
  }
    

  return fp_mcmc;
}

/* initialize everything - species tree, gene trees, locus structures etc.
   NOTE: *ALL* parameters of this function are output parameters, therefore
   do not concentrate on them when reading this function - they are filled
   at the end of the routine */
static FILE * init(stree_t ** ptr_stree,
                   gtree_t *** ptr_gtree,
                   locus_t *** ptr_locus,
                   double ** ptr_pjump,
                   unsigned long * ptr_curstep,
                   long * ptr_ft_round,
                   long * ptr_dparam_count,
                   long * ptr_ft_round_rj,
                   double * ptr_pjump_rj,
                   long * ptr_ft_round_spr,
                   long * ptr_pjump_slider,
                   double * ptr_mean_logl,
                   double * ptr_mean_root_age,
                   double * ptr_mean_root_theta,
                   stree_t ** ptr_sclone, 
                   gtree_t *** ptr_gclones)
{
  long i,j;
  long msa_count;
  double logl,logpr;
  double logl_sum = 0;
  double logpr_sum = 0;
  double * pjump;
  FILE * fp_mcmc;
  stree_t * stree;
  msa_t ** msa_list;
  gtree_t ** gtree;
  locus_t ** locus;

  /* method 10 specific variables */
  long dparam_count = 0;

  /* method 01 specific variables */
  stree_t * sclone = NULL;
  gtree_t ** gclones = NULL;
  
  /* load species tree */
  stree = load_tree();
  printf(" Done\n");

  /* parse the phylip file */

  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  assert(fd);

  printf("Parsing phylip file...");
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);
  printf(" Done\n");

  phylip_close(fd);

  /* set global variable with number of loci, if not set */
  if (!opt_locus_count)
    opt_locus_count = msa_count;

  /* remove ambiguous sites */
  if (opt_cleandata)
  {
    printf("Removing sites containing ambiguous characters...");
    for (i = 0; i < msa_count; ++i)
      if (!msa_remove_ambiguous(msa_list[i]))
        fatal("All sites in locus %d contain ambiguous characters",i);
    printf(" Done\n");
  }
  else
  {
    for (i = 0; i < msa_count; ++i)
      msa_count_ambiguous_sites(msa_list[i], pll_map_amb);
  }

  /* compress it */
  unsigned int ** weights = (unsigned int **)xmalloc(msa_count *
                                                     sizeof(unsigned int *));
  for (i = 0; i < msa_count; ++i)
  {
    msa_list[i]->original_length = msa_list[i]->length;
    weights[i] = compress_site_patterns(msa_list[i]->sequence,
                                        pll_map_nt,
                                        msa_list[i]->count,
                                        &(msa_list[i]->length),
                                        COMPRESS_JC69);
  }
  msa_summary(msa_list,msa_count);

  #if 0
  /* print the alignments */
  for (i = 0; i < msa_count; ++i)
    msa_print(msa_list[i]);
  #endif

  /* parse map file */
  printf("Parsing map file...");
  list_t * map_list = yy_parse_map(opt_mapfile);
  printf(" Done\n");
  #if 0
  maplist_print(map_list);
  #endif

  if (!(fp_mcmc = fopen(opt_mcmcfile, "w")))
    fatal("Cannot open file %s for writing...");

  /* TODO: PLACE DIPLOID CODE HERE */
  /* mapping from A2 -> A3 if diploid sequences used */
  unsigned long ** mapping = NULL;
  unsigned long ** resolution_count = NULL;
  int * unphased_length = NULL;

  if (opt_diploid)
  {
    /* store length of alignment A1 */
    unphased_length = (int *)xmalloc((size_t)msa_count * sizeof(int));
    for (i = 0; i < msa_count; ++i)
      unphased_length[i] = msa_list[i]->length;

    /* compute and replace msa_list with alignments A3. resolution_count
       contains the number of resolved sites in A2 for each site in A1,
       i.e. resolution_count[0][3] contains the number of resolved sites in A2
       for the fourth site of locus 0 */
    resolution_count = diploid_resolve(stree,
                                       msa_list,
                                       map_list,
                                       weights,
                                       msa_count,
                                       pll_map_nt);

    for (i = 0; i < msa_count; ++i)
      msa_print(msa_list[i]);

    /* TODO: KEEP WEIGHTS */
    //for (i = 0; i < msa_count; ++i) free(weights[i]);

    mapping = (unsigned long **)xmalloc((size_t)msa_count *
                                        sizeof(unsigned long *));

    for (i = 0; i < msa_count; ++i)
    {
      /* compress again for JC69 and get mappings */
      mapping[i] = compress_site_patterns_diploid(msa_list[i]->sequence,
                                                  pll_map_nt,
                                                  msa_list[i]->count,
                                                  &(msa_list[i]->length),
                                                  COMPRESS_JC69);
    }
  }

  if (opt_method == METHOD_10)          /* species delimitation */
  {
    long dmodels_count = delimitations_init(stree);
    printf("Number of delimitation models: %ld\n", dmodels_count);
  }

  /* initialize species tree (tau + theta) */
  stree_init(stree,msa_list,map_list,msa_count);
  stree_show_pptable(stree);

  double * lrate = (double *)xmalloc((size_t)opt_locus_count*sizeof(long));
  for (i = 0; i < opt_locus_count; ++i)
    lrate[i] = 1;

  if (opt_est_locusrate == MUTRATE_ESTIMATE)
  {
    double mean = 0;
    for (i = 0; i < opt_locus_count; ++i)
    {
      lrate[i] = 0.8 + 0.4*legacy_rndu();
      mean += lrate[i];
    }

    mean /= opt_locus_count;

    for (i = 0; i < opt_locus_count; ++i)
      lrate[i] /= mean;

    for (i = 0; i < opt_locus_count; ++i)
      printf("locusrate %ld: %f\n", i, lrate[i]);
  }
  else if (opt_est_locusrate == MUTRATE_FROMFILE)
  {
    parsefile_locusrates(lrate);

    double mean = 0;
    for (i = 0; i < opt_locus_count; ++i)
      mean += lrate[i];

    mean /= opt_locus_count;
      
    for (i = 0; i < opt_locus_count; ++i)
      lrate[i] /= mean;

    opt_est_locusrate = 0;
  }

  /* TODO CALL HERE */
  /* We must first link tip sequences (gene tips) to populations */
  if (opt_method == METHOD_10)          /* species delimitation */
    stree_rootdist(stree,map_list,msa_list,weights);

  gtree = gtree_init(stree,msa_list,map_list,msa_count);

  /* the below two lines are specific to method 01 and they generate
     space for cloning the species and gene trees */
  if (opt_method == METHOD_01)          /* species tree inference */
  {
    sclone = stree_clone_init(stree);
    gclones = (gtree_t **)xmalloc((size_t)msa_count*sizeof(gtree_t *));
    for (i = 0; i < msa_count; ++i)
      gclones[i] = gtree_clone_init(gtree[i], sclone);
  }

  locus = (locus_t **)xcalloc((size_t)msa_count, sizeof(locus_t *));

  /* Check that only first 32 bits of opt_arch are used */
  assert(opt_arch < (1l << 32)-1);

  gtree_update_branch_lengths(gtree, msa_count);
  for (i = 0; i < msa_count; ++i)
  {
    msa_t * msa = msa_list[i];
    double frequencies[4] = {0.25, 0.25, 0.25, 0.25};
    unsigned int pmatrix_count = gtree[i]->edge_count;

    /* if species tree inference or locusrate enabled, activate twice as many
       transition probability matrices */
    if (opt_method == METHOD_01 || opt_est_locusrate == MUTRATE_ESTIMATE)
      pmatrix_count *= 2;               /* double to account for cloned */

    /* TODO: In the future we can allocate double amount of p-matrices
       for the other methods as well in order to speedup rollback when
       rejecting proposals */

    /* create the locus structure */
    locus[i] = locus_create(gtree[i]->tip_count,        /* # tip sequence */
                            2*gtree[i]->inner_count,    /* # CLV vectors */
                            4,                          /* # states */
                            msa->length,                /* sequence length */
                            rate_matrices,              /* subst matrices (1) */
                            pmatrix_count,              /* # prob matrices */
                            1,                          /* # rate categories */
                            0,                          /* # scale buffers */
                            (unsigned int)opt_arch);    /* attributes */

    /* set frequencies for model with index 0 */
    pll_set_frequencies(locus[i],0,frequencies);

    if (opt_diploid)
    {
      for (j = 0; j < (long)(stree->tip_count); ++j)
        if (stree->nodes[j]->diploid)
        {
          locus[i]->diploid = 1;
          break;
        }
    }

    //if (opt_est_locusrate)
    //{
      /* TODO with more complex mixture models where rate_matrices > 1 we need
         to revisit this */
      assert(rate_matrices == 1);
      pll_set_mut_rates(locus[i],lrate+i);
    //}

    /* set pattern weights and free the weights array */
    if (locus[i]->diploid)
    {
      /* TODO: 1) pattern_weights_sum is not updated here, but it is not used in
         the program, perhaps remove.
         2) pattern_weights is allocated in locus_create with a size msa->length
            equal to length of A3, but in reality we only need |A1| storage
            space. Free and reallocate here  */

      locus[i]->diploid_mapping = mapping[i];
      locus[i]->diploid_resolution_count = resolution_count[i];
      /* since PLL does not support diploid sequences we make a small hack */
      memcpy(locus[i]->pattern_weights, weights[i], unphased_length[i]);
      free(weights[i]);
      locus[i]->likelihood_vector = (double *)xmalloc((size_t)(msa->length) *
                                                      sizeof(double));
      locus[i]->unphased_length = unphased_length[i];
    }
    else
    {
      pll_set_pattern_weights(locus[i], weights[i]);
      free(weights[i]);
    }

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
    logpr_sum += logpr;
  }

  /* deallocate unnecessary arrays */
  free(lrate);
  if (opt_diploid)
  {
    free(mapping);
    free(unphased_length);
    free(resolution_count);
  }

  printf("\nInitial MSC density and log-likelihood of observing data:\n");
  printf("log-P0 = %f   log-L0 = %f\n\n", logpr_sum, logl_sum);

  /* free weights array */
  free(weights);

  if (opt_method == METHOD_10)          /* species delimitation */
    rj_init(gtree,stree,msa_count);

  /* initialize pjump and finetune rounds */
  if (opt_est_locusrate)
    pjump = (double *)xcalloc(PROP_COUNT+1, sizeof(double));
  else
    pjump = (double *)xcalloc(PROP_COUNT, sizeof(double));

  /* TODO: Method 10 has a commented call to 'delimit_resetpriors()' */
  //delimit_resetpriors();

  /* if method 00 or 01 print corresponding header line in MCMC file */
  if (opt_method == METHOD_01)
    mcmc_printinitial(fp_mcmc,stree);
  else
    mcmc_printheader(fp_mcmc,stree);

  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);

  /* TODO: This is only for species delimitation */
  if (opt_method == METHOD_10)          /* species delimitation */
  {
    dparam_count = 0;
    for (i = 0; i < (int)(stree->tip_count + stree->inner_count); ++i)
    {
      if (stree->nodes[i]->theta > 0)
        dparam_count++;
      if (stree->nodes[i]->tau > 0)
        dparam_count++;
    }
  }

  *ptr_stree = stree;
  *ptr_gtree = gtree;
  *ptr_locus = locus;
  *ptr_pjump = pjump;

  *ptr_curstep = 0;
  *ptr_ft_round = 0;

  /* species delimitation relevant */
  *ptr_ft_round_rj = 0;
  *ptr_dparam_count = dparam_count;
  *ptr_pjump_rj = 0;

  /* species tree inference relevant */
  *ptr_ft_round_spr = 0;
  *ptr_pjump_slider = 0;
  *ptr_mean_logl = 0;
  *ptr_mean_root_age = 0;
  *ptr_mean_root_theta = 0;

  *ptr_sclone = sclone;
  *ptr_gclones = gclones;

  /* deallocate maplist */
  list_clear(map_list,map_dealloc);
  free(map_list);

  /* deallocate alignments */
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  return fp_mcmc;

}

void cmd_run()
{
  /* common variables for all methods */
  long i,j;
  long ft_round;
  double logl_sum = 0;
  double * pjump;
  FILE * fp_mcmc;
  stree_t * stree;
  gtree_t ** gtree;
  locus_t ** locus;

  /* method 10 specific variables */
  long dparam_count = 0;
  long ft_round_rj;
  double pjump_rj = 0;

  /* method 01 specific variables */
  long ft_round_spr = 0;
  long pjump_slider;
  long printk = opt_samplefreq * opt_samples;
  double mean_logl = 0;
  double mean_root_age = 0;
  double mean_root_theta = 0;
  stree_t * sclone;
  gtree_t ** gclones;

  unsigned long curstep = 0;


  if (opt_resume)
    fp_mcmc = resume(&stree,
                     &gtree,
                     &locus,
                     &pjump,
                     &curstep,
                     &ft_round,
                     &dparam_count,
                     &ft_round_rj,
                     &pjump_rj,
                     &ft_round_spr,
                     &pjump_slider,
                     &mean_logl,
                     &mean_root_age,
                     &mean_root_theta,
                     &sclone, 
                     &gclones);
  else
    fp_mcmc = init(&stree,
                   &gtree,
                   &locus,
                   &pjump,
                   &curstep,
                   &ft_round,
                   &dparam_count,
                   &ft_round_rj,
                   &pjump_rj,
                   &ft_round_spr,
                   &pjump_slider,
                   &mean_logl,
                   &mean_root_age,
                   &mean_root_theta,
                   &sclone, 
                   &gclones);

  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);

  /* start of MCMC loop */
  for (i = curstep-opt_burnin; i < opt_samples*opt_samplefreq; ++i)
  {
    /* update progress bar */
    if (!opt_quiet)
      progress_update(curstep);

    /* reset finetune parameters */
    if (i == 0 || (opt_finetune_reset && opt_burnin >= 200 && i < 0 &&
                   ft_round >= 100 && i%(opt_burnin/4)==0))
    {
      int pjump_size = PROP_COUNT + (opt_est_locusrate == 1);
      if (opt_finetune_reset && opt_burnin >= 200)
        reset_finetune(pjump);
      for (j = 0; j < pjump_size; ++j)
        pjump[j] = 0;

      /* reset pjump and number of steps since last finetune reset to zero */
      ft_round = 0;
      memset(pjump,0,pjump_size*sizeof(double));

      if (opt_method == METHOD_10)      /* species delimitation */
      {
        pjump_rj = 0;
        ft_round_rj = 0;
      }
      if (opt_method == METHOD_01)      /* species tree inference */
      {
        ft_round_spr = 0;
        pjump_slider = 0;

        mean_logl = 0;
        mean_root_theta = 0;
        mean_root_age = 0;
      }
    }
    
    ++ft_round;

    /* perform proposals sequentially */   
    double ratio;

    /* propose delimitation through merging/splitting of nodes */
    if (opt_method == METHOD_10)        /* species delimitation */
    {
      if (legacy_rndu() < 0.5)
        j = prop_split(gtree,stree,locus,0.5,&dparam_count);
      else
        j = prop_join(gtree,stree,locus,0.5,&dparam_count);

      if (j != 2)
      {
        ft_round_rj++;
        pjump_rj += j;
      }
    }

    /* propose species tree topology using SPR */
    if (opt_method == METHOD_01)        /* species tree inference */
    {
      if (legacy_rndu() > 0)   /* bpp4 compatible results (RNG to next state) */
      {
        if (stree_propose_spr(&stree, &gtree, &sclone, &gclones, locus))
        {
          /* accepted */

          /* swap the pointers of species tree and gene tree list with cloned */
          SWAP(stree,sclone);
          SWAP(gtree,gclones);

          stree_label(stree);

          pjump_slider++;
        }

        ft_round_spr++;
      }
    }

    /* propose gene tree ages */
    ratio = gtree_propose_ages(locus, gtree, stree);
    pjump[0] = (pjump[0]*(ft_round-1) + ratio) / (double)ft_round;

    /* propose gene tree topologies using SPR */
    ratio = gtree_propose_spr(locus,gtree,stree);
    pjump[1] = (pjump[1]*(ft_round-1) + ratio) / (double)ft_round;

    /* propose population sizes on species tree */
    ratio = stree_propose_theta(gtree,stree);
    pjump[2] = (pjump[2]*(ft_round-1) + ratio) / (double)ft_round;

    /* propose species tree taus */
    ratio = stree_propose_tau(gtree,stree,locus);
    pjump[3] = (pjump[3]*(ft_round-1) + ratio) / (double)ft_round;

    /* mixing step */
    ratio = proposal_mixing(gtree,stree,locus);
    pjump[4] = (pjump[4]*(ft_round-1) + ratio) / (double)ft_round;

    if (opt_est_locusrate)
    {
      ratio = prop_locusrate(gtree,stree,locus);
      pjump[5] = (pjump[5]*(ft_round-1) + ratio) / (double)ft_round;
    }

    /* log sample into file (dparam_count is only used in method 10) */
    if (i >= 0 && (i+1)%opt_samplefreq == 0)
      mcmc_logsample(fp_mcmc,i+1,stree,gtree,dparam_count);

    if (opt_method == METHOD_01)        /* species tree inference */
    {
      /* update stats for printing on screen */
      mean_root_theta = (mean_root_theta*(ft_round-1) + stree->root->theta)/ft_round;
      mean_root_age = (mean_root_age*(ft_round-1) + stree->root->tau)/ft_round;
      for (logl_sum = 0, j = 0; j < opt_locus_count; ++j)
        logl_sum += gtree[j]->logl;
      mean_logl = (mean_logl * (ft_round-1) + logl_sum / opt_bfbeta)/ft_round;
    }

    /* TODO: print on screen */

    if (opt_method == METHOD_01)        /* species tree inference */
    {
      if (printk <= 500 || (i+1) % (printk / 200) == 0)
      {
        printf("\r%3.0f%%", (i + 1.499) / printk * 100.);
        for (j = 0; j < 5 + (opt_est_locusrate == 1); ++j)
          printf(" %4.2f", pjump[j]);
        printf(" ");

        printf(" %5.3f  %6.4f  %6.4f", ft_round_spr ? 
                                         (double)pjump_slider / ft_round_spr : 0.,
                                       mean_root_theta,
                                       mean_root_age);

        if (opt_usedata)
          printf(" %8.4f", mean_logl);

        if (printk >= 50 && (i+1) % (printk / 20) == 0)
          printf("\n");
      }
    }

    curstep++;

    if (opt_checkpoint)
    {
      if (((long)curstep == opt_checkpoint_initial) ||
          (opt_checkpoint_step && ((long)curstep > opt_checkpoint_initial) &&
           (((long)curstep-opt_checkpoint_initial) % opt_checkpoint_step == 0)))
      {
        checkpoint_dump(stree,
                        gtree,
                        locus,
                        pjump,
                        curstep,
                        ft_round,
                        ftell(fp_mcmc),
                        dparam_count,
                        ft_round_rj,
                        pjump_rj,
                        ft_round_spr,
                        pjump_slider,
                        mean_logl,
                        mean_root_age,
                        mean_root_theta);
      }
    }

  }

  progress_done();

  free(pjump);

  fclose(fp_mcmc);

  /* print summary using the MCMC file */
  if (opt_method == METHOD_10)          /* species delimitation */
  {
    delimit_summary(stree);
    delimitations_fini();
    rj_fini();
  }

  for (i = 0; i < opt_locus_count; ++i)
    locus_destroy(locus[i]);
  free(locus);

  /* deallocate gene trees */
  for (i = 0; i < opt_locus_count; ++i)
    gtree_destroy(gtree[i],NULL);
  free(gtree);

  /* if species tree inference, deallocate cloned gene trees */
  if (opt_method == METHOD_01)          /* species tree inference */
  {
    for (i = 0; i < opt_locus_count; ++i)
      gtree_destroy(gclones[i],NULL);
    free(gclones);
  }

  gtree_fini(opt_locus_count);

  if (opt_method == METHOD_00)
    allfixed_summary(stree);

  unsigned int species_count = 0;
  char ** species_names = NULL;
  if (opt_method == METHOD_01)
  {
    /* order of species */
    species_names = (char **)xmalloc(stree->tip_count * sizeof(char *));
    species_count = stree->tip_count;

    for (i = 0; i < (long)(stree->tip_count); ++i)
      species_names[i] = xstrdup(stree->nodes[i]->label);
  }

  /* deallocate tree */
  stree_destroy(stree,NULL);
  if (opt_method == METHOD_01)          /* species tree inference */
    stree_destroy(sclone,NULL);         /* destroy cloned species tree */
    
  stree_fini();

  if (opt_diploid)
    free(opt_diploid);

  /* summary for method 01 */
  if (opt_method == METHOD_01)          /* species tree inference */
  {
    assert(species_count > 0);

    stree_summary(species_names,species_count);

    /* cleanup */
    for (i = 0; i < species_count; ++i)
      free(species_names[i]);
    free(species_names);
  }
}
