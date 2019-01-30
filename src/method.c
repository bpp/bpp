/*
    Copyright (C) 2016-2019 Tomas Flouri and Ziheng Yang

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

#define PI  3.1415926535897932384626433832795
#define PROP_COUNT 5

/* maximum number of theta/tau to output on screen during MCMC */
#define MAX_THETA_OUTPUT        3
#define MAX_TAU_OUTPUT          3

const static int rate_matrices = 1;
const static long thread_index_zero = 0;

static double pj_optimum = 0.3;
static thread_data_t td;

static stree_t * load_tree_or_network(void)
{
  stree_t * stree;

  printf("Using seed: %ld\n", opt_seed);
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

static void reset_finetune(double * pjump, double * pjump_phi)
{
  int i;
  int extra;
  
  extra = (opt_est_locusrate || opt_est_heredity);

  fprintf(stdout, "\nCurrent Pjump:    ");

  for (i = 0; i < PROP_COUNT + extra; ++i)
    fprintf(stdout, " %8.5f", pjump[i]);
  if (opt_msci)
    fprintf(stdout, " %8.5f", *pjump_phi);
  fprintf(stdout, "\n");

  fprintf(stdout, "Current finetune: ");
  fprintf(stdout, " %8.5f", opt_finetune_gtage);
  fprintf(stdout, " %8.5f", opt_finetune_gtspr);
  fprintf(stdout, " %8.5f", opt_finetune_theta);
  fprintf(stdout, " %8.5f", opt_finetune_tau);
  fprintf(stdout, " %8.5f", opt_finetune_mix);
  if (extra)
    fprintf(stdout, " %8.5f", opt_finetune_locusrate);
  if (opt_msci)
    fprintf(stdout, " %8.5f\n", opt_finetune_phi);
  else
    fprintf(stdout, "\n");


  reset_finetune_onestep(pjump+0,&opt_finetune_gtage);
  reset_finetune_onestep(pjump+1,&opt_finetune_gtspr);
  reset_finetune_onestep(pjump+2,&opt_finetune_theta);
  reset_finetune_onestep(pjump+3,&opt_finetune_tau);
  reset_finetune_onestep(pjump+4,&opt_finetune_mix);

  if (extra)
    reset_finetune_onestep(pjump+5,&opt_finetune_locusrate);
  if (opt_msci)
    reset_finetune_onestep(pjump_phi, &opt_finetune_phi);

  fprintf(stdout, "New finetune:     ");
  fprintf(stdout, " %8.5f", opt_finetune_gtage);
  fprintf(stdout, " %8.5f", opt_finetune_gtspr);
  fprintf(stdout, " %8.5f", opt_finetune_theta);
  fprintf(stdout, " %8.5f", opt_finetune_tau);
  fprintf(stdout, " %8.5f", opt_finetune_mix);
  if (extra)
    fprintf(stdout, " %8.5f", opt_finetune_locusrate);
  if (opt_msci)
    fprintf(stdout, " %8.5f\n", opt_finetune_phi);
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
      if (opt_est_theta && node->theta > 0)
        xasprintf(&s, " #%f: %f", node->theta, node->parent->tau - node->tau);
      else
        xasprintf(&s, ": %f", node->parent->tau - node->tau);
    }
    else
    {
      if (opt_est_theta && node->theta > 0)
        xasprintf(&s, " #%f", node->theta);
      else
      {
        /* replacement to keep GCC warning from showing up when executing
           the following:
           
           xasprintf(&s, "");
        */

        s = (char *)xmalloc(sizeof(char));
        *s = 0;

      }
    }
      
  }
  else
  {
    if (opt_est_theta && node->theta > 0)
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
  unsigned int snodes_total;
  
  if (opt_msci)
    snodes_total = stree->tip_count + stree->inner_count + stree->hybrid_count;
  else
    snodes_total = stree->tip_count + stree->inner_count;

  if (opt_method == METHOD_10)          /* species delimitation */
    fprintf(fp, "Gen\tnp\ttree");
  else
    fprintf(fp, "Gen");

  /* TODO: If number of species > 10 do not print labels */

  if (stree->tip_count > 10)
    print_labels = 0;

  /* 1. Print thetas */
  if (opt_est_theta)
  {
    for (i = 0; i < snodes_total; ++i)
    {
      /* TODO: Is the 'has_theta' check also necessary ? */
      if (stree->nodes[i]->theta >= 0)
      {
        if (print_labels)
          fprintf(fp, "\ttheta_%d%s", i+1, stree->nodes[i]->label);
        else
          fprintf(fp, "\ttheta_%d", i+1);
      }
    }
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

  if (opt_msci)
  {
    unsigned int offset=stree->tip_count+stree->inner_count;
    for (i = 0; i < stree->hybrid_count; ++i)
      fprintf(fp, "\tphi_%s", stree->nodes[offset+i]->hybrid->label);
  }

  /* 3. Print mutation rate for each locus */
  if (opt_est_locusrate && opt_print_locusrate)
  {
    for (i = 0; i < opt_locus_count; ++i)
      fprintf(fp, "\trate_L%d", i+1);
  }

  /* 4. Print mutation rate for each locus */
  if (opt_est_heredity && opt_print_hscalars)
  {
    for (i = 0; i < opt_locus_count; ++i)
      fprintf(fp, "\theredity_L%d", i+1);
  }

  /* 5. Print log likelihood */
  if (opt_usedata)
    fprintf(fp, "\tlnL\n"); 
  else
    fprintf(fp, "\n");
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
                           locus_t ** locus,
                           long dparam_count,
                           long ndspecies)
{
  unsigned int i;
  unsigned int snodes_total;
  
  if (opt_msci)
    snodes_total = stree->tip_count + stree->inner_count + stree->hybrid_count;
  else
    snodes_total = stree->tip_count + stree->inner_count;

  if (opt_method == METHOD_01)          /* species tree inference */
  {
    char * newick = stree_export_newick(stree->root, cb_serialize_branch);
    fprintf(fp, "%s\n", newick);
    free(newick);
    return;
  }

  if (opt_method == METHOD_11)    /* species tree inference and delimitation */
  {
    char * newick = stree_export_newick(stree->root, cb_serialize_branch);
    fprintf(fp, "%s %ld\n", newick, ndspecies);
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

  /* TODO: Combine the next two loops? */

  /* first print thetas for tips */
  if (opt_est_theta)
  {
    for (i = 0; i < stree->tip_count; ++i)
      if (stree->nodes[i]->theta >= 0)
        fprintf(fp, "\t%.5g", stree->nodes[i]->theta);
  }

  /* then for inner nodes */
  if (opt_est_theta)
  {
    /* TODO: Is the 'has_theta' check also necessary ? */
    for (i = stree->tip_count; i < snodes_total; ++i)
      if (stree->nodes[i]->theta >= 0)
        fprintf(fp, "\t%.5g", stree->nodes[i]->theta);
  }

  /* 2. Print taus for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->tau)
      fprintf(fp, "\t%.5g", stree->nodes[i]->tau);

  /* 2a. Print phi for hybridization nodes */
  if (opt_msci)
  {
    unsigned int offset=stree->tip_count+stree->inner_count;
    for (i = 0; i < stree->hybrid_count; ++i)
      fprintf(fp, "\t%.5g", stree->nodes[offset+i]->hybrid->hphi);
  }

  /* 3. Print mutation rate for each locus */
  if (opt_est_locusrate && opt_print_locusrate)
  {
    for (i = 0; i < opt_locus_count; ++i)
      fprintf(fp, "\t%.5g", locus[i]->mut_rates[0]);
  }

  /* 4. Print mutation rate for each locus */
  if (opt_est_heredity && opt_print_hscalars)
  {
    for (i = 0; i < opt_locus_count; ++i)
      fprintf(fp, "\t%.5g", locus[i]->heredity[0]);
  }

  /* 22.6.2018 - Testing gene tree node age proposal for MSCi */
  //fprintf(fp, "\t%f\t%s", gtree[0]->root->time, gtree[0]->root->pop->label);

  /* 5. print log-likelihood if usedata=1 */
  if (opt_usedata)
  {
    double logl = 0;

    for (i = 0; i < stree->locus_count; ++i)
      logl += gtree[i]->logl;

    fprintf(fp, "\t%.3f\n", logl/opt_bfbeta);
  }
  else
    fprintf(fp, "\n");
}

static void print_gtree(FILE ** fp, gtree_t ** gtree)
{
  long i;

  for (i = 0; i < opt_locus_count; ++i)
  {
    char * newick = gtree_export_newick(gtree[i]->root,NULL);
    fprintf(fp[i], "%s\n", newick);
    free(newick);
  }
}

static FILE * resume(stree_t ** ptr_stree,
                     gtree_t *** ptr_gtree,
                     locus_t *** ptr_locus,
                     double ** ptr_pjump,
                     unsigned long * ptr_curstep,
                     long * ptr_ft_round,
                     long * ptr_ndspecies,

                     long * ptr_dparam_count,
                     double ** ptr_posterior,
                     double ** ptr_pspecies,
                     long * ptr_ft_round_rj,
                     double * ptr_pjump_rj,

                     long * ptr_ft_round_spr,
                     long * ptr_pjump_slider,
                     double * ptr_mean_logl,
                     double ** ptr_mean_tau,
                     double ** ptr_mean_theta,
                     long * ptr_mean_tau_count,
                     long * ptr_mean_theta_count,
                     stree_t ** ptr_sclone, 
                     gtree_t *** ptr_gclones,
                     FILE *** ptr_fp_gtree,
                     FILE ** ptr_fp_out)
{
  long i,j;
  FILE * fp_mcmc;
  FILE * fp_out;
  long mcmc_offset;
  long out_offset;
  long * gtree_offset;
  char ** gtree_files = NULL;

  if (sizeof(BYTE) != 1)
    fatal("Checkpoint does not work on systems with sizeof(char) <> 1");

  /* load data from checkpoint file */
  checkpoint_load(ptr_gtree,
                  ptr_locus,
                  ptr_stree,
                  ptr_pjump,
                  ptr_curstep,
                  ptr_ft_round,
                  ptr_ndspecies,
                  &mcmc_offset,
                  &out_offset,
                  &gtree_offset,
                  ptr_dparam_count,
                  ptr_posterior,
                  ptr_pspecies,
                  ptr_ft_round_rj,
                  ptr_pjump_rj,
                  ptr_ft_round_spr,
                  ptr_pjump_slider,
                  ptr_mean_logl,
                  ptr_mean_tau,
                  ptr_mean_theta,
                  ptr_mean_tau_count,
                  ptr_mean_theta_count);

  /* truncate MCMC file to specific offset */
  checkpoint_truncate(opt_mcmcfile, mcmc_offset);

  /* truncate output file to specific offset */
  checkpoint_truncate(opt_outfile, out_offset);

  /* truncate gene tree files if available */
  if (opt_print_genetrees)
  {
    assert(gtree_offset);
    gtree_files = (char **)xmalloc((size_t)opt_locus_count*sizeof(char *));

    for (i = 0; i < opt_locus_count; ++i)
    {
      char * s = NULL;
      xasprintf(&s, "%s.gtree.L%ld", opt_outfile, i+1);
      gtree_files[i] = s;
      checkpoint_truncate(s,gtree_offset[i]);
    }
    free(gtree_offset);
  }

  gtree_t ** gtree = *ptr_gtree;
  stree_t  * stree = *ptr_stree;

  gtree_alloc_internals(gtree,opt_locus_count);
  reset_gene_leaves_count(stree,gtree);
  stree_reset_pptable(stree);


  /* compute MSC density */
  locus_t ** locus = *ptr_locus;
  if (opt_est_theta)
  {
    for (i = 0; i < opt_locus_count; ++i)
      gtree[i]->logpr = gtree_logprob(stree,
                                      locus[i]->heredity[0],
                                      i,
                                      thread_index_zero);
  }
  else
  {
    //assert(0);
  }

  /* set old_pop to NULL */
  for (i = 0; i < opt_locus_count; ++i)
    for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      gtree[i]->nodes[j]->old_pop = NULL;

  for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
    for (i = 0; i < opt_threads; ++i)
      stree->nodes[j]->mark[i] = 0;


  /* set method */
  if (!opt_est_stree && !opt_est_delimit)
    opt_method = METHOD_00;
  else if (!opt_est_stree)
    opt_method = METHOD_10;
  else if (!opt_est_delimit)
    opt_method = METHOD_01;
  else
    opt_method = METHOD_11;

  /* open truncated MCMC file for appending */
  if (!(fp_mcmc = fopen(opt_mcmcfile, "a")))
    fatal("Cannot open file %s for appending...", opt_mcmcfile);
  if (!(fp_out = fopen(opt_outfile, "a")))
    fatal("Cannot open file %s for appending...", opt_outfile);
  *ptr_fp_out = fp_out;

  /* open potential truncated gene trees files for appending */
  *ptr_fp_gtree = NULL;
  if (opt_print_genetrees)
  {
    FILE ** fp_gtree = (FILE **)xmalloc((size_t)opt_locus_count*sizeof(FILE *));
    for (i = 0; i < opt_locus_count; ++i)
    {
      if (!(fp_gtree[i] = fopen(gtree_files[i], "a")))
        fatal("Cannot open file %s for appending...", gtree_files[i]);
      free(gtree_files[i]);
    }
    free(gtree_files);
    *ptr_fp_gtree = fp_gtree;
  }

  /* if we are infering the species tree, then create another cloned copy of the
     species tree and gene trees */
  if (opt_est_stree)
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
     4) Find delimitation string and return index using delimit_getindex()
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
  else if (opt_method == METHOD_11)
  {
    rj_init(gtree,stree,opt_locus_count);
    if (opt_method == METHOD_11)
      partition_fast(opt_max_species_count);
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
                   double ** ptr_posterior,
                   long * ptr_ft_round_rj,
                   double * ptr_pjump_rj,
                   long * ptr_ft_round_spr,
                   long * ptr_pjump_slider,
                   double * ptr_mean_logl,
                   stree_t ** ptr_sclone, 
                   gtree_t *** ptr_gclones,
                   FILE *** ptr_fp_gtree,
                   FILE ** ptr_fp_out)
{
  long i,j;
  long msa_count;
  double logl,logpr;
  double logl_sum = 0;
  double logpr_sum = 0;
  double * pjump;
  FILE * fp_mcmc = NULL;
  FILE * fp_out;
  list_t * map_list = NULL;
  stree_t * stree;
  FILE ** fp_gtree;
  msa_t ** msa_list;
  gtree_t ** gtree;
  locus_t ** locus;

  /* method 10 specific variables */
  long dparam_count = 0;

  /* method 01 specific variables */
  stree_t * sclone = NULL;
  gtree_t ** gclones = NULL;

  /* load species tree */
  stree = load_tree_or_network();
  printf(" Done\n");

  /* Show network */
  if (opt_msci)
  {
    if (opt_finetune_phi == -1)
      fatal("Missing finetune value for phi parameter");

    print_network_table(stree);
  }

  /* parse the phylip file */
  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  assert(fd);

  printf("Parsing phylip file...");
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);
  printf(" Done\n");

  phylip_close(fd);
  if (opt_locus_count > msa_count)
    fatal("Expected %ld loci but found only %ld", opt_locus_count, msa_count);

  /* set global variable with number of loci, if not set */
  if (!opt_locus_count)
    opt_locus_count = msa_count;

  /* check for locusrate and one locus */
  if (opt_locus_count == 1 && opt_est_locusrate)
    fatal("Cannot use option 'locusrate' with only one locus");

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

  /* parse map file */
  if (stree->tip_count > 1)
  {
    printf("Parsing map file...");
    map_list = yy_parse_map(opt_mapfile);
    printf(" Done\n");
  }
  #if 0
  maplist_print(map_list);
  #endif

  if (!opt_onlysummary)
  {
    if (!(fp_mcmc = fopen(opt_mcmcfile, "w")))
      fatal("Cannot open file %s for writing...");
  }

  if (!(fp_out = fopen(opt_outfile, "w")))
    fatal("Cannot open file %s for writing...");
  *ptr_fp_out = fp_out;

  /* print compressed alignmens in output file */
  fprintf(fp_out, "COMPRESSED ALIGNMENTS\n\n");

  /* print the alignments */
  msa_print_phylip(fp_out,msa_list,msa_count);

  *ptr_fp_gtree = NULL;

  /* if print gtree */
  if (opt_print_genetrees)
  {
    fp_gtree = (FILE **)xmalloc((size_t)opt_locus_count*sizeof(FILE *));
    for (i = 0; i < opt_locus_count; ++i)
    {
      char * s = NULL;
      xasprintf(&s, "%s.gtree.L%ld", opt_outfile, i+1);
      fp_gtree[i] = xopen(s,"w");
      free(s);
    }
  }
  else
    fp_gtree = NULL;

  *ptr_fp_gtree = fp_gtree;

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
    fprintf(fp_out, "COMPRESSED ALIGNMENTS AFTER PHASING OF DIPLOID SEQUENCES\n\n");
    msa_print_phylip(fp_out,msa_list,msa_count);

  }

  /* allocate TLS mark variables on stree as they are used in delimitations_init
     TODO: Move this allocation somewhere else */
  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
    stree->nodes[i]->mark = (int *)xcalloc((size_t)opt_threads,sizeof(int));

  if (opt_method == METHOD_10)          /* species delimitation */
  {
    assert(opt_msci == 0);

    long dmodels_count = delimitations_init(stree);
    printf("Number of delimitation models: %ld\n", dmodels_count);

    *ptr_posterior = (double *)xcalloc((size_t)dmodels_count,sizeof(double));
  }

  /* initialize species tree (tau + theta) */
  stree_init(stree,msa_list,map_list,msa_count,fp_out);

  stree_show_pptable(stree, BPP_FALSE);

  if (opt_method == METHOD_11)
  {
    assert(opt_msci == 0);
    partition_fast(stree->tip_count);
  }

  /* allocate arrays for locus mutation rate and heredity scalars */
  double * locusrate = (double*)xmalloc((size_t)opt_locus_count*sizeof(double));
  double * heredity = (double *)xmalloc((size_t)opt_locus_count*sizeof(double));

  /* initialize both to 1 in case no options regarding them were given */
  for (i = 0; i < opt_locus_count; ++i)
    locusrate[i] = heredity[i] = 1;

  /* initialize heredity scalars if estimation was selected */
  if (opt_est_heredity == HEREDITY_ESTIMATE)
  {
    for (i = 0; i < opt_locus_count; ++i)
      heredity[i] = opt_heredity_alpha /
                    opt_heredity_beta*(0.8 + 0.4*legacy_rndu(0));

    /* TODO: Perhaps we can avoid the check every 100-th term by using the log
       of heredity scaler from the beginning. E.g. if this loop is replaced by
       
       for (j = 0; j < opt_locus_count; ++j)
         logpr -= log(locus[j]->heredity[0]);

       then we only need to add and subtract the two corresponding heredity
       multipliers (the old and new)
    */
    if (!opt_est_theta)
    {
      double hfactor = 0;
      double y = 1;
      for (i = 0; i < opt_locus_count; ++i)
      {
        y *= heredity[i];
        if ((i+1) % 100 == 0)
        {
          hfactor -= log(y);
          y = 1;
        }
      }
      hfactor -= log(y);
      stree->notheta_hfactor = hfactor;
    }
  }
  else if (opt_est_heredity == HEREDITY_FROMFILE)
  {
    long errcontext = 0;
    int rc = parsefile_doubles(opt_heredity_filename,
                               opt_locus_count,
                               heredity,
                               &errcontext);
    if (rc == ERROR_PARSE_MORETHANEXPECTED)
      fatal("File %s contains more heredity scalers than number of loci (%ld)",
            opt_heredity_filename, opt_locus_count);
    else if (rc == ERROR_PARSE_LESSTHANEXPECTED)
      fatal("File %s contains less heredity scalers (%ld) than number of loci (%ld)",
            opt_heredity_filename, errcontext, opt_locus_count);
    else if (rc == ERROR_PARSE_INCORRECTFORMAT)
      fatal("Incorrect format of file %s at line %ld",
            opt_heredity_filename, errcontext);

    /* disable estimation of heredity scalars */
    opt_est_heredity = 0;
  }

  /* initialize locus mutation rates if estimation was selected */
  if (opt_est_locusrate == MUTRATE_ESTIMATE)
  {
    double mean = 0;
    for (i = 0; i < opt_locus_count; ++i)
    {
      locusrate[i] = 0.8 + 0.4*legacy_rndu(thread_index_zero);
      mean += locusrate[i];
    }

    mean /= opt_locus_count;

    for (i = 0; i < opt_locus_count; ++i)
      locusrate[i] /= mean;
  }
  else if (opt_est_locusrate == MUTRATE_FROMFILE)
  {
    long errcontext = 0;
    int rc = parsefile_doubles(opt_locusrate_filename,
                               opt_locus_count,
                               locusrate,
                               &errcontext);
    if (rc == ERROR_PARSE_MORETHANEXPECTED)
      fatal("File %s contains more rates than number of loci (%ld)",
            opt_locusrate_filename, opt_locus_count);
    else if (rc == ERROR_PARSE_LESSTHANEXPECTED)
      fatal("File %s contains less rates (%ld) than number of loci (%ld)",
            opt_locusrate_filename, errcontext, opt_locus_count);
    else if (rc == ERROR_PARSE_INCORRECTFORMAT)
      fatal("Incorrect format of file %s at line %ld",
            opt_locusrate_filename, errcontext);

    double mean = 0;
    for (i = 0; i < opt_locus_count; ++i)
      mean += locusrate[i];

    mean /= opt_locus_count;
      
    for (i = 0; i < opt_locus_count; ++i)
      locusrate[i] /= mean;

    /* disable estimation of mutation rates */
    opt_est_locusrate = 0;
  }

  /* We must first link tip sequences (gene tips) to populations */
  if (opt_est_delimit)          /* species delimitation */
  {
    assert(opt_msci == 0);
    stree_rootdist(stree,map_list,msa_list,weights);
  }

  gtree = gtree_init(stree,msa_list,map_list,msa_count);

  /* the below two lines are specific to method 01 and they generate
     space for cloning the species and gene trees */
  if (opt_est_stree)            /* species tree inference */
  {
    assert(opt_msci == 0);
    sclone = stree_clone_init(stree);
    gclones = (gtree_t **)xmalloc((size_t)msa_count*sizeof(gtree_t *));
    for (i = 0; i < msa_count; ++i)
      gclones[i] = gtree_clone_init(gtree[i], sclone);
  }

  locus = (locus_t **)xcalloc((size_t)msa_count, sizeof(locus_t *));

  /* Check that only first 32 bits of opt_arch are used */

#ifdef _WIN32
  assert(opt_arch < ((long long)1 << 32) - 1);
#else
  assert(opt_arch < (1L << 32) - 1);
#endif

  /* ensure that heredity / locusrate estimation is now set to either 0 or 1 */
  assert(opt_est_locusrate >= 0 && opt_est_locusrate <= 1);
  assert(opt_est_heredity  >= 0 && opt_est_heredity  <= 1);

  gtree_update_branch_lengths(gtree, msa_count);

  for (i = 0; i < msa_count; ++i)
  {
    msa_t * msa = msa_list[i];
    double frequencies[4] = {0.25, 0.25, 0.25, 0.25};
    unsigned int pmatrix_count = gtree[i]->edge_count;

    unsigned int scale_buffers = opt_scaling ?
                                   2*gtree[i]->inner_count : 0;

    /* activate twice as many transition probability matrices (for reverting in
       locusrate, species tree SPR and mixing proposals)  */
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
                            scale_buffers,              /* # scale buffers */
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

    /* TODO with more complex mixture models where rate_matrices > 1 we need
       to revisit this */
    assert(rate_matrices == 1);
    locus_set_mut_rates(locus[i],locusrate+i);
    locus_set_heredity_scalers(locus[i],heredity+i);

    /* set pattern weights and free the weights array */
    if (locus[i]->diploid)
    {
      /* TODO: 1) pattern_weights_sum is not updated here, but it is not used in
         the program, perhaps remove.
         2) pattern_weights is allocated in locus_create with a size msa->length
            equal to length of A3, but in reality we only need |A1| storage
            space. Free and reallocate here. *UPDATE* Actually |A1| may be larger
            than |A3| !! */

      free(locus[i]->pattern_weights);
      locus[i]->pattern_weights = (unsigned int *)xmalloc((size_t)
                                     (unphased_length[i])*sizeof(unsigned int));

      locus[i]->diploid_mapping = mapping[i];
      locus[i]->diploid_resolution_count = resolution_count[i];
      /* since PLL does not support diploid sequences we make a small hack */
      memcpy(locus[i]->pattern_weights,
             weights[i],
             unphased_length[i]*sizeof(unsigned int));
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
    if (opt_est_theta)
    {
      logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index_zero);
      gtree[i]->logpr = logpr;
      logpr_sum += logpr;
    }
    else
    {
      for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
        logpr_sum += gtree_update_logprob_contrib(stree->nodes[j],
                                                  locus[i]->heredity[0],
                                                  i,
                                                  thread_index_zero);
    }
  }
  if (!opt_est_theta)
  {
    logpr_sum = 0;
    for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      logpr_sum += stree->nodes[j]->notheta_logpr_contrib;
    stree->notheta_logpr += logpr_sum;
    stree->notheta_old_logpr = 0;

    logpr_sum = stree->notheta_logpr;
  }

  /* deallocate unnecessary arrays */
  free(locusrate);
  free(heredity);
  if (opt_diploid)
  {
    free(mapping);
    free(unphased_length);
    free(resolution_count);
  }

  #if 0
  debug_print_network_node_attribs(stree);
  #endif

  printf("\nInitial MSC density and log-likelihood of observing data:\n");
  printf("log-P0 = %f   log-L0 = %f\n\n", logpr_sum, logl_sum);

  /* free weights array */
  free(weights);

  if (opt_est_delimit)          /* species delimitation */
    rj_init(gtree,stree,msa_count);

  /* initialize pjump and finetune rounds */
  if (opt_est_locusrate || opt_est_heredity)
    pjump = (double *)xcalloc(PROP_COUNT+1, sizeof(double));
  else
    pjump = (double *)xcalloc(PROP_COUNT, sizeof(double));

  /* TODO: Method 10 has a commented call to 'delimit_resetpriors()' */
  //delimit_resetpriors();

  /* if method 00 or 01 print corresponding header line in MCMC file */
  if (!opt_onlysummary)
  {
    if (opt_method == METHOD_01)
      mcmc_printinitial(fp_mcmc,stree);
    else
    {
      if (opt_method != METHOD_11)
        mcmc_printheader(fp_mcmc,stree);
    }
  }

  #if 0
  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);
  #endif

  /* TODO: This is only for species delimitation */
  if (opt_est_delimit)          /* species delimitation */
  {
    dparam_count = 0;
    for (i = 0; i < (int)(stree->tip_count + stree->inner_count); ++i)
    {
      if (opt_est_theta && stree->nodes[i]->theta > 0)
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

  *ptr_sclone = sclone;
  *ptr_gclones = gclones;

  /* deallocate maplist */
  if (stree->tip_count > 1)
  {
    list_clear(map_list,map_dealloc);
    free(map_list);
  }

  /* deallocate alignments */
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  return fp_mcmc;

}

void cmd_run()
{
  /* common variables for all methods */
  long i,j,k;
  long ft_round;
  double logl_sum = 0;
  double * pjump;
  FILE * fp_mcmc;
  FILE * fp_out;
  stree_t * stree;
  FILE ** fp_gtree = NULL;
  gtree_t ** gtree;
  locus_t ** locus;
  long * gtree_offset = NULL;   /* for checkpointing when printing gene trees */
  double ratio;
  long ndspecies;


  /* method 10 specific variables */
  long dparam_count = 0;
  long ft_round_rj;
  double pjump_rj = 0;
  double * posterior = NULL;
  double * pspecies = NULL;

  /* method 01 specific variables */
  long ft_round_spr = 0;
  long pjump_slider;
  long printk;// = opt_samplefreq * opt_samples;
  double mean_logl = 0;

  double pjump_phi = 0;

  double * mean_tau = NULL;
  double * mean_theta = NULL;

  long mean_theta_count;
  long mean_tau_count;

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
                     &ndspecies,
                     &dparam_count,
                     &posterior,
                     &pspecies,
                     &ft_round_rj,
                     &pjump_rj,
                     &ft_round_spr,
                     &pjump_slider,
                     &mean_logl,
                     &mean_tau,
                     &mean_theta,
                     &mean_tau_count,
                     &mean_theta_count,
                     &sclone, 
                     &gclones,
                     &fp_gtree,
                     &fp_out);
  else
  {
    fp_mcmc = init(&stree,
                   &gtree,
                   &locus,
                   &pjump,
                   &curstep,
                   &ft_round,
                   &dparam_count,
                   &posterior,
                   &ft_round_rj,
                   &pjump_rj,
                   &ft_round_spr,
                   &pjump_slider,
                   &mean_logl,
                   &sclone, 
                   &gclones,
                   &fp_gtree,
                   &fp_out);

    /* allocate mean_tau and mean_theta */
    if (opt_est_theta)
      mean_theta = (double *)xcalloc(MAX_THETA_OUTPUT,sizeof(double));
    mean_tau   = (double *)xcalloc(MAX_TAU_OUTPUT,sizeof(double));

    /* count number of delimited species with current species tree */
    ndspecies = 1;
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
      if (stree->nodes[i]->tau > 0)
        ++ndspecies;
    
    opt_max_species_count = stree->tip_count;

    /* TODO: pspecies should be in the checkpoint */
    if (opt_method == METHOD_11)
      pspecies = (double *)xcalloc((size_t)(stree->tip_count),sizeof(double));

  }

  if (opt_checkpoint && opt_print_genetrees)
    gtree_offset = (long *)xmalloc((size_t)opt_locus_count*sizeof(long));
  if (opt_exp_randomize)
    fprintf(stdout, "[EXPERIMENTAL] - Randomize nodes order on gtree SPR\n");
  if (opt_rev_gspr)
    fprintf(stdout, "[EXPERIMENTAL] - Revolutionary gene tree SPR algorithm\n");
  if (opt_revolutionary_spr_method)
    fprintf(stdout, "[EXPERIMENTAL] - Revolutionary species tree SPR algorithm\n");

  if (opt_threads > 1)
  {
    threads_init();
    memset(&td,0,sizeof(td));
  }

  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);

  printk = opt_samplefreq * opt_samples;

  /* check if summary only was requested (no MCMC) and initialize counter
     for MCMC loop appropriately */
  if (opt_onlysummary)
    i = opt_samples*opt_samplefreq;
  else
    i = curstep - opt_burnin;


  /* 9.10.2018 - Testing gene tree node age proposal for MSCi **************** */
#if (defined DEBUG_MSCi)
  /* two bidirections */
  int dbg_locus = 0;
  double dbg_left[4] = {0};
  double dbg_right[4] = {0};
  double dbg_tmean = 0;
  double dbg_ppop[11] = { 0 };
  printf("rndu status: %d\n", get_legacy_rndu_status());
#endif
#if (0)
  /* one bidirection */
  int dbg_locus = 0;
  double dbg_left[2] = {0};
  double dbg_right[2] = {0};
  double dbg_tmean = 0;
  double dbg_ppop[11] = { 0 };
  printf("rndu status: %d\n", get_legacy_rndu_status());
#endif
#if (0)
  int dbg_locus = 1;
  double dbg_left[2] = {0};
  double dbg_right[2] = {0};
  double dbg_tmean = 0;
  double dbg_ppop[11] = { 0 };
  printf("rndu status: %d\n", get_legacy_rndu_status());
#endif

  /* *** start of MCMC loop *** */
  for (; i < opt_samples*opt_samplefreq; ++i)
  {
    /* update progress bar */
    if (!opt_quiet)
      progress_update(curstep);

    /* reset finetune parameters */
    if (i == 0 || (opt_finetune_reset && opt_burnin >= 200 && i < 0 &&
                   ft_round >= 100 && i%(opt_burnin/4)==0))
    {
      int pjump_size = PROP_COUNT + (opt_est_locusrate || opt_est_heredity);

/* 9.10.2018 - Testing gene tree node age proposal for MSCi **************** */
#if (defined DEBUG_MSCi)
      /* two bidirections */
      dbg_tmean = 0;
      memset(dbg_left, 0, 4 * sizeof(double));
      memset(dbg_right, 0, 4 * sizeof(double));
      memset(dbg_ppop, 0, 11 * sizeof(double));
#endif      
#if (0)
      /* one bidirection */
      dbg_tmean = 0;
      memset(dbg_left, 0, 2 * sizeof(double));
      memset(dbg_right, 0, 2 * sizeof(double));
      memset(dbg_ppop, 0, 11 * sizeof(double));
#endif      

      if (opt_finetune_reset && opt_burnin >= 200)
        reset_finetune(pjump,&pjump_phi);
      for (j = 0; j < pjump_size; ++j)
        pjump[j] = 0;

      /* reset pjump and number of steps since last finetune reset to zero */
      ft_round = 0;
      memset(pjump, 0, pjump_size * sizeof(double));

      if (opt_msci)
        pjump_phi = 0;

      if (opt_est_delimit)
      {
        pjump_rj = 0;
        ft_round_rj = 0;
      }
      if (opt_method == METHOD_10)      /* species delimitation */
        memset(posterior,0,delimitation_getparam_count()*sizeof(double));

      if (opt_est_stree)
      {
        ft_round_spr = 0;
        pjump_slider = 0;
      }
      if (opt_method == METHOD_11)      /* species tree inference + delimitation */
        memset(pspecies,0,stree->tip_count*sizeof(double));
    }
    
    ++ft_round;

    /* propose delimitation through merging/splitting of nodes */
    if (opt_est_delimit)        /* species delimitation */
    {
      if (legacy_rndu(thread_index_zero) < 0.5)
        j = prop_split(gtree,stree,locus,0.5,&dparam_count,&ndspecies);
      else
        j = prop_join(gtree,stree,locus,0.5,&dparam_count,&ndspecies);

      if (j != 2)
      {
        ft_round_rj++;
        pjump_rj += j;
      }
    }

    /* propose species tree topology using SPR */
    if (ndspecies > 2 && (opt_est_stree))
    {
      if (legacy_rndu(thread_index_zero) > 0)   /* bpp4 compatible results (RNG to next state) */
      {
        long ret;
        ret = stree_propose_spr(&stree, &gtree, &sclone, &gclones, locus);
        if (ret == 1)
        {
          /* accepted */
          /* swap the pointers of species tree and gene tree list with cloned */
          SWAP(stree,sclone);
          SWAP(gtree,gclones);
          stree_label(stree);
          pjump_slider++;
        }
        if (ret != 2)
          ft_round_spr++;
      }
    }

    /* perform proposals sequentially */   

    /* propose gene tree ages */
    if (opt_threads == 1)
      ratio = gtree_propose_ages_serial(locus, gtree, stree);
    else
    {
      td.locus = locus; td.gtree = gtree; td.stree = stree;
      threads_wakeup(THREAD_WORK_GTAGE,&td);
      ratio = td.accepted ? ((double)(td.accepted)/td.proposals) : 0;

    }
    pjump[0] = (pjump[0]*(ft_round-1) + ratio) / (double)ft_round;

    /* propose gene tree topologies using SPR */
    if (opt_threads == 1)
      ratio = gtree_propose_spr_serial(locus,gtree,stree);
    else
    {
      td.locus = locus; td.gtree = gtree; td.stree = stree;
      threads_wakeup(THREAD_WORK_GTSPR,&td);
      ratio = td.accepted ? ((double)(td.accepted)/td.proposals) : 0;
    }
    pjump[1] = (pjump[1]*(ft_round-1) + ratio) / (double)ft_round;

/* 9.10.2018 - Testing gene tree node age proposal for MSCi **************** */
#if (defined DEBUG_MSCi)
    /* two bidirections */
    int aindex = 0, cindex = gtree[dbg_locus]->tip_count - 1;
    int dbg_lca = 0, dbg_lcapop = 0;

    /* safety check */
    for (j = 0; j < gtree[dbg_locus]->tip_count + gtree[dbg_locus]->inner_count; ++j)
      assert(gtree[dbg_locus]->nodes[j]->mark == 0);

    /* mark all nodes on the path from a to root */
    gnode_t * dbg_tmp = gtree[dbg_locus]->nodes[aindex];
    while (dbg_tmp)
    {
       dbg_tmp->mark = 1;
       dbg_tmp = dbg_tmp->parent;
    }
    /* find first marked node on the path from c to root */
    dbg_tmp = gtree[dbg_locus]->nodes[cindex];
    while (!dbg_tmp->mark)
       dbg_tmp = dbg_tmp->parent;

    /* clear marks */
    for (j = 0; j < gtree[dbg_locus]->inner_count; ++j)
       gtree[dbg_locus]->nodes[gtree[dbg_locus]->tip_count + j]->mark = 0;
    gtree[dbg_locus]->nodes[aindex]->mark = gtree[dbg_locus]->nodes[cindex]->mark = 0;

    dbg_lca = dbg_tmp->node_index;
    assert(dbg_lca == 2);
    dbg_lcapop = gtree[dbg_locus]->nodes[dbg_lca]->pop->node_index;
    dbg_tmean = (dbg_tmean*(ft_round - 1) + gtree[dbg_locus]->nodes[dbg_lca]->time) / ft_round;
    dbg_ppop[dbg_lcapop]++;

    /* count X */
    if (gtree[dbg_locus]->nodes[aindex]->hpath[1] == BPP_HPATH_LEFT)
       dbg_left[1]++;
    else if (gtree[dbg_locus]->nodes[aindex]->hpath[1] == BPP_HPATH_RIGHT)
       dbg_right[1]++;
    /* count Y */
    if (gtree[dbg_locus]->nodes[cindex]->hpath[2] == BPP_HPATH_LEFT)
       dbg_left[2]++;
    else if (gtree[dbg_locus]->nodes[cindex]->hpath[2] == BPP_HPATH_RIGHT)
       dbg_right[2]++;

    /* count S */
    if (gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_LEFT || gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_LEFT)
       dbg_left[0]++;
    if (gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_RIGHT || gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_RIGHT)
       dbg_right[0]++;
    if (gtree[dbg_locus]->nodes[2]->hpath[0] == BPP_HPATH_LEFT)
    {
      assert(gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_NONE);
      assert(gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_NONE);
      dbg_left[0]++;
    }
    else if (gtree[dbg_locus]->nodes[2]->hpath[0] == BPP_HPATH_RIGHT)
    {
      assert(gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_NONE);
      assert(gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_NONE);
      dbg_right[0]++;
    }
    #if 0
    if (gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_LEFT)
       dbg_left[0]++;
    else if (gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_RIGHT)
       dbg_right[0]++;
    #endif

    /* count T */
    if (gtree[dbg_locus]->nodes[aindex]->hpath[3] == BPP_HPATH_LEFT || gtree[dbg_locus]->nodes[cindex]->hpath[3] == BPP_HPATH_LEFT)
       dbg_left[3]++;
    if (gtree[dbg_locus]->nodes[aindex]->hpath[3] == BPP_HPATH_RIGHT || gtree[dbg_locus]->nodes[cindex]->hpath[3] == BPP_HPATH_RIGHT)
       dbg_right[3]++;
    if (gtree[dbg_locus]->nodes[2]->hpath[3] == BPP_HPATH_LEFT)
    {
      assert(gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_NONE);
      assert(gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_NONE);
      dbg_left[3]++;
    }
    else if (gtree[dbg_locus]->nodes[2]->hpath[3] == BPP_HPATH_RIGHT)
    {
      assert(gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_NONE);
      assert(gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_NONE);
      dbg_right[3]++;
    }
    #if 0
    if (gtree[dbg_locus]->nodes[cindex]->hpath[3] == BPP_HPATH_LEFT)
       dbg_left[3]++;
    else if (gtree[dbg_locus]->nodes[cindex]->hpath[3] == BPP_HPATH_RIGHT)
       dbg_right[3]++;
    #endif
#endif
#if (0)
    /* one bidirection */
    int aindex = 0, cindex = gtree[dbg_locus]->tip_count - 1;
    int dbg_lca = 0, dbg_lcapop = 0;

    /* safety check */
    for (j = 0; j < gtree[dbg_locus]->tip_count + gtree[dbg_locus]->inner_count; ++j)
      assert(gtree[dbg_locus]->nodes[j]->mark == 0);

    /* mark all nodes on the path from a to root */
    gnode_t * dbg_tmp = gtree[dbg_locus]->nodes[aindex];
    while (dbg_tmp)
    {
       dbg_tmp->mark = 1;
       dbg_tmp = dbg_tmp->parent;
    }
    /* find first marked node on the path from c to root */
    dbg_tmp = gtree[dbg_locus]->nodes[cindex];
    while (!dbg_tmp->mark)
       dbg_tmp = dbg_tmp->parent;

    /* clear marks */
    for (j = 0; j < gtree[dbg_locus]->inner_count; ++j)
       gtree[dbg_locus]->nodes[gtree[dbg_locus]->tip_count + j]->mark = 0;
    gtree[dbg_locus]->nodes[aindex]->mark = gtree[dbg_locus]->nodes[cindex]->mark = 0;

    dbg_lca = dbg_tmp->node_index;
    assert(dbg_lca == 2);
    dbg_lcapop = gtree[dbg_locus]->nodes[dbg_lca]->pop->node_index;
    dbg_tmean = (dbg_tmean*(ft_round - 1) + gtree[dbg_locus]->nodes[dbg_lca]->time) / ft_round;
    dbg_ppop[dbg_lcapop]++;

    if (gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_LEFT)
       dbg_left[0]++;
    else if (gtree[dbg_locus]->nodes[aindex]->hpath[0] == BPP_HPATH_RIGHT)
       dbg_right[0]++;
    if (gtree[dbg_locus]->nodes[cindex]->hpath[1] == BPP_HPATH_LEFT)
       dbg_left[1]++;
    else if (gtree[dbg_locus]->nodes[cindex]->hpath[1] == BPP_HPATH_RIGHT)
       dbg_right[1]++;
#endif
#if (0)
    int aindex = 0, cindex = gtree[dbg_locus]->tip_count - 1;
    int dbg_lca = 0, dbg_lcapop = 0;

    /* safety check */
    for (j = 0; j < gtree[dbg_locus]->tip_count + gtree[dbg_locus]->inner_count; ++j)
      assert(gtree[dbg_locus]->nodes[j]->mark == 0);

    /* mark all nodes on the path from a to root */
    gnode_t * dbg_tmp = gtree[dbg_locus]->nodes[aindex];
    while (dbg_tmp)
    {
       dbg_tmp->mark = 1;
       dbg_tmp = dbg_tmp->parent;
    }
    /* find first marked node on the path from c to root */
    dbg_tmp = gtree[dbg_locus]->nodes[cindex];
    while (!dbg_tmp->mark)
       dbg_tmp = dbg_tmp->parent;

    /* clear marks */
    for (j = 0; j < gtree[dbg_locus]->inner_count; ++j)
       gtree[dbg_locus]->nodes[gtree[dbg_locus]->tip_count + j]->mark = 0;
    gtree[dbg_locus]->nodes[aindex]->mark = gtree[dbg_locus]->nodes[cindex]->mark = 0;

    dbg_lca = dbg_tmp->node_index;
    dbg_lcapop = gtree[dbg_locus]->nodes[dbg_lca]->pop->node_index;
    dbg_tmean = (dbg_tmean*(ft_round - 1) + gtree[dbg_locus]->nodes[dbg_lca]->time) / ft_round;
    dbg_ppop[dbg_lcapop]++;

    if (gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_LEFT)
       dbg_left[0]++;
    else if (gtree[dbg_locus]->nodes[cindex]->hpath[0] == BPP_HPATH_RIGHT)
       dbg_right[0]++;
    if (gtree[dbg_locus]->nodes[cindex]->hpath[1] == BPP_HPATH_LEFT)
       dbg_left[1]++;
    else if (gtree[dbg_locus]->nodes[cindex]->hpath[1] == BPP_HPATH_RIGHT)
       dbg_right[1]++;
#endif

    /* propose population sizes on species tree */
    if (opt_est_theta)
    {
      ratio = stree_propose_theta(gtree,locus,stree);
      pjump[2] = (pjump[2]*(ft_round-1) + ratio) / (double)ft_round;
    }

    /* propose species tree taus */
    if (stree->tip_count > 1 && stree->root->tau > 0)
    {
      ratio = stree_propose_tau(gtree,stree,locus);
      pjump[3] = (pjump[3]*(ft_round-1) + ratio) / (double)ft_round;
    }

    /* mixing step */
    ratio = proposal_mixing(gtree,stree,locus);
    pjump[4] = (pjump[4]*(ft_round-1) + ratio) / (double)ft_round;

    if (opt_est_locusrate || opt_est_heredity)
    {
      ratio = prop_locusrate_and_heredity(gtree,stree,locus,thread_index_zero);
      pjump[5] = (pjump[5]*(ft_round-1) + ratio) / (double)ft_round;
    }

    /* phi proposal */
    if (opt_msci)
    {
      ratio = stree_propose_phi(stree,gtree);
      pjump_phi = (pjump_phi*(ft_round-1) + ratio) / (double)ft_round;
    }

    /* log sample into file (dparam_count is only used in method 10) */
    if (i >= 0 && (i+1)%opt_samplefreq == 0)
    {
      mcmc_logsample(fp_mcmc,i+1,stree,gtree,locus,dparam_count,ndspecies);
      if ((i + 1) % (opt_samplefreq*10) == 0)
         fflush(fp_mcmc);

      if (opt_print_genetrees)
        print_gtree(fp_gtree,gtree);
    }

    if (opt_method == METHOD_10)
      ++posterior[delimitation_getcurindex()];

    /* posterior on number of species for displaying on screen */
    if (opt_method == METHOD_11)
    {
      pspecies[ndspecies-1] = (pspecies[ndspecies-1]*(ft_round-1)+1)/ft_round;
      for (j = 0; j < stree->tip_count; ++j)
        if (j != ndspecies-1)
          pspecies[j] = (pspecies[j] * (ft_round-1) + 0) / ft_round;
    }

      /* update stats for printing on screen */
      
      /* if species delimitation or species inference, we only compute mean
         theta and tau for the root node, otherwise, if method 00 then we
         compute at most MAX_THETA_OUTPUT mean thetas and at most MAX_TAU_OUTPUT
         taus */
    if (opt_method != METHOD_00)    /* species tree inference or delimitation */
    {
      if (opt_est_theta)
        mean_theta[0] = (mean_theta[0]*(ft_round-1)+stree->root->theta)/ft_round;

      mean_tau[0] = (mean_tau[0]*(ft_round-1)+stree->root->tau)/ft_round;

      mean_theta_count = 1;
      mean_tau_count = 1;
    }
    else
    {
      /* compute mean thetas */

      /* 1. calculate number of thetas to print */
      long max_param_count = MIN(stree->tip_count+stree->inner_count,
                                 MAX_THETA_OUTPUT);

      /* 2. calculate means */
      k = 0;
      if (opt_est_theta)
      {
        for (j=0; j < stree->tip_count+stree->inner_count; ++j)
        {
          if (stree->nodes[j]->theta < 0) continue;

          mean_theta[k] = (mean_theta[k]*(ft_round-1)+stree->nodes[j]->theta) /
                          ft_round;
          if (++k == max_param_count) break;
        }
      }
      mean_theta_count = k;

      /* compute mean taus */

      /* 1. calculate number of taus to print */
      max_param_count = MIN(stree->tip_count+stree->inner_count,
                            MAX_TAU_OUTPUT);

      /* 2. calculate means */
      k = 0;
      for (j = stree->tip_count; j < stree->tip_count+stree->inner_count; ++j)
      {
        if (stree->nodes[j]->tau == 0) continue;
        
        mean_tau[k] = (mean_tau[k]*(ft_round-1) + stree->nodes[j]->tau)/ft_round;

        if (++k == max_param_count) break;
      }
      mean_tau_count = k;

    }

    /* compute mean log-L */
    if (opt_usedata)
    {
      for (logl_sum = 0, j = 0; j < opt_locus_count; ++j)
        logl_sum += gtree[j]->logl;
      mean_logl = (mean_logl * (ft_round-1) + logl_sum / opt_bfbeta)/ft_round;
    }

    /* print MCMC status on screen */
    if (printk <= 500 || (i+1) % (printk / 200) == 0)
    {
      printf("\r%3.0f%%", (i + 1.499) / printk * 100.);
      for (j = 0; j < 5 + (opt_est_locusrate || opt_est_heredity); ++j)
        printf(" %4.2f", pjump[j]);
      if (opt_msci)
        printf(" %4.2f", pjump_phi);
      printf(" ");

      if (opt_method == METHOD_01)
        printf(" %5.4f ", ft_round_spr ? (double)pjump_slider/ft_round_spr : 0.);

      /* species delimitation specific output */
      if (opt_method == METHOD_10)
      {
        long bmodel = 0;
        for (j = 1; j < delimitation_getparam_count(); ++j)
          if (posterior[j] > posterior[bmodel]) bmodel = j;

        
        /* Number of parameters, pjump/finetune round ratio, current
           delimitation binary string, and support value of 'most-supported
           delimitation delimitation' */
         printf(" %2ld %6.4f %s", dparam_count,
                                  (ft_round_rj ? pjump_rj / ft_round_rj : 0),
                                  delimitation_getparam_string());
         printf(" P[%ld]=%6.4f",
                //posterior[bmodel],
                bmodel+1,
                //delimitation_getcurindex()+1,
                posterior[bmodel] / ft_round);

      }

      if (opt_method == METHOD_11)
      {
        long ndspeciesbest = 0;
        for (j = 1; j < stree->tip_count; ++j)
          if (pspecies[j] > pspecies[ndspeciesbest]) ndspeciesbest = j;
       
        printf(" %2ld %2ld %6.4f %6.4f P(%ld)=%6.4f",
               ndspecies,
               dparam_count,
               ft_round_rj  ? pjump_rj  / ft_round_rj : 0.,
               ft_round_spr ? pjump_slider / (double)ft_round_spr : 0.,
               ndspeciesbest + 1,
               pspecies[ndspeciesbest]);
      }

      if (opt_est_theta)
      {
        for (j = 0; j < mean_theta_count; ++j)
          printf(" %6.4f", mean_theta[j]);
        printf(" ");
      }

      for (j = 0; j < mean_tau_count; ++j)
        printf(" %6.4f", mean_tau[j]);
      printf(" ");
      
      double logpr_sum = 0;
      if (opt_est_theta)
        for (j = 0; j < opt_locus_count; ++j)
          logpr_sum += gtree[j]->logpr;
      else
        logpr_sum = stree->notheta_logpr;
      printf(" %8.5f", logpr_sum);

      if (opt_usedata)
        printf(" %8.5f", mean_logl);

/* 9.10.2018 - Testing gene tree node age proposal for MSCi **************** */
#if (defined DEBUG_MSCi)
      /* two bidirections */
      printf(" %8.5f ", dbg_tmean);
      for (j = 0; j < 11; ++j)
         printf(" %6.4f (%s) ", dbg_ppop[j] / ft_round, stree->nodes[j]->label);
      /* 9.10.2018 - *********************/
      printf("gam(X,Y) : %.6f %.6f ", dbg_left[1] / (dbg_left[1] + dbg_right[1]), 
                                    dbg_left[2] / (dbg_left[2] + dbg_right[2]));
      printf("gam(S,T) : %.6f %.6f ", dbg_left[0] / (dbg_left[0] + dbg_right[0]), 
                                    dbg_left[3] / (dbg_left[3] + dbg_right[3]));
#endif
#if (0)
      /* one bidirection */
      printf(" %8.5f ", dbg_tmean);
      for (j = 0; j < 7; ++j)
         printf(" %6.4f (%s) ", dbg_ppop[j] / ft_round, stree->nodes[j]->label);
      /* 9.10.2018 - *********************/
      printf("c1 left: %.6f %.6f ", dbg_left[0] / (dbg_left[0] + dbg_right[0]), 
                                    dbg_left[1] / (dbg_left[1] + dbg_right[1]));
#endif
#if (0)
      printf(" %8.5f ", dbg_tmean);
      for (j = 2; j < 11; ++j)
         printf(" %6.4f (%s) ", dbg_ppop[j] / ft_round, stree->nodes[j]->label);
      /* 9.10.2018 - *********************/
      printf("c1 left: %.6f %.6f ", dbg_left[0] / (dbg_left[0] + dbg_right[0]), 
                                    dbg_left[1] / (dbg_left[1] + dbg_right[1]));
#endif
      if (printk >= 50 && (i+1) % (printk / 20) == 0)
        printf("\n");
    }

    curstep++;

    /* Create a checkpoint file... */
    if (opt_checkpoint)
    {
      if (((long)curstep == opt_checkpoint_initial) ||
          (opt_checkpoint_step && ((long)curstep > opt_checkpoint_initial) &&
           (((long)curstep-opt_checkpoint_initial) % opt_checkpoint_step == 0)))
      {

        /* if gene tree printing is enabled get current file offsets */
        if (opt_print_genetrees)
          for (j = 0; j < opt_locus_count; ++j)
            gtree_offset[j] = ftell(fp_gtree[j]);

        checkpoint_dump(stree,
                        gtree,
                        locus,
                        pjump,
                        curstep,
                        ft_round,
                        ndspecies,
                        ftell(fp_mcmc),
                        ftell(fp_out),
                        gtree_offset,
                        dparam_count,
                        posterior,
                        pspecies,
                        opt_est_delimit ?
                          delimitation_getparam_count() : 0,
                        ft_round_rj,
                        pjump_rj,
                        ft_round_spr,
                        pjump_slider,
                        mean_logl,
                        mean_tau,
                        mean_theta,
                        mean_tau_count,
                        mean_theta_count);
      }
    }

  }

  progress_done();

  free(pjump);

  if (opt_threads > 1)
    threads_exit();

  if (opt_bfbeta != 1 && !opt_onlysummary)
  {
    fprintf(stdout, "\nBFbeta = %8.6f  E_b(lnf(X)) = %9.4f\n\n", opt_bfbeta, mean_logl);
  }

  /* close mcmc file */
  if (!opt_onlysummary)
    fclose(fp_mcmc);

  if (opt_onlysummary)
  {
    /* read file and correctly set opt_samples */
    opt_samples = getlinecount(opt_mcmcfile);
    if (opt_samples)
    {
      if ((opt_method == METHOD_00) || (opt_method == METHOD_10))
        --opt_samples;
    }

    if (opt_samples == 0)
      fatal("Found %ld samples in file %s",opt_samples,opt_mcmcfile);
    else
      fprintf(stdout,"Read %ld samples from file %s\n",opt_samples,opt_mcmcfile);
  }

  if (opt_print_genetrees)
  {
    for (i = 0; i < opt_locus_count; ++i)
      fclose(fp_gtree[i]);
    free(fp_gtree);
    if (opt_checkpoint)
      free(gtree_offset);
  }

  /* print summary using the MCMC file */
  if (opt_method == METHOD_10)          /* species delimitation */
  {
    delimit_summary(fp_out, stree);
    delimitations_fini();
    rj_fini();
    free(posterior);
  }
  else if (opt_method == METHOD_11)
  {
    mixed_summary(fp_out);
    delimitations_fini();
    rj_fini();
    free(pspecies);

    /* TODO: FINISH IT */
    //assert(0);
  }

  for (i = 0; i < opt_locus_count; ++i)
    locus_destroy(locus[i]);
  free(locus);

  /* deallocate gene trees */
  for (i = 0; i < opt_locus_count; ++i)
    gtree_destroy(gtree[i],NULL);
  free(gtree);

  /* if species tree inference, deallocate cloned gene trees */
  if (opt_est_stree)          /* species tree inference */
  {
    for (i = 0; i < opt_locus_count; ++i)
      gtree_destroy(gclones[i],NULL);
    free(gclones);
  }

  gtree_fini(opt_locus_count);

  if (opt_method == METHOD_00)
    allfixed_summary(fp_out,stree);

  /* TODO: Perhaps we do not need 'species_count' as it should be equivalent
     to 'ndspecies' for the A01 method */
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
  if (opt_est_stree)          /* species tree inference */
    stree_destroy(sclone,NULL);         /* destroy cloned species tree */
    
  stree_fini();

  /* summary for method 01 */
  if (opt_method == METHOD_01)          /* species tree inference */
  {
    assert(species_count > 0);

    stree_summary(fp_out,species_names,species_count);

    /* cleanup */
    for (i = 0; i < species_count; ++i)
      free(species_names[i]);
    free(species_names);
  }

  if (opt_est_theta)
    free(mean_theta);
  free(mean_tau);

  if (opt_diploid)
    free(opt_diploid);

  /* close output file */
  fclose(fp_out);
}
