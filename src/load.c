/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#define PROP_COUNT 5
#define GTR_PROP_COUNT 3
#define CLOCK_PROP_COUNT 5

#define LOAD(x,n,fp) (fread((void *)(x),sizeof(*(x)),n,fp) == (size_t)(n))

size_t chk_size_int;
size_t chk_size_long;
size_t chk_size_double;

static BYTE dummy[256] = {0};

static stree_t * stree;
static gtree_t ** gtree;
static locus_t ** locus;

static void alloc_gtree()
{
  long i,j;

  /* initialize gene trees */
  gtree = (gtree_t **)xmalloc((size_t)opt_locus_count * sizeof(gtree_t *));

  /* allocate gene tree structures */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree[i] = (gtree_t *)xmalloc(sizeof(gtree_t));
    gtree_t * gt = gtree[i];

    gt->tip_count = 0;

    /* get number of gene tree tips by looking in the number of incoming
       sequences in the species tree tip nodes */
    for (j = 0; j < stree->tip_count; ++j)
      gt->tip_count += stree->nodes[j]->seqin_count[i];

    gt->inner_count = gt->tip_count-1;
    gt->edge_count = 2*gt->tip_count-2;
    gt->root = NULL;

    gt->nodes = (gnode_t **)xmalloc((size_t)(2*gt->tip_count-1) *
                                    sizeof(gnode_t *));
    size_t alloc_size = MAX(4,gtree[i]->tip_count + gtree[i]->inner_count);
    /* TODO: Change to xmalloc for the first-touch numa policy */
    gt->travbuffer = (gnode_t **)xcalloc(alloc_size, sizeof(gnode_t *));
    gt->migcount = NULL;
    gt->migpops = NULL;
    gt->rb_linked = NULL;
    
    for (j = 0; j < gt->tip_count + gt->inner_count; ++j)
    {
      gt->nodes[j] = (gnode_t *)xmalloc(sizeof(gnode_t));

      gt->nodes[j]->node_index = j;
      gt->nodes[j]->mi = NULL;

      if (stree->hybrid_count)
        gt->nodes[j]->hpath = (int *)xmalloc((size_t)(stree->hybrid_count) *
                                             sizeof(int));
      else
        gt->nodes[j]->hpath = NULL;

    }
  }
}

void load_chk_header(FILE * fp)
{
  long version_major;
  long version_minor;
  long version_patch;
  long version_chkp;
  BYTE magic[BPP_MAGIC_BYTES];
  BYTE buffer[16];

  if (!LOAD(magic,BPP_MAGIC_BYTES,fp))
    fatal(" Magic header mismatch");

  if (!LOAD(buffer,16,fp))
    fatal("Cannot read version information");


  version_major = (buffer[0]  <<  0) | 
                  (buffer[1]  <<  8) |
                  (buffer[2]  << 16) |
                  (buffer[3]  << 24);

  version_minor = (buffer[4]  <<  0) | 
                  (buffer[5]  <<  8) |
                  (buffer[6]  << 16) |
                  (buffer[7]  << 24);
 
  version_patch = (buffer[8]  <<  0) | 
                  (buffer[9]  <<  8) |
                  (buffer[10] << 16) |
                  (buffer[11] << 24);

  version_chkp  = (buffer[12] <<  0) | 
                  (buffer[13] <<  8) |
                  (buffer[14] << 16) |
                  (buffer[15] << 24);

  #if 0
  printf(" Magic: %c%c%c%c\n", magic[0],magic[1],magic[2],magic[3]);
  printf(" Version: %ld.%ld.%ld\n", version_major, version_minor, version_patch);
  printf(" CHKP version: %ld\n", version_chkp);
  #endif

  if (memcmp(magic,BPP_MAGIC,BPP_MAGIC_BYTES))
    fatal("File %s is not a BPP checkpoint file...", opt_resume);

  if ((version_major != VERSION_MAJOR) || (version_minor != VERSION_MINOR) || (version_patch != VERSION_PATCH))
    fatal("Incompatible CHKP: Checkpoint file version %ld, BPP version %ld",
          version_chkp, VERSION_CHKP);

  if (!LOAD(buffer,3,fp))
    fatal("Cannot read data type sizes");

  chk_size_int    = buffer[0];
  chk_size_long   = buffer[1];
  chk_size_double = buffer[2];

  /* TODO: Note that MSVC supports the %z only from verison 2015 (cl v19) */
  #if 0
  printf(" Type sizes: int(%zu) long(%zu) double(%zu)\n",
         chk_size_int, chk_size_long, chk_size_double);
  #endif

  if (sizeof(int) != chk_size_int)
    fatal("Mismatching int size");
  if (sizeof(long) != chk_size_long)
    fatal("Mismatching long size");
  if (sizeof(double) != chk_size_double)
    fatal("Mismatching double size");

  unsigned int sections;
  unsigned long size_section;

  if (!LOAD(&opt_threads,1,fp))
    fatal("Cannot read number of threads");
  if (!LOAD(&opt_threads_start,1,fp))
    fatal("Cannot read first thread slot");
  if (!LOAD(&opt_threads_step,1,fp))
    fatal("Cannot read thread stepping");

  /* Pin master thread for NUMA first touch policy */
  if (opt_threads > 1)
    threads_pin_master();

  
  unsigned int * rng = (unsigned int *)xmalloc((size_t)opt_threads *
                                               sizeof(unsigned int));
  if (!LOAD(rng,opt_threads,fp))
    fatal("Cannot read RNG states");
  set_legacy_rndu_array(rng);

  if (!LOAD(&sections,1,fp))
    fatal("Cannot read number of sections");
  #if 0
  printf(" Sections: %u\n", sections);
  #endif

  if (!LOAD(&size_section,1,fp))
    fatal("Cannot read number of sections");
  #if 0
  printf("   Section 1: %ld bytes\n\n", size_section);
  #endif
}

int load_string(FILE * fp, char ** buffer)
{
  size_t alloc = 0;
  size_t size = 0;
  size_t increment = 100;
  int c = 0;

  char * s = (char *)xmalloc(increment*sizeof(char));
  alloc = increment;

  while ((c=fgetc(fp)) && c != EOF)
  {
    s[size++] = (char)c;

    if (size == alloc)
    {
      char * temp = (char *)xmalloc((alloc+increment)*sizeof(char));
      memcpy(temp,s,size*sizeof(char));
      free(s);
      s = temp;
    }
  }
  s[size] = 0;

  char * temp = (char *)xmalloc((size+1)*sizeof(char));
  memcpy(temp,s,(size+1)*sizeof(char));
  free(s);

  *buffer = temp;

  if ((c == EOF) || size == 0) return 0;

  return 1;
}

static void print_filepaths()
{
  fprintf(stdout,
          "Control file     : %s\n"
          "Sequence file    : %s\n"
          "Map file         : %s\n"
          "MCMC sample file : %s\n"
          "Output file      : %s\n\n",
          opt_cfile, opt_msafile, opt_mapfile, opt_mcmcfile, opt_outfile);
}

static void load_chk_section_1(FILE * fp,
                               double ** pjump,
                               unsigned long * curstep,
                               long * ft_round,
                               long * ndspecies,
                               long * mcmc_offset,
                               long * out_offset,
                               long ** gtree_offset,
                               long ** rates_offset,
                               long * dparam_count,
                               double ** posterior,
                               double ** pspecies,
                               long * ft_round_rj,
                               double * pjump_rj,
                               long * ft_round_spr,
                               long * ft_round_snl,
                               long * pjump_spr,
                               long * pjump_snl,
                               double * mean_logl,
                               long ** mean_mrate_row,
                               long ** mean_mrate_col,
                               long ** mean_mrate_round,
                               double ** mean_mrate,
                               double ** mean_tau,
                               double ** mean_theta,
                               double ** mean_phi,
                               long * mean_mrate_count,
                               long * mean_tau_count,
                               long * mean_theta_count,
                               long * mean_phi_count,
                               int * prec_logpg,
                               int * prec_logl)
{
  long i;
  long total_nodes;
  char ** labels;

  *gtree_offset = NULL;
  *rates_offset = NULL;

  if (!LOAD(&opt_seed,1,fp))
    fatal("Cannot read seed");

  /* read control file name */
  if (!load_string(fp,&opt_cfile))
    fatal("Cannot read name of control file");

  /* read MSA filename */
  if (!load_string(fp,&opt_msafile))
    fatal("Cannot read name of MSA file");

  /* read constraint filename */
  if (!LOAD(&opt_constraint_count,1,fp))
    fatal("Cannot read number of constraints");
  if (opt_constraint_count)
  {
    if (!load_string(fp,&opt_constraintfile))
      fatal("Cannot read name of constraint file");
    printf(" Constraints file: %s\n", opt_constraintfile);
  }

  long mapfile_present;
  opt_mapfile = NULL;
  if (!LOAD(&mapfile_present,1,fp))
    fatal("Cannot read mapfile information");

  if (mapfile_present)
  {
    /* read imap filename */
    if (!load_string(fp,&opt_mapfile))
      fatal("Cannot read name of map file");
  }

  /* read output filename */
  if (!load_string(fp,&opt_outfile))
    fatal("Cannot read name of output file");

  /* read output filename */
  if (!load_string(fp,&opt_mcmcfile))
    fatal("Cannot read name of mcmc file");

  print_filepaths();

  /* read checkpoint info */
  if (!LOAD(&opt_checkpoint,1,fp))
    fatal("Cannot read 'checkpoint' flag");
  if (!LOAD(&opt_checkpoint_current,1,fp))
    fatal("Cannot read 'checkpoint' status");
  if (!LOAD(&opt_checkpoint_initial,1,fp))
    fatal("Cannot read 'checkpoint' tag initial value");
  if (!LOAD(&opt_checkpoint_step,1,fp))
    fatal("Cannot read 'checkpoint' tag step value");

  /* read network info */
  if (!LOAD(&opt_msci,1,fp))
    fatal("Cannot read species network flag");

  if (!LOAD(&opt_migration,1,fp))
    fatal("Cannot read migration flag");

  /* read method info */
  if (!LOAD(&opt_method,1,fp))
    fatal("Cannot read species network flag");

  /* read speciesdelimitation */
  if (!LOAD(&opt_est_delimit,1,fp))
    fatal("Cannot read 'speciesdelimitation' tag");

  if (!LOAD(&opt_rjmcmc_method,1,fp))
    fatal("Cannot read 'speciesdelimitation' tag");

  if (opt_est_delimit && opt_rjmcmc_method != 0 && opt_rjmcmc_method != 1)
    fatal("rj-MCMC method can be either 0 or 1, but found %ld", opt_rjmcmc_method);

  if (opt_rjmcmc_method == 0)
  {
    if (!LOAD(&opt_rjmcmc_epsilon,1,fp))
      fatal("Cannot read 'speciesdelimitation' tag");
    if (!LOAD(dummy,sizeof(double),fp))
      fatal("Cannot read 'speciesdelimitation' tag");
    #if 0
    if (opt_est_delimit)
      printf(" Speciesdelimitation: %ld %ld %f\n", opt_est_delimit, opt_rjmcmc_method, opt_rjmcmc_epsilon);
    else
      printf(" Speciesdelimitation: Disabled\n");
    #endif
  }
  else
  {
    if (!LOAD(&opt_rjmcmc_alpha,1,fp))
      fatal("Cannot read 'speciesdelimitation' tag");
    if (!LOAD(&opt_rjmcmc_mean,1,fp))
      fatal("Cannot read 'speciesdelimitation' tag");
    #if 0
    if (opt_est_delimit)
      printf(" Speciesdelimitation: %ld %ld %f %f\n", opt_est_delimit,
             opt_rjmcmc_method, opt_rjmcmc_alpha, opt_rjmcmc_mean);
    else
      printf(" Speciesdelimitation: Disabled\n");
    #endif
  }

  /* read speciestree */
  if (!LOAD(&opt_est_stree,1,fp))
    fatal("Cannot read 'speciestree' tag");
  #if 0
  printf(" Speciestree: %ld\n", opt_est_stree);
  #endif

  if (!LOAD(dummy,3*sizeof(double),fp))
    fatal("Cannot read 'speciestree' tag");

  /* read speciesmodelprior */
  if (!LOAD(&opt_delimit_prior,1,fp))
    fatal("Cannot read 'speciesmodelprior' tag");
  #if 0
  printf(" Speciesmodelprior: %ld\n", opt_delimit_prior);
  #endif

  /* read species&tree */
  unsigned int stree_tip_count;
  unsigned int stree_inner_count;
  unsigned int stree_hybrid_count;
  unsigned int stree_edge_count;
  if (!LOAD(&stree_tip_count,1,fp))
    fatal("Cannot read 'species&tree' tag");
  if (!LOAD(&stree_inner_count,1,fp))
    fatal("Cannot read 'species&tree' tag");
  if (!LOAD(&stree_hybrid_count,1,fp))
    fatal("Cannot read 'species&tree' tag");
  if (!LOAD(&stree_edge_count,1,fp))
    fatal("Cannot read 'species&tree' tag");

  labels = (char **)xmalloc((size_t)(stree_tip_count+stree_hybrid_count) *
                            sizeof(char *));
  for (i = 0; i < stree_tip_count; ++i)
  {
    if (!load_string(fp,labels+i))
      fatal("Cannot read species labels for 'species&tree' tag");
  }
  for (i = 0; i < stree_hybrid_count; ++i)
  {
    if (!load_string(fp,labels+stree_tip_count+i))
      fatal("Cannot read species labels for 'species&tree' tag");
  }
  #if 0
  printf(" Species&tree: %u species (", stree_tip_count);
  for (i = 0; i < stree_tip_count; ++i)
    printf(" %s", labels[i]);
  printf(")\n");
  if (stree_hybrid_count)
  {
    printf("Hybridizations:\n");
    for (i = 0; i < stree_hybrid_count; ++i)
      printf(" %s\n", labels[stree_tip_count+i]);
    printf("\n");
  }
  #endif

  /* read usedata, cleandata and nloci */
  if (!LOAD(&opt_usedata,1,fp))
    fatal("Cannot read 'usedata' tag");
  #if 0
  printf(" usedata: %ld\n", opt_usedata);
  #endif
  if (!LOAD(&opt_cleandata,1,fp))
    fatal("Cannot read 'cleandata' tag");
  #if 0
  printf(" cleandata: %ld\n", opt_cleandata);
  #endif
  if (!LOAD(&opt_locus_count,1,fp))
    fatal("Cannot read 'nloci' tag");
  #if 0
  printf(" nloci: %ld\n", opt_locus_count);
  #endif

  /* load print flags */
  if (!LOAD(&opt_print_samples,1,fp))
    fatal("Cannot read print flags");
  if (!LOAD(&opt_print_locusrate,1,fp))
    fatal("Cannot read print flags");
  if (!LOAD(&opt_print_hscalars,1,fp))
    fatal("Cannot read print flags");
  if (!LOAD(&opt_print_genetrees,1,fp))
    fatal("Cannot read print flags");
  if (!LOAD(&opt_print_rates,1,fp))
    fatal("Cannot read rates flags");
  if (!LOAD(&opt_print_qmatrix,1,fp))
    fatal("Cannot read qmatrix flag");
  if (!LOAD(&opt_print_locusfile,1,fp))
    fatal("Cannot read locusfile flag");
  if (opt_print_samples == 0)
    fatal("Corrupted checkfpoint file, opt_print_samples=0");

  /* read theta prior */
  if (!LOAD(&opt_theta_dist,1,fp))
    fatal("Cannot read type of theta prior");
  if (!LOAD(&opt_theta_alpha,1,fp))
    fatal("Cannot read alpha of 'theta' tag");
  if (!LOAD(&opt_theta_beta,1,fp))
    fatal("Cannot read beta 'theta' tag");
  if (!LOAD(&opt_theta_p,1,fp))
    fatal("Cannot read p of 'theta' tag");
  if (!LOAD(&opt_theta_q,1,fp))
    fatal("Cannot read q of 'theta' tag");
  if (!LOAD(&opt_theta_min,1,fp))
    fatal("Cannot read min of 'theta' tag");
  if (!LOAD(&opt_theta_max,1,fp))
    fatal("Cannot read max of 'theta' tag");
  if (!LOAD(&opt_est_theta,1,fp))
    fatal("Cannot read est 'theta' tag");
  if (!LOAD(&opt_linkedtheta,1,fp))
    fatal("Cannot read linked theta tag");
  #if 0
  printf(" theta: %f %f %ld\n", opt_theta_alpha, opt_theta_beta, opt_est_theta);
  #endif
  
  /* read tau prior */
  if (!LOAD(&opt_tau_dist,1,fp))
    fatal("Cannot read type of tau prior");
  if (!LOAD(&opt_tau_alpha,1,fp))
    fatal("Cannot read alpha of 'theta' tag");
  if (!LOAD(&opt_tau_beta,1,fp))
    fatal("Cannot read beta 'theta' tag");
  #if 0
  printf(" tau: %f %f\n", opt_tau_alpha, opt_tau_beta);
  #endif

  /* laod phi prior */
  if (!LOAD(&opt_phi_alpha,1,fp))
    fatal("Cannot read alpha of 'phiprior' tag"); 
  if (!LOAD(&opt_phi_beta,1,fp))
    fatal("Cannot read beta of 'phiprior' tag"); 

  /* laod migration rates prior */
  if (!LOAD(&opt_mig_alpha,1,fp))
    fatal("Cannot read alpha of 'migprior' tag"); 
  if (!LOAD(&opt_mig_beta,1,fp))
    fatal("Cannot read beta of 'migprior' tag"); 

  if (!LOAD(&opt_model,1,fp))
    fatal("Cannot read substitution model information");

  if (!LOAD(&opt_alpha_cats,1,fp))
    fatal("Cannot read number of gamma categories");
  if (!LOAD(&opt_alpha_alpha,1,fp))
    fatal("Cannot read alpha parameter of gamma rates");
  if (!LOAD(&opt_alpha_beta,1,fp))
    fatal("Cannot read beta parameter of gamma rates");

  /* load locus rate estimation flag */
  if (!LOAD(&opt_est_locusrate,1,fp))
    fatal("Cannot read locusrate tag"); 

  if (!LOAD(&opt_est_mubar,1,fp))
    fatal("Cannot read mubar estimation flag");

  /* load heredity scalers estimation flag */
  if (!LOAD(&opt_est_heredity,1,fp))
    fatal("Cannot read heredity tag"); 

  /* load heredity scalers alpha */
  if (!LOAD(&opt_heredity_alpha,1,fp))
    fatal("Cannot read heredity alpha"); 

  /* load heredity scalers beta */
  if (!LOAD(&opt_heredity_beta,1,fp))
    fatal("Cannot read heredity beta"); 

  /* load clock and locusrate info */
  if (!LOAD(&opt_clock,1,fp))
    fatal("Cannot read 'clock' tag");
  if (!LOAD(&opt_mubar_alpha,1,fp))
    fatal("Cannot read 'mubar_alpha'");
  if (!LOAD(&opt_mubar_beta,1,fp))
    fatal("Cannot read 'mubar_beta'");
  if (!LOAD(&opt_mui_alpha,1,fp))
    fatal("Cannot read 'mui_alpha'");
  if (!LOAD(&opt_vbar_alpha,1,fp))
    fatal("Cannot read 'vbar_alpha'");
  if (!LOAD(&opt_vbar_beta,1,fp))
    fatal("Cannot read 'vbar_beta'");
  if (!LOAD(&opt_vi_alpha,1,fp))
    fatal("Cannot read 'vi_alpha'");
  if (!LOAD(&opt_rate_prior,1,fp))
    fatal("Cannot read rate prior");
  if (!LOAD(&opt_locusrate_prior,1,fp))
    fatal("Cannot read locus rate prior");

  /* read finetune */
  if (!LOAD(&opt_finetune_reset,1,fp))
    fatal("Cannot read 'finetune' tag");
  if (!LOAD(&opt_finetune_migrates,1,fp))
    fatal("Cannot read migration rates finetune parameter");
  if (!LOAD(&opt_finetune_phi,1,fp))
    fatal("Cannot read gene tree phi finetune parameter");
  if (!LOAD(&opt_finetune_gtage,1,fp))
    fatal("Cannot read gene tree age finetune parameter");
  if (!LOAD(&opt_finetune_gtspr,1,fp))
    fatal("Cannot read gene tree SPR finetune parameter");
  if (!LOAD(&opt_finetune_theta,1,fp))
    fatal("Cannot read species tree theta finetune parameter");
  if (!LOAD(&opt_finetune_tau,1,fp))
    fatal("Cannot read species tree tau finetune parameter");
  if (!LOAD(&opt_finetune_mix,1,fp))
    fatal("Cannot read species mixing step finetune parameter");
  if (!LOAD(&opt_finetune_locusrate,1,fp))
    fatal("Cannot read species locusrate/heredity finetune parameter");
  if (!LOAD(&opt_finetune_qrates,1,fp))
    fatal("Cannot read qmatrix rates finetune parameter");
  if (!LOAD(&opt_finetune_freqs,1,fp))
    fatal("Cannot read base frequencies finetune parameter");
  if (!LOAD(&opt_finetune_alpha,1,fp))
    fatal("Cannot read alpha value for gamma rate variation");
  if (!LOAD(&opt_finetune_mubar,1,fp))
    fatal("Cannot read mubar finetune parameter");
  if (!LOAD(&opt_finetune_mui,1,fp))
    fatal("Cannot read mui finetune parameter");
  if (!LOAD(&opt_finetune_nubar,1,fp))
    fatal("Cannot read nubar finetune parameter");
  if (!LOAD(&opt_finetune_nui,1,fp))
    fatal("Cannot read nui finetune parameter");
  if (!LOAD(&opt_finetune_branchrate,1,fp))
    fatal("Cannot read branchrate finetune parameter");

  if (!LOAD(&opt_max_species_count,1,fp))
    fatal("Cannot read max number of species");

  if (!LOAD(&opt_prob_snl,1,fp))
    fatal("Cannot read frequency for SNL");
  if (!LOAD(&opt_prob_snl_shrink,1,fp))
    fatal("Cannot read frequency for SNL shrink");
  if (!LOAD(&opt_snl_lambda_expand,1,fp))
    fatal("Cannot read lambda for SNL expand");
  if (!LOAD(&opt_snl_lambda_shrink,1,fp))
    fatal("Cannot read lambda for SNL shrink");

  if (!LOAD(&opt_pseudop_exist,1,fp))
    fatal("Cannot read information on pseudo priors");
  if (!LOAD(&opt_mig_vrates_exist,1,fp))
    fatal("Cannot read information on variable migration rates");

  #if 0
  printf(" Current finetune: %ld: %f %f %f %f %f",
         opt_finetune_reset, opt_finetune_gtage, opt_finetune_gtspr,
         opt_finetune_theta, opt_finetune_tau,opt_finetune_mix);
  if (opt_est_locusrate == MUTRATE_ESTIMATE || opt_est_heredity)
    printf(" %f\n", opt_finetune_locusrate);
  else
    printf("\n");
  #endif

  double stree_locusrate_mubar = 0;
  double stree_locusrate_nubar = 0;
  double stree_nui_sum = 0;

  if (!LOAD(&stree_locusrate_mubar,1,fp))
    fatal("Cannot read locusrate_mubar value");
  if (!LOAD(&stree_locusrate_nubar,1,fp))
    fatal("Cannot read locusrate_nubar value");
  if (!LOAD(&stree_nui_sum,1,fp))
    fatal("Cannot read nui_sum value");

  /* read diploid */
  opt_diploid = (long *)xmalloc((size_t)stree_tip_count*sizeof(char *));
  if (!LOAD(opt_diploid,stree_tip_count,fp))
    fatal("Cannot read 'diploid' tag");

  for  (i = 0; i < stree_tip_count; ++i)
    if (opt_diploid[i]) break;

  if (i == stree_tip_count)  /* no diploids */
  {
    free(opt_diploid);
    opt_diploid = NULL;
    #if 0
    printf(" Diploid: None found\n");
    #endif
  }
  #if 0
  else
  {
    printf(" Diploid:");
    for (i = 0; i < stree_tip_count; ++i)
      printf(" %ld", opt_diploid[i]);
    printf("\n");
  }
  #endif
  if (!LOAD(&opt_diploid_size,1,fp))
    fatal("Cannot load diploid size variable");

  /* read MCMC run info  */
  if (!LOAD(&opt_burnin,1,fp))
    fatal("Cannot read 'burnin' tag");
  if (!LOAD(&opt_samplefreq,1,fp))
    fatal("Cannot read 'sampfreq' tag");
  if (!LOAD(&opt_samples,1,fp))
    fatal("Cannot read 'nsample' tag");
  if (!LOAD(curstep,1,fp))
    fatal("Cannot read current MCMC step");
  if (!LOAD(ft_round,1,fp))
    fatal("Cannot read current finetune round");
  if (!LOAD(ndspecies,1,fp))
    fatal("Cannot read number of delimited species");

  size_t pjump_size = PROP_COUNT + 1+1 + GTR_PROP_COUNT + CLOCK_PROP_COUNT + !!opt_migration;
  *pjump = (double *)xmalloc(pjump_size*sizeof(double));

  if (!LOAD(*pjump,pjump_size,fp))
    fatal("Cannot read pjump");

  if (!LOAD(mcmc_offset,1,fp))
    fatal("Cannot read MCMC file offset");

  if (!LOAD(out_offset,1,fp))
    fatal("Cannot read output file offset");

  if (!LOAD(&opt_bfbeta,1,fp))
    fatal("Cannot read bfbeta");

  if (opt_print_genetrees)
  {
    *gtree_offset = (long *)xmalloc((size_t)opt_locus_count*sizeof(long));
    if (!LOAD(*gtree_offset,opt_locus_count,fp))
      fatal("Cannot read gtree file offsets");
  }

  if (opt_print_locusfile)
  {
    *rates_offset = (long *)xmalloc((size_t)opt_locus_count*sizeof(long));
    if (!LOAD(*rates_offset,opt_locus_count,fp))
      fatal("Cannot read rates files offsets");
  }

  if (!LOAD(dparam_count,1,fp))
    fatal("Cannot read dparam_count");

  long dmodels_count;
  if (!LOAD(&dmodels_count,1,fp))
    fatal("Cannot read dmodels_count");

  *posterior = NULL;
  if (dmodels_count)
  {
    *posterior = (double *)xmalloc((size_t)dmodels_count*sizeof(double));
    if (!LOAD(*posterior,dmodels_count,fp))
      fatal("Cannot read posterior");
  }

  *pspecies = NULL;
  if (opt_est_stree && opt_est_delimit)
  {
    *pspecies = (double *)xmalloc((size_t)opt_max_species_count *
                                  sizeof(double));
    if (!LOAD(*pspecies,opt_max_species_count,fp))
      fatal("Cannot read pspecies");
  }

  if (!LOAD(ft_round_rj,1,fp))
    fatal("Cannot read RJ finetune round"); 

  if (!LOAD(pjump_rj,1,fp))
    fatal("Cannot read RJ pjump"); 

  if (!LOAD(ft_round_spr,1,fp))
    fatal("Cannot read finetune round for species tree SPR"); 

  if (!LOAD(ft_round_snl, 1, fp))
    fatal("Cannot read finetune round for species tree SNL");

  if (!LOAD(pjump_spr, 1, fp))
    fatal("Cannot read species tree SPR pjump");

  if (!LOAD(pjump_snl,1,fp))
    fatal("Cannot read species tree SNL pjump"); 

  if (!LOAD(mean_logl,1,fp))
    fatal("Cannot read mean logl"); 

  *mean_mrate = NULL;
  *mean_mrate_row = NULL;
  *mean_mrate_col = NULL;
  *mean_mrate_round = NULL;
  *mean_mrate_count = 0;
  if (opt_migration)
  {
    if (!LOAD(mean_mrate_count,1,fp))
      fatal("Cannot read number of mean migration rates");

    *mean_mrate = (double *)xmalloc((size_t)(*mean_mrate_count)*sizeof(double));
    *mean_mrate_row = (long *)xmalloc((size_t)(*mean_mrate_count)*sizeof(long));
    *mean_mrate_col = (long *)xmalloc((size_t)(*mean_mrate_count)*sizeof(long));
    *mean_mrate_round = (long *)xmalloc((size_t)(*mean_mrate_count)*sizeof(long));

    if (!LOAD(*mean_mrate, *mean_mrate_count, fp))
      fatal("Cannot read mean migration rate values");

    if (!LOAD(*mean_mrate_row, *mean_mrate_count, fp))
      fatal("Cannot read mean migration rate row values");

    if (!LOAD(*mean_mrate_col, *mean_mrate_count, fp))
      fatal("Cannot read mean migration rate col values");

    if (!LOAD(*mean_mrate_round, *mean_mrate_count, fp))
      fatal("Cannot read mean migration rate round values");

  }
  if (!LOAD(mean_tau_count,1,fp))
    fatal("Cannot read number of mean taus"); 

  if (opt_est_theta)
  {
    if (!LOAD(mean_theta_count,1,fp))
      fatal("Cannot read number of mean thetas"); 
  }

  if (!LOAD(mean_phi_count,1,fp))
    fatal("Cannot read number of mean phis");

  *mean_tau   = (double *)xmalloc((size_t)(*mean_tau_count)*sizeof(double));

  if (opt_est_theta)
  {
    *mean_theta = (double *)xmalloc((size_t)(*mean_theta_count)*sizeof(double));
  }

  if (!LOAD(*mean_tau,*mean_tau_count,fp))
    fatal("Cannot read mean tau values"); 

  if (opt_est_theta)
  {
    if (!LOAD(*mean_theta,*mean_theta_count,fp))
      fatal("Cannot read mean theta values"); 
  }

  *mean_phi = NULL;
  if (*mean_phi_count)
  {
    *mean_phi = (double *)xmalloc((size_t)(*mean_phi_count)*sizeof(double));
    
    if (!LOAD(*mean_phi,*mean_phi_count,fp))
      fatal("Cannot load mean phi value");
  }


  if (!LOAD(prec_logpg,1,fp))
    fatal("Cannot read logPG digits precision");

  if (!LOAD(prec_logl,1,fp))
    fatal("Cannot read logL digits precision");

  if (!LOAD(&opt_load_balance,1,fp))
    fatal("Cannot read load balance scheme");

  if (opt_threads > 1)
  {
    thread_info_t * ti = (thread_info_t *)xmalloc((size_t)opt_threads *
                                                  sizeof(thread_info_t));
    for (i = 0; i < opt_threads; ++i)
    {
      thread_info_t * tip = ti+i;
      if (!LOAD(&(tip->locus_first),1,fp))
        fatal("Cannot load thread_info");
      if (!LOAD(&(tip->locus_count),1,fp))
        fatal("Cannot load thread_info");
    }
    threads_set_ti(ti);
  }

  if (opt_migration)
  {
    assert(stree_hybrid_count == 0);
    unsigned int total_nodes = stree_tip_count+stree_inner_count;

    opt_migration_matrix = (long **)xmalloc((size_t)total_nodes*sizeof(long *));
    for (i = 0; i < total_nodes; ++i)
    {
     opt_migration_matrix[i] = (long*)xmalloc((size_t)total_nodes*sizeof(long));
     if (!LOAD(opt_migration_matrix[i],total_nodes,fp))
       fatal("Cannot load migration matrix");
    }

    opt_mig_bitmatrix = (long **)xmalloc((size_t)total_nodes * sizeof(long *));
    for (i = 0; i < total_nodes; ++i)
    {
      opt_mig_bitmatrix[i] = (long *)xmalloc((size_t)total_nodes*sizeof(long));
      if (!LOAD(opt_mig_bitmatrix[i],total_nodes,fp))
        fatal("Cannot load migration bitmatrix");
    }

    opt_mig_specs = (migspec_t *)xcalloc((size_t)opt_migration,sizeof(migspec_t));
    for (i = 0; i < opt_migration; ++i)
    {
      migspec_t * spec = opt_mig_specs+i;
      if (!load_string(fp,&(spec->source)))
        fatal("Cannot read list of migrations (migration %ld - source)", i+1);
      if (!load_string(fp,&(spec->target)))
        fatal("Cannot read list of migrations (migration %ld - target)", i+1);
      if (!LOAD(&(spec->si), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - si)", i+1);
      if (!LOAD(&(spec->ti), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - ti)", i+1);
      if (!LOAD(&(spec->am), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - am)", i+1);
      if (!LOAD(&(spec->alpha), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - alpha)", i+1);
      if (!LOAD(&(spec->beta), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - beta)", i+1);
      if (!LOAD(&(spec->pseudo_a), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - pseudo_a)", i+1);
      if (!LOAD(&(spec->pseudo_b), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - pseudo_b)", i+1);
      if (!LOAD(&(spec->params), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - params)", i+1);
      if (!LOAD(&(spec->M), 1, fp))
        fatal("Cannot read list of migrations (migration %ld - M)", i+1);
      if (spec->params == 1 || spec->params == 3 || spec->params == 5)
      {
        spec->Mi = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
        if (!LOAD(spec->Mi, opt_locus_count, fp))
          fatal("Cannot read list of migrations (migration %ld - Mi)", i+1);
      }
    }
  }

  #if 0
  fprintf(stdout, " Burnin: %ld\n", opt_burnin);
  fprintf(stdout, " Sampfreq: %ld\n", opt_samplefreq);
  fprintf(stdout, " Nsample: %ld\n", opt_samples);
  fprintf(stdout, " Next step: %ld\n\n", *curstep);
  #endif

  /* populate species tree */
  //stree = (stree_t *)xcalloc(1,sizeof(stree_t));
  stree = (stree_t *)xcalloc(1,sizeof(stree_t));
  
  stree->tip_count = stree_tip_count;
  stree->inner_count = stree_inner_count;
  stree->hybrid_count = stree_hybrid_count;
  stree->edge_count = stree_edge_count;
  stree->locus_count = opt_locus_count;
  stree->locusrate_mubar = stree_locusrate_mubar;
  stree->locusrate_nubar = stree_locusrate_nubar;
  stree->nui_sum = stree_nui_sum;

  stree->mi_tbuffer = NULL;
  if (opt_migration)
    stree->mi_tbuffer = (miginfo_t **)xcalloc((size_t)opt_threads,
                                              sizeof(miginfo_t *));

  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;
  stree->nodes = (snode_t **)xmalloc((size_t)total_nodes*sizeof(snode_t *));
  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i] = (snode_t *)xcalloc(1,sizeof(snode_t));

  stree->td = (snode_t **)xmalloc((size_t)total_nodes * sizeof(snode_t *));

  /* set tip labels */
  for (i = 0; i < stree->tip_count; ++i)
    stree->nodes[i]->label = labels[i];

  unsigned int hoffset = stree->tip_count + stree->inner_count;
  for (i = 0; i < stree->hybrid_count; ++i)
    stree->nodes[hoffset+i]->label = labels[stree->tip_count+i];
  free(labels);

  for (i = 0; i < stree->inner_count; ++i)
    stree->nodes[stree->tip_count+i]->label = NULL;

  for (i = 0; i < total_nodes; ++i)
  {
    stree->nodes[i]->node_index = i;
    stree->nodes[i]->data = NULL;
    stree->nodes[i]->hybrid = NULL;
  }

  /* allocate coalescent events */
  for (i = 0; i < total_nodes; ++i)
  {
    snode_t * node = stree->nodes[i];

    node->event_count = (int *)xmalloc((size_t)opt_locus_count * sizeof(int));
    node->seqin_count = (int *)xmalloc((size_t)opt_locus_count * sizeof(int));
    node->gene_leaves = (unsigned int *)xmalloc((size_t)opt_locus_count *
                                                sizeof(unsigned int));
    node->logpr_contrib = (double *)xmalloc((size_t)opt_locus_count *
                                            sizeof(double));
    node->old_logpr_contrib = (double *)xmalloc((size_t)opt_locus_count *
                                                sizeof(double));
    node->event = (dlist_t **)xmalloc((size_t)opt_locus_count *
                                      sizeof(dlist_t *));

    node->t2h = NULL;
    node->old_t2h = NULL;
    node->hphi_sum = 0;
    node->notheta_phi_contrib = NULL;
    node->notheta_old_phi_contrib = NULL;
    if (!opt_est_theta)
    {
      node->t2h = (double *)xcalloc((size_t)opt_locus_count,sizeof(double));
      node->old_t2h = (double *)xcalloc((size_t)opt_locus_count,sizeof(double));
      node->t2h_sum = 0;
      node->event_count_sum = 0;
    }

    node->brate = NULL;  /* will be allocated during reading */
  }

  /* allocate hx */
  if (opt_msci)
  {
    for (i = 0; i < total_nodes; ++i)
    {
      stree->nodes[i]->hx = (long *)xcalloc((size_t)opt_threads,sizeof(long));
    }
  }

  /* allocate re-entrant marks */
  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->mark = (int *)xcalloc((size_t)opt_threads,sizeof(int));

  stree->pptable = (int**)xcalloc((size_t)total_nodes,sizeof(int *));
  for (i = 0; i < total_nodes; ++i)
    stree->pptable[i] = (int *)xcalloc((size_t)total_nodes,sizeof(int));
}

static int cb_ascint(const void * a, const void * b)
{
  const unsigned int * x = a;
  const unsigned int * y = b;

  if (*x > *y) return 1;
  return -1;
}

void load_chk_section_2(FILE * fp)
{
  unsigned int total_nodes;
  unsigned int hoffset;
  long i,j,k;
  unsigned int * hindices; 

  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;
  hoffset = stree->tip_count + stree->inner_count;

  /* reset parent nodes to NULL */
  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->parent = NULL;

  /* reset child nodes to NULL for tips */
  for (i = 0; i < stree->tip_count; ++i)
    stree->nodes[i]->left = stree->nodes[i]->right = NULL;

  hindices = (unsigned int *)xcalloc((size_t)(stree->hybrid_count),
                                     sizeof(unsigned int));

  /* load node indices of hybridization events */
  if (!LOAD(hindices,stree->hybrid_count,fp))
    fatal("Cannot read species tree topology");

  /* set hybridization links */
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    stree->nodes[hoffset+i]->hybrid = stree->nodes[hindices[i]];
    stree->nodes[hindices[i]]->hybrid = stree->nodes[hoffset+i];
    stree->nodes[hindices[i]]->label = xstrdup(stree->nodes[hoffset+i]->label);
  }

  /* sort hindices in ascending order */
  qsort(hindices,(size_t)stree->hybrid_count,sizeof(unsigned int), cb_ascint);

  /* set left child for each inner node according to checkpoint data */
  for (i = 0; i < stree->inner_count; ++i)
  {
    unsigned int left_child_index;
    if (!LOAD(&left_child_index,1,fp))
      fatal("Cannot read species tree topology");
    stree->nodes[stree->tip_count+i]->left = stree->nodes[left_child_index];
    stree->nodes[left_child_index]->parent = stree->nodes[stree->tip_count+i];
  }

  /* set right child for each inner node according to checkpoint data */
  for (i = 0; i < stree->inner_count; ++i)
  {
    unsigned int valid;
    unsigned int right_child_index;

    if (!LOAD(&valid,1,fp))
      fatal("Cannot read species tree topology");

    if (!valid) continue;

    if (!LOAD(&right_child_index,1,fp))
      fatal("Cannot read species tree topology");

    stree->nodes[stree->tip_count+i]->right = stree->nodes[right_child_index];
    stree->nodes[right_child_index]->parent = stree->nodes[stree->tip_count+i];
  }

  free(hindices);

  /* now set the root node */
  if (stree->tip_count == 1)
  {
    stree->root = stree->nodes[0];
  }
  else
  {
    unsigned nullparent_count = 0;
    for (i = 0; i < stree->inner_count; ++i)
    {
      if (!stree->nodes[stree->tip_count+i]->parent)
      {
        stree->root = stree->nodes[stree->tip_count+i];
        nullparent_count++;
      }
    }
    if (nullparent_count != 1)
      fatal("Erroneous species tree structure");
  }


  /* now check species tree consistency */
  for (i = 0; i < stree->tip_count; ++i)
    assert((stree->nodes[i]->left == stree->nodes[i]->right) &&
           (stree->nodes[i]->left == NULL));

  /* a lot of assertions */
  for (i = 0; i < stree->inner_count; ++i)
  {
    unsigned int index = stree->tip_count + i;

    if (stree->nodes[index]->left && stree->nodes[index]->right)
      assert(stree->nodes[index]->left != stree->nodes[index]->right);

    if (stree->nodes[index]->parent)
    {
      assert(stree->nodes[index]->left  != stree->nodes[index]->parent);
      assert(stree->nodes[index]->right != stree->nodes[index]->parent);
    }
    assert(stree->nodes[index]->left  != stree->nodes[index]);
    assert(stree->nodes[index]->right != stree->nodes[index]);

    if (stree->nodes[index]->left)
      assert(stree->nodes[index]->left->parent == stree->nodes[index]);
    if (stree->nodes[index]->right)
      assert(stree->nodes[index]->right->parent == stree->nodes[index]);

  }

  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->length = 0;

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    if (!LOAD(&(stree->nodes[hoffset+i]->hybrid->hphi),1,fp))
      fatal("Cannot read genetic contribution (phi) for node %ld", i);

    stree->nodes[hoffset+i]->hphi = 1 - stree->nodes[hoffset+i]->hybrid->hphi;
  }

  for (i = 0; i < total_nodes; ++i)
  {
    if (!LOAD(&(stree->nodes[i]->has_phi),1,fp))
      fatal("Cannot read has_phi for node %ld", i);
  }

  for (i = 0; i < total_nodes; ++i)
  {
    if (!LOAD(&(stree->nodes[i]->htau),1,fp))
      fatal("Cannot read parent htau for node %ld", i);
  }

  for (i = 0; i < total_nodes; ++i)
  {
    if (!LOAD(&(stree->nodes[i]->prop_tau),1,fp))
      fatal("Cannot read parent prop_tau for node %ld", i);
  }

  for (i = 0; i < total_nodes; ++i)
  {
    long ltheta_index;

    if (!LOAD(&ltheta_index,1,fp))
      fatal("Cannot read linked theta index for node %ld", i);     
    if (ltheta_index == -1)
      stree->nodes[i]->linked_theta = NULL;
    else
      stree->nodes[i]->linked_theta = stree->nodes[ltheta_index];
  }

  stree_label(stree);

  #if 0
  char * newick = stree_export_newick(stree->root,NULL);
  printf("Current species tree: %s\n", newick);
  free(newick);
  #endif

  /* read thetas */
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(&(stree->nodes[i]->theta),1,fp))
      fatal("Cannot read species nodes theta");
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(&(stree->nodes[i]->has_theta),1,fp))
      fatal("Cannot read species nodes has_theta");

  /* read taus */
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(&(stree->nodes[i]->tau),1,fp))
      fatal("Cannot read species nodes tau");

  /* read support */
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(&(stree->nodes[i]->support),1,fp))
      fatal("Cannot read species nodes support values");

  /* read constraints */
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(&(stree->nodes[i]->constraint),1,fp))
      fatal("Cannot read species nodes constraints");

  /* read constraints line numbers*/
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(&(stree->nodes[i]->constraint_lineno),1,fp))
      fatal("Cannot read species nodes constraints line numbers");

  /* read number of coalescent events */
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(stree->nodes[i]->event_count,opt_locus_count,fp))
        fatal("Cannot read species event counts");

  /* read branch rates */
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    unsigned int valid;

    for (i = 0; i < total_nodes; ++i)
    {
      if (!LOAD(&valid,1,fp))
        fatal("Cannot read branch rates");

      if (!valid) continue;

      snode_t * node = stree->nodes[i];
      node->brate = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
      if (!LOAD(node->brate,opt_locus_count,fp))
        fatal("Cannot read branch rates");
    }
  }

//  /* read MSC density contributions */
//  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
//    if (fread((void *)&(stree->nodes[i]->logpr_contrib),sizeof(double),1,fp) != 1)
//      fatal("Cannot read species MSC density contribution");
  
  if (!opt_est_theta)
  {
    if (!LOAD(&(stree->notheta_logpr),1,fp))
      fatal("Cannot read MSC density");

    if (!LOAD(&(stree->notheta_hfactor),1,fp))
      fatal("Cannot read heredity multiplier precomputed density contribution");

    if (!LOAD(&(stree->notheta_sfactor),1,fp))
      fatal("Cannot read sequence count precomputed density contribution");

    stree->notheta_old_logpr = 0;

    for (i = 0; i < total_nodes; ++i)
    {
      if (!LOAD(stree->nodes[i]->t2h,opt_locus_count,fp))
        fatal("Cannot read per-locus t2h contributions");

      if (!LOAD(&(stree->nodes[i]->t2h_sum),1,fp))
        fatal("Cannot read t2h sum");

      if (!LOAD(&(stree->nodes[i]->event_count_sum),1,fp))
        fatal("Cannot read coalescent events sum");

      if (!LOAD(&(stree->nodes[i]->notheta_logpr_contrib),1,fp))
        fatal("Cannot read MSC density contribution for node");

      stree->nodes[i]->notheta_old_logpr_contrib = 0;
    }
    if (opt_msci)
    {
      for (i = 0; i < stree->hybrid_count; ++i)
      {
        unsigned int index = stree->tip_count+stree->inner_count;
        snode_t * x = stree->nodes[index+i];

        x->notheta_phi_contrib = (double *)xmalloc((size_t)opt_locus_count *
                                                   sizeof(double));
        x->notheta_old_phi_contrib = (double *)xcalloc((size_t)opt_locus_count,
                                                       sizeof(double));
        x->hybrid->notheta_phi_contrib = (double *)xmalloc((size_t)opt_locus_count *
                                                           sizeof(double));
        x->hybrid->notheta_old_phi_contrib = (double *)xcalloc((size_t)opt_locus_count,
                                                               sizeof(double));
        if (!LOAD(x->notheta_phi_contrib, opt_locus_count, fp))
          fatal("Cannot read per-locus phi contributions");
        if (!LOAD(x->hybrid->notheta_phi_contrib, opt_locus_count, fp))
          fatal("Cannot read per-locus phi contributions");
        if (!LOAD(&(x->hphi_sum),1,fp))
          fatal("Cannot read hphi sum");
        if (!LOAD(&(x->hybrid->hphi_sum),1,fp))
          fatal("Cannot read hphi sum");
      }
    }
  }

  /* read root_tau */
  if (!LOAD(&(stree->root_age),1,fp))
    fatal("Cannot read species root tau");

  /* read number of incoming sequences for each node node */
  for (i = 0; i < total_nodes; ++i)
    if (!LOAD(stree->nodes[i]->seqin_count,opt_locus_count,fp))
        fatal("Cannot read incoming sequence counts");

  alloc_gtree();
  unsigned int max_tips = 0;
  for (i = 0; i < opt_locus_count; ++i)
    if (gtree[i]->tip_count > max_tips)
      max_tips = gtree[i]->tip_count;

  unsigned int gtree_inner_sum = 0;
  for (i = 0; i < opt_locus_count; ++i)
    gtree_inner_sum += gtree[i]->inner_count;
  
  
  long * locus_seqcount = (long *)xmalloc((size_t)(opt_locus_count) * sizeof(long));
  for (i = 0; i < (unsigned int)opt_locus_count; ++i)
    locus_seqcount[i] = gtree[i]->tip_count;
  stree_alloc_internals(stree,locus_seqcount,gtree_inner_sum,opt_locus_count);
  free(locus_seqcount);

  /* read event indices for each node */
  unsigned int * buffer = (unsigned int *)xmalloc((size_t)(max_tips*2-1) *
                                                  sizeof(unsigned int));
  for (i = 0; i < total_nodes; ++i)
  {
    snode_t * snode = stree->nodes[i];

    for (j = 0; j < opt_locus_count; ++j)
    {
      snode->event[j] = dlist_create();

      if (!LOAD(buffer,snode->event_count[j],fp))
        fatal("Cannot read coalescent events");

      for (k = 0; k < snode->event_count[j]; ++k)
      {
        gnode_t * gt_node = gtree[j]->nodes[buffer[k]];

        dlist_item_t * dlitem = dlist_append(snode->event[j], (void *)gt_node);

        gt_node->event = dlitem;
      }
    }
  }
  free(buffer);

  if (opt_migration)
  {
    for (i = 0; i < total_nodes; ++i)
    {
      snode_t * x = stree->nodes[i];
      x->migevent_count = (long *)xcalloc((size_t)opt_locus_count,sizeof(long));
      
      if (!LOAD(x->migevent_count,opt_locus_count,fp))
        fatal("Cannot read migevent_count");

    }
    for (i = 0; i < total_nodes; ++i)
    {
      snode_t * x = stree->nodes[i];
      x->migbuffer = (migbuffer_t *)xcalloc((size_t)(stree->inner_count),
                                            sizeof(migbuffer_t));
      x->mig_source = (dlist_t **)xcalloc((size_t)opt_locus_count,
                                          sizeof(dlist_t *));
      x->mig_target = (dlist_t **)xcalloc((size_t)opt_locus_count,
                                          sizeof(dlist_t *));
      for (j = 0; j < opt_locus_count; ++j)
      {
        x->mig_source[j] = dlist_create();
        x->mig_target[j] = dlist_create();
      }
      if (!LOAD(&(x->mb_count), 1, fp))
        fatal("Cannot read migbuffer count");

      if (!LOAD(x->migbuffer, x->mb_count, fp))
        fatal("Cannot load node migbuffers");
    }
  }
}

static void load_gene_tree(FILE * fp, long index)
{
  long i,j;
  unsigned int gtree_tip_count = 0;

  gtree_t * gt = gtree[index];
  

  /* get number of tips */
  for (i = 0; i < stree->tip_count; ++i)
    gtree_tip_count += stree->nodes[i]->seqin_count[index]; 

  gt->tip_count = gtree_tip_count;
  gt->inner_count = gtree_tip_count-1;
  gt->edge_count = gt->tip_count + gt->inner_count - 1;

  /* read gene tree labels */
  char ** labels = (char **)xmalloc((size_t)gtree_tip_count*sizeof(char *));
  for (i = 0; i < gtree_tip_count; ++i)
  {
    if (!load_string(fp,labels+i))
      fatal("Cannot read gene tree %ld labels", index);
  }
  for (i = 0; i < gtree_tip_count; ++i)
  {
    gt->nodes[i]->label = labels[i];
  }

  free(labels);

  /* reset parent nodes to NULL */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    gt->nodes[i]->parent = NULL;

  /* reset data element */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    gt->nodes[i]->data = NULL;

  /* reset child nodes to NULL for tips */
  for (i = 0; i < gt->tip_count; ++i)
    gt->nodes[i]->left = gt->nodes[i]->right = NULL;

  /* set labels of inner nodes to NULL */
  for (i = 0; i < gt->inner_count; ++i)
    gt->nodes[gt->tip_count+i]->label = NULL;

  /* set left child for each inner node according to checkpoint data */
  for (i = 0; i < gt->inner_count; ++i)
  {
    unsigned int left_child_index;
    if (!LOAD(&left_child_index,1,fp))
      fatal("Cannot read species tree topology");
    gt->nodes[gt->tip_count+i]->left = gt->nodes[left_child_index];
    gt->nodes[left_child_index]->parent = gt->nodes[gt->tip_count+i];
  }

  /* set right child for each inner node according to checkpoint data */
  for (i = 0; i < gt->inner_count; ++i)
  {
    unsigned int right_child_index;
    if (!LOAD(&right_child_index,1,fp))
      fatal("Cannot read species tree topology");
    gt->nodes[gt->tip_count+i]->right = gt->nodes[right_child_index];
    gt->nodes[right_child_index]->parent = gt->nodes[gt->tip_count+i];
  }

  /* now set the root node */
  unsigned nullparent_count = 0;
  for (i = 0; i < gt->inner_count; ++i)
  {
    if (!gt->nodes[gt->tip_count+i]->parent)
    {
      gt->root = gt->nodes[gt->tip_count+i];
      nullparent_count++;
    }
  }
  if (nullparent_count != 1)
    fatal("Erroneous gene tree %ld structure", index);

  /* now check gene tree consistency */
  for (i = 0; i < gt->inner_count; ++i)
    assert((gt->nodes[i]->left == gt->nodes[i]->right) &&
           (gt->nodes[i]->left == NULL));

  for (i = 0; i < gt->inner_count; ++i)
  {
    unsigned int index = gt->tip_count + i;

    assert(gt->nodes[index]->left != gt->nodes[index]->right);
    if (gt->nodes[index]->parent)
      assert(gt->nodes[index]->left != gt->nodes[index]->parent);
    assert(gt->nodes[index]->left != gt->nodes[index]);

    assert(gt->nodes[index]->left->parent == gt->nodes[index]);
    assert(gt->nodes[index]->right->parent == gt->nodes[index]);
  }

  /* load branch lengths - TODO: Candidate for removal */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(&(gt->nodes[i]->length),1,fp))
      fatal("Cannot read gene tree branch lengths");

  /* load ages */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(&(gt->nodes[i]->time),1,fp))
      fatal("Cannot read gene tree node ages");

  /* load population index (corresponding species tree node index) */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
  {
    unsigned int pop_index;
    if (!LOAD(&(pop_index),1,fp))
      fatal("Cannot read gene tree population indices");
    gt->nodes[i]->pop = stree->nodes[pop_index];
  }

  /* load CLV indices */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(&(gt->nodes[i]->clv_index),1,fp))
      fatal("Cannot read gene tree clv indices");

  /* load scaler indices */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(&(gt->nodes[i]->scaler_index),1,fp))
      fatal("Cannot read gene tree scaler indices");

  /* load pmatrix indices */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(&(gt->nodes[i]->pmatrix_index),1,fp))
      fatal("Cannot read gene tree pmatrix indices");

  /* load mark - TODO: Candidate for removal */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(&(gt->nodes[i]->mark),1,fp))
      fatal("Cannot read gene tree marks");

  /* load hpath */
  for (i = 0; i < gt->tip_count + gt->inner_count; ++i)
    if (!LOAD(gt->nodes[i]->hpath,stree->hybrid_count,fp))
      fatal("Cannot read gene tree path flags");

  if (!LOAD(&(gt->rate_mui),1,fp))
      fatal("Cannot read gene tree mu_%ld", index);

  /* relaxed clock nu_i and logprior */
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    if (!LOAD(&(gt->rate_nui),1,fp))
      fatal("Cannot read gene tree nu_%ld", index);
    if (!LOAD(&(gt->lnprior_rates),1,fp))
      fatal("Cannot read gene tree %ld rate prior", index);
  }

  if (!LOAD(&(gt->original_index),1,fp))
    fatal("Cannot read gene tree original index");

  if (!LOAD(&(gt->msa_index),1,fp))
    fatal("Cannot read gene tree msa index");

  if (opt_migration)
  {
    /* mi structure */
    for (i = 0; i < gt->tip_count+gt->inner_count; ++i)
    {
      gnode_t * x = gt->nodes[i];
      long mi_count = 0;
      double mi_time = 0;
      unsigned int src_node_index;
      unsigned int tgt_node_index;
      

      x->mi = NULL;

      if (!LOAD(&mi_count,1,fp))
        fatal("Cannot load number of migration events on branch");

      if (!mi_count) continue;

      miginfo_check_and_extend(&(x->mi), mi_count);

      for (j = 0; j < mi_count; ++j)
      {
        if (!LOAD(&mi_time,1,fp))
          fatal("Cannot load migration event time");

        if (!LOAD(&src_node_index,1,fp))
          fatal("Cannot load migration event source index");

        if (!LOAD(&tgt_node_index,1,fp))
          fatal("Cannot load migration event target index");

        miginfo_append(&(x->mi),
                       stree->nodes[src_node_index],
                       stree->nodes[tgt_node_index],
                       mi_time,
                       gt->msa_index);
      }
      assert(x->mi->count == mi_count);
    }

    /*  migcount */
    long total_snodes = stree->tip_count+stree->inner_count;

    void * mem = xmalloc((size_t)(total_snodes*total_snodes)*sizeof(long) +
                         (size_t)total_snodes * sizeof(long *));
    gt->migcount = (long **)mem;
    gt->migcount[0] = (long *)(gt->migcount+total_snodes);
    memset(gt->migcount[0],0,total_snodes*sizeof(long));
    for (i = 1; i < total_snodes; ++i)
    {
      gt->migcount[i] = (long *)(gt->migcount[i-1] + total_snodes);
      memset(gt->migcount[i],0,total_snodes*sizeof(long));
    }
    for (i = 0; i < total_snodes; ++i)
    {
      if (!LOAD(gt->migcount[i],total_snodes,fp))
        fatal("Cannot load migcounts");
    }

    gt->migpops = (snode_t **)xcalloc((size_t)total_snodes, sizeof(snode_t *));
    gt->rb_linked = NULL;
    if (opt_exp_imrb)
      gt->rb_linked = (snode_t **)xmalloc((size_t)(total_snodes+1) *
                                          sizeof(snode_t *));
    gt->rb_lcount = 0;
  }
}

static void load_chk_section_3(FILE * fp, long msa_count)
{
  long i;

  for (i = 0; i < msa_count; ++i)
  {
    load_gene_tree(fp,i);
  }
}

static void load_locus(FILE * fp, long index)
{
  long i;

  /* load number of sites */
  unsigned int sites;
  unsigned int states;
  unsigned int rate_cats;
  unsigned int rate_matrices;
  unsigned int prob_matrices;
  unsigned int scale_buffers;
  unsigned int attributes;
  unsigned int dtype;
  unsigned int model;
  size_t span;

  gtree_t * gt = gtree[index];

  if (!LOAD(&dtype,1,fp))
    fatal("Cannot read data type");

  if (!LOAD(&model,1,fp))
    fatal("Cannot read substitution model");

  if (!LOAD(&sites,1,fp))
    fatal("Cannot read number of sites");

  /* load number of states */
  if (!LOAD(&states,1,fp))
    fatal("Cannot read number of states");

  /* load number of rate categories */
  if (!LOAD(&rate_cats,1,fp))
    fatal("Cannot read number of rate categories");

  /* load number of rate matrices */
  if (!LOAD(&rate_matrices,1,fp))
    fatal("Cannot read number of rate matrices");

  /* load number of prob matrices */
  if (!LOAD(&prob_matrices,1,fp))
    fatal("Cannot read number of rate matrices");

    /* load number of scale buffers */
  if (!LOAD(&scale_buffers,1,fp))
    fatal("Cannot read number of scale buffers");

  /* TODO: Store and load opt_scaling value instead */
  if (scale_buffers)
    opt_scaling = 1;

  /* load attributes */
  if (!LOAD(&attributes,1,fp))
    fatal("Cannot read attributes");

  locus[index] = locus_create(dtype,
                              model,
                              gt->tip_count,
                              2*gt->inner_count,
                              states,
                              sites,
                              rate_matrices,
                              prob_matrices,
                              rate_cats,
                              scale_buffers,
                              attributes);
  
  /* set frequencies for model with index 0 */
  //locus_set_frequencies_and_rates(locus[index]);

  for (i = 0; i < rate_matrices; ++i)
    locus[index]->eigen_decomp_valid[i] = 0;

  /* load pattern weights sum */
  if (!LOAD(&(locus[index]->pattern_weights_sum),1,fp))
    fatal("Cannot read pattern weights sum");
  
  /* load alpha */
  if (!LOAD(&(locus[index]->rates_alpha),1,fp))
    fatal("Cannot read alpha value");

  /* load qrates param count */
  if (!LOAD(&(locus[index]->qrates_param_count),1,fp))
    fatal("Cannot read qrates param count");

  /* load freqs param count */
  if (!LOAD(&(locus[index]->freqs_param_count),1,fp))
    fatal("Cannot read freqs param count");

  /* load category rates */
  if (!LOAD(locus[index]->rates,rate_cats,fp))
    fatal("Cannot read category rates");

  /* load base frequencies */
  for (i = 0; i < rate_matrices; ++i)
    if (!LOAD(locus[index]->frequencies[i],states,fp))
      fatal("Cannot read base frequencies");

  /* load qmatrix rates */
  for (i = 0; i < rate_matrices; ++i)
    if (!LOAD(locus[index]->subst_params[i],((states-1)*states)/2,fp))
      fatal("Cannot read qmatrix rates");

  /* load param indices */
  if (!LOAD(locus[index]->param_indices,locus[index]->rate_cats,fp))
    fatal("Cannot read param indices");

  /* load heredity scalars */
  if (!LOAD(locus[index]->heredity,locus[index]->rate_matrices,fp))
    fatal("Cannot read heredity scalars");

  /* load diploid */
  if (!LOAD(&(locus[index]->diploid),1,fp))
    fatal("Cannot read locus %ld diploid", index);
  
  /* TODO with more complex mixture models where rate_matrices > 1 we need
       to revisit this */
  assert(rate_matrices == 1);

  if (locus[index]->diploid)
  {
    size_t sites_a2 = 0;
    /* TODO: locus->pattern_weights is allocated in locis_create with a size
       equal to length of A3, but in reality we only need |A1| storage space.
       Free and reallocate here with 'unphased_length' */

    /* load original diploid number of sites */
    if (!LOAD(&(locus[index]->unphased_length),1,fp))
      fatal("Cannot read locus %ld unphased length", index);
    
    locus[index]->diploid_resolution_count = (unsigned long *)xmalloc((size_t)
                                             (locus[index]->unphased_length) *
                                             sizeof(unsigned long));
    locus[index]->likelihood_vector = (double *)xmalloc((size_t)
                                      (locus[index]->sites)*sizeof(double));

    /* load diploid resolution count */
    if (!LOAD(locus[index]->diploid_resolution_count,locus[index]->unphased_length,fp))
      fatal("Cannot read locus %ld diploid resolution count", index);

    /* load diploid mapping A1 -> A3 */
    for (i = 0; i < locus[index]->unphased_length; ++i)
      sites_a2 += locus[index]->diploid_resolution_count[i];
    locus[index]->diploid_mapping = (unsigned long *)xmalloc(sites_a2 *
                                              sizeof(unsigned long));
    if (!LOAD(locus[index]->diploid_mapping,sites_a2,fp))
      fatal("Cannot read locus %ld diploid mapping", index);


    /* load pattern weights for original diploid A1 alignment */
    free(locus[index]->pattern_weights);
    locus[index]->pattern_weights = (unsigned int *)xmalloc((size_t)
                                    (locus[index]->unphased_length) *
                                    sizeof(unsigned int));
    if (!LOAD(locus[index]->pattern_weights,locus[index]->unphased_length,fp))
      fatal("Cannot read pattern weights");
  }
  else
  {
    /* load pattern weights */
    if (!LOAD(locus[index]->pattern_weights,locus[index]->sites,fp))
      fatal("Cannot read pattern weights");
  }
    

  /* load tip CLVs */
  for (i = 0; i < gt->tip_count; ++i)
  {
    unsigned int clv_index = gt->nodes[i]->clv_index;
    span = locus[index]->sites * locus[index]->states * locus[index]->rate_cats;

    if (!LOAD(locus[index]->clv[clv_index],span,fp))
      fatal("Cannot read gene tree %ld tip CLV", index);
  }

  if (!LOAD(&(locus[index]->original_index),1,fp))
    fatal("Cannot read locus original index");
}

void load_chk_section_4(FILE * fp)
{
  long i;

  locus = (locus_t **)xmalloc((size_t)opt_locus_count * sizeof(locus_t));
  for (i = 0; i < opt_locus_count; ++i)
    load_locus(fp,i);

}

int checkpoint_load(gtree_t *** gtreep,
                    locus_t *** locusp,
                    stree_t ** streep,
                    double ** pjump,
                    unsigned long * curstep,
                    long * ft_round,
                    long * ndspecies,
                    long * mcmc_offset,
                    long * out_offset,
                    long ** gtree_offset,
                    long ** rates_offset,
                    long * dparam_count,
                    double ** posterior,
                    double ** pspecies,
                    long * ft_round_rj,
                    double * pjump_rj,
                    long * ft_round_spr,
                    long * ft_round_snl,
                    long * pjump_spr,
                    long * pjump_snl,
                    double * mean_logl,
                    long ** mean_mrate_row,
                    long ** mean_mrate_col,
                    long ** mean_mrate_round,
                    double ** mean_mrate,
                    double ** mean_tau,
                    double ** mean_theta,
                    double ** mean_phi,
                    long * mean_mrate_count,
                    long * mean_tau_count,
                    long * mean_theta_count,
                    long * mean_phi_count,
                    int * prec_logpg,
                    int * prec_logl)
{
  long i,j,k;
  FILE * fp;

  assert(opt_resume);

  fprintf(stdout, "\nResuming from checkpoint file %s\n\n", opt_resume);
  fp = fopen(opt_resume,"r");
  if (!fp)
    fatal("Cannot open checkpoint file %s", opt_resume);

  /* read header */
  #if 0
  fprintf(stdout,"HEADER:\n");
  #endif
  load_chk_header(fp);

  /* load section 1 */
  #if 0
  fprintf(stdout,"SECTION 1:\n");
  #endif

  load_chk_section_1(fp,
                     pjump,
                     curstep,
                     ft_round,
                     ndspecies,
                     mcmc_offset,
                     out_offset,
                     gtree_offset,
                     rates_offset,
                     dparam_count,
                     posterior,
                     pspecies,
                     ft_round_rj,
                     pjump_rj,
                     ft_round_spr,
                     ft_round_snl,
                     pjump_spr,
                     pjump_snl,
                     mean_logl,
                     mean_mrate_row,
                     mean_mrate_col,
                     mean_mrate_round,
                     mean_mrate,
                     mean_tau,
                     mean_theta,
                     mean_phi,
                     mean_mrate_count,
                     mean_tau_count,
                     mean_theta_count,
                     mean_phi_count,
                     prec_logpg,
                     prec_logl);

  /* load section 2 */
  load_chk_section_2(fp);

  /* initialize gene trees */
//  gtree = init_gtrees(opt_locus_count);

  /* load section 3 */
  load_chk_section_3(fp,opt_locus_count);

  /* load section 4 */
  load_chk_section_4(fp);

  /* TODO: set tip sequences, charmap etc when using tipchars */

  /* if migration then population migcount_sum */
  if (opt_migration)
  {
    unsigned int nodes_count = stree->tip_count + stree->inner_count;

    void * mem = xmalloc((size_t)(nodes_count*nodes_count)*sizeof(long) +
                         (size_t)nodes_count*sizeof(long *));
    stree->migcount_sum = (long **)mem;
    stree->migcount_sum[0] = (long *)(stree->migcount_sum+nodes_count);
    for (i = 1; i < nodes_count; ++i)
    {
      stree->migcount_sum[i] = (long *)(stree->migcount_sum[i-1] + nodes_count);
    }
    for (i = 0; i < nodes_count; ++i)
    {
      for (j = 0; j < nodes_count; ++j)
      {
        stree->migcount_sum[i][j] = 0;
        for (k = 0; k < opt_locus_count; ++k)
          stree->migcount_sum[i][j] = gtree[k]->migcount[i][j];
      }
    }
  }

  /* update pmatrices and CLVs */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_reset_leaves(gtree[i]->root);
    locus_update_all_matrices(locus[i],gtree[i],stree,i);
    locus_update_all_partials(locus[i],gtree[i]);

    gtree[i]->logl = locus_root_loglikelihood(locus[i],
                                              gtree[i]->root,
                                              locus[i]->param_indices,
                                              NULL);
  }

  #if 0
  logl = locus_root_loglikelihood(locus[i],
                                  gtree[i]->root,
                                  locus->param_indices,
                                  NULL);
  #endif
  fclose(fp);


  *streep = stree;
  *gtreep = gtree;
  *locusp = locus;

  #if 0
  printf("opt_est_locusrate: %ld\n", opt_est_locusrate);
  printf("opt_locusrate_prior: %ld\n", opt_locusrate_prior);
  #endif

  return 1;
}

void checkpoint_truncate(const char * filename, long offset)
{
  FILE * fp;
  
  if (!(fp = fopen(filename, "a")))
    fatal("Cannot open file %s for reading...", filename);

  if (xtruncate(fileno(fp),offset))
    fatal("Cannot truncate file %s to %ld bytes...", filename, offset);
  
  fclose(fp);
}
