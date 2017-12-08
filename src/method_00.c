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
  stree_t * stree;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing species tree...");

  assert(opt_streefile || opt_streenewick);
  assert(!(opt_streefile && opt_streenewick));

  if (opt_streefile)
    stree = stree_parse_newick(opt_streefile);
  else
    stree = stree_parse_newick_string(opt_streenewick);

  if (!stree)
    fatal("Error while reading species tree");

  return stree;
}

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

#if 0
static char * cb_serialize_none(const snode_t * snode)
{
  if (!snode->left) return xstrdup(snode->label);

  return xstrdup("");
}
#endif

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static long get_long(const char * line, long * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%ld%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static long get_double(const char * line, double * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%lf%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  return *x > *y;
}

static double eff_ict(double * x, long n, double mean, double stdev)
{
  /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
     sequence method */

  long i,j;
  double tint = 1;
  double rho, rho0 = 0;

  if (stdev/(fabs(mean)+1) < 1E-9)
  {
   tint = n;
  }
  else
  {
    for (i = 1; i < n-10; ++i)
    {
      rho = 0;
      for (j = 0; j < n - i; ++j)
        rho += x[j]*x[i+j];

      rho /= (n-1);

      if (i > 10 && rho+rho0 < 0)
        break;

      tint += rho*2;
      rho0 = rho;
    }
  }

  return 1/tint;
}

static void hpd_interval(double * x,
                         long n,
                         double * ltail,
                         double * rtail,
                         double alpha)
{
  long lrow = (long)(n*alpha/2);
  long urow = (long)(n*(1-alpha/2));
  long diffrow = urow - lrow;
  long l,r;

  long left = lrow;
  

  double w = x[urow] - x[lrow];

  *ltail = x[lrow];
  *rtail = x[urow];

  if (n <= 2) return;

  for (l=0,r=l+diffrow; r < n; l++,r++)
  {
    if (x[r] - x[l] < w)
    {
      left = l;
      w = x[r] - x[l];
    }
  }

  *ltail = x[left];
  *rtail = x[left + diffrow];
}

static void method_summary(stree_t * stree)
{
  long i,j,count;
  long sample_num;
  long rc = 0;
  FILE * fp;

  /* TODO: pretty-fy output */

  fp = xopen(opt_mcmcfile,"r");

  /* skip line containing header */
  getnextline(fp);
  assert(strlen(line) > 4);
  printf("          %s\n", line+4);

  /* compute number of columns in the file */
  long col_count = 0;

  /* compute number of theta parameters */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->theta >= 0)
      col_count++;

  /* compute number of tau parameters */
  for (i = 0; i < stree->inner_count; ++i)
    if (stree->nodes[stree->tip_count+i]->tau)
      col_count++;

  /* add one more for log-L */
  col_count++;

  /* allocate storage matrix */
  double ** matrix = (double **)xmalloc((size_t)col_count * sizeof(double *));
  for (i = 0; i < col_count; ++i)
    matrix[i] = (double *)xmalloc((size_t)opt_samples * sizeof(double));

  double * mean = (double *)xmalloc((size_t)col_count * sizeof(double));

  double * hpd025 = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * hpd975 = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * tint = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * stdev = (double *)xmalloc((size_t)col_count * sizeof(double));

  long line_count = 0;

  /* read data line by line and store in matrix */
  while (getnextline(fp))
  {
    double x;
    char * p = line;

    /* skip sample number */
    count = get_long(p,&sample_num);
    if (!count) goto l_unwind;

    p += count;

    /* read remaining elements of current row */

    for (i = 0; i < col_count; ++i)
    {
      count = get_double(p,&x);
      if (!count) goto l_unwind;

      p += count;

      matrix[i][line_count] = x;
    }

    line_count++;
  }

  /* compute means */
  printf("mean    ");
  for (i = 0; i < col_count; ++i)
  {
    double sum = 0;
    for (j = 0; j < opt_samples; ++j)
      sum += matrix[i][j];

    mean[i] = sum/opt_samples;
    printf("  %f", mean[i]);
  }
  printf("\n");

  /* compute standard deviation */
  for (i = 0; i < col_count; ++i)
  {
    double sd = 0;
    for (j = 0; j < opt_samples; ++j)
      sd += (matrix[i][j]-mean[i]) * (matrix[i][j]-mean[i]);

    stdev[i] = sqrt(sd/(opt_samples-1));
  }
  
  /* compute tint */
  for (i = 0; i < col_count; ++i)
    tint[i] = eff_ict(matrix[i],opt_samples,mean[i],stdev[i]);

  /* compute and print medians */
  printf("median  ");
  long median_line = opt_samples / 2;

  for (i = 0; i < col_count; ++i)
  {
    qsort(matrix[i], opt_samples, sizeof(double), cb_cmp_double);

    double median = matrix[i][median_line];
    if ((opt_samples & 1) == 0)
    {
      median += matrix[i][median_line-1];
      median /= 2;
    }

    printf("  %f", median);
  }
  printf("\n");

  /* print standard deviation */
  printf("S.D     ");
  for (i = 0; i < col_count; ++i)
  {
    printf("  %f", stdev[i]);
  }
  printf("\n");

  /* print minimum values */
  printf("min     ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][0]);
  printf("\n");

  /* print maximum values */
  printf("max     ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][opt_samples-1]);
  printf("\n");

  /* print line at 2.5% of matrix */
  printf("2.5%%    ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][(long)(opt_samples*.025)]);
  printf("\n");

  /* print line at 97.5% of matrix */
  printf("97.5%%   ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][(long)(opt_samples*.975)]);
  printf("\n");

  /* compute and print HPD 2.5% and 97.5% */
  for (i = 0; i < col_count; ++i)
    hpd_interval(matrix[i],opt_samples,hpd025+i,hpd975+i,0.05);

  /* print 2.5% HPD */
  printf("2.5%%HPD ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", hpd025[i]);
  printf("\n");

  /* print 97.5% HPD */
  printf("97.5%%HPD");
  for (i = 0; i < col_count; ++i)
    printf("  %f", hpd975[i]);
  printf("\n");

  /* print ESS */
  printf("ESS*    ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", opt_samples/tint[i]);
  printf("\n");
    
  /* print Eff */
  printf("Eff*    ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", 1/tint[i]);
  printf("\n");

  /* success */
  rc = 1;

l_unwind:
  for (i = 0; i < col_count; ++i)
    free(matrix[i]);
  free(matrix);

  free(mean);
  free(hpd025);
  free(hpd975);
  free(tint);
  free(stdev);

  fclose(fp);

  if (!rc)
    fatal("Error while reading/summarizing %s", opt_mcmcfile);

  return;
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
  unsigned int i;
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

void cmd_a00()
{
  int i,j;
  int msa_count;
  double logl,logpr;
  double logl_sum = 0;
  double logpr_sum = 0;
  msa_t ** msa_list;
  FILE * fp_mcmc;

  if (opt_samples < 1)
    fatal("--samples must be a positive integer greater than zero");

  if (opt_burnin < 0)
    fatal("--burnin must be a positive integer");

  /* load species tree */
  stree_t * stree = load_tree();
  printf(" Done\n");

  /* parse the phylip file */

  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  assert(fd);

  printf("Parsing phylip file...");
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);
  printf(" Done\n");

  phylip_close(fd);

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
    int ol = msa_list[i]->length;
    msa_list[i]->original_length = ol;
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

  /* initialize species tree (tau + theta) */
  stree_init(stree,msa_list,map_list,msa_count);
  stree_show_pptable(stree);

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
                            PLL_ATTRIB_ARCH_AVX);       /* attributes */

    /* set frequencies for model with index 0 */
    pll_set_frequencies(locus[i],0,frequencies);

    if (opt_diploid)
    {
      for (j = 0; j < (int)(stree->tip_count); ++j)
        if (stree->nodes[j]->diploid)
        {
          locus[i]->diploid = 1;
          break;
        }
    }

    /* set pattern weights and free the weights array */
    if (locus[i]->diploid)
    {
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

  /* start of MCMC loop */

  double * pjump = (double *)xcalloc(PROP_COUNT, sizeof(double));
  long ft_round = 0;

  /* print header in mcmc file */
  mcmc_printheader(fp_mcmc,stree);

  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);
  unsigned long curstep = 0;

  /* start of MCMC loop */
  for (i = -opt_burnin; i < opt_samples*opt_samplefreq; ++i)
  {
    
    /* update progress bar */
    if (!opt_quiet)
      progress_update(curstep);

    /* reset finetune parameters */
    if (i == 0 || (opt_finetune_reset && opt_burnin >= 200 && i < 0 &&
                   ft_round >= 100 && i%(opt_burnin/4)==0))
    {
      if (opt_finetune_reset && opt_burnin >= 200)
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
    if (i >= 0 && (i+1)%opt_samplefreq == 0)
    {
      mcmc_logsample(fp_mcmc,i+1,stree,gtree);
    }

    /* TODO: print on screen */
    
    curstep++;
  }

  progress_done();

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

  method_summary(stree);
  /* deallocate tree */
  stree_destroy(stree,NULL);
  stree_fini();

  /* deallocate alignments */
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  if (opt_diploid)
    free(opt_diploid);

  /* deallocate maplist */
  list_clear(map_list,map_dealloc);
  free(map_list);
}
