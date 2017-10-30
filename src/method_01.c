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
    fprintf(stdout, "Parsing tree file...\n");

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

static char * cb_serialize_none(const snode_t * snode)
{
  if (!snode->left) return xstrdup(snode->label);

  return xstrdup("");
}

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

static void strip_attributes(char * s)
{
  char * p = s;

  while (*s)
  {
    if (*s == '#' || *s == ':')
    {
      while (*s && *s != ',' && *s != ')' && *s != ';')
        ++s;
    }
    else if (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n')
      ++s;
    else
      *p++ = *s++;
  }
  *p = 0;
}

static int cb_strcmp(const void * a, const void * b)
{
  const char ** pa = (const char **)a;
  const char ** pb = (const char **)b;

  return strcmp(*pa,*pb);
}

static void stree_sort_recursive(snode_t * node)
{
  if (!node->left)
    return;

  stree_sort_recursive(node->left);
  stree_sort_recursive(node->right);

  if (strcmp(node->left->label, node->right->label) > 0)
    SWAP(node->left,node->right);

  node->label = (char *)xmalloc(strlen(node->left->label) +
                                strlen(node->right->label) + 1);

  /* concatenate species labels */
  node->label[0] = 0;
  strcat(node->label,node->left->label);
  strcat(node->label,node->right->label);
}

static void stree_sort(stree_t * stree)
{
  stree_sort_recursive(stree->root); 
}

struct freqtable_s
{
  long pos;
  long count;
};

static int cb_ftcmp(const void * a, const void * b)
{
  const struct freqtable_s * pa = (const struct freqtable_s *)a;
  const struct freqtable_s * pb = (const struct freqtable_s *)b;

  if (pa->count < pb->count) return 1;
  else if (pa->count > pb->count) return -1;

  return 0;
}

static void stree_summary(char ** species_names, long species_count)
{
  long i,line_count = 0;
  FILE * fp;
  char ** treelist;

  /* allocate space for holding all species tree samples */
  treelist = (char **)xmalloc((opt_samples+1)*sizeof(char *));

  /* open mcmc file */
  #ifndef DEBUG_MAJORITY
  fp = xopen(opt_mcmcfile,"r");
  #else
  fp = xopen("test.txt","r");
  #endif

  splits_init(4096,4096,species_names,species_count);

  /* read each line from the file, and strip all thetas and branch lengths
     such that only the tree topology and tip names remain, and store them
     in treelist */
  while (getnextline(fp))
  {
    strip_attributes(line);
    stree_t * t = stree_parse_newick_string(line);
    stree_sort(t);
    treelist[line_count++] = stree_export_newick(t->root,cb_serialize_none);

    splits_update(t);
    stree_destroy(t,NULL);
  }


  printf("Species in order:\n");
  for (i = 0; i < species_count; ++i)
  {
    printf(" %3ld. %s\n", i+1, species_names[i]);
  }
  printf("\n");

  qsort(treelist,line_count,sizeof(char *), cb_strcmp);

  long distinct = 1;
  long * uniquepos = (long *)xmalloc(line_count*sizeof(long));
  uniquepos[0] = 0;
  for (i = 1; i < line_count; ++i)
  {
    if (strcmp(treelist[i],treelist[i-1]))
    {
      uniquepos[distinct++] = i;
    }
  }

  struct freqtable_s * ft = (struct freqtable_s *)xmalloc(distinct*sizeof(struct freqtable_s));
  for (i = 0; i < distinct; ++i)
  {
    ft[i].pos   = uniquepos[i];
    ft[i].count = (i == distinct - 1) ?
                    line_count - uniquepos[i] : uniquepos[i+1] - uniquepos[i];
  }

  qsort(ft, distinct, sizeof(struct freqtable_s), cb_ftcmp);
  printf("(A) Best trees in the sample (%ld distinct trees in all)\n", distinct);
  double cdf = 0;
  for (i = 0; i < distinct; ++i)
  {
    double pdf = ft[i].count / (double)line_count;
    cdf += pdf;
    printf(" %8ld %8.5f %8.5f %s\n",
           ft[i].count, pdf, cdf, treelist[ft[i].pos]);
  }

  splits_finalize(line_count,species_names);

  fprintf(stdout, "\n(D) Best tree (or trees from the mastertree file) "
          "with support values\n");
  for (i = 0; i < distinct; ++i)
  {
    if (i && ft[i].count != ft[i-1].count)
      break;

    //printf("%s\n", treelist[ft[i].pos]);
    print_stree_with_support(treelist[ft[i].pos],ft[i].count,line_count);
  }

  summary_dealloc_hashtables();
  free(uniquepos);
  free(ft);
 
  /* deallocate list of trees */
  for (i = 0; i < line_count; ++i)
    free(treelist[i]);
  free(treelist);

  fclose(fp);
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

static void mcmc_printinitial(FILE * fp, stree_t * stree)
{
  char * newick = stree_export_newick(stree->root, cb_serialize_branch);
  fprintf(fp, "%s\n", newick);
  free(newick);
}

static void mcmc_logsample(FILE * fp,
                           int step,
                           stree_t * stree,
                           gtree_t ** gtree)
{
#if 0
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
  #endif
  char * newick = stree_export_newick(stree->root, cb_serialize_branch);
  fprintf(fp, "%s\n", newick);
  free(newick);
}

void cmd_a01()
{
  long i,j;
  int msa_count;
  long ft_round_spr = 0;
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

  /* parse the phylip file */
  printf("Parsed tree\n");

  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  assert(fd);

  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);

  phylip_close(fd);

  /* remove ambiguous sites */
  if (opt_cleandata)
  {
    for (i = 0; i < msa_count; ++i)
      if (!msa_remove_ambiguous(msa_list[i]))
        fatal("All sites in locus %d contain ambiguous characters",i);
  }

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
    printf("Locus %ld: original length %d, after compression %d\n", i, ol, msa_list[i]->length);
  }

  #if 0
  /* print the alignments */
  for (i = 0; i < msa_count; ++i)
    msa_print(msa_list[i]);
  #endif


  /* parse map file */
  printf("Parsing map file...\n");
  list_t * map_list = yy_parse_map(opt_mapfile);
  maplist_print(map_list);

  if (!(fp_mcmc = fopen(opt_mcmcfile, "w")))
    fatal("Cannot open file %s for writing...");

  /* initialize species tree (tau + theta) */
  stree_init(stree,msa_list,map_list,msa_count);

  gtree_t ** gtree = gtree_init(stree,msa_list,map_list,msa_count);

  /* the below two lines are specific to method 01 and they generate
     space for cloning the species and gene trees */
  
  stree_t * sclone = stree_clone_init(stree);
  gtree_t ** gclones = (gtree_t **)xmalloc((size_t)msa_count*sizeof(gtree_t *));
  for (i = 0; i < msa_count; ++i)
    gclones[i] = gtree_clone_init(gtree[i], sclone);


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
                            2*gtree[i]->edge_count,     /* # prob matrices */
                            1,                          /* # rate categories */
                            0,                          /* # scale buffers */
                            opt_arch);                  /* attributes */

    /* set frequencies for model with index 0 */
    pll_set_frequencies(locus[i],0,frequencies);

    /* set pattern weights and free the weights array */
    pll_set_pattern_weights(locus[i], weights[i]);
    free(weights[i]);


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

  printf("logL0 = %f   logP0 = %f\n", logl_sum, logpr_sum);

  /* free weights array */
  free(weights);

  /* start of MCMC loop */

  double * pjump = (double *)xcalloc(PROP_COUNT, sizeof(double));
  long ft_round = 0;
  long pjump_slider = 0;

  /* print header in mcmc file */
  mcmc_printinitial(fp_mcmc,stree);

  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);
  unsigned long curstep = 0;

  /* start of MCMC loop */
  #ifndef DEBUG_MAJORITY
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
      ft_round_spr = 0;
      pjump_slider = 0;
      memset(pjump,0,PROP_COUNT*sizeof(double));
    }

    ++ft_round;

    /* perform proposals sequentially */   
    double ratio;

    if (legacy_rndu() > 0)
    {
      long accepted = stree_propose_spr(&stree, &gtree, &sclone, &gclones, locus);
      if (accepted)
      {
        /* swap the pointers of species tree and gene tree list with cloned */
        SWAP(stree,sclone);
        SWAP(gtree,gclones);

        stree_label(stree);

        pjump_slider++;
      }

      ft_round_spr++;
    }

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
  #endif

  progress_done();

  free(pjump);

  fclose(fp_mcmc);

  for (i = 0; i < msa_count; ++i)
    locus_destroy(locus[i]);
  free(locus);

  /* deallocate gene trees */
  for (i = 0; i < msa_count; ++i)
  {
    gtree_destroy(gtree[i],NULL);
    gtree_destroy(gclones[i],NULL);
  }
  free(gtree);
  free(gclones);
  gtree_fini(msa_count);

  /* order of species */
  char ** species_names = (char **)xmalloc(stree->tip_count * sizeof(char *));
  unsigned int species_count = stree->tip_count;

  for (i = 0; i < stree->tip_count; ++i)
    species_names[i] = xstrdup(stree->nodes[i]->label);

  /* deallocate tree */
  stree_destroy(stree,NULL);
  stree_destroy(sclone,NULL);
  stree_fini();

  /* deallocate alignments */
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  /* deallocate maplist */
  list_clear(map_list,map_dealloc);
  free(map_list);

  stree_summary(species_names,species_count);

  for (i = 0; i < species_count; ++i)
    free(species_names[i]);
  free(species_names);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}
