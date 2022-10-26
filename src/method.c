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

//#define DEBUG_GTR

//#define CHECK_LOGL
//#define CHECK_LOGPR
//#define CHECK_LNPRIOR

/* maximum number of theta/tau to output on screen during MCMC */
#define MAX_MRATE_OUTPUT        4
#define MAX_THETA_OUTPUT        3
#define MAX_TAU_OUTPUT          3
#define MAX_PHI_OUTPUT          4

const static int rate_matrices = 1;
const static long thread_index_zero = 0;

static double pj_optimum = 0.3;
static thread_data_t td;
static time_t time_start;

static long enabled_prop_qrates = 0;
static long enabled_prop_freqs  = 0;
static long enabled_prop_alpha  = 0;

static int enabled_lrht  = 0;
static int enabled_hrdt  = 0;
static int enabled_mui   = 0;
static int enabled_mubar = 0;
static int enabled_nubar = 0;

static const char * template_ratesfile = "locus.%d.sample.txt";

static int prec_logl =  8;
static int prec_logpr = 8;

static int prec_ft = 6;

static void timer_start()
{
  time_start = time(NULL);
}

static void timer_print(const char * prefix, const char * suffix, FILE * fp)
{
  time_t t;
  long h,m,s;
  
  t = time(NULL) - time_start;

  h = (long)t / 3600;
  m = (long)(t % 3600) / 60;
  s = (long)(t - (t/60)*60);
  if (h)
  {
    fprintf(stdout, "%s%ld:%02ld:%02ld%s", prefix, h, m, s, suffix);
    fprintf(fp, "%s%ld:%02ld:%02ld%s", prefix, h, m, s, suffix);
  }
  else
  {
    fprintf(stdout, "%s%ld:%02ld%s", prefix, m, s, suffix);
    fprintf(fp, "%s%ld:%02ld%s", prefix, m, s, suffix);
  }

}

static void init_outfile(FILE * fp)
{
  struct tm * lt = NULL;
  char buffer[256];

  time_t t = time(NULL);  
  lt = localtime(&t);
  assert(lt);

  /* strftime(buffer, 256, "%a %b %d %T %Y", lt); */
  strftime(buffer, 256, "%c", lt);

  fprintf(fp, "Analysis started at: %s\n", buffer);
  fprintf(fp, "Using BPP version: %d.%d.%d\n",
          VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
  fprintf(fp, "Command: %s\n\n", cmdline);

  if (opt_seed > 0)
    fprintf(fp, "Seed: %ld (fixed by user)\n", opt_seed);
  else
  {
    unsigned int rseed = get_legacy_rndu_status(0);
    fprintf(fp, "Seed: %d (randomly generated)\n", rseed);
  }
}

static void print_mcmc_headerline(FILE * fp,
                                  stree_t * stree,
                                  gtree_t ** gtree,
                                  long * mean_mrate_row,
                                  long * mean_mrate_col,
                                  long mean_mrate_count)
{
  long i,j,k;
  long mean_theta_count = 0;
  long mean_tau_count = 0;
  long mean_phi_count = 0;

  long linewidth = 0;
  long ap_width = 0;
  long delim_count = 0;

  char * s = NULL;
  long * mean_theta_index = NULL;

  double logl = 0;
  double logpr = 0;
  if (opt_est_theta)
    for (i = 0; i < opt_locus_count; ++i)
      logpr += gtree[i]->logpr;
  else
    logpr = stree->notheta_logpr;

  for (i = 0; i < opt_locus_count; ++i)
    logl += gtree[i]->logl;

  /* compute digits for log-L */
  xasprintf(&s,"%8.5f", logl);
  long len_logl = strlen(s);
  free(s);

  xasprintf(&s,"%8.5f", logpr);
  long len_logpr = strlen(s);
  free(s);

  prec_logl  = (len_logl+1 > prec_logl) ? len_logl+1 : prec_logl;
  prec_logpr = (len_logpr+1 > prec_logpr) ? len_logpr+1 : prec_logpr;

  if (opt_method != METHOD_00)    /* species tree inference or delimitation */
  {
    mean_theta_count = 1;
    mean_tau_count = 1;
    mean_theta_index = (long *)xcalloc(1,sizeof(long));
    mean_theta_index[0] = stree->root->node_index;

    assert(!opt_msci);
  }
  else
  {
    /* compute mean thetas */

    /* 1a. calculate number of thetas to print */
    long max_param_count = MIN(stree->tip_count+stree->inner_count,
                               MAX_THETA_OUTPUT);
    mean_theta_index = (long *)xcalloc((size_t)max_param_count,sizeof(long));

    /* 1b. calculate number of mean thetas */
    k = 0;
    if (opt_est_theta)
    {
      for (j=0; j < stree->tip_count+stree->inner_count; ++j)
      {
        if (stree->nodes[j]->theta < 0) continue;

        mean_theta_index[k] = j;
        if (++k == max_param_count) break;
      }
    }
    mean_theta_count = k;

    /* compute mean taus */

    /* 2a. calculate number of taus to print */
    max_param_count = MIN(stree->tip_count+stree->inner_count,
                          MAX_TAU_OUTPUT);

    /* 2b. calculate means */
    k = 0;
    for (j = stree->tip_count; j < stree->tip_count+stree->inner_count; ++j)
    {
      if (stree->nodes[j]->tau == 0) continue;
      
      if (++k == max_param_count) break;
    }
    mean_tau_count = k;

    /* compute mean phis */
    mean_phi_count = MIN(stree->hybrid_count, MAX_PHI_OUTPUT);
  }

  /* print legend */
  fprintf(fp, "\n");
  fprintf(fp, "-*- Terms index -*-\n\n");
  fprintf(fp, "  Prgs: progress of MCMC run (negative progress means burnin)\n");
  fprintf(fp, "  Gage: gene-tree age proposal\n");
  fprintf(fp, "  Gspr: gene-tree SPR proposal\n");
  fprintf(fp, "  thet: species tree theta proposal\n");
  fprintf(fp, "   tau: species tree tau proposal\n");
  fprintf(fp, "   mix: mixing proposal\n");
  if (opt_migration)
    fprintf(fp, "  mrte: migration rates proposal\n");
  if (enabled_hrdt)
    fprintf(fp, "  hrdt: heredity proposal\n");
  if (enabled_lrht)
    fprintf(fp, "  lrht: locus rate and heredity proposal\n");
  if (enabled_mui)
    fprintf(fp, "  mu_i: locus rate (mu_i) proposal\n");
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    fprintf(fp, "  nu_i: locus rate variance (nu_i) proposal\n");
    fprintf(fp, "  brte: per-locus species tree Branch rate proposal\n");
  }
  if (enabled_mubar)
    fprintf(fp, "  mubr: average locus rate (mu_bar) proposal\n");
  if (enabled_nubar)
    fprintf(fp, "  nubr: average locus rate variance (nu_bar) proposal\n");
  if (opt_msci)
    fprintf(fp, "   phi: MSCi phi parameter proposal\n");
  if (enabled_prop_freqs)
    fprintf(fp, "    pi: base frequencies proposal\n");
  if (enabled_prop_qrates)
    fprintf(fp, "  qmat: instantaneous substitution rates proposal\n");
  if (enabled_prop_alpha)
    fprintf(fp, "  alfa: discretized gamma (rate variation among sites) alpha proposal\n");
  if (opt_method == METHOD_01)
    fprintf(fp, " stree: Sspr (%6.4f) & Ssnl (%6.4f)\n", 1 - opt_prob_snl, opt_prob_snl);
  if (opt_method == METHOD_10)
    fprintf(fp, "    rj: reversible-jump split/merge proposal\n");
  if (opt_method == METHOD_11)
  {
    fprintf(fp, "  stree proposal:  %6.4f for Ssnl, %6.4f for Sspr\n", opt_prob_snl, 1 - opt_prob_snl);
    fprintf(fp, "    sp: number of delimited species\n");
  }
  if (opt_method == METHOD_10 || opt_method == METHOD_11)
  {
    fprintf(fp, "    np: number of parameters\n");
  }
  if (opt_method == METHOD_10)
  {
    fprintf(fp, "   del: delimitation model\n");
  }
  if (opt_method == METHOD_10 || opt_method == METHOD_11)
  {
    fprintf(fp, "  mldp: most likely delimitation and probability\n");
  }
  if (mean_theta_count == 1)
    fprintf(fp, "mthet1: root node mean theta\n");
  else
  {
    for (i = 0; i < mean_theta_count; ++i)
      fprintf(fp, "mthet%ld: mean theta of node %ld\n", i+1,mean_theta_index[i]);
  }
  if (mean_tau_count == 1)
    fprintf(fp, " mtau1: root node mean tau\n");
  else
  {
    for (i = 0; i < mean_tau_count; ++i)
      fprintf(fp," mtau%ld: mean tau of node %ld\n", i+1,stree->tip_count+i);
  }
  if (opt_migration)
  {
    for (i = 0; i < mean_mrate_count; ++i)
    {
      long j = mean_mrate_row[i];
      long k = mean_mrate_col[i];

      fprintf(fp,
              " mmrt%ld: mean migration rate %s -> %s\n",
              i+1,
              stree->nodes[j]->label,
              stree->nodes[k]->label);
    }
  }
  if (opt_msci)
  {
    for (i = 0; i < mean_phi_count; ++i)
    {
      snode_t * tmpnode = stree->nodes[stree->tip_count+stree->inner_count+i];

      #if 0
      /* old code before introduction of has_phi */
      if (!node_is_bidirection(tmpnode))
      {
        /* hybridization node */

        /* if main node htau==0 and mirror node htau==1 then that is the only case
           we use the phi from the main node */
        if (tmpnode->hybrid->htau == 0 && tmpnode->htau == 1)
          tmpnode = tmpnode->hybrid;
      }
      #else
      /* new correct code */
      if (!tmpnode->has_phi)
        tmpnode = tmpnode->hybrid;

      #endif

      fprintf(fp,
              " mphi%ld: mean of phi_%s : %s -> %s\n",
              i+1,
              tmpnode->label,
              tmpnode->parent->label,
              tmpnode->label);
    }
  }
  if (opt_msci)
    fprintf(fp, "log-PG: log-probability of gene trees (MSCi)\n");
  else
    fprintf(fp, "log-PG: log-probability of gene trees (MSC)\n");
  fprintf(fp, " log-L: mean log-L of observing data\n");
  fprintf(fp,"\n");

  ap_width += 5*5;
  ap_width += enabled_hrdt ? 5 : 0;
  ap_width += enabled_lrht ? 5 : 0;
  ap_width += enabled_mui ? 5 : 0;
  ap_width += opt_clock != BPP_CLOCK_GLOBAL ? 10 : 0;
  ap_width += enabled_mubar ? 5 : 0;
  ap_width += enabled_nubar ? 5 : 0;
  ap_width += opt_msci ? 5 : 0;
  ap_width += enabled_prop_freqs ? 5 : 0;
  ap_width += enabled_prop_qrates ? 5 : 0;
  ap_width += enabled_prop_alpha ? 5 : 0;
  ap_width += opt_method == METHOD_01 ? 7*(1+!!opt_prob_snl) : 0;
  ap_width += opt_method == METHOD_10 ? 7 : 0;
  ap_width += opt_method == METHOD_11 ? 7*(2+!!opt_prob_snl) : 0;
  ap_width += opt_migration ? 5 : 0;
  ap_width += 1;
  
  if (opt_method == METHOD_10)
  {
    delim_count = delimitation_getparam_count();
  }
  else if (opt_method == METHOD_11)
    delim_count = delimitations_count(stree);

  char * ap_title = xstrdup("Acceptance proportions");
  long titlelen = strlen(ap_title);
  long prefix = (ap_width - titlelen)/2;
  long suffix = ap_width - titlelen  - prefix;
  fprintf(fp,"     |");
  for (i = 0; i < prefix; ++i)
    fprintf(fp," ");
  fprintf(fp,"%s",ap_title);
  for (i = 0; i < suffix; ++i)
    fprintf(fp," ");
  fprintf(fp,"|\n");
  free(ap_title);

  fprintf(fp,"Prgs |");     linewidth += 6;
  fprintf(fp," Gage");      linewidth += 5;
  fprintf(fp," Gspr");      linewidth += 5;
  fprintf(fp," thet");      linewidth += 5;
  fprintf(fp,"  tau");      linewidth += 5;
  fprintf(fp,"  mix");      linewidth += 5;
  if (enabled_hrdt)
  {
    fprintf(fp," hrdt");    linewidth += 5;
  }
  if (enabled_lrht)
  {
    fprintf(fp," lrht");    linewidth += 5;
  }
  if (enabled_mui)
  {
    fprintf(fp," mu_i");    linewidth += 5;
  }
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    fprintf(fp," nu_i");    linewidth += 5;
    fprintf(fp," brte");    linewidth += 5;
  }
  if (opt_migration)
  {
    fprintf(fp," mrte");    linewidth += 5;
  }
  if (enabled_mubar)
  {
    fprintf(fp," mubr");    linewidth += 5;
  }
  if (enabled_nubar)
  {
    fprintf(fp," nubr");    linewidth += 5;
  }
  if (opt_msci)
  {
    fprintf(fp,"  phi");    linewidth += 5;
  }
  if (enabled_prop_freqs)
  {
    fprintf(fp,"   pi");    linewidth += 5;
  }
  if (enabled_prop_qrates)
  {
    fprintf(fp," qmat");    linewidth += 5;
  }
  if (enabled_prop_alpha)
  {
    fprintf(fp," alfa");    linewidth += 5;
  }

  if (opt_method == METHOD_01)
  {
    fprintf(fp, "   Sspr");  linewidth += 7;
    if (opt_prob_snl)
    {
      fprintf(fp, "   Ssnl");  linewidth += 6;
    }
  }
  else if (opt_method == METHOD_10)
  {
    fprintf(fp,"     rj");  linewidth += 7;
  }
  else if (opt_method == METHOD_11)
  {
    fprintf(fp, "   Sspr");  linewidth += 7;
    if (opt_prob_snl)
    {
      fprintf(fp, "   Ssnl");  linewidth += 6;
    }
    fprintf(fp, "     rj");  linewidth += 7;
  }
  fprintf(fp," |");         linewidth += 2;

    /* TODO */
  if (opt_method == METHOD_10)
  {
    fprintf(fp," np");   linewidth += 3;
    fprintf(fp," %*s", stree->inner_count > 3 ? stree->inner_count : 3, "del");
    linewidth += stree->inner_count > 4 ? stree->inner_count+1 : 4;

    /* TODO */
    int digits = delim_count ? (int)floor(log10(abs(delim_count))) + 1 : 1;
    fprintf(fp," %*s", 4+6+digits, "mldp");     linewidth += 1+10+digits;

  }
  else if (opt_method == METHOD_11)
  {
    int digits = delim_count ? (int)floor(log10(abs(delim_count))) + 1 : 1;
    fprintf(fp," sp np %*s",4+6+digits,"mldp"); linewidth += 6+1+10+digits;

    /* TODO */
  }

  if (opt_est_theta)
  {
    for (i = 0; i < mean_theta_count; ++i)
    {
      fprintf(fp," mthet%ld", i+1);         linewidth += 7;
    }
    fprintf(fp," ");          linewidth += 1;
  }
  for (i = 0; i < mean_tau_count; ++i)
  {
    fprintf(fp,"  mtau%ld", i+1);           linewidth += 7;
  }

  if (opt_migration)
  {
    fprintf(fp, " ");                       linewidth += 1;
    for (i = 0; i < mean_mrate_count; ++i)
    {
      fprintf(fp," mmrte%ld", i+1);         linewidth += 7;
    }
  }

  if (opt_msci)
  {
    fprintf(fp, " ");                       linewidth += 1;
    for (i = 0; i < mean_phi_count; ++i)
    {
      fprintf(fp,"  mphi%ld",i+1);          linewidth += 7;
    }
  }

  fprintf(fp," ");          linewidth += 1;
  fprintf(fp," %*s", prec_logpr,"log-PG");  linewidth += prec_logpr+1;
  if (opt_usedata)
  {
    fprintf(fp," %*s", prec_logl,"log-L");  linewidth += prec_logl+1;
  }
  fprintf(fp,"\n");

  /* TODO */
  for (i = 0; i < linewidth; ++i)
    fprintf(fp,"-");
  fprintf(fp,"\n");
  
  if (mean_theta_index)
    free(mean_theta_index);
}

static stree_t * load_tree_or_network(void)
{
  stree_t * stree;

  if (opt_seed > 0)
    printf("Seed: %ld (fixed by user)\n", opt_seed);
  else
  {
    unsigned int rseed = get_legacy_rndu_status(0);
    printf("Seed: %d (randomly generated)\n", rseed);
  }

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing species tree...");

  assert(opt_streenewick);

  stree = bpp_parse_newick_string(opt_streenewick);
  if (!stree)
    fatal("Error while reading species tree");


  return stree;
}

static void reset_finetune_onestep(double pjump, double * param)
{
  double maxstep = 99;

  if (pjump < 0.001)
    *param /= 100;
  else if (pjump > 0.999)
    *param = MIN(maxstep, *param * 100);
  else
  {
    *param *= tan(BPP_PI/2*(pjump)) / tan(BPP_PI/2*pj_optimum);
    *param = MIN(maxstep, *param);
  }
  
}

static void reset_finetune(FILE * fp_out, double * pjump)
{
  int i,j;
  int extra;
  int spacing = 1;
  int empty = 4;  /* for entries not available */
  FILE * fp[2];

  fp[0] = stdout; fp[1] = fp_out;
  
  /* TODO: Adjust such that each move has a different pjump and steplength */
  extra = (opt_est_heredity == HEREDITY_ESTIMATE ||
           (opt_est_locusrate == MUTRATE_ESTIMATE &&
            opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR));

  for (j = 0; j < 2; ++j)
  {
    fprintf(fp[j], "\n                   ");
    fprintf(fp[j],  "%*s", prec_ft+spacing, "Gage");      /*  0 */
    fprintf(fp[j], " %*s", prec_ft+spacing, "Gspr");      /*  1 */
    fprintf(fp[j], " %*s", prec_ft+spacing, "thet");      /*  2 */
    fprintf(fp[j], " %*s", prec_ft+spacing, "tau");       /*  3 */
    fprintf(fp[j], " %*s", prec_ft+spacing, "mix");       /*  4 */
    if (extra)
      fprintf(fp[j], " %*s", prec_ft+spacing, "lrht");    /*  5 */
    else
      fprintf(fp[j], " %*s", empty, "lrht");
    if (opt_msci)
      fprintf(fp[j], " %*s", prec_ft+spacing, "phi");     /*  6 */
    else
      fprintf(fp[j], " %*s", empty, "phi");
    if (enabled_prop_freqs)
      fprintf(fp[j], " %*s", prec_ft+spacing, "pi");      /*  7 */
    else
      fprintf(fp[j], " %*s", empty, "pi");
    if (enabled_prop_qrates)
      fprintf(fp[j], " %*s", prec_ft+spacing, "qmat");    /*  8 */
    else
      fprintf(fp[j], " %*s", empty, "qmat");
    if (enabled_prop_alpha)
      fprintf(fp[j], " %*s", prec_ft+spacing, "alfa");    /*  9 */
    else
      fprintf(fp[j], " %*s", empty, "alfa");
    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
        opt_est_mubar)
      fprintf(fp[j], " %*s", prec_ft+spacing, "mubr");    /* 10 */
    else
      fprintf(fp[j], " %*s", empty, "mubr");
    if (opt_clock != BPP_CLOCK_GLOBAL &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      fprintf(fp[j], " %*s", prec_ft+spacing, "nubr");    /* 11 */
    else
      fprintf(fp[j], " %*s", empty, "nubr");
    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
         opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
      fprintf(fp[j], " %*s", prec_ft+spacing, "mu_i");    /* 12 */
    else
      fprintf(fp[j], " %*s", empty, "mu_i");
    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      fprintf(fp[j], " %*s", prec_ft+spacing, "nu_i");    /* 13 */
      fprintf(fp[j], " %*s", prec_ft+spacing, "brte");    /* 14 */
    }
    else
    {
      fprintf(fp[j], " %*s", empty, "nu_i");
      fprintf(fp[j], " %*s", empty, "brte");
    }
    if (opt_migration)
      fprintf(fp[j], " %*s", prec_ft+spacing, "mrte");    /* 15 */
    else
      fprintf(fp[j], " %*s", empty, "mrte");
  }
  for (j = 0; j < 2; ++j)
  {
    fprintf(fp[j], "\nCurrent Pjump:    ");

    for (i = 0; i < PROP_COUNT; ++i)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[i]);

    /* mu_i with GammaDir or Heredity scalars */
    if (extra)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_LRHT_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* phi pjump */
    if (opt_msci)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_PHI_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    if (enabled_prop_freqs)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_FREQS_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_qrates)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_QRATES_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_alpha)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_ALPHA_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* mubar pjump */
    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
        opt_est_mubar)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_MUBAR_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    
    /* nubar pjump */
    if (opt_clock != BPP_CLOCK_GLOBAL &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_NUBAR_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");


    /* mu_i with conditional iid */
    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
         opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
      fprintf(fp[j], " %*.5f", prec_ft+spacing, pjump[BPP_MOVE_MUI_INDEX]);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* nu_i and branch rates pjump */
    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      fprintf(fp[j], " %*.5f",prec_ft+spacing,pjump[BPP_MOVE_NUI_INDEX]);
      fprintf(fp[j], " %*.5f",prec_ft+spacing,pjump[BPP_MOVE_BRANCHRATE_INDEX]);
    }
    else
    {
      fprintf(fp[j], " %*s", empty, "- ");
      fprintf(fp[j], " %*s", empty, "- ");
    }

    if (opt_migration)
    {
      fprintf(fp[j], " %*.5f",prec_ft+spacing,pjump[BPP_MOVE_MRATE_INDEX]);
    }
    else
    {
      fprintf(fp[j], " %*s", empty, "- ");
    }


    fprintf(fp[j], "\n");

    //printf("finetune_mui: %.15f", opt_finetune_mui);
    fprintf(fp[j], "Current finetune:");
//    fprintf(fp[j], " %8.5f", opt_finetune_gtage);
    fprintf(fp[j], " %*.5f", prec_ft+spacing+1, opt_finetune_gtage);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_gtspr);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_theta);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_tau);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_mix);

    /* mu_i with GammaDir or Heredity scalars */
    if (extra)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_locusrate);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* phi */
    if (opt_msci)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_phi);
    else
      fprintf(fp[j], " %*s", empty, "- ");


    if (enabled_prop_freqs)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_freqs);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_qrates)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_qrates);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_alpha)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_alpha);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* mubar */
    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
        opt_est_mubar)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_mubar);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* nubar */
    if (opt_clock != BPP_CLOCK_GLOBAL &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_nubar);
    else
      fprintf(fp[j], " %*s", empty, "- ");


    /* mu_i with conditional iid */
    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
         opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_mui);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    /* nu_i and branch rates */
    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_nui);
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_branchrate);
    }
    else
    {
      fprintf(fp[j], " %*s", empty, "- ");
      fprintf(fp[j], " %*s", empty, "- ");
    }

    if (opt_migration)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_migrates);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    fprintf(fp[j], "\n");
  }


  reset_finetune_onestep(pjump[BPP_MOVE_GTAGE_INDEX],&opt_finetune_gtage);
  reset_finetune_onestep(pjump[BPP_MOVE_GTSPR_INDEX],&opt_finetune_gtspr);
  reset_finetune_onestep(pjump[BPP_MOVE_THETA_INDEX],&opt_finetune_theta);
  reset_finetune_onestep(pjump[BPP_MOVE_TAU_INDEX],&opt_finetune_tau);
  reset_finetune_onestep(pjump[BPP_MOVE_MIX_INDEX],&opt_finetune_mix);

  if (extra)
    reset_finetune_onestep(pjump[BPP_MOVE_LRHT_INDEX],&opt_finetune_locusrate);
  if (opt_msci)
    reset_finetune_onestep(pjump[BPP_MOVE_PHI_INDEX], &opt_finetune_phi);
  if (enabled_prop_freqs)
    reset_finetune_onestep(pjump[BPP_MOVE_FREQS_INDEX], &opt_finetune_freqs);
  if (enabled_prop_qrates)
    reset_finetune_onestep(pjump[BPP_MOVE_QRATES_INDEX], &opt_finetune_qrates);
  if (enabled_prop_alpha)
    reset_finetune_onestep(pjump[BPP_MOVE_ALPHA_INDEX], &opt_finetune_alpha);
  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
      opt_est_mubar)
    reset_finetune_onestep(pjump[BPP_MOVE_MUBAR_INDEX], &opt_finetune_mubar);
  if (opt_clock != BPP_CLOCK_GLOBAL &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    reset_finetune_onestep(pjump[BPP_MOVE_NUBAR_INDEX], &opt_finetune_nubar);
  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
       opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
    reset_finetune_onestep(pjump[BPP_MOVE_MUI_INDEX], &opt_finetune_mui);
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    reset_finetune_onestep(pjump[BPP_MOVE_NUI_INDEX], &opt_finetune_nui);
    reset_finetune_onestep(pjump[BPP_MOVE_BRANCHRATE_INDEX], &opt_finetune_branchrate);
  }
  if (opt_migration)
    reset_finetune_onestep(pjump[BPP_MOVE_MRATE_INDEX], &opt_finetune_migrates);

  for (j = 0; j < 2; ++j)
  {
    fprintf(fp[j], "New finetune:    ");
    fprintf(fp[j], " %*.5f", prec_ft+spacing+1, opt_finetune_gtage);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_gtspr);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_theta);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_tau);
    fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_mix);

    if (extra)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_locusrate);
    else
      fprintf(fp[j], " %*s", empty,  "- ");

    if (opt_msci)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_phi);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_freqs)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_freqs);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_qrates)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_qrates);
    else
      fprintf(fp[j], " %*s", empty, "- ");
    if (enabled_prop_alpha)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_alpha);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
        opt_est_mubar)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_mubar);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    if (opt_clock != BPP_CLOCK_GLOBAL &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_nubar);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
         opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_mui);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_nui);
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_branchrate);
    }
    else
    {
      fprintf(fp[j], " %*s", empty, "- ");
      fprintf(fp[j], " %*s", empty, "- ");
    }
    if (opt_migration)
      fprintf(fp[j], " %*.5f", prec_ft+spacing, opt_finetune_migrates);
    else
      fprintf(fp[j], " %*s", empty, "- ");

    fprintf(fp[j], "\n");
  }
  fprintf(fp[0],"\n");
  fprintf(fp_out, "\n");
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

static void status_print_pjump(FILE * fp,
                               double * pjump,
                               long ft_round_spr,
                               long ft_round_snl,
                               long pjump_spr,
                               long pjump_snl,
                               double mean_pjump_rj)
{
  long j;
  int extra = (opt_est_heredity == HEREDITY_ESTIMATE ||
               (opt_est_locusrate == MUTRATE_ESTIMATE &&
                opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR));

  for (j = 0; j < PROP_COUNT; ++j)
    fprintf(fp, " %4.2f", pjump[j]);
  if (extra)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_LRHT_INDEX]);
  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
       opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_MUI_INDEX]);
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_NUI_INDEX]);
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_BRANCHRATE_INDEX]);
  }
  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
      opt_est_mubar)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_MUBAR_INDEX]);
  if (opt_clock != BPP_CLOCK_GLOBAL &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_NUBAR_INDEX]);

  if (opt_msci)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_PHI_INDEX]);
  if (enabled_prop_freqs)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_FREQS_INDEX]);
  if (enabled_prop_qrates)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_QRATES_INDEX]);
  if (enabled_prop_alpha)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_ALPHA_INDEX]);
  if (opt_migration)
    fprintf(fp, " %4.2f", pjump[BPP_MOVE_MRATE_INDEX]);

  /* print pjump for species tree SPR */
  if (opt_method == METHOD_01)
  {
    fprintf(fp, " %5.4f", ft_round_spr ? (double)pjump_spr / ft_round_spr : 0.);
    if (opt_prob_snl)
      fprintf(fp, " %5.4f", ft_round_snl ? (double)pjump_snl / ft_round_snl : 0.);
  }
  else if (opt_method == METHOD_10)
    fprintf(fp," %5.4f", mean_pjump_rj);
  else if (opt_method == METHOD_11)
  {
    fprintf(fp, " %5.4f", ft_round_spr ? (double)pjump_spr / ft_round_spr : 0.);
    if (opt_prob_snl)
      fprintf(fp, " %5.4f", ft_round_snl ? (double)pjump_snl / ft_round_snl : 0.);
    fprintf(fp," %5.4f", mean_pjump_rj);
  }
  fprintf(fp, "  ");
}

static void mcmc_printheader(FILE * fp, stree_t * stree)
{
  int print_labels = 1;
  unsigned int i,j;
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
    {
      #if 0
      /* old code before the introduction of has_phi */

      if (node_is_bidirection(stree->nodes[offset+i]))
        fprintf(fp, "\tphi_%s", stree->nodes[offset+i]->label);
      else
      {
        /* hybridization node */

        /* if main node htau==0 and mirror node htau==1 then that is the only
           case we use the phi from the main node */
        snode_t * tmpnode = stree->nodes[offset+i];
        if (tmpnode->hybrid->htau == 0 && tmpnode->htau == 1)
          tmpnode = tmpnode->hybrid;

        fprintf(fp,
                "\tphi_%s<-%s",
                tmpnode->label,
                tmpnode->parent->label);
      }
      #else

      /* new correct code */

      snode_t * tmpnode = stree->nodes[offset+i];
      if (!tmpnode->has_phi)
        tmpnode = tmpnode->hybrid;
      fprintf(fp,
              "\tphi_%s<-%s",
              tmpnode->label,
              tmpnode->parent->label);
      #endif

    }
  }

  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL &&
      opt_est_mubar)
    fprintf(fp, "\tmu_bar");
    
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      fprintf(fp, "\tnu_bar");
    else
      fprintf(fp, "\tnu");
  }

  if (opt_migration)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        if (opt_mig_bitmatrix[i][j])
          fprintf(fp,
                  "\tM_%s->%s",
                  stree->nodes[i]->label,
                  stree->nodes[j]->label);
  }

  /* 5. Print log likelihood */
  if (opt_usedata)
    fprintf(fp, "\tlnL\n"); 
  else
    fprintf(fp, "\n");

}

static void mcmc_printheader_rates(FILE ** fp_locus,
                                   stree_t * stree,
                                   locus_t ** locus)
{
  int tab_required = 0;
  long i,j;
  unsigned int total_nodes;

  /* labels in TCAG order */
  char * flabels[4] = {"pi_T", "pi_C", "pi_A", "pi_G"};
  char * rlabels[6] = {"a(TC)", "b(TA)", "c(TG)", "d(CA)", "e(CG)", "f(AG)"};

  assert(!opt_msci || (opt_msci && !opt_est_stree));
  assert(fp_locus);

  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  for (i = 0; i < opt_locus_count; ++i)
  {
    tab_required = 0;
    /* print heredity scalars header */
    if (opt_est_heredity == HEREDITY_ESTIMATE && opt_print_hscalars)
    {
      fprintf(fp_locus[i],
              "%sheredity_L%d",
              tab_required ? "\t" : "", locus[i]->original_index+1);
      tab_required = 1;
    }

    /* print mu_i header */
    if (opt_est_locusrate == MUTRATE_ESTIMATE && opt_print_locusrate)
    {
      fprintf(fp_locus[i],
              "%smu_%d",
              tab_required ? "\t" : "", locus[i]->original_index+1);
      tab_required = 1;
    }
    /* print nu_i header */
    if (opt_clock != BPP_CLOCK_GLOBAL && opt_print_rates)
    {
      fprintf(fp_locus[i],
              "%snu_%d",
              tab_required ? "\t" : "", locus[i]->original_index+1);
      tab_required = 1;
    }
    /* print species tree branch rates header */
    if (opt_clock != BPP_CLOCK_GLOBAL && opt_print_rates)
    {
      fprintf(fp_locus[i],
              "%sr_%s",
              tab_required ? "\t" : "", stree->nodes[0]->label);
      tab_required = 1;
      for (j = 1; j < total_nodes; ++j)
        if (stree->nodes[j]->brate)
          fprintf(fp_locus[i], "\tr_%s", stree->nodes[j]->label);
    }

    /* print qmatrix parameters header */
    if (opt_print_qmatrix)
    {
      if (locus[i]->model == BPP_DNA_MODEL_GTR)
      {
        for (j = 0; j < 6; ++j)
        {
          fprintf(fp_locus[i], "%s%s", tab_required ? "\t" : "", rlabels[j]);
          tab_required = 1;
        }
        for (j = 0; j < 4; ++j)
        {
          fprintf(fp_locus[i], "%s%s", tab_required ? "\t" : "", flabels[j]);
          tab_required = 1;
        }
      }
      else if (locus[i]->model == BPP_DNA_MODEL_K80)
      {
        fprintf(fp_locus[i],
                "%skappa",
                tab_required ? "\t" : "");
        tab_required = 1;
      }
      else if (locus[i]->model == BPP_DNA_MODEL_F81)
      {
        for (j = 0; j < 4; ++j)
        {
          fprintf(fp_locus[i], "%s%s", tab_required ? "\t" : "", flabels[j]);
          tab_required = 1;
        }
      }
      else if (locus[i]->model == BPP_DNA_MODEL_HKY)
      {
        fprintf(fp_locus[i], "%skappa", tab_required ? "\t" : "");
        tab_required = 1;
        for (j = 0; j < 4; ++j)
          fprintf(fp_locus[i], "\t%s", flabels[j]);
      }
      else if (locus[i]->model == BPP_DNA_MODEL_F84)
      {
        fprintf(fp_locus[i], "%skappa", tab_required ? "\t" : "");
        tab_required = 1;
        for (j = 0; j < 4; ++j)
          fprintf(fp_locus[i], "\t%s", flabels[j]);
      }
      else if (locus[i]->model == BPP_DNA_MODEL_T92)
      {
        fprintf(fp_locus[i],
                "%skappa\tpi_GC",
                tab_required ? "\t" : "");
        tab_required = 1;
      }
      else if (locus[i]->model == BPP_DNA_MODEL_TN93)
      {
        fprintf(fp_locus[i], "%skappa1\tkappa2", tab_required ? "\t" : "");
        tab_required = 1;
        for (j = 0; j < 4; ++j)
          fprintf(fp_locus[i], "\t%s", flabels[j]);
      }
      else
      {
        /* TODO: Implement other models */
        assert(locus[i]->model == BPP_DNA_MODEL_JC69);
      }

      if (opt_alpha_cats > 1)
      {
        fprintf(fp_locus[i],
                "%salpha",
                tab_required ? "\t" : "");
        tab_required = 1;
      }
    }
    if (tab_required)
      fprintf(fp_locus[i], "\n");
  }
}

static void mcmc_printinitial(FILE * fp, stree_t * stree)
{
  char * newick = stree_export_newick(stree->root, cb_serialize_branch);
  fprintf(fp, "%s\n", newick);
  free(newick);
}

static void print_rates(FILE ** fp_locus,
                        stree_t * stree,
                        gtree_t ** gtree,
                        locus_t ** locus) 
{
  int tab_required = 0;
  long i,j;
  unsigned int total_nodes;

  /* indices in TCAG order */
  long findices[4] = {3,1,0,2};
  long rindices[6] = {4, 2, 5, 0, 3, 1};

  assert(!opt_msci || (opt_msci && !opt_est_stree));

  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  for (i = 0; i < opt_locus_count; ++i)
  {
    tab_required = 0;
    /* print heredity scalars */
    if (opt_est_heredity == HEREDITY_ESTIMATE && opt_print_hscalars)
    {
      fprintf(fp_locus[i],
              "%s%.6f",
              tab_required ? "\t" : "", locus[i]->heredity[0]);
      tab_required = 1;
    }

    /* print mu_i and nu_i */
    if (opt_est_locusrate == MUTRATE_ESTIMATE && opt_print_locusrate)
    {
      fprintf(fp_locus[i],
              "%s%.6f",
              tab_required ? "\t" : "", gtree[i]->rate_mui);
      tab_required = 1;
    }
    if (opt_clock != BPP_CLOCK_GLOBAL && opt_print_rates)
    {
      fprintf(fp_locus[i],
              "%s%.6f",
              tab_required ? "\t" : "", gtree[i]->rate_nui);
      tab_required = 1;
    }

    /* print r_i */
    if (opt_clock != BPP_CLOCK_GLOBAL && opt_print_rates)
    {
      /* first one is tip, it always have a branch rate */
      fprintf(fp_locus[i],
              "%s%.6f",
              tab_required ? "\t" : "", stree->nodes[0]->brate[i]);
      tab_required = 1;
      for (j = 1; j < total_nodes; ++j)
        if (stree->nodes[j]->brate)
          fprintf(fp_locus[i], "\t%.6f", stree->nodes[j]->brate[i]);
    }

    if (opt_print_qmatrix)
    {
      if (locus[i]->model == BPP_DNA_MODEL_GTR)
      {
        for (j = 0; j < 6; ++j)
        {
          fprintf(fp_locus[i],
                  "%s%.6f",
                  tab_required ? "\t" : "",
                  locus[i]->subst_params[0][rindices[j]]);
          tab_required = 1;
        }
        for (j = 0; j < 4; ++j)
        {
          fprintf(fp_locus[i],
                  "%s%.6f",
                  tab_required ? "\t" : "",
                  locus[i]->frequencies[0][findices[j]]);
          tab_required = 1;
        }
      }
      else if (locus[i]->model == BPP_DNA_MODEL_K80)
      {
        fprintf(fp_locus[i],
                "%s%.6f",
                tab_required ? "\t" : "",
                locus[i]->subst_params[0][0]/locus[i]->subst_params[0][1]);
        tab_required = 1;
      }
      else if (locus[i]->model == BPP_DNA_MODEL_F81)
      {
        for (j = 0; j < 4; ++j)
        {
          fprintf(fp_locus[i],
                  "%s%.6f",
                  tab_required ? "\t" : "",
                  locus[i]->frequencies[0][findices[j]]);
          tab_required = 1;
        }
      }
      else if (locus[i]->model == BPP_DNA_MODEL_HKY)
      {
        fprintf(fp_locus[i],
                "%s%.6f",
                tab_required ? "\t" : "",
                locus[i]->subst_params[0][0] / locus[i]->subst_params[0][1]);
        tab_required = 1;
        for (j = 0; j < 4; ++j)
          fprintf(fp_locus[i], "\t%.6f", locus[i]->frequencies[0][findices[j]]);
      }
      else if (locus[i]->model == BPP_DNA_MODEL_F84)
      {
        fprintf(fp_locus[i],
                "%s%.6f",
                tab_required ? "\t" : "",
                locus[i]->subst_params[0][0] / locus[i]->subst_params[0][1]);
        tab_required = 1;
        for (j = 0; j < 4; ++j)
          fprintf(fp_locus[i], "\t%.6f", locus[i]->frequencies[0][findices[j]]);
      }
      else if (locus[i]->model == BPP_DNA_MODEL_T92)
      {
        fprintf(fp_locus[i],
                "%s%.6f\t%.6f",
                tab_required ? "\t" : "",
                locus[i]->subst_params[0][0]/locus[i]->subst_params[0][1],
                locus[i]->frequencies[0][1]+locus[i]->frequencies[0][2]);
        tab_required = 1;
      }
      else if (locus[i]->model == BPP_DNA_MODEL_TN93)
      {
        fprintf(fp_locus[i],
                "%s%.6f\t%.6f",
                tab_required ? "\t" : "",
                locus[i]->subst_params[0][0]/locus[i]->subst_params[0][2],
                locus[i]->subst_params[0][1]/locus[i]->subst_params[0][2]);
        tab_required = 1;
        for (j = 0; j < 4; ++j)
          fprintf(fp_locus[i], "\t%.6f", locus[i]->frequencies[0][findices[j]]);
      }
      else
      {
        /* TODO: Implement other models */
        assert(locus[i]->model == BPP_DNA_MODEL_JC69);
      }

      if (opt_alpha_cats > 1)
      {
        fprintf(fp_locus[i],
                "%s%f",
                tab_required ? "\t" : "", locus[i]->rates_alpha);
        tab_required = 1;
      }
        
    }
    if (tab_required)
      fprintf(fp_locus[i], "\n");
  }
}

static void mcmc_logsample(FILE * fp,
                           int step,
                           stree_t * stree,
                           gtree_t ** gtree,
                           locus_t ** locus,
                           long dparam_count,
                           long ndspecies)
{
  unsigned int i,j;
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
        fprintf(fp, "\t%.6f", stree->nodes[i]->theta);
  }

  /* then for inner nodes */
  if (opt_est_theta)
  {
    /* TODO: Is the 'has_theta' check also necessary ? */
    for (i = stree->tip_count; i < snodes_total; ++i)
      if (stree->nodes[i]->theta >= 0)
        fprintf(fp, "\t%.6f", stree->nodes[i]->theta);
  }

  /* 2. Print taus for inner nodes */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->tau)
      fprintf(fp, "\t%.6f", stree->nodes[i]->tau);

  /* 2a. Print phi for hybridization nodes */
  if (opt_msci)
  {
    unsigned int offset=stree->tip_count+stree->inner_count;

    for (i = 0; i < stree->hybrid_count; ++i)
    {
      #if 0

      /* old code prior the introduction of has_phi */

      snode_t * tmpnode = stree->nodes[offset+i];
      if (!node_is_bidirection(tmpnode))
      {
        /* hybridization node */

        /* if main node htau==0 and mirror node htau==1 then that is the only
           case we use the phi from the main node */
        if (tmpnode->hybrid->htau == 0 && tmpnode->htau == 1)
          tmpnode = tmpnode->hybrid;
      }
      #else

      /* mew correct code */

      snode_t * tmpnode = stree->nodes[offset+i];
      if (!tmpnode->has_phi)
        tmpnode = tmpnode->hybrid;
      #endif

      fprintf(fp, "\t%.6f", tmpnode->hphi);
    }
  }

  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      opt_est_mubar &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    fprintf(fp,"\t%.6f",stree->locusrate_mubar);
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      fprintf(fp,"\t%.6f", stree->locusrate_nubar);
    else
      fprintf(fp,"\t%.6f", stree->nui_sum / opt_locus_count);
  }

  if (opt_migration)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        if (opt_mig_bitmatrix[i][j])
          fprintf(fp, "\t%.6f", opt_migration_matrix[i][j]);
  }

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

static void print_gtree(FILE ** fp, stree_t * stree, gtree_t ** gtree)
{
  long i,j;
  double tl;

  /* TODO: For IM, branch lengths are incorrect */
  assert(!opt_migration);

  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_update_branchlengths(stree,gtree[i]);
    for (tl = 0, j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      if (!gtree[i]->nodes[j]->parent) continue;

      tl += gtree[i]->nodes[j]->length;
    }

    char * newick = gtree_export_newick(gtree[i]->root,NULL);
    fprintf(fp[i], "%s [TH=%.6f, TL=%.6f]\n", newick, gtree[i]->root->time, tl);
    free(newick);
  }
}

static void empirical_base_freqs_dna(msa_t * msa,
                                     unsigned int * weights,
                                     const unsigned int * pll_map)
{
  int i,j,k;
  int index;
  unsigned int code;
  const int states = 4;         /* ACGT */


  assert(pll_map == pll_map_nt);
  assert(pll_map['A'] == 1 &&
         pll_map['C'] == 2 &&
         pll_map['G'] == 4 &&
         pll_map['T'] == 8);

  for (i = 0; i < msa->count; ++i)
  {
    for (j = 0; j < msa->length; ++j)
    {
      code  = pll_map[(int)msa->sequence[i][j]]; 
      index = 0;

      for (k = 0; k < states; ++k)
      {
        if (code & 1)
          msa->freqs[index] += weights[j];
        code >>= 1;
        index++;
      }
    }
  }

  double sum = 0;
  for (i = 0; i < states; ++i)
    sum += msa->freqs[i];
  for (i = 0; i < states; ++i)
    msa->freqs[i] /= sum;
}

static void compute_base_freqs(msa_t * msa, unsigned int * weights, const unsigned int * pll_map)
{
  long states;

  if (msa->dtype == BPP_DATA_DNA)
  {
    states = 4;
    msa->freqs = (double *)xcalloc((size_t)states,sizeof(double));

    switch (msa->model)
    {
      case BPP_DNA_MODEL_JC69:
      case BPP_DNA_MODEL_K80:
        msa->freqs[0] = msa->freqs[1] = msa->freqs[2] = msa->freqs[3] = 0.25;
        break;
      case BPP_DNA_MODEL_F81:
      case BPP_DNA_MODEL_HKY:
      case BPP_DNA_MODEL_T92:
      case BPP_DNA_MODEL_TN93:
      case BPP_DNA_MODEL_F84:
      case BPP_DNA_MODEL_GTR:
        empirical_base_freqs_dna(msa,weights,pll_map);
        break;
      default:
        assert(0);
    }
  }
  else if (msa->dtype == BPP_DATA_AA)
  {
    states = 20;
    msa->freqs = (double *)xcalloc((size_t)states,sizeof(double));

    switch (msa->model)
    {
      case BPP_AA_MODEL_DAYHOFF:
        memcpy(msa->freqs,pll_aa_freqs_dayhoff,states*sizeof(double));
        break;
      case BPP_AA_MODEL_LG:
        memcpy(msa->freqs,pll_aa_freqs_lg,states*sizeof(double));
        break;
      case BPP_AA_MODEL_DCMUT:
        memcpy(msa->freqs,pll_aa_freqs_dcmut,states*sizeof(double));
        break;
      case BPP_AA_MODEL_JTT:
        memcpy(msa->freqs,pll_aa_freqs_jtt,states*sizeof(double));
        break;
      case BPP_AA_MODEL_MTREV:
        memcpy(msa->freqs,pll_aa_freqs_mtrev,states*sizeof(double));
        break;
      case BPP_AA_MODEL_WAG:
        memcpy(msa->freqs,pll_aa_freqs_wag,states*sizeof(double));
        break;
      case BPP_AA_MODEL_RTREV:
        memcpy(msa->freqs,pll_aa_freqs_rtrev,states*sizeof(double));
        break;
      case BPP_AA_MODEL_CPREV:
        memcpy(msa->freqs,pll_aa_freqs_cprev,states*sizeof(double));
        break;
      case BPP_AA_MODEL_VT:
        memcpy(msa->freqs,pll_aa_freqs_vt,states*sizeof(double));
        break;
      case BPP_AA_MODEL_BLOSUM62:
        memcpy(msa->freqs,pll_aa_freqs_blosum62,states*sizeof(double));
        break;
      case BPP_AA_MODEL_MTMAM:
        memcpy(msa->freqs,pll_aa_freqs_mtmam,states*sizeof(double));
        break;
      case BPP_AA_MODEL_MTART:
        memcpy(msa->freqs,pll_aa_freqs_mtart,states*sizeof(double));
        break;
      case BPP_AA_MODEL_MTZOA:
        memcpy(msa->freqs,pll_aa_freqs_mtzoa,states*sizeof(double));
        break;
      case BPP_AA_MODEL_PMB:
        memcpy(msa->freqs,pll_aa_freqs_pmb,states*sizeof(double));
        break;
      case BPP_AA_MODEL_HIVB:
        memcpy(msa->freqs,pll_aa_freqs_hivb,states*sizeof(double));
        break;
      case BPP_AA_MODEL_HIVW:
        memcpy(msa->freqs,pll_aa_freqs_hivw,states*sizeof(double));
        break;
      case BPP_AA_MODEL_JTTDCMUT:
        memcpy(msa->freqs,pll_aa_freqs_jttdcmut,states*sizeof(double));
        break;
      case BPP_AA_MODEL_FLU:
        memcpy(msa->freqs,pll_aa_freqs_flu,states*sizeof(double));
        break;
      case BPP_AA_MODEL_STMTREV:
        memcpy(msa->freqs,pll_aa_freqs_stmtrev,states*sizeof(double));
        break;
      default:
        assert(0);
    }
  }
}

static void create_mig_bitmatrix(stree_t * stree)
{
  long int i,j,s,t;

  assert(!opt_msci);
  assert(opt_migration);

  unsigned int total_nodes = stree->tip_count + stree->inner_count;

  opt_mig_bitmatrix = (long **)xmalloc((size_t)total_nodes*sizeof(long *));
  for (i = 0; i < total_nodes; ++i)
    opt_mig_bitmatrix[i] = (long *)xcalloc((size_t)total_nodes,sizeof(long));
  opt_migration_matrix = (double **)xmalloc((size_t)total_nodes*sizeof(double *));
  for (i = 0; i < total_nodes; ++i)
    opt_migration_matrix[i] = (double *)xcalloc((size_t)total_nodes,sizeof(double));

  /* go through source-target pairs */
  for (i = 0; i < opt_migration; ++i)
  {
    s = t = -1;
    for (j = 0; j < total_nodes; ++j)
    {
      if (s == -1 && !strcmp(stree->nodes[j]->label,opt_mig_source[i]))
        s = j;
      if (t == -1 && !strcmp(stree->nodes[j]->label,opt_mig_target[i]))
        t = j;
    }
    if (s == -1 || t == -1)
      fatal("Cannot create migration band for pair %s %s",
            opt_mig_source[i], opt_mig_target[i]);
    if (s == t)
      fatal("Cannot create migration from one species to itself (species: %s)",
            stree->nodes[s]->label);

    opt_mig_bitmatrix[s][t] = 1;

    opt_migration_matrix[s][t] = opt_mig_alpha / opt_mig_beta;

  }

  for (i = 0; i < opt_migration; ++i)
  {
    free(opt_mig_target[i]);
    free(opt_mig_source[i]);
  }
  free(opt_mig_target);
  free(opt_mig_source);
  opt_mig_target = opt_mig_source = NULL;
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
                     long * ptr_ft_round_snl,
                     long * ptr_pjump_spr,
                     long * ptr_pjump_snl,

                     double * ptr_mean_logl,
                     long ** ptr_mrate_row,
                     long ** ptr_mrate_col,
                     long ** ptr_mrate_round,
                     double ** ptr_mrate,
                     double ** ptr_mean_tau,
                     double ** ptr_mean_theta,
                     double ** ptr_mean_phi,
                     long * ptr_mrate_count,
                     long * ptr_mean_tau_count,
                     long * ptr_mean_theta_count,
                     long * ptr_mean_phi_count,
                     stree_t ** ptr_sclone, 
                     gtree_t *** ptr_gclones,
                     FILE *** ptr_fp_gtree,
                     FILE *** ptr_fp_locus,
                     FILE ** ptr_fp_out)
{
  long i,j;
  FILE * fp_mcmc;
  FILE * fp_out;
  long mcmc_offset;
  long out_offset;
  long * gtree_offset;
  long * rates_offset;
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
                  &rates_offset,
                  ptr_dparam_count,
                  ptr_posterior,
                  ptr_pspecies,
                  ptr_ft_round_rj,
                  ptr_pjump_rj,
                  ptr_ft_round_spr,
                  ptr_ft_round_snl,
                  ptr_pjump_spr,
                  ptr_pjump_snl,
                  ptr_mean_logl,
                  ptr_mrate_row,
                  ptr_mrate_col,
                  ptr_mrate_round,
                  ptr_mrate,
                  ptr_mean_tau,
                  ptr_mean_theta,
                  ptr_mean_phi,
                  ptr_mrate_count,
                  ptr_mean_tau_count,
                  ptr_mean_theta_count,
                  ptr_mean_phi_count,
                  &prec_logpr,
                  &prec_logl);

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
      xasprintf(&s, "%s.gtree.L%d", opt_outfile, (*ptr_gtree)[i]->original_index+1);
      gtree_files[i] = s;
      checkpoint_truncate(s,gtree_offset[i]);
    }
    free(gtree_offset);
  }

  /* truncate rate files if available */
  if (opt_print_locusfile)
  {
    assert(rates_offset);
    
    for (i = 0; i < opt_locus_count; ++i)
    {
      char * s = NULL;
      xasprintf(&s,template_ratesfile,(*ptr_gtree)[i]->original_index+1);
      checkpoint_truncate(s,rates_offset[i]);
      free(s);
    }
    free(rates_offset);
  }

  gtree_t ** gtree = *ptr_gtree;
  stree_t  * stree = *ptr_stree;

  //ANNA
  if (opt_datefile) 
	  fatal("Check pointing is not yet implemented for tip dating");
  gtree_alloc_internals(gtree,opt_locus_count,stree->inner_count, 0);
  reset_gene_leaves_count(stree,gtree);
  stree_reset_pptable(stree);


  /* compute MSC density */
  locus_t ** locus = *ptr_locus;
  if (opt_est_theta)
  {
    if (opt_migration)
    {
      for (i = 0; i < opt_locus_count; ++i)
        gtree[i]->logpr = gtree_logprob_mig(stree,
                                            gtree[i],
                                            locus[i]->heredity[0],
                                            i,
                                            thread_index_zero);
    }
    else
    {
      for (i = 0; i < opt_locus_count; ++i)
        gtree[i]->logpr = gtree_logprob(stree,
                                        locus[i]->heredity[0],
                                        i,
                                        thread_index_zero);
    }
  }
  else
  {
    if (opt_migration)
      fatal("Integrating out thetas with IM model not implemented yet");

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

  /* open potential truncated rate files for appending */
  *ptr_fp_locus = NULL;
  if (opt_print_locusfile)
  {
    FILE ** fp_locus = (FILE **)xmalloc((size_t)opt_locus_count*sizeof(FILE *));
    for (i = 0; i < opt_locus_count; ++i)
    {
      char * s = NULL;
      xasprintf(&s, template_ratesfile, gtree[i]->original_index+1);
      if (!(fp_locus[i] = fopen(s, "a")))
        fatal("Cannot open file %s for appending...", s);
      free(s);
    }
    *ptr_fp_locus = fp_locus;
  }

  /* if we are infering the species tree or gene flow, then create another
     cloned copy of the species tree and gene trees */
  if (opt_est_stree || opt_migration)
  {
    assert(opt_msci == 0);
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
                   long* ptr_ft_round_spr,
                   long* ptr_ft_round_snl,
                   long* ptr_pjump_spr,
                   long* ptr_pjump_snl,
                   double * ptr_mean_logl,
                   stree_t ** ptr_sclone, 
                   gtree_t *** ptr_gclones,
                   FILE *** ptr_fp_gtree,
                   FILE *** ptr_fp_locus,
                   FILE ** ptr_fp_out)
{
  long i,j;
  long msa_count;
  long pindex;
  double logl,logpr;
  double logl_sum = 0;
  double logpr_sum = 0;
  double * pjump;
  list_t * map_list = NULL;
  list_t * date_list = NULL;
  stree_t * stree;
  const unsigned int * pll_map;
  FILE * fp_mcmc = NULL;
  FILE * fp_out;
  FILE ** fp_gtree;
  FILE ** fp_locus = NULL;
  msa_t ** msa_list;
  gtree_t ** gtree;
  locus_t ** locus;

  const long thread_index = 0;

  /* method 10 specific variables */
  long dparam_count = 0;

  /* method 01 specific variables */
  stree_t * sclone = NULL;
  gtree_t ** gclones = NULL;

  if (!(fp_out = fopen(opt_outfile, "w")))
    fatal("Cannot open file %s for writing...");
  *ptr_fp_out = fp_out;
  init_outfile(fp_out);


  /* load species tree */
  stree = load_tree_or_network();
  printf(" Done\n");

  if (opt_msci && opt_migration)
    fatal(BPP_ERROR " Cannot use isolation with migration (IM) and "
          "introgression (MSci) models together.");

  /* Show network */
  if (opt_msci)
  {
    if (opt_finetune_phi == -1)
      fatal("Missing finetune value for phi parameter");
    if (opt_clock == BPP_CLOCK_CORR)
      fatal("MSCi model with auto-correlated relaxed clock is not currently implemented.");

    print_network_table(stree,fp_out);
    print_network_table(stree,stdout);
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

  /* set the data type for each multiple sequence alignment */
  if (opt_partition_list)
  {
    assert(opt_partition_count > 0);
    if (opt_partition_list[opt_partition_count-1]->end != msa_count)
      fatal("Partition file %s differs in number of partitions (%ld) "
            "to the specified number of loci (%ld)",
            opt_partition_file, opt_partition_list[opt_partition_count-1]->end,
            msa_count);
  }

  assert(BPP_DNA_MODEL_CUSTOM > BPP_DNA_MODEL_MAX &&
         BPP_DNA_MODEL_CUSTOM < BPP_AA_MODEL_MIN);
  assert (BPP_DNA_MODEL_MAX < BPP_AA_MODEL_MIN);


  if (opt_partition_file)
  {
    /* we have a partition file that specifies the substitution model for
       each locus */
    assert(opt_model == BPP_DNA_MODEL_CUSTOM);
    assert(opt_partition_list);

    pindex = 0;
    for (i = 0; i < opt_locus_count; ++i)
    {
      assert((i+1) >= opt_partition_list[pindex]->start);
      if ((i+1) > opt_partition_list[pindex]->end)
        ++pindex;
      assert((i+1) >= opt_partition_list[pindex]->start &&
             (i+1) <= opt_partition_list[pindex]->end);

      msa_list[i]->dtype = opt_partition_list[pindex]->dtype;
      msa_list[i]->model = opt_partition_list[pindex]->model;
    }

    /* deallocate partition list */
    for (i = 0; i < opt_partition_count; ++i)
      free(opt_partition_list[i]);
    free(opt_partition_list);
    opt_partition_list = NULL;
  }
  else
  {
    /* all loci have the same substitution model */

    int dtype,model;

    assert(opt_model != BPP_DNA_MODEL_CUSTOM);

    if ((opt_model >= BPP_DNA_MODEL_MIN) && (opt_model <= BPP_DNA_MODEL_MAX))
      dtype = BPP_DATA_DNA;
    else if ((opt_model >= BPP_AA_MODEL_MIN) && (opt_model <= BPP_AA_MODEL_MAX))
      dtype = BPP_DATA_AA;
    else
      fatal("Internal error when selecting substitution model for all loci");

    model = opt_model;
    for (i = 0; i < opt_locus_count; ++i)
    {
      msa_list[i]->dtype = dtype;
      msa_list[i]->model = model;
    }
  }

  /* remove missing sequences */
  for (i = 0; i < msa_count; ++i)
  {
    int deleted = msa_remove_missing_sequences(msa_list[i]);
    if (deleted == -1)
      fatal("[ERROR]: Locus %ld contains missing sequences only.\n"
            "Please remove the locus and restart the analysis.\n");

    if (deleted)
    {
      fprintf(stdout,
              "[WARNING]: Removing %d missing sequences from locus %ld\n",
              deleted, i);
      fprintf(fp_out,
              "[WARNING]: Removing %d missing sequences from locus %ld\n",
              deleted, i);
    }
    msa_list[i]->original_index = i;
  }

  /* remove ambiguous sites */
  if (opt_cleandata)
  {
    printf("Removing sites containing ambiguous characters...");
    for (i = 0; i < msa_count; ++i)
    {
      if (msa_list[i]->dtype == BPP_DATA_AA)
        continue;

      if (!msa_remove_ambiguous(msa_list[i]))
        fatal("All sites in locus %d contain ambiguous characters",i);
    }
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
    int compress_method;

    if (msa_list[i]->dtype == BPP_DATA_DNA)
    {
      pll_map = pll_map_nt;
      if (msa_list[i]->model == BPP_DNA_MODEL_JC69)
        compress_method = COMPRESS_JC69;
      else if (msa_list[i]->model == BPP_DNA_MODEL_GTR)
        compress_method = COMPRESS_GENERAL;
      else
      {
        /* TODO: Custom compression routines for the various models */
        compress_method = COMPRESS_GENERAL;
      }
    }
    else if (msa_list[i]->dtype == BPP_DATA_AA)
    {
      pll_map = pll_map_aa;
      compress_method = COMPRESS_GENERAL;
    }
    else
      assert(0);

    msa_list[i]->freqs = NULL;

    /* NOTE: Original length is the length after opt_cleandata is applied */
    msa_list[i]->original_length = msa_list[i]->length;
    weights[i] = compress_site_patterns(msa_list[i]->sequence,
                                        pll_map,
                                        msa_list[i]->count,
                                        &(msa_list[i]->length),
                                        compress_method);

    /* compute base frequencies */
    compute_base_freqs(msa_list[i], weights[i], pll_map);
  }

  if (opt_diploid)
  {
    fprintf(stdout, "\nSummary of alignments *before* phasing sequences:");
    fprintf(fp_out, "\nSummary of alignments *before* phasing sequences:");
  }
  msa_summary(stdout, msa_list,msa_count);
  msa_summary(fp_out, msa_list,msa_count);

  /* parse map file */
  if (stree->tip_count > 1)
  {
    printf("Parsing map file...");
    if (!(map_list = parse_mapfile(opt_mapfile)))
      fatal("Failed parsing map file %s", opt_mapfile);
    else
      printf(" Done\n");
  }

  if (opt_datefile) {
	  //ANNA
        printf("Parsing date file...");
        if (! (date_list = parse_date_mapfile(opt_datefile)))
                fatal("Failed to parse date file %s", opt_datefile);
        else
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

  /* print compressed alignmens in output file */
  fprintf(fp_out, "COMPRESSED ALIGNMENTS\n\n");

  /* print the alignments */
  msa_print_phylip(fp_out,msa_list,msa_count, weights);

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
                                       msa_count);

    /* TODO: KEEP WEIGHTS */
    //for (i = 0; i < msa_count; ++i) free(weights[i]);

    mapping = (unsigned long **)xmalloc((size_t)msa_count *
                                        sizeof(unsigned long *));

    /* allocate temporary array for storing pattern weights for alignment A3 */
    unsigned int ** tmpwgt = (unsigned int **)xmalloc((size_t)(msa_count) *
                                                      sizeof(unsigned int *));
    for (i = 0; i < msa_count; ++i)
    {
      int compress_method;

      assert(msa_list[i]->dtype == BPP_DATA_DNA);
      if (msa_list[i]->dtype == BPP_DATA_DNA)
      {
        pll_map = pll_map_nt;
        if (msa_list[i]->model == BPP_DNA_MODEL_JC69)
          compress_method = COMPRESS_JC69;
        else if (msa_list[i]->model == BPP_DNA_MODEL_GTR)
          compress_method = COMPRESS_GENERAL;
        else
        {
          /* TODO: Custom compression routines for the various models */
          compress_method = COMPRESS_GENERAL;
        }
      }
      else if (msa_list[i]->dtype == BPP_DATA_AA)
      {
        pll_map = pll_map_aa;
        compress_method = COMPRESS_GENERAL;
      }
      else
        assert(0);

      /* compress again for JC69 and get mappings */
      mapping[i] = compress_site_patterns_diploid(msa_list[i]->sequence,
                                                  pll_map,
                                                  msa_list[i]->count,
                                                  &(msa_list[i]->length),
                                                  tmpwgt+i,
                                                  compress_method);
    }
    fprintf(fp_out, "COMPRESSED ALIGNMENTS AFTER PHASING OF DIPLOID SEQUENCES\n\n");
    msa_print_phylip(fp_out,msa_list,msa_count,tmpwgt);

    /* deallocate temporary pattern weights */
    for (i = 0; i < msa_count; ++i)
      free(tmpwgt[i]);
    free(tmpwgt);

    fprintf(stdout, "\nSummary of alignments *after* phasing sequences:");
    fprintf(fp_out, "\nSummary of alignments *after* phasing sequences:");
    msa_summary(stdout, msa_list,msa_count);
    msa_summary(fp_out, msa_list,msa_count);
  }

  /* Pin master thread for NUMA first policy touch
     TODO: Perhaps move this to an earlier point */
  if (opt_threads > 1)
  {
    long * indices = threads_load_balance(msa_list);
    if (indices)
    {
      assert(msa_count == opt_locus_count);

      unsigned long ** tmp_mapping = NULL;
      unsigned long ** tmp_rescount = NULL;
      unsigned int ** tmp_weights = NULL;
      int * tmp_unphased_length = NULL;

      /* allocate temporary arrays */
      if (opt_diploid)
      {
        tmp_mapping = (unsigned long **)xmalloc((size_t)msa_count *
                                                sizeof(unsigned int long *));
        tmp_rescount = (unsigned long **)xmalloc((size_t)msa_count *
                                                 sizeof(unsigned int long *));
        tmp_unphased_length = (int *)xmalloc((size_t)msa_count * sizeof(int));
      }
      tmp_weights = (unsigned int **)xmalloc((size_t)msa_count *
                                             sizeof(unsigned int *));

      /* reorder other arrays to the new positions of msa_list */
      for (i = 0; i < opt_locus_count; ++i)
      {
        if (opt_diploid)
        {
          tmp_mapping[i]         = mapping[indices[i]];
          tmp_rescount[i]        = resolution_count[indices[i]];
          tmp_unphased_length[i] = unphased_length[indices[i]];
        }
        tmp_weights[i]         = weights[indices[i]];
      }

      if (opt_diploid)
      {
        memmove(mapping,tmp_mapping,opt_locus_count * sizeof(unsigned long *));
        memmove(resolution_count,tmp_rescount,opt_locus_count * sizeof(unsigned long *));
        memmove(unphased_length,tmp_unphased_length,opt_locus_count*sizeof(int));
      }
      memmove(weights,tmp_weights,opt_locus_count*sizeof(unsigned int *));

      if (opt_diploid)
      {
        free(tmp_mapping);
        free(tmp_rescount);
        free(tmp_unphased_length);
      }
      free(tmp_weights);

      for (i = 0; i < opt_locus_count; ++i)
        msa_list[i]->original_index = indices[i];
      free(indices);
    }
    threads_pin_master();
  }

  /* gene tree and locus output files */
  *ptr_fp_gtree = NULL;

  /* if print gtree */
  if (opt_print_genetrees)
  {
    fp_gtree = (FILE **)xmalloc((size_t)opt_locus_count*sizeof(FILE *));
    for (i = 0; i < opt_locus_count; ++i)
    {
      char * s = NULL;
      xasprintf(&s, "%s.gtree.L%d", opt_outfile, msa_list[i]->original_index+1);
      fp_gtree[i] = xopen(s,"w");
      free(s);
    }
  }
  else
    fp_gtree = NULL;

  *ptr_fp_gtree = fp_gtree;

  /* if print rates */
  *ptr_fp_locus = NULL;
  if (opt_print_locusfile)
  {
    fp_locus = (FILE **)xmalloc((size_t)opt_locus_count*sizeof(FILE *));
    for (i = 0; i < opt_locus_count; ++i)
    {
      char * s = NULL;
      xasprintf(&s, template_ratesfile, msa_list[i]->original_index+1);
      if (!(fp_locus[i] = fopen(s, "w")))
        fatal("Cannot open file %s for appending...", s);
      free(s);
    }
    *ptr_fp_locus = fp_locus;
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

  /* initialize migration (presence/absence) matrix */
  if (opt_migration)
  {
    stree_label(stree);
    create_mig_bitmatrix(stree);
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

    if (!opt_est_theta)
    {
      double hfactor = 0;
      for (i = 0; i < opt_locus_count; ++i)
        hfactor -= (msa_list[i]->count-1)*log(heredity[i]);
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
    if (!opt_est_theta)
    {
      double hfactor = 0;
      for (i = 0; i < opt_locus_count; ++i)
        hfactor -= (msa_list[i]->count-1)*log(heredity[i]);
      stree->notheta_hfactor = hfactor;
    }
  }

  /* initialize locus mutation rates if estimation was selected */
  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR ||
       opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR))
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
    opt_est_locusrate = MUTRATE_CONSTANT;
  }

  /* We must first link tip sequences (gene tips) to populations */
  if (opt_est_delimit)          /* species delimitation */
  {
    assert(opt_msci == 0);
    stree_rootdist(stree,map_list,msa_list,weights);
  }

  if (opt_migration)
    stree_update_mig_subpops(stree, thread_index);

  gtree = gtree_init(stree,msa_list,map_list,date_list,msa_count);
  free(date_list);
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree[i]->original_index = msa_list[i]->original_index;
    gtree[i]->msa_index = i;
  }



  /* Generate space for cloning the species and gene trees (used for species
     tree inference or inference with migration) */
  if (opt_est_stree || opt_migration)
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

  /* ensure that heredity is now set to either 0 or 1, and locusrate is not set
     to MUTATE_FROMFILE */
  assert(opt_est_locusrate >= MUTRATE_CONSTANT && opt_est_locusrate != MUTRATE_FROMFILE);

  gtree_update_branch_lengths(gtree, msa_count);

    /* set universal mean/variance */
  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
  {
    if (opt_est_mubar)
      stree->locusrate_mubar = opt_mubar_alpha / opt_mubar_beta;
    else
      stree->locusrate_mubar = 1;
  }
  if (opt_clock != BPP_CLOCK_GLOBAL)
    stree->locusrate_nubar = opt_vbar_alpha / opt_vbar_beta;

  stree->nui_sum = 0;

  if (opt_msci && !opt_est_theta)
  {
    for (j = 0; j < stree->tip_count + stree->inner_count+stree->hybrid_count; ++j)
    {
      snode_t * x = stree->nodes[j];
      if (x->hybrid)
      {
        x->notheta_old_phi_contrib = (double *)xcalloc((size_t)opt_locus_count,
                                                       sizeof(double));
        x->notheta_phi_contrib = (double *)xcalloc((size_t)opt_locus_count,
                                                   sizeof(double));
      }
      else
        x->notheta_phi_contrib = x->notheta_old_phi_contrib = NULL;
    }
  }
  else
  {
    for (j = 0; j < stree->tip_count + stree->inner_count+stree->hybrid_count; ++j)
    {
      snode_t * x = stree->nodes[j];
      x->notheta_old_phi_contrib = NULL;
      x->notheta_phi_contrib = NULL;
    }
  }

  for (i = 0, pindex=0; i < msa_count; ++i)
  {
    int states = 0;
    unsigned int pmatrix_count = gtree[i]->edge_count;
    msa_t * msa = msa_list[i];
    unsigned int scale_buffers = opt_scaling ? 2*gtree[i]->inner_count : 0;

    /* activate twice as many transition probability matrices (for reverting in
       locusrate, species tree SPR and mixing proposals)  */
    pmatrix_count *= 2;               /* double to account for cloned */

    /* TODO: In the future we can allocate double amount of p-matrices
       for the other methods as well in order to speedup rollback when
       rejecting proposals */

    if (msa_list[i]->dtype == BPP_DATA_DNA)
    {
      //assert(msa_list[i]->model == BPP_DNA_MODEL_JC69);
      states = 4;
      pll_map = pll_map_nt;
    }
    else if (msa_list[i]->dtype == BPP_DATA_AA)
    {
      states = 20;
      pll_map = pll_map_aa;
    }
    else
      fatal("Internal error when setting states for locus %ld", i);

    /* create the locus structure */
    locus[i] = locus_create((unsigned int)(msa_list[i]->dtype),        /* data type */
                            (unsigned int)(msa_list[i]->model),        /* subst model */
                            gtree[i]->tip_count,        /* # tip sequence */
                            2*gtree[i]->inner_count,    /* # CLV vectors */
                            states,                     /* # states */
                            msa->length,                /* sequence length */
                            rate_matrices,              /* subst matrices (1) */
                            pmatrix_count,              /* # prob matrices */
                            opt_alpha_cats,             /* # rate categories */
                            scale_buffers,              /* # scale buffers */
                            (unsigned int)opt_arch);    /* attributes */

    locus[i]->original_index = msa_list[i]->original_index;
    /* set frequencies and substitution rates */
    /* TODO: For GTR perhaps set to empirical frequencies */
    locus_set_frequencies_and_rates(locus[i]);

    if (opt_diploid)
    {
      for (j = 0; j < (long)(stree->tip_count); ++j)
        if (stree->nodes[j]->diploid)
        {
          locus[i]->diploid = 1;
          break;
        }
    }

    /* set rate of evolution and heredity scalar for each locus */
    gtree[i]->rate_mui = locusrate[i];
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
      pll_set_tip_states(locus[i], j, pll_map, msa_list[i]->sequence[j]);

    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    {
      long thread_index = 0;

      gtree[i]->rate_mui  = stree->locusrate_mubar;
      if (!opt_debug_rates)
        gtree[i]->rate_mui *= (.9 + .2*legacy_rndu(thread_index));
    }

    if (opt_clock != BPP_CLOCK_GLOBAL)
    {
      long thread_index = 0;

      gtree[i]->rate_nui = stree->locusrate_nubar;

      /* TODO: Remove if condition after debugging */
      if (!opt_debug_rates)
        gtree[i]->rate_nui *= (.9 + .2*legacy_rndu(thread_index));

      for (j = 0; j < stree->tip_count+stree->inner_count+stree->hybrid_count; ++j)
      {
        snode_t * node = stree->nodes[j];

        /* Mirror nodes in bidirectional introgression */
        if (opt_msci && node->hybrid)
        {
          if (node_is_hybridization(node) && !node->htau) continue;
          if (node_is_bidirection(node) && node_is_mirror(node)) continue;
        }

        /* set values around Log-Normal mean */
        /* TODO: Keep only else after debugging */
        if (opt_debug_rates)
          node->brate[i] = gtree[i]->rate_mui;
        else
        {
          if (opt_clock == BPP_CLOCK_CORR && !node->parent)  /* root node */
            node->brate[i] = gtree[i]->rate_mui;
          else
            node->brate[i] = gtree[i]->rate_mui *
                             (.9+.2*legacy_rndu(thread_index));
        }
      }
      gtree[i]->lnprior_rates   = lnprior_rates(gtree[i],stree,i);

      stree->nui_sum += gtree[i]->rate_nui;
    }

    /* compute the conditional probabilities for each inner node */
    locus_update_matrices(locus[i],
                          gtree[i],
                          gtree[i]->nodes,
                          stree,
                          i,
                          gtree[i]->edge_count);
    locus_update_partials(locus[i],
                          gtree[i]->nodes+gtree[i]->tip_count,
                          gtree[i]->inner_count);

    /* now that we computed the CLVs, calculate the log-likelihood for the
       current gene tree */
    logl = locus_root_loglikelihood(locus[i],
                                    gtree[i]->root,
                                    locus[i]->param_indices,
                                    NULL);
    logl_sum += logl;
    if (isinf(logl))
      fatal("\n[ERROR] log-L for locus %d is -inf.\n"
            "Please run BPP with numerical scaling. This is enabled by adding the line:\n"
            "\n  scaling = 1\n\nto the control file", i+1);

    /* store current log-likelihood in each gene tree structure */
    gtree[i]->logl = logl;
    if (opt_est_theta)
    {
      if (opt_migration)
        logpr = gtree_logprob_mig(stree,
                                  gtree[i],
                                  locus[i]->heredity[0],
                                  i,
                                  thread_index_zero);
      else
        logpr = gtree_logprob(stree,locus[i]->heredity[0],i,thread_index_zero);
      gtree[i]->logpr = logpr;
      logpr_sum += logpr;
    }
    else
    {
      if (opt_migration)
        fatal("Integrating out thetas not implemented yet for IM model");

      for (j = 0; j < stree->tip_count + stree->inner_count+stree->hybrid_count; ++j)
        logpr_sum += gtree_update_logprob_contrib(stree->nodes[j],
                                                  locus[i]->heredity[0],
                                                  i,
                                                  thread_index_zero);
    }
  }
  if (!opt_est_theta)
  {
    logpr_sum = 0;
    for (j = 0; j < stree->tip_count + stree->inner_count+stree->hybrid_count; ++j)
      logpr_sum += stree->nodes[j]->notheta_logpr_contrib;
    stree->notheta_logpr += logpr_sum;
    stree->notheta_old_logpr = 0;

    logpr_sum = stree->notheta_logpr;
  }

  /* Reading constraints file */
  if (opt_constraintfile)
  {
    fprintf(stdout, "Reading constraint file %s\n", opt_constraintfile);
    fprintf(fp_out, "Reading constraint file %s\n", opt_constraintfile);
    parse_and_set_constraints(stree, fp_out);
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

  fprintf(stdout,"\nInitial MSC density and log-likelihood of observing data:\n");
  fprintf(stdout,"log-PG0 = %f   log-L0 = %f\n\n", logpr_sum, logl_sum);
  fprintf(fp_out,"\nInitial MSC density and log-likelihood of observing data:\n");
  fprintf(fp_out,"log-PG0 = %f   log-L0 = %f\n\n", logpr_sum, logl_sum);

  if (isinf(logl_sum))
    fatal("\n[ERROR] The sum of log-L for all loci is -inf.\n"
          "Please run BPP with numerical scaling. This is enabled by adding the line:\n"
          "\n  scaling = 1\n\nto the control file");

  #if 0
  debug_print_stree(stree);
  assert(0);
  #endif
  /* free weights array */
  free(weights);

  if (opt_est_delimit)          /* species delimitation */
    rj_init(gtree,stree,msa_count);

  /* initialize pjump and finetune rounds */

  int pjump_size = PROP_COUNT + 1+1 + GTR_PROP_COUNT + CLOCK_PROP_COUNT + !!opt_migration;
  pjump = (double *)xcalloc(pjump_size,sizeof(double));

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

    if (opt_print_locusfile)
      mcmc_printheader_rates(fp_locus,stree,locus);
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
  *ptr_ft_round_snl = 0;
  *ptr_pjump_spr = 0;
  *ptr_pjump_snl = 0;
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

#ifdef CHECK_LOGL
#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))
#define SWAP_PMAT_INDEX(e,i) (((e)+(i))%((e)<<1))
#define SWAP_SCALER_INDEX(n,i) (((n)+((i)-1))%(2*(n)-2)) 
#define GET_HINDEX(t,p) (((node_is_mirror((p)) ? \
                          (p)->node_index : (p)->hybrid->node_index)) - \
                        ((t)->tip_count+(t)->inner_count))
static void all_partials_recursive(gnode_t * node,
                                   unsigned int * trav_size,
                                   gnode_t ** outbuffer)
{
  if (!node->left)
    return;

  all_partials_recursive(node->left,  trav_size, outbuffer);
  all_partials_recursive(node->right, trav_size, outbuffer);

  outbuffer[*trav_size] = node;
  *trav_size = *trav_size + 1;
}

static void gtree_all_partials(gnode_t * root,
                               gnode_t ** travbuffer,
                               unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return;

  all_partials_recursive(root, trav_size, travbuffer);
}

static double debug_full_lh(stree_t * stree, gtree_t * gtree, locus_t * locus, long msa_index)
{
  unsigned int j,k;
  gnode_t ** buffer;
  double logl = 0;

  buffer = (gnode_t **)xcalloc((size_t)(gtree->tip_count+gtree->inner_count),
                               sizeof(gnode_t *));
  k=0;
  for (j = 0; j < gtree->tip_count + gtree->inner_count; ++j)
  {
    gnode_t * tmp = gtree->nodes[j];
    if (tmp->parent)
    {
      tmp->pmatrix_index = SWAP_PMAT_INDEX(gtree->edge_count,
                                           tmp->pmatrix_index);
      buffer[k++] = tmp;
    }
  }
  
  
  locus_update_matrices(locus,gtree,buffer,stree,msa_index,k);
  
  gtree_all_partials(gtree->root,buffer,&k);
  for (j = 0; j < k; ++j)
  {
    buffer[j]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                          buffer[j]->clv_index);
    if (opt_scaling)
      buffer[j]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                  buffer[j]->scaler_index);
  }
  
  locus_update_partials(locus,buffer,k);
  
  /* compute log-likelihood */
  assert(!gtree->root->parent);
  logl = locus_root_loglikelihood(locus,
                                  gtree->root,
                                  locus->param_indices,
                                  NULL);
  //printf("LOGL: %f    OLDD LOGL: %f", gtree->logl, old_logl);  
  free(buffer);
  return logl;
}

static void check_logl(stree_t * stree, gtree_t ** gtree, locus_t ** locus, long iter, const char * move)
{
  long i;
  double logl, old_logl;

  for (i = 0; i < opt_locus_count; ++i)
  {
    logl = debug_full_lh(stree,gtree[i],locus[i],i);

    old_logl = gtree[i]->logl;

    if (fabs(logl - old_logl) > 1e-9)
    {
      printf("[FATAL iter %ld locus %ld] LOGL: %f   OLD_LOGL: %f\n", iter, i, logl, old_logl);
      
      fatal("Invalid logl iter: %ld locus: %ld move: %s    correct logl: %f   wrong logl: %f", iter, i, move, logl, old_logl);
    }
  }
}

#endif
#ifdef CHECK_LNPRIOR
static void check_lnprior(stree_t * stree, gtree_t ** gtree, long iter, const char * move)
{
  long i;

  if (opt_clock == BPP_CLOCK_GLOBAL) return;
  /* */
  for (i = 0; i < opt_locus_count; ++i)
  {
    double debug_new_prior_rates = lnprior_rates(gtree[i], stree, i);
    if (fabs(debug_new_prior_rates - gtree[i]->lnprior_rates) > 1e-5)
    {
      printf("[FATAL iter %ld locus %ld] lnPrior: %f   OLD lnPrior: %f\n", iter, i, debug_new_prior_rates, gtree[i]->lnprior_rates);
      fatal("Invalid lnprior iter: %ld locus: %ld move: %s    correct lnprior: %f   wrong lnprior: %f", iter, i, move, debug_new_prior_rates, gtree[i]->lnprior_rates);
    }
  }
}
#endif
static void fill_mean_mrate_indices(stree_t * stree, long * row, long * col, long count)
{
  long i,j,k = 0;
  long total_nodes = stree->tip_count + stree->inner_count;

  for (i = 0; i < total_nodes && k < count; ++i)
    for (j = 0; j < total_nodes && k < count; ++j)
      if (opt_mig_bitmatrix[i][j])
      {
        row[k] = i;
        col[k] = j;
        ++k;
      }
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
  FILE ** fp_locus = NULL;
  gtree_t ** gtree;
  locus_t ** locus;
  long * gtree_offset = NULL;   /* for checkpointing when printing gene trees */
  long * rates_offset = NULL;
  double ratio;
  long ndspecies;

  /* method 10 specific variables */
  long dparam_count = 0;
  long ft_round_rj;
  double pjump_rj = 0;
  double * posterior = NULL;
  double * pspecies = NULL;

  /* method 01 specific variables */
  long ft_round_spr = 0, ft_round_snl = 0;
  long pjump_spr, pjump_snl;
  long printk;// = opt_samplefreq * opt_samples;
  double mean_logl = 0;
  double dbg_mean_rate = 0;
  #ifdef DEBUG_GTR
  double mean_freqa = 0;
  double mean_freqc = 0;
  double mean_freqg = 0;
  double mean_freqt = 0;

  double mean_ratea = 0;
  double mean_rateb = 0;
  double mean_ratec = 0;
  double mean_rated = 0;
  double mean_ratee = 0;
  double mean_ratef = 0;

  double mean_alpha0 = 0;
  #endif

  long * mean_mrate_row = NULL;
  long * mean_mrate_col = NULL;
  long * mean_mrate_round = NULL;
  double * mean_mrate = NULL;
  double * mean_tau = NULL;
  double * mean_theta = NULL;
  double * mean_phi = NULL;

  long mean_mrate_count = 0;
  long mean_theta_count = 0;
  long mean_tau_count = 0;
  long mean_phi_count = 0;

  stree_t * sclone;
  gtree_t ** gclones;

  unsigned long curstep = 0;

  printf("\nStarting timer..\n");
  timer_start();

  /* TODO: Decide whether to have a different option for printing rates */
  if (opt_clock != BPP_CLOCK_GLOBAL)
    opt_print_rates = opt_print_locusrate;

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
                     &ft_round_snl,
                     &pjump_spr,
                     &pjump_snl,
                     &mean_logl,
                     &mean_mrate_row,
                     &mean_mrate_col,
                     &mean_mrate_round,
                     &mean_mrate,
                     &mean_tau,
                     &mean_theta,
                     &mean_phi,
                     &mean_mrate_count,
                     &mean_tau_count,
                     &mean_theta_count,
                     &mean_phi_count,
                     &sclone, 
                     &gclones,
                     &fp_gtree,
                     &fp_locus,
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
                   &ft_round_snl,
                   &pjump_spr,
                   &pjump_snl,
                   &mean_logl,
                   &sclone, 
                   &gclones,
                   &fp_gtree,
                   &fp_locus,
                   &fp_out);

    /* allocate mean_mrate, mean_tau, mean_theta */
    if (opt_migration)
    {
      long mig_size    = MIN(MAX_MRATE_OUTPUT,opt_migration);
      mean_mrate_row   = (long *)xcalloc((size_t)mig_size,sizeof(long));
      mean_mrate_col   = (long *)xcalloc((size_t)mig_size,sizeof(long));
      mean_mrate_round = (long *)xcalloc((size_t)mig_size,sizeof(long));
      mean_mrate       = (double *)xcalloc((size_t)mig_size, sizeof(double));

      fill_mean_mrate_indices(stree, mean_mrate_row, mean_mrate_col, mig_size);

      mean_mrate_count = mig_size;
    }
    if (opt_est_theta)
      mean_theta = (double *)xcalloc(MAX_THETA_OUTPUT,sizeof(double));
    mean_tau   = (double *)xcalloc(MAX_TAU_OUTPUT,sizeof(double));

    if (opt_msci)
      mean_phi = (double *)xcalloc(MAX_PHI_OUTPUT,sizeof(double));

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

  if (opt_migration && opt_msci)
    fatal("Cannot use both MSci and IM models together");
  if (opt_migration && opt_clock != BPP_CLOCK_GLOBAL)
    fatal("Cannot use IM model with relaxed clock models");

  if (opt_checkpoint && opt_print_genetrees)
    gtree_offset = (long *)xmalloc((size_t)opt_locus_count*sizeof(long));
  if (opt_checkpoint && opt_print_locusfile)
    rates_offset = (long *)xmalloc((size_t)opt_locus_count*sizeof(long));
  if (opt_exp_randomize)
    fprintf(stdout, "[EXPERIMENTAL] - Randomize nodes order on gtree SPR\n");
  if (opt_rev_gspr)
    fprintf(stdout, "[EXPERIMENTAL] - Revolutionary gene tree SPR algorithm\n");
  if (opt_revolutionary_spr_method)
    fprintf(stdout, "[EXPERIMENTAL] - Revolutionary species tree SPR algorithm\n");
  if (opt_exp_theta)
    fprintf(stdout, "[EXPERIMENTAL] - Theta proposal using a sliding window log(theta)\n");
  if (opt_exp_imrb)
    fprintf(stdout, "[EXPERIMENTAL] - New improved IM rubberband algorithm\n");

  assert(opt_theta_move==BPP_THETA_SLIDE || opt_theta_move==BPP_THETA_GIBBS);
  fprintf(stdout, "Theta proposal: %s\n",
          (opt_theta_move == BPP_THETA_SLIDE) ?
            "sliding window" : "gibbs sampler");

  if (opt_linkedtheta == BPP_LINKEDTHETA_NONE)
    fprintf(stdout, "Linked thetas: none\n");
  else if (opt_linkedtheta == BPP_LINKEDTHETA_ALL)
    fprintf(stdout, "Linked thetas: all nodes\n");
  else if (opt_linkedtheta == BPP_LINKEDTHETA_INNER)
    fprintf(stdout, "Linked thetas: inner nodes\n");
  else if (opt_linkedtheta == BPP_LINKEDTHETA_MSCI)
    fprintf(stdout, "Linked thetas: introgression nodes\n");

  /* enable proposals */
  if (opt_model != BPP_DNA_MODEL_JC69)
  {
    enabled_prop_qrates = 0;
    enabled_prop_freqs  = 0;
    for (i = 0; i < opt_locus_count; ++i)
    {
      if (locus[i]->qrates_param_count)
        enabled_prop_qrates = 1;
      if (locus[i]->freqs_param_count)
        enabled_prop_freqs = 1;
    }
  }

  /* enable/disable variables for mcmc headerline */
  if (opt_est_heredity == HEREDITY_ESTIMATE &&
      (opt_est_locusrate == MUTRATE_ESTIMATE &&
       opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR))
    enabled_lrht = 1;
  else if (opt_est_heredity == HEREDITY_ESTIMATE)
    enabled_hrdt = 1;
  else if (opt_est_locusrate == MUTRATE_ESTIMATE &&
           opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR)
    enabled_mui = 1;

  assert(!(enabled_hrdt && enabled_lrht));

  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
       opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
  {
    assert(!enabled_mui);
    enabled_mui = 1;
  }

  if (opt_est_locusrate == MUTRATE_ESTIMATE &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL && opt_est_mubar)
    enabled_mubar = 1;
  if (opt_clock != BPP_CLOCK_GLOBAL &&
      opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    enabled_nubar = 1;

  assert(!(enabled_mui && enabled_lrht));


  if (opt_threads > 1)
  {
    threads_lb_stats(locus, fp_out);
    threads_init();
    memset(&td,0,sizeof(td));
  }

  /* flush all open files */
  fflush(NULL);
  if (!opt_resume)
    timer_print("", " taken to read and process data..\nRestarting timer...\n\n",
                fp_out);
  timer_start();

  #if 0
  unsigned long total_steps = opt_samples * opt_samplefreq + opt_burnin;
  progress_init("Running MCMC...", total_steps);
  #endif
  if (opt_prob_snl && opt_clock != BPP_CLOCK_GLOBAL)
  {
    if (!opt_debug_full)
    {
      fprintf(stdout, "[DEBUG] SNL and relaxed clock requires full "
                      "recomputation of log-L (turning on)\n");
      opt_debug_full = 1;                 
    }
  }
  if (opt_debug_full)
    fprintf(stdout, "[DEBUG] Full recomputation of gene tree probabilities and "
                    "log-likelihood in SNL/SPR moves\n");
  if (opt_prob_snl && opt_est_stree)
  {
    fprintf(stdout, "[DEBUG] Using SNL move\n");
    fprintf(stdout,
            "[DEBUG] lambda_expand = %f   lambda_shrink = %f\n",
            opt_snl_lambda_expand, opt_snl_lambda_shrink);
    fprintf(stdout,
            "[DEBUG] SHRINK move proportion: %f\n",
            opt_prob_snl_shrink);
  }
  if (opt_exp_sim)
  {
    fprintf(stdout,
            "[DEBUG] Using experimental gene tree SPR move with simulation\n");
  }

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (locus[i]->rate_cats > 1)
      enabled_prop_alpha = 1;
  }

  printk = opt_samplefreq * opt_samples;

  /* check if summary only was requested (no MCMC) and initialize counter
     for MCMC loop appropriately */
  if (opt_onlysummary)
    i = opt_samples*opt_samplefreq;
  else
    i = curstep - opt_burnin;

  /* TODO: Delete after debugging */
  if (opt_debug_rates)
  {
    opt_seed = 1;
    legacy_init();
  }

  int delim_digit_count = 0;
  if (opt_method == METHOD_10)
    delim_digit_count = delimitation_getparam_count() ?
                      (int)floor(log10(abs(delimitation_getparam_count())))+1 : 1;
  else if (opt_method == METHOD_11)
  {
    long dels = delimitations_count(stree);
    delim_digit_count = dels ? (int)floor(log10(abs(dels)))+1 : 1;
  }

  print_mcmc_headerline(stdout,
                        stree,
                        gtree,
                        mean_mrate_row,
                        mean_mrate_col,
                        mean_mrate_count);
  if (!opt_resume)
  {
    print_mcmc_headerline(fp_out,
                          stree,
                          gtree,
                          mean_mrate_row,
                          mean_mrate_col,
                          mean_mrate_count);
    
    if (!opt_debug_start)
      opt_debug_start = 1;
    if (!opt_debug_end)
      opt_debug_end = opt_burnin+opt_samples*opt_samplefreq;
  }

  FILE * fp_debug = stdout;

  //Anna set taus and thetas (wait jk the taus need to be set before the gene trees are simulated
/*	stree->nodes[0]->theta =  .004; 
 	stree->nodes[1]->theta =  .004; 
 	stree->nodes[2]->theta =  .005; 
 	stree->nodes[3]->theta =  .006; 
 	stree->nodes[4]->theta =  .0035; */

	stree->nodes[0]->theta =  .007; 
 	stree->nodes[1]->theta =  .004; 
 	stree->nodes[2]->theta =  .005; 
 	stree->nodes[3]->theta =  .006; 
 	stree->nodes[4]->theta =  .008; 

  /* *** start of MCMC loop *** */
  for ( ; i < opt_samples*opt_samplefreq; ++i)
  {
    long print_newline = 0;
    #if 0
    /* update progress bar */
    if (!opt_quiet)
      progress_update(curstep);
    #endif

    ++opt_debug_counter;

    /* reset finetune parameters */
    if (i == 0 || (opt_finetune_reset && opt_burnin >= 200 && i < 0 &&
                   ft_round >= 100 && i%(opt_burnin/4)==0))
    {
      int pjump_size = PROP_COUNT + 1+1 + GTR_PROP_COUNT + CLOCK_PROP_COUNT + !!opt_migration;

      if (opt_finetune_reset && opt_burnin >= 200)
      {
        /* If last MCMC status line did not end with a newline, print one */
        if (!print_newline)
          fprintf(stdout, "\n");
        reset_finetune(fp_out, pjump);
      }

      /* reset pjump and number of steps since last finetune reset to zero */
      ft_round = 0;
      memset(pjump, 0, pjump_size * sizeof(double));

      if (opt_est_delimit)
      {
        pjump_rj = 0;
        ft_round_rj = 0;
      }
      if (opt_method == METHOD_10)      /* species delimitation */
        memset(posterior,0,delimitation_getparam_count()*sizeof(double));

      if (opt_migration)
        for (j = 0; j < mean_mrate_count; ++j)
          mean_mrate_round[j] = 0;

      if (opt_est_stree)
      {
        ft_round_spr = ft_round_snl = 0;
        pjump_spr = pjump_snl = 0;
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
        long stree_snl = 0;
        if (opt_prob_snl == 0)      stree_snl = 0;
        else if (opt_prob_snl == 1) stree_snl = 1;
        else                        stree_snl = (legacy_rndu(thread_index_zero) < opt_prob_snl);
        if (stree_snl==0) {
          ret = stree_propose_spr(&stree, &gtree, &sclone, &gclones, locus);
          ft_round_spr++;
          if (ret == 1) pjump_spr++;
        }
        else {
          ret = stree_propose_stree_snl(&stree, &gtree, &sclone, &gclones, locus);
          ft_round_snl++;
          if (ret==1) pjump_snl++;
        }
        if (ret == 1)
        {
          /* accepted */
          /* swap the pointers of species tree and gene tree list with cloned */
          SWAP(stree,sclone);
          SWAP(gtree,gclones);
          stree_label(stree);
        }
        if (opt_debug_bruce)
          debug_bruce(stree,gtree,stree_snl == 0 ? "SSPR" : "SNL", i, fp_debug);
      }
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "SSPR");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "SSPR");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "SSPR");
      #endif
    }

    /* perform proposals sequentially */   

      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "GAGE-1");
      #endif

    /* propose gene tree ages */
    /* Note: call serial version when thetas are integrated out */
    if (!opt_est_theta || opt_threads == 1)
      ratio = gtree_propose_ages_serial(locus, gtree, stree);
    else
    {
      td.locus = locus; td.gtree = gtree; td.stree = stree;
      threads_wakeup(THREAD_WORK_GTAGE,&td);
      ratio = td.accepted ? ((double)(td.accepted)/td.proposals) : 0;
    }
    pjump[BPP_MOVE_GTAGE_INDEX] = (pjump[BPP_MOVE_GTAGE_INDEX]*(ft_round-1)+ratio) /
                                  (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "GAGE");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "GAGE");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "GAGE");
      #endif
      if (opt_debug_bruce)
        debug_bruce(stree,gtree,"GAGE", i, fp_debug);

    /* propose migration ages */
    if (opt_migration)
    {
      ratio = gtree_propose_migevent_ages_serial(locus,gtree,stree);
    }

    /* propose gene tree topologies using SPR */
    /* Note: call serial version when thetas are integrated out */
    if (!opt_est_theta || opt_threads == 1)
      ratio = gtree_propose_spr_serial(locus,gtree,stree);
    else
    {
      td.locus = locus; td.gtree = gtree; td.stree = stree;
      threads_wakeup(THREAD_WORK_GTSPR,&td);
      ratio = td.accepted ? ((double)(td.accepted)/td.proposals) : 0;
    }
    pjump[BPP_MOVE_GTSPR_INDEX] = (pjump[BPP_MOVE_GTSPR_INDEX]*(ft_round-1)+ratio) /
                                  (double)ft_round;

      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "GSPR");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "GSPR");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "GSPR");
      #endif
      if (opt_debug_bruce)
        debug_bruce(stree,gtree,"GSPR", i, fp_debug);

    /* propose population sizes on species tree */
    if (opt_est_theta)
    {
      //ratio = stree_propose_theta(gtree,locus,stree);
      //pjump[BPP_MOVE_THETA_INDEX] = (pjump[BPP_MOVE_THETA_INDEX]*(ft_round-1)+ratio) /
                                    (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "THETA");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "THETA");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "THETA");
      #endif
    }

    /* propose species tree taus */
   /* if (stree->tip_count > 1 && stree->root->tau > 0)
    {
      if (opt_migration)
        ratio = stree_propose_tau_mig(&stree, &gtree, &sclone, &gclones, locus);
      else
        ratio = stree_propose_tau(gtree,stree,locus);
      pjump[BPP_MOVE_TAU_INDEX] = (pjump[BPP_MOVE_TAU_INDEX]*(ft_round-1)+ratio) /
                                  (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "TAU");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "TAU");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "TAU");
      #endif
    } */

    /* propose migration rates */
    if (opt_migration)
    {
      ratio = prop_migrates(stree,gtree,locus);
      pjump[BPP_MOVE_MRATE_INDEX] = (pjump[BPP_MOVE_MRATE_INDEX]*(ft_round-1)+ratio) /
                                  (double)ft_round;

    }

    /* mixing step */
    /*ratio = proposal_mixing(gtree, stree, locus);
    pjump[BPP_MOVE_MIX_INDEX] = (pjump[BPP_MOVE_MIX_INDEX] * (ft_round - 1) + ratio) /
                                (double)ft_round;
				*/
    #ifdef CHECK_LOGL
    check_logl(stree, gtree, locus, i, "MIXING");
    #endif
    #ifdef CHECK_LOGPR
    debug_validate_logpg(stree, gtree, locus, "MIXING");
    #endif
    #ifdef CHECK_LNPRIOR
    check_lnprior(stree, gtree, i, "MIXING");
    #endif

    if ((opt_est_locusrate == MUTRATE_ESTIMATE &&
         opt_locusrate_prior == BPP_LOCRATE_PRIOR_DIR) ||
         opt_est_heredity == HEREDITY_ESTIMATE)
    {
      ratio = prop_locusrate_and_heredity(gtree,stree,locus,thread_index_zero);
      pjump[BPP_MOVE_LRHT_INDEX] = (pjump[BPP_MOVE_LRHT_INDEX]*(ft_round-1)+ratio) /
                                   (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "LRHT");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "LRHT");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "LRHT");
      #endif
    }

    /* phi proposal */
    if (opt_msci)
    {
      ratio = stree_propose_phi(stree,gtree);
      pjump[BPP_MOVE_PHI_INDEX] = (pjump[BPP_MOVE_PHI_INDEX]*(ft_round-1)+ratio) /
                                  (double)ft_round;
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "PHI");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "PHI");
      #endif
    }

    if (enabled_prop_freqs)
    {
      if (opt_threads == 1)
        ratio = locus_propose_freqs_serial(stree,locus,gtree);
      else
      {
        td.locus = locus; td.gtree = gtree;
        threads_wakeup(THREAD_WORK_FREQS,&td);
        ratio = td.proposals ? ((double)(td.accepted)/td.proposals) : 0;
      }
      pjump[BPP_MOVE_FREQS_INDEX] = (pjump[BPP_MOVE_FREQS_INDEX]*(ft_round-1)+ratio) /
                                    (double)ft_round;
    }

    if (enabled_prop_qrates)
    {
      if (opt_threads == 1)
        ratio = locus_propose_qrates_serial(stree,locus,gtree);
      else
      {
        td.locus = locus; td.gtree = gtree;
        threads_wakeup(THREAD_WORK_RATES,&td);
        ratio = td.proposals ? ((double)(td.accepted)/td.proposals) : 0;
      }
      pjump[BPP_MOVE_QRATES_INDEX] = (pjump[BPP_MOVE_QRATES_INDEX]*(ft_round-1)+ratio) /
                                    (double)ft_round;
    }

    if (enabled_prop_alpha)
    {
      if (opt_threads == 1)
        ratio = locus_propose_alpha_serial(stree,locus,gtree);
      else
      {
        td.locus = locus; td.gtree = gtree;
        threads_wakeup(THREAD_WORK_ALPHA,&td);
        ratio = td.accepted ? ((double)(td.accepted)/td.proposals) : 0;
      }
      pjump[BPP_MOVE_ALPHA_INDEX] = (pjump[BPP_MOVE_ALPHA_INDEX]*(ft_round-1)+ratio) /
                                    (double)ft_round;
    }

    /* TODO: Delete after debugging */
    if (opt_debug_rates && i==0)
    {
      opt_seed = 1;
      legacy_init();
    }

    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL ||
         opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR))
    {
      ratio = prop_locusrate_mui(gtree,stree,locus,thread_index_zero);
      pjump[BPP_MOVE_MUI_INDEX] = (pjump[BPP_MOVE_MUI_INDEX]*(ft_round-1)+ratio) /
                                  (double)ft_round;

      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "MUI");
      #endif
      #ifdef CHECK_LOGPR
      debug_validate_logpg(stree, gtree, locus, "MUI");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "MUI");
      #endif

      if (opt_est_mubar)
      {
        ratio = prop_locusrate_mubar(stree,gtree);
        pjump[BPP_MOVE_MUBAR_INDEX] = (pjump[BPP_MOVE_MUBAR_INDEX]*(ft_round-1)+ratio) /
                                       (double)ft_round;
        #ifdef CHECK_LOGL
        check_logl(stree, gtree, locus, i, "MUBAR");
        #endif
        #ifdef CHECK_LNPRIOR
        check_lnprior(stree, gtree, i, "MUBAR");
        #endif
      }
    }

    if (opt_clock != BPP_CLOCK_GLOBAL)
    {

      ratio = prop_locusrate_nui(gtree,stree,locus,thread_index_zero);
      pjump[BPP_MOVE_NUI_INDEX] = (pjump[BPP_MOVE_NUI_INDEX]*(ft_round-1)+ratio) /
                                        (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "NUI");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "NUI");
      #endif

      if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
      {
        ratio = prop_locusrate_nubar(stree,gtree);
        pjump[BPP_MOVE_NUBAR_INDEX] = (pjump[BPP_MOVE_NUBAR_INDEX]*(ft_round-1)+ratio) /
                                          (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "NUBAR");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "NUBAR");
      #endif
      }

      if (opt_threads == 1)
        ratio = prop_branch_rates_serial(gtree,stree,locus);
      else
      {
        td.locus = locus; td.gtree = gtree;
        threads_wakeup(THREAD_WORK_BRATE,&td);
        ratio = td.proposals ? ((double)(td.accepted)/td.proposals) : 0;
      }
      pjump[BPP_MOVE_BRANCHRATE_INDEX] = (pjump[BPP_MOVE_BRANCHRATE_INDEX]*(ft_round-1)+ratio) /
                                         (double)ft_round;
      #ifdef CHECK_LOGL
      check_logl(stree, gtree, locus, i, "BRATE");
      #endif
      #ifdef CHECK_LNPRIOR
      check_lnprior(stree, gtree, i, "BRATE");
      #endif
    }

    /* log sample into file (dparam_count is only used in method 10) */
    if ((i + 1) % (opt_samplefreq*5) == 0)
       fflush(NULL);
    if (i >= 0 && (i+1)%opt_samplefreq == 0)
    {
      mcmc_logsample(fp_mcmc,i+1,stree,gtree,locus,dparam_count,ndspecies);

      /* log gene trees */
      if (opt_print_genetrees)
        print_gtree(fp_gtree,stree,gtree);

      /* log rates */
      if (opt_print_locusfile)
        print_rates(fp_locus,stree,gtree,locus);
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
      
    if (opt_migration)
    {
      for (j = 0; j < mean_mrate_count; ++j)
      {
        double mv = opt_migration_matrix[mean_mrate_row[j]][mean_mrate_col[j]];
        if (migration_valid(stree,
                            stree->nodes[mean_mrate_row[j]],
                            stree->nodes[mean_mrate_col[j]]))
        {
          ++mean_mrate_round[j];
          mean_mrate[j] = (mean_mrate[j]*(mean_mrate_round[j]-1) + mv) / mean_mrate_round[j];
        }
      }
    }
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
      max_param_count = MIN(stree->tip_count+stree->inner_count,MAX_TAU_OUTPUT);

      /* 2. calculate means */
      k = 0;
      for (j = stree->tip_count; j < stree->tip_count+stree->inner_count; ++j)
      {
        if (stree->nodes[j]->tau == 0) continue;
        
        mean_tau[k] = (mean_tau[k]*(ft_round-1) + stree->nodes[j]->tau)/ft_round;

        if (++k == max_param_count) break;
      }
      mean_tau_count = k;
      if (opt_msci)
      {
        mean_phi_count = MIN(stree->hybrid_count, MAX_PHI_OUTPUT);

        for (j = 0; j < mean_phi_count; ++j)
        {
          snode_t * tmpnode = stree->nodes[stree->tip_count+stree->inner_count+j];
          if (!tmpnode->has_phi)
            tmpnode = tmpnode->hybrid;

          mean_phi[j] = (mean_phi[j]*(ft_round-1) + tmpnode->hphi)/ft_round;
        }
      }

      #ifdef DEBUG_GTR
      if (fabs(locus[0]->frequencies[0][0] + locus[0]->frequencies[0][1] + locus[0]->frequencies[0][2] + locus[0]->frequencies[0][3]-1) >= 1e-10)
      {
        printf ("%f %f %f %f = %f\n", locus[0]->frequencies[0][0], locus[0]->frequencies[0][1], locus[0]->frequencies[0][2], locus[0]->frequencies[0][3], 
                                 locus[0]->frequencies[0][0]+ locus[0]->frequencies[0][1]+ locus[0]->frequencies[0][2]+ locus[0]->frequencies[0][3]);
        assert(0);
      }

      mean_freqa = (mean_freqa*(ft_round-1) + locus[0]->frequencies[0][0])/ft_round;
      mean_freqc = (mean_freqc*(ft_round-1) + locus[0]->frequencies[0][1])/ft_round;
      mean_freqg = (mean_freqg*(ft_round-1) + locus[0]->frequencies[0][2])/ft_round;
      mean_freqt = (mean_freqt*(ft_round-1) + locus[0]->frequencies[0][3])/ft_round;

      mean_ratea = (mean_ratea*(ft_round-1) + locus[0]->subst_params[0][0])/ft_round;
      mean_rateb = (mean_rateb*(ft_round-1) + locus[0]->subst_params[0][1])/ft_round;
      mean_ratec = (mean_ratec*(ft_round-1) + locus[0]->subst_params[0][2])/ft_round;
      mean_rated = (mean_rated*(ft_round-1) + locus[0]->subst_params[0][3])/ft_round;
      mean_ratee = (mean_ratee*(ft_round-1) + locus[0]->subst_params[0][4])/ft_round;
      mean_ratef = (mean_ratef*(ft_round-1) + locus[0]->subst_params[0][5])/ft_round;

      mean_alpha0 = (mean_alpha0*(ft_round-1) + locus[0]->rates_alpha)/ft_round;
      #endif
    }
    dbg_mean_rate = (dbg_mean_rate*(ft_round-1) + gtree[0]->rate_mui)/ft_round;

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
      print_newline = (printk >= 50 && (i+1) % (printk / 20) == 0);

      /* print progress percentage */
      printf("\r%4.0f%% ", (i + 1.499) / printk * 100.);

      double mean_pjump_rj = 0;
      if (opt_method == METHOD_10)
        mean_pjump_rj = ft_round_rj ? pjump_rj / ft_round_rj : 0;

      /* print pjumps */
      status_print_pjump(stdout, pjump, ft_round_spr, ft_round_snl, pjump_spr, pjump_snl, mean_pjump_rj);
      if (print_newline)
      {
        fprintf(fp_out, "%4.0f%% ", (i + 1.499) / printk * 100.);
        status_print_pjump(fp_out,pjump, ft_round_spr, ft_round_snl, pjump_spr, pjump_snl, mean_pjump_rj);
      }

      /* species delimitation specific output */
      if (opt_method == METHOD_10)
      {
        long bmodel = 0;
        for (j = 1; j < delimitation_getparam_count(); ++j)
          if (posterior[j] > posterior[bmodel]) bmodel = j;

        
        /* Number of parameters, pjump/finetune round ratio, current
           delimitation binary string, and support value of 'most-supported
           delimitation delimitation' */
         printf(" %2ld %s", dparam_count,
                            delimitation_getparam_string());
         char * stmp = NULL;
         xasprintf(&stmp,"P[%ld]=%6.4f", bmodel+1, posterior[bmodel]/ft_round);
         printf(" %*s",4+6+delim_digit_count, stmp);
         if (print_newline)
         {
           fprintf(fp_out, " %2ld %s", dparam_count, delimitation_getparam_string());
           fprintf(fp_out, " %*s", 4+6+delim_digit_count, stmp);
         }
         free(stmp);

      }

      if (opt_method == METHOD_11)
      {
        long ndspeciesbest = 0;
        for (j = 1; j < stree->tip_count; ++j)
          if (pspecies[j] > pspecies[ndspeciesbest]) ndspeciesbest = j;
       
        char * stmp = NULL;
        xasprintf(&stmp,"P[%ld]=%6.4f", ndspeciesbest+1,pspecies[ndspeciesbest]);
        printf(" %2ld %2ld", ndspecies, dparam_count);
        printf(" %*s", 4+6+delim_digit_count,stmp);
        if (print_newline)
        {
          fprintf(fp_out, " %2ld %2ld", ndspecies, dparam_count);
          fprintf(fp_out, " %*s", 4+6+delim_digit_count,stmp);
        }
        free(stmp);
      }

      if (opt_est_theta)
      {
        for (j = 0; j < mean_theta_count; ++j)
          printf(" %6.4f", mean_theta[j]);
        printf(" ");
        if (print_newline)
        {
          for (j = 0; j < mean_theta_count; ++j)
            fprintf(fp_out, " %6.4f", mean_theta[j]);
          fprintf(fp_out, " ");
        }
      }

      for (j = 0; j < mean_tau_count; ++j)
        printf(" %6.4f", mean_tau[j]);
      printf(" ");
      if (print_newline)
      {
        for (j = 0; j < mean_tau_count; ++j)
          fprintf(fp_out, " %6.4f", mean_tau[j]);
        fprintf(fp_out, " ");
      }

      if (opt_migration)
      {
        for (j = 0; j < mean_mrate_count; ++j)
          printf(" %6.4f", mean_mrate[j]);
        printf(" ");
        if (print_newline)
        {
          for (j = 0; j < mean_mrate_count; ++j)
            fprintf(fp_out, " %6.4f", mean_mrate[j]);
          fprintf(fp_out, " ");
        }
      }

      if (opt_msci)
      {
        for (j = 0; j < mean_phi_count; ++j)
          printf(" %6.4f", mean_phi[j]);
        printf(" ");
        if (print_newline)
        {
          for (j = 0; j < mean_phi_count; ++j)
            fprintf(fp_out, " %6.4f", mean_phi[j]);
          fprintf(fp_out, " ");
        }
      }

      #if 0
      if (opt_est_locusrate == MUTRATE_ESTIMATE)
        printf(" MU_0=%6.4f ", dbg_mean_rate);
      #endif

      #ifdef DEBUG_GTR
        printf(" ( %6.4f %6.4f %6.4f %6.4f ) ", mean_freqa,mean_freqc,mean_freqg,mean_freqt);
        printf(" ( %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ) ", mean_ratea,mean_rateb,mean_ratec,mean_rated,mean_ratee,mean_ratef);
        printf(" ( %6.4f ) ", mean_alpha0);
      #endif

      double logpr_sum = 0;
      if (opt_est_theta)
        for (j = 0; j < opt_locus_count; ++j)
          logpr_sum += gtree[j]->logpr;
      else
        logpr_sum = stree->notheta_logpr;
      printf(" %*.5f", prec_logpr, logpr_sum);
      if (print_newline)
        fprintf(fp_out, " %*.5f", prec_logpr, logpr_sum);

      if (opt_usedata)
      {
        printf(" %*.5f", prec_logl, mean_logl);
        if (print_newline)
          fprintf(fp_out, " %*.5f", prec_logl, mean_logl);
      }

      if (print_newline)
      {
        timer_print("  ","\n", fp_out);

        if (isinf(mean_logl))
          fatal("\n[ERROR] The mean log-L over loci is -inf.\n"
                "Please run BPP with numerical scaling. This is enabled by adding the line:\n"
                "\n  scaling = 1\n\nto the control file");
      }
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

        /* if relaxed clock is enabled get offsets for rates files */
        if (opt_print_locusfile)
          for (j = 0; j < opt_locus_count; ++j)
            rates_offset[j] = ftell(fp_locus[j]);

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
                        rates_offset,
                        dparam_count,
                        posterior,
                        pspecies,
                        opt_est_delimit ?
                          delimitation_getparam_count() : 0,
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
                        prec_logpr,
                        prec_logl);
      }
    }
    if (opt_debug_abort == opt_debug_counter)
      fatal("[DBG] Aborting debugging (reached step %ld)", opt_debug_abort);
  }
  if (!opt_onlysummary)
    timer_print("\n", " spent in MCMC\n\n", fp_out);

  #if 0
  progress_done();
  #endif

  free(pjump);

  if (opt_threads > 1)
    threads_exit();

  if (opt_bfbeta != 1 && !opt_onlysummary)
  {
    fprintf(stdout, "\nBFbeta = %8.6f  E_b(lnf(X)) = %9.4f\n\n", opt_bfbeta, mean_logl);
    fprintf(fp_out, "\nBFbeta = %8.6f  E_b(lnf(X)) = %9.4f\n\n", opt_bfbeta, mean_logl);
  }

  /* close mcmc file */
  if (!opt_onlysummary)
    fclose(fp_mcmc);

  /* close files containing rates sample for each locus */
  if (opt_print_locusfile)
  {
    for (i = 0; i < opt_locus_count; ++i)
      fclose(fp_locus[i]);
    free(fp_locus);
    if (opt_checkpoint)
      free(rates_offset);
  }

  /* TODO: Delete after debugging */
  if (opt_debug_rates)
  {
    char * s;
    char * cmd;

    /* execute ds on rates sample of first locus */
    xasprintf(&s, template_ratesfile, 1);
    xasprintf(&cmd, "/a/c/ds %s", s);
    system(cmd);
    free(s);
    free(cmd);

    if (opt_locus_count > 1)
    {
      /* execute ds on rates sample of last locus */
      xasprintf(&s, template_ratesfile, opt_locus_count-1);
      xasprintf(&cmd, "/a/c/ds %s", s);
      system(cmd);
      free(s);
      free(cmd);
    }
  }

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
    mixed_summary(fp_out,stree->tip_count);
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
  if (opt_est_stree || opt_migration)          /* species tree inference */
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
  if (opt_migration)
  {
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    {
      free(opt_mig_bitmatrix[i]);
      free(opt_migration_matrix[i]);
    }
    free(opt_mig_bitmatrix);
    free(opt_migration_matrix);
  }


  /* deallocate tree */
  stree_destroy(stree,NULL);
  if (opt_est_stree || opt_migration)
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

    #if 0
    if (opt_prob_snl && !opt_onlysummary)
    {
      long opt_debug_sum = opt_debug_expand_count +
                           opt_debug_expshr_count +
                           opt_debug_shrink_count;
      printf("[DEBUG] SHRINK: %ld (%.2f%%) EXPAND: %ld (%.2f%%) EXPSHR: %ld (%.2f%%)\n",
             opt_debug_shrink_count, (opt_debug_shrink_count / (double)opt_debug_sum)*100,
             opt_debug_expand_count, (opt_debug_expand_count / (double)opt_debug_sum)*100,
             opt_debug_expshr_count, (opt_debug_expshr_count / (double)opt_debug_sum)*100);
    }
    #endif
  }

  if (opt_est_theta)
    free(mean_theta);
  free(mean_tau);
  if (opt_msci)
    free(mean_phi);

  if (opt_migration)
  {
    free(mean_mrate);
    free(mean_mrate_row);
    free(mean_mrate_col);
    free(mean_mrate_round);
  }

  if (opt_diploid)
    free(opt_diploid);

  if (opt_partition_file)
    free(opt_partition_file);

  /* close output file */
  fclose(fp_out);
}
