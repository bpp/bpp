/*
    Copyright (C) 2016-2018 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#ifdef _MSC_VER
#include "getopt_win.h"
#endif

static char * progname;
static char progheader[80];
char * cmdline;

/* global error message buffer */
__THREAD int bpp_errno;
__THREAD char bpp_errmsg[200] = {0};

/* options */
long opt_arch;
long opt_burnin;
long opt_checkpoint;
long opt_checkpoint_current;
long opt_checkpoint_initial;
long opt_checkpoint_step;
long opt_cleandata;
long opt_debug;
long opt_delimit_prior;
long opt_diploid_size;
long opt_est_delimit;
long opt_est_heredity;
long opt_est_locusrate;
long opt_est_stree;
long opt_est_theta;
long opt_experimental_method;
long opt_experimental_debug;
long opt_finetune_reset;
long opt_help;
long opt_max_species_count;
long opt_method;
long opt_locus_count;
long opt_print_genetrees;
long opt_print_hscalars;
long opt_print_locusrate;
long opt_print_samples;
long opt_quiet;
long opt_rjmcmc_method;
long opt_samplefreq;
long opt_samples;
long opt_seed;
long opt_usedata;
long opt_version;
double opt_bfbeta;
double opt_finetune_gtage;
double opt_finetune_gtspr;
double opt_finetune_locusrate;
double opt_finetune_mix;
double opt_finetune_tau;
double opt_finetune_theta;
double opt_heredity_alpha;
double opt_heredity_beta;
double opt_locusrate_alpha;
double opt_rjmcmc_alpha;
double opt_rjmcmc_epsilon;
double opt_rjmcmc_mean;
double opt_tau_alpha;
double opt_tau_beta;
double opt_theta_alpha;
double opt_theta_beta;
long * opt_diploid;
long * opt_sp_seqcount;
char * opt_cfile;
char * opt_heredity_filename;
char * opt_locusrate_filename;
char * opt_mapfile;
char * opt_msafile;
char * opt_mcmcfile;
char * opt_outfile;
char * opt_reorder;
char * opt_resume;
char * opt_streenewick;

long mmx_present;
long sse_present;
long sse2_present;
long sse3_present;
long ssse3_present;
long sse41_present;
long sse42_present;
long popcnt_present;
long avx_present;
long avx2_present;
long altivec_present;

static struct option long_options[] =
{
  {"help",       no_argument,       0, 0 },  /* 0 */
  {"version",    no_argument,       0, 0 },  /* 1 */
  {"quiet",      no_argument,       0, 0 },  /* 2 */
  {"cfile",      required_argument, 0, 0 },  /* 3 */
  {"arch",       required_argument, 0, 0 },  /* 4 */
  {"exp_method", required_argument, 0, 0 },  /* 5 */
  {"exp_debug",  no_argument,       0, 0 },  /* 6 */
  {"resume",     required_argument, 0, 0 },  /* 7 */
  { 0, 0, 0, 0 }
};

static long args_getlong(char * arg)
{
  int len = 0;
  long temp;

  int ret = sscanf(arg, "%ld%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

#if 0
static double args_getdouble(char * arg)
{
  int len = 0;
  double temp = 0;
  int ret = sscanf(arg, "%lf%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}
#endif

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;

  /* set defaults */

  progname = argv[0];

  opt_arch = -1;
  opt_bfbeta = 1;
  opt_burnin = 100;
  opt_cfile = NULL;
  opt_checkpoint = 0;
  opt_checkpoint_initial = 0;
  opt_checkpoint_current = 0;
  opt_checkpoint_step = 0;
  opt_cleandata = 0;
  opt_debug = 0;
  opt_est_delimit = 0;
  opt_delimit_prior = BPP_SPECIES_PRIOR_UNIFORM;
  opt_diploid = NULL;
  opt_diploid_size = 0;
  opt_est_heredity = 0;
  opt_est_locusrate = 0;
  opt_est_theta = 1;
  opt_experimental_debug = 0;
  opt_experimental_method = 0;
  opt_finetune_gtage = 5;
  opt_finetune_gtspr = 0.001;
  opt_finetune_mix   = 0.3;
  opt_finetune_locusrate = 0.33;
  opt_finetune_reset = 0;
  opt_finetune_tau   = 0.001;
  opt_finetune_theta = 0.001;
  opt_help = 0;
  opt_heredity_alpha = 0;
  opt_heredity_beta = 0;
  opt_heredity_filename = NULL;
  opt_locusrate_alpha = 0;
  opt_locusrate_filename = NULL;
  opt_locus_count = 0;
  opt_mapfile = NULL;
  opt_mcmcfile = NULL;
  opt_method = -1;
  opt_msafile = NULL;
  opt_outfile = NULL;
  opt_print_genetrees = 0;
  opt_print_hscalars = 0;
  opt_print_locusrate = 0;
  opt_print_samples = 1;
  opt_quiet = 0;
  opt_resume = NULL;
  opt_rjmcmc_alpha = -1;
  opt_rjmcmc_epsilon = -1;
  opt_rjmcmc_mean = -1;
  opt_rjmcmc_method = -1;
  opt_samplefreq = 10;
  opt_samples = 0;
  opt_seed = (long)time(NULL);
  opt_est_stree = 0;
  opt_streenewick = NULL;
  opt_tau_alpha = 0;
  opt_tau_beta = 0;
  opt_theta_alpha = 0;
  opt_theta_beta = 0;
  opt_usedata = 1;
  opt_version = 0;
  opt_max_species_count = 0;
  opt_sp_seqcount = NULL;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        opt_quiet = 1;
        break;

      case 3:
        opt_cfile = xstrdup(optarg);
        break;

      case 4:
        if (!strcasecmp(optarg,"cpu"))
          opt_arch = PLL_ATTRIB_ARCH_CPU;
        else if (!strcasecmp(optarg,"sse"))
          opt_arch = PLL_ATTRIB_ARCH_SSE;
        else if (!strcasecmp(optarg,"avx"))
          opt_arch = PLL_ATTRIB_ARCH_AVX;
        else if (!strcasecmp(optarg,"avx2"))
          opt_arch = PLL_ATTRIB_ARCH_AVX2;
        else
          fatal("Invalid instruction set (%s)", optarg);
        break;

      case 5:
        opt_experimental_method = args_getlong(optarg);
        break;

      case 6:
        opt_experimental_debug = 1;
        break;

      case 7:
        opt_resume = optarg;
        break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  if (opt_cfile)
    load_cfile();

  /* check for number of independent commands selected */
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
  if (opt_cfile)
    commands++;
  if (opt_resume)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  /* if no command specified, turn on --help */
  if (!commands)
  {
    opt_help = 1;
    return;
  }
}

static void dealloc_switches()
{
  if (opt_cfile) free(opt_cfile);
  if (opt_mapfile) free(opt_mapfile);
  if (opt_mcmcfile) free(opt_mcmcfile);
  if (opt_msafile) free(opt_msafile);
  if (opt_outfile) free(opt_outfile);
  if (opt_reorder) free(opt_reorder);
  if (opt_sp_seqcount) free(opt_sp_seqcount);
  if (opt_streenewick) free(opt_streenewick);
}

void cmd_help()
{
  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help             display help information\n"
          "  --version          display version information\n"
          "  --quiet            only output warnings and fatal errors to stderr\n"
          "  --cfile FILENAME   run analysis for the specified control file\n"
          "  --resume FILENAME  resume analysis from a specified checkpoint file\n"
          "  --arch SIMD        force specific vector instruction set (default: auto)\n"
          "\n"
         );

  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           arch_get_cores());
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/xflouris/bpp\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  cpu_features_detect();
  cpu_features_show();
  if (!opt_version && !opt_help)
    cpu_setarch();

  /* intiialize random number generators */
  legacy_init();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_version)
  {
    ;
  }
  else if (opt_resume || opt_cfile)
  {
    cmd_run();
  }

  dealloc_switches();
  free(cmdline);
  return (0);
}
