/*
    Copyright (C) 2016 Tomas Flouri and Ziheng Yang

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

#ifdef _MSC_VER
#include "getopt.h"
#endif

static char * progname;
static char progheader[80];
char * cmdline;

/* global error message buffer */
__THREAD int bpp_errno;
__THREAD char bpp_errmsg[200] = {0};

/* number of mandatory options for the user to input */
static const char mandatory_options_count = 3;
static const char * mandatory_options_list = " --stree_file --output_file --mcmc_file";

/* options */
long opt_help;
long opt_version;
long opt_quiet;
long opt_seed;
long opt_stree;
long opt_arch;
long opt_delimit;
long opt_delimit_prior;
long opt_cleandata;
long opt_debug;
long opt_est_theta;
long opt_samples;
long opt_samplefreq;
long opt_burnin;
long opt_finetune_reset;
long opt_rjmcmc_method;
long opt_usedata;
double opt_tau_alpha;
double opt_tau_beta;
double opt_theta_alpha;
double opt_theta_beta;
double opt_finetune_gtage;
double opt_finetune_gtspr;
double opt_finetune_theta;
double opt_finetune_tau;
double opt_finetune_mix;
double opt_rjmcmc_alpha;
double opt_rjmcmc_mean;
double opt_rjmcmc_epsilon;
char * opt_cfile;
char * opt_mapfile;
char * opt_msafile;
char * opt_mapfile;
char * opt_mcmcfile;
char * opt_outfile;
char * opt_reorder;
char * opt_streefile;
char * opt_streenewick;

long opt_debugflag = 0;

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
  {"help",            no_argument,       0, 0 },  /*  0 */
  {"version",         no_argument,       0, 0 },  /*  1 */
  {"quiet",           no_argument,       0, 0 },  /*  2 */
  {"stree_file",      required_argument, 0, 0 },  /*  3 */
  {"map_file",        required_argument, 0, 0 },  /*  4 */
  {"output_file",     required_argument, 0, 0 },  /*  5 */
  {"msa_file",        required_argument, 0, 0 },  /*  6 */
  {"seed",            required_argument, 0, 0 },  /*  7 */
  {"stree",           no_argument,       0, 0 },  /*  8 */
  {"delimit",         no_argument,       0, 0 },  /*  9 */
  {"tauprior",        required_argument, 0, 0 },  /* 10 */
  {"thetaprior",      required_argument, 0, 0 },  /* 11 */
  {"cleandata",       no_argument,       0, 0 },  /* 12 */
  {"debug",           no_argument,       0, 0 },  /* 13 */
  {"samples",         required_argument, 0, 0 },  /* 14 */
  {"samplefreq",      required_argument, 0, 0 },  /* 15 */
  {"burnin",          required_argument, 0, 0 },  /* 16 */
  {"finetune_reset",  no_argument,       0, 0 },  /* 17 */
  {"finetune_params", required_argument, 0, 0 },  /* 18 */
  {"mcmc_file",       required_argument, 0, 0 },  /* 19 */
  {"reorder",         required_argument, 0, 0 },  /* 20 */
  {"delimit_prior",   required_argument, 0, 0 },  /* 21 */
  {"rjmcmc_alpha",    required_argument, 0, 0 },  /* 22 */
  {"rjmcmc_mean",     required_argument, 0, 0 },  /* 23 */
  {"rjmcmc_epsilon",  required_argument, 0, 0 },  /* 24 */
  {"cfile",           required_argument, 0, 0 },  /* 25 */
  {"nodata",          no_argument,       0, 0 },  /* 26 */
  {"arch",            required_argument, 0, 0 },  /* 27 */
  { 0, 0, 0, 0 }
};

static int args_getgammaprior(char * arg, double * a, double * b)
{
  int len = 0;

  int ret = sscanf(arg, "%lf,%lf%n", a, b, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    return 0;
  return 1;
}

static int args_getftparams(char * arg)
{
  int len = 0;

  int ret = sscanf(arg,
                   "%lf,%lf,%lf,%lf,%lf%n", 
                   &opt_finetune_gtage,
                   &opt_finetune_gtspr,
                   &opt_finetune_theta,
                   &opt_finetune_tau,
                   &opt_finetune_mix,
                   &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    return 0;
  return 1;
}

static long args_getlong(char * arg)
{
  int len = 0;
  long temp;

  int ret = sscanf(arg, "%ld%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

static double args_getdouble(char * arg)
{
  int len = 0;
  double temp = 0;
  int ret = sscanf(arg, "%lf%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  int mand_options = 0;

  /* set defaults */

  progname = argv[0];

  opt_burnin = 100;
  opt_cleandata = 0;
  opt_debug = 0;
  opt_est_theta = 1;
  opt_delimit = 0;
  opt_delimit_prior = BPP_DELIMIT_PRIOR_UNIFORM;
  opt_finetune_gtage = 5;
  opt_finetune_gtspr = 0.001;
  opt_finetune_mix   = 0.3;
  opt_finetune_reset = 0;
  opt_finetune_tau   = 0.001;
  opt_finetune_theta = 0.001;
  opt_help = 0;
  opt_quiet = 0;
  opt_samplefreq = 10;
  opt_samples = 10000;
  opt_stree = 0;
  opt_tau_alpha = 0;
  opt_tau_beta = 0;
  opt_theta_alpha = 0;
  opt_theta_beta = 0;
  opt_version = 0;
  opt_cfile = NULL;
  opt_msafile = NULL;
  opt_mapfile = NULL;
  opt_mcmcfile = NULL;
  opt_outfile = NULL;
  opt_seed = (long)time(NULL);
  opt_streefile = NULL;
  opt_streenewick = NULL;
  opt_rjmcmc_alpha = -1;
  opt_rjmcmc_mean = -1;
  opt_rjmcmc_epsilon = -1;
  opt_rjmcmc_method = -1;
  opt_usedata = 1;
  opt_arch = PLL_ATTRIB_ARCH_CPU;

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
        opt_streefile = xstrdup(optarg);
        break;

      case 4:
        opt_mapfile = xstrdup(optarg);
        break;

      case 5:
        opt_outfile = xstrdup(optarg);
        break;

      case 6:
        opt_msafile = xstrdup(optarg);
        break;

      case 7:
        opt_seed = args_getlong(optarg);
        break;

      case 8:
        opt_stree = 1;
        break;

      case 9:
        opt_delimit = 1;
        break;

      case 10:
        if (!args_getgammaprior(optarg, &opt_tau_alpha, &opt_tau_beta))
          fatal("Illegal format for --tauprior");
        break;

      case 11:
        if (!args_getgammaprior(optarg, &opt_theta_alpha, &opt_theta_beta))
          fatal("Illegal format for --thetaprior");
        break;

      case 12:
        opt_cleandata = 1;
        break;

      case 13:
        opt_debug = 1;
        break;

      case 14:
        opt_samples = args_getlong(optarg);
        break;

      case 15:
        opt_samplefreq = args_getlong(optarg);
        break;

      case 16:
        opt_burnin = args_getlong(optarg);
        break;

      case 17:
        opt_finetune_reset = 1;
        break;

      case 18:
        if (!args_getftparams(optarg))
          fatal("Illegal format for --finetune_params");
        break;

      case 19:
        opt_mcmcfile = xstrdup(optarg);
        break;

      case 20:
        opt_reorder = xstrdup(optarg);
        break;

      case 21:
        if (!strcmp(optarg,"dirichlet"))
          opt_delimit_prior = BPP_DELIMIT_PRIOR_DIRICHLET;
        else if (!strcmp(optarg,"uniform"))
          opt_delimit_prior = BPP_DELIMIT_PRIOR_UNIFORM;
        else
          fatal("Unknown species delimitation prior: %s", optarg);
        break;

      case 22:
        if (opt_rjmcmc_method == 0)
          fatal("Cannot use --rjmcmc_alpha with --rjmcmc_epsilon");
        opt_rjmcmc_alpha = args_getdouble(optarg);
        opt_rjmcmc_method = 1;
        break;

      case 23:
        if (opt_rjmcmc_method == 0)
          fatal("Cannot use --rjmcmc_mean with --rjmcmc_epsilon");
        opt_rjmcmc_mean = args_getdouble(optarg);
        opt_rjmcmc_method = 1;
        break;

      case 24:
        if (opt_rjmcmc_method == 1)
          fatal("Cannot use --rjmcmc_epsilon with --rjmcmc_alpha/--rjmcmc_mean");
        opt_rjmcmc_epsilon = args_getdouble(optarg);
        opt_rjmcmc_method = 0;
        break;

      case 25:
        opt_cfile = xstrdup(optarg);
        break;

      case 26:
        opt_usedata = 0;
        break;

      case 27:
        if (!strcmp(optarg,"cpu"))
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
        
      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  if (opt_cfile)
    load_cfile();

  /* check for mandatory options */
  if (opt_streefile || opt_streenewick)
    mand_options++;
  if (opt_outfile)
    mand_options++;
  if (opt_mcmcfile)
    mand_options++;

  /* check for number of independent commands selected */
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
//  if (opt_stree || opt_delimit)
//    commands++;
  if (opt_streefile)
    commands++;
  if (opt_cfile)
    commands++;

  if (opt_streefile && opt_streenewick)
    fatal("Cannot use --stree when using a control file (--cfile)");

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  /* if no command specified, turn on --help */
  if (!commands)
  {
    opt_help = 1;
    return;
  }
  /* check for mandatory options */
  if (!opt_version && !opt_help)
    if (mand_options != mandatory_options_count && !opt_cfile)
      fatal("Mandatory options are:\n\n%s", mandatory_options_list);
}

static void dealloc_switches()
{
  if (opt_streefile) free(opt_streefile);
  if (opt_mapfile) free(opt_mapfile);
  if (opt_outfile) free(opt_outfile);
  if (opt_msafile) free(opt_msafile);
  if (opt_mcmcfile) free(opt_mcmcfile);
  if (opt_reorder) free(opt_reorder);
  if (opt_cfile) free(opt_cfile);
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
          "  --help                    display help information\n"
          "  --version                 display version information\n"
          "  --quiet                   only output warnings and fatal errors to stderr\n"
          "  --samples INT             total number of MCMC samples to log (default: 10000)\n"
          "  --samplefreq INT          log every INT sample (default: 10)\n"
          "  --burnin INT              discard first INT MCMC samples (default: 100)\n"
          "  --seed INT                seed for pseudo-random number generator\n"
          "  --stree                   estimate species tree (not implemented)\n"
          "  --delimit                 estimate species delimitation (not implemented)\n"
          "  --cleandata               remove sites containing ambiguous characters\n"
          "  --finetune_reset          reset finetune steps during MCMC\n"
          "  --finetune_params STRING  specify fine-tuning parameters\n"
          "  --tauprior REAL,REAL      specify prior for species tree ages\n"
          "  --thetaprior REAL,REAL    specify prior for population size\n"
          "  --reorder STRING          reorder sequence of species tips\n"
          "\n"
          "Input and output options:\n"
          "  --stree_file FILENAME     species tree file in newick format\n"
          "  --msa_file FILENAME       PHYLIP file containing loci\n"
          "  --map_file FILENAME       file containing mapping of sequences to species\n"
          "  --output_file FILENAME    output file name\n"
          "  --mcmc_file FILENAME      output file containing logged MCMC samples\n"
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

  /* intiialize random number generators */
  //srand48(opt_seed);
  legacy_init();

  printf("opt_stree = %ld\n", opt_stree);
  printf("opt_delimit = %ld\n", opt_delimit);

  if (opt_help)
  {
    cmd_help();
  }
  else if (!opt_stree && !opt_delimit)
  {
    cmd_a00();
  }
  else if (!opt_stree)
  {
    cmd_a10();
  }
  else if (!opt_delimit)
  {
    printf("HERE!!!\n");
    cmd_a01();
  }

  dealloc_switches();
  free(cmdline);
  return (0);
}
