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

static char * progname;
static char progheader[80];
char * cmdline;

/* global error message buffer */
__thread int bpp_errno;
__thread char bpp_errmsg[200] = {0};

/* number of mandatory options for the user to input */
static const char mandatory_options_count = 2;
static const char * mandatory_options_list = " --stree_file --output_file";

/* options */
long opt_help;
long opt_version;
long opt_quiet;
long opt_mcmc_steps;
long opt_mcmc_rate;
long opt_mcmc_burnin;
long opt_seed;
long opt_stree;
long opt_delimit;
long opt_cleandata;
double opt_tau_alpha;
double opt_tau_beta;
double opt_theta_alpha;
double opt_theta_beta;
char * opt_streefile;
char * opt_mapfile;
char * opt_outfile;
char * opt_msafile;
char * opt_mapfile;

static struct option long_options[] =
{
  {"help",        no_argument,       0, 0 },  /*  0 */
  {"version",     no_argument,       0, 0 },  /*  1 */
  {"quiet",       no_argument,       0, 0 },  /*  2 */
  {"stree_file",  required_argument, 0, 0 },  /*  3 */
  {"map_file",    required_argument, 0, 0 },  /*  4 */
  {"output_file", required_argument, 0, 0 },  /*  5 */
  {"msa_file",    required_argument, 0, 0 },  /*  6 */
  {"seed",        required_argument, 0, 0 },  /*  7 */
  {"mcmc_rate",   required_argument, 0, 0 },  /*  8 */
  {"mcmc_burnin", required_argument, 0, 0 },  /*  9 */
  {"mcmc_steps",  required_argument, 0, 0 },  /* 10 */
  {"stree",       no_argument,       0, 0 },  /* 11 */
  {"delimit",     no_argument,       0, 0 },  /* 12 */
  {"tauprior",    required_argument, 0, 0 },  /* 13 */
  {"thetaprior",  required_argument, 0, 0 },  /* 14 */
  {"cleandata",   no_argument,       0, 0 },  /* 15 */
  {"map_file",    required_argument, 0, 0 },  /* 16 */
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

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  int mand_options = 0;

  /* set defaults */

  progname = argv[0];

  opt_help = 0;
  opt_version = 0;
  opt_quiet = 0;
  opt_streefile = NULL;
  opt_outfile = NULL;
  opt_msafile = NULL;
  opt_mapfile = NULL;
  opt_mcmc_steps = 100000;
  opt_mcmc_rate = 10000;
  opt_mcmc_burnin = 10000;
  opt_seed = (long)time(NULL);
  opt_stree = 0;
  opt_delimit = 0;
  opt_tau_alpha = 0;
  opt_tau_beta = 0;
  opt_theta_alpha = 0;
  opt_theta_beta = 0;
  opt_cleandata = 0;

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
        opt_streefile = optarg;
        break;

      case 4:
        opt_mapfile = optarg;
        break;

      case 5:
        opt_outfile = optarg;
        break;

      case 6:
        opt_msafile = optarg;
        break;

      case 7:
        opt_seed = atol(optarg);
        break;

      case 8:
        opt_mcmc_rate = atol(optarg);
        break;

      case 9:
        opt_mcmc_burnin = atol(optarg);
        break;

      case 10:
        opt_mcmc_steps = atol(optarg);
        break;

      case 11:
        opt_stree = 1;
        break;

      case 12:
        opt_delimit = 1;
        break;

      case 13:
        if (!args_getgammaprior(optarg, &opt_tau_alpha, &opt_tau_beta))
          fatal("Illegal format for --tauprior");
        break;

      case 14:
        if (!args_getgammaprior(optarg, &opt_theta_alpha, &opt_theta_beta))
          fatal("Illegal format for --thetaprior");
        break;

      case 15:
        opt_cleandata = 1;
        break;

      case 16:
        opt_mapfile = optarg;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for mandatory options */
  if (opt_streefile)
    mand_options++;
  if (opt_outfile)
    mand_options++;

  /* check for number of independent commands selected */
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
  if (opt_stree || opt_delimit)
    commands++;
  if (opt_streefile)
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
  /* check for mandatory options */
  if (!opt_version && !opt_help)
    if (mand_options != mandatory_options_count)
      fatal("Mandatory options are:\n\n%s", mandatory_options_list);

}

void cmd_help()
{
  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help                  display help information.\n"
          "  --version               display version information.\n"
          "  --quiet                 only output warnings and fatal errors to stderr.\n"
          "  --mcmc_steps INT        Perform INT MCMC steps (default: 100000).\n"
          "  --mcmc_rate INT         Sample every INT step (default: 1000).\n"
          "  --mcmc_burnin INT       discard the first INT MCMC steps (default: 10000).\n"
          "  --seed INT              Seed for pseudo-random number generator.\n"
          "  --stree                 estimate species tree.\n"
          "  --delimit               estimates species delimitation.\n"
          "\n"
          "Input and output options:\n"
          "  --stree_file FILENAME   species tree file in newick format.\n"
          "  --msa_file FILENAMES    comma separated list of loci MSAs.\n"
          "  --map_file FILENAME     map file connecting  "
          "  --output_file FILENAME  output file name.\n"
         );
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
           sysconf(_SC_NPROCESSORS_ONLN));
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

  srand48(opt_seed);

  if (opt_help)
  {
    cmd_help();
  }
  else if (!opt_stree && !opt_delimit)
  {
    cmd_a00();
  }

  free(cmdline);
  return (0);
}
