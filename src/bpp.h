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

#include <assert.h>
#include <search.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#if (!defined(WINNT) && !defined(WIN32) && !defined(WIN64))
#include <sys/resource.h>
#endif

/* constants */

#define PROG_NAME "bpp"
#define PROG_VERSION "v0.0.0"

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

#define LINEALLOC 2048

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

typedef struct species_s
{
  double theta;
  double tau;
  char * label;
  int * perloci_seqcount;
} species_t;

typedef struct gtree_data_s
{
  species_t * species;
  double time;
} gtree_data_t;

typedef struct rtree_s
{
  char * label;
  double length;
  struct rtree_s * left;
  struct rtree_s * right;
  struct rtree_s * parent;
  int leaves;

  species_t * species_data;
  gtree_data_t * gtree_data;

  int node_index;

} rtree_t;

typedef struct alignment_s
{
  int seq_count;
  int seq_len;

  char ** seq;
  char ** label;
  unsigned int * weight;

} alignment_t;

typedef struct locus_s
{
  unsigned int tips_count;
  unsigned int clv_count;
  unsigned int pmatrix_count;
  unsigned int states;
  unsigned int sites;

  double ** clv;
  double ** pmatrix;
  unsigned int * pattern_weights;
  unsigned char ** tipstates;
} locus_t;

typedef struct pll_fasta
{
  FILE * fp;
  char line[LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} pll_fasta_t;

typedef struct map_s
{
  char * individual;
  char * species;
  int lineno;
} map_t;

typedef struct list_s
{
  void * data;
  struct list_s * next;
} list_t;

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* options */

extern long opt_help;
extern long opt_version;
extern long opt_quiet;
extern long opt_mcmc_steps;
extern long opt_mcmc_rate;
extern long opt_mcmc_burnin;
extern long opt_seed;
extern long opt_stree;
extern long opt_delimit;
extern long opt_cleandata;
extern double opt_tau_alpha;
extern double opt_tau_beta;
extern double opt_theta_alpha;
extern double opt_theta_beta;
extern char * opt_streefile;
extern char * opt_mapfile;
extern char * opt_outfile;
extern char * opt_msafile;
extern char * opt_mapfile;
extern char * cmdline;

/* common data */

extern char errmsg[200];

extern int pll_errno;
extern const unsigned int pll_map_nt[256];
extern const unsigned int pll_map_fasta[256];

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

/* functions in util.c */

void fatal(const char * format, ...) __attribute__ ((noreturn));
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned int progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
void show_rusage(void);
FILE * xopen(const char * filename, const char * mode);

/* functions in bpp.c */

void args_init(int argc, char ** argv);
void cmd_help(void);
void getentirecommandline(int argc, char * argv[]);
void fillheader(void);
void show_header(void);
void cmd_ml(void);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);
void rtree_destroy(rtree_t * root);

/* functions in parse_phylip.y */

alignment_t ** yy_parse_phylip(const char * filename, int * loci_count);

/* functions in rtree.c */

void rtree_show_ascii(rtree_t * tree);
char * rtree_export_newick(rtree_t * root);
int rtree_query_tipnodes(rtree_t * root, rtree_t ** node_list);
int rtree_query_innernodes(rtree_t * root, rtree_t ** node_list);
void rtree_reset_info(rtree_t * root);
void rtree_print_tips(rtree_t * node, FILE * out);
int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   struct drand48_data * rstate,
                   rtree_t ** outbuffer);
rtree_t * rtree_clone(rtree_t * node, rtree_t * parent);
int rtree_traverse_postorder(rtree_t * root,
                             int (*cbtrav)(rtree_t *),
                             rtree_t ** outbuffer);
rtree_t ** rtree_tipstring_nodes(rtree_t * root,
                                 char * tipstring,
                                 unsigned int * tiplist_count);

/* functions in arch.c */

unsigned long arch_get_memused(void);
unsigned long arch_get_memtotal(void);

/* functions in alignment.c */

void alignment_print(alignment_t * alignment);

void alignment_destroy(alignment_t * alignment);

/* functions in parse_map.y */

list_t * yy_parse_map(const char * filename);

/* functions in mapping.c */

void maplist_print(list_t * map_list);
void maplist_destroy(list_t * map_list);
