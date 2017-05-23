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

#define BPP_FAILURE  0
#define BPP_SUCCESS  1

#define LINEALLOC 2048
#define ASCII_SIZE 256

#define RTREE_SHOW_LABEL                1
#define RTREE_SHOW_BRANCH_LENGTH        2

#define TREE_TRAVERSE_POSTORDER         1
#define TREE_TRAVERSE_PREORDER          2

/* error codes */

#define ERROR_PHYLIP_SYNTAX            106
#define ERROR_PHYLIP_LONGSEQ           107
#define ERROR_PHYLIP_NONALIGNED        108
#define ERROR_PHYLIP_ILLEGALCHAR       109
#define ERROR_PHYLIP_UNPRINTABLECHAR   110

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

typedef struct rnode_s
{
  char * label;
  double length;
  struct rnode_s * left;
  struct rnode_s * right;
  struct rnode_s * parent;
  unsigned int leaves;

  species_t * species_data;
  gtree_data_t * gtree_data;

  unsigned int index;
} rnode_t;

typedef struct rtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  rnode_t ** nodes;

  rnode_t * root;
} rtree_t;

typedef struct msa_s
{
  int count;
  int length;

  char ** sequence;
  char ** label;

} msa_t;

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

/* Simple structure for handling PHYLIP parsing */

typedef struct fasta
{
  FILE * fp;
  char line[LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} fasta_t;


/* Simple structure for handling PHYLIP parsing */

typedef struct phylip_s
{
  FILE * fp;
  char * line;
  size_t line_size;
  size_t line_maxsize;
  char buffer[LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} phylip_t;

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
#define SWAP(x,y) do { __typeof__ (x) _t = x; x = y; y = _t; } while(0)

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

extern __thread int bpp_errno;
extern __thread char bpp_errmsg[200];

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
void * xcalloc(size_t nmemb, size_t size);
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

/* functions in phylip.c */

phylip_t * phylip_open(const char * filename,
                       const unsigned int * map);

int phylip_rewind(phylip_t * fd);

void phylip_close(phylip_t * fd);

msa_t * phylip_parse_interleaved(phylip_t * fd);

msa_t * phylip_parse_sequential(phylip_t * fd);

msa_t ** phylip_parse_multisequential(phylip_t * fd, int * count);

/* functions in rtree.c */

void rtree_show_ascii(const rnode_t * root, int options);
char * rtree_export_newick(const rnode_t * root,
                           char * (*cb_serialize)(const rnode_t *));
int rtree_traverse(rnode_t * root,
                   int traversal,
                   int (*cbtrav)(rnode_t *),
                   rnode_t ** outbuffer,
                   unsigned int * trav_size);
rtree_t ** rtree_tipstring_nodes(rtree_t * root,
                                 char * tipstring,
                                 unsigned int * tiplist_count);

/* functions in arch.c */

unsigned long arch_get_memused(void);
unsigned long arch_get_memtotal(void);

/* functions in msa.c */

void msa_print(msa_t * msa);

void msa_destroy(msa_t * msa);

/* functions in parse_map.y */

list_t * yy_parse_map(const char * filename);

/* functions in mapping.c */

void maplist_print(list_t * map_list);
void maplist_destroy(list_t * map_list);
