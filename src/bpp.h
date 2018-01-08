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

#include <assert.h>
#include <fcntl.h>
#include <search.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/stat.h>
#include <stdint.h>

#ifdef _MSC_VER
#include <pmmintrin.h>
#include <immintrin.h>
#else
#include <x86intrin.h>
#endif

#ifndef _MSC_VER
#include <getopt.h>
#endif

#ifndef _MSC_VER
#include <sys/time.h>
#include <unistd.h>
#endif

/* platform specific */

#if (defined(__BORLANDC__) || defined(_MSC_VER))
#define __THREAD __declspec(thread)
#else
#define __THREAD __thread
#endif

#ifdef _MSC_VER
#define PLL_ALIGN_HEADER(X) __declspec(align(X))
#define PLL_ALIGN_FOOTER(X)
#else
#define PLL_ALIGN_HEADER(X)
#define PLL_ALIGN_FOOTER(X) __attribute__((aligned(X)))
#endif

#ifndef _MSC_VER
#define xasprintf asprintf
#endif

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

/* constants */

#define PROG_NAME "bpp"

#define PLL_STRING(x) #x
#define PLL_C2S(x) PLL_STRING(x)


#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_PATCH 0

/* checkpoint version */
#define VERSION_CHKP 1

#define PROG_VERSION "v" PLL_C2S(VERSION_MAJOR) "." PLL_C2S(VERSION_MINOR) "." \
        PLL_C2S(VERSION_PATCH)

#define BPP_MAGIC_BYTES 4
#define BPP_MAGIC "BPPX"

#ifdef __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#else
#error "Big endian ppc64 CPUs not supported"
#endif

#else

#define PROG_CPU "x86_64"

#endif

#ifdef __APPLE__
#define PROG_OS "osx"
#include <sys/resource.h>
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#define PROG_OS "linux"
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef _WIN32
#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>
#endif

#define PROG_ARCH PROG_OS "_" PROG_CPU

#define BPP_FAILURE  0
#define BPP_SUCCESS  1

#define LINEALLOC 2048
#define ASCII_SIZE 256

#define RTREE_SHOW_LABEL                1
#define RTREE_SHOW_BRANCH_LENGTH        2

#define TREE_TRAVERSE_POSTORDER         1
#define TREE_TRAVERSE_PREORDER          2

#define COMPRESS_GENERAL                1
#define COMPRESS_JC69                   2

#define BPP_DELIMIT_PRIOR_DIRICHLET     0
#define BPP_DELIMIT_PRIOR_UNIFORM       1

/* libpll related definitions */

#define PLL_ALIGNMENT_CPU   8
#define PLL_ALIGNMENT_SSE  16
#define PLL_ALIGNMENT_AVX  32

#define PLL_ATTRIB_ARCH_CPU            0
#define PLL_ATTRIB_ARCH_SSE       (1 << 0)
#define PLL_ATTRIB_PATTERN_TIP    (1 << 4)

#define PLL_ATTRIB_ARCH_AVX       (1 << 1)
#define PLL_ATTRIB_ARCH_AVX2      (1 << 2)
#define PLL_ATTRIB_ARCH_AVX512    (1 << 3)
#define PLL_ATTRIB_ARCH_MASK         0xF

#define PLL_ATTRIB_PATTERN_TIP    (1 << 4)

#define PLL_ATTRIB_RATE_SCALERS   (1 << 9)

#define PLL_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define PLL_SCALE_THRESHOLD (1.0/PLL_SCALE_FACTOR)
#define PLL_SCALE_FACTOR_SQRT 340282366920938463463374607431768211456.0 /* 2**128 */
#define PLL_SCALE_THRESHOLD_SQRT (1.0/PLL_SCALE_FACTOR_SQRT)
#define PLL_SCALE_BUFFER_NONE -1

#define PLL_MISC_EPSILON 1e-8

/* error codes */

#define ERROR_PHYLIP_SYNTAX            106
#define ERROR_PHYLIP_LONGSEQ           107
#define ERROR_PHYLIP_NONALIGNED        108
#define ERROR_PHYLIP_ILLEGALCHAR       109
#define ERROR_PHYLIP_UNPRINTABLECHAR   110

/* available methods */

#define METHOD_00       0
#define METHOD_01       1
#define METHOD_10       2
#define METHOD_11       3

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

typedef struct dlist_item_s
{
  void * data;
  struct dlist_item_s * prev;
  struct dlist_item_s * next;
} dlist_item_t;

typedef struct dlist_s
{
  dlist_item_t * head;
  dlist_item_t * tail;
} dlist_t;

typedef struct snode_s
{
  char * label;
  double length;
  double theta;
  double tau;
  double old_tau;
  double old_theta;
  struct snode_s * left;
  struct snode_s * right;
  struct snode_s * parent;
  unsigned int leaves;
  unsigned int * gene_leaves;
  int mark;

  void * data;

  /* clade support when delimiting species */
  double support;

  /* branch weight when inferring species tree */
  double weight;

  /* list of per-locus coalescent events */
  dlist_t ** event;

  int * event_count;

  /* number of lineages coming in the population */
  int * seqin_count;
  double * logpr_contrib;
  double * old_logpr_contrib;
  //unsigned int * seqin_count;
  //unsigned int * seqout_count;

  //unsigned int ** seqin_indices;
  unsigned long * bitmask;

  unsigned int node_index;
  unsigned int diploid;
} snode_t;

typedef struct stree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  unsigned int locus_count;

  snode_t ** nodes;

  int ** pptable;

  snode_t * root;

  double root_age;

//  int mark;
} stree_t;

typedef struct gnode_s
{
  char * label;
  double length;
  double time;
  double old_time;
  struct gnode_s * left;
  struct gnode_s * right;
  struct gnode_s * parent;
  unsigned int leaves;

  void * data;

  snode_t * pop;
  snode_t * old_pop;

  /* pointer to the dlist item this node is wrapped into */
  dlist_item_t * event;

  unsigned int node_index;
  unsigned int clv_valid;

  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;

  int mark;

} gnode_t;

typedef struct gtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  gnode_t ** nodes;
  gnode_t * root;

  /* auxiliary space for traversals */
  double logl;
  double logpr;
  double old_logpr;
  double old_logl;

} gtree_t;

typedef struct msa_s
{
  int count;
  int length;

  char ** sequence;
  char ** label;

  int amb_sites_count;
  int original_length;

} msa_t;

typedef struct locus_s
{
  unsigned int tips;
  unsigned int clv_buffers;
  unsigned int states;
  unsigned int sites;
  unsigned int rate_matrices;
  unsigned int prob_matrices;
  unsigned int rate_cats;
  unsigned int scale_buffers;
  unsigned int attributes;

  /* vectorization parameters */
  size_t alignment;
  unsigned int states_padded;

  double ** clv;
  double ** pmatrix;
  double * rates;
  double * rate_weights;
  double ** subst_params;
  unsigned int ** scale_buffer;
  double ** frequencies;
  unsigned int * pattern_weights;
  unsigned int pattern_weights_sum;

  int * eigen_decomp_valid;
  double ** eigenvecs;
  double ** inv_eigenvecs;
  double ** eigenvals;

  /* tip-tip precomputation data */
  unsigned int maxstates;
  unsigned char ** tipchars;
  unsigned char * charmap;
  double * ttlookup;
  unsigned int * tipmap;

  /* diploid related */
  int diploid;
  unsigned long * diploid_mapping;
  unsigned long * diploid_resolution_count;
  double * likelihood_vector;
  int unphased_length;


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

typedef struct mapping_s
{
  char * individual;
  char * species;
  int lineno;
} mapping_t;

typedef struct list_item_s
{
  void * data;
  struct list_item_s * next;
} list_item_t;

typedef struct list_s
{
  list_item_t * head;
  list_item_t * tail;
  long count;
} list_t;

typedef struct ht_item_s
{
  unsigned long key;
  void * value;
} ht_item_t;

typedef struct hashtable_s
{
  unsigned long table_size;
  unsigned long entries_count;
  list_t ** entries;
} hashtable_t;

typedef struct pair_s
{
  char * label;
  void * data;
} pair_t;

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef _MSC_VER
#define SWAP(x,y) do                                                  \
  {                                                                   \
    size_t s = MAX(sizeof(x),sizeof(y));                              \
    unsigned char * temp = (unsigned char *)malloc(s*sizeof(char));   \
    memcpy(temp,&y,s);                                                \
    memcpy(&y,&x,s);                                                  \
    memcpy(&x,temp,s);                                                \
    free(temp);                                                       \
  } while(0)
#else
#define SWAP(x,y) do { __typeof__ (x) _t = x; x = y; y = _t; } while(0)
#endif

#ifdef _MSC_VER
#define PLL_POPCOUNT pll_popcount
#define PLL_POPCOUNTL pll_popcount64
#define PLL_CTZ pll_ctz
#else
#define PLL_POPCOUNT __builtin_popcount
#define PLL_POPCOUNTL __builtin_popcountl
#define PLL_CTZ __builtin_ctz
#endif

#define legacy_rndexp(mean) (-(mean)*log(legacy_rndu()))

#define FLAG_AGE_UPDATE                 1 
#define FLAG_POP_UPDATE                 2
#define FLAG_BRANCH_UPDATE              4
#define FLAG_MISC                       8
#define FLAG_PARTIAL_UPDATE           128


/* options */

extern long opt_help;
extern long opt_version;
extern long opt_quiet;
extern long opt_seed;
extern long opt_stree;
extern long opt_arch;
extern long opt_delimit;
extern long opt_delimit_prior;
extern long opt_cleandata;
extern long opt_debug;
extern long opt_est_theta;
extern long opt_samples;
extern long opt_samplefreq;
extern long opt_burnin;
extern long opt_finetune_reset;
extern long opt_rjmcmc_method;
extern long opt_usedata;
extern long opt_nloci;
extern long opt_experimental_method;
extern long opt_experimental_debug;
extern long opt_diploid_size;
extern long opt_checkpoint_initial;
extern long opt_checkpoint_step;
extern long opt_checkpoint_current;
extern long opt_method;
extern double opt_bfbeta;
extern double opt_tau_alpha;
extern double opt_tau_beta;
extern double opt_theta_alpha;
extern double opt_theta_beta;
extern double opt_finetune_gtage;
extern double opt_finetune_gtspr;
extern double opt_finetune_theta;
extern double opt_finetune_tau;
extern double opt_finetune_mix;
extern double opt_rjmcmc_alpha;
extern double opt_rjmcmc_mean;
extern double opt_rjmcmc_epsilon;
extern long * opt_diploid;
extern char * opt_cfile;
extern char * opt_mapfile;
extern char * opt_msafile;
extern char * opt_mapfile;
extern char * opt_mcmcfile;
extern char * opt_reorder;
extern char * opt_outfile;
extern char * opt_streefile;
extern char * opt_streenewick;
extern char * opt_resume;
extern char * cmdline;

/* common data */

extern __THREAD int bpp_errno;
extern __THREAD char bpp_errmsg[200];

extern const unsigned int pll_map_nt[256];
extern const unsigned int pll_map_fasta[256];
extern const unsigned int pll_map_amb[256];
extern const unsigned int pll_map_validjc69[16];

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
extern long altivec_present;

/* functions in util.c */

#ifdef _MSC_VER
__declspec(noreturn) void fatal(const char * format, ...);
int xasprintf(char ** strp, const char * fmt, ...);
#else
void fatal(const char * format, ...) __attribute__ ((noreturn));
#endif
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xcalloc(size_t nmemb, size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
FILE * xopen(const char * filename, const char * mode);
void * pll_aligned_alloc(size_t size, size_t alignment);
void pll_aligned_free(void * ptr);

/* functions in bpp.c */

void args_init(int argc, char ** argv);

void cmd_help(void);

void getentirecommandline(int argc, char * argv[]);

void fillheader(void);

void show_header(void);

void cmd_ml(void);

/* functions in parse_stree.y */

stree_t * stree_parse_newick(const char * filename);

stree_t * stree_parse_newick_string(const char * s);

void stree_destroy(stree_t * tree,
                   void (*cb_destroy)(void *));

/* functions in phylip.c */

phylip_t * phylip_open(const char * filename,
                       const unsigned int * map);

int phylip_rewind(phylip_t * fd);

void phylip_close(phylip_t * fd);

msa_t * phylip_parse_interleaved(phylip_t * fd);

msa_t * phylip_parse_sequential(phylip_t * fd);

msa_t ** phylip_parse_multisequential(phylip_t * fd, int * count);

/* functions in stree.c */

void stree_show_ascii(const snode_t * root, int options);

char * stree_export_newick(const snode_t * root,
                           char * (*cb_serialize)(const snode_t *));

int stree_traverse(snode_t * root,
                   int traversal,
                   int (*cbtrav)(snode_t *),
                   snode_t ** outbuffer,
                   unsigned int * trav_size);

stree_t ** stree_tipstring_nodes(stree_t * root,
                                 char * tipstring,
                                 unsigned int * tiplist_count);

hashtable_t * species_hash(stree_t * tree);

hashtable_t * maplist_hash(list_t * maplist, hashtable_t * sht);

double stree_propose_theta(gtree_t ** gtree, stree_t * stree);

double stree_propose_tau(gtree_t ** gtree, stree_t * stree, locus_t ** loci);

void stree_fini(void);

void stree_rootdist(stree_t * stree,
                    list_t * maplist,
                    msa_t ** msalist,
                    unsigned int ** weights);

void fill_seqin_counts(stree_t * stree, int msa_index);

void reset_gene_leaves_count(stree_t * stree);

void stree_reset_pptable(stree_t * stree);

long stree_propose_spr(stree_t ** streeptr,
                       gtree_t *** gtree_list_ptr,
                       stree_t ** scloneptr,
                       gtree_t *** gclonesptr,
                       locus_t ** loci);

gtree_t * gtree_clone_init(gtree_t * gtree, stree_t * stree);

stree_t * stree_clone_init(stree_t * stree);

void stree_label(stree_t * stree);

void stree_show_pptable(stree_t * stree);

void stree_init(stree_t * stree, msa_t ** msa, list_t * maplist, int msa_count);

void stree_alloc_internals(stree_t * stree,
                           unsigned int gtree_inner_sum,
                           long msa_count);


/* functions in arch.c */

unsigned long arch_get_memused(void);

unsigned long arch_get_memtotal(void);

long arch_get_cores(void);

/* functions in msa.c */

void msa_print(msa_t * msa);

void msa_destroy(msa_t * msa);

int msa_remove_ambiguous(msa_t * msa);

void msa_summary(msa_t ** msa_list, int msa_count);

void msa_count_ambiguous_sites(msa_t * msa, const unsigned int * map);

/* functions in parse_map.y */

list_t * yy_parse_map(const char * filename);

/* functions in mapping.c */

void maplist_print(list_t * map_list);

void map_dealloc(void * data);

hashtable_t * maplist_hash(list_t * maplist, hashtable_t * sht);

/* functions in list.c */

void list_append(list_t * list, void * data);

void list_prepend(list_t * list, void * data);

void list_clear(list_t * list, void (*cb_dealloc)(void *));

/* functions in dlist.c */

dlist_item_t * dlist_append(dlist_t * dlist, void * data);
void dlist_clear(dlist_t * dlist, void (*cb_dealloc)(void *));
void dlist_item_remove(dlist_item_t * item);
void dlist_item_append(dlist_t * dlist, dlist_item_t * item);
void dlist_item_prepend(dlist_t * dlist, dlist_item_t * item);
dlist_t * dlist_create();
void dlist_destroy(dlist_t *);

/* functions in hash.c */

void * hashtable_find(hashtable_t * ht,
                      void * x,
                      unsigned long hash,
                      int (*cb_cmp)(void *, void *));

hashtable_t * hashtable_create(unsigned long items_count);

int hashtable_strcmp(void * x, void * y);

int hashtable_ptrcmp(void * x, void * y);

unsigned long hash_djb2a(char * s);

unsigned long hash_fnv(char * s);

int hashtable_insert(hashtable_t * ht,
                     void * x,
                     unsigned long hash,
                     int (*cb_cmp)(void *, void *));

void hashtable_insert_force(hashtable_t * ht,
                            void * x,
                            unsigned long hash);

void hashtable_destroy(hashtable_t * ht, void (*cb_dealloc)(void *));

int cb_cmp_pairlabel(void * a, void * b);

/* functions in random.c */

double legacy_rndu(void);
double legacy_rnd_symmetrical(void);
void legacy_init(void);
double legacy_rndbeta (double p, double q);
double legacy_rndgamma (double a);
unsigned int get_legacy_rndu_status(void);
void set_legacy_rndu_status(unsigned int x);

/* functions in gtree.c */

void gtree_alloc_internals(gtree_t ** gtree, long msa_count);

gtree_t ** gtree_init(stree_t * stree,
                      msa_t ** msalist,
                      list_t * maplist,
                      int msa_count);

char * gtree_export_newick(const gnode_t * root,
                           char * (*cb_serialize)(const gnode_t *));

void gtree_destroy(gtree_t * tree, void (*cb_destroy)(void *));

int gtree_traverse(gnode_t * root,
                   int traversal,
                   int (*cbtrav)(gnode_t *),
                   gnode_t ** outbuffer,
                   unsigned int * trav_size);

void gtree_update_branch_lengths(gtree_t ** gtree_list, int msa_count);

double gtree_propose_ages(locus_t ** locus, gtree_t ** gtree, stree_t * stree);

void gtree_reset_leaves(gnode_t * node);

void gtree_fini(int msa_count);

double gtree_logprob(stree_t * stree, long msa_index);
double gtree_update_logprob_contrib(snode_t * snode, long msa_index);
double gtree_propose_spr(locus_t ** locus, gtree_t ** gtree, stree_t * stree);
double reflect(double t, double minage, double maxage);
gnode_t ** gtree_return_partials(gnode_t * root,
                                 unsigned int msa_index,
                                 unsigned int * trav_size);
void unlink_event(gnode_t * node, int msa_index);

/* functions in prop_mixing.c */

long proposal_mixing(gtree_t ** gtree, stree_t * stree, locus_t ** locus);

/* functions in prop_rj.c */

long prop_split(gtree_t ** gtree,
                stree_t * stree,
                locus_t ** locus,
                double pr_split,
                long * param_count);

long prop_join(gtree_t ** gtree,
               stree_t * stree,
               locus_t ** locus,
               double pr_split,
               long * param_count);

void rj_init(gtree_t ** gtreelist, stree_t * stree, unsigned int count);

void rj_fini();

double lnprior_species_model(stree_t * stree); /* TODO: Move function */

/* functions in locus.c */

locus_t * locus_create(unsigned int tips,
                       unsigned int clv_buffers,
                       unsigned int states,
                       unsigned int sites,
                       unsigned int rate_matrices,
                       unsigned int prob_matrices,
                       unsigned int rate_cats,
                       unsigned int scale_buffers,
                       unsigned int attributes);

void locus_destroy(locus_t * locus);

int pll_set_tip_states(locus_t * locus,
                       unsigned int tip_index,
                       const unsigned int * map,
                       const char * sequence);

int pll_set_tip_clv(locus_t * locus,
                    unsigned int tip_index,
                    const double * clv,
                    int padding);

void pll_set_frequencies(locus_t * locus,
                         unsigned int freqs_index,
                         const double * frequencies);

void locus_update_partials(locus_t * locus, gnode_t ** traversal, unsigned int count);

void locus_update_all_partials(locus_t * locus, gtree_t * gtree);

void pll_set_pattern_weights(locus_t * locus,
                             const unsigned int * pattern_weights);

void locus_update_matrices_jc69(locus_t * locus,
                                gnode_t ** traversal,
                                unsigned int count);

double locus_root_loglikelihood(locus_t * locus,
                                gnode_t * root,
                                const unsigned int * freqs_indices,
                                double * persite_lnl);

/* functions in compress.c */

unsigned int * compress_site_patterns(char ** sequence,
                                      const unsigned int * map,
                                      int count,
                                      int * length,
                                      int attrib);

unsigned long * compress_site_patterns_diploid(char ** sequence,
                                               const unsigned int * map,
                                               int count,
                                               int * length,
                                               int attrib);
/* functions in allfixed.c */

void allfixed_summary(stree_t * stree);

/* functions in summary.c */

void bipartitions_init(char ** species, long species_count);

void bipartitions_update(stree_t * stree);

void bipartitions_finalize(size_t trees_count, char ** species);

void print_stree_with_support(const char * treestr, size_t freq, size_t trees_count);

void summary_dealloc_hashtables(void);

void stree_summary(char ** species_names, long species_count);

/* functions in hardware.c */

void cpu_features_show(void);

void cpu_features_detect(void);

void cpu_setarch(void);

#ifdef _MSC_VER
int pll_ctz(unsigned int x);
unsigned int pll_popcount(unsigned int x);
unsigned int pll_popcount64(unsigned long x);
#endif

/* functions in diploid.c */

unsigned long ** diploid_resolve(stree_t * stree,
                                 msa_t ** msa_list,
                                 list_t * maplist,
                                 unsigned int ** weights,
                                 int msa_count,
                                 const unsigned int * map);

/* functions in dump.c */

int checkpoint_dump(stree_t * stree,
                    gtree_t ** gtree_list,
                    locus_t ** locus_list,
                    double * pjump,
                    long curstep,
                    long ft_round,
                    long mcmc_offset);

/* functions in load.c */

int checkpoint_load(gtree_t *** gtreep, locus_t *** locusp, stree_t ** streep, double ** pjump, long * curstep, long * ft_round, long * mcmc_offset);
void checkpoint_truncate(long mcmc_offset);

/* functions in core_partials.c */

void pll_core_update_partial_tt_4x4(unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const unsigned char * right_tipchars,
                                    const double * lookup,
                                    unsigned int attrib);

void pll_core_update_partial_tt(unsigned int states,
                                unsigned int sites,
                                unsigned int rate_cats,
                                double * parent_clv,
                                unsigned int * parent_scaler,
                                const unsigned char * left_tipchars,
                                const unsigned char * right_tipchars,
                                const unsigned int * tipmap,
                                unsigned int tipmap_size,
                                const double * lookup,
                                unsigned int attrib);

void pll_core_update_partial_ti_4x4(unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * right_scaler,
                                    unsigned int attrib);

void pll_core_update_partial_ti(unsigned int states,
                                unsigned int sites,
                                unsigned int rate_cats,
                                double * parent_clv,
                                unsigned int * parent_scaler,
                                const unsigned char * left_tipchars,
                                const double * right_clv,
                                const double * left_matrix,
                                const double * right_matrix,
                                const unsigned int * right_scaler,
                                const unsigned int * tipmap,
                                unsigned int tipmap_size,
                                unsigned int attrib);

void pll_core_update_partial_ii(unsigned int states,
                                unsigned int sites,
                                unsigned int rate_cats,
                                double * parent_clv,
                                unsigned int * parent_scaler,
                                const double * left_clv,
                                const double * right_clv,
                                const double * left_matrix,
                                const double * right_matrix,
                                const unsigned int * left_scaler,
                                const unsigned int * right_scaler,
                                unsigned int attrib);

void pll_core_create_lookup_4x4(unsigned int rate_cats,
                                double * lookup,
                                const double * left_matrix,
                                const double * right_matrix);

void pll_core_create_lookup(unsigned int states,
                            unsigned int rate_cats,
                            double * lookup,
                            const double * left_matrix,
                            const double * right_matrix,
                            const unsigned int * tipmap,
                            unsigned int tipmap_size,
                            unsigned int attrib);

/* functions in core_pmatrix.c */

int pll_core_update_pmatrix(double ** pmatrix,
                            unsigned int states,
                            unsigned int rate_cats,
                            const double * rates,
                            const double * branch_lengths,
                            const unsigned int * matrix_indices,
                            const unsigned int * params_indices,
                            double * const * eigenvals,
                            double * const * eigenvecs,
                            double * const * inv_eigenvecs,
                            unsigned int count,
                            unsigned int attrib);

int pll_core_update_pmatrix_4x4_jc69(double ** pmatrix,
                                     unsigned int states,
                                     unsigned int rate_cats,
                                     const double * rates,
                                     const double * branch_lengths,
                                     const unsigned int * matrix_indices,
                                     const unsigned int * params_indices,
                                     unsigned int count,
                                     unsigned int attrib);

/* functions in core_likelihood.c */

double pll_core_root_loglikelihood(unsigned int states,
                                   unsigned int sites,
                                   unsigned int rate_cats,
                                   const double * clv,
                                   const unsigned int * scaler,
                                   double * const * frequencies,
                                   const double * rate_weights,
                                   const unsigned int * pattern_weights,
                                   const unsigned int * freqs_indices,
                                   double * persite_lnl,
                                   unsigned int attrib);

void pll_core_root_likelihood_vector(unsigned int states,
                                     unsigned int sites,
                                     unsigned int rate_cats,
                                     const double * clv,
                                     const unsigned int * scaler,
                                     double * const * frequencies,
                                     const double * rate_weights,
                                     const unsigned int * pattern_weights,
                                     const unsigned int * freqs_indices,
                                     double * persite_lnl,
                                     unsigned int attrib);
/* functions in output.c */

void pll_show_pmatrix(const locus_t * locus,
                                 unsigned int index,
                                 unsigned int float_precision);

void pll_show_clv(const locus_t * locus,
                             unsigned int clv_index,
                             int scaler_index,
                             unsigned int float_precision);

/* functions in experimental.c */

void experimental_tselect_logl(gnode_t * mnode,
                               gnode_t ** target_list,
                               long target_count,
                               locus_t * locus,
                               double * weights);


/* functions in method.c */

void cmd_run(void);

///* functions in method_00.c */
//
//void cmd_a00(void);
//void cmd_a00_chk(void);
//
///* functions in method_10.c */
//
//void cmd_a10(void);
//
///* functions in method_01.c */
//
//void cmd_a01(void);

/* functions in delimit.c */

long delimitations_count(stree_t * stree);

long delimitations_init(stree_t * stree);

long histories(stree_t * stree);

void delimitations_fini(void);

void delimitation_set(stree_t * stree, long index);

long delimitation_getparam_count(void);

char * delimitation_getparam_string();

long delimit_getindex(stree_t * stree);
long delimit_getindexfromstring(char * model);

void delimit_setindex(long index);

void delimit_resetpriors(void);

void delimit_summary(stree_t * stree);

/* functions in cfile.c */

void load_cfile(void);

/* functions in core_partials_sse.c */

#ifdef HAVE_SSE3
void pll_core_create_lookup_sse(unsigned int states,
                                unsigned int rate_cats,
                                double * ttlookup,
                                const double * left_matrix,
                                const double * right_matrix,
                                const unsigned int * tipmap,
                                unsigned int tipmap_size);

void pll_core_create_lookup_4x4_sse(unsigned int rate_cats,
                                    double * lookup,
                                    const double * left_matrix,
                                    const double * right_matrix);

void pll_core_update_partial_tt_sse(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const unsigned char * right_tipchars,
                                    const double * lookup,
                                    unsigned int tipstates_count,
                                    unsigned int attrib);

void pll_core_update_partial_tt_4x4_sse(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const unsigned char * left_tipchars,
                                        const unsigned char * right_tipchars,
                                        const double * lookup,
                                        unsigned int attrib);

void pll_core_update_partial_ti_sse(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * right_scaler,
                                    const unsigned int * tipmap,
                                    unsigned int tipmap_size,
                                    unsigned int attrib);


void pll_core_update_partial_ti_4x4_sse(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const unsigned char * left_tipchar,
                                        const double * right_clv,
                                        const double * left_matrix,
                                        const double * right_matrix,
                                        const unsigned int * right_scaler,
                                        unsigned int attrib);

void pll_core_update_partial_ii_sse(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const double * left_clv,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * left_scaler,
                                    const unsigned int * right_scaler,
                                    unsigned int attrib);

void pll_core_update_partial_ii_4x4_sse(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const double * left_clv,
                                        const double * right_clv,
                                        const double * left_matrix,
                                        const double * right_matrix,
                                        const unsigned int * left_scaler,
                                        const unsigned int * right_scaler,
                                        unsigned int attrib);

/* functions in core_likelihood_sse.c */


double pll_core_root_loglikelihood_sse(unsigned int states,
                                       unsigned int sites,
                                       unsigned int rate_cats,
                                       const double * clv,
                                       const unsigned int * scaler,
                                       double * const * frequencies,
                                       const double * rate_weights,
                                       const unsigned int * pattern_weights,
                                       const unsigned int * freqs_indices,
                                       double * persite_lnl);

double pll_core_root_loglikelihood_4x4_sse(unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * clv,
                                           const unsigned int * scaler,
                                           double * const * frequencies,
                                           const double * rate_weights,
                                           const unsigned int * pattern_weights,
                                           const unsigned int * freqs_indices,
                                           double * persite_lnl);

void pll_core_root_likelihood_vec_sse(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      const double * clv,
                                      const unsigned int * scaler,
                                      double * const * frequencies,
                                      const double * rate_weights,
                                      const unsigned int * pattern_weights,
                                      const unsigned int * freqs_indices,
                                      double * persite_lh);

void pll_core_root_likelihood_vec_4x4_sse(unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * clv,
                                          const unsigned int * scaler,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const unsigned int * freqs_indices,
                                          double * persite_lh);
#endif

/* functions in core_partials_avx.c */

#ifdef HAVE_AVX
void pll_core_create_lookup_avx(unsigned int states,
                                unsigned int rate_cats,
                                double * lookup,
                                const double * left_matrix,
                                const double * right_matrix,
                                const unsigned int * tipmap,
                                unsigned int tipmap_size);

void pll_core_create_lookup_4x4_avx(unsigned int rate_cats,
                                    double * lookup,
                                    const double * left_matrix,
                                    const double * right_matrix);

void pll_core_create_lookup_20x20_avx(unsigned int rate_cats,
                                      double * ttlookup,
                                      const double * left_matrix,
                                      const double * right_matrix,
                                      const unsigned int * tipmap,
                                      unsigned int tipmap_size);

void pll_core_update_partial_tt_avx(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const unsigned char * right_tipchars,
                                    const double * lookup,
                                    unsigned int tipstates_count,
                                    unsigned int attrib);

void pll_core_update_partial_tt_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const unsigned char * left_tipchars,
                                        const unsigned char * right_tipchars,
                                        const double * lookup,
                                        unsigned int attrib);

void pll_core_update_partial_ti_avx(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const unsigned char * left_tipchars,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * right_scaler,
                                    const unsigned int * tipmap,
                                    unsigned int tipmap_size,
                                    unsigned int attrib);

void pll_core_update_partial_ti_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const unsigned char * left_tipchar,
                                        const double * right_clv,
                                        const double * left_matrix,
                                        const double * right_matrix,
                                        const unsigned int * right_scaler,
                                        unsigned int attrib);

void pll_core_update_partial_ti_20x20_avx(unsigned int sites,
                                          unsigned int rate_cats,
                                          double * parent_clv,
                                          unsigned int * parent_scaler,
                                          const unsigned char * left_tipchar,
                                          const double * right_clv,
                                          const double * left_matrix,
                                          const double * right_matrix,
                                          const unsigned int * right_scaler,
                                          const unsigned int * tipmap,
                                          unsigned int tipmap_size,
                                          unsigned int attrib);

void pll_core_update_partial_ii_avx(unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_cats,
                                    double * parent_clv,
                                    unsigned int * parent_scaler,
                                    const double * left_clv,
                                    const double * right_clv,
                                    const double * left_matrix,
                                    const double * right_matrix,
                                    const unsigned int * left_scaler,
                                    const unsigned int * right_scaler,
                                    unsigned int attrib);

void pll_core_update_partial_ii_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double * parent_clv,
                                        unsigned int * parent_scaler,
                                        const double * left_clv,
                                        const double * right_clv,
                                        const double * left_matrix,
                                        const double * right_matrix,
                                        const unsigned int * left_scaler,
                                        const unsigned int * right_scaler,
                                        unsigned int attrib);

/* functions in core_likelihood_avx.c */


double pll_core_root_loglikelihood_avx(unsigned int states,
                                       unsigned int sites,
                                       unsigned int rate_cats,
                                       const double * clv,
                                       const unsigned int * scaler,
                                       double * const * frequencies,
                                       const double * rate_weights,
                                       const unsigned int * pattern_weights,
                                       const unsigned int * freqs_indices,
                                       double * persite_lnl);

double pll_core_root_loglikelihood_4x4_avx(unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * clv,
                                           const unsigned int * scaler,
                                           double * const * frequencies,
                                           const double * rate_weights,
                                           const unsigned int * pattern_weights,
                                           const unsigned int * freqs_indices,
                                           double * persite_lnl);

void pll_core_root_likelihood_vec_avx(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      const double * clv,
                                      const unsigned int * scaler,
                                      double * const * frequencies,
                                      const double * rate_weights,
                                      const unsigned int * pattern_weights,
                                      const unsigned int * freqs_indices,
                                      double * persite_lh);

void pll_core_root_likelihood_vec_4x4_avx(unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * clv,
                                          const unsigned int * scaler,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const unsigned int * freqs_indices,
                                          double * persite_lh);
#endif


#ifdef HAVE_AVX2

/* functions in core_partials_avx2.c */

void pll_core_update_partial_ti_avx2(unsigned int states,
                                     unsigned int sites,
                                     unsigned int rate_cats,
                                     double * parent_clv,
                                     unsigned int * parent_scaler,
                                     const unsigned char * left_tipchars,
                                     const double * right_clv,
                                     const double * left_matrix,
                                     const double * right_matrix,
                                     const unsigned int * right_scaler,
                                     const unsigned int * tipmap,
                                     unsigned int tipmap_size,
                                     unsigned int attrib);


void pll_core_update_partial_ti_20x20_avx2(unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchar,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * right_scaler,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size,
                                           unsigned int attrib);

void pll_core_update_partial_ii_avx2(unsigned int states,
                                     unsigned int sites,
                                     unsigned int rate_cats,
                                     double * parent_clv,
                                     unsigned int * parent_scaler,
                                     const double * left_clv,
                                     const double * right_clv,
                                     const double * left_matrix,
                                     const double * right_matrix,
                                     const unsigned int * left_scaler,
                                     const unsigned int * right_scaler,
                                     unsigned int attrib);

/* functions in core_likelihood_avx2.c */

double pll_core_root_loglikelihood_avx2(unsigned int states,
                                        unsigned int sites,
                                        unsigned int rate_cats,
                                        const double * clv,
                                        const unsigned int * scaler,
                                        double * const * frequencies,
                                        const double * rate_weights,
                                        const unsigned int * pattern_weights,
                                        const unsigned int * freqs_indices,
                                        double * persite_lnl);

void pll_core_root_likelihood_vec_avx2(unsigned int states,
                                       unsigned int sites,
                                       unsigned int rate_cats,
                                       const double * clv,
                                       const unsigned int * scaler,
                                       double * const * frequencies,
                                       const double * rate_weights,
                                       const unsigned int * pattern_weights,
                                       const unsigned int * freqs_indices,
                                       double * persite_lh);
#endif
