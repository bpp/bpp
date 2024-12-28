/*
    Copyright (C) 2016-2024 Tomas Flouri, Bruce Rannala and Ziheng Yang

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
#include <inttypes.h>
#include <ctype.h>
#include <pthread.h>


#define PROG_NAME "bpp"

#if defined(__x86_64__) || defined(_M_AMD64)

#define PROG_CPU "x86_64"
#ifdef _MSC_VER
#include <pmmintrin.h>
#include <immintrin.h>
#else
#include <x86intrin.h>
#endif

#elif __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#include <altivec.h>
#else
#error "Big endian ppc64 CPUs not supported"
#endif

#elif __aarch64__

#define PROG_CPU "aarch64"
#include <arm_neon.h>
#undef HAVE_SSE3
#undef HAVE_AVX
#undef HAVE_AVX2

#else

#error Unknown architecture (not ppc64le, aarch64 or x86_64)

#endif

#ifdef __APPLE__
#define PROG_OS "macos"
#include <sys/resource.h>
#include <sys/sysctl.h>

#elif __linux__
#define PROG_OS "linux"
#include <sys/resource.h>
#include <sys/sysinfo.h>

#elif _WIN32
#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>

#elif __FreeBSD__
#define PROG_OS "freebsd"
#include <sys/resource.h>
#include <sys/sysctl.h>

#elif __NetBSD__
#define PROG_OS "netbsd"
#include <sys/resource.h>

#else

#define PROG_OS "unknown"
#include <sys/sysinfo.h>
#include <sys/resource.h>

#endif

#define PROG_ARCH PROG_OS "_" PROG_CPU

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

/* debugging switches */

#if 0
#define DEBUG_THREADS
#define DEBUG_THREADS_COUNT 4
#endif

/* constants */


#define PLL_STRING(x) #x
#define PLL_C2S(x) PLL_STRING(x)


#define VERSION_MAJOR 4
#define VERSION_MINOR 8
#define VERSION_PATCH 3

#define PVER_SHA1 "079bb48a9dc5c766bcc28e7cf766e47f5b4edf08"

/* checkpoint version */
#define VERSION_CHKP 1

#define PROG_VERSION "v" PLL_C2S(VERSION_MAJOR) "." PLL_C2S(VERSION_MINOR) "." \
        PLL_C2S(VERSION_PATCH)

#define BPP_MAGIC_BYTES 4
#define BPP_MAGIC "BPPX"

#define BPP_FALSE 0
#define BPP_TRUE  1

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

#define BPP_SPECIES_PRIOR_MIN           0
#define BPP_SPECIES_PRIOR_LH            0
#define BPP_SPECIES_PRIOR_UNIFORM       1
#define BPP_SPECIES_PRIOR_SLH           2
#define BPP_SPECIES_PRIOR_SUNIFORM      3
#define BPP_SPECIES_PRIOR_MAX           3

#define BPP_CLOCK_MIN                   1
#define BPP_CLOCK_GLOBAL                1
#define BPP_CLOCK_IND                   2
#define BPP_CLOCK_CORR                  3
#define BPP_CLOCK_MAX                   3

#define BPP_HPATH_NONE                  0
#define BPP_HPATH_LEFT                  1
#define BPP_HPATH_RIGHT                 2

#define BPP_DATA_DNA                    0
#define BPP_DATA_AA                     1

#define BPP_DNA_MODEL_MIN               0
#define BPP_DNA_MODEL_DEFAULT           0
#define BPP_DNA_MODEL_JC69              0
#define BPP_DNA_MODEL_K80               1
#define BPP_DNA_MODEL_F81               2
#define BPP_DNA_MODEL_HKY               3
#define BPP_DNA_MODEL_T92               4
#define BPP_DNA_MODEL_TN93              5
#define BPP_DNA_MODEL_F84               6
#define BPP_DNA_MODEL_GTR               7
#define BPP_DNA_MODEL_MAX               7

#define BPP_DNA_MODEL_CUSTOM            8

#define BPP_AA_MODEL_MIN                9
#define BPP_AA_MODEL_DAYHOFF            9
#define BPP_AA_MODEL_LG                10
#define BPP_AA_MODEL_DCMUT             11
#define BPP_AA_MODEL_JTT               12
#define BPP_AA_MODEL_MTREV             13
#define BPP_AA_MODEL_WAG               14
#define BPP_AA_MODEL_RTREV             15
#define BPP_AA_MODEL_CPREV             16
#define BPP_AA_MODEL_VT                17
#define BPP_AA_MODEL_BLOSUM62          18
#define BPP_AA_MODEL_MTMAM             19
#define BPP_AA_MODEL_MTART             20
#define BPP_AA_MODEL_MTZOA             21
#define BPP_AA_MODEL_PMB               22
#define BPP_AA_MODEL_HIVB              23
#define BPP_AA_MODEL_HIVW              24
#define BPP_AA_MODEL_JTTDCMUT          25
#define BPP_AA_MODEL_FLU               26
#define BPP_AA_MODEL_STMTREV           27
#define BPP_AA_MODEL_MAX               27

/* update the following table in bpp.c if any of the numbers above change */
extern const char * global_model_strings[28];
extern const char * global_freqs_strings[28];

/* Note: the XXX_PRIOR_DIR is after PRIOR_MAX as it's not available as input,
   but instead BPP sets it if mubar is fixed to 1 */
#define BPP_LOCRATE_PRIOR_MIN           0
#define BPP_LOCRATE_PRIOR_GAMMADIR      0
#define BPP_LOCRATE_PRIOR_HIERARCHICAL  1
#define BPP_LOCRATE_PRIOR_NONE  	2
#define BPP_LOCRATE_PRIOR_MAX           1
#define BPP_LOCRATE_PRIOR_DIR           2

#define BPP_BRATE_PRIOR_MIN             0
#define BPP_BRATE_PRIOR_LOGNORMAL       0
#define BPP_BRATE_PRIOR_GAMMA           1
#define BPP_BRATE_PRIOR_MAX             1

#define BPP_TAU_PRIOR_MIN               1
#define BPP_TAU_PRIOR_GAMMA             1
#define BPP_TAU_PRIOR_INVGAMMA          2
#define BPP_TAU_PRIOR_MAX               2


#define BPP_THETA_PRIOR_MIN             1
#define BPP_THETA_PRIOR_INVGAMMA        1
#define BPP_THETA_PRIOR_GAMMA           2
#define BPP_THETA_PRIOR_MAX             2
#define BPP_THETA_PROP_MG_INVG          1
#define BPP_THETA_PROP_MG_GAMMA         2

#define BPP_THETA_MOVE_NONE            -1
#define BPP_THETA_MOVE_MIN              0
#define BPP_THETA_MOVE_SLIDE            0
#define BPP_THETA_MOVE_GIBBS            1
#define BPP_THETA_MOVE_MAX              1

#define BPP_PHI_MOVE_SLIDE              0
#define BPP_PHI_MOVE_GIBBS              1

#define BPP_MRATE_SLIDE                 0
#define BPP_MRATE_GIBBS                 1

#define BPP_LINKEDTHETA_NONE            0
#define BPP_LINKEDTHETA_ALL             1  /* model M0 */
#define BPP_LINKEDTHETA_INNER           2  /* model M1 */
#define BPP_LINKEDTHETA_MSCI            3  /* model M2 */
#define BPP_LINKEDTHETA_MSCM            4  /* model M3 */


#define BPP_LB_NONE                     0
#define BPP_LB_ZIGZAG                   1

#define BPP_PI  3.1415926535897932384626433832795

#define THREAD_WORK_GTAGE               1
#define THREAD_WORK_GTSPR               2
#define THREAD_WORK_TAU                 3
#define THREAD_WORK_TAU_MIG             4
#define THREAD_WORK_MIXING              5
#define THREAD_WORK_ALPHA               6
#define THREAD_WORK_RATES               7
#define THREAD_WORK_FREQS               8
#define THREAD_WORK_BRATE               9

#define BPP_MOVE_INDEX_MIN              0
#define BPP_MOVE_GTAGE_INDEX            0
#define BPP_MOVE_GTSPR_INDEX            1
#define BPP_MOVE_THETA_INDEX            2
#define BPP_MOVE_TAU_INDEX              3
#define BPP_MOVE_MIX_INDEX              4
#define BPP_MOVE_LRHT_INDEX             5
#define BPP_MOVE_PHI_INDEX              6
#define BPP_MOVE_FREQS_INDEX            7
#define BPP_MOVE_QRATES_INDEX           8
#define BPP_MOVE_ALPHA_INDEX            9
#define BPP_MOVE_MUBAR_INDEX            10
#define BPP_MOVE_NUBAR_INDEX            11
#define BPP_MOVE_MUI_INDEX              12
#define BPP_MOVE_NUI_INDEX              13
#define BPP_MOVE_BRANCHRATE_INDEX       14
#define BPP_MOVE_MRATE_INDEX            15
#define BPP_MOVE_MIGVR_INDEX            16
#define BPP_MOVE_INDEX_MAX              16

#define BPP_MSCIDEFS_TREE               1
#define BPP_MSCIDEFS_DEFINE             2
#define BPP_MSCIDEFS_HYBRID             3
#define BPP_MSCIDEFS_BIDIR              4
#define BPP_MSCIDEFS_SHOWBL             5

#define BPP_CONSTDEFS_CONSTRAINT        1
#define BPP_CONSTDEFS_DEFINE            2
#define BPP_CONSTDEFS_OUTGROUP          3

#define BPP_OUTGROUP_NONE               0
#define BPP_OUTGROUP_FULL               1
#define BPP_OUTGROUP_PARTIAL            2

/* migbuffer events due to migration */
#define EVENT_COAL        0
#define EVENT_MIG_SOURCE  1
#define EVENT_MIG_TARGET  2
#define EVENT_TAU         3
#define EVENT_SAMPLE      4


/* libpll related definitions */

#define PLL_ALIGNMENT_CPU               8
#define PLL_ALIGNMENT_NEON             16
#define PLL_ALIGNMENT_SSE              16
#define PLL_ALIGNMENT_AVX              32

#define PLL_ATTRIB_ARCH_CPU            0
#define PLL_ATTRIB_ARCH_SSE       (1 << 0)
#define PLL_ATTRIB_ARCH_AVX       (1 << 1)
#define PLL_ATTRIB_ARCH_AVX2      (1 << 2)
#define PLL_ATTRIB_ARCH_AVX512    (1 << 3)
#define PLL_ATTRIB_ARCH_NEON      (1 << 4)
#define PLL_ATTRIB_ARCH_MASK         0x1F

#define PLL_ATTRIB_PATTERN_TIP    (1 << 5)

#define PLL_ATTRIB_RATE_SCALERS   (1 << 9)

#define PLL_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define PLL_SCALE_THRESHOLD (1.0/PLL_SCALE_FACTOR)
#define PLL_SCALE_FACTOR_SQRT 340282366920938463463374607431768211456.0 /* 2**128 */
#define PLL_SCALE_THRESHOLD_SQRT (1.0/PLL_SCALE_FACTOR_SQRT)
#define PLL_SCALE_BUFFER_NONE -1

#define PLL_MISC_EPSILON 1e-8

#define PLL_GAMMA_RATES_MEAN             0
#define PLL_GAMMA_RATES_MEDIAN           1

/* error codes */

#define ERROR_PHYLIP_SYNTAX            106
#define ERROR_PHYLIP_LONGSEQ           107
#define ERROR_PHYLIP_NONALIGNED        108
#define ERROR_PHYLIP_ILLEGALCHAR       109
#define ERROR_PHYLIP_UNPRINTABLECHAR   110
#define ERROR_PARSE_MORETHANEXPECTED   111
#define ERROR_PARSE_LESSTHANEXPECTED   112
#define ERROR_PARSE_INCORRECTFORMAT    113

/* available methods */

#define METHOD_00       0
#define METHOD_01       1
#define METHOD_10       2
#define METHOD_11       3

/* other */
#define MUTRATE_CONSTANT                0
#define MUTRATE_ESTIMATE                1
#define MUTRATE_FROMFILE                2
#define MUTRATE_ONLY                    3
#define HEREDITY_ESTIMATE               1
#define HEREDITY_FROMFILE               2

/* Handle memory (de)allocation of dlist_item_t within the miginfo_t structure.
   NOOP (no operation, do not allocate), ALLOC (allocate) FREE (deallocate) */
#define MI_DLI_NOOP                  0
#define MI_DLI_ALLOC                 1
#define MI_DLI_FREE                  2

/* ANSI color codes */

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_WHITE   "\x1b[1;37m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define BPP_ERROR  "[" ANSI_COLOR_RED "ERROR" ANSI_COLOR_RESET "]"

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

typedef struct migspec_s
{
  char * source;
  char * target;
  unsigned int si;
  unsigned int ti;

  double am; /* shape param for variable migration rates among loci */
  double alpha;
  double beta;
  double pseudo_a;
  double pseudo_b;
  long params;

  double M;
  double old_M;
  double * Mi;

  int mark;             /* used in rubberband to mark entry as affected */

  char * outfile;       /* used only to store rates when simulating  */
} migspec_t;

typedef struct migbuffer_s
{
  double time;
  long type;
  double * mrsum; /* per locus total migration rate (variable mig rates) */
  long active_count; /* number of active elements (1 or nloci) */
  struct snode_s ** donors;
  long donors_count;
} migbuffer_t;

typedef struct snode_s
{
  char * label;
  char * attrib;
  double length;
  double rate;          /* used for simulations (MCcoal) */
  double theta;
  double tau;
  double old_tau;
  double old_theta;
  struct snode_s * left;
  struct snode_s * right;
  struct snode_s * parent;
  unsigned int leaves;
  unsigned int * gene_leaves;
  int prop_tau;
  int * mark;

  void * data;

  /* clade support when delimiting species */
  double support;

  /* branch weight when inferring species tree */
  double weight;

  /* list of per-locus coalescent events */
  dlist_t ** coalevent;

  int * coal_count;

  /* branch rate (per locus)*/
  double * brate;

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

  /* TODO: This is a temporary fix for the Split (rj-MCMC) function
     to indicate the tip nodes that do not have a theta assigned to them,
     because the max number of lineages at every locus ending at the particular
     tip is 1. This should be implemented for inner nodes as well, but for
     keeping the compatibility with bpp4, all inner nodes have thetas. */
  long has_theta;

  /* no theta related variables */
  double t2h_sum;                   /* t2h sum for all loci */
  double hphi_sum;                  /* hphi sum for all loci */
  long coal_count_sum;              /* sum of coalencent events count */
  double notheta_logpr_contrib;     /* MSC density contribution from pop */
  double notheta_old_logpr_contrib; /* storage space for rollback */
  double * notheta_phi_contrib;     /* per-locus phi contribution for hybrid */
  double * notheta_old_phi_contrib; /* per-locus phi contribution for hybrid */

  /* constraints */
  long constraint;
  long constraint_lineno;

  /* introgression */
  long has_phi;                     /* has the phi parameter (1: yes, 0: no) */
  double hphi;                      /* genetic contribution */
  long htau;                        /* tau parameter (1: yes, 0: no) */
  struct snode_s * hybrid;          /* linked hybridization node */
  long * hx;                        /* sum of events count and seqin_count (per msa) */

  /* migration events associated with population across loci */
  long * migevent_count;
  migbuffer_t * migbuffer;
  long mb_mrsum_isarray;
  long mb_count;
  dlist_t ** mig_source;
  dlist_t ** mig_target;

  /* tip dating */
  int * epoch_count; /*Number of sampling epochs */
  double ** tip_date; /* Date at time of epoch*/
  int ** date_count;  /* Number of sequences sampled at the epoch */

  /* linked theta model */
  struct snode_s * linked_theta;

  /* independent theta step lengths */
  long theta_step_index;

  /* total coalescent waiting time for pop j at locus i */
  double * old_C2ji;
  double * C2ji;  /* total coal waiting time x2 in current pop (j) at locus i */

  long flag;
} snode_t;

typedef struct migevent_s
{
  double time;
  double old_time;
  snode_t * source;
  snode_t * target;

  dlist_item_t * di_src;
  dlist_item_t * di_tgt;
} migevent_t;

typedef struct miginfo_s
{
  long alloc_size;
  long count;

  migevent_t * me;
} miginfo_t;

typedef struct stree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;
  unsigned int hybrid_count;

  unsigned int locus_count;

  snode_t ** nodes;
  snode_t ** td;

  int ** pptable;

  snode_t * root;

  double root_age;

  /* no theta related variables */
  double notheta_logpr;      /* precomputed MSC density over all loci */
  double notheta_old_logpr;  /* storage space for rollback */
  double notheta_hfactor;    /* precomputed heredity term over all loci */
  double notheta_sfactor;    /* sum of per-locus (sequences-1)*sqrt(2) */

//  int mark;
  /* mean rate and variance across loci */
  double locusrate_mubar;
  double locusrate_nubar;
  double nui_sum;

  /* migration related elements */
  miginfo_t ** mi_tbuffer;
  long ** migcount_sum;      /* migrations across loci */
  double *** Wsji;
  double *** old_Wsji;

  double * u_constraint;
  double * l_constraint;

} stree_t;

typedef struct mutation_s
{
  char state; 
  int site;
  double time; 
  snode_t * pop;
  
} mutation_t;

typedef struct gnode_s
{
  char * label;
  double length;
  double time;
  double time_fixed;
  double old_time;
  struct gnode_s * left;
  struct gnode_s * right;
  struct gnode_s * parent;
  unsigned int leaves;

  void * data;

  snode_t * pop;
  snode_t * old_pop;

  /* pointer to the dlist item this node is wrapped into */
  dlist_item_t * coalevent;

  unsigned int node_index;
  unsigned int clv_valid;

  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;

  int mark;

  /* Array of flags describing the direction the lineage originating from the
     current node towards its parent takes on each  hybridization node.
     
     Assume x = stree->tip_count and y = stree->inner_count, then entry i of
     hpath corresponds to the hybridization node at stree->nodes[x+y+i]. The
     possible values are BPP_HPATH_NONE, BPP_HPATH_LEFT, BPP_HPATH_RIGHT */
  int * hpath;

  /* migration related structures */
  miginfo_t * mi;

} gnode_t;

typedef struct gtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  int original_index;
  long msa_index;

  gnode_t ** nodes;
  gnode_t * root;

  gnode_t ** travbuffer;

  /* auxiliary space for traversals */
  double logl;
  double logpr;
  double old_logpr;
  double old_logl;

  /* locus rate */
  double rate_mui;
  double rate_nui;
  double lnprior_rates;
  double old_lnprior_rates;

  /* migration counts between populations */
  long ** migcount;

  /* buffer to store feasible populations where migrants come from during simulation */
  snode_t ** migpops;

  /* per-locus affected snodes in new rubberband */
  snode_t ** rb_linked;
  long rb_lcount;

  /* per locus migration rate */
  double locus_mig_rate;

} gtree_t;


/* multifurcating tree structure */
typedef struct node_s
{
  char * label;
  char * attr;
  double length;
  double theta;
  struct node_s ** children;
  struct node_s * parent;
  int children_count;
  int mark;
  int leaves;
  long node_index;
  double tau;

  void * data;

} node_t;

typedef struct ntree_s
{
  int tip_count;
  int inner_count;
  node_t * root;
  node_t ** leaves;
  node_t ** inner;
} ntree_t;

typedef struct msa_s
{
  int count;
  int length;

  char ** sequence;
  char ** label;

  int amb_sites_count;
  int original_length;

  double * freqs;

  int dtype;
  int model;
  int original_index;

} msa_t;

/* used for creating FigTree.tre */
typedef struct nodepinfo_s 
{
  /* 95% HPD CI */
  double lo;
  double hi;

  /* mean age */
  double age;

  /* mean theta */
  double theta;
} nodepinfo_t; 


typedef struct partition_s
{
  long start;
  long end;
  long dtype;
  long model;
} partition_t;

typedef struct locus_s
{
  unsigned int tips;
  unsigned int clv_buffers;
  unsigned int states;
  unsigned int sites;
  unsigned int rate_matrices;
  unsigned int prob_matrices;
  unsigned int rate_cats;
  unsigned int qrates_param_count;
  unsigned int freqs_param_count;
  unsigned int scale_buffers;
  unsigned int attributes;
  unsigned int model;
  unsigned int dtype;

  /* vectorization parameters */
  size_t alignment;
  unsigned int states_padded;

  double rates_alpha;
  double ** clv;
  double ** pmatrix;
  double * rates;
  double * rate_weights;
  double ** subst_params;
  unsigned int ** scale_buffer;
  double ** frequencies;
  double * heredity;
  unsigned int * pattern_weights;
  unsigned int pattern_weights_sum;

  int * eigen_decomp_valid;
  double ** eigenvecs;
  double ** inv_eigenvecs;
  double ** eigenvals;

  /* index of frequency/qmatrix values set to use for computing the pmatrix
     for each rate category */
  unsigned int * param_indices;

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

  int original_index;

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

typedef struct mappingDate_s
{
  char * individual;
  double date;
  int lineno;
} mappingDate_t;

typedef struct mscidefs_s
{
  long type;
  long lineno;
  char * node1_1;
  char * node1_2;
  char * node2_1;
  char * node2_2;
  char * label1;
  char * label2;
  long has_tau1;
  long has_tau2;
  double phi1;
  double phi2;
} mscidefs_t;

typedef struct constdefs_s
{
  long type;
  long lineno;
  char * arg1;
  char * arg2;
} constdefs_t;

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

typedef struct thread_data_s
{
  /* contains common data that are passed to all threads, and variables that
     are filled with return values upon work completion (all threads finish) as
     part of a 'threads reduce' operation */

  /* arguments for all proposals */
  locus_t ** locus;
  gtree_t ** gtree;
  stree_t * stree;

  /* arguments for tau proposal */
  snode_t * snode;
  double oldage;
  double minage;
  double maxage;
  double minfactor;
  double maxfactor;
  snode_t ** affected;
  unsigned int paffected_count;

  /* arguments for mixing proposal */
  double c;

  /* return values for gene tree age/spr moves */
  long proposals;
  long accepted;

  /* return values for tau proposal */
  unsigned int count_above;
  unsigned int count_below;
  double logl_diff;
  double logpr_diff;

  /* return values for mixing proposal */
  double lnacceptance;

  /* rejection for tau rubberband with migration */
  long mig_reject;

} thread_data_t;

typedef struct thread_info_s
{
  pthread_t thread;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  
  /* type of work (0 no work) */
  volatile int work;

  /* thread parameters */
  long locus_first;
  long locus_count;

  thread_data_t td;

} thread_info_t;


/* macros */

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
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
#define xtruncate _chsize
#else
#define PLL_POPCOUNT __builtin_popcount
#define PLL_POPCOUNTL __builtin_popcountl
#define PLL_CTZ __builtin_ctz
#define xtruncate ftruncate
#endif

#define legacy_rndexp(index,mean) (-(mean)*log(legacy_rndu(index)))
#define QuantileGamma(prob,alpha,beta) QuantileChi2(prob,2.0*(alpha))/(2.0*(beta))


#define FLAG_AGE_UPDATE                 1 
#define FLAG_POP_UPDATE                 2
#define FLAG_BRANCH_UPDATE              4
#define FLAG_MISC                       8
#define FLAG_PRUNED                    16
#define FLAG_MIGRATE                   32
#define FLAG_PARTIAL_UPDATE           128
#define FLAG_RED_LEFT                 256
#define FLAG_RED_RIGHT                512
#define FLAG_SIMULATE                1024

#define SN_AFFECT                       (1 << 0)
#define BPP_SN_AFFECT_A7                (1 << 0)
#define BPP_SN_AFFECT_A8                <1 << 1)
#define BPP_SN_UPDATE_KC2               (1 << 2)
#define BPP_SN_UPDATE_PG                (1 << 3)
#define SN_MIG_UPDATE                   (1 << 4)

/* options */

extern long opt_alpha_cats;
extern long opt_arch;
extern long opt_basefreqs_fixed;
extern long opt_bfd_points;
extern long opt_burnin;
extern long opt_checkpoint;
extern long opt_checkpoint_current;
extern long opt_checkpoint_initial;
extern long opt_checkpoint_step;
extern long opt_cleandata;
extern long opt_clock;
extern long opt_comply;
extern long opt_constraint_count;
extern long opt_corepin;
extern long opt_debug;
extern long opt_debug_abort;
extern long opt_debug_br;
extern long opt_debug_bruce;
extern long opt_debug_counter;
extern long opt_debug_end;
extern long opt_debug_expand_count;
extern long opt_debug_expshr_count;
extern long opt_debug_full;
extern long opt_debug_gage;
extern long opt_debug_gspr;
extern long opt_debug_hs;
extern long opt_debug_migration;
extern long opt_debug_mix;
extern long opt_debug_mui;
extern long opt_debug_parser;
extern long opt_debug_rates;
extern long opt_debug_rj;
extern long opt_debug_shrink_count;
extern long opt_debug_sim;
extern long opt_debug_snl;
extern long opt_debug_sspr;
extern long opt_debug_start;
extern long opt_debug_tau;
extern long opt_debug_theta;
extern long opt_delimit_prior;
extern long opt_diploid_size;
extern long opt_est_heredity;
extern long opt_est_delimit;
extern long opt_est_locusrate;
extern long opt_est_geneflow;
extern long opt_est_mubar;
extern long opt_est_stree;
extern long opt_est_theta;
extern long opt_exp_imrb;
extern long opt_exp_randomize;
extern long opt_exp_theta;
extern long opt_exp_sim;
extern long opt_finetune_reset;
extern long opt_finetune_theta_count;
extern long opt_finetune_theta_mode;
extern long opt_help;
extern long opt_linkedtheta;
extern long opt_load_balance;
extern long opt_locusrate_prior;
extern long opt_locus_count;
extern long opt_locus_simlen;
extern long opt_max_species_count;
extern long opt_method;
extern long opt_migration;
extern long opt_migration_count;
extern long opt_mig_vrates_exist;
extern long opt_mix_theta_update;
extern long opt_mix_w_update;
extern long opt_model;
extern long opt_mrate_move;
extern long opt_msci;
extern long opt_onlysummary;
extern long opt_partition_count;
extern long opt_print_a1b1;
extern long opt_print_genetrees;
extern long opt_print_locus;
extern long opt_print_hscalars;
extern long opt_print_locusfile;
extern long opt_print_locusrate;
extern long opt_print_qmatrix;
extern long opt_print_rates;
extern long opt_print_samples;
extern long opt_pseudop_exist;
extern long opt_qrates_fixed;
extern long opt_quiet;
extern long opt_rate_prior;
extern long opt_rb_w_update;
extern long opt_rb_theta_update;
extern long opt_revolutionary_spr_method;
extern long opt_revolutionary_spr_debug;
extern long opt_rev_gspr;
extern long opt_rjmcmc_method;
extern long opt_samplefreq;
extern long opt_samples;
extern long opt_scaling;
extern long opt_seed;
extern long  opt_simulate_read_depth;
extern long opt_siterate_cats;
extern long opt_siterate_fixed;
extern long opt_tau_dist;
extern long opt_theta_gibbs_showall_eps;
extern long opt_theta_prior;
extern long opt_theta_prop;
extern long opt_threads;
extern long opt_threads_start;
extern long opt_threads_step;
extern long opt_usedata;
extern long opt_usedata_fix_gtree;
extern long opt_version;
extern long opt_extend;
extern double opt_alpha_alpha;
extern double opt_alpha_beta;
extern double opt_bfbeta;
extern double opt_finetune_alpha;
extern double opt_finetune_branchrate;
extern double opt_finetune_freqs;
extern double opt_finetune_gtage;
extern double opt_finetune_gtspr;
extern double opt_finetune_locusrate;
extern double opt_finetune_mix;
extern double opt_finetune_migrates;
extern double opt_finetune_mig_Mi;
extern double opt_finetune_mubar;
extern double opt_finetune_mui;
extern double opt_finetune_phi;
extern double opt_finetune_qrates;
extern double opt_finetune_nubar;
extern double opt_finetune_nui;
extern double opt_finetune_tau;
extern double opt_heredity_alpha;
extern double opt_heredity_beta;
extern double opt_snl_lambda_expand;
extern double opt_snl_lambda_shrink;
extern double opt_locusrate_mubar;      /* used only in simulation */
extern double opt_mig_alpha;
extern double opt_mig_beta;
extern double opt_mubar_alpha;
extern double opt_mubar_beta;
extern double opt_mui_alpha;
extern double opt_phi_alpha;
extern double opt_phi_beta;
extern double opt_phi_slide_prob;
extern double opt_prob_snl;            /* probability for SNL move, with 1 - opt_prob_snl for SPR */
extern double opt_prob_snl_shrink;
extern double opt_pseudo_alpha;
extern double opt_pseudo_beta;
extern double opt_rjmcmc_alpha;
extern double opt_rjmcmc_mean;
extern double opt_rjmcmc_epsilon;
extern double opt_simulate_base_err;
extern double opt_simulate_a_samples;
extern double opt_simulate_a_sites;
extern double opt_siterate_alpha;
extern double opt_siterate_beta;
extern double opt_tau_alpha;
extern double opt_tau_beta;
extern double opt_theta_alpha;
extern double opt_theta_beta;
extern double opt_theta_max;
extern double opt_theta_min;
extern double opt_theta_p;
extern double opt_theta_slide_prob;
extern double opt_theta_q;
extern double opt_vbar_alpha;
extern double opt_vbar_beta;
extern double opt_clock_vbar;
extern double opt_vi_alpha;
extern long * opt_diploid;
extern long * opt_finetune_theta_mask;
extern long * opt_print_locus_num;
extern long * opt_sp_seqcount;
extern char * opt_bfdriver;
extern char * cmdline;
extern char * opt_a1b1file;
extern char * opt_cfile;
extern char * opt_concatfile;
extern char * opt_constraintfile;
extern char * opt_datefile;
extern char * opt_heredity_filename;
extern char * opt_jobname;
extern char * opt_mapfile;
extern char * opt_mcmcfile;
extern char * opt_modelparafile;
extern char * opt_msafile;
extern char * opt_mscifile;
extern char * opt_locusrate_filename;
extern char * opt_partition_file;
extern char * opt_reorder;
extern char * opt_resume;
extern char * opt_seqDates;
extern char * opt_simulate;
extern char * opt_streenewick;
extern char * opt_treefile;
extern double * opt_basefreqs_params;
extern double * opt_finetune_theta;
extern double * opt_qrates_params;
extern migspec_t * opt_mig_specs;
extern long ** opt_migration_matrix;
extern long ** opt_mig_bitmatrix;
extern double ** opt_migration_events;
extern partition_t ** opt_partition_list;
extern int  opt_seqAncestral;

/* common data */

extern __THREAD int bpp_errno;
extern __THREAD char bpp_errmsg[200];

extern const unsigned int pll_map_nt[256];
extern const unsigned int pll_map_nt_tcag[256];
extern const unsigned int pll_map_aa[256];
extern const unsigned int pll_map_fasta[256];
extern const unsigned int pll_map_amb[256];
extern const unsigned int pll_map_validjc69[16];
extern const unsigned int bpp_tolower_table[256];
extern const unsigned int pll_map_nt_missing[256];
extern const unsigned int pll_map_aa_missing[256];

extern const double pll_aa_rates_dayhoff[190];
extern const double pll_aa_rates_lg[190];
extern const double pll_aa_rates_dcmut[190];
extern const double pll_aa_rates_jtt[190];
extern const double pll_aa_rates_mtrev[190];
extern const double pll_aa_rates_wag[190];
extern const double pll_aa_rates_rtrev[190];
extern const double pll_aa_rates_cprev[190];
extern const double pll_aa_rates_vt[190];
extern const double pll_aa_rates_blosum62[190];
extern const double pll_aa_rates_mtmam[190];
extern const double pll_aa_rates_mtart[190];
extern const double pll_aa_rates_mtzoa[190];
extern const double pll_aa_rates_pmb[190];
extern const double pll_aa_rates_hivb[190];
extern const double pll_aa_rates_hivw[190];
extern const double pll_aa_rates_jttdcmut[190];
extern const double pll_aa_rates_flu[190];
extern const double pll_aa_rates_stmtrev[190];

extern const double pll_aa_freqs_dayhoff[20];
extern const double pll_aa_freqs_lg[20];
extern const double pll_aa_freqs_dcmut[20];
extern const double pll_aa_freqs_jtt[20];
extern const double pll_aa_freqs_mtrev[20];
extern const double pll_aa_freqs_wag[20];
extern const double pll_aa_freqs_rtrev[20];
extern const double pll_aa_freqs_cprev[20];
extern const double pll_aa_freqs_vt[20];
extern const double pll_aa_freqs_blosum62[20];
extern const double pll_aa_freqs_mtmam[20];
extern const double pll_aa_freqs_mtart[20];
extern const double pll_aa_freqs_mtzoa[20];
extern const double pll_aa_freqs_pmb[20];
extern const double pll_aa_freqs_hivb[20];
extern const double pll_aa_freqs_hivw[20];
extern const double pll_aa_freqs_jttdcmut[20];
extern const double pll_aa_freqs_flu[20];
extern const double pll_aa_freqs_stmtrev[20];

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
extern long neon_present;

extern double ** global_sortbuffer_r;
extern migbuffer_t ** global_migbuffer_r;

/* pjumps */
extern double g_pj_gage;
extern double g_pj_gspr;
extern double g_pj_tau;
extern double g_pj_mix;
extern double g_pj_lrht;
extern double g_pj_phi_gibbs;
extern double g_pj_phi_slide;
extern double g_pj_freqs;
extern double g_pj_qmat;
extern double g_pj_alpha;
extern double g_pj_mubar;
extern double g_pj_nubar;
extern double g_pj_mui;
extern double g_pj_nui;
extern double g_pj_brate;
extern double g_pj_mrate;
extern double g_pj_migvr;
extern double * g_pj_theta_slide;
extern double * g_pj_theta_gibbs;

extern long g_pj_sspr;
extern long g_pj_ssnl;
extern double g_pj_rj;

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
int xtolower(int c);

/* functions in bpp.c */

void args_init(int argc, char ** argv);

void cmd_help(void);

void getentirecommandline(int argc, char * argv[]);

void fillheader(void);

void show_header(void);

void cmd_ml(void);

/* functions in phylip.c */

phylip_t * phylip_open(const char * filename,
                       const unsigned int * map);

int phylip_rewind(phylip_t * fd);

void phylip_close(phylip_t * fd);

msa_t * phylip_parse_interleaved(phylip_t * fd);

msa_t * phylip_parse_sequential(phylip_t * fd);

msa_t ** phylip_parse_multisequential(phylip_t * fd, long * count);

/* functions in rtree.c */

void stree_show_ascii(const snode_t * root, int options);

char * stree_export_newick(const snode_t * root, char * (*cb_serialize)(const snode_t *));

char* msci_export_newick(const snode_t* root, char* (*cb_serialize)(const snode_t*));

int stree_traverse(snode_t * root,
                   int traversal,
                   int (*cbtrav)(snode_t *),
                   snode_t ** outbuffer,
                   unsigned int * trav_size);

char * cb_stree_print_tau(const snode_t * node);
char * cb_stree_print_none(const snode_t * node);
char * cb_stree_print_node_tau(const snode_t * node);
char * cb_gtree_print_none(const gnode_t * node);
char * cb_gtree_print_age(const gnode_t * node);
char * cb_gtree_print_node_age(const gnode_t * node);
char * cb_gtree_print_length(const gnode_t* node);

/* functions in stree.c */
stree_t ** stree_tipstring_nodes(stree_t * root,
                                 char * tipstring,
                                 unsigned int * tiplist_count);

hashtable_t * species_hash(stree_t * tree);

hashtable_t * maplist_hash(list_t * maplist, hashtable_t * sht);

void stree_propose_theta(gtree_t ** gtree,
                         locus_t ** locus,
                         stree_t * stree,
                         double * acceptvec_gibbs,
                         double * acceptvec_slide,
                         long * acceptvec_movetype);

hashtable_t * datelist_hash(list_t * datelist);

double stree_propose_tau(gtree_t ** gtree, stree_t * stree, locus_t ** loci);
double stree_propose_tau_mig(stree_t ** streeptr,
                             gtree_t *** gtreeptr,
                             stree_t ** scloneptr,
                             gtree_t *** gcloneptr,
                             locus_t ** loci);

void stree_propose_phi(stree_t * stree,
                       gtree_t ** gtree,
                       long * acceptvec,
                       long * movecount,
                       long mcmc_step,
                       FILE * fp_a1b1);

void stree_fini(void);

void stree_rootdist(stree_t * stree,
                    list_t * maplist,
                    msa_t ** msalist,
                    unsigned int ** weights);

void fill_seqin_counts(stree_t * stree, gtree_t * gtree, int msa_index);

void reset_gene_leaves_count(stree_t * stree, gtree_t ** gtree);

void stree_reset_pptable(stree_t * stree);

long stree_propose_spr(stree_t ** streeptr,
                       gtree_t *** gtree_list_ptr,
                       stree_t ** scloneptr,
                       gtree_t *** gclonesptr,
                       locus_t ** loci);
long stree_propose_stree_snl(stree_t ** streeptr,
                             gtree_t *** gtree_list_ptr,
                             stree_t ** scloneptr,
                             gtree_t *** gclonesptr,
                             locus_t ** loci);

gtree_t * gtree_clone_init(gtree_t * gtree, stree_t * stree);

stree_t * stree_clone_init(stree_t * stree);

void stree_label(stree_t * stree);

void stree_show_pptable(stree_t * stree, int show_taus_and_thetas);

void stree_init(stree_t * stree,
                msa_t ** msa,
                list_t * maplist,
                int msa_count,
		int * tau_ctl,
                FILE * fp_out);

void stree_init_pptable(stree_t * stree);

void stree_alloc_internals(stree_t * stree,
                           long * locus_seqcount,
                           unsigned int gtree_inner_sum,
                           long msa_count);

int node_is_bidirection(const snode_t * node);
int node_is_mirror(const snode_t * node);
int node_is_hybridization(const snode_t * node);

void print_network_table(stree_t * stree, FILE * fp);
void debug_print_network_node_attribs(stree_t * stree);

void propose_tau_update_gtrees(locus_t ** loci,
                               gtree_t ** gtree,
                               stree_t * stree,
                               snode_t * snode,
                               double oldage,
                               double minage,
                               double maxage,
                               double minfactor,
                               double maxfactor,
                               long locus_start,
                               long locus_count,
                               snode_t ** affected,
                               unsigned int paffected_count,
                               unsigned int * ret_count_above,
                               unsigned int * ret_count_below,
                               double * ret_logl_diff,
                               double * ret_logpr_diff,
                               long thread_index);
void propose_tau_update_gtrees_mig(locus_t ** loci,
                                   gtree_t ** gtree,
                                   stree_t * stree,
                                   snode_t * snode,
                                   double oldage,
                                   double minage,
                                   double maxage,
                                   double minfactor,
                                   double maxfactor,
                                   long locus_start,
                                   long locus_count,
                                   snode_t ** affected,
                                   unsigned int paffected_count,
                                   long * ret_mig_reject,
                                   unsigned int * ret_count_above,
                                   unsigned int * ret_count_below,
                                   double * ret_logl_diff,
                                   double * ret_logpr_diff,
                                   long thread_index);

double lnprior_rates(gtree_t * gtree, stree_t * stree, long msa_index);

void stree_reset_leaves(stree_t * stree);

double prop_migrates(stree_t * stree, gtree_t ** gtree, locus_t ** locus);

double prop_mig_vrates(stree_t * stree, gtree_t ** gtree, locus_t ** locus);

void stree_update_mig_subpops(stree_t * stree, long msa_index);
void stree_update_mig_subpops_single(snode_t * x,
                                     snode_t * y,
                                     double oldM);
void stree_update_mig_subpops_single_vrates(snode_t * x,
                                            snode_t * y,
                                            long msa_index,
                                            double oldMi);

long migration_valid(snode_t * from, snode_t * to);

long stree_migration_rj(gtree_t *** gtreeptr,
                            gtree_t *** gcloneptr,
                            stree_t ** streeptr,
                            stree_t ** scloneptr,
                            locus_t ** locus);
long stree_migration_flip_wrapper(gtree_t *** gtreeptr,
                                  gtree_t *** gcloneptr,
                                  stree_t ** streeptr,
                                  stree_t ** scloneptr,
                                  locus_t ** locus);
                        

void stree_init_tau_recursive_constraint(stree_t * stree,
                                     snode_t * node,
                                     double prop,
                                     long thread_index,
                                     double *u_constraint,
                                     double *l_constraint);


double prop_tipDate_muGtree(gtree_t ** gtree,
                          stree_t * stree,
                          locus_t ** locus,
                          long thread_index); 

double find_maxMuGtree(stree_t * stree);
double prop_mu_updateCoal(gtree_t * gtree, stree_t * stree, double rateMultiplier, double new_mui);
void reset_mu_coal(gtree_t * gtree);

int get_gamma_conditional_approx(double a, double b, long k, double T,
                                 double * a1, double * b1);


/* functions in arch.c */

uint64_t arch_get_memused(void);

uint64_t arch_get_memtotal(void);

long arch_get_cores(void);

/* functions in msa.c */

void msa_print_phylip(FILE * fp,
                      msa_t ** msa,
                      long count,
                      unsigned int ** weights);

void msa_destroy(msa_t * msa);

int msa_remove_ambiguous(msa_t * msa);

void msa_summary(FILE * fp, msa_t ** msa_list, int msa_count);

void msa_count_ambiguous_sites(msa_t * msa, const unsigned int * map);

int msa_remove_missing_sequences(msa_t * msa);

/* functions in mapping.c */

void maplist_print(list_t * map_list);

void map_dealloc(void * data);

void mapDate_dealloc(void * data);

hashtable_t * maplist_hash(list_t * maplist, hashtable_t * sht);


/* functions in list.c */

void list_append(list_t * list, void * data);

void list_prepend(list_t * list, void * data);

void list_clear(list_t * list, void (*cb_dealloc)(void *));

long list_reposition_tail(list_t * list, list_item_t * item);
long list_delitem(list_t * list, list_item_t * item, void (*cb_dealloc)(void *));

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

double legacy_rndu(long index);
double legacy_rnd_symmetrical(long index);
void legacy_init(void);
void legacy_fini(void);
double legacy_rndbeta(long index, double p, double q);
double legacy_rndgamma(long index, double a);
unsigned int get_legacy_rndu_status(long index);
void set_legacy_rndu_status(long index, unsigned int x);
void legacy_rnddirichlet(long index, double * output, double * alpha, long k);
long legacy_rndpoisson(long index, double m);
unsigned int * get_legacy_rndu_array(void);
void set_legacy_rndu_array(unsigned int * x);
double rndNormal(long index);
int MultiNomialAlias(long index, int n, int ncat, double* F, int* L, int* nobs);
int MultiNomialAliasSetTable(int ncat, double* prob, double* F, int* L);
long legacy_rndbinomial(long index, int n, double p);

/* functions in gamma.c */

int pll_compute_gamma_cats(double alpha,
                           double beta,
                           unsigned int categories,
                           double * output_rates,
                           int rates_mode);

/* functions in gtree.c */

void gtree_alloc_internals(gtree_t ** gtree,
                           long msa_count,
                           unsigned int stree_inner_count, 
			   int totEpochs);

gtree_t ** gtree_init(stree_t * stree,
                      msa_t ** msalist,
                      list_t * maplist,
		      list_t * datelist,
		      int tau_ctl,
                      int msa_count);
void gtree_simulate_init(stree_t * stree, list_t * maplist);
void gtree_simulate_fini(void);

char * gtree_export_newick(const gnode_t * root,
                           char * (*cb_serialize)(const gnode_t *));

char * gtree_export_migration(const gnode_t * root);

void gtree_destroy(gtree_t * tree, void (*cb_destroy)(void *));

int gtree_traverse(gnode_t * root,
                   int traversal,
                   int (*cbtrav)(gnode_t *),
                   gnode_t ** outbuffer,
                   unsigned int * trav_size);

void gtree_update_branch_lengths(gtree_t ** gtree_list, int msa_count);

double gtree_propose_ages_serial(locus_t ** locus,
                                 gtree_t ** gtree,
                                 stree_t * stree);

void gtree_propose_ages_parallel(locus_t ** locus,
                                 gtree_t ** gtree,
                                 stree_t * stree,
                                 long locus_start,
                                 long locus_count,
                                 long thread_index,
                                 long * p_proposal_count,
                                 long * p_accepted);

void gtree_propose_spr_parallel(locus_t ** locus,
                                gtree_t ** gtree,
                                stree_t * stree,
                                long locus_start,
                                long locus_count,
                                long thread_index,
                                long * p_proposal_count,
                                long * p_accepted);

void gtree_reset_leaves(gnode_t * node);

void gtree_fini();

double gtree_logprob_mig(stree_t * stree,
                         gtree_t * gtree,
                         double heredity,
                         long msa_index,
                         long thread_index);
double gtree_logprob(stree_t * stree,
                     double heredity,
                     long msa_index,
                     long thread_index);
void gtree_update_branchlengths(stree_t * stree, gtree_t * gtree);

//double gtree_logprob_notheta(stree_t * stree);
double gtree_update_logprob_contrib_mig(snode_t * snode,
                                        stree_t * stree,
                                        gtree_t * gtree,
                                        double heredity,
                                        long msa_index,
                                        long thread_index);
double gtree_update_logprob_contrib(snode_t * snode,
                                    double heredity,
                                    long msa_index,
                                    long thread_index);
void gtree_update_C2j(snode_t * snode,
                      double heredity,
                      long msa_index,
                      long thread_index);
double update_logpg_contrib(stree_t * stree, snode_t * snode);
void logprob_revert_C2j(snode_t * snode, long msa_index);
void logprob_revert_contribs(snode_t * snode);
double gtree_propose_spr_serial(locus_t ** locus,
                                gtree_t ** gtree,
                                stree_t * stree);

double reflect(double t, double minage, double maxage, long thread_index);
void gtree_return_partials(gnode_t * root,
                           gnode_t ** trav,
                           unsigned int * trav_size);
void unlink_event(gnode_t * node, int msa_index);

double prop_locusrate_and_heredity(gtree_t ** gtree,
                                   stree_t * stree,
                                   locus_t ** locus,
                                   long thread_index);

gtree_t * gtree_simulate(stree_t * stree, msa_t * msa, int msa_index, 
                        mappingDate_t ** tipDateArray,
                        int tipDateArrayLen, 
			int tau_ctl);

double prop_branch_rates_serial(gtree_t ** gtree,
                                stree_t * stree,
                                locus_t ** locus);
void prop_branch_rates_parallel(gtree_t ** gtree,
                                stree_t * stree,
                                locus_t ** locus,
                                long locus_start,
                                long locus_count,
                                long thread_index,
                                long * p_proposal_count,
                                long * p_accepted);

long prop_locusrate_mubar(stree_t * stree, gtree_t ** gtree);
long prop_locusrate_nubar(stree_t * stree, gtree_t ** gtree);

double prop_locusrate_nui(gtree_t ** gtree,
                          stree_t * stree,
                          locus_t ** locus,
                          long thread_index);

double prop_locusrate_mui(gtree_t ** gtree,
                          stree_t * stree,
                          locus_t ** locus,
                          long thread_index);

void update_tau_constraint_recursive_to_root(stree_t * stree,
                                                snode_t * node,
                                                double * constraint);

void update_tau_constraint_recursive_to_tip(stree_t * stree,
                                                snode_t * node,
                                                double * constraint);

double gtree_propose_migevent_ages_serial(locus_t ** locus,
                                          gtree_t ** gtree,
                                          stree_t * stree);

void migbuffer_check_and_realloc(long thread_index, size_t alloc_required);
int cb_migbuf_asctime(const void * x, const void * y);

/* functions in prop_mixing.c */

long proposal_mixing(gtree_t ** gtree, stree_t * stree, locus_t ** locus);
void prop_mixing_update_gtrees(locus_t ** locus,
                               gtree_t ** gtree,
                               stree_t * stree,
                               long locus_start,
                               long locus_count,
                               double c,
                               long thread_index,
                               double * ret_lnacceptance);
                               #if 0
                               double * ret_logpr);
                               #endif

/* functions in prop_rj.c */

long prop_split(gtree_t ** gtree,
                stree_t * stree,
                locus_t ** locus,
                double pr_split,
                long * param_count,
                long * ndspecies);

long prop_join(gtree_t ** gtree,
               stree_t * stree,
               locus_t ** locus,
               double pr_split,
               long * param_count,
               long * ndspecies);

void rj_init(gtree_t ** gtreelist, stree_t * stree, unsigned int count);

void rj_fini();

/* functions in prop_gamma.c */

double locus_propose_alpha_serial(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree);
void locus_propose_alpha_parallel(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree,
                                  long locus_start,
                                  long locus_count,
                                  long thread_index,
                                  long * p_proposal_count,
                                  long * p_accepted);

/* functions in locus.c */

locus_t * locus_create(unsigned int dtype,
                       unsigned int model,
                       unsigned int tips,
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

void pll_set_subst_params(locus_t * locus,
                          unsigned int params_index,
                          const double * params);

void pll_set_frequencies(locus_t * locus,
                         unsigned int freqs_index,
                         const double * frequencies);

void locus_set_frequencies_and_rates(locus_t * locus);
void pll_set_category_rates(locus_t * locus, const double * rates);
void locus_set_heredity_scalers(locus_t * locus, const double * heredity);

void locus_update_partials(locus_t * locus, gnode_t ** traversal, unsigned int count);

void locus_update_all_partials(locus_t * locus, gtree_t * gtree);

void pll_set_pattern_weights(locus_t * locus,
                             const unsigned int * pattern_weights);

void locus_update_matrices(locus_t * locus,
                           gtree_t * gtree,
                           gnode_t ** traversal,
                           stree_t * stree,
                           long msa_index,
                           unsigned int count);

void locus_update_all_matrices(locus_t * locus,
                               gtree_t * gtree,
                               stree_t * stree,
                               long msa_index);

double locus_root_loglikelihood(locus_t * locus,
                                gnode_t * root,
                                const unsigned int * freqs_indices,
                                double * persite_lnl);

double locus_propose_qrates_serial(stree_t * stree,
                                   locus_t ** locus,
                                   gtree_t ** gtree);
void locus_propose_qrates_parallel(stree_t * stree,
                                   locus_t ** locus,
                                   gtree_t ** gtree,
                                   long locus_start,
                                   long locus_count,
                                   long thread_index,
                                   long * p_proposal_count,
                                   long * p_accepted);

double locus_propose_freqs_serial(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree);
void locus_propose_freqs_parallel(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree,
                                  long locus_start,
                                  long locus_count,
                                  long thread_index,
                                  long * p_proposal_count,
                                  long * p_accepted);

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
                                               unsigned int ** wptr,
                                               int attrib);
/* functions in allfixed.c */

void allfixed_summary(FILE * fp_out, stree_t * stree);
double eff_ict(double * y, long n, double mean, double stdev, double * rho1);

/* functions in summary.c */

void bipartitions_init(char ** species, long species_count);

void bipartitions_update(stree_t * stree);

void summary_dealloc_hashtables(void);

void stree_summary(FILE * fp_out, char ** species_names, long species_count);

long getlinecount(const char * filename);

/* functions in summary11.c */

void mixed_summary(FILE * fp_out, unsigned int sp_count);

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
                                 list_t * datelist,
                                 unsigned int ** weights,
                                 int msa_count);

/* functions in dump.c */

int checkpoint_dump(stree_t * stree,
                    gtree_t ** gtree_list,
                    locus_t ** locus_list,
                    unsigned long curstep,
                    long ft_round,
                    long ndspecies,
                    long mcmc_offset,
                    long out_offset,
                    long a1b1_offset,
                    long * gtree_offset,
                    long * mig_offset,
                    long * rates_offset,
                    long * migcount_offset,
                    long dparam_count,
                    double * posterior,
                    double * pspecies,
                    long dmodels_count,
                    long ft_round_rj,
                    long ft_round_spr,
                    long ft_round_snl,
                    double mean_logl,
                    long * mean_mrate_row,
                    long * mean_mrate_col,
                    long * mean_mrate_round,
                    double * mean_mrate,
                    double * mean_tau,
                    double * mean_theta,
                    double * mean_phi,
                    long mean_mrate_count,
                    long mean_tau_count,
                    long mean_theta_count,
                    long mean_phi_count,
                    int prec_logpg,
                    int prec_logl, 
		    int * printLocusIndex);

/* functions in load.c */

int checkpoint_load(gtree_t *** gtreep,
                    locus_t *** locusp,
                    stree_t ** streep,
                    unsigned long * curstep,
                    long * ft_round,
                    long * ndspecies,
                    long * mcmc_offset,
                    long * out_offset,
                    long * a1b1_offset,
                    long ** gtree_offset,
                    long ** mig_offset,
                    long ** rates_offset,
                    long ** migcount_offset,
                    long * dparam_count,
                    double ** posterior,
                    double ** pspecies,
                    long * ft_round_rj,
                    long * ft_round_spr,
                    long * ft_round_snl,
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
                    int * prec_logl,
		    int ** ptr_printLocusIndex);

void checkpoint_truncate(const char * filename, long mcmc_offset);

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

void bpp_core_update_pmatrix(locus_t * locus,
                             gtree_t * gtree,
                             gnode_t ** traversal,
                             stree_t * stree,
                             long msa_index,
                             unsigned int count);

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

void pll_update_eigen(double * eigenvecs,
                      double * inv_eigenvecs,
                      double * eigenvals,
                      double * freqs,
                      double * subst_params,
                      long states,
                      long states_padded);

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

void pll_show_eigendecomp(const locus_t * locus, unsigned int float_precision);

/* functions in revolutionary.c */

void revolutionary_spr_tselect_logl(gnode_t * mnode,
                               gnode_t ** target_list,
                               long target_count,
                               locus_t * locus,
                               double * weights);

void rev_spr_tselect(gnode_t * mnode,
                     double t,
                     gnode_t ** targets,
                     unsigned int target_count,
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
long delimitation_getcurindex();

long delimit_getindex(stree_t * stree);
long delimit_getindexfromstring(char * model);

void delimit_setindex(long index);

void delimit_resetpriors(void);

void delimit_summary(FILE * fp, stree_t * stree);

double lnprior_species_model(stree_t * stree); /* TODO: Move function */

void partition_fast(long n);

double * getpriorA11(void);

/* functions in cfile.c */

void load_cfile(void);
int parsefile_doubles(const char * filename,
                      long n,
                      double * outbuffer,
                      long * errcontext);
long parse_printlocus(const char * line, long * lcount);

#ifdef HAVE_NEON

/* functions in core_partials_neon.c */

void pll_core_create_lookup_neon(unsigned int states,
                                 unsigned int rate_cats,
                                 double * ttlookup,
                                 const double * left_matrix,
                                 const double * right_matrix,
                                 const unsigned int * tipmap,
                                 unsigned int tipmap_size);

void pll_core_create_lookup_4x4_neon(unsigned int rate_cats,
                                     double * lookup,
                                     const double * left_matrix,
                                     const double * right_matrix);

void pll_core_update_partial_tt_neon(unsigned int states,
                                     unsigned int sites,
                                     unsigned int rate_cats,
                                     double * parent_clv,
                                     unsigned int * parent_scaler,
                                     const unsigned char * left_tipchars,
                                     const unsigned char * right_tipchars,
                                     const double * lookup,
                                     unsigned int tipstates_count,
                                     unsigned int attrib);

void pll_core_update_partial_tt_4x4_neon(unsigned int sites,
                                         unsigned int rate_cats,
                                         double * parent_clv,
                                         unsigned int * parent_scaler,
                                         const unsigned char * left_tipchars,
                                         const unsigned char * right_tipchars,
                                         const double * lookup,
                                         unsigned int attrib);

void pll_core_update_partial_ti_neon(unsigned int states,
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


void pll_core_update_partial_ti_4x4_neon(unsigned int sites,
                                         unsigned int rate_cats,
                                         double * parent_clv,
                                         unsigned int * parent_scaler,
                                         const unsigned char * left_tipchar,
                                         const double * right_clv,
                                         const double * left_matrix,
                                         const double * right_matrix,
                                         const unsigned int * right_scaler,
                                         unsigned int attrib);

void pll_core_update_partial_ii_neon(unsigned int states,
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

void pll_core_update_partial_ii_4x4_neon(unsigned int sites,
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
/* functions in core_likelihood_neon.c */


double pll_core_root_loglikelihood_neon(unsigned int states,
                                        unsigned int sites,
                                        unsigned int rate_cats,
                                        const double * clv,
                                        const unsigned int * scaler,
                                        double * const * frequencies,
                                        const double * rate_weights,
                                        const unsigned int * pattern_weights,
                                        const unsigned int * freqs_indices,
                                        double * persite_lnl);

double pll_core_root_loglikelihood_4x4_neon(unsigned int sites,
                                            unsigned int rate_cats,
                                            const double * clv,
                                            const unsigned int * scaler,
                                            double * const * frequencies,
                                            const double * rate_weights,
                                            const unsigned int * pattern_weights,
                                            const unsigned int * freqs_indices,
                                            double * persite_lnl);

void pll_core_root_likelihood_vec_neon(unsigned int states,
                                       unsigned int sites,
                                       unsigned int rate_cats,
                                       const double * clv,
                                       const unsigned int * scaler,
                                       double * const * frequencies,
                                       const double * rate_weights,
                                       const unsigned int * pattern_weights,
                                       const unsigned int * freqs_indices,
                                       double * persite_lh);

void pll_core_root_likelihood_vec_4x4_neon(unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * clv,
                                           const unsigned int * scaler,
                                           double * const * frequencies,
                                           const double * rate_weights,
                                           const unsigned int * pattern_weights,
                                           const unsigned int * freqs_indices,
                                           double * persite_lh);
#endif

#ifdef HAVE_SSE3

/* functions in core_partials_sse.c */

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

/* functions in cfile_sim.c */

void load_cfile_sim(void);

/* functions in simulate.c */

void cmd_simulate(void);
double QuantileChi2(double prob, double v);

/* functions in threads.c */

long * threads_load_balance(msa_t ** msa_list);
void threads_lb_stats(locus_t ** locus, FILE * fp_out);
void threads_init(void);
void threads_wakeup(int work_type, thread_data_t * tp);
void threads_exit(void);
void threads_pin_master(void);
thread_info_t * threads_ti(void);
void threads_set_ti(thread_info_t * tip);

/* functions in treeparse.c */

int ntree_check_rbinary(ntree_t * tree);
int ntree_check_ubinary(ntree_t * tree);
void stree_destroy(stree_t * tree, void (*cb_destroy)(void *));
void ntree_destroy(ntree_t * tree, void (*cb_destroy)(void *));
char * ntree_export_newick(ntree_t * tree, long print_bl);
stree_t * bpp_parse_newick_string(const char * line);
ntree_t * bpp_parse_newick_string_ntree(const char * line);
ntree_t * ntree_wraptree(node_t * root, int tip_count, int inner_count);
stree_t * stree_from_ntree(ntree_t * ntree);
mappingDate_t ** prepTipDatesInfer(stree_t * stree, list_t *  dateList, 
		 double mu_bar);

/* functions in parsemap.c */

list_t * parse_mapfile(const char * mapfile);
list_t * parse_date_mapfile(const char * mapfile);

/* functions in msci_gen.c */

void cmd_msci_create(void);

/* functions in constraint.c */
void parse_and_set_constraints(stree_t * stree, FILE * fp_out);
void cmd_comply();

/* functions in debug.c */
void debug_linked_notheta(stree_t * stree,
                          gtree_t * gtree,
                          double bpp_logPG,
                          const char * move,
                          int only_pg,
                          long lmodel);
void debug_linked_notheta3(stree_t * stree, gtree_t ** gtree, double bpp_logPG, const char * move, int only_pg, long lmodel);
void debug_print_wsji(stree_t * stree,
                      int prec,
                      const char * prefix_msg,
                      FILE * fp);
void debug_print_migmatrix(stree_t * stree);
void debug_print_migrations(stree_t * stree);
void debug_print_bitmatrix(stree_t * stree);
void debug_print_gtree(gtree_t * gtree);
void debug_print_gtree_detailed(gtree_t * gtree);
void debug_print_stree(stree_t * stree);
void debug_check_relations(stree_t * stree, gtree_t * gtree, long msa_index);
void debug_check_leaves(gtree_t ** gtree);
void debug_validate_logpg(stree_t * stree,
                          gtree_t ** gtree,
                          locus_t ** locus,
                          const char * move);
void debug_snl_stage1(stree_t * stree,
                      gtree_t ** gtree_list,
                      snode_t * y,
                      snode_t * c,
                      snode_t * a,
                      long movetype,
                      int downwards,
                      double tauy,
                      double newtauy);
void debug_snl_stage2(stree_t * stree,
                      gtree_t ** gtree,
                      double logpr_notheta,
                      double r,
                      double lnacceptance);
void debug_bruce(stree_t * stree,
                 gtree_t ** gtree,
                 const char * move,
                 long iter,
                 FILE * fp_out);
void debug_migration_internals(stree_t * stree, gtree_t * gtree, long msa_index);
void debug_consistency(stree_t * stree, gtree_t ** gtree_list, const char * msg);

/* functions in miginfo.c */
miginfo_t * miginfo_create(size_t alloc_size, int alloc_dlist);
miginfo_t * miginfo_create_default();
void miginfo_destroy(miginfo_t * mi, long msa_index, int dealloc_dlist);
miginfo_t * miginfo_extend(miginfo_t * mi, size_t units);
void miginfo_append(miginfo_t ** miptr, snode_t * s, snode_t * t, double time, long msa_index);
void miginfo_clean(miginfo_t * mi, long msa_index);
void miginfo_clone(miginfo_t * mi,
                   miginfo_t ** clonebuf,
                   stree_t * clone_stree,
                   long msa_index);
void miginfo_move(migevent_t * me, miginfo_t ** miptr);
void migevent_link(migevent_t * me, long msa_index);
void migevent_unlink(migevent_t * me, long msa_index);
void miginfo_check_and_extend(miginfo_t ** miptr, size_t extra);

/* functions in lswitch.c */
void lswitch(stree_t * stree,
             const char * header,
             double ** matrix,
             long col_count,
             FILE * fp_outfile);

/* functions in ming2.c */
int ming2(FILE *fout, double *f, double(*fun)(double x[], int n),
          int(*dfun)(double x[], double *f, double dx[], int n),
          double x[], double xb[][2], double space[], double e, int n);

/* functions in bfdriver.c */
void cmd_bfdriver();

/* functions in visual.c */
void stree_export_pdf(const stree_t * stree);

/* functions in a1b1.c */
void conditional_to_marginal(double * ai_full,
                             double * bi_full,
                             long nsamples,
                             long nbins,
                             long cond_dist,
                             double tail,
                             double * out_mean,
                             double * out_sd,
                             double * out_et025,
                             double * out_et975,
                             double * out_hpd025,
                             double * out_hpd975,
                             double * out_c,
                             double * out_effu,
                             double * out_effy);

void conditional_to_marginal_M(double * ai1_full,
                               double * bi1_full,
                               double * ai2_full,
                               double * bi2_full,
                               int nsamples,
                               int nbins,
                               int cond_dist1,
                               double tail,
                               double * out_wmean,
                               double * out_wsd,
                               double * out_mean,
                               double * out_sd,
                               double * out_et025,
                               double * out_et975,
                               double * out_hpd025,
                               double * out_hpd975,
                               double * out_c,
                               double * out_effu,
                               double * out_effy);
