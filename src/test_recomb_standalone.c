/*
    Standalone unit tests for recombination module
    These tests do not require linking with the full BPP codebase.
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============== Minimal type definitions ============== */

typedef struct gnode_s gnode_t;
typedef struct gtree_s gtree_t;
typedef struct snode_s snode_t;
typedef struct stree_s stree_t;
typedef struct locus_s locus_t;

struct gnode_s
{
  unsigned int node_index;
  double time;
  gnode_t * parent;
  gnode_t * left;
  gnode_t * right;
};

struct gtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;
  gnode_t ** nodes;
  gnode_t * root;
  unsigned int msa_index;
};

struct snode_s
{
  double theta;
  double tau;
};

struct stree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  snode_t ** nodes;
};

struct locus_s
{
  unsigned int sites;
  int has_recombination;
  double * heredity;
  unsigned int * pattern_weights;
  unsigned int * param_indices;
  void * arg;  /* Will be cast to arg_t */
};

/* Recombination structures */
typedef struct breakpoint_s
{
  unsigned int position;
  unsigned int lineage;
  double recomb_time;
  unsigned int target_pop;
  double coal_time;
  unsigned int target_node;
} breakpoint_t;

typedef struct block_s
{
  unsigned int start;
  unsigned int end;
  gtree_t * local_tree;
  double logl;
} block_t;

typedef struct arg_s
{
  unsigned int num_breakpoints;
  breakpoint_t * breakpoints;
  unsigned int max_breakpoints;
  unsigned int num_blocks;
  block_t * blocks;
  gtree_t * base_tree;
  double recomb_rate;
  double log_recomb_prior;
} arg_t;

/* ============== Stub utility functions ============== */

static void fatal(const char * format, ...)
{
  fprintf(stderr, "FATAL: %s\n", format);
  exit(1);
}

static void * xmalloc(size_t size)
{
  void * ptr = malloc(size);
  if (!ptr) fatal("Out of memory");
  return ptr;
}

static void * xrealloc(void * ptr, size_t size)
{
  void * new_ptr = realloc(ptr, size);
  if (!new_ptr) fatal("Out of memory");
  return new_ptr;
}

/* ============== Recombination functions (from recomb.c) ============== */

static void local_tree_destroy(gtree_t * tree)
{
  if (tree) free(tree);
}

arg_t * arg_create(unsigned int max_bp, gtree_t * base, unsigned int sites)
{
  arg_t * arg;

  arg = (arg_t *)xmalloc(sizeof(arg_t));

  arg->num_breakpoints = 0;
  arg->max_breakpoints = max_bp;
  arg->breakpoints = (breakpoint_t *)xmalloc(max_bp * sizeof(breakpoint_t));

  arg->num_blocks = 1;
  arg->blocks = (block_t *)xmalloc((max_bp + 1) * sizeof(block_t));

  arg->blocks[0].start = 0;
  arg->blocks[0].end = sites;
  arg->blocks[0].local_tree = base;
  arg->blocks[0].logl = 0.0;

  arg->base_tree = base;
  arg->recomb_rate = 0.001;
  arg->log_recomb_prior = 0.0;

  return arg;
}

void arg_destroy(arg_t * arg)
{
  unsigned int i;

  if (!arg) return;

  for (i = 1; i < arg->num_blocks; i++)
  {
    if (arg->blocks[i].local_tree &&
        arg->blocks[i].local_tree != arg->base_tree)
    {
      local_tree_destroy(arg->blocks[i].local_tree);
    }
  }

  if (arg->breakpoints) free(arg->breakpoints);
  if (arg->blocks) free(arg->blocks);
  free(arg);
}

void arg_reset(arg_t * arg, unsigned int sites)
{
  unsigned int i;

  if (!arg) return;

  for (i = 1; i < arg->num_blocks; i++)
  {
    if (arg->blocks[i].local_tree &&
        arg->blocks[i].local_tree != arg->base_tree)
    {
      local_tree_destroy(arg->blocks[i].local_tree);
    }
  }

  arg->num_breakpoints = 0;
  arg->num_blocks = 1;

  arg->blocks[0].start = 0;
  if (sites > 0)
    arg->blocks[0].end = sites;
  arg->blocks[0].local_tree = arg->base_tree;
  arg->blocks[0].logl = 0.0;
}

double arg_log_prior(arg_t * arg, stree_t * stree, locus_t * locus)
{
  double logpr = 0.0;
  unsigned int L;
  double rho = arg->recomb_rate;
  unsigned int k = arg->num_breakpoints;

  if (locus)
    L = locus->sites;
  else
    L = arg->blocks[0].end;

  double lambda = rho * (double)L;
  if (k > 0)
    logpr += (double)k * log(lambda) - lambda - lgamma((double)k + 1.0);
  else
    logpr -= lambda;

  if (k > 0)
    logpr -= (double)k * log((double)L);

  arg->log_recomb_prior = logpr;
  return logpr;
}

void insert_breakpoint(arg_t * arg, breakpoint_t * bp)
{
  unsigned int i, j;

  if (arg->num_breakpoints >= arg->max_breakpoints)
  {
    arg->max_breakpoints *= 2;
    arg->breakpoints = (breakpoint_t *)xrealloc(arg->breakpoints,
                        arg->max_breakpoints * sizeof(breakpoint_t));
    arg->blocks = (block_t *)xrealloc(arg->blocks,
                   (arg->max_breakpoints + 1) * sizeof(block_t));
  }

  for (i = 0; i < arg->num_breakpoints; i++)
  {
    if (arg->breakpoints[i].position > bp->position)
      break;
  }

  for (j = arg->num_breakpoints; j > i; j--)
  {
    arg->breakpoints[j] = arg->breakpoints[j - 1];
  }

  arg->breakpoints[i] = *bp;
  arg->num_breakpoints++;
}

void remove_breakpoint(arg_t * arg, unsigned int idx)
{
  unsigned int i;

  if (idx >= arg->num_breakpoints)
    return;

  for (i = idx; i < arg->num_breakpoints - 1; i++)
  {
    arg->breakpoints[i] = arg->breakpoints[i + 1];
  }

  arg->num_breakpoints--;
}

void update_blocks(arg_t * arg, locus_t * locus)
{
  unsigned int i;
  unsigned int sites;

  if (locus)
    sites = locus->sites;
  else
    sites = arg->blocks[arg->num_blocks - 1].end;

  arg->num_blocks = arg->num_breakpoints + 1;

  arg->blocks[0].start = 0;

  for (i = 0; i < arg->num_breakpoints; i++)
  {
    arg->blocks[i].end = arg->breakpoints[i].position;
    arg->blocks[i + 1].start = arg->breakpoints[i].position;
  }

  arg->blocks[arg->num_blocks - 1].end = sites;
  arg->blocks[0].local_tree = arg->base_tree;
}

void merge_blocks(arg_t * arg, unsigned int idx)
{
  unsigned int i;

  if (idx >= arg->num_blocks - 1)
    return;

  if (arg->blocks[idx + 1].local_tree &&
      arg->blocks[idx + 1].local_tree != arg->base_tree)
  {
    local_tree_destroy(arg->blocks[idx + 1].local_tree);
  }

  arg->blocks[idx].end = arg->blocks[idx + 1].end;

  for (i = idx + 1; i < arg->num_blocks - 1; i++)
  {
    arg->blocks[i] = arg->blocks[i + 1];
  }

  arg->num_blocks--;
}

/* ============== Test framework ============== */

#define TEST_ASSERT(cond, msg) do { \
  if (!(cond)) { \
    fprintf(stderr, "FAIL: %s (line %d): %s\n", __func__, __LINE__, msg); \
    return 1; \
  } \
} while(0)

#define TEST_PASS() do { \
  fprintf(stderr, "PASS: %s\n", __func__); \
  return 0; \
} while(0)

static int tests_run = 0;
static int tests_passed = 0;

/* ============== Test cases ============== */

static int test_arg_create_destroy(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;

  arg = arg_create(max_bp, NULL, sites);

  TEST_ASSERT(arg != NULL, "arg_create returned NULL");
  TEST_ASSERT(arg->max_breakpoints == max_bp, "max_breakpoints not set correctly");
  TEST_ASSERT(arg->num_breakpoints == 0, "initial num_breakpoints should be 0");
  TEST_ASSERT(arg->num_blocks == 1, "initial num_blocks should be 1");
  TEST_ASSERT(arg->breakpoints != NULL, "breakpoints array not allocated");
  TEST_ASSERT(arg->blocks != NULL, "blocks array not allocated");
  TEST_ASSERT(arg->blocks[0].start == 0, "first block start should be 0");
  TEST_ASSERT(arg->blocks[0].end == sites, "first block end should be sites");
  TEST_ASSERT(arg->recomb_rate >= 0, "recomb_rate should be non-negative");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_arg_reset(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  arg->num_breakpoints = 3;
  arg->num_blocks = 4;
  arg->breakpoints[0].position = 100;
  arg->breakpoints[1].position = 500;
  arg->breakpoints[2].position = 800;

  arg_reset(arg, sites);

  TEST_ASSERT(arg->num_breakpoints == 0, "num_breakpoints should be 0 after reset");
  TEST_ASSERT(arg->num_blocks == 1, "num_blocks should be 1 after reset");
  TEST_ASSERT(arg->blocks[0].start == 0, "block start should be 0 after reset");
  TEST_ASSERT(arg->blocks[0].end == sites, "block end should be sites after reset");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_breakpoint_insertion(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  breakpoint_t bp;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  bp.position = 500;
  bp.lineage = 0;
  bp.recomb_time = 0.001;
  bp.target_pop = 0;
  bp.coal_time = 0.002;
  bp.target_node = 1;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 1, "should have 1 breakpoint");

  bp.position = 200;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 2, "should have 2 breakpoints");
  TEST_ASSERT(arg->breakpoints[0].position == 200, "first bp should be at 200");
  TEST_ASSERT(arg->breakpoints[1].position == 500, "second bp should be at 500");

  bp.position = 800;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 3, "should have 3 breakpoints");
  TEST_ASSERT(arg->breakpoints[2].position == 800, "third bp should be at 800");

  bp.position = 350;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 4, "should have 4 breakpoints");
  TEST_ASSERT(arg->breakpoints[0].position == 200, "positions should be sorted: 200");
  TEST_ASSERT(arg->breakpoints[1].position == 350, "positions should be sorted: 350");
  TEST_ASSERT(arg->breakpoints[2].position == 500, "positions should be sorted: 500");
  TEST_ASSERT(arg->breakpoints[3].position == 800, "positions should be sorted: 800");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_breakpoint_removal(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  breakpoint_t bp;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  bp.position = 200;
  bp.lineage = 0;
  bp.recomb_time = 0.001;
  bp.target_pop = 0;
  bp.coal_time = 0.002;
  bp.target_node = 1;
  insert_breakpoint(arg, &bp);

  bp.position = 500;
  insert_breakpoint(arg, &bp);

  bp.position = 800;
  insert_breakpoint(arg, &bp);

  update_blocks(arg, NULL);
  TEST_ASSERT(arg->num_blocks == 4, "should have 4 blocks");
  TEST_ASSERT(arg->blocks[0].start == 0, "block 0 start");
  TEST_ASSERT(arg->blocks[0].end == 200, "block 0 end");
  TEST_ASSERT(arg->blocks[1].start == 200, "block 1 start");
  TEST_ASSERT(arg->blocks[1].end == 500, "block 1 end");
  TEST_ASSERT(arg->blocks[2].start == 500, "block 2 start");
  TEST_ASSERT(arg->blocks[2].end == 800, "block 2 end");
  TEST_ASSERT(arg->blocks[3].start == 800, "block 3 start");
  TEST_ASSERT(arg->blocks[3].end == 1000, "block 3 end");

  remove_breakpoint(arg, 1);
  TEST_ASSERT(arg->num_breakpoints == 2, "should have 2 breakpoints after removal");
  TEST_ASSERT(arg->breakpoints[0].position == 200, "first bp should be at 200");
  TEST_ASSERT(arg->breakpoints[1].position == 800, "second bp should be at 800");

  merge_blocks(arg, 1);
  TEST_ASSERT(arg->num_blocks == 3, "should have 3 blocks after merge");
  TEST_ASSERT(arg->blocks[1].start == 200, "merged block start");
  TEST_ASSERT(arg->blocks[1].end == 800, "merged block end");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_smc_prior_zero_breakpoints(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  double log_prior;
  double rho = 0.001;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  arg->recomb_rate = rho;

  log_prior = arg_log_prior(arg, NULL, NULL);

  double expected = -rho * sites;
  TEST_ASSERT(fabs(log_prior - expected) < 1e-10,
              "SMC prior with k=0 should be -rho*L");

  arg->recomb_rate = 0.01;
  log_prior = arg_log_prior(arg, NULL, NULL);
  expected = -0.01 * sites;
  TEST_ASSERT(fabs(log_prior - expected) < 1e-10,
              "SMC prior should scale with rho");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_smc_prior_with_breakpoints(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  double log_prior;
  double rho = 0.001;
  breakpoint_t bp;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  arg->recomb_rate = rho;

  bp.position = 200;
  bp.lineage = 0;
  bp.recomb_time = 0.001;
  bp.target_pop = 0;
  bp.coal_time = 0.002;
  bp.target_node = 1;
  insert_breakpoint(arg, &bp);

  bp.position = 500;
  insert_breakpoint(arg, &bp);

  log_prior = arg_log_prior(arg, NULL, NULL);

  double lambda = rho * sites;
  double expected_poisson = 2 * log(lambda) - lambda - lgamma(3);
  double expected_uniform = -2 * log(sites);
  double expected = expected_poisson + expected_uniform;

  TEST_ASSERT(fabs(log_prior - expected) < 1e-10,
              "SMC prior with k=2 not matching expected");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_block_structure(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 500;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");
  TEST_ASSERT(arg->num_blocks == 1, "should have 1 block with k=0");
  TEST_ASSERT(arg->blocks[0].start == 0, "single block should cover all sites");
  TEST_ASSERT(arg->blocks[0].end == sites, "single block should cover all sites");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_max_breakpoints_reallocation(void)
{
  arg_t * arg;
  unsigned int max_bp = 5;
  unsigned int sites = 1000;
  breakpoint_t bp;
  unsigned int i;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");
  TEST_ASSERT(arg->max_breakpoints == 5, "max_breakpoints should be 5");

  for (i = 0; i < 5; i++)
  {
    bp.position = (i + 1) * 100;
    bp.lineage = 0;
    bp.recomb_time = 0.001;
    bp.target_pop = 0;
    bp.coal_time = 0.002;
    bp.target_node = 1;
    insert_breakpoint(arg, &bp);
  }

  TEST_ASSERT(arg->num_breakpoints == 5, "should have 5 breakpoints");

  /* Add one more - should trigger reallocation */
  bp.position = 600;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 6, "should have 6 breakpoints after realloc");
  TEST_ASSERT(arg->max_breakpoints == 10, "max_breakpoints should have doubled");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_small_sequence(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 10;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");
  TEST_ASSERT(arg->blocks[0].end == 10, "block end should be 10");

  arg_destroy(arg);

  arg = arg_create(max_bp, NULL, 1);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL for single site");
  TEST_ASSERT(arg->blocks[0].start == 0, "single site: block start");
  TEST_ASSERT(arg->blocks[0].end == 1, "single site: block end");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_insertion_removal_sequence(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  breakpoint_t bp;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  bp.lineage = 0;
  bp.recomb_time = 0.001;
  bp.target_pop = 0;
  bp.coal_time = 0.002;
  bp.target_node = 1;

  bp.position = 500;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 1, "insert 1");

  bp.position = 250;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 2, "insert 2");

  bp.position = 750;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 3, "insert 3");

  remove_breakpoint(arg, 1);
  TEST_ASSERT(arg->num_breakpoints == 2, "after remove");
  TEST_ASSERT(arg->breakpoints[0].position == 250, "remaining bp 0");
  TEST_ASSERT(arg->breakpoints[1].position == 750, "remaining bp 1");

  bp.position = 400;
  insert_breakpoint(arg, &bp);
  TEST_ASSERT(arg->num_breakpoints == 3, "insert again");
  TEST_ASSERT(arg->breakpoints[0].position == 250, "sorted bp 0");
  TEST_ASSERT(arg->breakpoints[1].position == 400, "sorted bp 1");
  TEST_ASSERT(arg->breakpoints[2].position == 750, "sorted bp 2");

  remove_breakpoint(arg, 0);
  TEST_ASSERT(arg->num_breakpoints == 2, "after remove first");
  TEST_ASSERT(arg->breakpoints[0].position == 400, "after remove: bp 0");
  TEST_ASSERT(arg->breakpoints[1].position == 750, "after remove: bp 1");

  remove_breakpoint(arg, 1);
  TEST_ASSERT(arg->num_breakpoints == 1, "after remove last");
  TEST_ASSERT(arg->breakpoints[0].position == 400, "single bp remaining");

  remove_breakpoint(arg, 0);
  TEST_ASSERT(arg->num_breakpoints == 0, "all removed");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_lgamma_values(void)
{
  double lg1, lg2, lg3, lg4;

  lg1 = lgamma(1);
  lg2 = lgamma(2);
  lg3 = lgamma(3);
  lg4 = lgamma(4);

  TEST_ASSERT(fabs(lg1) < 1e-10, "lgamma(1) should be 0");
  TEST_ASSERT(fabs(lg2) < 1e-10, "lgamma(2) should be 0");
  TEST_ASSERT(fabs(lg3 - log(2)) < 1e-10, "lgamma(3) should be log(2)");
  TEST_ASSERT(fabs(lg4 - log(6)) < 1e-10, "lgamma(4) should be log(6)");

  TEST_PASS();
}

static int test_poisson_probability(void)
{
  double lambda, logp;
  int k;

  lambda = 1.0;
  k = 0;
  logp = k * log(lambda) - lambda - lgamma(k + 1);
  TEST_ASSERT(fabs(logp - (-1.0)) < 1e-10, "Poisson(0|1)");

  k = 1;
  logp = k * log(lambda) - lambda - lgamma(k + 1);
  TEST_ASSERT(fabs(logp - (-1.0)) < 1e-10, "Poisson(1|1)");

  lambda = 2.0;
  k = 2;
  logp = k * log(lambda) - lambda - lgamma(k + 1);
  double expected = log(2.0) - 2.0;
  TEST_ASSERT(fabs(logp - expected) < 1e-10, "Poisson(2|2)");

  TEST_PASS();
}

static int test_block_merge_edge_cases(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  breakpoint_t bp;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  /* Add two breakpoints */
  bp.position = 300;
  bp.lineage = 0;
  bp.recomb_time = 0.001;
  bp.target_pop = 0;
  bp.coal_time = 0.002;
  bp.target_node = 1;
  insert_breakpoint(arg, &bp);

  bp.position = 700;
  insert_breakpoint(arg, &bp);

  update_blocks(arg, NULL);
  TEST_ASSERT(arg->num_blocks == 3, "should have 3 blocks");

  /* Try merging at invalid index (beyond num_blocks) */
  merge_blocks(arg, 5);  /* Should do nothing */
  TEST_ASSERT(arg->num_blocks == 3, "num_blocks unchanged after invalid merge");

  /* Merge last valid block */
  merge_blocks(arg, 1);
  TEST_ASSERT(arg->num_blocks == 2, "should have 2 blocks after merge");
  TEST_ASSERT(arg->blocks[0].start == 0, "block 0 start");
  TEST_ASSERT(arg->blocks[0].end == 300, "block 0 end");
  TEST_ASSERT(arg->blocks[1].start == 300, "block 1 start");
  TEST_ASSERT(arg->blocks[1].end == 1000, "block 1 end");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_duplicate_position_handling(void)
{
  arg_t * arg;
  unsigned int max_bp = 50;
  unsigned int sites = 1000;
  breakpoint_t bp;

  arg = arg_create(max_bp, NULL, sites);
  TEST_ASSERT(arg != NULL, "arg_create returned NULL");

  /* Add breakpoints at same position (implementation should allow this) */
  bp.position = 500;
  bp.lineage = 0;
  bp.recomb_time = 0.001;
  bp.target_pop = 0;
  bp.coal_time = 0.002;
  bp.target_node = 1;
  insert_breakpoint(arg, &bp);

  bp.lineage = 1;  /* Different lineage, same position */
  bp.recomb_time = 0.0015;
  insert_breakpoint(arg, &bp);

  TEST_ASSERT(arg->num_breakpoints == 2, "should have 2 breakpoints at same position");
  TEST_ASSERT(arg->breakpoints[0].position == 500, "bp 0 at 500");
  TEST_ASSERT(arg->breakpoints[1].position == 500, "bp 1 at 500");

  arg_destroy(arg);

  TEST_PASS();
}

static int test_prior_with_different_rho_values(void)
{
  arg_t * arg;
  unsigned int sites = 1000;
  double log_prior;

  arg = arg_create(50, NULL, sites);

  /* Test with very small rho */
  arg->recomb_rate = 1e-6;
  log_prior = arg_log_prior(arg, NULL, NULL);
  TEST_ASSERT(fabs(log_prior - (-0.001)) < 1e-8, "prior with small rho");

  /* Test with larger rho */
  arg->recomb_rate = 0.1;
  log_prior = arg_log_prior(arg, NULL, NULL);
  TEST_ASSERT(fabs(log_prior - (-100.0)) < 1e-8, "prior with large rho");

  arg_destroy(arg);

  TEST_PASS();
}

/* ============== Run all tests ============== */

static void run_test(int (*test_func)(void), const char * name)
{
  (void)name;  /* Unused, test function prints its own name */
  tests_run++;
  if (test_func() == 0)
    tests_passed++;
}

int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;

  fprintf(stderr, "\n========================================\n");
  fprintf(stderr, "BPP Recombination Unit Tests (Standalone)\n");
  fprintf(stderr, "========================================\n\n");

  run_test(test_arg_create_destroy, "ARG create/destroy");
  run_test(test_arg_reset, "ARG reset");
  run_test(test_breakpoint_insertion, "Breakpoint insertion");
  run_test(test_breakpoint_removal, "Breakpoint removal");
  run_test(test_smc_prior_zero_breakpoints, "SMC prior (k=0)");
  run_test(test_smc_prior_with_breakpoints, "SMC prior (k>0)");
  run_test(test_block_structure, "Block structure");
  run_test(test_max_breakpoints_reallocation, "Max breakpoints reallocation");
  run_test(test_small_sequence, "Small sequence");
  run_test(test_insertion_removal_sequence, "Insert/remove sequence");
  run_test(test_lgamma_values, "lgamma values");
  run_test(test_poisson_probability, "Poisson probability");
  run_test(test_block_merge_edge_cases, "Block merge edge cases");
  run_test(test_duplicate_position_handling, "Duplicate position handling");
  run_test(test_prior_with_different_rho_values, "Prior with different rho");

  fprintf(stderr, "\n========================================\n");
  fprintf(stderr, "Results: %d/%d tests passed\n", tests_passed, tests_run);
  fprintf(stderr, "========================================\n\n");

  return (tests_passed == tests_run) ? 0 : 1;
}
