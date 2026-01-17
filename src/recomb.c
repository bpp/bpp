/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

/*
 * Recombination module implementing within-locus recombination using
 * the Sequentially Markov Coalescent (SMC) approximation.
 *
 * Under the SMC, a recombination event at position x creates a new lineage
 * that coalesces with the marginal tree at x-1. The local tree at position x
 * differs from x-1 only in where the recombinant lineage attaches.
 */

/* Forward declarations for static helper functions */
static void local_tree_destroy(gtree_t * tree);

/**
 * Create an ARG structure with initial allocation.
 *
 * @param max_bp   Maximum number of breakpoints to allocate space for
 * @param base     Backbone tree (leftmost block)
 * @param sites    Number of sites in the locus
 * @return         Newly allocated ARG structure
 */
arg_t * arg_create(unsigned int max_bp, gtree_t * base, unsigned int sites)
{
  arg_t * arg;

  arg = (arg_t *)xmalloc(sizeof(arg_t));

  arg->num_breakpoints = 0;
  arg->max_breakpoints = max_bp;
  arg->breakpoints = (breakpoint_t *)xmalloc(max_bp * sizeof(breakpoint_t));

  /* Initially one block covering entire sequence */
  arg->num_blocks = 1;
  arg->blocks = (block_t *)xmalloc((max_bp + 1) * sizeof(block_t));

  /* Initialize the single block */
  arg->blocks[0].start = 0;
  arg->blocks[0].end = sites;
  arg->blocks[0].local_tree = base;  /* Use the base tree for block 0 */
  arg->blocks[0].logl = 0.0;

  arg->base_tree = base;
  arg->recomb_rate = 0.001;  /* Default rho value */
  arg->log_recomb_prior = 0.0;

  return arg;
}

/**
 * Destroy an ARG structure and free all memory.
 *
 * @param arg  ARG structure to destroy
 */
void arg_destroy(arg_t * arg)
{
  unsigned int i;

  if (!arg) return;

  /* Free local trees for blocks other than the first (base tree) */
  for (i = 1; i < arg->num_blocks; i++)
  {
    if (arg->blocks[i].local_tree &&
        arg->blocks[i].local_tree != arg->base_tree)
    {
      local_tree_destroy(arg->blocks[i].local_tree);
    }
  }

  if (arg->breakpoints)
    free(arg->breakpoints);
  if (arg->blocks)
    free(arg->blocks);

  free(arg);
}

/**
 * Reset ARG to have zero breakpoints.
 *
 * @param arg    ARG structure to reset
 * @param sites  Number of sites (0 to use existing blocks[0].end)
 */
void arg_reset(arg_t * arg, unsigned int sites)
{
  unsigned int i;

  if (!arg) return;

  /* Free local trees for blocks other than the first */
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

  /* Reset first block to cover entire sequence */
  arg->blocks[0].start = 0;
  if (sites > 0)
    arg->blocks[0].end = sites;
  arg->blocks[0].local_tree = arg->base_tree;
  arg->blocks[0].logl = 0.0;
}

/**
 * Compute the log prior probability of the ARG under the SMC model.
 *
 * The prior includes:
 * - Poisson(rho*L) on the number of breakpoints
 * - Uniform positions: -k*log(L)
 * - MSC prior for each local tree (if stree and locus provided)
 * - SMC transition probabilities (if stree provided)
 *
 * @param arg    ARG structure
 * @param stree  Species tree (NULL for basic prior only)
 * @param locus  Locus data (NULL for basic prior only)
 * @return       Log prior probability
 */
double arg_log_prior(arg_t * arg, stree_t * stree, locus_t * locus)
{
  double logpr = 0.0;
  unsigned int L;
  double rho = arg->recomb_rate;
  unsigned int k = arg->num_breakpoints;
  unsigned int b, i;

  /* Get number of sites - from locus if available, otherwise from first block */
  if (locus)
    L = locus->sites;
  else
    L = arg->blocks[0].end;  /* Assume single block covers all sites */

  /* Poisson(rho*L) on number of breakpoints */
  double lambda = rho * (double)L;
  if (k > 0)
    logpr += (double)k * log(lambda) - lambda - lgamma((double)k + 1.0);
  else
    logpr -= lambda;  /* k=0: just -lambda term */

  /* Uniform positions: -k*log(L) */
  if (k > 0)
    logpr -= (double)k * log((double)L);

  /* MSC prior for each local tree (only if stree and locus provided) */
  if (stree && locus)
  {
    for (b = 0; b < arg->num_blocks; b++)
    {
      if (arg->blocks[b].local_tree)
      {
        logpr += gtree_logprob(stree, locus->heredity[0],
                               arg->blocks[b].local_tree->msa_index, 0);
      }
    }
  }

  /* SMC transition probabilities (only if stree provided) */
  if (stree)
  {
    for (i = 0; i < k; i++)
    {
      logpr += smc_transition_logprob(arg, &arg->breakpoints[i], stree);
    }
  }

  arg->log_recomb_prior = logpr;
  return logpr;
}

/**
 * Compute the SMC transition probability for a single breakpoint.
 *
 * Under the SMC, when a lineage is cut at time t_recomb, it must reattach
 * to a lineage that exists in the marginal tree at the previous position.
 *
 * @param arg    ARG structure
 * @param bp     Breakpoint
 * @param stree  Species tree
 * @return       Log probability of the SMC transition
 */
double smc_transition_logprob(arg_t * arg, breakpoint_t * bp, stree_t * stree)
{
  double logpr = 0.0;
  double dt;
  snode_t * pop;

  /*
   * Under SMC, after detachment at time t_recomb, the lineage enters
   * a coalescent process. The probability depends on:
   * 1. The number of available lineages at time t_recomb
   * 2. The waiting time until coalescence
   * 3. The population parameters (theta)
   */

  /* Time from recombination to re-coalescence */
  dt = bp->coal_time - bp->recomb_time;
  if (dt < 0)
    return -INFINITY;  /* Invalid: coalescence before recombination */

  /* Get the population where re-coalescence occurs */
  if (bp->target_pop < stree->tip_count + stree->inner_count)
  {
    pop = stree->nodes[bp->target_pop];

    /* Coalescent rate in this population */
    if (pop->theta > 0)
    {
      /* Simple exponential waiting time model */
      /* Rate = 1/theta for a single lineage pair */
      double rate = 1.0 / pop->theta;
      logpr = log(rate) - rate * dt;
    }
  }

  return logpr;
}

/**
 * Construct a local tree by applying a breakpoint's recombination event
 * to a base tree.
 *
 * Under the SMC model, a recombination at site position bp->position causes
 * the lineage bp->lineage to detach at time bp->recomb_time and reattach
 * to bp->target_node at time bp->coal_time.
 *
 * This follows BPP's SPR pattern: the pruned parent node is reused as
 * the coalescence point for the regraft.
 *
 * @param local  Output local tree (pre-allocated with nodes)
 * @param base   Base tree to copy from
 * @param bp     Breakpoint defining the recombination event
 * @param stree  Species tree (for population references)
 */
void construct_local_tree(gtree_t * local, gtree_t * base, breakpoint_t * bp,
                          stree_t * stree)
{
  gnode_t * node;
  gnode_t * target;
  gnode_t * parent;

  if (!local || !base || !bp)
    return;

  /* Step 1: Clone the base tree to the local tree */
  gtree_clone_for_recomb(local, base, stree);

  /* Step 2: Get the lineage that recombines (using index from breakpoint) */
  if (bp->lineage >= local->tip_count + local->inner_count)
    return;  /* Invalid lineage index */
  node = local->nodes[bp->lineage];

  /* Step 3: Get the target node where lineage reattaches */
  if (bp->target_node >= local->tip_count + local->inner_count)
    return;  /* Invalid target index */
  target = local->nodes[bp->target_node];

  /* Step 4: Prune the recombinant lineage
   * This detaches node's parent from the tree and returns it for reuse */
  parent = gtree_prune_for_smc(local, node, bp->recomb_time);
  if (!parent)
    return;  /* Prune failed (e.g., node is root) */

  /* Step 5: Regraft onto target branch at coalescence time
   * Uses parent as the new internal coalescence node */
  gtree_regraft_for_smc(local, node, parent, target, bp->coal_time);

  /* Step 6: Reset CLV indices to mark all partials as invalid */
  gtree_reset_clv_indices(local);
}

/**
 * Insert a breakpoint into the ARG, maintaining sorted order by position.
 *
 * @param arg  ARG structure
 * @param bp   Breakpoint to insert
 */
void insert_breakpoint(arg_t * arg, breakpoint_t * bp)
{
  unsigned int i, j;

  if (arg->num_breakpoints >= arg->max_breakpoints)
  {
    /* Reallocate if needed */
    arg->max_breakpoints *= 2;
    arg->breakpoints = (breakpoint_t *)xrealloc(arg->breakpoints,
                        arg->max_breakpoints * sizeof(breakpoint_t));
    arg->blocks = (block_t *)xrealloc(arg->blocks,
                   (arg->max_breakpoints + 1) * sizeof(block_t));
  }

  /* Find insertion position to maintain sorted order */
  for (i = 0; i < arg->num_breakpoints; i++)
  {
    if (arg->breakpoints[i].position > bp->position)
      break;
  }

  /* Shift existing breakpoints */
  for (j = arg->num_breakpoints; j > i; j--)
  {
    arg->breakpoints[j] = arg->breakpoints[j - 1];
  }

  /* Insert new breakpoint */
  arg->breakpoints[i] = *bp;
  arg->num_breakpoints++;
}

/**
 * Remove a breakpoint from the ARG.
 *
 * @param arg  ARG structure
 * @param idx  Index of breakpoint to remove
 */
void remove_breakpoint(arg_t * arg, unsigned int idx)
{
  unsigned int i;

  if (idx >= arg->num_breakpoints)
    return;

  /* Shift remaining breakpoints */
  for (i = idx; i < arg->num_breakpoints - 1; i++)
  {
    arg->breakpoints[i] = arg->breakpoints[i + 1];
  }

  arg->num_breakpoints--;
}

/**
 * Update block boundaries after breakpoint changes.
 *
 * @param arg    ARG structure
 * @param locus  Locus data (NULL to use existing block end for sites)
 */
void update_blocks(arg_t * arg, locus_t * locus)
{
  unsigned int i;
  unsigned int sites;

  /* Get sites from locus if available, otherwise use current last block end */
  if (locus)
    sites = locus->sites;
  else
    sites = arg->blocks[arg->num_blocks - 1].end;

  arg->num_blocks = arg->num_breakpoints + 1;

  /* Set block boundaries based on breakpoint positions */
  arg->blocks[0].start = 0;

  for (i = 0; i < arg->num_breakpoints; i++)
  {
    arg->blocks[i].end = arg->breakpoints[i].position;
    arg->blocks[i + 1].start = arg->breakpoints[i].position;
  }

  arg->blocks[arg->num_blocks - 1].end = sites;

  /* First block always uses base tree */
  arg->blocks[0].local_tree = arg->base_tree;

  /*
   * TODO: For proper SMC implementation, subsequent blocks should have
   * distinct local trees constructed via construct_local_tree().
   * This requires allocating gtree structures with CLVs and matrices.
   *
   * For now, use the base tree for all blocks. This means the likelihood
   * is computed correctly (sum of per-site likelihoods) but without the
   * topology changes that recombination should introduce. The prior still
   * accounts for the breakpoint positions and recombination events.
   */
  for (i = 1; i < arg->num_blocks; i++)
  {
    arg->blocks[i].local_tree = arg->base_tree;
  }
}

/**
 * Merge two adjacent blocks after removing a breakpoint.
 *
 * @param arg  ARG structure
 * @param idx  Index of the breakpoint being removed (between blocks idx and idx+1)
 */
void merge_blocks(arg_t * arg, unsigned int idx)
{
  unsigned int i;

  if (idx >= arg->num_blocks - 1)
    return;

  /* Free the local tree of the block being merged away */
  if (arg->blocks[idx + 1].local_tree &&
      arg->blocks[idx + 1].local_tree != arg->base_tree)
  {
    local_tree_destroy(arg->blocks[idx + 1].local_tree);
  }

  /* Extend the current block to cover the merged block */
  arg->blocks[idx].end = arg->blocks[idx + 1].end;

  /* Shift remaining blocks */
  for (i = idx + 1; i < arg->num_blocks - 1; i++)
  {
    arg->blocks[i] = arg->blocks[i + 1];
  }

  arg->num_blocks--;
}

/**
 * Compute the log-likelihood across all blocks.
 *
 * @param locus  Locus data
 * @param gtree  Gene tree (used for base tree reference)
 * @return       Total log-likelihood
 */
double locus_root_loglikelihood_blocks(locus_t * locus, gtree_t * gtree)
{
  /* If no recombination, use standard likelihood computation */
  if (!locus->has_recombination || !locus->arg)
  {
    return locus_root_loglikelihood(locus, gtree->root,
                                    locus->param_indices, NULL);
  }

  /*
   * TODO: Proper block-wise likelihood requires handling the mapping between
   * original sequence positions (breakpoint coordinates) and compressed site
   * patterns (persite_lnl indices). For now, use the full sequence likelihood.
   *
   * With this simplification:
   * - The likelihood is computed correctly using the base tree
   * - The prior (Poisson on k, SMC transitions) still affects proposal acceptance
   * - Actual topology changes from recombination are not reflected in likelihood
   *
   * To implement proper block-wise likelihood:
   * 1. Track mapping from original positions to pattern indices
   * 2. Construct distinct local trees for each block
   * 3. Compute partials only for patterns in each block's range
   */
  return locus_root_loglikelihood(locus, gtree->root,
                                  locus->param_indices, NULL);
}

/**
 * Compute the log-likelihood for a single block (site range).
 *
 * @param locus  Locus data
 * @param gtree  Gene tree for this block
 * @param start  Start site (inclusive)
 * @param end    End site (exclusive)
 * @return       Log-likelihood for this block
 */
double compute_block_likelihood(locus_t * locus, gtree_t * gtree,
                                unsigned int start, unsigned int end)
{
  double logl = 0.0;
  unsigned int i;

  if (!gtree || !locus || start >= end)
    return 0.0;

  /* If block covers entire sequence, just use standard likelihood */
  if (start == 0 && end >= locus->sites)
  {
    return locus_root_loglikelihood(locus, gtree->root,
                                    locus->param_indices, NULL);
  }

  /*
   * For partial blocks, sum the per-site log-likelihoods for sites in [start, end)
   * Note: persite_lnl[i] already includes pattern_weights[i], so just sum directly
   *
   * TODO: Optimize by computing CLVs only for affected sites
   */

  /* Update all matrices and partials for this tree */
  locus_update_all_matrices(locus, gtree, NULL, gtree->msa_index);
  locus_update_all_partials(locus, gtree);

  /* Get per-site log-likelihoods (already weighted by pattern_weights) */
  double * persite_lnl = (double *)xmalloc(locus->sites * sizeof(double));

  locus_root_loglikelihood(locus, gtree->root, locus->param_indices, persite_lnl);

  /* Sum log-likelihoods for sites in this block (no extra weighting needed) */
  for (i = start; i < end && i < locus->sites; i++)
  {
    logl += persite_lnl[i];
  }

  free(persite_lnl);
  return logl;
}

/* ============== Static helper functions ============== */

/**
 * Destroy a local tree.
 */
static void local_tree_destroy(gtree_t * tree)
{
  if (tree)
  {
    gtree_destroy(tree, NULL);
  }
}
