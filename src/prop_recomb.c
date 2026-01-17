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
 * MCMC proposals for within-locus recombination using the SMC approximation.
 *
 * Proposals:
 * 1. Add breakpoint (reversible-jump)
 * 2. Remove breakpoint (reversible-jump)
 * 3. Move breakpoint position
 * 4. Modify recombination event (resample times and target)
 * 5. Update recombination rate (rho)
 */

/* Finetune parameters (can be made configurable later) */
static double finetune_breakpoint_window = 10.0;  /* window for moving breakpoints */
static double finetune_rho = 0.1;                 /* window for rho proposals */

/* Forward declarations */
static int position_is_valid(arg_t * arg, unsigned int pos, unsigned int sites);
static unsigned int sample_active_lineage(gtree_t * gtree, double time,
                                          long thread_index);
static double sample_recomb_time(gnode_t * node, long thread_index);
static void sample_recoalescence(gtree_t * gtree, stree_t * stree,
                                 double recomb_time, unsigned int lineage,
                                 unsigned int * target_pop,
                                 double * coal_time,
                                 unsigned int * target_node,
                                 long thread_index);

/**
 * Proposal 1: Add a new breakpoint (reversible-jump).
 *
 * Algorithm:
 * 1. Choose random position (not existing breakpoint)
 * 2. Choose random active lineage at that position
 * 3. Sample recombination time along lineage
 * 4. Sample re-coalescence per SMC (time and target)
 * 5. Insert breakpoint, rebuild blocks
 * 6. Hastings ratio: log(k) - log(L - k)
 *
 * @param locus         Locus data
 * @param gtree         Gene tree
 * @param stree         Species tree
 * @param thread_index  Thread index for random number generation
 * @return              1 if accepted, 0 if rejected
 */
long propose_add_breakpoint(locus_t * locus, gtree_t * gtree, stree_t * stree,
                            long thread_index)
{
  arg_t * arg;
  breakpoint_t new_bp;
  unsigned int sites;
  unsigned int k;
  double lnacceptance;
  double old_logl, new_logl;
  double old_prior, new_prior;

  if (!locus->has_recombination || !locus->arg)
    return 0;

  arg = locus->arg;
  sites = locus->original_sites;
  k = arg->num_breakpoints;

  /* Check if we can add more breakpoints */
  if (k >= arg->max_breakpoints)
    return 0;

  /* Check if there are valid positions left */
  if (k >= sites - 1)
    return 0;

  /* Sample a new breakpoint position (1 to sites-1, not existing) */
  unsigned int attempts = 0;
  unsigned int max_attempts = 100;
  do {
    new_bp.position = 1 + (unsigned int)(legacy_rndu(thread_index) * (sites - 1));
    attempts++;
  } while (!position_is_valid(arg, new_bp.position, sites) &&
           attempts < max_attempts);

  if (attempts >= max_attempts)
    return 0;

  /* Sample which lineage recombines */
  new_bp.lineage = sample_active_lineage(gtree, 0.0, thread_index);
  if (new_bp.lineage == (unsigned int)-1)
    return 0;

  gnode_t * recomb_node = gtree->nodes[new_bp.lineage];

  /* Sample recombination time along this lineage */
  new_bp.recomb_time = sample_recomb_time(recomb_node, thread_index);

  /* Sample re-coalescence under SMC */
  sample_recoalescence(gtree, stree, new_bp.recomb_time, new_bp.lineage,
                       &new_bp.target_pop, &new_bp.coal_time, &new_bp.target_node,
                       thread_index);

  /* Store old values */
  old_logl = gtree->logl;
  old_prior = arg->log_recomb_prior;

  /* Insert the breakpoint */
  insert_breakpoint(arg, &new_bp);
  update_blocks(arg, locus);

  /* Compute new likelihood and prior */
  new_logl = locus_root_loglikelihood_blocks(locus, gtree);
  new_prior = arg_log_prior(arg, stree, locus);

  /* Hastings ratio for reversible-jump:
   * Forward: propose add, choose position, choose lineage, sample times
   * Backward: propose remove, choose which breakpoint to remove
   * HR = P(remove) * (1/k+1) / [P(add) * (1/(L-k)) * q(times)]
   *
   * If P(add) = P(remove) = 0.5 and q(times) = prior, simplifies to:
   * log(HR) = log(L-k) - log(k+1)
   */
  lnacceptance = new_logl - old_logl;
  lnacceptance += new_prior - old_prior;
  lnacceptance += log((double)(sites - k)) - log((double)(k + 1));

  /* Accept or reject */
  if (lnacceptance >= 0 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* Accept */
    gtree->logl = new_logl;
    return 1;
  }
  else
  {
    /* Reject - restore old state */
    remove_breakpoint(arg, arg->num_breakpoints - 1);
    update_blocks(arg, locus);
    arg->log_recomb_prior = old_prior;
    return 0;
  }
}

/**
 * Proposal 2: Remove an existing breakpoint (reversible-jump).
 *
 * Algorithm:
 * 1. Choose random breakpoint to remove
 * 2. Remove and merge adjacent blocks
 * 3. Hastings ratio: log(L - k) - log(k + 1)
 *
 * @param locus         Locus data
 * @param gtree         Gene tree
 * @param stree         Species tree
 * @param thread_index  Thread index for random number generation
 * @return              1 if accepted, 0 if rejected
 */
long propose_remove_breakpoint(locus_t * locus, gtree_t * gtree, stree_t * stree,
                               long thread_index)
{
  arg_t * arg;
  breakpoint_t removed_bp;
  unsigned int sites;
  unsigned int k;
  unsigned int remove_idx;
  double lnacceptance;
  double old_logl, new_logl;
  double old_prior, new_prior;

  if (!locus->has_recombination || !locus->arg)
    return 0;

  arg = locus->arg;
  sites = locus->original_sites;
  k = arg->num_breakpoints;

  /* Can't remove if there are no breakpoints */
  if (k == 0)
    return 0;

  /* Choose random breakpoint to remove */
  remove_idx = (unsigned int)(legacy_rndu(thread_index) * k);
  removed_bp = arg->breakpoints[remove_idx];

  /* Store old values */
  old_logl = gtree->logl;
  old_prior = arg->log_recomb_prior;

  /* Remove the breakpoint */
  remove_breakpoint(arg, remove_idx);
  merge_blocks(arg, remove_idx);
  update_blocks(arg, locus);

  /* Compute new likelihood and prior */
  new_logl = locus_root_loglikelihood_blocks(locus, gtree);
  new_prior = arg_log_prior(arg, stree, locus);

  /* Hastings ratio (reverse of add):
   * log(HR) = log(k) - log(L-k+1)
   */
  lnacceptance = new_logl - old_logl;
  lnacceptance += new_prior - old_prior;
  lnacceptance += log((double)k) - log((double)(sites - k + 1));

  /* Accept or reject */
  if (lnacceptance >= 0 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* Accept */
    gtree->logl = new_logl;
    return 1;
  }
  else
  {
    /* Reject - restore old state */
    insert_breakpoint(arg, &removed_bp);
    update_blocks(arg, locus);
    arg->log_recomb_prior = old_prior;
    return 0;
  }
}

/**
 * Proposal 3: Move a breakpoint to a new position.
 *
 * Algorithm:
 * 1. Choose random breakpoint
 * 2. Sliding window to new position
 * 3. Rebuild blocks
 * 4. Symmetric proposal (HR = 0)
 *
 * @param locus         Locus data
 * @param gtree         Gene tree
 * @param stree         Species tree
 * @param thread_index  Thread index for random number generation
 * @return              1 if accepted, 0 if rejected
 */
long propose_move_breakpoint(locus_t * locus, gtree_t * gtree, stree_t * stree,
                             long thread_index)
{
  arg_t * arg;
  unsigned int sites;
  unsigned int k;
  unsigned int move_idx;
  unsigned int old_pos, new_pos;
  int delta;
  double lnacceptance;
  double old_logl, new_logl;
  double old_prior, new_prior;

  if (!locus->has_recombination || !locus->arg)
    return 0;

  arg = locus->arg;
  sites = locus->original_sites;
  k = arg->num_breakpoints;

  /* Need at least one breakpoint */
  if (k == 0)
    return 0;

  /* Choose random breakpoint to move */
  move_idx = (unsigned int)(legacy_rndu(thread_index) * k);
  old_pos = arg->breakpoints[move_idx].position;

  /* Propose new position using sliding window */
  delta = (int)((legacy_rndu(thread_index) - 0.5) * 2.0 * finetune_breakpoint_window);
  new_pos = (unsigned int)((int)old_pos + delta);

  /* Reflect at boundaries */
  if ((int)new_pos < 1)
    new_pos = 2 - new_pos;
  if (new_pos >= sites)
    new_pos = 2 * (sites - 1) - new_pos;

  /* Check if new position is valid */
  arg->breakpoints[move_idx].position = new_pos;
  if (!position_is_valid(arg, new_pos, sites))
  {
    arg->breakpoints[move_idx].position = old_pos;
    return 0;
  }

  /* Store old values */
  old_logl = gtree->logl;
  old_prior = arg->log_recomb_prior;

  /* Update blocks with new position */
  update_blocks(arg, locus);

  /* Compute new likelihood and prior */
  new_logl = locus_root_loglikelihood_blocks(locus, gtree);
  new_prior = arg_log_prior(arg, stree, locus);

  /* Symmetric proposal: HR = 0 */
  lnacceptance = new_logl - old_logl;
  lnacceptance += new_prior - old_prior;

  /* Accept or reject */
  if (lnacceptance >= 0 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* Accept */
    gtree->logl = new_logl;
    return 1;
  }
  else
  {
    /* Reject - restore old position */
    arg->breakpoints[move_idx].position = old_pos;
    update_blocks(arg, locus);
    arg->log_recomb_prior = old_prior;
    return 0;
  }
}

/**
 * Proposal 4: Modify a recombination event (resample times and target).
 *
 * Algorithm:
 * 1. Choose random breakpoint
 * 2. Resample recomb_time, coal_time, target_node
 * 3. Update local tree
 *
 * @param locus         Locus data
 * @param gtree         Gene tree
 * @param stree         Species tree
 * @param thread_index  Thread index for random number generation
 * @return              1 if accepted, 0 if rejected
 */
long propose_modify_recomb_event(locus_t * locus, gtree_t * gtree, stree_t * stree,
                                 long thread_index)
{
  arg_t * arg;
  unsigned int k;
  unsigned int modify_idx;
  breakpoint_t old_bp, new_bp;
  double lnacceptance;
  double old_logl, new_logl;
  double old_prior, new_prior;

  if (!locus->has_recombination || !locus->arg)
    return 0;

  arg = locus->arg;
  k = arg->num_breakpoints;

  /* Need at least one breakpoint */
  if (k == 0)
    return 0;

  /* Choose random breakpoint to modify */
  modify_idx = (unsigned int)(legacy_rndu(thread_index) * k);
  old_bp = arg->breakpoints[modify_idx];
  new_bp = old_bp;

  /* Resample the recombination event */
  gnode_t * recomb_node = gtree->nodes[new_bp.lineage];
  new_bp.recomb_time = sample_recomb_time(recomb_node, thread_index);

  sample_recoalescence(gtree, stree, new_bp.recomb_time, new_bp.lineage,
                       &new_bp.target_pop, &new_bp.coal_time, &new_bp.target_node,
                       thread_index);

  /* Store old values */
  old_logl = gtree->logl;
  old_prior = arg->log_recomb_prior;

  /* Apply new breakpoint */
  arg->breakpoints[modify_idx] = new_bp;
  update_blocks(arg, locus);

  /* Compute new likelihood and prior */
  new_logl = locus_root_loglikelihood_blocks(locus, gtree);
  new_prior = arg_log_prior(arg, stree, locus);

  /* Prior ratio serves as proposal ratio for Gibbs-like sampling */
  lnacceptance = new_logl - old_logl;
  lnacceptance += new_prior - old_prior;

  /* Accept or reject */
  if (lnacceptance >= 0 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* Accept */
    gtree->logl = new_logl;
    return 1;
  }
  else
  {
    /* Reject - restore old breakpoint */
    arg->breakpoints[modify_idx] = old_bp;
    update_blocks(arg, locus);
    arg->log_recomb_prior = old_prior;
    return 0;
  }
}

/**
 * Proposal 5: Update recombination rate (rho) across all loci.
 *
 * Algorithm:
 * 1. Log-normal sliding window: rho' = rho * exp(c*(U - 0.5))
 * 2. Prior ratio from gamma prior
 * 3. Likelihood ratio from Poisson on k
 *
 * @param locus         Array of locus data
 * @param locus_count   Number of loci
 * @param stree         Species tree
 * @param thread_index  Thread index for random number generation
 * @return              1 if accepted, 0 if rejected
 */
long propose_rho(locus_t ** locus, long locus_count, stree_t * stree,
                 long thread_index)
{
  long i;
  double old_rho, new_rho;
  double lnacceptance = 0.0;
  double old_prior_sum = 0.0;
  double new_prior_sum = 0.0;
  arg_t * arg;

  /* Get current rho from first locus with recombination */
  old_rho = 0.001;  /* Default */
  for (i = 0; i < locus_count; i++)
  {
    if (locus[i]->has_recombination && locus[i]->arg)
    {
      old_rho = locus[i]->arg->recomb_rate;
      break;
    }
  }

  /* Propose new rho using log-normal sliding window */
  double c = finetune_rho;
  new_rho = old_rho * exp(c * (legacy_rndu(thread_index) - 0.5));

  /* Reject if rho is too small or too large */
  if (new_rho < 1e-10 || new_rho > 1.0)
    return 0;

  /* Compute prior and Poisson ratios */
  for (i = 0; i < locus_count; i++)
  {
    if (!locus[i]->has_recombination || !locus[i]->arg)
      continue;

    arg = locus[i]->arg;
    unsigned int sites = locus[i]->sites;
    unsigned int k = arg->num_breakpoints;

    /* Old Poisson log-likelihood */
    double old_lambda = old_rho * (double)sites;
    old_prior_sum += (double)k * log(old_lambda) - old_lambda;

    /* New Poisson log-likelihood */
    double new_lambda = new_rho * (double)sites;
    new_prior_sum += (double)k * log(new_lambda) - new_lambda;
  }

  /* Gamma prior on rho: Gamma(alpha, beta) */
  double alpha = opt_rho_alpha;
  double beta = opt_rho_beta;

  old_prior_sum += (alpha - 1.0) * log(old_rho) - beta * old_rho;
  new_prior_sum += (alpha - 1.0) * log(new_rho) - beta * new_rho;

  /* Jacobian for log-normal proposal */
  lnacceptance = new_prior_sum - old_prior_sum;
  lnacceptance += log(new_rho) - log(old_rho);  /* Jacobian */

  /* Accept or reject */
  if (lnacceptance >= 0 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    /* Accept - update rho in all ARGs */
    for (i = 0; i < locus_count; i++)
    {
      if (locus[i]->has_recombination && locus[i]->arg)
      {
        locus[i]->arg->recomb_rate = new_rho;
      }
    }
    return 1;
  }
  else
  {
    return 0;
  }
}

/* ============== Static helper functions ============== */

/**
 * Check if a position is valid for a new breakpoint.
 * Position must be in [1, sites-1] and not already used.
 */
static int position_is_valid(arg_t * arg, unsigned int pos, unsigned int sites)
{
  unsigned int i;

  if (pos < 1 || pos >= sites)
    return 0;

  for (i = 0; i < arg->num_breakpoints; i++)
  {
    if (arg->breakpoints[i].position == pos)
      return 0;
  }

  return 1;
}

/**
 * Sample an active lineage at a given time.
 * Returns the node index, or -1 if no lineage found.
 */
static unsigned int sample_active_lineage(gtree_t * gtree, double time,
                                          long thread_index)
{
  unsigned int count = 0;
  unsigned int * active;
  unsigned int i;
  unsigned int selected;
  unsigned int total_nodes = gtree->tip_count + gtree->inner_count;

  active = (unsigned int *)xmalloc(total_nodes * sizeof(unsigned int));

  /* Find lineages active at the given time */
  for (i = 0; i < total_nodes; i++)
  {
    gnode_t * node = gtree->nodes[i];

    /* A lineage is active if:
     * - node's time <= time (lineage has started)
     * - parent's time > time or no parent (lineage hasn't ended)
     */
    if (node->time <= time)
    {
      if (!node->parent || node->parent->time > time)
      {
        active[count++] = i;
      }
    }
  }

  if (count == 0)
  {
    free(active);
    return (unsigned int)-1;
  }

  /* Select random active lineage */
  selected = active[(unsigned int)(legacy_rndu(thread_index) * count)];
  free(active);

  return selected;
}

/**
 * Sample a recombination time along a lineage.
 * Time is uniform between node's time and parent's time.
 */
static double sample_recomb_time(gnode_t * node, long thread_index)
{
  double min_time, max_time;

  min_time = node->time;
  max_time = node->parent ? node->parent->time : min_time + 1.0;

  return min_time + legacy_rndu(thread_index) * (max_time - min_time);
}

/**
 * Sample re-coalescence under SMC.
 * After detachment at recomb_time, the lineage coalesces with the
 * marginal tree at the previous position.
 */
static void sample_recoalescence(gtree_t * gtree, stree_t * stree,
                                 double recomb_time, unsigned int lineage,
                                 unsigned int * target_pop,
                                 double * coal_time,
                                 unsigned int * target_node,
                                 long thread_index)
{
  unsigned int i;
  unsigned int total_nodes = gtree->tip_count + gtree->inner_count;
  gnode_t * node = gtree->nodes[lineage];
  snode_t * pop;

  /* Find the population at recombination time */
  pop = node->pop;
  while (pop->parent && pop->parent->tau <= recomb_time)
  {
    pop = pop->parent;
  }

  *target_pop = pop->node_index;

  /* Sample coalescence time (exponential waiting time) */
  double theta = pop->theta > 0 ? pop->theta : 0.001;
  double rate = 1.0 / theta;
  *coal_time = recomb_time + legacy_rndexp(thread_index, 1.0 / rate);

  /* Sample target node (lineage to coalesce with) */
  /* Find lineages active at coal_time */
  unsigned int * active = (unsigned int *)xmalloc(total_nodes * sizeof(unsigned int));
  unsigned int count = 0;

  for (i = 0; i < total_nodes; i++)
  {
    if (i == lineage)
      continue;  /* Can't coalesce with self */

    gnode_t * other = gtree->nodes[i];
    if (other->time <= *coal_time)
    {
      if (!other->parent || other->parent->time > *coal_time)
      {
        active[count++] = i;
      }
    }
  }

  if (count > 0)
  {
    *target_node = active[(unsigned int)(legacy_rndu(thread_index) * count)];
  }
  else
  {
    /* No available target - use root */
    *target_node = gtree->root->node_index;
  }

  free(active);
}
