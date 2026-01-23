# rjMCMC Species Delimitation with Migration (MSC-M)

## Overview

This document describes the implementation of MSC-M (Multispecies Coalescent with Migration) support for species delimitation in BPP, along with analysis of acceptance rate issues and planned improvements.

**Date**: January 2026
**Branch**: delimitation
**Related files**: `src/prop_rj.c`, `src/miginfo.c`, `src/stree.c`

---

## Part 1: Implemented Changes

### 1.1 Problem Statement

Prior to these changes, the species delimitation proposals (`prop_split()` and `prop_join()` in `src/prop_rj.c`) ignored migration events entirely. When populations were split or joined:

- **SPLIT**: Migration events with the split population as source/target were not reassigned to child populations
- **JOIN**: Migration events with child populations as source/target were not reassigned to the parent; migrations between the two children (L↔R) could become invalid self-migrations

This caused silent data corruption and incorrect MCMC behavior when running species delimitation under the MSC-M model.

### 1.2 Implementation Strategy

**Approach**: Reject-on-conflict with migration reassignment

1. For JOIN: Reject proposals if any L↔R migrations exist (they would become self-migrations)
2. For SPLIT/JOIN: Reassign migration source/target similar to how coalescent events are reassigned
3. Maintain all invariants for dual storage consistency and count tracking
4. Follow patterns from `propose_tau_mig()` in `stree.c`

### 1.3 Code Changes

#### 1.3.1 New Data Structures (`src/prop_rj.c`, lines 47-58)

```c
/* Storage for migration rollback */
typedef struct mig_rollback_s
{
  migevent_t * me;
  snode_t * old_source;
  snode_t * old_target;
  long msa_index;
} mig_rollback_t;

static mig_rollback_t * mig_rb;
static long mig_rb_count;
static long mig_rb_alloc;
```

**Purpose**: Store migration event modifications for rollback on proposal rejection.

#### 1.3.2 New Helper Functions

| Function | Lines | Purpose |
|----------|-------|---------|
| `mig_rb_reset()` | 140-143 | Reset rollback storage at start of each proposal |
| `mig_rb_store()` | 146-162 | Store a migration event's old state for potential rollback |
| `check_mig_join_validity()` | 166-200 | Check if any L↔R migrations exist (reject join if so) |
| `update_migs_split()` | 205-283 | Reassign migrations from parent to children during split |
| `update_migs_join()` | 286-357 | Reassign migrations from children to parent during join |
| `revert_migs()` | 360-403 | Restore all migrations to original state on rejection |

#### 1.3.3 Function Details

**`check_mig_join_validity()`**
```c
static int check_mig_join_validity(stree_t * stree,
                                   gtree_t ** gtree,
                                   snode_t * snode)
```
- Iterates through all gene trees and all migration events
- Returns 1 if any migration has (source=L, target=R) or (source=R, target=L)
- Called early in `prop_join()` to reject invalid proposals

**`update_migs_split()`**
```c
static void update_migs_split(stree_t * stree,
                              gtree_t * gtree,
                              snode_t * snode,
                              double tau_new,
                              long msa_index)
```
- For each migration event with time < tau_new involving snode:
  - Unlink from current population lists via `migevent_unlink()`
  - Determine child (L or R) using gene tree node ancestry marks
  - Update `me->source` and/or `me->target` to child
  - Relink to new population lists via `migevent_link()`
  - Update `snode->migevent_count[]` and `gtree->migcount[][]`
  - Store in rollback buffer

**`update_migs_join()`**
```c
static void update_migs_join(stree_t * stree,
                             gtree_t * gtree,
                             snode_t * snode,
                             long msa_index)
```
- For each migration event involving L or R as source/target:
  - Unlink, update to parent snode, relink
  - Update counts and store for rollback

**`revert_migs()`**
```c
static void revert_migs(gtree_t ** gtree)
```
- Iterates through rollback buffer in reverse
- Restores original source/target for each modified migration
- Restores all count arrays

#### 1.3.4 Modifications to `prop_split()` (lines 667-1226)

1. **Line 689-690**: Reset migration rollback at function start
   ```c
   if (opt_migration)
     mig_rb_reset();
   ```

2. **Line 734-735**: Update migration bands after tau change
   ```c
   if (opt_migration)
     stree_update_mig_subpops(stree, thread_index);
   ```

3. **Line 849-850**: Reassign migrations after rubber_proportional
   ```c
   if (opt_migration)
     update_migs_split(stree, gtree[i], node, tau_new, i);
   ```

4. **Line 1077-1078**: Revert migrations on rejection (before tau revert)
   ```c
   if (opt_migration)
     revert_migs(gtree);
   ```

5. **Line 1085-1086**: Restore migration bands after tau revert
   ```c
   if (opt_migration)
     stree_update_mig_subpops(stree, thread_index);
   ```

#### 1.3.5 Modifications to `prop_join()` (lines 1228-1640)

1. **Line 1248-1249**: Reset migration rollback at function start
   ```c
   if (opt_migration)
     mig_rb_reset();
   ```

2. **Line 1286-1287**: Early rejection for L↔R migrations
   ```c
   if (opt_migration && check_mig_join_validity(stree, gtree, node))
     return 2;  /* reject - L<->R migrations exist */
   ```

3. **Line 1367-1368**: Update migration bands after tau=0
   ```c
   if (opt_migration)
     stree_update_mig_subpops(stree, thread_index);
   ```

4. **Line 1413-1414**: Reassign migrations after rubber_proportional
   ```c
   if (opt_migration)
     update_migs_join(stree, gtree[i], node, i);
   ```

5. **Line 1571-1572**: Revert migrations on rejection
   ```c
   if (opt_migration)
     revert_migs(gtree);
   ```

6. **Line 1579-1580**: Restore migration bands after tau revert
   ```c
   if (opt_migration)
     stree_update_mig_subpops(stree, thread_index);
   ```

#### 1.3.6 Initialization and Cleanup

**`rj_init()`** (lines 108-111):
```c
/* initialize migration rollback storage */
mig_rb = NULL;
mig_rb_count = 0;
mig_rb_alloc = 0;
```

**`rj_fini()`** (lines 121-123):
```c
/* free migration rollback storage */
if (mig_rb)
  free(mig_rb);
```

### 1.4 Key Invariants Maintained

1. **Dual storage consistency**: `gnode->mi` events match `snode->mig_source[]` and `snode->mig_target[]` lists
2. **Count consistency**: `snode->migevent_count[]` matches actual count of events involving that population
3. **Direction counts**: `gtree->migcount[src][tgt]` matches actual migration directions
4. **Time ordering**: Migration events remain sorted by time within each gene tree node

### 1.5 Testing

Test case created in `test/mig_delimit_test/`:
- Simulated data with 2 populations, M=10 bidirectional migration, 50 loci, 10 individuals per population
- Control files: `sim_veryhighmig.ctl` (simulation), `infer_delimit_mig.ctl` (inference)
- Results: P[2 species] = 1.0000, migration rates recovered correctly (W1≈7.0, W2≈10.1, true=10)

---

## Part 2: Analysis of Acceptance Rate Issues

### 2.1 Observed Behavior

In test runs with moderate-to-high migration rates, the rj (split/join) acceptance rate is 0.0000. This occurs because:

1. Join proposals are immediately rejected when L↔R migrations exist
2. With high migration rates, L↔R migrations are almost always present
3. This creates asymmetry: splits are possible, joins are not

### 2.2 Root Causes

#### 2.2.1 Hard Rejection of L↔R Migrations (Critical)

```c
if (opt_migration && check_mig_join_validity(stree, gtree, node))
  return 2;  /* reject - L<->R migrations exist */
```

**Impact**: With M=10 and 50 loci, there are typically many L↔R migration events across loci, making joins essentially impossible.

**Mathematical issue**: L↔R migrations become self-migrations after join (source=target), which are undefined in the model.

#### 2.2.2 Fixed Beta(2,8) Parameters for Tau Proposal (Moderate)

```c
double pbetatau = 2;
double qbetatau = 8;
```

Mean = 2/(2+8) = 0.2, so tau_new ≈ 0.2 × tau_upper on average.

**Impact**: Proposals are concentrated in a narrow range; poor mixing when true tau differs significantly from 0.2 × tau_upper.

#### 2.2.3 No Migration Rate Adjustment (Moderate)

When populations split, new migration bands are created (e.g., L↔R, L↔external, R↔external), but no migration rates are proposed for these new bands. The reverse happens for joins.

**Impact**: Likelihood may drop substantially due to mismatched rates.

#### 2.2.4 Theta Proposal Assumptions (Minor)

Both theta proposal methods assume child effective sizes are similar to parent:
- Method 0: `child_theta = parent_theta * exp(epsilon * (U - 0.5))`
- Method 1: `child_theta ~ Gamma(alpha, alpha/(mean * parent_theta))`

**Impact**: Poor proposals when child populations have very different sizes.

---

## Part 3: Future Improvements

### 3.1 Priority 1: Migration-Aware Join Proposals

**Goal**: Allow join proposals to succeed even when L↔R migrations exist by properly removing them.

#### 3.1.1 Approach A: Migration Removal with Hastings Correction

Instead of hard rejection, compute the probability of the L↔R migration configuration and include in acceptance ratio.

```c
/* Proposed new function */
static double handle_migs_join_with_removal(stree_t *stree,
                                            gtree_t **gtree,
                                            snode_t *snode,
                                            double *lnacceptance)
{
  snode_t *left = snode->left;
  snode_t *right = snode->right;
  long removed_count = 0;

  for (i = 0; i < stree->locus_count; ++i)
  {
    for each gnode in gtree[i]:
      for each migration me in gnode->mi:
        if ((me->source == left && me->target == right) ||
            (me->source == right && me->target == left))
        {
          /* Store for potential restoration */
          store_for_removal(me, i);

          /* Remove migration event */
          migevent_unlink(me, i);
          remove_from_gnode_mi(gnode, me);

          /* Adjust counts */
          update_counts_for_removal(me, gtree[i], i);

          removed_count++;
        }
  }

  /* Hastings ratio: probability of proposing to add these back in split */
  /* This requires proposing migration events in split as well */
  *lnacceptance -= removed_count * log(proposal_prob_per_migration);

  /* Prior ratio: prior probability of these migrations */
  *lnacceptance -= compute_mig_prior_contribution(removed_migrations);

  return removed_count;
}
```

**Complexity**: Requires symmetric treatment in `prop_split()` to propose adding migrations.

#### 3.1.2 Approach B: Probabilistic Soft Rejection

Weight the acceptance probability by the cost of L↔R migrations rather than hard rejection:

```c
/* Instead of return 2 */
if (opt_migration)
{
  long lr_mig_count = count_lr_migrations(stree, gtree, node);
  if (lr_mig_count > 0)
  {
    /* Penalize but don't reject outright */
    /* Penalty reflects the "impossibility" of these migrations */
    lnacceptance -= lr_mig_count * LARGE_PENALTY;  /* e.g., 100.0 */
  }
}
```

**Simplicity**: Easy to implement but theoretically unsatisfying.

#### 3.1.3 Approach C: Migration Rerouting (Complex)

When L↔R migrations exist, reroute them through an external population if available:

```
Before join: L → R migration
After join:  (AB) → external → (AB) [invalid]
         OR  (AB) → external [if external accepts]
```

**Complexity**: Requires topological analysis and may not always be possible.

**Recommendation**: Start with Approach A for theoretical correctness.

### 3.2 Priority 2: Tunable Tau Proposal Parameters

**Goal**: Allow users to adjust the Beta distribution parameters for tau proposals.

#### 3.2.1 Control File Option

Add new option to control file:
```
rjtau = 2 8       # Beta(p, q) parameters for tau proposal
```

#### 3.2.2 Code Changes

In `src/bpp.h`:
```c
extern double opt_rjtau_alpha;  /* Beta p parameter, default 2 */
extern double opt_rjtau_beta;   /* Beta q parameter, default 8 */
```

In `src/prop_rj.c`:
```c
/* Replace hardcoded values */
double pbetatau = opt_rjtau_alpha > 0 ? opt_rjtau_alpha : 2.0;
double qbetatau = opt_rjtau_beta > 0 ? opt_rjtau_beta : 8.0;
```

#### 3.2.3 Adaptive Alternative

Track acceptance rates and adjust parameters automatically:
```c
static double adapt_beta_params(double current_p, double current_q,
                                double acceptance_rate, double target_rate)
{
  /* Increase variance if acceptance too low */
  /* Decrease variance if acceptance too high */
  if (acceptance_rate < target_rate * 0.5)
  {
    /* Proposals too narrow, increase variance */
    /* Move p and q closer to 1 */
  }
  else if (acceptance_rate > target_rate * 1.5)
  {
    /* Proposals too wide, decrease variance */
  }
}
```

### 3.3 Priority 3: Migration Rate Proposals for Split/Join

**Goal**: Propose migration rates for new bands during split, account for removed bands during join.

#### 3.3.1 Split: Propose New Migration Rates

When population P splits into L and R:
- Old bands: P↔X for each external population X
- New bands: L↔X, R↔X, L↔R

```c
static void propose_mig_rates_split(stree_t *stree, snode_t *node,
                                    double *lnacceptance)
{
  snode_t *left = node->left;
  snode_t *right = node->right;

  /* For each new migration band, propose rate from prior */
  for each new_band in {L↔R, L↔externals, R↔externals}:
  {
    double new_rate = sample_from_gamma(opt_mig_alpha, opt_mig_beta);
    set_migration_rate(new_band, new_rate);

    /* Jacobian: proposal density */
    *lnacceptance -= log_pdf_gamma(new_rate, opt_mig_alpha, opt_mig_beta);

    /* Prior contribution already handled by lnprior_species_model() */
  }
}
```

#### 3.3.2 Join: Account for Removed Rates

```c
static void account_mig_rates_join(stree_t *stree, snode_t *node,
                                   double *lnacceptance)
{
  /* For each band being removed, add proposal density to acceptance */
  for each removed_band in {L↔R, L↔externals, R↔externals}:
  {
    double old_rate = get_migration_rate(removed_band);

    /* Reverse Jacobian */
    *lnacceptance += log_pdf_gamma(old_rate, opt_mig_alpha, opt_mig_beta);
  }
}
```

### 3.4 Priority 4: Improved Theta Proposals

**Goal**: Better theta proposals that don't assume child sizes similar to parent.

#### 3.4.1 Data-Informed Proposals

Use sufficient statistics from gene tree coalescent times to inform theta proposals:

```c
static double propose_theta_from_data(snode_t *child, gtree_t **gtree,
                                      long locus_count)
{
  /* Compute mean coalescent rate in child population across loci */
  double sum_coal_times = 0;
  long coal_count = 0;

  for (i = 0; i < locus_count; ++i)
  {
    sum_coal_times += sum_of_coalescent_times_in_pop(child, gtree[i]);
    coal_count += child->coal_count[i];
  }

  /* Propose theta based on observed coalescent density */
  double estimated_theta = sum_coal_times / coal_count;  /* rough estimate */
  return sample_around(estimated_theta, epsilon);
}
```

#### 3.4.2 Wider Independent Proposals

Propose child thetas independently from the prior with wider variance:

```c
/* Method 2: Independent prior proposals */
if (opt_rjmcmc_method == 2)
{
  node->left->theta = sample_from_prior(opt_theta_alpha, opt_theta_beta);
  node->right->theta = sample_from_prior(opt_theta_alpha, opt_theta_beta);

  /* Jacobian is simply 1/prior_density for each */
  thetafactor /= pdf_prior(node->left->theta);
  thetafactor /= pdf_prior(node->right->theta);
}
```

### 3.5 Priority 5: Diagnostic Output

**Goal**: Help users understand why proposals are being rejected.

#### 3.5.1 Rejection Reason Tracking

```c
/* Add to prop_rj.c */
static long rj_reject_lr_mig = 0;      /* L<->R migration rejection */
static long rj_reject_theta_bounds = 0; /* theta outside bounds */
static long rj_reject_likelihood = 0;   /* likelihood ratio rejection */
static long rj_reject_tau_bounds = 0;   /* tau >= tau_upper */

/* In prop_join() */
if (opt_migration && check_mig_join_validity(stree, gtree, node))
{
  rj_reject_lr_mig++;
  return 2;
}

/* Report periodically */
void rj_print_diagnostics()
{
  printf("rjMCMC rejection reasons:\n");
  printf("  L<->R migrations: %ld\n", rj_reject_lr_mig);
  printf("  Theta bounds:     %ld\n", rj_reject_theta_bounds);
  printf("  Likelihood:       %ld\n", rj_reject_likelihood);
  printf("  Tau bounds:       %ld\n", rj_reject_tau_bounds);
}
```

#### 3.5.2 Debug Output Option

Add `opt_debug_rj_detailed` flag for verbose output:

```c
if (opt_debug_rj_detailed)
{
  printf("[rj-split] node=%s tau_new=%.6f tau_upper=%.6f\n",
         node->label, tau_new, tau_upper);
  printf("[rj-split] theta_L=%.6f theta_R=%.6f (parent=%.6f)\n",
         node->left->theta, node->right->theta, node->theta);
  printf("[rj-split] lnacceptance=%.4f (prior=%.4f, likelihood=%.4f)\n",
         lnacceptance, prior_contrib, likelihood_contrib);
}
```

---

## Part 4: Implementation Roadmap

### Phase 1: Diagnostics and Tuning (Low Risk)
- [ ] Add rejection reason tracking
- [ ] Add `opt_debug_rj_detailed` flag
- [ ] Make Beta parameters tunable via control file
- [ ] Add warning when rj acceptance rate is 0

### Phase 2: Migration Rate Proposals (Medium Risk)
- [ ] Implement `propose_mig_rates_split()`
- [ ] Implement `account_mig_rates_join()`
- [ ] Test with varying migration band configurations
- [ ] Verify detailed balance

### Phase 3: Migration-Aware Joins (High Complexity)
- [ ] Design L↔R migration removal mechanism
- [ ] Implement symmetric migration addition in split
- [ ] Implement proper Hastings ratio calculation
- [ ] Extensive testing with known true models
- [ ] Verify reversibility and detailed balance

### Phase 4: Advanced Proposals (Future)
- [ ] Data-informed theta proposals
- [ ] Adaptive tau proposal parameters
- [ ] Alternative proposal distributions

---

## Part 5: References

1. Yang Z, Rannala B (2010) Bayesian species delimitation using multilocus sequence data. PNAS 107:9264-9269. (YR2010 - original rjMCMC algorithm)

2. Flouri T, Jiao X, Rannala B, Yang Z (2018) Species tree inference with BPP using genomic sequences and the multispecies coalescent. Molecular Biology and Evolution 35:2585-2593.

3. Flouri T, Jiao X, Rannala B, Yang Z (2020) A Bayesian implementation of the multispecies coalescent model with introgression for phylogenomic analysis. Molecular Biology and Evolution 37:1211-1223. (MSC-I/MSC-M implementation)

---

## Appendix: Test Commands

```bash
# Simulate data with high migration
cd /home/bruce/repos/bpp/test/mig_delimit_test
../../src/bpp --simulate sim_veryhighmig.ctl

# Run species delimitation with migration
../../src/bpp --cfile infer_delimit_mig.ctl

# Check results
tail -20 delimit_mig.txt
```
