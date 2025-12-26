# Piecewise-Constant Population Size (Theta) with Variable Divergence Time (Tau)

## Summary of Theory, Implementation, and Preliminary Results

---

## 1. Introduction and Motivation

In Bayesian phylogenetic inference under the multispecies coalescent (MSC), population size parameters (θ) are typically assumed constant through time for each population. However, population sizes may vary across different time periods. This implementation extends BPP to allow **piecewise-constant theta values** across user-defined time intervals.

A key challenge arises when the divergence time (τ) is estimated rather than fixed: different populations exist during different time intervals, and the set of "real" parameters changes depending on where τ falls. We address this using the **Carlin-Chib product-space formulation**, which avoids the complexity of reversible-jump MCMC by maintaining all parameters in a single augmented state space.

---

## 2. Theoretical Framework

### 2.1 Piecewise-Constant Theta Model

For a species tree with K global time intervals defined by boundaries:
```
0 = t₀ < t₁ < t₂ < ... < t_K = ∞
```

Each population has K theta parameters, one for each interval:
```
θ_pop = (θ_pop[1], θ_pop[2], ..., θ_pop[K])
```

The MSC density for coalescent events in interval k uses θ_pop[k].

### 2.2 The Variable Tau Problem

Consider a 2-species case with species A, B and ancestor AB:
- Species A and B exist from time 0 to τ_AB
- Ancestor AB exists from time τ_AB to ∞

When τ is **fixed**, determining which intervals each population occupies is straightforward. When τ is **estimated**, the model dimension effectively changes as τ crosses interval boundaries.

**Example with K=3 intervals [0, t₁), [t₁, t₂), [t₂, ∞):**

| τ location | A real intervals | B real intervals | AB real intervals | Total real params |
|------------|------------------|------------------|-------------------|-------------------|
| τ ∈ [0, t₁) | [0, t₁) only | [0, t₁) only | [0, t₁), [t₁, t₂), [t₂, ∞) | 1+1+3 = 5 |
| τ ∈ [t₁, t₂) | [0, t₁), [t₁, t₂) | [0, t₁), [t₁, t₂) | [t₁, t₂), [t₂, ∞) | 2+2+2 = 6 |
| τ ∈ [t₂, ∞) | [0, t₁), [t₁, t₂), [t₂, ∞) | [0, t₁), [t₁, t₂), [t₂, ∞) | [t₂, ∞) only | 3+3+1 = 7 |

### 2.3 Carlin-Chib Product-Space Formulation

Instead of using reversible-jump MCMC (which requires carefully designed dimension-matching proposals), we use the **Carlin-Chib approach**:

1. **Product Space**: Maintain ALL 3K theta parameters in the state space at all times
2. **Real vs Pseudo Parameters**:
   - **Real parameters**: Contribute to the likelihood (population exists in that interval)
   - **Pseudo parameters**: Do not affect likelihood (population doesn't exist in that interval)
3. **Pseudo-prior**: Pseudo parameters are sampled from the prior distribution
4. **Model Indicator**: k = which interval contains τ

**Key insight**: The posterior over the product space correctly marginalizes to give the correct within-model posteriors, and the model posterior probabilities P(M_k | data) are estimated by counting samples.

### 2.4 MCMC Updates

**Theta Gibbs Sampling (modified):**
```
For each interval i in population:
    if θ[i] is REAL:
        Sample from conjugate posterior: θ[i] ~ InvGamma(α + n_coal, β + C₂h)
    else (PSEUDO):
        Sample from prior: θ[i] ~ InvGamma(α, β)
```

**Tau Proposals (modified):**
```
After accepting new τ:
    Determine new model k' = which interval contains τ
    Update real/pseudo status for all populations
    (No Jacobian needed - product space handles dimension change)
```

---

## 3. Implementation Details

### 3.1 New Control File Options

```
theta_intervals = K t0 t1 t2 ... tK
```
- K: number of intervals
- t0, t1, ..., tK: interval boundaries (t0 should be 0, tK typically 1.0 or large)

```
theta_variable_tau = 1
```
- Enables Carlin-Chib model switching when τ is estimated
- When disabled (default), all intervals are treated as real

### 3.2 Data Structure Extensions

**In `snode_t`:**
```c
int * theta_intv_real;        /* [K] boolean: 1=real, 0=pseudo */
int theta_intv_first_real;    /* first real interval index */
int theta_intv_last_real;     /* last real interval index */
```

**In `stree_t`:**
```c
long * model_count;           /* [K] samples per model */
int current_model_k;          /* current model indicator (which interval contains root τ) */
```

### 3.3 Key Functions Added/Modified

| File | Function | Purpose |
|------|----------|---------|
| stree.c | `get_model_k()` | Determine which interval contains τ |
| stree.c | `update_theta_real_pseudo_status()` | Set real/pseudo flags for a node |
| stree.c | `update_all_theta_real_pseudo_status()` | Update all nodes after τ change |
| stree.c | `propose_theta_intv_gibbs()` | Modified to sample pseudo from prior |
| stree.c | `propose_tau()` | Updated to call status update on acceptance |
| gtree.c | `gtree_update_logprob_contrib()` | Skip pseudo intervals in likelihood |
| method.c | MCMC loop | Model counting and display |
| method.c | Summary output | Model posterior probabilities |

### 3.4 Output

**During MCMC (screen display):**
```
Mk:pp - current model and its posterior probability
Example: M2:0.97 means τ is in interval 2, with 97% posterior probability
```

**End of run:**
```
Model Posterior Probabilities (tau interval):
  Model 1 [0.000000, 0.002000): P = 0.0000 (samples = 0)
  Model 2 [0.002000, 0.005000): P = 1.0000 (samples = 20000)
  Model 3 [0.005000, 1.000000): P = 0.0000 (samples = 0)
```

---

## 4. Preliminary Analyses

### 4.1 Test 1: Single-Species Validation

**Purpose:** Verify theta_variable_tau=1 works correctly when there's only one population (no model switching should occur).

**Setup:**
- 1 species (Human), 64 sequences
- 50 loci
- theta_intervals = 3: [0, 0.0001), [0.0001, 0.0005), [0.0005, 1.0)
- theta_variable_tau = 1

**Result:** All samples in Model 1 (M1:1.00), as expected for single species with no divergence time.

### 4.2 Test 2: Two-Species Basic Test

**Purpose:** Verify implementation with 2-species data where τ is estimated.

**Simulation Parameters:**
- 2 species (A, B) with ancestor AB
- 10 sequences per species, 20 loci, 3000 bp each
- True values: τ_AB = 0.002, θ_A = 0.001, θ_B = 0.002, θ_AB = 0.003

**Inference Setup:**
- theta_intervals = 3: [0, 0.001), [0.001, 0.005), [0.005, 1.0)
- theta_variable_tau = 1
- Priors: θ ~ InvGamma(3, 0.002), τ ~ InvGamma(3, 0.004)
- MCMC: 4000 burnin, 10000 samples

**Results:**

| Parameter | True | Estimated |
|-----------|------|-----------|
| τ_AB | 0.002 | 0.00195 |
| θ_A | 0.001 | 0.00095 |
| θ_B | 0.002 | 0.00195 |
| θ_AB | 0.003 | 0.00373 |

**Model Posteriors:**
- Model 1 [0.000, 0.001): P = 0.028
- Model 2 [0.001, 0.005): P = 0.972 ✓ (true τ=0.002 is in this interval)
- Model 3 [0.005, 1.000): P = 0.000

### 4.3 Test 3: Discordant Theta Values - Consistency Check

**Purpose:** Test with more distinct theta values and verify consistency across independent runs.

**Simulation Parameters:**
- 2 species (A, B) with ancestor AB
- 20 sequences per species, 50 loci, 2000 bp each
- True values: τ_AB = 0.003, θ_A = 0.002, θ_B = 0.004, θ_AB = 0.008

**Inference Setup:**
- theta_intervals = 3: [0, 0.002), [0.002, 0.005), [0.005, 1.0)
- theta_variable_tau = 1
- Priors: θ ~ InvGamma(3, 0.008), τ ~ InvGamma(3, 0.006)
- MCMC: 8000 burnin, 20000 samples
- Two independent runs with different seeds (11111 and 99999)

**Results:**

| Parameter | True | Run 1 (seed=11111) | Run 2 (seed=99999) |
|-----------|------|-------------------|-------------------|
| τ_AB | 0.003 | 0.0033 | 0.0032 |
| θ_A | 0.002 | 0.00191 | 0.00191 |
| θ_B | 0.004 | 0.00403 | 0.00403 |
| θ_AB (intv 2) | 0.008 | 0.00669 | 0.00731 |

**Model Posteriors (both runs):**
- Model 1 [0.000, 0.002): P = 0.0000
- Model 2 [0.002, 0.005): P = 1.0000 ✓ (true τ=0.003 is in this interval)
- Model 3 [0.005, 1.000): P = 0.0000

**Conclusions:**
1. Both runs show excellent consistency in parameter estimates
2. Model 2 correctly receives 100% posterior probability
3. θ_A and θ_B estimates are very close to true values
4. θ_AB is somewhat underestimated but within reasonable range
5. The implementation produces reproducible results

---

## 5. Control File Examples

### 5.1 Simulation Control File
```
seed = 54321

species&tree = 2 A B
                  20 20
                  (A #0.002, B #0.004): 0.003 #0.008;

imapfile = imap.txt
seqfile = simulated.phy
loci&length = 50 2000
clock = 1
locusrate = 0
model = 0
```

### 5.2 Inference Control File
```
seed = 11111

seqfile = simulated.phy
Imapfile = imap.txt
jobname = output

speciesdelimitation = 0
speciestree = 0
species&tree = 2 A B
               20 20
               (A, B): 0.003;

theta_intervals = 3 0 0.002 0.005 1.0
theta_variable_tau = 1

thetaprior = invgamma 3 0.008
tauprior = invgamma 3 0.006

phase = 0
nloci = 50
usedata = 1
cleandata = 0

finetune = 1
print = 1 0 0 0
burnin = 8000
sampfreq = 2
nsample = 20000
```

---

## 6. Limitations and Future Work

### Current Limitations:
1. **2-species only**: Current implementation focuses on 2-species case (A, B, AB)
2. **Single τ**: Only handles single divergence time (root τ)
3. **No interval-specific simulation**: Simulation doesn't yet support different θ per interval

### Future Extensions:
1. Extend to arbitrary species trees with multiple divergence times
2. Each internal node τ would define its own model indicator
3. Support for interval-specific theta in simulation
4. Model-conditional summary statistics for each interval

---

## 7. References

1. Carlin, B.P. and Chib, S. (1995). Bayesian model choice via Markov chain Monte Carlo methods. *Journal of the Royal Statistical Society: Series B*, 57(3), 473-484.

2. Flouri, T., et al. (2018). A Bayesian implementation of the multispecies coalescent model with introgression for phylogenomic analysis. *Molecular Biology and Evolution*, 35(7), 1630-1640.

---

## 8. Files

### Source Code (committed):
- `src/bpp.h` - Data structure definitions
- `src/bpp.c` - Global variables
- `src/stree.c` - Tree operations and theta/tau proposals
- `src/gtree.c` - Gene tree likelihood calculations
- `src/method.c` - MCMC loop and output
- `src/cfile.c` - Control file parsing
- `src/cfile_sim.c` - Simulation control file parsing
- `src/treeparse.c` - Tree parsing and memory management
- `src/allfixed.c` - Fixed-parameter analysis

### Test Files (in bpp_popsize_test/):
- `sim_discordant.ctl` - Simulation control file
- `imap_discordant.txt` - Individual-to-species mapping
- `infer_discordant_run1.ctl` - Inference run 1
- `infer_discordant_run2.ctl` - Inference run 2
- `sim_discordant.phy` - Simulated sequence data
- `infer_discordant_run1.txt` - Run 1 output
- `infer_discordant_run2.txt` - Run 2 output
