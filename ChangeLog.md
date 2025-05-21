# Change Log
All notable changes to `bpp` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [4.8.5] - 2025-05-21
## Added
 - clock 4 (lineage rate model)
## Fixed
 - processing of A1B1 crash due to missing columns (some thetas that cannot be
   estimated, i.e. there is one or no sequence at the species) are still
   required for converting W to M, but they were skipped.
 - disabled PDF visualization for the single population coalescent model

## [4.8.4] - 2025-02-11
## Added
 - usedata=2 for using fixed gene trees as data instead of sequences
 - Simulation of sequences with genotyping errors with mixed diploid and haploid
   sequences
 - Command-line option --keep-label
## Changed
 - Number of digits in a1b1 file
 - Printing inf for Effu and c if Effy=1
 - Alingments from compressed site patterns are now written in separate file
## Fixed
 - Crash due to mismatching header line in a1b1 file
 - Order of sampling W values in a1b1 file
 - Theta pjump calculation for mixture of gibbs and sliding window sampling by
   introducing move counters for each eps
 - Incorrect number of steplengths when using --theta_mode 3 under MSC-I model 3
 - Heredity scaler estimation when combined with gibbs sampler for thetas
 - Incorrect typecast in diploid code when tipdating

## [4.8.2] - 2024-11-28
## Changed
 - thetaprior: always estimate theta with invgamma prior unless 'int' specified
 - invgamma theta prior requires that a>2
## Fixed
 - duplicate logging of a1,b1 for theta and w under certain conditions
 - logging a1,b1 for phi even when sliding window is used
 - SNL move assertion fails with thetas integrated (called old density function)

## [4.8.1] - 2024-11-25
### Changed
 - a1b1 summary handles missing data
 - per-locus affected lists for extended rubberband algorithm
 - finetune option to use dictionary-like syntax
### Added
 - logging of a1,b1 parameters of theta conditionals when thetas integrated out
### Fixed
 - -nan for th1 caused by having one sequence per species and  differnt step
   lengths for inner and tip nodes
 - logging a1 and b1 for W parameters
 - theta distribution flag when summarizing conditional a1b1 file
 - when running --summary the conditional a1b1 file was being overwritten

## [4.8.0] - 2024-11-12
### Changed
 - Re-parameterization of the MSC-M model (M to w)
 - Summary statistics output format (transposed, aligned, parameters enumerated)
 - Summary for both w and M (MSC-M model)
 - Prior for HKY and K80
 - Only summarize unique theta parameters
 - Removed unnecessary sorting in several proposals
 - Replaced outfile and mcmcfile with a single jobname option
 - All output filenames now begin with string defined in jobname option
 - SeedUsed no longer created when invoking BPP without parameters
### Added
 - Gibbs move for sampling theta from conditional distribution
 - Gibbs move for sampling phi from conditional distribution
 - Gibbs move for sampling w from conditional distribution
 - Linked theta model for MSC-M
 - PDF visualization of species tree (MSC and MSC-M models)
 - Posterior summaries for parameters whose conditionals are tractable or can
   be approximated
 - Ported likelihood calculations for aarch64 neon instructions (apple silicon)
 - Option --extend to extend runs by additional MCMC steps
 - Options --phi-slide-prob and --theta-slide-prob to control the frequency of
   phi and theta sliding window moves
### Fixed
 - Option of integrating out thetas with linked theta models
 - Theta model "linked-msci" when a branch is broken into three or more segments
 - Illegal memory accesses when constructing faketree under linked theta models
 - Fixed efficiency calculation routine, added rho1 in output


## [4.7.0] - 2023-11-04
### Changed
 - Example data
 - Extended rubberband move is now the default move for MSCM
### Added
 - Pseudopriors for migration rates
 - Option to specify a separate prior for each migration rate
 - Per-locus variable migration rates
 - Reversible-jump MCMC move for adding and removing migration rates
 - Move to flip the direction of a migration rate
 - Option to change the probability of using a gibbs sampler vs sliding window
   proposal for theta moves (default: 0.2 sliding window)
 - Option for having one stepsize per theta parameter
 - Option for having one stepsize for tip nodes and one for inner nodes
 - Option for disabling core pinning
### Fixed
 - Memory leaks
 - Invalid access of migbuffer structure in theta gibbs proposal
 - Binary mode when opening checkpoint files

## [4.6.2] - 2022-11-03
### Changed
 - Default theta proposal to the original sliding window move
 - Default number of thetas on monitor output set back to 3
### Added
 - Gibbs sampler for thetas with t1/t4 conditions
### Fixed
 - Typos

## [4.6.1] - 2022-07-26
### Changed
 - Number of thetas,taus,phis and mrates on screen output
### Added
 - Gibbs sampler for migration rates
### Fixed
 - Priors on frogs example control files

## [4.6.0] - 2022-07-16
### Changed
 - Disabled twin towers unidentifiability algorithm when not using data
 - Printouts from label switching algorithm are now included in output file
 - Gibbs sampler is now the default move for thetas
### Fixed
 - Crash when resuming from checkpoint due to uninitialized IM variables
 - Checkpointing when using linked theta models
 - Copying of linked_theta attribute when cloning snodes -- fixes linked theta
   models under IM
### Added
 - Consistency checks for tau definitions in --simulate
 - Error when using MSCi and IM together
 - Possibility of displaying multiple phi means on screen output
 - Gibbs sampler for thetas under MSC,MSCi,IM with gamma prior

## [4.5.2] - 2022-06-27
### Fixed
 - Phi proposal for assymetric beta priors -- phi proposed on the branch the 
   parameter resides on
### Added
 - Gibbs sampler for population size parameters under the IM model with inverse
   gamma prior
 - Parallelization for freqs,qrates,alpha,brte moves when thetas integrated out
 - 'defphi' keyword in species tree attributes to define the branch on which
   the phi parameter is placed on

## [4.5.1] - 2022-05-22
### Changed
 - Migration rate symbol from small m to big M in MCMC sample file
 - Output seed value when seed<0
 - Qrates in simulation option now sum to 1
 - Order of qrates and base frequencies in the output file and added an
   indicative label for each substitution rate and base frequency
### Fixed
 - Edge conditions for phi label switching
 - Enabled --summary option  (previously disabled)
 - Bug that corrupted opt_tau_alpha when reading control file
 - Gene tree migration event move to save proposed time
 - Segfault in msci generator when adding time-travelling BDIs (issue #158)
 - Race condition due to migcount_sum
 - Memory leaks, improved reallocation of migbuffer, issues in miginfo_extend
 - Gene tree branch lengths in print_gene_trees output when using relaxed clocks
### Added
 - Indexing of migration events (improves performance)
 - Option --bfdriver (issue #129)
 - New rubberband algorithm under IM
 - Linked theta models (none,all,inner,msci)
 - Gene tree height and length in gene tree output file
 - Gibbs sampler for population sizes under the MSC and MSCi models with
   inverse gamma prior

## [4.5.0] - 2021-12-21
### Changed
 - Proposal kernel to bactrial laplace
 - Testbed generation hash
### Added
 - Isolation with migration model


## [4.4.1] - 2021-12-13
### Changed
 - Shortened MSci model information table
### Fixed
 - Crash when resuming from a checkpoint with SNL enabled
 - Assignment of mean phi values for model B on summarized species tree/network
 - Position of theta and tau in the header of the last output table
### Added
 - Option --summary for only summarizing the MCMC file
 - Algorithms CoG0,CoGN,BetaGamma for resolving identifiability issues
   associated with BDI events

## [4.4.0] - 2021-06-29
### Changed
 - Phi parameters (MSCi model) for non-BDI events always correspond to
   horizontal edges. If both parental edges of a hybridization node are
   horizontal or both are non-horizontal, then the phi parameter always
   corresponds to the mirror node.
 - Phi parameters (MSci model) for BDI events always correspond to the
   horizontal branches
 - MSci generator syntax accepts source branch first, then target branch
 - MSci generator prints phis on horizontal branches for BDI events
 - Numbering of loci in load balancing stats output starts from 1 instead of 0
 - Removed checks for consecutive unary nodes which fixes the rejection of some
   correctly constructed MSCi newick formats
### Fixed
 - Failing assertion caused by trying to update a non-existent root branch
   length for circle/triangle/square nodes in species tree spr with relaxed
   clock
 - Error message in simulation for locusrate syntax
 - MRCA function (for MSCi) in gene tree age proposal to find the MRCA
   population of two daughter nodes (a and b), such that the MRCA is also a
   descendant of the population of the parent of a and b. This fixes a failing
   assertion
 - Rule in MSCi creator that prevented from having multiple bidirections on
   sister species
 - Header of modelparafile when using relaxed clock models (printing of species
   tree branch rate labels)
 - Removed SNL branch rate swap strategies when species tree root changes
 - Removed branch rates swapping from sSPR and SNL when move is NNI

### Added
 - Zig-zag load balancing scheme
 - Removal of sequences that contain only missing data
 - Gene tree topology change move using simulation (option --exp_sim)
 - Gamma prior for tau
 - Gamma and Beta for theta

## [4.3.8] - 2020-12-02
### Changed
 - Initialization of phi values
### Fixed
 - Loading correct beta value for bayes factor computation when resuming from a
   checkpoint

## [4.3.7] - 2020-11-26
### Added
 - Added option --snl_noswap to test the strategy of overwriting new root rates
   when using correlated clock (SNL move)

## [4.3.6] - 2020-11-25
### Added
 - New mapping of branch rates on species tree for moves that correspond to NNI
### Changed
 - 'diploid' to 'phase' in frogs example control file

## [4.3.5] - 2020-10-29
### Added
 - Error message for number of threads
 - Code to generate 'faketree' (MSCi graph model represented as binary tree)
 - Swapping of branch rates (correlated clock) when SNL changes root node 
### Changed
 - Renamed 'diploid' to 'phase' in sumulation code
 - Updated README.md file for anopheles example

### Fixed
 - Simulation code to set thetas for species with 1 sequence and phase flag set

## [4.3.4] - 2020-09-07
### Changed
 - Removed command-line options for reject and repeat strategy (SNL move) and
   kept the reject strategy
 - All SNL command-line options were removed and made part of the control file
   (speciestree option)
 - Debugging function for validating log-PG
### Fixed
 - Acceptance proportion (SNL move)
 - Notheta option when using the MSCi model
 - Acceptance rate (pjump) output of SNL

## [4.3.3] - 2020-08-21
### Added
 - Options --snl_le and --snl_ls for setting lambdas (SNL move)
 - Printing of shrink proportion
### Fixed
 - Rejection and repeat options (SNL move)
 - Acceptance proportion for expand with downwards path (SNL move)

## [4.3.2] - 2020-08-17
### Added
 - Option for checking log-likelihood after each move
 - Debugging system
 - Options for repeated and rejection sampling for SNL move
### Fixed
 - Bug in locusrate proposal for correlated clock model using Dirichlet and when
   mubar is fixed to 1
 - Rules for detecting long branches in SSPR when using relaxed clocks
 - Swapping of p-matrices in mu_i and nu_i proposals
 - Acceptance ratio for SNL move
 - Memory leak associated with SNL move

## [4.3.1] - 2020-08-11
### Added
 - Added NOT keyword for defining topological constraints
 - Added '--comply' switch for checking compatibility of constraints against
   a tree
 - Added an initial implementation of the snakes and ladders move
### Changed
 - Simplified constraints definition (multiple compatible constraints may be
   specified)


## [4.3.0] - 2020-07-07
### Fixed
- Notation for simulating data with only one species
- Simulations with relaxed clock and log-normal distribution
- Mapping of sites from A2 to A3 in diploid compression which was causing the
  program to crash
- Crashes when user does not specify thetas in simulations. Error messages are
  now printed
### Added
- Specification of topological constraints for species tree and of outgroup

## [4.2.9] - 2020-04-28
### Fixed
 - Acceptance ratio for nu_i proposal with Dirichlet prior
 - nubar appears in summary statistics for Dirichlet prior
 - Parsing of custom model in control file
 - Memory leaks and segfault when checkpointing with molecular clock, GTR and
   per-locus file printing

## [4.2.8] - 2020-04-19
### Fixed
 - Conditions for enabling nubar estimation
 - Recomputation of log-L for correlated clock when proposing mu_i (full 
 - Computation of log prior ratio for branch rates with correlated clock and
   log-normal distribution
 - Branch rates prior for correlated clock and gamma distribution
 - Swapping of pmatrices when proposing mu_i
 - Deallocation problem when phylip alignment contains less sequences than
   specified in header
 - Acceptance ratio for correlated clock and log-normal distribution when
   proposing nu_i (extra variance term)

## [4.2.7] - 2020-04-04
### Changed
 - Species numbering in pptable starts from 1 instead of 0
 - Reduced spacing in finetune adjustment output
 - Finetune option now accepts dashes as entries indicating to use default the
   default value for that step length
### Fixed
 - Log prior ratio for correlated clock with gamma distribution
 - Update branch rate prior in mixing move when using correlated clock and gamma
   distribution

## [4.2.6] - 2020-03-28
### Fixed
 - Made arch option CPU to be case insensitive
 - Computation of branch lengths for gamma rates heterogeneity

## [4.2.5] - 2020-03-26
### Changed
 - Clock prior specification in (--simulate) from 0 and 1 to DIR and IID
### Fixed
 - Bug in compressing patterns when GTR and diploid sequences; compression was
   still done as if JC69 was used, instead of compressing only unique patterns
 - Nucleotide map was used to translate the characters of AA data when using
   partitioned analyses
 - pmatrix computation did not account for gamma rates heterogeneity

## [4.2.4] - 2020-03-20
### Added
 - Checks for -inf log-L and prompt for enabling scaling to prevent numerical
   underflow
### Changed
 - Clock and locusrate prior specification from 0 and 1 to DIR and IID


## [4.2.3] - 2020-03-14
### Added
 - Header line on current pjump finetune and new finetune
### Changed
 - Rearranged order of finetune arguments in control file and included alpha,
   freqs,qmat steplengths
### Fixed
 - Corrected typecasting of returned pjump for mu_i and nu_i proposals

## [4.2.2] - 2020-03-12
### Added
 - Compatible relaxed clock models in the BPP simulator
 - Screen output is now printed in output file as well
 - Correlated relaxed clock models
 - Parallelized branch rate proposal
 - Implemented MSci notation generator from a species tree and a list of
   edge (option --msci-create)
### Changed
 - New parser for Imap file (not using flex/bison anymore)
 - mubar and vbar no longer printed in finetune output if not used
 - Clock arguments for specifying branch rates distribution from 0 and 1 to
   LN and G
 - Renamed option diploid to phase
 - Prettified screen output
### Fixed
 - Inconsistency between current finetune and pjump in output
 - MSci model concerning nodes participating in two hybridization events with
   conflicting htau causing the program to enter an infinite loop before
   starting

## [4.2.1] - 2020-01-31
### Added
 - Additional print flag for printing qmatrix parameters, frequencies and alpha
   value into a locus specific file
 - Custom pmatrix computation for different nucleotide models and output of
   their in locus-specific files
### Changed
 - Implemented new general tree parser (not using flex/bison anymore)
 - Starting values for phi are not needed anymore in inference mode, and are
   used as starting values if specified
 - Updated Anopheles example documentation
 - Locus rates (mu_i) and heredity scalars are now printed in locus specific
   files
### Fixed
 - branch rates prior computation
 - Qrates and frequencies proposal acceptance ratio
 - Branch length computation for relaxed clock

## [4.2.0] - 2019-10-29
### Added
 - Rates proposal for nucleotide GTR model
 - Computation of empirical frequencies
 - Summary table for loci
 - Option alphaprior=X Y Z and roposal for alpha parameter (site rate variation)
 - Yeast dataset
 - Relaxed clock (iid rates model) following gamma or lognormal
 - Relaxed clock (Gamma-Dir model) following gamma or lognormal
### Changed
 - Finetune option now accepts arbitrary number of parameters (non-specified
   are set to default values)
 - Format of locusrate option
### Fixed
 - Checkpoint file now stores phi pjump
 - Parsing of species tree with one species and hybridization when simulating

## [4.1.4] - 2019-05-22
### Added
 - All protein models are now activated and can be used with the 'model' option
 - Option to set a different model for each locus by using model = custom
### Changed
 - MCMC output file columns are now aligned
 - 'threads' option can take two optional arguments indicating a) starting core
   index to pin threads to, and b) a stride.
 - Removed progress indicator
 - Reading seed from /dev/urandom when simulating data
 - Finetuning results are stored also in output file
### Fixed
  - Checkpointing now works when bidirectional introgressions are present

## [4.1.3] - 2019-03-15
### Added
 - Printing starting time/date of analysis, BPP version and command-line
   arguments at the beginning of the output file
 - Printing pattern weights for compressed alignments in output file, and
   sequneces are now printed in aligned form
 - Added check that Inverse-Gamma priors are indeed used. If prior mean is
   greater than 1 then BPP complains.

## [4.1.2] - 2019-02-22
### Added
 - Pinning threads to cores on linux systems. Improved multithread performance

## [4.1.1] - 2019-01-30
### Added
 - Parallelized mixing and tau proposals
 - Experimental option for randomizing nodes order before gene tree SPR move
   (switch --exp\_random)
 - Revolutionary gene tree spr move (switch --rev\_gspr)
### Changed
 - Changed 'gammaprior' to 'phiprior' and 'gamma' to 'phi' in MSci model

## [4.1.0] - 2019-01-15
## Added
 - Parallelized gene tree SPR and gene tree age proposal
## Fixed
 - Minor bug when cleandata=1 and locus has no ambiguous characters

## [4.0.7] - 2019-01-06
### Added
 - Code for generating full data and random resolution files when using diploid
   option in simulations
### Changed
 - Random number generator
 - Format of diploid sequence labels in simulations code

## [4.0.6] - 2019-01-04
### Fixed
 - Variable names in QuantileChi2() to conform with MSVC compiler
 - Assignment of thetas and taus to hybridization nodes in simulations code

## [4.0.5] - 2019-01-02
### Added
 - Simulations (MCcoal) code via the --simulate switch
 - Added multispecies coalescent introgression (MSci) model

## [4.0.4] - 2018-10-18
### Fixed
  - Fixed minor bug in gene tree spr move

## [4.0.3] - 2018-10-04
### Fixed
  - Numerical boundary problems in reflect
  - Usage of exptm1 function in JC69 pmatrix computation code for preventing
    numerical problems
### Added
  - Estimation of thetas (method A00) for one species only

## [4.0.2] - 2018-09-14
### Fixed
 - When resuming from a checkpoint, additional future expected checkpoints were
   no longer created. This is fixed now
## [4.0.1] - 2018-08-30
### Added
 - Scaling option in control file that prevents numerical underflow when
   calculating partial likelihoods
### Fixed
 - Fixed qsort callback function in A00 summary to adher to BSD qsort
