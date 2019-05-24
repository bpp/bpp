# Change Log
All notable changes to `bpp` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
