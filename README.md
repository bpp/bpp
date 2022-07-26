# bpp

[![Build Status](https://travis-ci.org/bpp/bpp.svg?branch=master)](https://travis-ci.org/bpp/bpp)
[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/DOI/10.1093/molbev/msy147.svg)](https://doi.org/10.1093/molbev/msy147)
[![Version](https://img.shields.io/badge/version-4.6.1-blue.svg)](https://github.com/bpp/bpp/releases/tag/v4.6.1)

## Introduction

The aim of this project is to implement a versatile high-performance version
of the BPP software. It should have the following properties:

* open-source code with an appropriate open-source license.
* 64-bit multi-threaded design that handles very large datasets.
* easy to use and well-documented.
* SIMD implementations of time-consuming parts.
* Linux, Mac and Microsoft Windows compatibility.

BPP currently implements four methods:

* Estimation of the parameters of species divergence times and population sizes
  under the multi-species coalescent (MSC) model when the species phylogeny is
  given (Rannala and Yang, 2003)

* Inference of the species tree when the assignments are given by the user
  (Rannala and Yang, 2017)

* Species delimitation using a user-specified guide tree (Yang and Rannala,
  2010; Rannala and Yang, 2013)

* Joint species delimitation and species tree estimation (Yang and Rannala 2014)

BPP can also accommodate variable mutation rates among loci (Burgess and Yang,
2008) and heredity multipliers (Hey and Nielsen, 2004).  Finally, BPP supports
diploid data. Phasing is done analytically as described by Gronau et al, 2011.


**Update 2019-01-02**

BPP now also implements the Multispecies-coalescent-with-introgression (MSci) model
(see Flouri et al, 2020), an extension of the multispecies coalescent model to
incorporate introgression/hybridization.
For more information on usage please see the [BPP manual](https://github.com/bpp/bpp/raw/master/bppDOC.pdf) and/or the [wiki documentation](https://github.com/bpp/bpp/wiki).

**Update 2020-04-20**

Since v4.3.8 BPP implements the following:
* Parallelized computation using POSIX threads.
* Nucleotide substitution models JC69, K80, F81, HKY, T92, TN93, F84 and GTR.
* 19 amino acid substition models (Dayhoff, LG, DCMUT, JTT, etc)
* Site rate variation (+Gamma)
* Partitioned analyses
* Relaxed clock (independent rates model and correlated clock) as described in (Zhu et al, 2015).
* Easy to use tool for creating MSci networks that accepts a species tree and a list of definitions
  specifying introgressions/hybridizations.
For information on using those options check the [BPP manual](https://github.com/bpp/bpp/raw/master/bppDOC.pdf)

## Download and install

**Binary distribution** Starting with version 4.1.3, binary distribution files
containing pre-compiled binaries will be available as part of each
[release](https://github.com/bpp/bpp/releases). The included executables are
statically compiled whenever possible such that no library dependencies are
necessary.

Binary distributions are provided for x86-64 systems running GNU/Linux, macOS
(version 10.13 or higher) and Windows (64-bit, version 7 or higher).

Download the appropriate executable for your system using the following
commands if you are using a Linux x86_64 system:

```bash
wget https://github.com/bpp/bpp/releases/download/v4.6.1/bpp-4.6.1-linux-x86_64.tar.gz
tar zxvf bpp-4.6.1-linux-x86_64.tar.gz
```

Or these commands if you using a Mac:

```bash
wget https://github.com/bpp/bpp/releases/download/v4.6.1/bpp-4.6.1-macos-x86_64.tar.gz
tar zxvf bpp-4.6.1-macos-x86_64.tar.gz
```

Or if you are using Windows, download and extract (unzip) the contents of this file:

```
https://github.com/bpp/bpp/releases/download/v4.6.1/bpp-4.6.1-win-x86_64.zip
```

Linux and Mac: You will now have the binary distribution in a folder called `bpp-4.6.1-linux-x86_64` or `bpp-4.6.1-macos-x86_64`. The binary file is located in the `bin` subfolder, i.e. `bin/bpp`. We recommend making a copy or a symbolic link to the binary in a folder included in your `$PATH`. 

Windows: You will now have the binary distribution in a folder called `bpp-4.6.1-win-x86_64`. The bpp executable is called `bpp.exe`.


**Compiling from source** You can either download the *source distribution* for a particular version or *clone the repository*.

**Source distribution** To download the source distribution from a [release](https://github.com/bpp/bpp/releases) and build the executable and documentation, use the following commands:

```
wget https://github.com/bpp/bpp/archive/v4.6.1.tar.gz
tar zxvf v4.6.1.tar.gz
cd bpp-4.6.1/src
make
```

**Cloning the repo** Instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below.

```
git clone https://github.com/bpp/bpp.git
cd bpp/src
make
```

Compiling BPP requires that your system has GCC version 4.7 or newer, as
[AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions)  and AVX-2
optimized functions are compiled even if your processor does not support them.
This is fine, as BPP will automatically select the right instruction set that
your processor supports at run-time. This means, you can compile on one system,
and run BPP on any other compatible system.

However, if your compiler is older than 4.7, you will get errors such as:

```bash
cc1: error: unrecognized command line option "-mavx2"
```
or

```bash
cc1: error: unrecognized command line option "-mavx"
```

If your compiler is GCC 4.6.x then you can compile BPP using:

```bash
make clean
make -e DISABLE_AVX2=1
```

In case your compiler is older than GCC 4.6 then compile using:

```bash
make -e DISABLE_AVX2=1 DISABLE_AVX=1
```

You can check your compiler version with:
```bash
gcc --version
```

## Running BPP

After creating the control file, one can run BPP as follows:

```bash
bpp --cfile [CONTROL-FILE]
```

If you would like to resume a checkpoint file, please run:

```bash
bpp --resume [CHECKPOINT-FILE]
```

If you would like to run the simulator (previously MCcoal), please run:

```bash
bpp --simulate [CONTROL-FILE]
```

If you would like to run the MSci network generator, please run:

```bash
bpp --msci-create [DEFS-FILE]

```

If you would like to only summarize the results of an analysis, please run:
```bash
bpp --summary [CONTROL-FILE]
```


For an example of a DEFS-FILE see the [MSci generator notes](https://github.com/bpp/bpp/releases/download/v4.6.1/msci-create.pdf)

More documentation regarding control files, will be available soon on the [wiki](https://github.com/bpp/bpp/wiki).

## Documentation

The most up-to-date documentation of BPP is [bppDOC.pdf](https://github.com/bpp/bpp/releases/download/v4.6.1/bppDOC.pdf) distribution together with BPP.

A tutorial on BPP was recently published as a book chapter:
[A Tutorial on the Use of BPP for Species Tree Estimation and Species Delimitation](https://hal.inria.fr/PGE/hal-02536475)

For information on the MSci model please read:

* Flouri T., Jiao X., Rannala B., Yang Z. (2020) **A Bayesian Implementation of the Multispecies Coalescent Model with Introgression for Phylogenomic Analysis.** *Molecular Biology and Evolution*, 37(4):1211-1223.

* Jiao X., Flouri T., Rannala B., Yang Z. (2020) **The Impact of Cross-Species Gene Flow on Species Tree Estimation.** *Systematic Biology* (in press). doi:[10.1093/sysbio/syaa001](https://doi.org/10.1093/sysbio/syaa001)

## Citing BPP

Please cite the following publication if you use BPP:

Flouri T., Jiao X., Rannala B., Yang Z. (2018) **Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent.** *Molecular Biology and Evolution*, 35(10):2585-2593.
doi:[10.1093/molbev/msy147](https://doi.org/10.1093/molbev/msy147)

Please note, citing the corresponding of the four underlying methods, may also be appropriate.

If you use the MSci model please also cite:

Flouri T., Jiao X., Rannala B., Yang Z. (2020) **A Bayesian Implementation of the Multispecies Coalescent Model with Introgression for Phylogenomic Analysis.** *Molecular Biology and Evolution*, 37(4):1211-1223.
doi:[10.1093/molbev/msz296](https://doi.org/10.1093/molbev/msz296)


## License and third party licenses

The code is currently licensed under the [GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).

## The team

* Tom&aacute;&scaron; Flouri
* Bruce Rannala
* Ziheng Yang

## Code

| File                       | Description                                                                       |
| -------------------------- | --------------------------------------------------------------------------------- |
| **arch.c**                 | Architecture specific code (Linux/Mac/Windows)                                    |
| **allfixed.c**             | Summary statistics for method A00 (fixed species tree)                            |
| **bpp.c**                  | Main file handling command-line parameters and executing selected methods         |
| **bpp.h**                  | BPP header file including function prototypes and data structures                 |
| **cfile.c**                | Functions for parsing the control file                                            |
| **cfile_sim.c**            | Functions for parsing the control file (simulation mode)                          |
| **compress.c**             | Functions for compressing multiple sequence alignments into site patterns         |
| **constraint.c**           | Functions for placing topological constraints on species tree                     |
| **core_likelihood.c**      | Core functions for evaluating the likelihood of a tree (non-vectorized)           |
| **core_likelihood_avx.c**  | Core functions for evaluating the likelihood of a tree (AVX version)              |
| **core_likelihood_avx2.c** | Core functions for evaluating the likelihood of a tree (AVX-2 version)            |
| **core_likelihood_sse.c**  | Core functions for evaluating the likelihood of a tree (SSE-3 version)            |
| **core_partials.c**        | Core functions for computing partial likelihoods (non-vectorized)                 |
| **core_partials_avx.c**    | Core functions for computing partial likelihoods (AVX version)                    |
| **core_partials_avx2.c**   | Core functions for computing partial likelihoods (AVX-2 version)                  |
| **core_partials_sse.c**    | Core functions for computing partial likelihoods (SSE-3 version)                  |
| **core_pmatrix.c**         | Core functions for constructing the transition probability matrix                 |
| **debug.c**                | Functions for debugging purposes                                                  |
| **delimit.c**              | Species delimitation auxiliary functions and summary statistics                   |
| **diploid.c**              | Functions for resolving/phasing diploid sequences                                 |
| **dlist.c**                | Functions for handling doubly linked-lists                                        |
| **dump.c**                 | Functions for dumping the MCMC state into a checkpoint file                       |
| **gamma.c**                | Functions for obtaining rates from a discretized Gamma distribution               |
| **gtree.c**                | Functions for setting and processing gene trees                                   |
| **hardware.c**             | Functions for hardware detection                                                  |
| **hash.c**                 | Hash table implementation and related functions                                   |
| **list.c**                 | Linked list implementation and related functions                                  |
| **load.c**                 | Functions for loading a checkpoint file                                           |
| **locus.c**                | Locus specific functions                                                          |
| **lswitch.c**              | Algorithms for resolving identifiability issues associated with BDI events        |
| **Makefile**               | Makefile                                                                          |
| **mapping.c**              | Functions for handling map files                                                  |
| **maps.c**                 | Character mapping arrays for converting sequences to the internal representation  |
| **method.c**               | Function containing the MCMC loop and calls to proposals                          |
| **msa.c**                  | Code for processing multiple sequence alignments                                  |
| **msci_gen.c**             | Functions for the MSci generator                                                  |
| **output.c**               | Auxiliary functions for printing pmatrices (to-be-renamed)                        |
| **parsemap.c**             | Functions for parsing map files                                                   |
| **phylip.c**               | Functions for parsing phylip files                                                |
| **prop_gamma.c**           | Functions for proposing site rates                                                |
| **prop_mixing.c**          | Functions for the mixing proposal                                                 |
| **prop_rj.c**              | Functions for the reversible-jumps MCMC proposals for species delimitation        |
| **random.c**               | Pseudo-random number generator functions                                          |
| **revolutionary.c**        | Experimental functions for new (r)evolutionary algorithms                         |
| **rtree.c**                | Species tree export functions (to-be-renamed).                                    |
| **simulate.c**             | Functions for the simulation program (MCcoal)                                     |
| **stree.c**                | Functions for setting and processing the species tree                             |
| **summary.c**              | Species tree inference summary related functions                                  | 
| **summary11.c**            | Functions for summarizing joint species tree inference and delimitation           |
| **threads.c**              | Functions for parallelizing computation using POSIX threads                       |
| **treeparse.c**            | Functions for parsing trees                                                       |
| **util.c**                 | Various common utility functions                                                  |

# Acknowledgements

Special thanks to:
 * Paul M. Hime
 * Jiayi Ji
 * Paschalia Kapli
 * Mario dos Reis Barros
 * Yuttapong Thawornwattana,

for testing and bug reports.

# References

* Burgess R., Yang Z. (2008)
**Estimation of hominoid ancestral population sizes under Bayesian coalescent models incorporating mutation rate variation and sequencing errors.**
*Molecular Biology and Evolution*, 25(9):1979-1994.
doi:[10.1093/molbev/msn148](http://dx.doi.org/10.1093/molbev/msn148)

* Flouri T., Carrasco FI, Darriba D., Aberer AJ, Nguyen LT, Minh BQ, Haeseler A., Stamatakis A. (2015)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2):356-362.
doi:[10.1093/sysbio/syu084](http://dx.doi.org/10.1093/sysbio/syu084)

* Flouri T., Jiao X., Rannala B., Yang Z. (2018)
**Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent.**
*Molecular Biology and Evolution*, 35(10):2585-2593.
doi:[10.1093/molbev/msy147](http://dx.doi.org/10.1093/molbev/msy147)

* Flouri T., Jiao X., Rannala B., Yang Z. (2020)
**A Bayesian Implementation of the Multispecies Coalescent Model with Introgression for Phylogenomic Analysis.**
*Molecular Biology and Evolution*, 37(4):1211-1223.
doi:[10.1093/molbev/msz296](https://doi.org/10.1093/molbev/msz296)

* Gronau I., Hubisz MJ, Gulko B., Danko CG, Siepel A. (2011)
**Bayesian inference of ancient human demography from individual genome sequences.**
*Nature Genetics*, 43(10):1031-1035.
doi:[10.1038/ng.937](http://dx.doi.org/10.1038/ng.937)

* Hey J., Nielsen R. (2004)
**Multilocus methods for estimating population sizes, migration rates and divergence time, with applications to the divergence of Drosophila pseudoobscura and D. persimilis.**
*Genetics*, 167(2):747-760.
doi:[10.1534/genetics.103.024182](http://dx.doi.org/10.1534/genetics.103.024182)

* Jiao X., Flouri T., Rannala B., Yang Z. (2020)
**The Impact of Cross-Species Gene Flow on Species Tree Estimation.**
*Systematic Biology* (in press).
doi:[10.1093/sysbio/syaa001](https://doi.org/10.1093/sysbio/syaa001)

* Rannala B., Yang Z. (2013)
**Improved reversible jump algorithms for Bayesian species delimitation.**
*Genetics*, 194:245-253.
doi:[10.1534/genetics.112.149039](http://dx.doi.org/10.1534/genetics.112.149039)

* Rannala B., Yang Z. (2017)
**Efficient Bayesian Species Tree Inference under the Multispecies Coalescent.**
*Systematic Biology*, 66(5):823-842.
doi:[0.1093/sysbio/syw119](http://dx.doi.org/0.1093/sysbio/syw119)

* Yang Z., Rannala B. (2003)
**Bayes Estimation of Species Divergence Times and Ancestral Population Sizes using DNA Sequences From Multiple Loci.**
*Genetics*, 164:1645-1656.
Available at: [http://www.genetics.org/content/164/4/1645.long](http://www.genetics.org/content/164/4/1645.long)

* Yang Z., Rannala B. (2010)
**Bayesian species delimitation using multilocus sequence data.**
*Proceedings of the National Academy of Sciences*, 107(20):9264-9269.
doi:[10.1073/pnas.0913022107](http://dx.doi.org/10.1073/pnas.0913022107)

* Yang Z., Rannala B. (2014)
**Unguided species delimitation using DNA sequence data from multiple loci.**
*Molecular Biology and Evolution*, 31(12):3125-3135.
doi:[10.1093/molbev/msu279](http://dx.doi.org/10.1093/molbev/msu279)
