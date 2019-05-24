# bpp

[![Build Status](https://travis-ci.org/bpp/bpp.svg?branch=master)](https://travis-ci.org/bpp/bpp)
[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/DOI/10.1093/molbev/msy147.svg)](https://doi.org/10.1093/molbev/msy147)
[![Version](https://img.shields.io/badge/version-4.1.4-blue.svg)](https://github.com/bpp/bpp/releases/tag/v4.1.4)

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

BPP now also implements the Multispecies-coalescent-with-introgression (MSci) model,
an extension of the multispecies coalescent model to incorporate introgression/hybridization.
For more information on usage please see the [wiki documentation](https://github.com/bpp/bpp/wiki).

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
wget https://github.com/bpp/bpp/releases/download/v4.1.4/bpp-4.1.4-linux-x86_64.tar.gz
tar zxvf bpp-4.1.4-linux-x86_64.tar.gz
```

Or these commands if you using a Mac:

```bash
wget https://github.com/bpp/bpp/releases/download/v4.1.4/bpp-4.1.4-macos-x86_64.tar.gz
tar zxvf bpp-4.1.4-macos-x86_64.tar.gz
```

Or if you are using Windows, download and extract (unzip) the contents of this file:

```
https://github.com/bpp/bpp/releases/download/v4.1.4/bpp-4.1.4-win-x86_64.zip
```

Linux and Mac: You will now have the binary distribution in a folder called `bpp-4.1.4-linux-x86_64` or `bpp-4.1.4-macos-x86_64`. The binary file is located in the `bin` subfolder, i.e. `bin/bpp`. We recommend making a copy or a symbolic link to the binary in a folder included in your `$PATH`. 

Windows: You will now have the binary distribution in a folder called `bpp-4.1.4-win-x86_64`. The bpp executable is called `bpp.exe`.


**Compiling from source** You can either download the *source distribution* for a particular version or *clone the repository* (see below). In both cases, you will need several packages installed on your system. Currently, BPP requires that [GNU Bison](http://www.gnu.org/software/bison/)
and [Flex](http://flex.sourceforge.net/) are installed on the target system.
On a Debian-based Linux system, the two packages can be installed using the
command

```bash
apt-get install flex bison
```

**Source distribution** To download the source distribution from a [release](https://github.com/bpp/bpp/releases) and build the executable and documentation, use the following commands:

```
wget https://github.com/bpp/bpp/archive/v4.1.4.tar.gz
tar zxvf v4.1.4.tar.gz
cd bpp-4.1.4/src
make
```

**Cloning the repo** Instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below.

```
git clone https://github.com/bpp/bpp.git
cd bpp
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

More documentation regarding control files, will be available soon on the [wiki](https://github.com/bpp/bpp/wiki).

## Citing BPP

Please cite the following publication if you use BPP:

Flouri T., Jiao X., Rannala B., Yang Z. (2018) Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent. Molecular Biology and Evolution (accepted manuscript).
doi:[10.1093/molbev/msy147](https://doi.org/10.1093/molbev/msy147)

Please note, citing the corresponding of the four underlying methods, may also be appropriate.

## License and third party licenses

The code is currently licensed under the [GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).

## Code

| File                       | Description                                                                       |
| -------------------------- | --------------------------------------------------------------------------------- |
| **arch.c**                 | Architecture specific code (Linux/Mac/Windows)                                    |
| **allfixed.c**             | Summary statistics for method A00 (fixed species tree)                            |
| **bpp.c**                  | Main file handling command-line parameters and executing selected methods         |
| **cfile.c**                | Functions for parsing the control file                                            |
| **compress.c**             | Functions for compressing multiple sequence alignments into site patterns         |
| **core_likelihood.c**      | Core functions for evaluating the likelihood of a tree (non-vectorized)           |
| **core_likelihood_avx.c**  | Core functions for evaluating the likelihood of a tree (AVX version)              |
| **core_likelihood_avx2.c** | Core functions for evaluating the likelihood of a tree (AVX-2 version)            |
| **core_likelihood_sse.c**  | Core functions for evaluating the likelihood of a tree (SSE-3 version)            |
| **core_partials.c**        | Core functions for computing partial likelihoods (non-vectorized)                 |
| **core_partials_avx.c**    | Core functions for computing partial likelihoods (AVX version)                    |
| **core_partials_avx2.c**   | Core functions for computing partial likelihoods (AVX-2 version)                  |
| **core_partials_sse.c**    | Core functions for computing partial likelihoods (SSE-3 version)                  |
| **core_pmatrix.c**         | Core functions for constructing the transition probability matrix                 |
| **delimit.c**              | Species delimitation auxiliary functions and summary statistics                   |
| **diploid.c**              | Functions for resolving/phasing diploid sequences                                 |
| **dlist.c**                | Functions for handling doubly linked-lists                                        |
| **dump.c**                 | Functions for dumping the MCMC state into a checkpoint file                       |
| **experimental.c**         | Experimental functions that are not yet production-ready                          |
| **gtree.c**                | Functions for setting and processing gene trees                                   |
| **hardware.c**             | Functions for hardware detection                                                  |
| **hash.c**                 | Hash table implementation and related functions                                   |
| **lex_map.l**              | Lexical analyzer for parsing map files                                            |
| **lex_rtree.l**            | Lexical analyzer for parsing newick rooted trees                                  |
| **list.c**                 | Linked list implementation and related functions                                  |
| **load.c**                 | Functions for loading a checkpoint file                                           |
| **locus.c**                | Locus specific functions                                                          |
| **Makefile**               | Makefile                                                                          |
| **mapping.c**              | Functions for handling map files                                                  |
| **maps.c**                 | Character mapping arrays for converting sequences to the internal representation  |
| **method.c**               | Function containing the MCMC loop and calls to proposals                          |
| **msa.c**                  | Code for processing multiple sequence alignments                                  |
| **output.c**               | Auxiliary functions for printing pmatrices (to-be-renamed)                        |
| **parse_map.y**            | Functions for parsing map files                                                   |
| **parse_rtree.y**          | Functions for parsing rooted trees in newick format                               |
| **phylip.c**               | Functions for parsing phylip files                                                |
| **random.c**               | Pseudo-random number generator functions                                          |
| **rtree.c**                | Species tree export functions (to-be-renamed).                                    |
| **stree.c**                | Functions for setting and processing the species tree                             |
| **summary.c**              | Species tree inference summary related functions                                  | 
| **util.c**                 | Various common utility functions                                                  |

# Acknowledgements

Special thanks to Yuttapong Thawornwattana and [Mario dos Reis Barros](http://www.sbcs.qmul.ac.uk/staff/mariodosreisbarros.html) for testing and bug reports.

# References

* Flouri T., Carrasco FI, Darriba D., Aberer AJ, Nguyen LT, Minh BQ, Haeseler A., Stamatakis A. (2015)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2):356-362.
doi:[10.1093/sysbio/syu084](http://dx.doi.org/10.1093/sysbio/syu084)

* Flouri T., Jiao X., Rannala B., Yang Z. (2018)
**Species Tree Inference with BPP using Genomic Sequences and the Multispecies Coalescent.**
*Molecular Biology and Evolution*, Accepted Manuscript.
doi:[10.1093/molbev/msy147](http://dx.doi.org/10.1093/molbev/msy147)

* Yang Z., Rannala B. (2003)
**Bayes Estimation of Species Divergence Times and Ancestral Population Sizes using DNA Sequences From Multiple Loci.**
*Genetics*, 164:1645-1656.
Available at: [http://www.genetics.org/content/164/4/1645.long](http://www.genetics.org/content/164/4/1645.long)

* Hey J., Nielsen R. (2004)
**Multilocus methods for estimating population sizes, migration rates and divergence time, with applications to the divergence of Drosophila pseudoobscura and D. persimilis.**
*Genetics*, 167(2):747-760.
doi:[10.1534/genetics.103.024182](http://dx.doi.org/10.1534/genetics.103.024182)

* Burgess R., Yang Z. (2008)
**Estimation of hominoid ancestral population sizes under Bayesian coalescent models incorporating mutation rate variation and sequencing errors.**
*Molecular Biology and Evolution*, 25(9):1979-1994.
doi:[10.1093/molbev/msn148](http://dx.doi.org/10.1093/molbev/msn148)

* Yang Z., Rannala B. (2010)
**Bayesian species delimitation using multilocus sequence data.**
*Proceedings of the National Academy of Sciences*, 107(20):9264-9269.
doi:[10.1073/pnas.0913022107](http://dx.doi.org/10.1073/pnas.0913022107)

* Gronau I., Hubisz MJ, Gulko B., Danko CG, Siepel A. (2011)
**Bayesian inference of ancient human demography from individual genome sequences.**
*Nature Genetics*, 43(10):1031-1035.
doi:[10.1038/ng.937](http://dx.doi.org/10.1038/ng.937)

* Rannala B., Yang Z. (2013)
**Improved reversible jump algorithms for Bayesian species delimitation.**
*Genetics*, 194:245-253.
doi:[10.1534/genetics.112.149039](http://dx.doi.org/10.1534/genetics.112.149039)

* Yang Z., Rannala B. (2014)
**Unguided species delimitation using DNA sequence data from multiple loci.**
*Molecular Biology and Evolution*, 31(12):3125-3135.
doi:[10.1093/molbev/msu279](http://dx.doi.org/10.1093/molbev/msu279)

* Rannala B., Yang Z. (2017)
**Efficient Bayesian Species Tree Inference under the Multispecies Coalescent.**
*Systematic Biology*, 66(5):823-842.
doi:[0.1093/sysbio/syw119](http://dx.doi.org/0.1093/sysbio/syw119)
