# bpp

[![Build Status](https://travis-ci.org/xflouris/bpp.svg?branch=master)](https://travis-ci.org/xflouris/bpp)
[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

The aim of this project is to implement a versatile high-performance version
of the BPP software. It should have the following properties:

* open-source code with an appropriate open-source license.
* 64-bit multi-threaded design that handles very large datasets.
* easy to use and well-documented.
* SIMD implementations of time-consuming parts.
* Linux, Mac and Microsoft Windows compatibility.

## Compilation instructions

Currently, BPP requires that [GNU Bison](http://www.gnu.org/software/bison/)
and [Flex](http://flex.sourceforge.net/) are installed on the target system.
On a Debian-based Linux system, the two packages can be installed using the
command

```bash
apt-get install flex bison
```

BPP can then be compiled using the provided `Makefile`

```bash
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


## License and third party licenses

The code is currently licensed under the [GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).

## Code

| File                  | Description                                                                       |
| --------------------- | --------------------------------------------------------------------------------- |
| **arch.c**            | Architecture specific code (Mac/Linux).                                           |
| **bpp.c**             | Main file handling command-line parameters and executing selected methods.        |
| **compress.c**        | Functions for compressing multiple sequence alignments into site patterns         |
| **gtree.c**           | Functions for setting and processing gene trees                                   |
| **hash.c**            | Hash table implementation and related functions                                   |
| **lex_map.l**         | Lexical analyzer for parsing map files.                                           |
| **lex_rtree.l**       | Lexical analyzer for parsing newick rooted trees.                                 |
| **likelihood_avx.c**  | AVX likelihood functions.                                                         |
| **likelihood_sse.c**  | SSE likelihood functions.                                                         |
| **likelihood.c**      | Likelihood related functions.                                                     |
| **list.c**            | Linked list implementation and related functions                                  |
| **locus.c**           | Locus specific functions.                                                         |
| **Makefile**          | Makefile                                                                          |
| **mapping.c**         | Functions for handling map files.                                                 |
| **maps.c**            | Character mapping arrays for converting sequences to the internal representation. |
| **msa.c**             | Code for processing multiple sequence alignments                                  |
| **parse_map.y**       | Functions for parsing map files.                                                  |
| **parse_rtree.y**     | Functions for parsing rooted trees in newick format.                              |
| **phylip.c**          | Functions for parsing phylip files.                                               |
| **random.c**          | Pseudo-random number generator functions                                          |
| **rtree.c**           | Rooted tree manipulation functions.                                               |
| **stree.c**           | Functions for setting and processing the species tree                             |
| **util.c**            | Various common utility functions.                                                 |

# Acknowledgements

Special thanks to Yuttapong Thawornwattana for testing and bug reports.

# References

* Flouri T., Carrasco FI, Darriba D., Aberer AJ, Nguyen LT, Minh BQ, Haeseler A., Stamatakis A. (2015)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2):356-362.
doi:[10.1093/sysbio/syu084](10.1093/sysbio/syu084)

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
doi:[10.1093/molbev/msu279](10.1093/molbev/msu279)
