This dataset represents the third block of 100 loci for the non-coding 2L1
region as used in (Flouri et al. 2019).

Below is the species tree (network) for the anopheles dataset with annotated introgressions 

![anopheles species network](https://raw.githubusercontent.com/xflouris/assets/master/bpp/anopheles/anopheles.png)

The newick string format specification for the network closely follows the notation from (Cardona et al. 2008):

```
 ((R,(Q)h[&phi=0.3,&tau-parent=no])g,(f[&tau-parent=yes,&phi=0.3],(((((G,C)b)f[&tau-parent=no],A)e,h[&tau-parent=yes])d,L)c)a)o;
```

The annotations in the square brackets **[...]** are used to define the type of
hybridization (admixture) events. Each annotation is relative to the
corresponding parent of the node defined on. For example, the first annotation
(reading from left to right) is defined on inner node **h** and relates to the
parent node **g**. The attribute `tau-parent=no` dictates that parent node
**g** is not to have its own &tau; parameter, but instead, share the tau
parameter with node 'h'.  The `phi` annotation indicates the value of `phi` in
case BPP is used for simulation purposes (`--simulate` switch), and is ignored
otherwise (although it must still be specified).

The next annotation is defined for node **f** and its parent node **a** (i.e,
the node after the closing parenthesis matching the opening parenthesis before
**f**). In this case, `tau-parent=yes` indicates that **a** is to have its own
&tau; parameter. Finally, the last two annotations are for nodes **f** (with
respect to parent node **e**) and **h** (with respect to node **d**).

For more information and a step-by-step guide on creating newick strings for
networks, please see section
[Specification of hybridization events with the MSci model](https://github.com/bpp/bpp/wiki/Specification-of-hybridization-events-with-the-MSci-model)
in the [BPP wiki](https://github.com/bpp/bpp/wiki).

### Running the dataset

Download and install the BPP version suitable for your operating system (see
[Download and install](https://github.com/bpp/bpp#download-and-install)).

Run BPP using the control file

```bash
bpp --cfile bpp.ctl
```
This assumes the JC mutation model and the molecular clock.  With the settings
`burnin = 32000`, `sampfreq = 2`, `nsample = 500000`, the run took about 1h10m
using 4 threads on an
[Intel Xeon Gold 6154 CPU](https://ark.intel.com/content/www/us/en/ark/products/120495/intel-xeon-gold-6154-processor-24-75m-cache-3-00-ghz.html)
CPU. The posterior means and 95% HPD CI are **0.176 (0.069, 0.288)** for
`phi_h` (&phi;<sub>R&rarr;Q</sub>), and **0.981 (0.946, 1.000)** for
`phi_f` (&phi;<sub>A&rarr;GC</sub>).  Estimates of other parameters
are in the excel file anopheles-noncoding-block3-MSci-estimates.xlsx.


For illustration, the control file bpp.gtr.ctl specifies the GTR mutation
model.

```bash
bpp --cfile bpp.gtr.ctl
```

With all other settings the same as for the JC model, the run took about 2h10m
using 4 threads on the same CPU.  The posterior means and 95% HPD CI are
**0.173 (0.069, 0.288)** for `phi_h` (&phi;<sub>R&rarr;Q</sub>), and 
**0.982 (0.949, 1.000)** for `phi_f` (&phi;<sub>A&rarr;GC</sub>).  The results
are virtually identical to those under JC.  Those species are very closely
related and the sequences are highly similar, so that JC does a good job at
multiple-hit correction, and the use of GTR is unnecessary.

### Estimates

#### Theta estimates

<table>
<tr><th>JC theta estimates</th><th>GTR theta estimates</th></tr>
<tr><td>

|   &theta;   | mean     | 2.5% HPD | 97.5% HPD |
|-------------|----------|----------|-----------|
| theta_1G    | 0.020353 | 0.009347 | 0.034987  |
| theta_2C    | 0.047352 | 0.013925 | 0.098182  |
| theta_3R    | 0.00353  | 0.002627 | 0.004482  |
| theta_4L    | 0.002372 | 0.001745 | 0.003022  |
| theta_5A    | 0.004086 | 0.003091 | 0.005206  |
| theta_6Q    | 0.007419 | 0.005723 | 0.009282  |
| theta_7o    | 0.009217 | 0.005397 | 0.01396   |
| theta_8g    | 0.036599 | 0.007    | 0.081634  |
| theta_10a   | 0.018444 | 0.003386 | 0.046029  |
| theta_11c   | 0.00868  | 0.004698 | 0.01328   |
| theta_12d   | 0.00516  | 0.002772 | 0.007984  |
| theta_13e   | 0.008921 | 0.006147 | 0.011845  |
| theta_15b   | 0.004749 | 0.003339 | 0.006245  |
| theta_16h   | 0.020236 | 0.003429 | 0.049355  |
| theta_17f   | 0.019664 | 0.003521 | 0.048952  |

</td><td>

|   &theta;   | mean      |  2.5% HPD | 97.5% HPD  |
|-------------|-----------|-----------|------------|
| theta_1G    | 0.02006   | 0.00941   | 0.0343645  |
| theta_2C    | 0.048687  | 0.0142785 | 0.102599   |
| theta_3R    | 0.0035275 | 0.0026305 | 0.0044685  |
| theta_4L    | 0.002368  | 0.0017355 | 0.0030165  |
| theta_5A    | 0.004078  | 0.0030815 | 0.005151   |
| theta_6Q    | 0.0075725 | 0.0058235 | 0.009477   |
| theta_7o    | 0.0089225 | 0.005873  | 0.0122755  |
| theta_8g    | 0.037665  | 0.007166  | 0.082657   |
| theta_10a   | 0.0188865 | 0.0034215 | 0.047591   |
| theta_11c   | 0.0083765 | 0.00481   | 0.0125065  |
| theta_12d   | 0.0051465 | 0.0028495 | 0.0079185  |
| theta_13e   | 0.0087985 | 0.0060765 | 0.0116845  |
| theta_15b   | 0.004638  | 0.0032715 | 0.006078   |
| theta_16h   | 0.0198265 | 0.0035675 | 0.0491005  |
| theta_17f   | 0.020095  | 0.00343   | 0.049422   |

</td></tr> </table>

## References

* Flouri T., Jiao X., Rannala B., Yang Z. (2019)
**A Bayesian multispecies coalescent model with introgression for comparative genomic analysis.**
*Under review*

* Cardona G., Rossello F., Valiente G. (2008)
**Extended Newick: it is time for a standard representation of phylogenetic networks.**
*BMC Bioinformatics*, 9:532.
doi:[10.1186/1471-2105-9-532](https://doi.org/10.1186/1471-2105-9-532)
