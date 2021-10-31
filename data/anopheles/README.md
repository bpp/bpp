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
| theta_1G    | 0.020797 | 0.008611 |  0.03795  |
| theta_2C    | 0.046268 | 0.012601 |  0.102933 |
| theta_3R    | 0.003237 | 0.002439 |  0.004095 |
| theta_4L    | 0.002276 | 0.001677 |  0.002904 |
| theta_5A    | 0.003871 | 0.002919 |  0.004911 |
| theta_6Q    | 0.006419 | 0.00495  |  0.00802  |
| theta_7o    | 0.016853 | 0.013098 |  0.020786 |
| theta_8g    | 0.054167 | 0.006109 |  0.130988 |
| theta_10a   | 0.01379  | 0.003627 |  0.03025  |
| theta_11c   | 0.016743 | 0.003748 |  0.03896  |
| theta_12d   | 0.007415 | 0.00428  |  0.01096  |
| theta_13e   | 0.007023 | 0.004091 |  0.01032  |
| theta_15b   | 0.005967 | 0.004181 |  0.007883 |
| theta_16h   | 0.018344 | 0.003228 |  0.045553 |
| theta_17f   | 0.019365 | 0.00361  |  0.047733 |

</td><td>

|   &theta;   | mean     | 2.5% HPD | 97.5% HPD |
|-------------|----------|----------|-----------|
| theta_1G    | 0.020498 | 0.008627 |  0.037276 |
| theta_2C    | 0.047439 | 0.012648 |  0.106778 |
| theta_3R    | 0.003222 | 0.002418 |  0.004063 |
| theta_4L    | 0.002269 | 0.00167  |  0.002881 |
| theta_5A    | 0.00386  | 0.002903 |  0.004881 |
| theta_6Q    | 0.006467 | 0.004962 |  0.008069 |
| theta_7o    | 0.016882 | 0.013169 |  0.020831 |
| theta_8g    | 0.058729 | 0.005864 |  0.1408   |
| theta_10a   | 0.013847 | 0.003712 |  0.029899 |
| theta_11c   | 0.016571 | 0.003772 |  0.038701 |
| theta_12d   | 0.007343 | 0.004283 |  0.01083  |
| theta_13e   | 0.006861 | 0.003983 |  0.010142 |
| theta_15b   | 0.005909 | 0.004122 |  0.007801 |
| theta_16h   | 0.0187   | 0.003194 |  0.046225 |
| theta_17f   | 0.020447 | 0.003473 |  0.049717 |

</td></tr> </table>

#### Tau estimates

<table>
<tr><th>JC tau estimates</th><th>GTR tau estimates</th></tr>
<tr><td>

|  &tau;  | mean     | 2.5% HPD | 97.5% HPD |
|---------|----------|----------|-----------|
| tau_7o  | 0.0145   | 0.013218 |  0.015808 |
| tau_8g  | 0.007343 | 0.006155 |  0.00849  |
| tau_9h  | 0.007343 | 0.006155 |  0.00849  |
| tau_10a | 0.013294 | 0.011698 |  0.014927 |
| tau_11c | 0.01231  | 0.011029 |  0.013554 |
| tau_12d | 0.007896 | 0.006978 |  0.008804 |
| tau_13e | 0.005295 | 0.004574 |  0.00601  |
| tau_14f | 0.005295 | 0.004574 |  0.00601  |
| tau_15b | 0.001976 | 0.001494 |  0.00246  |

</td><td>

|  &tau;  | mean     | 2.5% HPD | 97.5% HPD |
|---------|----------|----------|-----------|
| tau_7o, | 0.014481 | 0.013191 |  0.015745 |
| tau_8g, | 0.007274 | 0.006074 |  0.008431 |
| tau_9h, | 0.007274 | 0.006074 |  0.008431 |
| tau_10a | 0.013261 | 0.011655 |  0.014908 |
| tau_11c | 0.012276 | 0.01101  |  0.01353  |
| tau_12d | 0.007852 | 0.006953 |  0.008759 |
| tau_13e | 0.005333 | 0.00461  |  0.006062 |
| tau_14f | 0.005333 | 0.00461  |  0.006062 |
| tau_15b | 0.001992 | 0.001503 |  0.002466 |

</td></tr> </table>

#### Phi estimates

<table>
<tr><th>JC phi estimates</th><th>GTR phi estimates</th></tr>
<tr><td>

|  &phi;  | mean     | 2.5% HPD | 97.5% HPD |
|---------|----------|----------|-----------|
|  phi_h  | 0.245163 | 0.064886 |  0.415232 |
|  phi_f  | 0.974894 | 0.929059 |  0.999989 |

</td><td>

|  &phi;  | mean     | 2.5% HPD | 97.5% HPD |
|---------|----------|----------|-----------|
|  phi_h  | 0.243097 | 0.064691 | 0.411921  |
|  phi_f  | 0.975208 | 0.929233 | 1         |

</td></tr> </table>



### Drawing the species tree with introgressions using FigTree

At the end of the MCMC run under the MSci model, BPP creates a file named
`FakeTree.tre`, which has branch lengths given by the posterior means of node ages (&tau;) and the 95% HPD CI 
as node bars.  It also includes the estimated population size parameters (&theta;) as branch/node labels. 

This tree is binary with "hybrid mirror nodes" converted into ghost tips, so the tree can be read into [FigTree](http://tree.bio.ed.ac.uk/software/figtree/).
You can  save the tree into a vector graphics file (instead of a pictcure file), and import it into a graphics program for editing. 
On Windows .wmf (Windows Metafile) and .emf (Enhanced Windows Metafile) formats seem to work well, so the procedure will be as follows.

First, open FakeTree.tre in FigTree.  Choose node bars, node labels, change curvature etc. as you like.  
Save the file as a graphics file (File - Export Graphics or Ctrl-Alt-E, choose file type .emf).

Second start MS Powerpoint, Insert-Picture from file, Ungroup (Ctrl-Shift-G), twice.  Texts and lines will be editable.

On a MAC, you can try to save the graphics into a sgv file.


Below is an example of the generated `FakeTree.tre`.

![anopheles species network](https://raw.githubusercontent.com/xflouris/assets/master/bpp/anopheles/anopheles-faketree.png)

Example `FakeTree.tre` file:

```

#NEXUS

BEGIN TREES;
        UTREE 1 = ((R: 0.006978, (Q: 0.006978, ghost_h :0.006978)h[&height_95%_HPD={0.00002400, 0.01377600}, theta=-1.0000000]: 0.000000)g[&height_95%_HPD={0.00002400, 0.01377600}, theta=0.0189369]: 0.014756, (f :0.018738, (((((G: 0.001570, C: 0.001570)b[&height_95%_HPD={0.00001200, 0.00376500}, theta=0.0147427]: 0.003392, ghost_f :0.004962)f[&height_95%_HPD={0.00108900, 0.00882000}, theta=-1.0000000]: 0.000000, A: 0.004962)e[&height_95%_HPD={0.00108900, 0.00882000}, theta=0.0108381]: 0.007289, h :0.012250)d[&height_95%_HPD={0.00580100, 0.01864700}, theta=0.0186650]: 0.003424, L: 0.015674)c[&height_95%_HPD={0.00920400, 0.02270300}, theta=0.0190814]: 0.003063)a[&height_95%_HPD={0.01148500, 0.02611200}, theta=0.0185349]: 0.002996)o[&height_95%_HPD={0.01503800, 0.02909800}, theta=0.0118272];
END;

```

## References

* Flouri T., Jiao X., Rannala B., Yang Z. (2020)
**A Bayesian Implementation of the Multispecies Coalescent Model with Introgression for Phylogenomic Analysis.**
*Molecular Biology and Evolution*, 37(4):1211-1223.
doi:[10.1093/molbev/msz296](https://doi.org/10.1093/molbev/msz296)

* Cardona G., Rossello F., Valiente G. (2008)
**Extended Newick: it is time for a standard representation of phylogenetic networks.**
*BMC Bioinformatics*, 9:532.
doi:[10.1186/1471-2105-9-532](https://doi.org/10.1186/1471-2105-9-532)
