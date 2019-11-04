This is the 106-locus data set from (Rokas et al. 2003) restricted to the five yeast species Scer, Spar, Smik, Skud and Sbay. 

| Species acronym | Full name                     |
| --------------- | ----------------------------- |
| **Scer**        | *Saccharomyces cerevisiae*    |
| **Spar**        | *Saccharomyces paradoxus*     |
| **Smik**        | *Saccharomyces mikatae*       |
| **Skud**        | *Saccharomyces kudriavzevii*  |
| **Sbay**        | *Saccharomyces bayanus*       |

Below is the species tree (network) for the yeast dataset with annotated introgressions

![yeast species network](https://raw.githubusercontent.com/xflouris/assets/master/bpp/yeast/yeast-5.png)

The newick string format specification for the network closely follows the notation from (Cardona et al. 2008):

```
((((Scer,Spar)A,Smik)B,(Skud,(Sbay)H[&phi=0.5,&tau-parent=no])D)C,H[&tau-parent=yes])R;
```

The annotations in the square brackets **[...]** are used to define the type of
hybridization (admixture) events. Each annotation is relative to the
corresponding parent of the node defined on. For example, the first annotation
(reading from left to right) is defined on inner node **H** and relates to the
parent node **D**. The attribute `tau-parent=no` dictates that parent node
**D** is not to have its own &tau; parameter, but instead, share the tau
parameter with node **H**.  The `phi` annotation indicates the value of `phi` in
case BPP is used for simulation purposes (`--simulate` switch), and is ignored
otherwise (although it must still be specified).

The second annotation is defined for node **H** and its parent node **R** (i.e,
the node after the closing parenthesis matching the opening parenthesis before
**H**). In this case, `tau-parent=yes` indicates that **R** is to have its own
&tau; parameter. 

For more information and a step-by-step guide on creating newick strings for
networks, please see section
[Specification of hybridization events with the MSci model](https://github.com/bpp/bpp/wiki/Specification-of-hybridization-events-with-the-MSci-model)
in the [BPP wiki](https://github.com/bpp/bpp/wiki).

### Running the dataset

Download and install the BPP version suitable for your operating system (see [Download and install](https://github.com/bpp/bpp#download-and-install)).

Run BPP using the control file

```bash
bpp --cfile Rokas2003-5species-bpp.ctl
```

This may take some minutes, depending on the processor.  You may also
experiment with the number of threads used to see which value gives you the
shortest running time.

The results should be comparable with Figure 22e in (Wen and Nakhleh, 2018).
For example, the posterior mean and median of the introgression probability
should be 0.70 with the 95% HPD interval to be (0.56,0.83).






## References

* Rokas A., Williams BL, King N., Carroll SB. (2003)
**Genome-scale approaches to resolving incongruence in molecular phylogenies.**
*Nature*, 425:798-804.
doi:[10.1038/nature020503](http://dx.doi.org/10.1038/nature020503)

* Cardona G., Rossello F., Valiente G. (2008)
**Extended Newick: it is time for a standard representation of phylogenetic networks.**
*BMC Bioinformatics*, 9:532.
doi:[10.1186/1471-2105-9-532](https://doi.org/10.1186/1471-2105-9-532)

* Wen D., Nakhleh L. (2018) 
**Coestimating reticulate phylogenies and gene trees from multilocus sequence data.**
*Systematic Biology*, 67(3):439-457.
doi:[10.1093/sysbio/syx085](http://dx.doi.org/10.1093/sysbio/syx085)
