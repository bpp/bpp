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

Download and install the BPP version suitable for your operating system (see [Download and install](https://github.com/bpp/bpp#download-and-install)).

Run BPP using the control file

```bash
bpp --cfile bpp.ctl
```

## References

* Flouri T., Jiao X., Rannala B., Yang Z. (2019)
**A Bayesian multispecies coalescent model with introgression for comparative genomic analysis.**
*Under review*

* Cardona G., Rossello F., Valiente G. (2008)
**Extended Newick: it is time for a standard representation of phylogenetic networks.**
*BMC Bioinformatics*, 9:532.
doi:[10.1186/1471-2105-9-532](https://doi.org/10.1186/1471-2105-9-532)
