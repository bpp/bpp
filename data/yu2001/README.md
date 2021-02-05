This is the dataset from (Yu et al. 2001) consisting of a sample of 61 human sequences. When analyzing
this one species dataset, the MSC model becomes the single-population coalescent (Kingman, 1982). There
is no need for an Imap file, or the need to tag the sequence names in the sequence file. There is no
need for a species tree, so the block for specifying species names and species tree looks like this:

```
species&tree = 1 H
                 61
```

A tau prior is also not specified since there is no tau in this case.


## References

* Kingman JFC. (1982)
**The coalescent.**
*Stochastic Processes and their Applications*, 13(3):235-248.
doi:[10.1016/0304-4149(82)90011-4](https://doi.org/10.1016/0304-4149(82)90011-4)

* Yu N., Zhao Z., Fu YX, Sambuughin N., Ramsay M., Jenkins T., Leskinen E., Patthy L., Jorde LB, Kuromori T., Li WH. (2001)
**Global Patterns of Human DNA Sequence Variation in a 10-kb Region on Chromosome 1.**
*Molecular Biology and Evolution*, 18(2):214-222.
doi:[10.1093/oxfordjournals.molbev.a003795](https://doi.org/10.1093/oxfordjournals.molbev.a003795)
