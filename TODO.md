TODOS:

* gwasim simulations still scale so-so, improve via paralellization or faster(?) large matrix operations,
* write a function where effect size is independent of maf,
* discuss whichh architectures do we want to study, how to specify allele effects (per gt. distributions, what with additivity?),
* do we want direct support for tools other than SKAT & GenABEL?
* what input formats do we need? Currently: VCF, GenABEL:gwaa_data-class objects and other, e.g. ped or Affymetrix but via GenABEL which is obsolete,
* we can save as tibble (although conversion of large datasets may take time), do we want some other formats,
* perhaps running more tests by others
