# QUANTA

QUANTA is a software package to distinguish transcriptionally silent genes and analyze their post-transcriptional regulation.
QUANTA uses kinetic modeling to compare total and polyA+ expression patterns of genes,
and quantitatively dissect the interconnected kinetics of their mRNA polyadenylation and degradation.

You can find the full description in:

[A comparative analysis across species of maternal mRNA regulation in embryos using QUANTA](https://www.ncbi.nlm.nih.gov/).
Tawil M, Alcalay D, Greenberg P, Har-Sheffer S, Fishman L, Rabani M (2024). Submitted.


## Package Description and Prerequisites

The first part of QUANTA analyzes RNA-Seq data to get FPKM counts per gene,
and is implemented in linux.

The second part of QUANTA performs kinetic analysis,
and runs with Matlab on your local machine.

We developed and tested QUANTA with Matlab R2023b. Matlab can be obtained and
installed from [Mathworks](https://www.mathworks.com/products/matlab.html).


## Example data for RNA-Seq analysis

TBD


## Example data for kinetic analysis

The example data directory (data_files) contains FPKM counts from two temporal zebrafish datasets 
(Medina et. al., 2021) that were already downloaded and quantified.

To run the example, simply type within matlab:

```
addpath <path to installation directory>;
analyze;
```

This will run the kinetics analysis on the zebrafish data in the example directory.
The results of this analysis will be available in a "results" directory.


## License

This project is licensed under the MIT License - see source files for details.

