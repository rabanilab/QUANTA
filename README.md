# QUANTA

QUANTA is a software package to distinguish transcriptionally silent genes and analyze their post-transcriptional regulation.
QUANTA uses kinetic modeling to compare total and polyA+ expression patterns of genes,
and quantitatively dissect the interconnected kinetics of their mRNA polyadenylation and degradation.

You can find the full description in:

[A comparative analysis across species of maternal mRNA regulation in embryos](https://www.ncbi.nlm.nih.gov/pubmed/29727622).
Tawil M, Alcalay D, Greenberg P, Har-Sheffer S, Fishman L, Rabani M (2024). Submitted.

## Installation

These instructions will get you a copy of the QUANTA analysis software.

The first part of QUANTA analyzes FASTQ files to get FPKM counts per gene,
and is implemented in linux.

The second part of QUANTA performs kinetic analysis,
and runs with Matlab on your local machine.

### Prerequisites

We developed and tested QUANTA with Matlab R2023b. Matlab can be obtained and
installed from [Mathworks](https://www.mathworks.com/products/matlab.html).

### Installing

Download the package source code from GitHub.



## Example data for FASTQ file analysis

TBD


## Example data for kinetic analysis

The example data directory (data_files) contains FPKM counts from multiple
temporal zebrafish datasets that were already downloaded and quantified.

Datasets:

```
Medina et. al., 2021 (polyA and total-RNA)
Meier et. al., 2018 (total-RNA)
Zhao et. al., 2017 (total-RNA)
Pauli et. al., 2011 (polyA-RNA)
Harvey et. al., 2013 (polyA-RNA)
Yang et. al., 2019 (polyA-RNA)
```

To run the example, simply type within matlab:

```
addpath <path to installation directory>;
analyze;
```

This will run the kinetics analysis on the zebrafish data in the example directory.
The results of this analysis will be available in a "results" directory.


## License

This project is licensed under the MIT License - see source files for details.

