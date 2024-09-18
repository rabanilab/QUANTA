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

To run the pre-processing you will need to make sure two external tools are installed on your machine:
1. samtools (can be downloaded and installed from [https://www.htslib.org/download](https://www.htslib.org/download/).
2. cufflinks (can be downloaded and installed from [https://cole-trapnell-lab.github.io/cufflinks](https://cole-trapnell-lab.github.io/cufflinks).

Add a bam file to the installation directory in the "bam_files" directory.

After installation, you can test the run by typing:
```
cd bam_analysis;
make test_run;
```

This will run the FPKM quantification distinguishing pre-mRNA and mRNA on zebrafish data in the bam_flies directory.



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

