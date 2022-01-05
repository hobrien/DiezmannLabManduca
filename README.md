# Tobacco Hornworm as a novel host model for the study of fungal virulence and drug efficacy

- Scripts and Workflows used for the analyses in [Lyons N, I Softley, A Balfour, C Williamson, HE O'Brien, AC Shetty, VM Bruno and S Diezmann. 2020. Tobacco Hornworm (_Manduca sexta_) caterpillars as a novel host model for the study of fungal virulence and drug efficacy. Virulence 11(1): 1075-1089](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7549948/)

- All R code can be found [here](https://github.com/hobrien/DiezmannLabManduca/blob/master/R/Results.Rmd)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to run RNAseq analyses and R code can be found [here](https://github.com/hobrien/DiezmannLabManduca/blob/master/Snakefile)
- RNAseq analysis uses [HISAT2](http://ccb.jhu.edu/software/hisat2/manual.shtml) for mapping, [htseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html) for quantification, and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for Differential Expression analysis
    - Differential Expression analysis uses a [modified](https://github.com/hobrien/SARTools) version of the [SARTools](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0157022) pipeline
- Survival analysis uses the [Survival](https://github.com/therneau/survival) and survival plots were made using [survminer](https://github.com/kassambara/survminer)
