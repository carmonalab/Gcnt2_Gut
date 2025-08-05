# GCNT2 Expression Profiles in Human and Mouse Gut Resident T Cells 🧬

## Description ✅ <a name="description"></a> 

This project relies on scRNA-seq  datasets from mouse and human guts to quantify _GCNT2_ expression across various cell types and distinct gastrointestinal conditions. 
The cellular focus of the analysis is CD8αα⁺ T cells (T⍺⍺) and CD8αβ⁺ tissue-resident memory T cells (Trm) although other gut resident cells are also investigated as a source of comparison.
Moreover, the study seeks to determine whether gastrointestinal conditions such as healthy, inflammatory bowel disease (IBD), or celiac disease influence _GCNT2_ expression per cell type. 
It is also interesting to compare how conserved _GCNT2_ expression profiles are between human and mouse guts per cell type. Ultimately, this project will suggest potential immmunoregulatory roles
for _GCNT2_ and hopefully elucidate how T⍺⍺ cells are activated.

## Table of Contents 📝 <a name="tof"></a>

-   [Description](#description)

-   [Table of Contents](#tof)

-   [R Environment](#renv)

-   [Processing](#processing)

## R Environment 🌲 <a name="renv"></a>

This project comes with its own R environment. In order to activate and be able to run all the comprised code, 
open the project through the **Gcnt2_Gut.Rproj** file and execute the following commands: 

```R
renv::activate()
renv::restore()
```
For ease, these commands have already been added to a code cell in **scRNAseq_data_processing_template.Rmd**, the
first file in the workflow. Hence, if you run this file first, as it is intended, there is no need to run these 
commands again elsewhere.

## Processing 📉 <a name="processing"></a>

A quality control procedure was applied to the human scRNA-seq dataset to single out healthy, normal cells and discard 
extracellular droplets and vesicles or abnormal and dying cells. The table below shows the filtering indicators 
that were used and their cutoffs.

| Indicator | Description | Cutoff |
|-----------|-------------|--------|
| Number of features | Quantifies the number of distinct genes identified within each scan | < 400 & ≥ 6000 |
| Number of counts | Quantifies the total number of mRNA molecules present in each scan | < 800 & ≥ 40000 |
| Ribosomal content | % of total mRNA molecules originating from ribosomes within each scan | < 0.5% & ≥ 60% |
| Mitochondrial content | % of total mRNA molecules originating from mitochondria within each scan | ≥ 20% |
| MALAT1 Expression | Number of MALAT1 mRNA reads within each scan | ≤ 2 |
| Log genes per UMI | log10(number of features) / log10(number of counts) | ≤ 0.6 |

The quality control pipeline is found in **scRNAseq_data_processing_template.Rmd**. The pipeline displays plots 
about various quality control metrics but also outputs the raw scRNA-seq dataset along with two processed scRNA-seq
datasets: the quality control end product and a lighter version.

For comparison statistics and graphs with the author's filtering procedure see **data_processing_comparison_author.Rmd** 

