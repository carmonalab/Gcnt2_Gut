# GCNT2 Expression Profiles in Human and Mouse Gut Resident T Cells 🧬

## Description ✅ <a name="description"></a> 

This project relies on the Salas A scRNA-seq dataset (GEO Accession Code: GSE214695) from human colonic mucosa
to quantify _GCNT2_ expression across various cell types and distinct gastrointestinal conditions. 
_GCNT2_ expression is investigated in broad cell types, such as epithelial cells, B cells, mast cells, etc.
The study also seeks to determine whether gastrointestinal conditions such as healthy condition (HC), 
Crohn's disease (CD), or ulcerative colitis (UC) influence _GCNT2_ expression per cell type. 

## Table of Contents 📝 <a name="tof"></a>

-   [Description](#description)

-   [Table of Contents](#tof)

-   [R Environment](#renv)

-   [Processing](#processing)

-   [Clustering](#clustering)

-   [_GCNT2_ Expression](#gcnt2_expression)

-   [Author](#author)

## R Environment 🌲 <a name="renv"></a>

This project comes with its own R environment. In order to activate and be able to run all the comprised code, 
open the project through the **Gcnt2_Gut.Rproj** file or execute the following commands: 

```R
renv::activate()
renv::restore()
```
For ease, these commands have already been added to a code cell in **scRNAseq_data_processing_template.Rmd**, the
first file in the workflow. Hence, if you run this file first, as it is intended, there is no need to run these 
commands again elsewhere.

## Processing 📉 <a name="processing"></a>

A quality control pipeline was applied to the raw human colonic mucosa scRNA-seq dataset to single out healthy, 
normal cells, and discard extracellular droplets and vesicles or abnormal and dying cells. The table below shows 
the filtering indicators that were used and their cutoffs.

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

## Clustering 🧑‍🧑‍🧒‍🧒 <a name="clustering"></a>

Following the Seurat Guided Clustering Tutorial, the cells of the human colonic mucosa scRNA-seq dataset were 
manually annotated in the **cell_clustering.Rmd** file. The scGate package was also used to aid with the manual 
annotation. The clustering resolution was low and resulted in 9 distinct cell types, displayed in the table below.

| Cell Types | 
|------------|
| CD4 T Cells |
| CD8 T Cells |
| B Cells |
| Plasma Cells |
| Mast Cells |
| MoMac |
| Epithelial Cells |
| Endothelial Cells |
| Fibroblasts |

Given that the all cells were clustered and annotated at once, the **batch_effects_verification.Rmd** file was
created to check that no batch effects had been ignored throughout the process.

CD4 T cells and CD8 T cells were also subdivided into subtypes in the **t_cell_projection.Rmd** file, using
ProjecTILs. However, given that  cancerous reference maps were used, which do not apply in HC, CD, and UC
gastrointestinal conditions, the results of this step shall be ignored.

## _GCNT2_ Expression 🧪 <a name="gcnt2_expression"></a>

_GCNT2_ expression was quantified per cell type, per gut condition, and per cell type per gut condition.
The quantification was done in three distinct ways. Firstly, the ratios of _GCNT2_ expressing cells were investigated.
Secondly, _GCNT2_ expression was studied using the normalized counts. Thirdly, a pseudobulk analysis of _GCNT2_
expression was conducted. All results can be found in the **gcnt2_expression.Rmd** file.

### Deeper Dive Into B Cells 🧫

Differential expression of _GCNT2_ was observed within the sub-population of B cells. Hence, the 
**gcnt2_bcell_expression.Rmd** is dedicated to studying _GCNT2_ expression in this sub-population. In this file, 
the Seurat Guided Clustering Tutorial was executed again, solely on the B cells, in attempt to create clusters
that outline differential _GCNT2_ expression.

## Author ✍️ <a name="author"></a>

Aleksandar Mihaylov <aleksandar.mihaylov@epfl.ch>













