<!-- README.md is generated from README.Rmd. Please edit that file -->



# n1pas

Welcome to the <span style="color:purple">N-of-1-*pathways*</span> Alternatively Spliced (N1PAS) R package github repo! Below is a short overview of the method, installation instructions, and basic usage. Each function included in the package has complete documentation for more details.

## Overview

The software transforms a pair of isoform-level RNA-seq data sets from a individual subject to a personal profile of pathway enrichment of alternatively spliced genes. The schematic below illustrates the N1PAS workflow:

![plot of chunk workflow_jpg](vignettes/workflow_diagram.jpg)

**Workflow of N1PAS**. A) Isoform-specific mRNA-Seq data are obtained from an individual. Gene-level distances between the two samples indicates the magnitude of alternative splicing. B) Gene-level distances are aggregated across the whole transcriptome and unsupervised clustering classifies genes as alternatively spliced genes (ASGs; blue) and not (red). The vertical axis shows the count of genes with Hellinger distances binned together to form the histogram. Note that since this is a univariate setting, 2-means simply finds a threshold Hellinger distance to classify the genes into two groups. Gene set (pathway) enrichment analysis is conducted by first C) computing the odds ratio (OR) to quantify the relative abundance of ASGs in the pathway versus genes not in the pathway (background). ORs are calculated for each pathway in the data-base to produce an empirical, subject-specific distribution (D). A local false discovery procedure provides uncertainty quantification (FDR) and classified pathways as enriched using a simple threshold at FDR < 20%. E) Results are tabulated to provide an individualized profile of alternative spliced enrichment.

## Installation


```r
## Install from github directly using the R package 'devtools'
# install.packages("devtools")
devtools::install_github("grizant/n1pas")
```

## Usage

Included with the package is one TCGA BLCA (Bladder cancer) patient's paired tumor-normal isoform RSEM counts.


```r
data(blca_patient_iso_kegg)
## Display a few genes for the first two pairs of transcriptomes.
head(blca_patient_iso_kegg)
#>            geneSymbol TCGA-BL-A13J-N TCGA-BL-A13J-T
#> uc002bgz.2      OR4F5         8.3372        15.9111
#> uc003wfr.3     OR4F16        60.7593        76.1321
#> uc011cbi.1      ISG15         0.0000         0.0000
#> uc001jfy.3       AGRN         0.0000         0.0000
#> uc001jji.2    B3GALT6         0.0000         0.0000
#> uc002cyw.2     TAS1R3         0.0000         0.0000
## Note that the rows are labeled by isoform ID
```

Note that the above RNA-seq data have been filtered to include only genes annotated the KEGG pathway database. Now we demostrate the basic usage of the main wrapper function `compute_n1pas` that completes the above workflow. Note that the mixture model fit in Panel D can be customized by passing `locfdr` specific parameters. For a small ontology, such as KEGG, custom parameters are provided. Simply set `small_ontology=TRUE` as in the example below. See the documentation via `?transform_gene_pathway` for more details.


```r
library(n1pas)
#> Welcome to the n1pas package!
#> Need help? Email the listserv mailing list nof1pathwayssupport@list.arizona.edu
#> Stackoverflow is a great place for general help: http://stackoverflow.com
## Retrieve the library location to access files needed for this tutorial
my_path <- find.package("n1pas")
annot_file <- file.path(my_path, "extdata/kegg_tb.txt")
desc_file <- file.path(my_path, "extdata/kegg.description_tb.txt")

iso_data <- blca_patient_iso_kegg

enrichment_profile <- compute_n1pas(iso_data = iso_data, annot_file = annot_file, desc_file = desc_file, iso_range = c(2,30), genes_range = c(15,500), small_ontology = TRUE, plot_locfdr = 1)
```

![plot of chunk example](vignettes/example-1.png)

```r

## Explore the top five dysregulated pathways for the first patient.
head(enrichment_profile, 15)
#>          pathway_score  fdr_value num_genes_annot num_genes_measured
#> hsa00670      4.566867 0.00000000              18                 15
#> hsa05211      2.354654 0.00000000              70                 46
#> hsa04540      2.164417 0.00000000              90                 58
#> hsa05213      2.121612 0.02662560              52                 39
#> hsa05218      2.072953 0.02794445              71                 47
#> hsa03020      2.022463 0.04838938              29                 20
#> hsa04115      1.961458 0.08247495              68                 46
#> hsa03440      1.928977 0.10835410              28                 18
#> hsa03022      1.866148 0.17434023              37                 21
#> hsa05110      1.856213 0.18781769              54                 29
#> hsa05219      1.856213 0.18781769              42                 29
#> hsa05221      1.826214 0.22851314              57                 40
#> hsa04962      1.820241 0.23661523              44                 24
#> hsa03040      1.778470 0.30825666             127                 68
#> hsa04530      1.761299 0.34170396             132                 85
#>          upper_fdr_threshold diff_splice_call
#> hsa00670            1.847233                1
#> hsa05211            1.847233                1
#> hsa04540            1.847233                1
#> hsa05213            1.847233                1
#> hsa05218            1.847233                1
#> hsa03020            1.847233                1
#> hsa04115            1.847233                1
#> hsa03440            1.847233                1
#> hsa03022            1.847233                1
#> hsa05110            1.847233                1
#> hsa05219            1.847233                1
#> hsa05221            1.847233                0
#> hsa04962            1.847233                0
#> hsa03040            1.847233                0
#> hsa04530            1.847233                0
#>                                      pathway_desc
#> hsa00670                One carbon pool by folate
#> hsa05211                     Renal cell carcinoma
#> hsa04540                             Gap junction
#> hsa05213                       Endometrial cancer
#> hsa05218                                 Melanoma
#> hsa03020                           RNA polymerase
#> hsa04115                    p53 signaling pathway
#> hsa03440                 Homologous recombination
#> hsa03022              Basal transcription factors
#> hsa05110                Vibrio cholerae infection
#> hsa05219                           Bladder cancer
#> hsa05221                   Acute myeloid leukemia
#> hsa04962 Vasopressin-regulated water reabsorption
#> hsa03040                              Spliceosome
#> hsa04530                           Tight junction
```

## Getting help

There are three main places to get help with n1pas:

1.  Read the function specific documentation included with the package. Including ways to customize.

2.  [Stack Overflow][so] is a great source of answers to common ggplot2
    questions. It is also a great place to get help, once you have
    created a reproducible example that illustrates your problem.
	
3.  Need help? Email the listserv mailing list nof1pathwayssupport@list.arizona.edu.
