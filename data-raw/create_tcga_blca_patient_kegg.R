## Create isoform-level RNA-seq data set of one TCGA BLCA patient to use in examples
## AG Schissler, 2 Apr 2019

## Isoform RNA-seq data for one TCGA BLCA patient
## Data already filtered to genes annotated in KEGG to reduce size
blca_patient_iso_kegg <- readRDS("inst/extdata/blca_patient_iso_kegg.rds")
## remove the isoformID column as the rows are labeled with this information
blca_patient_iso_kegg <- blca_patient_iso_kegg[,-1]
usethis::use_data(blca_patient_iso_kegg, overwrite = TRUE)

## Derived from manuscript files:
## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_kegg_data.RData")
## blca_patient_iso_kegg <- blca_iso_kegg_data[,1:4]
## saveRDS(blca_patient_iso_kegg, file = "/Users/alfred/n1pas/inst/extdata/blca_patient_iso_kegg.rds")
