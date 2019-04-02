#' Normal-tumor paired isoform-level RNA-seq for a BLCA TCGA patient.
#'
#' An isoform-level RNA-seq data set containing read counts for a pair of transcriptomes
#' from a bladder cancer (BLCA) patients.
#' All gene expression counts have been adjusted for multi-read assignment by
#' RSEM (RNA-Seq by Expectation-Maximization). Otherwise, data are unnormalized.
#' Platform used was Illumina RNA-seq V2. Downloaded Dec 2014.
#' Rows are labeled by isoform ID and further annotated to HGNC gene symbols.
#' Genes were filtered to include only those in the KEGG pathway ontology.
#'
#' @format A data frame with 18823 rows and 3 variables:
#' \describe{
#' \item{row.names}{isoform ID}
#'   \item{geneSymbol}{HUGO gene symbol corresponding to the isoform ID}
#' \item{TCGA-BL-A13J-N}{The counts from a AC.A2FM non-tumor bladder tissue, N=normal}
#'   \item{TCGA-BL-A13J-T}{The counts from AC.A2FM tumor bladder tissue, T=tumor}
#' }
#' @source \url{https://tcga-data.nci.nih.gov}
"blca_patient_iso_kegg"
