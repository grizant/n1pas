#' 8 paired RNA-seq samples from 4 BRCA_TCGA patients.
#'
#' A RNA-seq data set containing read counts for 12 pairs of transcriptomes
#' from 12 breast cancer (BRCA) patients.
#' All gene expression counts have been adjusted for multi-read assignment by
#' RSEM (RNA-Seq by Expectation-Maximization). Otherwise, data are unnormalized.
#' Platform used was Illumina RNA-seq V2. Downloaded Dec 2014.
#' Rows are labeled by HGNC gene symbol.
#'
#' @format A data frame with 20501 rows and 24 variables:
#' \describe{
#' \item{TCGA.A7.A0DB.N}{The counts from A7.A0DB non-tumor breast tissue, N=normal}
#'   \item{TCGA.A7.A0DB.T}{The counts from A7.A0DB tumor breast tissue, T=tumor}
#' \item{TCGA.AC.A2FM.N}{The counts from AC.A2FM non-tumor breast tissue, N=normal}
#'   \item{TCGA.AC.A2FM.T}{The counts from AC.A2FM tumor breast tissue, T=tumor}
#' \item{TCGA.BH.A1EU.N}{The counts from BH.A1EU non-tumor breast tissue, N=normal}
#'   \item{TCGA.BH.A1EU.T}{The counts from BH.A1EU tumor breast tissue, T=tumor}
#'  \item{TCGA.E2.A153.N}{The counts from E2.A153 non-tumor breast tissue, N=normal}
#'   \item{TCGA.E2.A153.T}{The counts from E2.A153 tumor breast tissue, T=tumor}
#'   
#' }
#' @source \url{https://tcga-data.nci.nih.gov}
"blca_patient_iso_kegg"