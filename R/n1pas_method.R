#' Compute Hellinger distance for gene-specific isoforms given pair of measurements
#'
#' \code{compute_hellinger} transforms a p_value vector to signed z-score vector.
#'
#' @details This function provides the Hellinger distance between two samples.
#' This distance is used to quantify differential isoform usage.
#' PROVIDE FORMULA
#' 
#' @param X a N x matrix of isoform counts, values restricted between 0 and 
#' @param direction a numeric vector indicating the direction of differentially expression
#' with 1 indicating higher expressed and -1 indicating lower expressed in the 'case' sample. 
#'
#' @return a numeric vector of signed z-scores
#'
#' @examples
#' set.seed(44)
#' N <- 100
#' p_values <- stats::runif(N)
#' directions <- sample(c(-1,1), size = N, replace = TRUE)
#' zscores <- compute_zscore(p_values, directions)
#' \dontrun{
#' stats::qqnorm(zscores); abline(0,1)
#' }
#' 
#' @export
compute_hellinger <- function(proportion_mat, only_expressed = T) {
    ## if there is more than one isoform compute the hellinger distance
    if (length(dim(proportion_mat)) == 2) {
        ## check that the gene is expressed in both samples
        if (only_expressed) {
            if (colSums(proportion_mat)[1] <= 0 | colSums(proportion_mat)[2] <= 0) {
                ## only alternatively spliced genes are of interest
                ## hel_dist <- 0
                hel_dist <- NA
            } else {
                hel_dist <- sqrt(0.5*sum((sqrt(proportion_mat[,1]) - sqrt(proportion_mat[,2]))^2))
            }
        } else {
            hel_dist <- sqrt(0.5*sum((sqrt(proportion_mat[,1]) - sqrt(proportion_mat[,2]))^2))
        }
    } else hel_dist <- NA
    ## return the value
    return(hel_dist)
}

#' Compute Hellinger distance for gene-specific isoforms given pair of measurements
#'
#' \code{compute_hellinger} transforms a p_value vector to signed z-score vector.
#'
#' @details This function provides the Hellinger distance between two samples.
#' This distance is used to quantify differential isoform usage.
#' PROVIDE FORMULA
#' 
#' @param X a N x matrix of isoform counts, values restricted between 0 and 
#' @param direction a numeric vector indicating the direction of differentially expression
#' with 1 indicating higher expressed and -1 indicating lower expressed in the 'case' sample. 
#'
#' @return a numeric vector of signed z-scores
#'
#' @examples
#' set.seed(44)
#' N <- 100
#' p_values <- stats::runif(N)
#' directions <- sample(c(-1,1), size = N, replace = TRUE)
#' zscores <- compute_zscore(p_values, directions)
#' \dontrun{
#' stats::qqnorm(zscores); abline(0,1)
#' }
#' 
#' @export
get_OR <- function(tmp_genes, dist_data, genes_range = c(15,500)) {
            measured_genes <- tmp_genes[tmp_genes %in% rownames(dist_data)]
            ## genes outside pathway
            not_genes <- rownames(dist_data)[!(rownames(dist_data) %in% measured_genes)]
            ## impose genes range
            if (length(measured_genes) >= genes_range[1] & length(measured_genes) <= genes_range[2]) {
                ## calculate counts
                (x11 <- sum(dist_data[measured_genes, "call"]))
                (x21 <- sum(dist_data[measured_genes, "call"] == 0))
                (x12 <- sum(dist_data[not_genes, "call"]))
                (x22 <- sum(dist_data[not_genes, "call"] == 0))
                ## sum(x11, x21, x12, x22) == nrow(dist_data)
                ## compute FET
                ## tmp_fet <- fisher.test(x = matrix( c(x11, x21, x12, x22), nrow = 2, ncol = 2))
                ## to_return <- tmp_fet$estimate
                ## names(to_return) <- NULL
                ## compute odds ratio manually
                (to_return <-  (x11/x21)/(x12/x22))
            } else to_return <- NA
            return(list(odds_ratio = to_return, num_genes = length(measured_genes)))
}

#' N-of-1-\emph{pathways} Mahalanobis distance (MD) Calculation
#'
#' \code{compute_nof1_md} computes N-of-1-\emph{pathways} Mahalanobis distance (MD)
#' statistic, direction, and P-value for a pathway (Schissler et al. 2015).
#'
#' @details This function enables gene set (pathway) testing for a single pathway while providing a context-meaningful effect size. In the original publication by Schissler et al. 2015, the P-value was produced via a nonparametric bootstrapping procedure. However, in later yet unpublished studies, it has been determined that a simple paired \emph{t}-test performs as well and much faster. Therefore the \emph{t}-test has been implemented. For a nonparametric approach use N-of-1-\emph{pathways} Wilcoxon for Gardeux et al. 2014.
#' 
#' @param baseline A vector (paired with case) containing all gene expression
#' @param case A vector (paired with baseline) containing all gene expression
#' @param ... Other arguments for stats::t.test.
#'
#' @return A list with the pathway_score (average MD score),
#' direction (up or down relative to the baseline),
#' and p_value. P_value is computed via a t.test.
#' 
#' @references
#' Gardeux, Vincent, et al. "'N-of-1-\emph{pathways}' unveils personal deregulated mechanisms from a single pair of RNA-Seq samples: towards precision medicine." Journal of the American Medical Informatics Association 21.6 (2014): 1015-1025.
#' 
#' Schissler, A. Grant, et al. "Dynamic changes of RNA-sequencing expression for precision medicine: N-of-1-\emph{pathways} Mahalanobis distance within pathways of single subjects predicts breast cancer survival." Bioinformatics 31.12 (2015): i293-i302.
#'
#' @seealso \code{link{compute_nof1_wilcoxon}}, \code{link{compute_nof1_pathways}}, \code{link{compute_multi_nof1_pathways}}
#'
#' @examples
#' ## compute N-of-1-pathways MD for a single pathway
#' set.seed(44)
#' N <- 10
#' gene1 <- stats::rchisq(N, 5)
#' gene2 <- stats::rchisq(N, 5)
#' md_res <- compute_nof1_md(gene1, gene2)
#' print(md_res)
#' 
#' @export
transform_iso_gene <- function(gene_list, expr_threshold = 0, ...) {

    ## 1. remove low expressing genes
    if (expr_threshold > 0) {
        print("removing low expressing genes")
        ## tmp_gene <- gene_list[[1]]
        keep_logic <- unlist(lapply(gene_list, function(tmp_gene) {
            ## sum the isoform level to find the gene-level expression
            (colSums(tmp_gene)[1] >= expr_threshold) & (colSums(tmp_gene)[2] >= expr_threshold)
        }))
        num_genes_removed <- length(keep_logic) - sum(keep_logic)
        prop_genes_removed <- round(num_genes_removed/length(keep_logic),3)
        print(paste("Removed", num_genes_removed, "genes, proportion of total is", prop_genes_removed))
        gene_list <- gene_list[keep_logic]
    }
    
    ## 2. transform into a distance per gene
    ## first compute the relative expression within a gene
    rel_list <- lapply(gene_list, function(tmp_gene) {
        apply(tmp_gene, 2, function(tmp_col) {
            ## check that there is no division by 0 
            if (sum(tmp_col) > 0) {
                tmp_col/sum(tmp_col)
            } else tmp_col
        })  
    })

    ## 3. Compute the Hellinger distance
    to_return <- unlist(lapply(rel_list, compute_hellinger))
    ## remove missing distances (due to non-expression or low expression in either sample)
    to_return <- to_return[!is.na(to_return)]

    ## return the gene-level metrics
    return(to_return)
}

#' Transforms gene-level distances into a pathway enrichment profiles
#'
#' \code{transform_gene_pathway} computes odds ratios of enrichment with alternatively spliced genes. Then fits a mixture model to identify enriched pathways.
#'
#' @details This function is the main workhorse of the package. It first computes odds ratios of enrichment then fits a mixture model via the \code{locfdr} package. For large ontologies (greater than 1000 gene sets), the defaults of locfdr are likely to be adequate. For small ontologies we have developed custom parameters (set small_ontology to TRUE). The fit can be customized by setting \code{custom_locfdr} to TRUE and \code{small_ontology} to FALSE. Please see locfdr for more details.
#' 
#' @param gene_dist A named vector of Hellinger distances. Names are HUGO gene symbols
#' @param annot_file The path and file name for the gene set annotations
#' @param desc_file The path and file name for gene set descriptions
#' @param genes_range a two-component vector to filter gene sets based on minimum and maximum number of genes. Default is (15, 500).
#' @param fdr_threshold A numeric atomic threshold for local fdr
#' @param plot_locfdr A numeric indicating whether the mixture model should be plotted. 0 = no, 1 = yes. See locfdr for more options and details.
#' @param small_ontology Logical. Set to TRUE to use small ontology customizations and FALSE to not.
#' @param custom_locfdr Logical. Set to TRUE to pass your own locfdr parameters and FALSE to use the defaults or \code{small_ontology} customizations.
#' @param ... Other arguments for n1pas, locfdr, and other package's functions
#'
#' @return pathway (gene set) enrichment profile organized in a data frame with 7 variables sorted by local fdr value.
#' \describe{
#'   \item{row.names}{rows are labeled by user-specified gene set ID, the "path_id"}
#'   \item{pathway_score}{numeric. N1PAS's effect size: an odds ratio quantifying pathawy enrichment of alternatively spliced genes}
#'   \item{fdr_value}{numeric. local false discovery rate. See Efron 2013}
#'   \item{num_genes_annot}{numeric. Number of genes annotated to the pathway from gene set input}
#'   \item{num_genes_measured}{numeric. Number of genes measured and expressed (in at least one sample) within the pathway.}
#'   \item{upper_fdr_threshold}{numeric. Odds ratios larger than this threshold are identified as enriched with alternatively spliced genes.}
#'   \item{diff_splice_call}{numeric. Binary values with 1 being pathway is enriched and 0 otherwise.}
#'  \item{pathway_desc}{character. Description of the pathway provided by the user.}
#' }
#'
#' @references
#' #' Efron, B. (2013). “Local False Discovery Rates,” in Large-Scale Inference doi:10.1017/cbo9780511761362.006
#' 
#' Schissler, A. Grant, et al. "A single-subject method to detect pathway enriched with alternatively spliced genes." Frontiers in Genetics, to appear, (2019)
#' Gardeux, Vincent, et al. "'N-of-1-\emph{pathways}' unveils personal deregulated mechanisms from a single pair of RNA-Seq samples: towards precision medicine." Journal of the American Medical Informatics Association 21.6 (2014): 1015-1025.
#' 
#' Schissler, A. Grant, et al. "Dynamic changes of RNA-sequencing expression for precision medicine: N-of-1-\emph{pathways} Mahalanobis distance within pathways of single subjects predicts breast cancer survival." Bioinformatics 31.12 (2015): i293-i302.
#'
#' @seealso \code{link{compute_nof1_wilcoxon}}, \code{link{compute_nof1_pathways}}, \code{link{compute_multi_nof1_pathways}}
#'
#' @export
transform_gene_pathway <- function(gene_dist, annot_file, desc_file = NULL, genes_range = c(15, 500), fdr_threshold = 0.2, plot_locfdr = 0, small_ontology = TRUE, custom_locfdr = FALSE, ...) {

    ## i. set locfdr parameters for small ontology 
    if (small_ontology) {
        bre = 8
        df = 4
        pct = 0
        pct0 = 1/4
        nulltype = 2
        use_sd_mlests = TRUE
    }
    
    ## 1. Structure gene set definitions
    ## read in gene set (pathway) annotation
    annot_data <- read.delim2(file = annot_file, stringsAsFactors = F)

    ## restructure into a list of pathways
    ## make path_id_name a variable, path_id_name <- "path_id"
    id_column <- 1 ## standard is that it will be the first column
    annot_list <- split(annot_data$symbol, annot_data[,id_column])

    ## store the number of genes annotated
    annot_lengths <- unlist(lapply(annot_list, function(tmp_genes) {length(unique(tmp_genes))}))
    ## count genes actually measured
    measured_lengths <- unlist(lapply(annot_list, function(tmp_genes) {
        tmp_genes <- as.character(tmp_genes)
        sum(tmp_genes %in% names(gene_dist))
    }))

    ## 2. Quantify enrichment using locfdr

    ## 2.1 Cluster genes
    gene_clusters <- stats::kmeans(x = gene_dist, 2)
    dist_data <- data.frame(gene_dist, call = 0)
        
    ## determine which cluster has larger center
    alt_center <- which.max(gene_clusters$centers)
    dist_data[names(gene_clusters$cluster), "call"] <- ifelse(gene_clusters$cluster == alt_center, 1, 0)

    ## compute odds ratios
    odds_ratio <- lapply(annot_list, get_OR, dist_data, genes_range=genes_range)

    ## count genes 
    num_genes <- unlist(lapply(odds_ratio, function(tmp_id) {tmp_id$num_genes}))
    odds_ratio <- unlist(lapply(odds_ratio, function(tmp_id) {tmp_id$odds_ratio}))
        
    ## detect outlying pathways above prior to fitting
    fil_odds <- odds_ratio[!is.na(odds_ratio)]
    outliers <- grDevices::boxplot.stats(fil_odds)$out
    if (length(outliers) > 0) fil_odds <- fil_odds[!(names(fil_odds) %in% names(outliers))]

    ## 2.2 Fit mixture model using locfdr package
    ## Suppress warnings to complete in parallel

    ## defaults are for a small ontology, such as KEGG

    ## allow a user to customize locfdr
    if (custom_locfdr) {
        tmp_locfdr <- locfdr::locfdr(zz = fil_odds, plot = plot_locfdr, ...)
    } else {
        if (small_ontology) {
            ## or use the small ontology customizations
            if (use_sd_mlests) {mlests <- c(1, stats::sd(fil_odds))}
            suppressWarnings(tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/bre), df = df, pct = pct, pct0 = pct0, nulltype = nulltype, plot = plot_locfdr, mlests = mlests))
        } else {
            ## or use the program defaults, suitable for large ontologies (> 1000)
            tmp_locfdr <- locfdr::locfdr(zz = fil_odds, plot = plot_locfdr)
        }
    }

    ## 2.3 use the fit to make decisions
    tmp_upper <- tmp_locfdr$z.2[2]
    ## then indicate the hits the hits
    tmp_call <- rep(0, length(odds_ratio))
    names(tmp_call) <- names(odds_ratio)
    if (!is.na(tmp_upper)) {
        odds_ratio[which(odds_ratio >= tmp_upper)]
        tmp_call[which(odds_ratio >= tmp_upper)] <- 1
    } else {
        if (length(outliers) > 0) {
            tmp_call[names(outliers)] <- 1
        }
    }

    ## 3. structure as a data frame
    to_return <- data.frame(pathway_score = odds_ratio, fdr_value = 1, num_genes_annot = annot_lengths, num_genes_measured = measured_lengths, row.names = names(odds_ratio), upper_fdr_threshold = tmp_upper, diff_splice_call = tmp_call, stringsAsFactors = F)
    ## add fdr_values
    to_return[names(fil_odds), "fdr_value"] <- tmp_locfdr$fdr
    ## and give p-value of 0 for outliers
    if (length(outliers) > 0) {to_return[names(outliers), "fdr_value"] <- 0}
    ## sort by lowest fdr value
    to_return <- to_return[order(to_return$pathway_score, decreasing = T),]

    ## 4. add pathway descriptions
    to_return$pathway_desc <- "N/A"
    if (!is.null(desc_file)) {
        desc_data <- utils::read.delim2(file = desc_file, stringsAsFactors = F)
        rownames(desc_data) <- desc_data[, id_column]
        ## check that scored pathways have annotation
        if (any(rownames(to_return) %in% desc_data$path_id)) {
            ## add to the output
            to_return$pathway_desc <- desc_data[rownames(to_return), "description"]
        } else {
            warning("Pathway ids in description data do not match ids in gene set definitions.")
        }
    }
    
    ## 5. return the pathway-level metrics
    return(to_return)
}

#' Execute N-of-1-\emph{pathways} Alternatively Spliced (N1PAS)
#'
#' \code{compute_n1pas} transforms a pair of isoform expression data sets
#' to a profile of pathways enriched with alternatively spliced genes.
#'
#' @details This function performs the entire N1PAS workflow for a pair of samples.
#' 
#' @param iso_data A data frame containing gene symbols, 'case' and 'baseline' isoform measurements. Each row corresponds to a unique isoform.
#' @param annot_file the path and file name for the gene set annotations
#' @param desc_file optionally provide the path and file name for gene set descriptions
#' @param iso_range a two-component vector to filter genes based on minimum and maximum number of isoforms. Default is (2, 30).
#' @param genes_range a two-component vector to filter gene sets based on minimum and maximum number of genes. Default is (15, 500).
#' @param ... Other arguments for other n1pas functions
#'
#' @return pathway (gene set) enrichment profile organized in a data frame with 7 variables sorted by local fdr value.
#' \describe{
#'   \item{row.names}{rows are labeled by user-specified gene set ID, the "path_id"}
#'   \item{pathway_score}{numeric. N1PAS's effect size: an odds ratio quantifying pathawy enrichment of alternatively spliced genes}
#'   \item{fdr_value}{numeric. local false discovery rate. See Efron 2013}
#'   \item{num_genes_annot}{numeric. Number of genes annotated to the pathway from gene set input}
#'   \item{num_genes_measured}{numeric. Number of genes measured and expressed (in at least one sample) within the pathway.}
#'   \item{upper_fdr_threshold}{numeric. Odds ratios larger than this threshold are identified as enriched with alternatively spliced genes.}
#'   \item{diff_splice_call}{numeric. Binary values with 1 being pathway is enriched and 0 otherwise.}
#'  \item{pathway_desc}{character. Description of the pathway provided by the user.}
#' }
#'
#' @note This function requires does not require library-size normalized data as the Hellinger distance is invariant under scaling.
#' 
#' @references
#' Efron, B. (2013). “Local False Discovery Rates,” in Large-Scale Inference doi:10.1017/cbo9780511761362.006
#' 
#' Schissler, A. Grant, et al. "A single-subject method to detect pathway enriched with alternatively spliced genes." Frontiers in Genetics, to appear, (2019)
#'
#' Gardeux, Vincent, et al. "'N-of-1-\emph{pathways}' unveils personal deregulated mechanisms from a single pair of RNA-Seq samples: towards precision medicine." Journal of the American Medical Informatics Association 21.6 (2014): 1015-1025.
#' 
#' Schissler, A. Grant, et al. "Dynamic changes of RNA-sequencing expression for precision medicine: N-of-1-\emph{pathways} Mahalanobis distance within pathways of single subjects predicts breast cancer survival." Bioinformatics 31.12 (2015): i293-i302.
#'
#'
#' @seealso \code{\link{compute_hellinger}} \code{\link{get_OR}} \code{link{transform_iso_gene}} \code{link{transform_gene_pathway}} \code{\link{blca_patient_iso_kegg}}
#'
#' @export
compute_n1pas <- function(iso_data, annot_file, desc_file = NULL, iso_range = c(2,30), genes_range = c(15,500), ...) {

    ######################################################################
    ## i. pre-processing
    ## restructure into a list of genes
    ## gene symbols should appear as the first column
    gene_list <- split(iso_data[,-1], iso_data[,1])

    ## impose min/max number of isoforms and at least one alternative isoform
    filter_logic <- unlist(lapply(gene_list, function(tmp_gene) {
        ## retrive number of isoforms
        tmp_num_iso <- dim(tmp_gene)[1]
        tmp_num_iso < iso_range[1] | tmp_num_iso > iso_range[2]
    }))

    genes_to_keep <- names(filter_logic)[!filter_logic]
    gene_list <- gene_list[genes_to_keep]
    
    ## 1. Find genewise distances
    gene_dist <- transform_iso_gene(gene_list = gene_list, ...)

    ## gene_dist <- transform_iso_gene(X = gene_list, method = gene_method)
 
    ## 2. Compute pathway-level metrics
    to_return <- transform_gene_pathway(gene_dist = gene_dist, annot_file = annot_file,
                                        desc_file = desc_file, method = pathway_method,
                                        genes_range = genes_range, ...)

    ## 3. Return scored and sorted pathways
    return(to_return)
}
