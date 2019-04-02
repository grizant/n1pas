## load utilities.R for reading
devtools::load_all()

## KEGG pathways
kegg_annot <- read_gene_set("inst/extdata/kegg_tb.txt")
usethis::use_data(kegg_annot, overwrite = TRUE)

## KEGG pathway descriptions to annotate output
kegg_desc <- suppressWarnings(read_gene_set("inst/extdata/kegg.description_tb.txt"))
rownames(kegg_desc) <- kegg_desc$path_id
usethis::use_data(kegg_desc, overwrite = TRUE)
