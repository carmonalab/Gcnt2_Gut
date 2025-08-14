library(dplyr)
library(ggplot2)
library(purrr)
# takes a single-cell (genes x cells) matrix with named rows (given gene names) and the path to Ensembl annotation (output of biomart with fields "Gene name", "Gene Synonym" and "Gene stable ID"
# Returns a matrix with filtered rows of genes with standard gene names. Of accepted synonyms were found, these are converted to standard gene names, as long as these are not redundant

standardizeGeneSymbols = function(matrix, EnsemblGeneTable=NULL, EnsemblGeneFile=NULL){
  
   #If file is given
  if (is.null(EnsemblGeneTable)) {
    if (is.null(EnsemblGeneFile)) {
       stop("Please provide EnsemblID table or file")
    }
    EnsemblGeneTable <- fread(EnsemblGeneFile)
  } 
     
  #Translate Ensembl IDs if necessary
  ens.format <- FALSE
  genes.in <- rownames(matrix)
  ngenes <- length(genes.in)
  
  ens.count <- length(intersect(genes.in, EnsemblGeneTable[["Gene stable ID"]]))
  gname.count <- length(intersect(genes.in, EnsemblGeneTable[["Gene name"]]))
  
  max <- max(ens.count, gname.count)
  if (max < length(genes.in)/2) {
    warning("Over 50% of genes in input object not found in reference gene table")
  }
  if (ens.count > gname.count) {
    ens.format <- TRUE
  }
  
  if (ens.count > gname.count) {  #Input object has Ensembl IDs
    genes.tr <- EnsemblGeneTable[["Gene name"]][match(genes.in, EnsemblGeneTable[["Gene stable ID"]])]
    names(genes.tr) <- genes.in
    
    genes.tr <- genes.tr[!is.na(genes.tr) & genes.tr != ""]
  } else {
    genes.tr <- genes.in
    names(genes.tr) <- genes.in
  }
  
  ###### 1. First match dictionary 
  geneRef_dict <- EnsemblGeneTable[["Gene name"]]
  names(geneRef_dict) <- EnsemblGeneTable[["Gene Synonym"]]
  geneRef_dict <- geneRef_dict[!is.null(names(geneRef_dict))]
  
  message(paste("Number of genes in input matrix:", ngenes))
  genesAllowList1 <- genes.tr[!is.na(genes.tr) & genes.tr != "" &
                                genes.tr %in% EnsemblGeneTable[["Gene name"]]] #keep genes with standard Gene name
  l <- length(genesAllowList1)
  
  message(sprintf("Number of genes with standard symbols: %i (%.2f%%)", l, l/ngenes*100))

  if (l < ngenes & !ens.format){
    message(paste("Examples of non-standard gene names:"))
    message(paste(head(genes.tr[ !genes.tr %in% EnsemblGeneTable[["Gene name"]] ])))
  }
  
  ###### 2. Search among synonyms
  genesAllowList2 <- genes.tr[!genes.tr %in% EnsemblGeneTable[["Gene name"]] & 
                                genes.tr %in% EnsemblGeneTable[["Gene Synonym"]]] # keep genes with accepted gene name synonym
  genesAllowList2.gn <- geneRef_dict[genesAllowList2] # translate gene synonym to standard gene name
  
  message(paste("Additional number of genes with accepted gene name synonym: ",length(genesAllowList2.gn)))
  
  #Names of genesAllowList contain IDs in matrix - elements contain the new names
  genesAllowList <- c(genesAllowList1,genesAllowList2.gn)
  
  ###### 3. Check for duplicates
  is.dup <- duplicated(genesAllowList)
  genesAllowList <- genesAllowList[!is.dup]
  message(sprintf("Number of duplicated gene name: %i (%.2f%%)", sum(is.dup), sum(is.dup)/ngenes*100))
  
  l <- length(genesAllowList)
  message(sprintf("Final number of genes: %i (%.2f%%)", l, l/ngenes*100))
  
  ###### 4. Subset matrix for allowed genes, and translate names
  rows.select <- rownames(matrix)[rownames(matrix) %in% names(genesAllowList)]
  matrix <- matrix[rows.select, ]
  rownames(matrix) <- genesAllowList[rows.select]
  
  return(matrix)
}

#Calculate variable genes, excluding genes from a blacklist (e.g. ribosomal), and poorly expressed genes
select.variable.genes = function(obj, nfeat=1500, blacklist=NULL, plot=TRUE, min.exp=0.01, max.exp=3){
  
  obj <- FindVariableFeatures(obj, nfeatures = nfeat)
  
  varfeat <- obj@assays$RNA@var.features
  
  if (!is.null(blacklist)) {
    removeGenes1 <- varfeat[varfeat %in% unlist(blacklist)]
    varfeat <- setdiff(varfeat, removeGenes1)
    
    print("Remove from variable genes:")
    print(removeGenes1)
  }
  #Also remove genes that are very poorly or always expressed (=not really variable genes)
  means <- apply(obj@assays$RNA@data[varfeat,], 1, mean)
  
  if (plot) {
    par(mfrow=c(1,2))
    hist(means, breaks=30)
    hist(log10(means), breaks=30)
    par(mfrow=c(1,1))
  }
  removeGenes2 <- names(means[means<min.exp | means>max.exp])
  
  varfeat <- setdiff(varfeat, removeGenes2 )
  print("Remove from variable genes:")
  print(removeGenes2)
  
  obj@assays$RNA@var.features <- varfeat
  
  return(obj)
} 


# Read an xlsx document with multiple tabs into a list of dataframes
read_excel_allsheets <- function(filename,
                                 tibble = FALSE,
                                 skip = 0,
                                 col_names = TRUE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  
  x <- lapply(sheets, function(X){
    readxl::read_excel(filename,
                       sheet = X,
                       skip = skip,
                       col_names = col_names)
          }
    )
  
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


#################################################################################
#################################################################################
# Plot functions
## QC metrics distribution
qc_plot <- function(metadata, #seurat object metadata
                    var, # QC variables to plot
                    cutoffs,  # list of cutoffs with same names as QC vars
                    fill = NULL, # groupin variable for density plots
                    show.legend = T,
                    split = F
                    ) {
  

  metadata %>%
    ggplot(aes(.data[[var]])) +
    {if(is.null(fill)){
      geom_density(fill = "grey50",
                   color = "grey30")
    } else {
      geom_density(aes(fill = .data[[fill]]),
                   color = "grey10",
                   alpha = 0.5,
                   show.legend = show.legend)
    }}+
    {if(max(metadata[[var]])>100){
          scale_x_log10()}} +
    {if(!all(is.na(cutoffs[[var]]))){
      geom_vline(xintercept = cutoffs[[var]],
                 color = "red",
                 linetype="dashed")}} +
    {if(split){
      facet_wrap(~ .data[[fill]])}} +
    theme_bw() +
    theme(aspect.ratio = 1)
}


#################################################################################
#################################################################################
#' Write a dgCMatrix to an h5 file similar to cellRanger format
# Source: https://github.com/AllenInstitute/scrattch.io/blob/master/R/write_10x.R
if (!require("rhdf5", quietly = TRUE)){
  BiocManager::install("rhdf5")
}


write_dgCMatrix_h5 <- function(mat,
                               cols_are = "gene_names",
                               h5_target,
                               ref_name = "mm10-1.2.0_premrna",
                               gene_ids = NULL) {
  
  #library(Matrix)
  
  if(grepl("gene",cols_are)) {
    mat <- Matrix::t(mat)
  }
  
  # Create target file
  rhdf5::h5createFile(h5_target)
  # Create data group
  rhdf5::h5createGroup(h5_target,
                       ref_name)
  
  # Store sample ids (barcodes) and gene names
  rhdf5::h5write(colnames(mat),
                 h5_target,
                 paste0("/",ref_name,"/barcodes"))
  rhdf5::h5write(rownames(mat),
                 h5_target,
                 paste0("/",ref_name,"/gene_names"))
  
  if(is.null(gene_ids)) {
    gene_ids <- rownames(mat)
  }
  
  rhdf5::h5write(gene_ids,
                 h5_target,
                 paste0("/",ref_name,"/gene"))
  
  # Store dimensions as shape
  rhdf5::h5write(dim(mat),
                 h5_target,
                 paste0("/",ref_name,"/shape"))
  
  # Store values from mat@x as data
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/data"),
                         dims = length(mat@x),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@x,
                 h5_target,
                 paste0("/",ref_name,"/data"))
  
  # Store row indices from mat@i as indices
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/indices"),
                         dims = length(mat@i),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@i,
                 h5_target,
                 paste0("/",ref_name,"/indices"))
  
  # Store column pointers from mat@p as indptr
  rhdf5::h5write(mat@p,
                 h5_target,
                 paste0("/",ref_name,"/indptr"))
  
}


cell_type_classification <- function(object, cell_type, ref_map, 
                                     ident_col = "scGate_multi", 
                                     filter.cells = FALSE){
  
  # Assuming that your Seurat object was scGated  
  Seurat::Idents(object) <- object@meta.data[[ident_col]]
  
  # Subsetting only relevant cell type
  cells_of_interest <- subset(object, idents = cell_type)
  
  # ProjecTILs classification using passed reference map
  data("Hs2Mm.convert.table", package = "ProjecTILs")
  cells_of_interest <- ProjecTILs::ProjecTILs.classifier(query = cells_of_interest, 
                                                         ref = ref_map, 
                                                         filter.cells	= filter.cells)
  
  return(cells_of_interest)
}

multiple.wilcox.tests <- function(groups_df, factor_col, num_col, pairs_list){
  group1 <- pairs_list[1]
  group2 <- pairs_list[2]
  
  data1 <- groups_df %>%
    dplyr::filter(.data[[factor_col]] == group1) %>%
    dplyr::pull(.data[[num_col]])  
  
  data2 <- groups_df %>%
    dplyr::filter(.data[[factor_col]] == group2) %>%
    dplyr::pull(.data[[num_col]]) 
  
  test <- wilcox.test(data1, data2)
  
  return(tibble(
    group1 = group1,
    group2 = group2,
    p.value = test$p.value
  ))
} 





