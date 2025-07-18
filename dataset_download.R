# Downloading human scRNA-seq human dataset
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE214695&format=file"
file_name <- "raw_scRNA_seq_human.tar"
file_path <- "/Users/aleksandarmihaylov/Documents/Carmona Lab/Gcnt2_Gut/datasets/"
options(timeout = 300)
download.file(url, paste(file_path, file_name, sep = ""), mode  = "wb")

# Downloading author's human cells annotation dataset
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE214695&format=file&file=GSE214695%5Fcell%5Fannotation%2Ecsv%2Egz"
file_name <- "cells_annotation_human.csv"
file_path <- "/Users/aleksandarmihaylov/Documents/Carmona Lab/Gcnt2_Gut/datasets/"
download.file(url, paste(file_path, file_name, sep = ""), mode  = "wb")

# Downloading human bulk RNA-seq dataset
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235236&format=file"
file_name <- "bulk_RNA_seq_human.tar"
file_path <- "/Users/aleksandarmihaylov/Documents/Carmona Lab/Gcnt2_Gut/datasets/"
download.file(url, paste(file_path, file_name, sep = ""), mode  = "wb")

