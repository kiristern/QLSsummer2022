setwd('Desktop/')
getwd()

# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }

# BiocManager::install("spatialLIBD")

## Load the package
library("spatialLIBD")

create_csr <- function(mat){
    csr_count <- summary(mat)
    csr_counts <- data.frame(gene = rownames(genexspot_count)[csr_count$i], 
                        spot = colnames(genexspot_count)[csr_count$j],
                        count = csr_count$x)
}

## Download the spot-level data
spe <- fetch_data(type = "spe")

## This is a SpatialExperiment object
spe

# investigate data
gene_meta <- rowData(spe)
write.csv(gene_meta, "projects/spatialLIBD/data/gene_meta.csv")

# counts matrix
genexspot <- assays(spe)
genexspot_count <- genexspot$counts
genexspot_logcount <- genexspot$logcounts

# create csr
csr_counts <- create_csr(genexspot_count)
head(csr_counts)
write.csv(csr_counts, "projects/spatialLIBD/data/spatialLIBD_csr_counts.csv")

csr_logcounts <- create_csr(genexspot_logcount)
head(csr_logcounts)
write.csv(csr_logcounts, file = "projects/spatialLIBD/data/spatialLIBD_csr_logcounts.csv")

# get col metadata
spot_meta <- colData(spe)
head(spot_meta)
write.csv(spot_meta, "projects/spatialLIBD/data/spatialLIBD_spot_meta.csv")

spot_st <- spatialCoords(spe)
head(spot_st)
write.csv(spot_st, "projects/spatialLIBD/data/spatialLIBD_spot_st.csv")


##### SingleCellExperiment data #####
sce <- fetch_data(type = 'sce')

# create csr
genexcell <- assays(sce)
genexcell_count <- genexcell$counts
genexcell_logcount <- genexcell$logcounts

# count matrix
csr_counts_cell <- create_csr(genexcell_count)
head(csr_counts_cell)
write.csv(csr_counts_cell, "projects/spatialLIBD/data/spatialLIBD_csr_counts_cell.csv")

csr_logcounts_cell <- create_csr(genexcell_logcount)
head(csr_logcounts_cell)
write.csv(csr_logcounts_cell, file = "projects/spatialLIBD/data/spatialLIBD_csr_logcounts_cell.csv")

# gene meta
gene_meta_sce <- rowData(sce)
write.csv(gene_meta_sce, "projects/spatialLIBD/data/gene_meta_sce.csv")

# get col metadata
cell_meta <- colData(sce)
head(cell_meta)
write.csv(cell_meta, "projects/spatialLIBD/data/spatialLIBD_cell_meta.csv")

cell_st <- spatialCoords(sce)
head(cell_st)
write.csv(cell_st, "projects/spatialLIBD/data/spatialLIBD_cell_st.csv")