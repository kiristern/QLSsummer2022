setwd('Desktop/')
getwd()

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("spatialLIBD")

## Load the package
library("spatialLIBD")

## Download the spot-level data
spe <- fetch_data(type = "spe")

## This is a SpatialExperiment object
spe

genexspot <- assays(spe)
genexspot_count <- genexspot$counts
genexspot_logcount <- genexspot$logcounts

create_csr <- function(mat){
    csr_count <- summary(mat)
    csr_counts <- data.frame(gene = rownames(genexspot_count)[csr_count$i], 
                        spot = colnames(genexspot_count)[csr_count$j],
                        count = csr_count$x)
}

# create csr
csr_counts <- create_csr(genexspot_count)
head(csr_counts)
write.csv(csr_counts, "Desktop/spatialLIBD_csr_counts.csv")

csr_logcounts <- create_csr(genexspot_logcount)
head(csr_logcounts)
write.csv(csr_logcounts, file = "Desktop/spatialLIBD_csr_logcounts.csv")

# get col metadata
spot_counts <- colData(spe)
head(spot_counts)
write.csv(spot_counts, "Desktop/spatialLIBD_spot_counts.csv")

spot_st <- spatialCoords(spe)
head(spot_st)
write.csv(spot_st, "Desktop/spatialLIBD_spot_st.csv")


##### SingleCellExperiment data #####
sce <- fetch_data(type = 'sce')

# create csr
genexcell <- assays(sce)
genexcell_count <- genexcell$counts
genexcell_logcount <- genexcell$logcounts

csr_counts_cell <- create_csr(genexcell_count)
head(csr_counts_cell)
write.csv(csr_counts_cell, "Desktop/spatialLIBD_csr_counts_cell.csv")

csr_logcounts_cell <- create_csr(genexcell_logcount)
head(csr_logcounts_cell)
write.csv(csr_logcounts_cell, file = "Desktop/spatialLIBD_csr_logcounts_cell.csv")

# get col metadata
cell_counts <- colData(sce)
head(cell_counts)
write.csv(cell_counts, "Desktop/spatialLIBD_cell_counts.csv")

cell_st <- spatialCoords(sce)
head(cell_st)
write.csv(cell_st, "Desktop/spatialLIBD_cell_st.csv")