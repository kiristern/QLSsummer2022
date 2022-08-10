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
head(gene_meta)
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
colnames(spot_meta)
write.csv(spot_meta, "projects/spatialLIBD/data/spatialLIBD_spot_meta.csv")

spot_st <- spatialCoords(spe)
head(spot_st)
write.csv(spot_st, "projects/spatialLIBD/data/spatialLIBD_spot_st.csv")


##### SingleCellExperiment data #####
sce <- fetch_data(type = 'sce')
sce

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

# layer-level data
sce_layer <- fetch_data(type="sce_layer")
sce_layer
head(rowData(sce_layer))
head(colData(sce_layer))



## Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()

# results from 'enrichment model' (t-stat: assess whether gene had higher expression in a given layer compared to the rest);
#              'pariwise model' (t-stat: assess whether gene had higher expression between one layer and another layer);
#              'anova model' (F-stat: expression variability across all layers)
modeling_results <- fetch_data("modeling_results", eh = ehub)

## list of modeling result tables
sapply(modeling_results, class)
sapply(modeling_results, dim)

sapply(modeling_results, function(x) {
    head(colnames(x))
})

## Convert modeling stats from wide to long format and extract significant genes
## This takes a few seconds to run
system.time(
    sig_genes <-
        sig_genes_extract_all(
            n = nrow(sce_layer),
            modeling_results = modeling_results,
            sce_layer = sce_layer
        )
)
head(sig_genes) # see http://research.libd.org/spatialLIBD/articles/spatialLIBD.html#extract-significant-genes for col details

## Explore the enrichment t-statistics derived from Tran et al's snRNA-seq
## DLPFC data
dim(tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer)

## Compute the correlation matrix of enrichment t-statistics between our data
## and Tran et al's snRNA-seq data
cor_stats_layer <- layer_stat_cor(
    tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer,
    modeling_results,
    "enrichment"
)


## Visualize the estimated number of cells per spot
vis_gene(
    spe = spe,
    sampleid = "151673",
    geneid = "cell_count"
)