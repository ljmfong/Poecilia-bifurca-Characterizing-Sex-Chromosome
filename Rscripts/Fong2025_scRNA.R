# Install Packages ---------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

install.packages('Seurat')

install.packages('grr')

install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)

install.packages('Matrix')

install.packages('dplyr')

install.packages('magrittr')

install.packages("tidyverse")

BiocManager::install("DESeq2")

install.packages("pheatmap")

BiocManager::install("apeglm")

install.packages("RColorBrewer")

BiocManager::install("scran")

BiocManager::install("scuttle")

install.packages("ashr")

BiocManager::install("scater")

# Load Libraries ------------
library(SingleCellExperiment)
library(Seurat)
library(Matrix.utils)
library(Matrix)
library(dplyr)
library(magrittr)
library(purrr)
library(stringr)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(apeglm)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scran)
library(scuttle)
library(ashr)
library(scater)
library(tidyverse)
# Load and view the data --------
load("gonad_umap_3_rename_clusters.RData")
View(gonad_umap_3_rename_clusters@meta.data) # View all the information
head(gonad_umap_3_rename_clusters) # Summary of information

metadata_gonad_unfilt <- gonad_umap_3_rename_clusters@meta.data # Subset the meta data
summary(metadata_gonad_unfilt) # Look at the data before filtering

# Output:
#seq_folder             nUMI            nGene         log10GenesPerUMI     mtUMI           mitoRatio      
#Length:28014       Min.   :   500   Min.   :  100.0   Min.   :0.4976   Min.   :   0.00   Min.   :0.00000  
#Class :character   1st Qu.:   613   1st Qu.:  346.0   1st Qu.:0.8680   1st Qu.:   3.00   1st Qu.:0.00346  
#Mode  :character   Median :   784   Median :  439.0   Median :0.9060   Median :  11.00   Median :0.01131  
#Mean   :  1723   Mean   :  610.1   Mean   :0.8949   Mean   :  53.96   Mean   :0.03436  
#3rd Qu.:  1507   3rd Qu.:  680.0   3rd Qu.:0.9307   3rd Qu.:  59.00   3rd Qu.:0.03510  
#Max.   :277960   Max.   :12656.0   Max.   :0.9700   Max.   :9381.00   Max.   :0.73017  

#cells              sample              sex               tissue            nCount_RNA      nFeature_RNA    
#Length:28014       Length:28014       Length:28014       Length:28014       Min.   :   500   Min.   :  100.0  
#Class :character   Class :character   Class :character   Class :character   1st Qu.:   613   1st Qu.:  346.0  
#Mode  :character   Mode  :character   Mode  :character   Mode  :character   Median :   784   Median :  439.0  
#Mean   :  1723   Mean   :  610.1  
#3rd Qu.:  1507   3rd Qu.:  680.0  
#Max.   :277960   Max.   :12656.0  

#nCount_SCT    nFeature_SCT    SCT_snn_res.0.5 seurat_clusters
#Min.   : 545   Min.   :  16.0   0      : 2967   0      : 2967  
#1st Qu.: 645   1st Qu.: 344.0   1      : 2092   1      : 2092  
#Median : 859   Median : 423.0   2      : 1737   2      : 1737  
#Mean   :1028   Mean   : 488.1   3      : 1732   3      : 1732  
#3rd Qu.:1271   3rd Qu.: 587.0   4      : 1583   4      : 1583  
#Max.   :5985   Max.   :1906.0   5      : 1499   5      : 1499  
#(Other):16404   (Other):16404  

# Visuals for different Quality metrics before filtering --------
# Read here: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# and here: https://bioconductor.github.io/CSAMA-labs/single-cell-rnaseq/singlecell_CSAMA2024.html#normalization
# and here: https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/norm.html
# and here: https://rdrr.io/bioc/scuttle/f/vignettes/overview.Rmd

# Visualize the number of cell counts per sample
metadata_gonad_unfilt %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  geom_hline(yintercept = 5000, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  scale_fill_brewer(palette = "PiYG") +  
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(ylim(0,8000))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Unique Cells per Sample")
#ggsave('plots/uniquecells_persample.png', width = 6, height = 8, dpi = 300)

# Visualize the number UMIs/transcripts per cell
metadata_gonad_unfilt %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of UMI counts (transcripts) per cell") +
  ylab("Log10(Cell Density)") +
  geom_vline(xintercept = 500, linetype = "dashed", color = "gray30", linewidth = 0.8)
#ggsave('plots/celldensity_persample.png', width = 6, height = 8, dpi = 300)

# Visualize the distribution of genes detected per cell via histogram
metadata_gonad_unfilt %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of genes detected per cell")
#ggsave('plots/genedensity_percell.png', width = 6, height = 8, dpi = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata_gonad_unfilt %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of genes detected per cell")
#ggsave('plots/genedensity_percell_bar.png', width = 6, height = 8, dpi = 300)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_gonad_unfilt %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)  +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("UMIs vs. Genes Detected  by Mitochondrial Ratio")
#ggsave('plots/umis_vsgenes.png', width = 8, height = 8, dpi = 300)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata_gonad_unfilt %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Mitochondrial gene expression detected per cell")
#ggsave('plots/mitochondrial_countsratio.png', width = 8, height = 8, dpi = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata_gonad_unfilt %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity of gene expression by genes/UMI")
#ggsave('plots/complexity_geneexpression.png', width = 8, height = 8, dpi = 300)


# Cell-level filter (by nUMI, nGene, log10GenesPerUMI, mitoRatio) ----------
filtered_seurat <- subset(x = gonad_umap_3_rename_clusters, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))


metadata_gonad_filter <- filtered_seurat@meta.data # Subset the meta data

summary(metadata_gonad_filter) # Look at the data after filtering

# Output:
#seq_folder             nUMI           nGene         log10GenesPerUMI     mtUMI           mitoRatio       
#Length:26123       Min.   :  500   Min.   :  250.0   Min.   :0.8001   Min.   :   0.00   Min.   :0.000000  
#Class :character   1st Qu.:  613   1st Qu.:  353.0   1st Qu.:0.8803   1st Qu.:   3.00   1st Qu.:0.003339  
#Mode  :character   Median :  778   Median :  452.0   Median :0.9078   Median :   9.00   Median :0.010152  
#Mean   : 1598   Mean   :  619.6   Mean   :0.9022   Mean   :  44.44   Mean   :0.027827  
#3rd Qu.: 1477   3rd Qu.:  694.0   3rd Qu.:0.9330   3rd Qu.:  49.00   3rd Qu.:0.029983  
#Max.   :81176   Max.   :10414.0   Max.   :0.9700   Max.   :1728.00   Max.   :0.199847  

#cells              sample              sex               tissue            nCount_RNA     nFeature_RNA    
#Length:26123       Length:26123       Length:26123       Length:26123       Min.   :  500   Min.   :  250.0  
#Class :character   Class :character   Class :character   Class :character   1st Qu.:  613   1st Qu.:  353.0  
#Mode  :character   Mode  :character   Mode  :character   Mode  :character   Median :  778   Median :  452.0  
#Mean   : 1598   Mean   :  619.6  
#3rd Qu.: 1477   3rd Qu.:  694.0  
#Max.   :81176   Max.   :10414.0  

#nCount_SCT    nFeature_SCT    SCT_snn_res.0.5 seurat_clusters          CellType    
#Min.   : 545   Min.   : 112.0   0      : 2955   0      : 2955   Spermatid    :17839  
#1st Qu.: 637   1st Qu.: 351.0   1      : 2090   1      : 2090   Spermatocytes: 2643  
#Median : 816   Median : 437.0   2      : 1731   2      : 1731   Granulosa    : 1473  
#Mean   : 992   Mean   : 502.8   3      : 1728   3      : 1728   Endocrine    : 1456  
#3rd Qu.:1192   3rd Qu.: 600.0   5      : 1497   5      : 1497   Macrophage   : 1290  
#Max.   :3006   Max.   :1878.0   4      : 1473   4      : 1473   Germ         :  929  
#(Other):14649   (Other):14649   (Other)      :  493  


# Visuals for different Quality metrics after filtering --------

# Visualize the number of cell counts per sample
metadata_gonad_filter %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  geom_hline(yintercept = 5000, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  scale_fill_brewer(palette = "PiYG") +  
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(ylim(0,8000))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Unique Cells per Sample")
ggsave('plots/uniquecells_persample_filtered.png', width = 6, height = 8, dpi = 300)

# Visualize the number UMIs/transcripts per cell
metadata_gonad_filter %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of UMI counts (transcripts) per cell") +
  ylab("Log10(Cell Density)") +
  geom_vline(xintercept = 500, linetype = "dashed", color = "gray30", linewidth = 0.8)
ggsave('plots/celldensity_persample_filtered.png', width = 6, height = 8, dpi = 300)

# Visualize the distribution of genes detected per cell via histogram
metadata_gonad_filter %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of genes detected per cell")
ggsave('plots/genedensity_percell_filtered.png', width = 6, height = 8, dpi = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata_gonad_filter %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of genes detected per cell")
ggsave('plots/genedensity_percell_bar_filtered.png', width = 6, height = 8, dpi = 300)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_gonad_filter %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)  +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("UMIs vs. Genes Detected  by Mitochondrial Ratio")
ggsave('plots/umis_vsgenes_filtered.png', width = 8, height = 8, dpi = 300)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata_gonad_filter %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Mitochondrial gene expression detected per cell")
ggsave('plots/mitochondrial_countsratio_filtered.png', width = 8, height = 8, dpi = 300)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata_gonad_filter %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity of gene expression by genes/UMI")
ggsave('plots/complexity_geneexpression_filtered.png', width = 8, height = 8, dpi = 300)


# Subset and organize filtered data --------

counts_gonad_filter <- filtered_seurat@assays$RNA@counts # Subset the count data

metadata_gonad_filter$CellType <- factor(filtered_seurat@active.ident) # Add info for cell type

gonad_subset <- subset(filtered_seurat, idents="Germ", invert=T) # Removes Germ cell types

counts_gonad <- gonad_subset@assays$RNA@counts # Subset the count data

metadata_gonad <- gonad_subset@meta.data # Subset the meta data

metadata_gonad$CellType <- factor(gonad_subset@active.ident) # Add info for cell type

sce_gonad <- SingleCellExperiment(assays=list(counts=counts_gonad), colData=metadata_gonad) # Put together counts and meta data into a named list
sce_gonad 

# Keep only male data
sce_gonad_male <- sce_gonad[, sce_gonad$sex == "Male"]

# Check that I only have male data
table(sce_gonad_male$sex) 

# Converting some metadata to factors for grouping, and putting into new column
# Comparing male-only data to female and male data 
colData(sce_gonad_male)$sample_id <- as.factor(colData(sce_gonad_male)$sample)
summary(sce_gonad_male$sample_id) # Only male data

colData(sce_gonad)$sample_id <- as.factor(colData(sce_gonad)$sample)
summary(sce_gonad$sample_id) # Data for both sexes


colData(sce_gonad_male)$sex_id <- as.factor(colData(sce_gonad_male)$sex)
summary(sce_gonad_male$sex_id) # Only male data

colData(sce_gonad)$sex_id <- as.factor(colData(sce_gonad)$sex)
summary(sce_gonad$sex_id) # Data for both sexes


colData(sce_gonad_male)$cluster_id <- as.factor(colData(sce_gonad_male)$CellType)
summary(sce_gonad_male$cluster_id) # Only male data, there should be less female-specific cells (Granulosa), but retain almost all Spermatid and Spermatocytes

colData(sce_gonad)$cluster_id <- as.factor(colData(sce_gonad)$CellType)
summary(sce_gonad$cluster_id) # Data for both sexes

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
kids_gonad <- purrr::set_names(levels(sce_gonad_male$cluster_id)) # List all possible cluster/ cell-type names
kids_gonad

nk_gonad <- length(kids_gonad) # How many clusters/ cell-types are there?
nk_gonad
# 7

sids_gonad <- purrr::set_names(levels(sce_gonad_male$sample_id)) # List number of samples - 3 M only
sids_gonad
 
ns_gonad <- length(sids_gonad) # How many samples are there?
ns_gonad
# 3

# Number of cells per sample
table(sce_gonad_male$sample_id)
# GonadM1 GonadM2 GonadM3 
#  7946    6468    6738  

# Turn named vector into a numeric vector
n_cells_gonad <- as.numeric(table(sce_gonad_male$sample_id))

# Reorder samples (rows) of the metadata to match the order of the sample names
m_gonad <- match(sids_gonad, sce_gonad_male$sample_id)
m_gonad
#     1  7947 14415

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_gonad <- data.frame(colData(sce_gonad_male)[m_gonad, ], n_cells_gonad, row.names=NULL) %>% select(-"cluster_id")

# Remove genes with low expression 
dim(sce_gonad_male) # Before filtering for low expression
# 21793 21152

sce_gonad_male <- sce_gonad_male[rowSums(counts(sce_gonad_male)>0)>0,] # This removes genes that have zero expression in all cells
dim(sce_gonad_male)	
# 19349 21152

sce_gonad_male <- sce_gonad_male[rowSums(counts(sce_gonad_male)>1)>=10,] # This removes genes that are not expressed (above 1 count) in at least 10 cells, filtering out low-expressed genes
dim(sce_gonad_male)
# 8358 21152


# Aggregates gene counts (pb_gonad) by summing across cells within each cluster and sample to create a pseudobulk expression matrix
groups_gonad <- colData(sce_gonad_male)[, c("cluster_id", "sample_id")]
pb_gonad <- aggregate.Matrix(t(counts(sce_gonad_male)), groupings=groups_gonad, fun="sum")

dim(pb_gonad)
#  21 8358

pb_gonad[1:21, 1:6]
#21 x 6 sparse Matrix of class "dgCMatrix"
#                           arl13b zgc:152904 timmdc1 sap18 nid2a rtraf
#Spermatid_GonadM4           69        417      26   863     9   452
#Spermatocytes_GonadM4       15         32      17   124     .    93
#Granulosa_GonadM4            1          2       1     8     .     2
#Macrophage_GonadM4           1          .       .     2     .     .
#Endocrine_GonadM4            2          9       1    23     .     4
#Monocyte_GonadM4             .          .       .     .     .     .
# ...

## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_gonad <- sapply(stringr::str_split(rownames(pb_gonad), pattern="_", n=2), '[', 1) 

splitf_gonad # Check that each cluster is present in all samples

pb_gonad <- split.data.frame(pb_gonad, factor(splitf_gonad)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

str(pb_gonad)

# Number of cells for each cell type and sample
options(width=100)
table(sce_gonad_male$cluster_id, sce_gonad_male$sample_id)

# GonadM4 GonadM5 GonadM6
#Spermatid          7310    4823    5630
#Spermatocytes       311    1347     984
#Granulosa            89      11      14
#Macrophage           18      38      13
#Endocrine           147     232      80
#Monocyte             11       2       4
#Hemoglobin rich      60      15      13



# Function to get sample IDs (which are the column names)
get_sample_ids_gonad <- function(x){
  pb_gonad[[x]] %>% colnames()
}

# Use function to get sample IDs for each cell type/ cluster ID, convert from list to vector
de_samples_gonad <- map(1:length(kids_gonad), get_sample_ids_gonad) %>% unlist()

# Use function to get sample IDs for each cell type/ cluster ID, keep as list
samples_list_gonad <- map(1:length(kids_gonad), get_sample_ids_gonad)

# Function to get cell type/cluster IDs (which are the column names)
get_cluster_ids_gonad <- function(x){
  rep(names(pb_gonad)[x], each=length(samples_list_gonad[[x]]))
}

# Use function to get cluster IDs for each sample, convert from list to vector
de_cluster_ids_gonad <- map(1:length(kids_gonad), get_cluster_ids_gonad) %>% unlist()

# Construct metadata into a dataframe, where we have both cluster ID and sample ID
gg_df_gonad <- data.frame(cluster_id=de_cluster_ids_gonad, sample_id=de_samples_gonad)

# Add in info from the ei_gonad dataset (including sample_id, sex, and number of cells), joining with `by = join_by(sample_id)`
gg_df_gonad <- left_join(gg_df_gonad, ei_gonad[, c("sample_id", "sex_id", "n_cells_gonad")])

# Organize data by keeping only relevant columns (if messy)
metadata_gonad <- gg_df_gonad %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_gonad)

# Convert the cluster_id to a factor (it is categorical data)
metadata_gonad$cluster_id <- as.factor(metadata_gonad$cluster_id)

# Extract the unique cell types/ cluster names as a vector
clusters_gonad <- levels(metadata_gonad$cluster_id)

# Normalization method using DESeq -------------
# This method uses a median-of-ratios (size factors) to normalize, it produces DESeq2 normalized counts. This 
# is robust to outliers, and is better with sparse data (but we removed lowly expressed genes).

# Function to normalize pseudobulk matrix using DESeq2
deseq_normalize_cluster <- function(cluster_name) {
  message("Normalizing: ", cluster_name)
  
  # Extract pseudobulk count matrix for the cluster
  counts <- as.matrix(pb_gonad[[cluster_name]])
  
  # Get matching metadata for these samples
  sample_ids <- colnames(counts)
  cluster_metadata <- metadata_gonad %>%
    filter(cluster_id == cluster_name, sample_id %in% sample_ids)
  
  # Ensure rownames of colData match colnames of count matrix
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  cluster_metadata <- cluster_metadata[sample_ids, , drop = FALSE]  # Match order
  
  # Create DESeq2 object with a dummy design since we won’t run DE
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = cluster_metadata,
                                design = ~1)  # no condition modeling
  
  # Normalize using size factor estimation
  dds <- estimateSizeFactors(dds)
  
  # Extract normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Save (optional)
  write.csv(norm_counts, paste0("DESeq2_normalized_", cluster_name, ".csv"), quote = FALSE)
  
  return(norm_counts)
}

# Apply to all clusters
deseq_norm_list <- map(clusters_gonad, deseq_normalize_cluster)
names(deseq_norm_list) <- clusters_gonad


# Normalize the pseudobulk data using TMM-normalization in edgeR (generates TMM-normalized CPM values) ---------
# This method will correct for composition bias (for example if one highly expressed gene dominates total counts)
# and library size differences and it is scaling the counts per million. 

# Function to TMM-normalize a single cluster's pseudobulk matrix
tmm_normalize_cluster <- function(cluster_name) {
  # Extract count matrix for that cluster
  counts <- as.matrix(pb_gonad[[cluster_name]])
  
  # Create DGEList object
  dge <- DGEList(counts = counts)
  
  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Compute CPM (normalized.lib.sizes = TRUE ensures TMM scaling is used)
  cpm_tmm <- cpm(dge, normalized.lib.sizes = TRUE)
  
  # Save to CSV (optional)
  write.csv(cpm_tmm, paste0("TMM_CPM_", cluster_name, ".csv"), quote = FALSE)
  
  # Return the normalized matrix (optional)
  return(cpm_tmm)
}

# Apply the function to all clusters in pb_gonad
tmm_cpm_list <- map(clusters_gonad, tmm_normalize_cluster)
names(tmm_cpm_list) <- clusters_gonad

# Visualize the distribution of the normalized counts --------------

# Load the data (either CPM values from edgeR or normalized counts from DESeq)
normdata <- read_csv("TMM_CPM_Spermatid.csv")

normdata_long <- normdata |>
  pivot_longer(cols = -1, names_to = "sample", values_to = "CPM") |>
  mutate(log_CPM = log1p(CPM))  # log1p for the zeros 

# Histogram of log-transformed (distribution of data)
ggplot(normdata_long, aes(x = log_CPM)) +
  geom_histogram(bins = 50, fill = "#277da1", color = "white", alpha = 0.8) +
  labs(
    title = "Log-transformed CPM Distribution",
    x = "log1p(CPM)",  # log1p(x) = log(x + 1)
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)

# QQ plots (normality of data)
ggplot(normdata_long, aes(sample = log_CPM)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~sample) +
  theme_minimal(base_size = 14) +
  labs(title = "QQ Plot of log(CPM) per Sample")

# Histogram of each sample (distribution of data)
ggplot(normdata_long, aes(x = log_CPM, fill = sample)) +
  geom_density(alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Density of log(CPM) Across Samples",
    x = "log1p(CPM)",
    y = "Density"
  )

