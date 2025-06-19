# R script for clustering NMDS projections of AD risk factor traits
# using MCLUST (model-based clustering) and DBSCAN (density-based clustering)

# 1. Load required libraries
if (!require("mclust")) install.packages("mclust", dependencies = TRUE)
if (!require("dbscan")) install.packages("dbscan", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("readr")) install.packages("readr", dependencies = TRUE)

library(mclust)
library(dbscan)
library(ggplot2)
library(readr)

# 2. Read in data
# - nmds_coords.tsv: NMDS coordinates (43 points x columns: Dim1, Dim2, Trait, IsAD)
# - significant_snp_matrix.tsv is available if you want to overlay later
nmds_coords <- read_tsv("nmds_coords.tsv", col_types = cols())

#filter traits
excluded_traits <- c("HIV_5531", "BRE_CANCER_174", "EBS_20002") 
nmds_coords <- subset(nmds_coords, !(Trait %in% excluded_traits))

# 3. Prepare data
coords <- nmds_coords[, c("Dim1", "Dim2")]
rownames(coords) <- nmds_coords$Trait
# Strip trailing underscore+digits from trait names for labels
nmds_coords$Label <- sub("_\\d+(\\.\\d+)?$", "", nmds_coords$Trait)
# Identify AD trait index
ad_trait <- "AD_7158"
ad_idx <- which(nmds_coords$Trait == ad_trait)

# 4. Model-based clustering with MCLUST
set.seed(42) 
mclust_res <- Mclust(coords)
nmds_coords$MclustCluster <- factor(mclust_res$classification)

# 5. Density-based clustering with DBSCAN
kNNdistplot(coords, k = 4)
abline(h = 0.2, lty = 2)  # adjust eps based on your elbow

db_res <- dbscan(coords, eps = 0.2, minPts = 4)
nmds_coords$DBSCANcluster <- factor(db_res$cluster)

# 6. Plotting function: discrete colors, AD highlighted, size-differentiated labels
plot_clusters <- function(df, cluster_col, title) {
  ggplot(df, aes_string(x = "Dim1", y = "Dim2", color = cluster_col)) +
    scale_color_brewer(palette = "Set2", na.value = "grey50", name = "Cluster", drop= TRUE, na.translate= FALSE) +
    geom_point(size = 4) +
    # AD outline
    geom_point(data = df[ad_idx, ], aes(x = Dim1, y = Dim2), 
               color = "darkblue", size = 8, shape = 1, stroke = 1.5, inherit.aes = FALSE, show.legend = FALSE) +
    # all other labels
    geom_text(data = df[-ad_idx, ], aes(label = Label), 
              size = 3, fontface = "plain", hjust = 0.5, vjust = -0.7, show.legend = FALSE) +
    # AD label bold and larger
    geom_text(data = df[ad_idx, ], aes(label = Label), 
              size = 5, fontface = "bold", hjust = 0.5, vjust = -0.7, show.legend = FALSE, colour = "darkblue") +
    labs(title = title) +
    theme_bw() +
    theme(legend.position = "right")
}

# 7. Generate plots
p_mclust <- plot_clusters(nmds_coords, "MclustCluster", "MCLUST Clustering of NMDS Traits")
print(p_mclust)

p_dbscan <- plot_clusters(nmds_coords, "DBSCANcluster", "DBSCAN Clustering of NMDS Traits")
print(p_dbscan)

# 8. Save results
write_tsv(nmds_coords, "nmds_with_clusters.tsv")
ggsave("mclust_plt.png", p_mclust, width = 6, height = 5, dpi = 300)
ggsave("dbscan_plt.png", p_dbscan, width = 6, height = 5, dpi = 300)

# Notes:
# - Trait labels now omit the numeric suffix (e.g. "DIAB_TYP1" not "DIAB_TYP1_3791").
# - All labels are size 3, but AD is size 5 and bold for clear emphasis.
# - Clusters are colored discretely using a qualitative Set2 palette.
# - Outline of AD is a large red hollow circle. Feel free to adjust sizes or colors.
# - Tune DBSCAN parameters (eps/minPts) by inspecting the kNN distance plot.


# 9. Filter for the cluster where AD is found, redo the projection and the clustering
# -- Load Required Packages ------------------------------------------------
library(smacof)
library(ggplot2)
library(ggrepel)
library(readr)
library(mclust)
library(dbscan)

# -- Load Full NMDS with Clusters -----------------------------------------
nmds_full <- read_tsv("nmds_with_clusters.tsv", col_types = cols())

# Filter for Cluster 1 Traits
cluster1_traits <- nmds_full$Trait[nmds_full$MclustCluster == 1]

# Load and Subset SNP Overlap Matrix --------------------------------------
mat <- as.matrix(read.delim("rescaled_significant_snp_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE))
mat_cluster1 <- mat[cluster1_traits, cluster1_traits]

# Avoid infinite distances by replacing zeros
eps <- 1e-3
mat_cluster1[mat_cluster1 == 0] <- eps
dist_mat <- 1 / mat_cluster1
diss <- as.dist(dist_mat)

# -- Run NMDS on Filtered Data --------------------------------------------
set.seed(42) 
nmds_cluster1 <- smacofSym(diss, ndim = 2, type = "ordinal", ties = "tertiary")
coords <- as.data.frame(nmds_cluster1$conf)
colnames(coords) <- c("Dim1", "Dim2")
coords$Trait <- rownames(coords)
coords$Label <- sub("_\\d+(\\.\\d+)?$", "", coords$Trait)
coords$IsAD <- coords$Trait == "AD_7158"

# -- Plot 1: Filtered NMDS Projection -------------------------------------
p1 <- ggplot(coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = IsAD), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick"), guide = FALSE) +
  geom_text_repel(
    data = subset(coords, !IsAD),
    aes(label = Label),
    size = 3, max.overlaps = 30, segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(coords, IsAD),
    aes(label = Label),
    size = 4, fontface = "bold", color = "firebrick",
    nudge_y = 0.1, segment.color = "firebrick"
  )+
  theme_bw(base_size = 14) +
  labs(
    title = "NMDS Projection of Cluster 1 Traits",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  )

ggsave("nmds_cluster1_Mclust_projection.png", p1, width = 6, height = 5, dpi = 300)

# -- Run MCLUST on Filtered NMDS ------------------------------------------
set.seed(42) 
mclust_res <- Mclust(coords[, c("Dim1", "Dim2")])
coords$MclustCluster <- factor(mclust_res$classification)

# -- Plot 2: MCLUST Clustering on Filtered NMDS ---------------------------
p2 <- ggplot(coords, aes(x = Dim1, y = Dim2, color = MclustCluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    data = subset(coords, !IsAD),
    aes(label = Label),
    size = 3, max.overlaps = 30, segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(coords, IsAD),
    aes(label = Label),
    size = 4, fontface = "bold", color = "firebrick",
    nudge_y = 0.1, segment.color = "firebrick"
  )+
  theme_bw(base_size = 14) +
  labs(
    title = "MCLUST Clustering of Cluster 1 Traits",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2",
    color = "Cluster"
  )

ggsave("nmds_cluster1_Mclust_mclust.png", p2, width = 6, height = 5, dpi = 300)

# -- Run DBSCAN on Filtered NMDS ------------------------------------------
# Estimate eps from kNN distance (optional tuning)
# kNNdistplot(coords[, c("Dim1", "Dim2")], k = 3); abline(h = 0.5, lty = 2)
db_res <- dbscan(coords[, c("Dim1", "Dim2")], eps = 0.5, minPts = 3)
coords$DBSCAN <- factor(ifelse(db_res$cluster == 0, "Noise", paste0("C", db_res$cluster)))

# -- Plot 3: DBSCAN Clustering on Filtered NMDS ---------------------------
p3 <- ggplot(coords, aes(x = Dim1, y = Dim2, color = DBSCAN)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    data = subset(coords, !IsAD),
    aes(label = Label),
    size = 3, max.overlaps = 30, segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(coords, IsAD),
    aes(label = Label),
    size = 4, fontface = "bold", color = "firebrick",
    nudge_y = 0.1, segment.color = "firebrick"
  )+
  theme_bw(base_size = 14) +
  labs(
    title = "DBSCAN Clustering of Cluster 1 Traits",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2",
    color = "Cluster"
  )

ggsave("nmds_cluster1_Mclust_dbscan.png", p3, width = 6, height = 5, dpi = 300)

