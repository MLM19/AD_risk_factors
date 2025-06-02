# nmds_and_clustering.R

# === Load required packages ===
library(smacof)
library(ggplot2)
library(ggrepel)
library(mclust)
library(dbscan)
library(readr)

# === Load rescaled matrix ===
mat <- as.matrix(read.delim(
  "rescaled_significant_snp_matrix.tsv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
))

# === Replace zeros with small value to avoid errors ===
eps <- 1e-3
mat[mat == 0] <- eps

# === Create dissimilarity matrix ===
dist_mat <- 1 / mat
diss <- as.dist(dist_mat)

# === NMDS Projection ===
set.seed(42)
nmds <- smacofSym(diss, ndim = 2, type = "ordinal", ties = "tertiary")

coords <- as.data.frame(nmds$conf)
colnames(coords) <- c("Dim1", "Dim2")
coords$Trait <- rownames(coords)
coords$Label <- sub("_\\d+$", "", coords$Trait)
coords$IsAD <- coords$Trait == "AD_7158"

# === Save NMDS coordinates ===
write.table(coords, file = "nmds_coords.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# === Plot NMDS ===
p <- ggplot(coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = IsAD), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick"), guide = FALSE) +
  geom_text_repel(data = subset(coords, !IsAD), aes(label = Label), size = 3, segment.size = 0.2) +
  geom_text_repel(data = subset(coords, IsAD), aes(label = Label), size = 4, fontface = "bold",
                  color = "firebrick", nudge_y = 0.1, segment.color = "firebrick") +
  theme_bw(base_size = 14) +
  labs(title = "Non-metric MDS of Trait Similarities",
       x = "MDS Dimension 1", y = "MDS Dimension 2")

ggsave("nmds_projection.png", plot = p, width = 6, height = 5, dpi = 300)

# === Clustering ===
# MCLUST
mclust_res <- Mclust(coords[, c("Dim1", "Dim2")])
coords$MclustCluster <- factor(mclust_res$classification)

# DBSCAN
db_res <- dbscan(coords[, c("Dim1", "Dim2")], eps = 0.2, minPts = 4)
coords$DBSCANcluster <- factor(db_res$cluster)

# Save combined data
write_tsv(coords, "nmds_with_clusters.tsv")

# MCLUST plot
p_mclust <- ggplot(coords, aes(x = Dim1, y = Dim2, color = MclustCluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3) +
  labs(title = "MCLUST Clustering", x = "Dim1", y = "Dim2") +
  theme_bw()

ggsave("mclust_clustering.png", plot = p_mclust, width = 6, height = 5, dpi = 300)

# DBSCAN plot
p_dbscan <- ggplot(coords, aes(x = Dim1, y = Dim2, color = DBSCANcluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3) +
  labs(title = "DBSCAN Clustering", x = "Dim1", y = "Dim2") +
  theme_bw()

ggsave("dbscan_clustering.png", plot = p_dbscan, width = 6, height = 5, dpi = 300)
