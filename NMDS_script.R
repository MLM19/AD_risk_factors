# -- INSTALL / LOAD PACKAGES --------------------------------------------
# Install once if you haven't already:
# install.packages(c("smacof", "ggplot2"))

library(smacof)    # smacofSym() for nonmetric MDS
library(ggplot2)   # for plotting
library(ggrepel)

# -- 1) LOAD YOUR MATRIX -----------------------------------------------
# CSV assumed: first column = rownames, header = trait names
mat <- as.matrix(
  read.delim(
    "significant_snp_matrix.tsv",
    header      = TRUE,
    row.names   = 1,
    check.names = FALSE
  )
)

# -- 2) HANDLE ZERO COUNTS ---------------------------------------------
# Choose a small epsilon (e.g. 1e-3) to replace zeros
eps <- 1e-3
mat[mat == 0] <- eps

# -- 3) CREATE DISTANCE MATRIX -----------------------------------------
# Invert counts -> larger counts = smaller distances
dist_mat <- 1 / mat

# Convert to 'dist' object for smacof
diss <- as.dist(dist_mat)

# -- 4) RUN NON‑METRIC MDS ----------------------------------------------
# ndim = 2 dimensions; you can increase trymax if convergence is an issue
nmds <- smacofSym(diss, ndim = 2, verbose = TRUE)

# Extract 2D coordinates
coords <- as.data.frame(nmds$conf)
colnames(coords) <- c("Dim1", "Dim2")
coords$Trait <- rownames(coords)
coords$IsAD <- coords$Trait == "AD_7158"

# -- 5) PLOT WITH ggplot2 ------------------------------------------------
# Highlight the "AD" trait in red; adjust trait name if different.
p <- ggplot(coords, aes(x = Dim1, y = Dim2)) +
  # all points, gray; AD in red
  geom_point(aes(color = IsAD), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick"),
                     guide = FALSE) +
  # non‑overlapping labels
  geom_text_repel(aes(label = Trait),
                  size = 3,
                  max.overlaps = 30,
                  segment.size = 0.2) +
  # emphasize AD label in bold red
  geom_text_repel(data = subset(coords, IsAD),
                  aes(label = Trait),
                  size = 4,
                  fontface = "bold",
                  color = "firebrick",
                  nudge_y = 0.1,
                  segment.color = "firebrick") +
  # clean
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title     = element_text(face = "bold"),
    plot.subtitle  = element_text(face = "italic")
  ) +
  labs(
    title    = "Non metric MDS of Trait Similarities",
    x        = "MDS Dimension 1",
    y        = "MDS Dimension 2"
  )

print(p + theme_bw())
# -- OPTIONAL: SAVE PLOT ------------------------------------------------
ggsave("nmds_smacof_plot_AD.png", width = 6, height = 5, dpi = 300)
