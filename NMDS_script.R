#non metric multidimensional scaling projection

#--LOAD PACKAGES --------------------------------------------
#Install once if you haven't already:
#install.packages(c("smacof", "ggplot2", "ggrepel"))

library(smacof)    #smacofSym() for nonmetric MDS
library(ggplot2)   #for plotting the projection
library(ggrepel)   #to avoid label overlap

#-- 1) LOAD MATRIX -----------------------------------------------
#first column = rownames, header = trait names
mat <- as.matrix(
  read.delim(
    "rescaled_significant_snp_matrix.tsv",
    header      = TRUE,
    row.names   = 1,
    check.names = FALSE
  )
)

#-- 2)Transformations for using distance matrix (dissimilarity) ------------
##Replace 0 with small number to avoid errors in transformations
eps <- 1e-3
mat[mat == 0] <- eps

#-- 2.1) CREATE DISTANCE MATRIX -----------------------------------------
##Invert counts -> larger counts = smaller distances
dist_mat <- 1 / mat

##Convert to 'dist' object for smacof
diss <- as.dist(dist_mat)

##It can be done using:
#diss <- sim2diss(mat, method = "reciprocal", to.dist = TRUE)

# -- 4) RUN NON‑METRIC MDS ----------------------------------------------
#Default: direct ratio scaling
#nmds <- smacofSym(diss, ndim = 2, verbose = TRUE)

#Ordinal MDS
#nmds <- smacofSym(diss, ndim = 2, type = "ordinal")

#ties? how the ties are handled: primary, secondary, tertiary
set.seed(42)
nmds <- smacofSym(diss, ndim = 2, type = "ordinal", ties = "tertiary")
#nmds <- smacofSym(diss, ndim = 2, ties = "primary")

#print stress to know how well the rep fits
message("NMDS stress: ", nmds$stress)

# Extract 2D coordinates
coords <- as.data.frame(nmds$conf)
colnames(coords) <- c("Dim1", "Dim2")
coords$Trait <- rownames(coords)
write.table(
  coords,
  file      = "nmds_coords.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
coords$IsAD <- coords$Trait == "AD_7158"
coords$Label <- sub("_\\d+(\\.\\d+)?$", "", coords$Trait)

# -- 5) PLOT WITH ggplot2 ------------------------------------------------
# Highlight the "AD" trait in red; adjust trait name if different.
p <- ggplot(coords, aes(x = Dim1, y = Dim2)) +
  #all points, gray; AD in red
  geom_point(aes(color = IsAD), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick"),
                     guide = FALSE) +
  # Only label non-AD traits here
  geom_text_repel(
    data = subset(coords, !IsAD),
    aes(label = Label),
    size = 3,
    max.overlaps = 30,
    segment.size = 0.2
  ) +
  # Label AD separately and in red/bold
  geom_text_repel(
    data = subset(coords, IsAD),
    aes(label = Label),
    size = 4,
    fontface = "bold",
    color = "firebrick",
    nudge_y = 0.1,
    segment.color = "firebrick"
  ) +
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
ggsave("nmds_ordinal_tertiary_AD.png", width = 6, height = 5, dpi = 300)
 


##Filtering out otliars (Change manually)
excluded_traits <- c("HIV_5531", "BRE_CANCER_174", "EBS_20002", "ASBEST_3282", "ALCOH_3727") 

coords_filtered <- subset(coords, !(Trait %in% excluded_traits))
coords_filtered$IsAD <- coords_filtered$Trait == "AD_7158"
coords_filtered$Label <- sub("_\\d+(\\.\\d+)?$", "", coords_filtered$Trait)

p <- ggplot(coords_filtered, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = IsAD), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick"), guide = FALSE) +
  geom_text_repel(
    data = subset(coords_filtered, !IsAD),
    aes(label = Label),
    size = 3,
    max.overlaps = 30,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(coords_filtered, IsAD),
    aes(label = Label),
    size = 4,
    fontface = "bold",
    color = "firebrick",
    nudge_y = 0.1,
    segment.color = "firebrick"
  ) +
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
ggsave("nmds_filtered2_ordinal_tertiary_AD.png", width = 6, height = 5, dpi = 300)
