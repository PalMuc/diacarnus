library(ggplot2)
library(dplyr)



eigenvectors <- read.table("~/Downloads/Diacarnus_pca_mind.eigenvec", quote="\"", comment.char="")
eigenvalues <- read.table("~/Downloads/Diacarnus_pca_mind.eigenval", quote="\"", comment.char="")
pop_info <- read.csv("~/Downloads/pop_info.txt", sep="")
# after your read.table() call
colnames(eigenvectors)[1:2] <- c("FID", "IID")

n_pcs <- ncol(eigenvectors) - 2
pc_names <- paste0("PC", seq_len(n_pcs))
colnames(eigenvectors)[3:ncol(eigenvectors)] <- pc_names


eigenvectors <- eigenvectors %>%
  left_join(pop_info, by = "IID")

eigenvectors$pop <- factor(eigenvectors$pop)

ggplot(eigenvectors, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 0.7) +  # use linewidth, not size
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  theme_minimal()
