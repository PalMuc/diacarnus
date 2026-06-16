setwd("/home/lvilleg/Diacarnus_0526/04_inStrain_46MAGS/MAG_compare_new/MAGs_compare_clustered_new/output")

cmp <- read.table("MAGs_compare_clustered_new_genomeWide_compare.tsv",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

strip_host <- function(x) sub("\\.MAGs.*$", "", x)

cmp$host1 <- strip_host(cmp$name1)
cmp$host2 <- strip_host(cmp$name2)

all_mags <- sort(unique(cmp$genome))
length(all_mags)  # just to see how many

for (mag_id in all_mags) {
  cmp_mag <- cmp[cmp$genome == mag_id, ]
  
  hosts <- sort(unique(c(cmp_mag$host1, cmp_mag$host2)))
  n <- length(hosts)
  if (n < 2) next  # skip weird cases
  
  mag_dist <- matrix(NA_real_, nrow = n, ncol = n,
                     dimnames = list(hosts, hosts))
  
  for (i in seq_len(nrow(cmp_mag))) {
    h1 <- cmp_mag$host1[i]
    h2 <- cmp_mag$host2[i]
    d  <- 1 - cmp_mag$popANI[i]
    mag_dist[h1, h2] <- d
    mag_dist[h2, h1] <- d
  }
  diag(mag_dist) <- 0
  
  safe_id <- gsub("[^A-Za-z0-9_.-]", "_", mag_id)
  out_file <- paste0("MAG_dist_", safe_id, ".tsv")
  write.table(mag_dist,
              file = out_file,
              sep = "\t", quote = FALSE, col.names = NA)
  cat("Wrote:", out_file, "\n")
}
