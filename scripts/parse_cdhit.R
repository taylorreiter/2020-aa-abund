library(readr)
library(dplyr)
library(tidyr)
# "outputs/cd-hit95/hu-genome19.cdhit95.faa.clstr"
clstr <- read_tsv(snakemake@input[["cluster"]],
                  col_names = c("cluster", "info"))

# If a row contains the word "Cluster", keep it's name the same.
# Otherwise, inherit the name from the previous row.
for(i in 1:nrow(clstr)){
  clstr$cluster[i] <- ifelse(grepl("Cluster", clstr$cluster[i]), 
                             clstr$cluster[i], clstr$cluster[i-1])
}

# remove rows that have no information associated with them
clstr <- clstr %>%
  filter(!is.na(info)) %>%
  mutate(cluster = gsub(">Cluster ", "", cluster)) %>%
  separate(info, into = c("length", "name", "tmp", "identity"), sep = " ") %>%
  mutate(length = gsub("aa,", "", length)) %>%
  mutate(name = gsub(">", "", name)) %>%
  select(cluster, length, name, identity)

clstr %>%
  group_by(cluster) %>%
  select(name) %>%
  group_walk(~ write_csv(.x, paste0(snakemake@output[[1]], .y$cluster, "_cluster.csv"), colnames = F))
