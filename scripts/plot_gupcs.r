library(tidyverse)
setwd("~/Desktop/ebame3/data/")
# Let's read the files
m_vs_u <- read_tsv("resultDB80_tophits.m8.gz", col_names = FALSE, progress = TRUE) %>%
  select(X1, X2, X3, X11) %>%
  rename(gene_callers_id = X1, unk_id = X2, id = X3, evalue = X11)
m_genes <- read_tsv("orf2mag.tsv.gz", col_names = TRUE, progress = TRUE) %>%
  select(gene_callers_id, contig)
mag_cdata <- read_tsv("TARA_MAGs_v3_metadata.txt", col_names = TRUE)

# Some data massaging
# Combine tables and clean MAG names
m_vs_u_mag <- m_vs_u %>%
  left_join(m_genes) %>%
  separate(unk_id, into = "uclass", extra = "drop", remove = FALSE) %>%
  tidyr::extract(contig, into = paste("V", 1:4, sep = ""), regex = "([[:alnum:]]+)_([[:alnum:]]+)_([[:alnum:]]+)_([[:alnum:]]+)") %>%
  unite(col = MAG, V1,V2,V3,V4, sep = "_", remove = TRUE)

mag_n_orfs <- m_genes %>%
  tidyr::extract(contig, into = paste("V", 1:4, sep = ""), regex = "([[:alnum:]]+)_([[:alnum:]]+)_([[:alnum:]]+)_([[:alnum:]]+)") %>%
  unite(col = MAG, V1,V2,V3,V4, sep = "_", remove = TRUE) %>%
  group_by(MAG) %>%
  count() %>%
  rename(n_orf = n)

# Count proportion of ORFs in each MAG
m_vs_u_mag_n <- m_vs_u_mag %>%
  group_by(uclass, MAG) %>%
  count() %>%
  ungroup() %>%
  tidyr::complete(uclass, MAG, fill = list(n = 0)) %>%
  left_join(mag_n_orfs) %>%
  mutate(prop = n/n_orf)

# Get the 25 more abundants
m_vs_u_mag_n_25 <- m_vs_u_mag_n %>%
  filter(uclass == "eupc") %>%
  dplyr::top_n(25) %>%
  .$MAG

# Plot them
ggplot(m_vs_u_mag_n %>% filter(MAG %in% m_vs_u_mag_n_25) , aes(MAG, n, fill = uclass)) +
  geom_bar(stat = "identity") +
  theme_light() +
  xlab("") +
  ylab("# hits") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate histograms
all <- ggplot(m_vs_u_mag_n, aes(prop)) +
  geom_histogram(alpha = 0.5, fill = "#2E3239") +
  theme_light() +
  ylab("Frequency") +
  xlab("Proportion of unknown ORFs") +
  scale_x_continuous(labels = scales::percent)

eupc <- ggplot(m_vs_u_mag_n %>% filter(uclass == "eupc"), aes(prop)) +
  geom_histogram(alpha = 0.5, fill = "#2E3239") +
  theme_light() +
  ylab("Frequency") +
  xlab("Proportion of env. unknowns ORFs") +
  scale_x_continuous(labels = scales::percent)

gupc <- ggplot(m_vs_u_mag_n %>% filter(uclass == "gupc"), aes(prop)) +
  geom_histogram(alpha = 0.5, fill = "#2E3239") +
  theme_light() +
  ylab("Frequency") +
  xlab("Proportion of gen. unknowns ORFs") +
  scale_x_continuous(labels = scales::percent)

# Plot them
ggpubr::ggarrange(all, gupc, eupc, nrow = 1, ncol = 3)

m_vs_u_mag %>%
  group_by(MAG, uclass) %>%
  count %>% View

