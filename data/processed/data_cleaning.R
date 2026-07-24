library(dplyr)

T_len <- 10; K <- 11
prev_cut <- 0.30; species_only <- TRUE

# read data
meta <- read.csv("hmp2_metadata_2018-08-20.csv", check.names = FALSE)
meta <- meta %>%filter(data_type == "metagenomics")

tax <- read.delim("taxonomic_profiles_3.tsv.gz", check.names = FALSE)
rownames(tax) <- tax[[1]]
tax <- tax[, -1]

tax <- as.data.frame(t(tax), stringsAsFactors = FALSE)

# clean sample ids
tax$sample_id_clean <- sub("_P$", "", sub("_profile$", "", rownames(tax)))
meta$sample_id_clean <- sub("_P$", "", meta$`External ID`)

# keep one row per sample
tax <- tax[!duplicated(tax$sample_id_clean), ]
meta <- meta[!duplicated(meta$sample_id_clean), ]

common_ids <- intersect(tax$sample_id_clean, meta$sample_id_clean)

tax <- tax[match(common_ids, tax$sample_id_clean), ]
meta <- meta[match(common_ids, meta$sample_id_clean), ]

stopifnot(all(tax$sample_id_clean == meta$sample_id_clean))

# choose subjects with at least 10 distinct weeks
subj_info <- meta %>%group_by(`Participant ID`) %>%
               summarise(n_sample = n(),
                        n_week = n_distinct(week_num),
                        .groups = "drop")

keep_subj <- subj_info$`Participant ID`[subj_info$n_week >= T_len]

meta_T <- meta %>%
  filter(`Participant ID` %in% keep_subj) %>%
  arrange(`Participant ID`, week_num) %>%
  group_by(`Participant ID`) %>%
  distinct(week_num, .keep_all = TRUE) %>%
  slice(1:T_len) %>%
  ungroup()

tax_T <- tax[match(meta_T$sample_id_clean, tax$sample_id_clean), ]
stopifnot(all(tax_T$sample_id_clean == meta_T$sample_id_clean))

# taxon table
tax_mat <- tax_T[, setdiff(colnames(tax_T), "sample_id_clean")]
tax_mat[] <- lapply(tax_mat, as.numeric)

prev <- colMeans(tax_mat > 0, na.rm = TRUE)
avg <- colMeans(tax_mat, na.rm = TRUE)

tax_stat <- data.frame(taxon = colnames(tax_mat),prevalence = prev,
                       mean_abund = avg,level = lengths(strsplit(colnames(tax_mat), "[|]")))

tax_stat <- tax_stat %>% filter(prevalence >= prev_cut, mean_abund > 0)

if (species_only) {
  tax_stat <- tax_stat %>% filter(level == 7)
}

tax_stat <- tax_stat %>%
  arrange(desc(prevalence), desc(mean_abund))

top_taxa <- head(tax_stat$taxon, K)

Y <- tax_mat[, top_taxa, drop = FALSE]

# reorder by subject and time
ord <- order(meta_T$`Participant ID`, meta_T$week_num)

meta_T <- meta_T[ord, ]
Y <- Y[ord, , drop = FALSE]

meta_T <- meta_T %>%
  group_by(`Participant ID`) %>%
  mutate(t_index = row_number()) %>%
  ungroup()

# make N x T x K array
ids <- unique(meta_T$`Participant ID`)
N <- length(ids)

Y_array <- array(NA_real_,dim = c(N, T_len, K),dimnames = list(ids, paste0("t", 1:T_len), top_taxa))

for (i in seq_along(ids)) {
  ii <- which(meta_T$`Participant ID` == ids[i])
  Y_array[i, , ] <- as.matrix(Y[ii, ])
}

c(N = N, T_len = T_len, K = K)
top_taxa