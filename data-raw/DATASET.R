## Preparation of `ecoli_bw25113` dataset and tn5 transposon insertion sites
# uses dplyr and readr

# data from https://biocyc.org/GCF_000750555/organism-summary
ecoli_bw25113 <- readr::read_tsv("./data-raw/All-genes-of-E.-coli-K-12-substr.-BW25113.txt")

# rename, sort by start and drop Product and Accession-1
ecoli_bw25113 <- ecoli_bw25113[, c(1,3,4)]
names(ecoli_bw25113) <- c("gene", "start", "end")
ecoli_bw25113 <- ecoli_bw25113 %>% arrange(start)

# Load final IS from Goodall et al. (2018),
is_count_neg <-
  read.table("./data-raw/Ecoli_BW25113_chlor-tn_position_depth_count_neg.txt")
is_count_pos <-
  read.table("./data-raw/Ecoli_BW25113_chlor-tn_position_depth_count_pos.txt")

# Combine + & - (in total 884,675 IS)
is_pos <- which(is_count_neg+is_count_pos > 0)

# make a less dense library
set.seed(1)
is_pos <- sample(is_pos, 100000, replace = F)

# save data
usethis::use_data(ecoli_bw25113, overwrite = TRUE)
usethis::use_data(is_pos, overwrite = TRUE)


