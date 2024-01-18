
# VCFGENIE ANALYSIS SCRIPT ------------------------------------------------

# Author: Chase W. Nelson
# Email: chase.nelson@nih.gov
# Date created: 2023-07-24
# Last updated: 2024-01-18
# Description: run manually line-by-line with input data available from github
# URL1: https://github.com/NCI-CGR/HPV_low_VAF_SNV_prediction
# URL2: https://github.com/chasewnelson/VCFgenie
# Citation: Mishra SK, Nelson CW, ... Yeager M, Mirabello L (2024) "Improved detection of low-frequency within-host variants from deep sequencing: A case study with human papillomavirus"

library(tidyverse)
library(scales)


# IMPORT Pre-VCFgenie Variants --------------------------------------------

# BEFORE VCFgenie
preVCFgenie <- read_tsv("initial_iSNVs_preVCFgenie_maskID.txt")
nrow(preVCFgenie)  # 9225

# create unique ID
preVCFgenie$ID_uniq <- paste0(preVCFgenie$sample, '_', preVCFgenie$POS, '_', preVCFgenie$REF, '_', preVCFgenie$ALT)
length(unique(preVCFgenie$ID_uniq))  # 9225 sample/POS/REF/ALT combos

# AFTER VCFenie, AC=0 for failing variants
postVCFgenie <- read_tsv("initial_iSNVs_postVCFgenie_maskID.txt")
nrow(postVCFgenie)  # 9225

# create unique ID
postVCFgenie$ID_uniq <- paste0(postVCFgenie$sample, '_', postVCFgenie$POS, '_', postVCFgenie$REF, '_', postVCFgenie$ALT)
length(unique(postVCFgenie$ID_uniq))  # 9225 sample/POS/REF/ALT combos


# PRE-VCFGENIE: Bookkeeping -----------------------------------------------

# Samples
length(unique(preVCFgenie$sample))  # 93, but including replicates
93 / 3  # 31 distinct samples

# BEFORE: AO, FAO; trim the replicate suffix e.g. *_C1
length(unique(str_replace(preVCFgenie$sample, pattern = "(sample\\d+)_.+", replacement = "\\1")))  # 31, QED

nrow(filter(preVCFgenie, AO == 0))  # 1
nrow(filter(preVCFgenie, FAO == 0))  # 70
nrow(filter(preVCFgenie, FDP < 100))  # 300
nrow(filter(preVCFgenie, FAO == 0, FDP < 100))  # 0 overlap

# Filter OUT 0 FAO (failing flow-correction)
preVCFgenie_FAOgt0 <- filter(preVCFgenie, FAO > 0)
nrow(preVCFgenie_FAOgt0)  # 9155
(nrow(preVCFgenie) - nrow(preVCFgenie_FAOgt0))  # 70 removed, QED


# Visualize Depth (FDP) ---------------------------------------------------

# Coverage comparison
(preVCFgenie_FDPbyRep <- preVCFgenie |>
    group_by(sample) |>
    summarise(
       median_FDP = median(FDP)
    ))

# rename
preVCFgenie_FDPbyRep <- preVCFgenie_FDPbyRep |>
   rename(sample_rep = sample)

# add sample
(preVCFgenie_FDPbyRep$sample <- str_replace(preVCFgenie_FDPbyRep$sample_rep, "_\\w+$", ""))

# add replicate of sample
(preVCFgenie_FDPbyRep$rep <- str_replace(preVCFgenie_FDPbyRep$sample_rep, "^\\w+_", ""))

# relocate columns
(preVCFgenie_FDPbyRep <- preVCFgenie_FDPbyRep |>
      relocate(sample_rep, sample, rep, median_FDP))

# arrange by FDP
(preVCFgenie_FDPbyRep <- arrange(preVCFgenie_FDPbyRep, sample, median_FDP))

# every sample have 3 replicates?
preVCFgenie_FDPbyRep |> 
   count(sample) |>
   filter(n != 3)
# none QED

# number replicates by FDP
preVCFgenie_FDPbyRep$rep_num <- rep(1:3, nrow(preVCFgenie_FDPbyRep) / 3)

# arrange by num_variants
(preVCFgenie_FDPbyRep <- arrange(preVCFgenie_FDPbyRep))

# factor by this order
preVCFgenie_FDPbyRep$sample <- factor(preVCFgenie_FDPbyRep$sample, levels = unique(preVCFgenie_FDPbyRep$sample))

# max coverage to scale second y axis
(max_median_FDP <- max(preVCFgenie_FDPbyRep$median_FDP))  # 1997


# POST-VCFGENIE: Bookkeeping ----------------------------------------------

# AFTER: AO, FAO
nrow(filter(postVCFgenie, AO == 0))  # 1 QED
nrow(filter(postVCFgenie, FAO == 0))  # 86, this is higher than 70 because VCFgenie reset some (86-70 == 16) to 0
nrow(filter(postVCFgenie, FDP < 100))  # 300 QED
nrow(filter(postVCFgenie, FDP < 100, FAO == 0))  # 3, VCFgenie set 3 with FDP < 100 to FAO = 0

# FILTER OUT 0 FAO
postVCFgenie_FAOgt0 <- filter(postVCFgenie, FAO > 0)
nrow(postVCFgenie_FAOgt0)  # 9139
9139 / 9155  # 99.8% of the full Initial Set PASS VCFgenie
(nrow(preVCFgenie) - nrow(postVCFgenie_FAOgt0))  # 86 removed, QED

# VCFGENIE ELIMINATED
(nrow(preVCFgenie_FAOgt0) - nrow(postVCFgenie_FAOgt0))  # 16
# ^after flow corrected quality values, VCFgenie eliminated 16


# Characterize Eliminated Variants ----------------------------------------

(VCFgenie_elim_IDs <- setdiff(preVCFgenie_FAOgt0$ID_uniq, postVCFgenie_FAOgt0$ID_uniq))

preVCFgenie_FAOgt0$pass_VCFgenie <- ! preVCFgenie_FAOgt0$ID_uniq %in% VCFgenie_elim_IDs
sum(! preVCFgenie_FAOgt0$pass_VCFgenie)  # 16, QED
(preVCFgenie_FAOgt0$pass_VCFgenie <- factor(preVCFgenie_FAOgt0$pass_VCFgenie, levels = c(FALSE, TRUE), labels = c('FAIL', 'PASS')))


# Confusion Matrices ------------------------------------------------------
# Supplementary Fig. S8

# PROPORTIONS of VCFgenie pass/fail that are true/false, i.e., CONFUSION MATRIX
preVCFgenie_FAOgt0$variant_status <- as.character(NA)
preVCFgenie_FAOgt0[near(preVCFgenie_FAOgt0$PERCENT_OVERLAP, 1/3, tol = .1), ]$variant_status <- "FALSE"
preVCFgenie_FAOgt0[near(preVCFgenie_FAOgt0$PERCENT_OVERLAP, 2/3, tol = .1), ]$variant_status <- "AMBIGUOUS"
preVCFgenie_FAOgt0[near(preVCFgenie_FAOgt0$PERCENT_OVERLAP, 3/3, tol = .1), ]$variant_status <- "TRUE"

# factor order
preVCFgenie_FAOgt0$variant_status <- factor(preVCFgenie_FAOgt0$variant_status, levels = c('TRUE', 'AMBIGUOUS', 'FALSE'))
preVCFgenie_FAOgt0$pass_VCFgenie <- factor(preVCFgenie_FAOgt0$pass_VCFgenie, levels = c('PASS', 'FAIL'))

preVCFgenie_FAOgt0 %>% group_by(variant_status, pass_VCFgenie) %>%
 summarise(
  count = n()
 )

# Summarize by replicate frequency
preVCFgenie_FAOgt0 |> 
   group_by(variant_status) |> 
   summarize(
      count = n(),
      .groups = "drop"
   ) |> 
   mutate(
      prop = count / sum(count)
   )

# Summarize by VCFgenie result
preVCFgenie_FAOgt0 |> 
   group_by(pass_VCFgenie) |> 
   summarize(
      count = n(),
      .groups = "drop"
   ) |> 
   mutate(
      prop = count / sum(count)
   )

# SUM OF ALL
nrow(preVCFgenie_FAOgt0)
# 9155

# SUMMARIZE
(preVCFgenie_FAOgt0_CONFUSION <- preVCFgenie_FAOgt0 %>% group_by(variant_status, pass_VCFgenie) %>%
  summarise(
   n = n()
  ))

sum(preVCFgenie_FAOgt0_CONFUSION$n)  # 9155
sum(filter(preVCFgenie_FAOgt0_CONFUSION, pass_VCFgenie == 'PASS')$n)  # 9139
9139/9155  # 0.9982523
sum(filter(preVCFgenie_FAOgt0_CONFUSION, pass_VCFgenie == 'FAIL')$n)  # 16
16/9155  # 0.001747679
sum(filter(preVCFgenie_FAOgt0_CONFUSION, variant_status == 'TRUE')$n)  # 4810
4810/9155  # 0.525396
sum(filter(preVCFgenie_FAOgt0_CONFUSION, variant_status == 'AMBIGUOUS')$n)  # 673
673/9155  # 0.07351174
sum(filter(preVCFgenie_FAOgt0_CONFUSION, variant_status == 'FALSE')$n)  # 3672
3672/9155  # 0.4010923


# Minor iSNVs -------------------------------------------------------------

# SUMMARIZE
(preVCFgenie_FAOgt0_CONFUSION_MINOR <- filter(preVCFgenie_FAOgt0, AF < 0.5) |> 
    group_by(variant_status, pass_VCFgenie) %>%
    summarise(
       n = n()
       )
 )

# Summarize by replicate frequency
preVCFgenie_FAOgt0 |> 
   filter(AF < 0.5) |> # <== MINOR
   group_by(variant_status) |> 
   summarize(
      count = n(),
      .groups = "drop"
   ) |> 
   mutate(
      prop = count / sum(count)
   )

# Summarize by VCFgenie result
preVCFgenie_FAOgt0 |> 
   filter(AF < 0.5) |> # <== MINOR
   group_by(pass_VCFgenie) |> 
   summarize(
      count = n(),
      .groups = "drop"
   ) |> 
   mutate(
      prop = count / sum(count)
   )

sum(preVCFgenie_FAOgt0_CONFUSION_MINOR$n)  # 4614
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MINOR, pass_VCFgenie == 'PASS')$n)  # 4598
4598 / 4614  # 0.9965323
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MINOR, pass_VCFgenie == 'FAIL')$n)  # 16
16 / 4614  # 0.003467707
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MINOR, variant_status == 'TRUE')$n)  # 469
469 / 4614  # 0.1016472
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MINOR, variant_status == 'AMBIGUOUS')$n)  # 550
550 / 4614  # 0.1192024
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MINOR, variant_status == 'FALSE')$n)  # 3595
3595 / 4614  # 0.7791504


# Major iSNVs -------------------------------------------------------------

# SUMMARIZE
(preVCFgenie_FAOgt0_CONFUSION_MAJOR <- filter(preVCFgenie_FAOgt0, AF >= 0.5) |> 
    group_by(variant_status, pass_VCFgenie) %>%
    summarise(
       n = n()
       )
 )

# Summarize by replicate frequency
preVCFgenie_FAOgt0 |> 
   filter(AF >= 0.5) |> # <== MAJOR
   group_by(variant_status) |> 
   summarize(
      count = n(),
      .groups = "drop"
   ) |> 
   mutate(
      prop = count / sum(count)
   )

# Summarize by VCFgenie result
preVCFgenie_FAOgt0 |> 
   filter(AF >= 0.5) |> # <== MAJOR
   group_by(pass_VCFgenie) |> 
   summarize(
      count = n(),
      .groups = "drop"
   ) |> 
   mutate(
      prop = count / sum(count)
   )

sum(preVCFgenie_FAOgt0_CONFUSION_MAJOR$n)  # 4541
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MAJOR, pass_VCFgenie == 'PASS')$n)  # 4541
4541 / 4541  # 1
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MAJOR, pass_VCFgenie == 'FAIL')$n)  # 0
0 / 4541  # 0
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MAJOR, variant_status == 'TRUE')$n)  # 4341
4341 / 4541  # 0.9559568
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MAJOR, variant_status == 'AMBIGUOUS')$n)  # 123
123 / 4541  # 0.02708654
sum(filter(preVCFgenie_FAOgt0_CONFUSION_MAJOR, variant_status == 'FALSE')$n)  # 77
77 / 4541  # 0.01695662


# Lowest-VAF Passing Variant ----------------------------------------------

# arrange and rank by increasing FDP
preVCFgenie_FAOgt0 <- arrange(preVCFgenie_FAOgt0, FDP)
preVCFgenie_FAOgt0$FDP_rank <- 1:nrow(preVCFgenie_FAOgt0)

# arrange and rank by increasing AF
preVCFgenie_FAOgt0$AF_check <- preVCFgenie_FAOgt0$FAO / preVCFgenie_FAOgt0$FDP
preVCFgenie_FAOgt0 <- arrange(preVCFgenie_FAOgt0, AF)
preVCFgenie_FAOgt0$AF_rank <- 1:nrow(preVCFgenie_FAOgt0)

# AF summary
summary(filter(preVCFgenie_FAOgt0, pass_VCFgenie == 'PASS')$AF)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003024 0.039013 0.400700 0.512718 1.000000 1.000000 

summary(filter(preVCFgenie_FAOgt0, pass_VCFgenie == 'FAIL')$AF)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005086 0.0081513 0.0132517 0.0119098 0.0156930 0.0263158

# FDP summary
summary(filter(preVCFgenie_FAOgt0, pass_VCFgenie == 'FAIL')$FDP)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 38.0   111.8   126.0   348.3   167.8  1966.0


# Pass vs. Fail Metrics ---------------------------------------------------
# Fig. 6

(VCFgenie_pass_summary <- preVCFgenie_FAOgt0 %>%
 group_by(pass_VCFgenie) %>%
 summarise(
  QUAL = mean(QUAL),
  AF = mean(AF),
  AO = mean(AO),
  DP = mean(DP),
  FAO = mean(FAO),
  FDP = mean(FDP),
  FRO = mean(FRO),
  FSAF = mean(FSAF),
  FSAR = mean(FSAR),
  FSRF = mean(FSRF),
  FSRR = mean(FSRR),
  FWDB = mean(FWDB),
  FXX = mean(FXX),
  HRUN = mean(HRUN),
  LEN = mean(LEN),
  MLLD = mean(MLLD),
  QD = mean(QD),
  RBI = mean(RBI),
  REFB = mean(REFB),
  REVB = mean(REVB),
  RO = mean(RO),
  SAF = mean(SAF),
  SAR = mean(SAR),
  SRF = mean(SRF),
  SRR = mean(SRR),
  SSEN = mean(SSEN),
  SSEP = mean(SSEP),
  SSSB = mean(SSSB),
  STB = mean(STB),
  STBP = mean(STBP),
  VARB = mean(VARB),
  GQ = mean(GQ),
  PERCENT_OVERLAP = mean(PERCENT_OVERLAP),
  COVERAGE = mean(COVERAGE)
 ))

# LONG
(preVCFgenie_FAOgt0_LONG <- preVCFgenie_FAOgt0 %>% 
  pivot_longer(cols = c('AF',
                        'AO',
                        'DP',
                        'FAO',
                        'FDP',
                        'FRO',
                        'FSAF',
                        'FSAR',
                        'FSRF',
                        'FSRR',
                        'FWDB',
                        'FXX',
                        'GQ',
                        'HRUN',
                        'MLLD',
                        'QD',
                        'QUAL',
                        'RBI',
                        'REFB',
                        'REVB',
                        'RO',
                        'SAF',
                        'SAR',
                        'SRF',
                        'SRR',
                        'SSSB',
                        'STB',
                        'STBP',
                        'VARB'), 
               names_to = 'feature',
               values_to = 'value') %>%
  select(sample, CHROM, POS, ID, REF, ALT, FILTER, pass_VCFgenie, feature, value)
)

# code good=higher values better, bad=higher values worse
good_features <- c('AF', 'AO', 'DP', 'FAO', 'FDP', 'FRO', 'FSAF', 'FSAR', 'FSRF', 'FSRR', 'GQ', 'MLLD', 'QD', 'QUAL', 'RO', 'SAF', 'SAR', 'SRF', 'SRR')
bad_features <- c('FWDB', 'FXX', 'HRUN', 'RBI', 'REFB', 'REVB', 'SSSB', 'STB', 'STBP', 'VARB')

preVCFgenie_FAOgt0_LONG$feature_type <- 'good'
preVCFgenie_FAOgt0_LONG[preVCFgenie_FAOgt0_LONG$feature %in% bad_features, ]$feature_type <- 'bad'


# Compare Feature Means & Medians -----------------------------------------

# determine which means / medians higher, etc.
(preVCFgenie_FAOgt0_LONG_SUMMARY <- preVCFgenie_FAOgt0_LONG %>% group_by(feature_type, feature, pass_VCFgenie) %>%
 summarise(count = n(),
           mean = mean(value),
           median = median(value)))


# MOLECULAR SPECTRUM of SEQUENCING ERRORS ---------------------------------

# form trinuc
(preVCFgenie_FAOgt0$nuc1 <- str_replace(string = preVCFgenie_FAOgt0$TRI_NUCLEOTIDE_CONTEXT, pattern = "(\\w).\\w", replacement = "\\1"))
(preVCFgenie_FAOgt0$nuc2 <- str_replace(string = preVCFgenie_FAOgt0$Mutation, pattern = "(\\w)>\\w", replacement = "\\1"))
(preVCFgenie_FAOgt0$nuc3 <- str_replace(string = preVCFgenie_FAOgt0$TRI_NUCLEOTIDE_CONTEXT, pattern = "\\w.(\\w)", replacement = "\\1"))

# from / to
(preVCFgenie_FAOgt0$from_trinuc <- paste0(preVCFgenie_FAOgt0$nuc1, preVCFgenie_FAOgt0$nuc2, preVCFgenie_FAOgt0$nuc3))
(preVCFgenie_FAOgt0$from <- str_replace(string = preVCFgenie_FAOgt0$Mutation, pattern = "(\\w)>\\w", replacement = "\\1"))
(preVCFgenie_FAOgt0$to <- str_replace(string = preVCFgenie_FAOgt0$Mutation, pattern = "\\w>(\\w)", replacement = "\\1"))

# Initialize nuc / trinuc variables

# Sequence content
nucs <- c('A', 'C', 'G', 'T')
trinucs <- c('AAA', 'AAC', 'AAG', 'AAT', 
             'ACA', 'ACC', 'ACG', 'ACT', 
             'AGA', 'AGC', 'AGG', 'AGT', 
             'ATA', 'ATC', 'ATG', 'ATT',
             'CAA', 'CAC', 'CAG', 'CAT', 
             'CCA', 'CCC', 'CCG', 'CCT', 
             'CGA', 'CGC', 'CGG', 'CGT', 
             'CTA', 'CTC', 'CTG', 'CTT', 
             'GAA', 'GAC', 'GAG', 'GAT', 
             'GCA', 'GCC', 'GCG', 'GCT', 
             'GGA', 'GGC', 'GGG', 'GGT', 
             'GTA', 'GTC', 'GTG', 'GTT', 
             'TAA', 'TAC', 'TAG', 'TAT', 
             'TCA', 'TCC', 'TCG', 'TCT', 
             'TGA', 'TGC', 'TGG', 'TGT', 
             'TTA', 'TTC', 'TTG', 'TTT')

# Mutations
mutation_to_mutationType <- c("T>G", "T>C", "T>A", "C>A", "C>G", "C>T", "C>T", "C>G", "C>A", "T>A", "T>C", "T>G")
names(mutation_to_mutationType) <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
mutations_used <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mutations_used_labels <- c("C%->%A", "C%->%G", "C%->%T", "T%->%A", "T%->%C", "T%->%G")

# WRANGLE
(preVCFgenie_FAOgt0_SUMMARY <- preVCFgenie_FAOgt0 %>% group_by(variant_status, from_trinuc, from, to) %>%
 summarise(
  count = n()
 ))

# ADD 0s for those that don't exist
nrow(preVCFgenie_FAOgt0_SUMMARY)  # 346, but want 3*(4*4*4)*3 = 576

for(this_variant_status in c('TRUE', 'AMBIGUOUS', 'FALSE')) {
 for (this_from_trinuc in trinucs) {
  
  # this_from is the middle of the trinuc
  this_from <- str_replace(string = this_from_trinuc, pattern = "\\w(\\w)\\w", replacement = "\\1")
  
  # this_to can only be those nucs that are NOT the from_to
  for (this_to in setdiff(nucs, this_from)) {
   # if there isn't already an entry, add the row with count of 0
   if (nrow(filter(preVCFgenie_FAOgt0_SUMMARY, variant_status == this_variant_status, from_trinuc == this_from_trinuc, from == this_from, to == this_to)) == 0) {
    preVCFgenie_FAOgt0_SUMMARY <- rbind(preVCFgenie_FAOgt0_SUMMARY,
                                        tibble(variant_status = this_variant_status,
                                               from_trinuc = this_from_trinuc,
                                               from = this_from,
                                               to = this_to,
                                               count = 0))
   }
  }
 }
}

preVCFgenie_FAOgt0_SUMMARY
nrow(preVCFgenie_FAOgt0_SUMMARY)  # 576, QED!

# Add mutation
preVCFgenie_FAOgt0_SUMMARY$mutation <- paste0(preVCFgenie_FAOgt0_SUMMARY$from, '>', preVCFgenie_FAOgt0_SUMMARY$to)

# Add mutation type
preVCFgenie_FAOgt0_SUMMARY$mutation_type <- preVCFgenie_FAOgt0_SUMMARY$mutation
preVCFgenie_FAOgt0_SUMMARY$mutation_type <- mutation_to_mutationType[preVCFgenie_FAOgt0_SUMMARY$mutation]
preVCFgenie_FAOgt0_SUMMARY[is.na(preVCFgenie_FAOgt0_SUMMARY$mutation_type), ]$mutation_type <- 'identity'

# Arrange columns
preVCFgenie_FAOgt0_SUMMARY <- dplyr::select(preVCFgenie_FAOgt0_SUMMARY, mutation_type, mutation, from_trinuc, from, to, everything())
preVCFgenie_FAOgt0_SUMMARY

# FACTOR: add mutation type labels
preVCFgenie_FAOgt0_SUMMARY$mutation_type <- factor(x = preVCFgenie_FAOgt0_SUMMARY$mutation_type,
                                        levels = mutations_used,
                                        labels = mutations_used_labels)

# FACTOR: variant status
preVCFgenie_FAOgt0_SUMMARY$variant_status <- factor(preVCFgenie_FAOgt0_SUMMARY$variant_status,
                                                    levels = c('TRUE', 'AMBIGUOUS', 'FALSE'),
                                                    labels = c('True', 'Ambiguous', 'False'))


# Count HPV16 Trinucs -----------------------------------------------------

# COUNT TRINUCLEOTIDES in the genome of HPV18
(HPV18REF_trinuc_counts <- read_tsv("HPV18REF_3mer_counts.tsv"))
names(HPV18REF_trinuc_counts) <- c('from_trinuc', 'from_trinuc_count')

# JOIN
preVCFgenie_FAOgt0_SUMMARY <- left_join(x = preVCFgenie_FAOgt0_SUMMARY, y = HPV18REF_trinuc_counts, by = 'from_trinuc')

# RATE
preVCFgenie_FAOgt0_SUMMARY$count_per_trinuc <- preVCFgenie_FAOgt0_SUMMARY$count / preVCFgenie_FAOgt0_SUMMARY$from_trinuc_count

# RATE SUMS PER STATUS
(preVCFgenie_FAOgt0_SUMMARY_StatusSums <- preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status) %>%
 summarise(
  status_rate_sum = sum(count_per_trinuc)
 ))

# JOIN SUMS
(preVCFgenie_FAOgt0_SUMMARY <- left_join(x = preVCFgenie_FAOgt0_SUMMARY, y = preVCFgenie_FAOgt0_SUMMARY_StatusSums, by = 'variant_status'))

# NORMALIZED RATE
preVCFgenie_FAOgt0_SUMMARY$norm_count_per_trinuc <- preVCFgenie_FAOgt0_SUMMARY$count_per_trinuc / preVCFgenie_FAOgt0_SUMMARY$status_rate_sum

# VERIFY that each STATUS' rate sum is equal to 1.0
preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status) %>%
 summarise(
  norm_count_per_trinuc_sum = sum(norm_count_per_trinuc)
 )
# QED


# Trinuc - Pool Strands ---------------------------------------------------

trinucs_to_woStrand <- c('TTT', 'GTT', 'CTT', 'ATT', 
                      'ACA', 'ACC', 'ACG', 'ACT', 
                      'TCT', 'GCT', 'CCT', 'ACT', 
                      'ATA', 'ATC', 'ATG', 'ATT',
                      'TTG', 'GTG', 'CTG', 'ATG', 
                      'CCA', 'CCC', 'CCG', 'CCT', 
                      'TCG', 'GCG', 'CCG', 'ACG', 
                      'CTA', 'CTC', 'CTG', 'CTT', 
                      'TTC', 'GTC', 'CTC', 'ATC', 
                      'GCA', 'GCC', 'GCG', 'GCT', 
                      'TCC', 'GCC', 'CCC', 'ACC', 
                      'GTA', 'GTC', 'GTG', 'GTT', 
                      'TTA', 'GTA', 'CTA', 'ATA', 
                      'TCA', 'TCC', 'TCG', 'TCT', 
                      'TCA', 'GCA', 'CCA', 'ACA', 
                      'TTA', 'TTC', 'TTG', 'TTT')
names(trinucs_to_woStrand) <- trinucs

preVCFgenie_FAOgt0_SUMMARY$from_trinuc_woStrand <- trinucs_to_woStrand[preVCFgenie_FAOgt0_SUMMARY$from_trinuc]

# re-summarize
(preVCFgenie_FAOgt0_SUMMARY_woStrand <- preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status, from_trinuc_woStrand, mutation_type) %>%
 summarise(
  norm_count_per_trinuc_sum = sum(norm_count_per_trinuc)
 ))

# Sample sizes
preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status) %>%
 summarize(
  sum = sum(count)
 )

# CONCORDANCE
(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE <- preVCFgenie_FAOgt0_SUMMARY_woStrand %>% 
 pivot_wider(names_from = variant_status,
             values_from = norm_count_per_trinuc_sum))

# direction of change
(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous <- preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous - preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True)
(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False <- preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$False - preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous)
(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_False <- preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$False - preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True)

# factor
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous_TYPE <- as.character(NA)
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE[preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous > 0, ]$True_to_Ambiguous_TYPE <- 'Increase'
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE[preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous < 0, ]$True_to_Ambiguous_TYPE <- 'Decrease'

preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False_TYPE <- as.character(NA)
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE[preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False > 0, ]$Ambiguous_to_False_TYPE <- 'Increase'
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE[preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False < 0, ]$Ambiguous_to_False_TYPE <- 'Decrease'

# Concordance of change
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$concordance <- as.character(NA)
preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE[! is.na(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous_TYPE) &
                                          ! is.na(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False_TYPE) &
                                          preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous_TYPE == preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False_TYPE, 
 ]$concordance <- 'Concordant'

preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE[! is.na(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous_TYPE) &
                                          ! is.na(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False_TYPE) &
                                          preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_Ambiguous_TYPE != preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$Ambiguous_to_False_TYPE, 
]$concordance <- 'Disconcordant'

preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE
# view(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE)

# summarise concordance
(concordance_SUMMARY <- preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE %>% group_by(concordance) %>%
 summarise(
  count = n()
 ))

(concordance_trinuc_SUMMARY <- preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE %>% group_by(mutation_type, concordance) %>%
  summarise(
   count = n()
  ))

# Largest increase or decrease
summary(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE$True_to_False)
arrange(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE, True_to_False)
arrange(preVCFgenie_FAOgt0_SUMMARY_woStrand_WIDE, -True_to_False)

# Re-summarize strand agnostic RAW COUNTS
(preVCFgenie_FAOgt0_SUMMARY_woStrand_COUNTS <- preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status, from_trinuc_woStrand, mutation_type) %>%
  summarise(
   count_sum = sum(count)
  ))


# TALLY - Nuc iSNV Type Counts --------------------------------------------

# re-summarize
(preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS <- preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status, mutation, mutation_type) %>%
  summarise(
   count_sum = sum(count)
  ))

# double check totals
(preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_SUMS <- preVCFgenie_FAOgt0_SUMMARY %>% group_by(variant_status) %>%
 summarise(
  n = sum(count)
 ))

# JOIN
(preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN <- left_join(x = preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS,
                                                         y = preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_SUMS,
                                                         by = "variant_status"))

# add prop
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$prop <- preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$count_sum / preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$n

# add cols
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$CI_lo <- as.numeric(NA)
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$CI_hi <- as.numeric(NA)
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$X_squared <- as.numeric(NA)
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$p_value <- as.numeric(NA)
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$prop_confirm <- as.numeric(NA)

# proportions and CIs
for(this_variant_status in unique(preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$variant_status)) {
 # this_variant_status <- 'True'
 this_data <- filter(preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN, variant_status == this_variant_status)
 this_n <- preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_SUMS[preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_SUMS$variant_status == this_variant_status, ]$n
 
 for(this_mutation in unique(this_data$mutation)) {
  # this_mutation <- 'A>C'
  this_x <- this_data[this_data$mutation == this_mutation, ]$count_sum
  this_prop_test <- prop.test(x = this_x, n = this_n, alternative = 'two.sided', conf.level = 0.95)
  
  # fill data into row
  preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN[preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$variant_status == this_variant_status &
                                              preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$mutation == this_mutation, ]$CI_lo <- this_prop_test$conf.int[1]
  preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN[preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$variant_status == this_variant_status &
                                              preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$mutation == this_mutation, ]$CI_hi <- this_prop_test$conf.int[2]
  preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN[preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$variant_status == this_variant_status &
                                              preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$mutation == this_mutation, ]$X_squared <- unname(this_prop_test$statistic)
  preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN[preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$variant_status == this_variant_status &
                                              preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$mutation == this_mutation, ]$p_value <- this_prop_test$p.value
  preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN[preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$variant_status == this_variant_status &
                                              preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$mutation == this_mutation, ]$prop_confirm <- unname(this_prop_test$estimate)
 }
}

# check
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN
sum(! near(preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$prop, preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$prop_confirm))  # 0 QED
sum(! (preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$CI_lo < preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$prop) &
     (preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$CI_hi > preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN$prop))  # 0 QED
# GREAT!

# summarize
preVCFgenie_FAOgt0_SUMMARY_nuc_COUNTS_JOIN %>% group_by(variant_status, mutation_type) %>%
 summarise(
  n = sum(count_sum)
 )


# Seq Content -------------------------------------------------------------
# SEQ CONTENT & add FEATURE NAME $ add FEATURE PROP

# nuc counts
(HPV18REF_nuc_counts <- read_tsv("HPV18REF_1mer_counts.tsv"))
HPV18REF_nuc_counts$feature <- 'nuc'
HPV18REF_nuc_counts$prop <- HPV18REF_nuc_counts$count / sum(HPV18REF_nuc_counts$count)
HPV18REF_nuc_counts$pct <- percent(HPV18REF_nuc_counts$prop, accuracy = .1)
HPV18REF_nuc_counts
nrow(HPV18REF_nuc_counts)  # 4

# CpG status
HPV18REF_nuc_counts$motif_status <- 'Other'  # 'CG-only'
HPV18REF_nuc_counts[HPV18REF_nuc_counts$kmer %in% c('A', 'T'), ]$motif_status <- 'Other'  # 'AT-only'
HPV18REF_nuc_counts

# dinuc counts
(HPV18REF_dinuc_counts <- read_tsv("HPV18REF_2mer_counts.tsv"))
HPV18REF_dinuc_counts$feature <- 'dinuc'
HPV18REF_dinuc_counts$prop <- HPV18REF_dinuc_counts$count / sum(HPV18REF_dinuc_counts$count)
HPV18REF_dinuc_counts$pct <- percent(HPV18REF_dinuc_counts$prop, accuracy = .1)
HPV18REF_dinuc_counts
nrow(HPV18REF_dinuc_counts)  # 16

# dinucs
CpG_dinucs <- c('CG')
TpC_dinucs <- c('TC', 'GA')

HPV18REF_dinuc_counts$motif_status <- "Other"
HPV18REF_dinuc_counts[HPV18REF_dinuc_counts$kmer %in% CpG_dinucs, ]$motif_status <- 'CpG'
HPV18REF_dinuc_counts[HPV18REF_dinuc_counts$kmer %in% TpC_dinucs, ]$motif_status <- 'TpC'
HPV18REF_dinuc_counts

# trinuc counts
(HPV18REF_trinuc_counts <- read_tsv("HPV18REF_3mer_counts.tsv"))
HPV18REF_trinuc_counts$feature <- 'trinuc'
HPV18REF_trinuc_counts$prop <- HPV18REF_trinuc_counts$count / sum(HPV18REF_trinuc_counts$count)
HPV18REF_trinuc_counts$pct <- percent(HPV18REF_trinuc_counts$prop, accuracy = .1)
HPV18REF_trinuc_counts
nrow(HPV18REF_trinuc_counts)  # 64

CpG_trinucs <- c('CGA', 'CGC', 'CGG', 'CGT', 'ACG', 'CCG', 'GCG', 'TCG')  # revcom CpG, so same on both strands
TpC_trinucs <- c('TCA', 'TCC', 'TCG', 'TCT', 'ATC', 'CTC', 'GTC', 'TTC',  # revcom CpA, next line
                 'TGA', 'GGA', 'CGA', 'AGA', 'GAT', 'GAG', 'GAC', 'GAA')

CpG_and_TpC_trinucs <- intersect(CpG_trinucs, TpC_trinucs)  # [1] "CGA" "TCG"

HPV18REF_trinuc_counts$motif_status <- "Other"
HPV18REF_trinuc_counts[HPV18REF_trinuc_counts$kmer %in% CpG_trinucs, ]$motif_status <- 'CpG'
HPV18REF_trinuc_counts[HPV18REF_trinuc_counts$kmer %in% TpC_trinucs, ]$motif_status <- 'TpC'
HPV18REF_trinuc_counts[HPV18REF_trinuc_counts$kmer %in% CpG_and_TpC_trinucs, ]$motif_status <- 'CpG and TpC'
HPV18REF_trinuc_counts


# Combine Nuc, Dinuc, Trinuc ----------------------------------------------

# COMBINE ALL
HPV18REF_kmer_counts <- rbind(HPV18REF_nuc_counts, HPV18REF_dinuc_counts, HPV18REF_trinuc_counts)
nrow(HPV18REF_kmer_counts)  # 84, QED

# FACTOR CpG status
HPV18REF_kmer_counts$motif_status <- factor(HPV18REF_kmer_counts$motif_status,
                                            levels = c('CpG', 'CpG and TpC', 'TpC', 'Other'))

# FACTOR feature
HPV18REF_kmer_counts$feature <- factor(HPV18REF_kmer_counts$feature,
                                       levels = c('nuc', 'dinuc', 'trinuc'))

# calc content
sum(filter(HPV18REF_kmer_counts, kmer %in% c('C', 'G'))$prop)  # G:C content = 0.4043528
filter(HPV18REF_kmer_counts, kmer == 'CG')$prop  # CpG = 0.02189409
filter(HPV18REF_kmer_counts, kmer == 'TC')$prop  # TpC = 0.03080448


# SFS ---------------------------------------------------------------------

preVCFgenie_FAOgt0

# Bin by frequency in 5% steps
preVCFgenie_FAOgt0$AF_binned <- cut(x = preVCFgenie_FAOgt0$AF, breaks = seq(0, 1, by = .05))
str(preVCFgenie_FAOgt0$AF_binned)
preVCFgenie_FAOgt0$AF_binned <- factor(preVCFgenie_FAOgt0$AF_binned,
                             levels = c('(0,0.05]', '(0.05,0.1]', '(0.1,0.15]', '(0.15,0.2]', '(0.2,0.25]', '(0.25,0.3]', '(0.3,0.35]', '(0.35,0.4]', '(0.4,0.45]', '(0.45,0.5]', '(0.5,0.55]', '(0.55,0.6]', '(0.6,0.65]', '(0.65,0.7]', '(0.7,0.75]', '(0.75,0.8]', '(0.8,0.85]', '(0.85,0.9]', '(0.9,0.95]', '(0.95,1]'),
                             labels = c('0-5', '5-10', '10-15', '15-20', '20-25', '25-30', '30-35', '35-40', '40-45', '45-50', '50-55', '55-60', '60-65', '65-70', '70-75', '75-80', '80-85', '85-90', '90-95', '95-100'))

# 10% steps
preVCFgenie_FAOgt0$AF_binned10 <- cut(x = preVCFgenie_FAOgt0$AF, breaks = seq(0, 1, by = .1))
str(preVCFgenie_FAOgt0$AF_binned10)
preVCFgenie_FAOgt0$AF_binned10 <- factor(preVCFgenie_FAOgt0$AF_binned10,
                               levels = c('(0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]', '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1]'),
                               labels = c('0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100'))

# CODE variant status
preVCFgenie_FAOgt0$variant_status <- as.character(NA)
preVCFgenie_FAOgt0[near(preVCFgenie_FAOgt0$PERCENT_OVERLAP, 1/3, tol = .1), ]$variant_status <- "Error"
preVCFgenie_FAOgt0[near(preVCFgenie_FAOgt0$PERCENT_OVERLAP, 2/3, tol = .1), ]$variant_status <- "Ambig"
preVCFgenie_FAOgt0[near(preVCFgenie_FAOgt0$PERCENT_OVERLAP, 3/3, tol = .1), ]$variant_status <- "True"

# FACTOR
preVCFgenie_FAOgt0$variant_status <- factor(preVCFgenie_FAOgt0$variant_status, levels = c('True', 'Ambig', 'Error'))
preVCFgenie_FAOgt0$pass_VCFgenie <- factor(preVCFgenie_FAOgt0$pass_VCFgenie, levels = c('PASS', 'FAIL'))

# Bin counts
preVCFgenie_FAOgt0 %>% group_by(AF_binned) %>%
 summarise(
  n = n()
 )


