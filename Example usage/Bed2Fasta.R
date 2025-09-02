## ===== CNV codon encoding from a single BED/TSV =====
rm(list = ls())

library(dplyr)
library(readr)
library(stringr)
library(seqinr)

#--------------------- Inputs & parameters --------------------------------------
setwd("set path to bed file")
bed_path    <- "cnv_simulated.bed"  # <-- put your path here
region_size <- 50                   # segment/bin size (bp)
individuals_per_pop <- 10
#--------------------- Load BED/TSV --------------------------------------------
# Expected columns: chr, start, end, copy_number, individual_id, population_id
cnv <- read_tsv(
  bed_path,
  col_types = cols(
    chr            = col_character(),
    start          = col_double(),
    end            = col_double(),
    copy_number    = col_double(),
    individual_id  = col_character(),
    population_id  = col_character()
  )
)

cnv <- cnv %>%
  mutate(chr = as.character(chr))

cnv <- cnv %>%
  mutate(copy_number = pmin(copy_number, 3))

# Basic sets
pop_levels <- unique(cnv$population_id)
pop_number <- length(pop_levels)

# Chromosomes present and their sizes inferred from max end
chr_levels <- cnv %>%
  group_by(chr) %>%
  summarise(chr_size = max(end, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr)

chro_number <- nrow(chr_levels)
chro_sizes  <- chr_levels$chr_size
names(chro_sizes) <- chr_levels$chr
#--------------------- Allocate data matrices -----------------------------------
data_pop     <- array(vector("list", pop_number * chro_number),
                      dim = c(pop_number, chro_number))

# For each population & chromosome, initialize matrix to CN=2
for (p in seq_len(pop_number)) {
  pop_id <- pop_levels[p]
  individual_id <- unique(cnv[cnv$population_id==pop_id, ]$individual_id)
  
  for (c in seq_len(chro_number)) {
    chr_name <- chr_levels$chr[c]
    segment_numbers <- round(chro_sizes[c] / region_size) + 1
    
    # Initialize matrix to 2 for all individuals
    mat <- matrix(2L, nrow = individuals_per_pop, ncol = segment_numbers)
    
    # Subset CNV calls for this (pop, chr)
    sub <- cnv %>%
      filter(population_id == pop_id, chr == chr_name)
    
    # If individual has no variant, mat row remains 2 by default
    if(nrow(sub) > 0){
      for (j in seq_len(nrow(sub))) {
        start_seg <- round(sub$start[j] / region_size) + 1
        end_seg   <- round(sub$end[j] / region_size)
        if (end_seg < start_seg) end_seg <- start_seg
        
        individual_index <- match( sub$individual_id[j], individual_id)
        mat[individual_index, start_seg:end_seg] <- as.integer(sub$copy_number[j])
      }
    } # sub
    data_pop[[p, c]] <- mat
  } # chr
} # pop
#--------------------- Extract biallelic sites ----------------------------------
# Keep sites where, across populations, each population has <3 unique copy numbers
# (i.e., at most 2 unique copy numbers per pop).

# Pre-allocate a working list for all pops
tmp_holder <- vector("list", pop_number)

for (c in seq_len(chro_number)) {
  n_seg <- ncol(data_pop[[1, c]])
  
  for (site_id in seq_len(n_seg)) {
    uniq_counts <- sapply(seq_len(pop_number), function(p) length(unique(data_pop[[p, c]][, site_id])))
    
    if (max(uniq_counts) < 3) {  # keep mono- and biallelic only
      for (p in seq_len(pop_number)) {
        tmp_holder[[p]] <- cbind(tmp_holder[[p]], data_pop[[p, c]][, site_id])
      }
    }
  }
} # chr
#--------------------- Codon encoding & FASTA output ----------------------------
codons <- c(
  "aaa","aac","aag","aat","aca","acc","acg","act","aga","agc","agg","agt",
  "ata","atc","atg","att","caa","cac","cag","cat","cca","ccc","ccg","cct",
  "cga","cgc","cgg","cgt","cta","ctc","ctg","ctt","gaa","gac","gag","gat",
  "gca","gcc","gcg","gct","gga","ggc","ggg","ggt","gta","gtc","gtg","gtt",
  "taa","tac","tag","tat","tca","tcc","tcg","tct","tga","tgc"
)

copy_number_alphabet <- c(0, 1, 2, 3)
siallelic_combination <- combn(copy_number_alphabet, 1)
biallelic_combination <- combn(copy_number_alphabet, 2)

# Initialize sequences
seq_by_pop <- lapply(seq_len(pop_number), function(p) character())

# Iterate over populations
for (p in seq_len(pop_number)) {
  n_sites <- ncol(tmp_holder[[p]])
  
  for (site_id in seq_len(n_sites)) {
    #site_id = 26
    pop_vec <- as.integer(tmp_holder[[p]][, site_id])
    alleles_tmp <- sort(unique(pop_vec))
    
    # Encode monallelic
    if (length(alleles_tmp) == 1) {
      codon_index <- which(siallelic_combination == alleles_tmp)
      
      # Encode biallelic
    } else if (length(alleles_tmp) == 2) {
      allele_index <- which(apply(biallelic_combination, 2, function(col) all(col == alleles_tmp)))
      count_1 <- sum(pop_vec == alleles_tmp[1])
      codon_index <- length(copy_number_alphabet) + (allele_index - 1) * (length(pop_vec) - 1) + count_1
    } else {
      next  # skip if more than 2 alleles (should not happen)
    }
    
    seq_by_pop[[p]] <- c(seq_by_pop[[p]], codons[codon_index])
  }
}

# Loop over populations
for (pop_id in seq_len(pop_number)) {
  # Skip if no codons
  if (length(seq_by_pop[[pop_id]]) == 0) next
  
  # Concatenate all codons for this population
  seq2fasta <- paste(seq_by_pop[[pop_id]], collapse = "")
  
  # Write a single FASTA file
  if(pop_id == 1){
    write.fasta(seq2fasta, paste0("pop_", pop_id), file.out = paste0("codon_encoded.fasta"), open = "w", nbchar = 60, as.string = FALSE)
  }else{
    write.fasta(seq2fasta, paste0("pop_", pop_id), file.out = paste0("codon_encoded.fasta"), open = "a", nbchar = 60, as.string = FALSE)
  }
}




















