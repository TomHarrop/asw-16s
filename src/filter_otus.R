#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(phyloseq)
library(data.table)

###########
# GLOBALS #
###########

phyloseq_file <- snakemake@input[["phyloseq"]]
filtered_phyloseq_file <- snakemake@output[["phyloseq"]]

########
# MAIN #
########

ps <- readRDS("output/040_phyloseq/ps.Rds")

#  
otu_wide <- data.table(otu_table(ps), keep.rownames = TRUE)
setnames(otu_wide, "rn", "my_taxid")
otu <- melt(otu_wide,
            id.vars = "my_taxid",
            variable.name = "indiv",
            value.name = "reads")
number_of_indivs <- otu[, length(unique(indiv))]


###########
# FILTERS #
###########

# (i) were present in at least one sample at a relative abundance > 1% of the
# total sequences of that sample or (ii) were present in at least 2% of samples
# at a relative abundance > 0.1% for a given sample, or (iii) were present in at
# least 5% of samples at any abundance level.

# otus at or above 3.5% abundance
overall_community_thr <- 0.035

# were present in at least 2% of samples at a relative abundance > 0.1% for a
# given sample
indivs_2pc <- as.integer(                   
    round(max(2, 0.02 * number_of_indivs), 0))

# relative abundance > 0.1%  for indivs_2pc
abundance_2pc_thr <- 0.001

# present in at least 5% of samples at any abundance level.
indivs_5pc <- as.integer(
    round(max(2, 0.05 * number_of_indivs), 0))


####################
# COUNT AND FILTER #
####################

# count total reads per sample
otu[, total_indiv_reads := sum(reads), by = indiv]

# count total reads per OTU
otu[, total_otu_reads := sum(reads), by = my_taxid]

# calculate fractions
otu[, relative_otu_abundance_in_indiv := reads / total_indiv_reads]

# run filters
pass_community_thr <- otu[, .(
    pass_community_thr = any(relative_otu_abundance_in_indiv >= overall_community_thr)),
    by = my_taxid]

pass_2pc_thr <- otu[, .(
    pass_2pc_thr = sum(relative_otu_abundance_in_indiv > abundance_2pc_thr) >= indivs_2pc),
    by = my_taxid]

pass_5pc_thr <- otu[, reads > 0, by = .(my_taxid, indiv)][
    , .(pass_5pc_thr = sum(V1) >= indivs_5pc), by = my_taxid]

all_filters <- merge(merge(pass_community_thr, pass_2pc_thr),
      pass_5pc_thr)
keep_taxids <- all_filters[, any(pass_community_thr, pass_2pc_thr, pass_5pc_thr),
            by = my_taxid][V1 == TRUE, unique(my_taxid)]

# generate count and abundance tables
filtered_otu <- otu[my_taxid %in% keep_taxids]

filtered_otu_table <- dcast(filtered_otu, my_taxid ~ indiv, value.var = "reads")

new_otu <- otu_table(as.matrix(data.frame(filtered_otu_table, row.names = "my_taxid")),
          taxa_are_rows = TRUE)

# modify the phyloseq object
otu_table(ps) <- new_otu
tax_table(ps) <- tax_table(ps)[keep_taxids]

# write output
saveRDS(ps, filtered_phyloseq_file)

# log
sessionInfo()
