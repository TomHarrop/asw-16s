#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

options(nwarnings = 10000)

library(data.table)
library(phyloseq)

#############
# FUNCTIONS #
#############


# Parse the metaphlan into a taxonomy. WARNING - bracken and kraken both return
# duplicate taxon names with the same taxon level for some taxa. I'm just
# arbitrarily picking the first one.
SplitTaxLine <- function(tax_line){
    # define prefixes
    prefix_resolver <- structure(list(
        prefix = c("d", "p", "c", "o", "f", "g", "s"), 
        level = c("domain", "phylum", "class", "order", "family", 
                  "genus", "species")),
        row.names = c(NA, -7L),
        class = c("data.table", 
                  "data.frame"))
    
    # split the taxline
    tax_levels <- data.table(unlist(strsplit(tax_line, "|", fixed = TRUE)))
    tax_dt <- tax_levels[, tstrsplit(V1,
                                     "_",
                                     fixed = TRUE,
                                     names = c("prefix", "taxon"))]
    
    # deal with duplicate levels here (just pick the first one)
    if(any(duplicated(tax_dt, by = "prefix"))) {
        warning(
            paste("Duplicated taxa for tax_line",
                  tax_line,
                  "Arbitrarily picking the first on",
                  sep = "\n"))
        tax_dt <- unique(tax_dt, by = "prefix")
    }
    # generate a wide table of levels
    tax_wide <- merge(prefix_resolver,
                      tax_dt,
                      all.x = TRUE)
    tax_wide[, tax_line := tax_line]
    my_tax <- dcast(tax_wide, tax_line ~ level, value.var = "taxon")
    setcolorder(my_tax, c("tax_line", prefix_resolver[, unique(level)]))
    return(my_tax)
}

###########
# GLOBALS #
###########

kraken_files <- snakemake@input[["mpa_files"]]
sample_catalog_file <- snakemake@input[["sample_catalog"]]
indiv_mapping_file <- snakemake@input[["indiv_mapping"]]
phyloseq_file <- snakemake@output[["phyloseq"]]
parsed_data_file <- snakemake@output[["parsed_data"]]

# dev
# kraken_files <- list.files("output/030_bracken",
#                            recursive = TRUE,
#                            full.names = TRUE,
#                            pattern = "bracken_report_mpa.txt")
# sample_catalog_file <- "data/sample_catalog.csv"
# indiv_mapping_file <- "data/ogbf_sample_info.csv"


########
# MAIN #
########

# read files
names(kraken_files) <- basename(dirname(kraken_files))
kraken_data_list <- lapply(kraken_files, fread, col.names = c("tax", "reads"))
kraken_data <- rbindlist(kraken_data_list, idcol = "indiv")

# generate a taxid for the analysis & parse the taxonomy (very slow)
unique_taxa <- unique(kraken_data[, .(tax)])
unique_taxa[, my_taxid := paste0("tax", 1:.N)]
full_taxonomy <- unique_taxa[, SplitTaxLine(tax), by = my_taxid]
warnings()

# add the tax info
kraken_data_with_tax <- merge(kraken_data,
                              full_taxonomy,
                              all.x = TRUE,
                              all.y = FALSE,
                              by.x = "tax",
                              by.y = "tax_line")

# map the sample names
indiv_mapping <- fread(indiv_mapping_file)
sample_catalog <- fread(sample_catalog_file)
indiv_mapping[, indiv := sprintf("%03d",
                                 as.integer(gsub(".*?([[:digit:]]+)$",
                                                 "\\1",
                                                 `Library ID`)))]
indiv_to_sampleinfo <- merge(indiv_mapping[, .(indiv, Individual = `Library name`)],
                             sample_catalog,
                             all.x = TRUE,
                             all.y = FALSE)
indiv_to_sampleinfo[is.na(Parasitoid) & !is.na(Location), Parasitoid := FALSE]

# merge teh sample info
kraken_data_with_info <- merge(kraken_data_with_tax,
                               indiv_to_sampleinfo[, .(indiv,
                                                       indiv_name = Individual,
                                                       Location,
                                                       Pasture,
                                                       Parasitoid)],
                               by = "indiv")

# subset at the genus level
kraken_level_subset <- unique(kraken_data_with_info[!is.na(genus)],
       by = c("indiv_name", "my_taxid"))


# generate phyloseq objects
sd_wide <- unique(kraken_level_subset[, .(indiv_name,
                                          Location,
                                          Pasture,
                                          Parasitoid)])
sample_data <- data.frame(sd_wide, row.names = "indiv_name")

otu_wide <- dcast(kraken_level_subset,
                  my_taxid ~ indiv_name,
                  value.var = "reads",
                  fill = 0)
otu <- otu_table(as.matrix(data.frame(otu_wide, row.names = "my_taxid")),
                 taxa_are_rows = TRUE)

tax_dt <- unique(kraken_level_subset[, .(my_taxid,
                        domain,
                        phylum,
                        class,
                        order,
                        family,
                        genus)])

tax <- tax_table(as.matrix(data.frame(tax_dt, row.names = "my_taxid")))

# combine
ps <- phyloseq(otu, sample_data, tax)

# write output
fwrite(kraken_data_with_info, parsed_data_file)
saveRDS(ps, phyloseq_file)

# log
sessionInfo()
