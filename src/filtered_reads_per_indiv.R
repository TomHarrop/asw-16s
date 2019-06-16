#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(phyloseq)
library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

phyloseq_file <- snakemake@input[["phyloseq"]]
plot_file <- snakemake@output[["plot"]]

########
# MAIN #
########

ps <- readRDS(phyloseq_file)

# dev
# ps <- readRDS("output/040_phyloseq/ps_filtered.Rds")

sd <- data.table(data.frame(sample_data(ps)), keep.rownames = TRUE)

pd_wide <- data.table(otu_table(ps), keep.rownames = TRUE)
pd <- melt(pd_wide, id.vars = "rn", variable.name = "indiv", value.name = "reads")

# total counts in passed OTUs
pd_sum <- pd[, .(`Reads per  individual` = sum(reads)), by = indiv]
pd_sum[, pop := gsub("[[:digit:]]+", "", indiv)]

# otus with reads
otu_sum <- pd[, .(`OTUs detected` = sum(reads > 0)), by = indiv]
otu_sum[, pop := gsub("[[:digit:]]+", "", indiv)]

# plot together
merge_sums <- merge(otu_sum, pd_sum)
comb_plot <- melt(merge_sums, id.vars = c("indiv", "pop"))


pd_para <-  merge(sd[, .(indiv = rn, Location, Pasture, Parasitoid)],
               comb_plot)
pd_para[, bpgrp := paste(pop, Parasitoid, sep = "_")]

gp <- ggplot(pd_para, aes(y = value, x = pop, colour = Parasitoid)) +
    theme(strip.placement = "outside",
          strip.background = element_blank()) +
    scale_colour_brewer(palette = "Set1") +
    facet_grid(variable ~ ., scales = "free_y",
               switch = "y") +
    xlab("Population") + ylab(NULL) +
    geom_boxplot(fill = NA,
                 aes(group = bpgrp),
                 width = 0.4,
                 outlier.size = -1,
                 colour = alpha("black", 0.5),
                 position = position_dodge(width = 0.75),
                 alpha = 0.25, show.legend = FALSE) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                               dodge.width = 0.75),
               size = 2,
               alpha = 0.75,
               shape = 16)

ggsave(plot_file,
       gp,
       device = cairo_pdf, 
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()
