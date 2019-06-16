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

# subset data
ps_pruned <- prune_taxa(taxa_sums(ps) > 0, ps)
sd <- data.table(data.frame(sample_data(ps_pruned)),
                 keep.rownames = TRUE)

# calculate richness
chao <- data.table(estimate_richness(ps_pruned,
                                     measures = "Chao1"),
                   keep.rownames = TRUE)

richness_pd <- merge(chao, sd, all.x = TRUE, all.y = FALSE)
richness_pd[, bpgrp := paste(Location, Parasitoid, sep = "_")]

# plot
gp <- ggplot(richness_pd, aes(y = Chao1,
                        x = Location,
                        colour = Parasitoid,
                        ymin = Chao1 - se.chao1,
                        ymax = Chao1 + se.chao1)) +
    theme(strip.placement = "outside",
          strip.background = element_blank()) +
    scale_colour_brewer(palette = "Set1") +
    xlab(NULL) + ylab("Alpha diversity (Chao)") +
    geom_boxplot(fill = NA,
                 aes(group = bpgrp),
                 width = 0.4,
                 outlier.size = -1,
                 colour = alpha("black", 0.5),
                 position = position_dodge(width = 0.75),
                 alpha = 0.25, show.legend = FALSE) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                               dodge.width = 0.75,
                                               seed = 14),
               size = 2,
               alpha = 0.75,
               shape = 16) +
    geom_linerange(position = position_jitterdodge(jitter.width = 0.2,
                                                   dodge.width = 0.75,
                                                   seed = 14),
                   alpha = 0.75)


ggsave(plot_file,
       gp,
       device = cairo_pdf, 
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()

