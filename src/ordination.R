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

# phyloseq_file <- "output/040_phyloseq/ps_filtered.Rds"

########
# MAIN #
########

ps <- readRDS(phyloseq_file)

# run ordination
dca <- ordinate(ps)
dca_dt <- data.table(dca$rproj, keep.rownames = TRUE)

# add sample data
sd <- data.table(data.frame(sample_data(ps)),
                 keep.rownames = TRUE)
pd_wide <- merge(dca_dt, sd)
pd <- melt(pd_wide, id.vars = c("rn", "Location", "Parasitoid", "Pasture"),
           variable.name = "ordinate")

# deal with poa, since it may affect the ordination
pd[, xlab := ifelse(Pasture == "Poa",
            paste0(Location,
                   " (", Pasture, ")"),
            Location)]
pd[, bpgrp := paste(xlab, Parasitoid, sep = "_")]

# plot
gp <- ggplot(pd, aes(y = value,
               x = xlab,
               colour = Parasitoid)) +
    theme(strip.placement = "outside",
          strip.background = element_blank()) +
    facet_wrap(~ ordinate, strip.position = "left") +
    scale_colour_brewer(palette = "Set1") +
    xlab(NULL) + ylab(NULL) +
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
               shape = 16) 

ggsave(plot_file,
       gp,
       device = cairo_pdf, 
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()



