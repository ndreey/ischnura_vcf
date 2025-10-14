#!/usr/bin/env Rscript

library(tidyverse)


lmiss <- read_delim("/cfs/klemming/projects/supr/snic2020-6-170/Sofie_N/Make_Primers/makeVFD_with_DToL_2023_11_03/Analysis_DToL_2023-11-03_17-50-13/11_VCF_wrangling/vcf-final/Stats/stats_downsam_50/SwD_vcftools_filt.no-miss_noX_noCHR13.txt.lmiss", delim = "\t")


# Summary statistics
summary_stats <- summary(lmiss$F_MISS)
numbAboveThreshold <- sum(lmiss$F_MISS > 0.9)

# Create density plot
lmiss %>% 
  ggplot(aes(x = F_MISS)) +
  geom_density(fill = "orchid", alpha = 0.5) +
  labs(title = "Density of Missingness (F_MISS)",
       x = "F_MISS",
       y = "Density") +
  annotate(
    "text",
    x = 0.4,  # You can adjust this for a better fit
    y = max(density(lmiss$F_MISS)$y) * 0.5,  # Dynamic Y position based on data density
    label = paste(
      " Min:", round(summary_stats[1], 2), "\n",
      "1st Qu.:", round(summary_stats[2], 2), "\n",
      "Median:", round(summary_stats[3], 2), "\n",
      "Mean:", round(summary_stats[4], 2), "\n",
      "3rd Qu.:", round(summary_stats[5], 2), "\n",
      "Max:", round(summary_stats[6], 2)
    ),
    size = 5,
    color = "black",
    hjust = 0
  ) +
  annotate(
    "text", x = 0.75, y = 4.5, label = paste("Number of variants with F_MISS > 0.9:", numbAboveThreshold),
    size = 4, color = "red", hjust = 0.5
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())

# Save plot
ggsave("missingness_density_plot.png", width = 8, height = 6)
