#!/usr/bin/env Rscript

library(tidyverse)
library(scales)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
stats_tsv <- args[1]
out_dir <- args[2]

# Create output directory if it doesn't exist
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Loading file: ", stats_tsv)
df <- read_delim(stats_tsv, delim = "\t", col_names = T)

# Make BaseQRankSum, MQRankSum, and ReadPosRankSum numeric
df <- df %>%
  mutate(BaseQRankSum = as.numeric(BaseQRankSum),
         MQRankSum = as.numeric(MQRankSum),
         ReadPosRankSum = as.numeric(ReadPosRankSum))

message("Plotting QUAL")
# Plot density plot per CHROM for QUAL
p_qual <- df %>%
  ggplot(aes(x = QUAL, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "QUAL_density_per_chr.png"), 
       p_qual, width = 8, height = 6)


message("Plotting AN")
# Plot density plot per CHROM for AN
p_AN <- df %>%
  ggplot(aes(x = AN, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AN_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting AC")
# Plot density plot per CHROM for AC
p_AC <- df %>%
  ggplot(aes(x = AC, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AC_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting AF")
# Plot density plot per CHROM for AF
p_AF <- df %>%
  ggplot(aes(x = AF, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AF_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting DP")
# Plot density plot per CHROM for DP
p_DP <- df %>%
  ggplot(aes(x = DP, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "DP_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting QD")
# Plot density plot per CHROM for QD
p_QD <- df %>%
  ggplot(aes(x = QD, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "QD_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting MQ")
# Plot density plot per CHROM for MQ
p_MQ <- df %>%
  ggplot(aes(x = MQ, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "MQ_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting FS")
# Plot density plot per CHROM for FS
p_FS <- df %>%
  ggplot(aes(x = FS, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "FS_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting SOR")
# Plot density plot per CHROM for SOR
p_SOR <- df %>%
  ggplot(aes(x = SOR, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "SOR_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting ExcessHEt")
# Plot density plot per CHROM for ExcessHet
p_ExcessHet <- df %>%
  ggplot(aes(x = ExcessHet, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "ExcessHet_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting InbreedingCoeff")
# Plot density plot per CHROM for InbreedingCoeff
p_InbreedingCoeff <- df %>%
  ggplot(aes(x = InbreedingCoeff, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "InbreedingCoeff_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting BaseQRankSum")
# Plot density plot per CHROM for BaseQRankSum
p_BaseQRankSum <- df %>% 
  filter(!is.na(BaseQRankSum)) %>%
  ggplot(aes(x = BaseQRankSum, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "BaseQRankSum_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting MQRankSum")
# Plot density plot per CHROM for MQRankSum
p_MQRankSum <- df %>% 
  filter(!is.na(MQRankSum)) %>%
  ggplot(aes(x = MQRankSum, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "MQRankSum_density_per_chr.png"),
       p_AN, width = 8, height = 6)

message("Plotting ReadPosRankSum")
# Plot density plot per CHROM for ReadPosRankSum
p_ReadPosRankSum <- df %>% 
  filter(!is.na(ReadPosRankSum)) %>%
  ggplot(aes(x = ReadPosRankSum, fill = CHROM)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "ReadPosRankSum_density_per_chr.png"),
       p_AN, width = 8, height = 6)

       
