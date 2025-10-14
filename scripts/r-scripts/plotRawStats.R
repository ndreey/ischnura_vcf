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
# Plot density plot  for QUAL
p_qual <- df %>%
  ggplot(aes(x = QUAL)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "QUAL_density.png"), 
       p_qual, width = 8, height = 6)


message("Plotting AN")
# Plot density plot  for AN
p_AN <- df %>%
  ggplot(aes(x = AN)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AN_density.png"),
       p_AN, width = 8, height = 6)

message("Plotting AC")
# Plot density plot  for AC
p_AC <- df %>%
  ggplot(aes(x = AC)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AC_density.png"),
       p_AC, width = 8, height = 6)

message("Plotting AF")
# Plot density plot  for AF
p_AF <- df %>%
  ggplot(aes(x = AF)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AF_density.png"),
       p_AF, width = 8, height = 6)

message("Plotting DP")
# Plot density plot  for DP
p_DP <- df %>%
  ggplot(aes(x = DP)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "DP_density.png"),
       p_DP, width = 8, height = 6)

message("Plotting QD")
# Plot density plot  for QD
p_QD <- df %>%
  ggplot(aes(x = QD)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "QD_density.png"),
       p_QD, width = 8, height = 6)

message("Plotting MQ")
# Plot density plot  for MQ
p_MQ <- df %>%
  ggplot(aes(x = MQ)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "MQ_density.png"),
       p_MQ, width = 8, height = 6)

message("Plotting FS")
# Plot density plot  for FS
p_FS <- df %>%
  ggplot(aes(x = FS)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "FS_density.png"),
       p_FS, width = 8, height = 6)

message("Plotting SOR")
# Plot density plot  for SOR
p_SOR <- df %>%
  ggplot(aes(x = SOR)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "SOR_density.png"),
       p_SOR, width = 8, height = 6)

message("Plotting ExcessHEt")
# Plot density plot  for ExcessHet
p_ExcessHet <- df %>%
  ggplot(aes(x = ExcessHet)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "ExcessHet_density.png"),
       p_ExcessHet, width = 8, height = 6)

message("Plotting InbreedingCoeff")
# Plot density plot  for InbreedingCoeff
p_InbreedingCoeff <- df %>%
  ggplot(aes(x = InbreedingCoeff)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "InbreedingCoeff_density.png"),
       p_InbreedingCoeff, width = 8, height = 6)

message("Plotting BaseQRankSum")
# Plot density plot  for BaseQRankSum
p_BaseQRankSum <- df %>% 
  filter(!is.na(BaseQRankSum)) %>%
  ggplot(aes(x = BaseQRankSum)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "BaseQRankSum_density.png"),
       p_BaseQRankSum, width = 8, height = 6)

message("Plotting MQRankSum")
# Plot density plot  for MQRankSum
p_MQRankSum <- df %>% 
  filter(!is.na(MQRankSum)) %>%
  ggplot(aes(x = MQRankSum)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "MQRankSum_density.png"),
       p_MQRankSum, width = 8, height = 6)

message("Plotting ReadPosRankSum")
# Plot density plot  for ReadPosRankSum
p_ReadPosRankSum <- df %>% 
  filter(!is.na(ReadPosRankSum)) %>%
  ggplot(aes(x = ReadPosRankSum)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "ReadPosRankSum_density.png"),
       p_ReadPosRankSum, width = 8, height = 6)

message("Plotting complete")       
