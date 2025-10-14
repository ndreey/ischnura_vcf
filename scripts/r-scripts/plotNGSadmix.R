#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggrepel)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plotNGSadmix.R <combined.tsv> <metadata.tsv> <out_dir>")
}
COMBINED <- args[1]   # e.g. nana.txt  (ID Q1 Q2 ... whitespace-delimited)
OUTDIR  <- args[2]   # e.g. doc/ischnura_metadata.tsv (must have column ID)

# Create output directory if it doesn't exist
output_dir <- OUTDIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===== Arguments =====
K <- 2
POP_ORDER <- c("Bunkeflostrand","Lomma","Borgeby","Genarp","Ilstorp","Vombs Vattenverk")
META_FILE <- "doc/ischnura_metadata.tsv"       # has ID, Population, etc.

# 1) read combined file
df <- read_delim(COMBINED, delim = "\t", col_names=F)
K <- ncol(df) - 1
colnames(df) <- c("ID", paste0("Q", seq_len(K)))

meta <- read_delim(META_FILE, delim = "\t", col_names = TRUE)

# ---- 3) attach IDs + metadata (wide form for ordering) ----
wide <- df %>%
  left_join(meta, by = "ID") %>%
  mutate(Population = factor(Population, levels = POP_ORDER))

# ---- 4) order individuals within each population by Q2 (desc) ----
SORT_CLUSTER <- "Q2"   # change to "Q1" etc. if you prefer
ordered_ids <- wide %>%
  arrange(Population) %>%
  group_by(Population) %>%
  arrange(desc(.data[[SORT_CLUSTER]]), .by_group = TRUE) %>%
  ungroup() %>%
  pull(ID) %>%
  unique()

# 5) long format
long_df <- wide %>%
  pivot_longer(starts_with("Q"), names_to = "Cluster", values_to = "Proportion") %>%
  mutate(ID = factor(ID, levels = rev(ordered_ids)))
# ---- 6) main plot (no grey facet strip) ----
p_admix <- ggplot(long_df, aes(ID, Proportion, fill = Cluster)) +
  geom_col(width = 0.99) +
  facet_grid(~ Population, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),  # removes grey rectangle
    strip.text = element_text(size = 14)
  )

label_plot <- ggdraw() + draw_label(paste0("K = ", K), angle = 270, size = 14)
final_plot <- plot_grid(p_admix, label_plot, ncol = 2, rel_widths = c(1, 0.03))

ggsave(file.path(output_dir, "admixture.png"), final_plot, 
    width = 14, height = 4, dpi = 300)