#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)
source("codes/figures/utils.R")

df <- read_tsv("data/tables/benchmark/benchmark_summary_primary.tsv", col_types = cols(.default = col_character()))

df_combined <- df %>%
  mutate(Optimal_Cluster_by_ECS = as.numeric(Optimal_Cluster_by_ECS)) %>%
  filter(method %in% c("Paper", "umap_run_in_julia")) %>%
  transmute(sample, sample_type, source = case_when(method == "Paper" ~ "scICE (original)", method == "umap_run_in_julia" ~ "scICE (reproduced)"), optimal = Optimal_Cluster_by_ECS, consistent_cluster)

df_ranges <- df_combined %>% rowwise() %>% mutate(ranges = list(parse_clusters_to_ranges(consistent_cluster))) %>% unnest(ranges) %>% ungroup()

plot_list <- list()
for (stype in unique(df_combined$sample_type)) {
  df_sub <- filter(df_combined, sample_type == stype)
  range_sub <- filter(df_ranges, sample_type == stype)
  sample_levels <- unique(df_sub$sample)

  df_sub <- df_sub %>% mutate(sample = factor(sample, levels = sample_levels), sample_num = as.numeric(sample), x_pos = ifelse(source == "scICE (original)", sample_num - 0.15, sample_num + 0.15))
  range_sub <- range_sub %>% mutate(sample = factor(sample, levels = sample_levels), sample_num = as.numeric(sample), x_pos = ifelse(source == "scICE (original)", sample_num - 0.15, sample_num + 0.15), end_adj = ifelse(start == end, start + 0.1, end), line_color = source)
  rect_data <- df_sub %>% distinct(sample, sample_num) %>% mutate(xmin = sample_num - 0.5, xmax = sample_num + 0.5)

  plot_list[[stype]] <- ggplot() +
    geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = NA, color = "#ADB6B6FF", linewidth = 0.6) +
    geom_segment(data = range_sub, aes(x = x_pos, xend = x_pos, y = start, yend = end_adj, color = line_color), linewidth = 0.7, arrow = arrow(type = "closed", length = unit(0.05, "inches"), ends = "both")) +
    geom_point(data = filter(df_sub, !is.na(optimal)), aes(x = x_pos, y = optimal), size = 2, color = "black") +
    scale_color_manual(values = c("scICE (original)" = "#F39B7FFF", "scICE (reproduced)" = "#4DBBD5FF"), name = "Method") +
    scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 2)) +
    scale_x_continuous(breaks = seq_along(sample_levels), labels = sample_levels, expand = expansion(mult = c(0, 0))) +
    labs(title = stype, y = "Number of Clusters", x = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), legend.position = "right")
}

wrap_plots(plot_list, ncol = 1, guides = "collect")
