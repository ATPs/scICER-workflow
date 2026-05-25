#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)
source("codes/figures/utils.R")

df <- read_tsv("data/tables/benchmark/benchmark_summary_workflow_graph.tsv", col_types = cols(.default = col_character()))

df_plot <- df %>%
  mutate(Optimal_Cluster_by_ECS = as.numeric(Optimal_Cluster_by_ECS)) %>%
  filter(method %in% c("RNA", "SCT", "sclens"), graph %in% c("snn", "knn", "umap")) %>%
  mutate(graph = factor(graph, levels = c("snn", "knn", "umap")), method = factor(method, levels = c("RNA", "SCT", "sclens")))

df_ranges <- df_plot %>% rowwise() %>% mutate(ranges = list(parse_clusters_to_ranges(consistent_cluster))) %>% unnest(ranges) %>% ungroup()
color_palette <- c("snn" = "#1F77B4", "knn" = "#FF7F0E", "umap" = "#2CA02C")
plot_list <- list()

for (current_method in unique(df_plot$method)) {
  for (current_type in unique(df_plot$sample_type)) {
    df_sub <- filter(df_plot, method == current_method, sample_type == current_type)
    range_sub <- filter(df_ranges, method == current_method, sample_type == current_type)
    current_samples <- unique(df_sub$sample)
    if (length(current_samples) == 0) next

    graph_types <- c("snn", "knn", "umap")
    df_sub <- df_sub %>% mutate(sample = factor(sample, levels = current_samples), sample_num = as.numeric(sample), x_pos = sample_num + (match(graph, graph_types) - 2) * 0.15)
    range_sub <- range_sub %>% mutate(sample = factor(sample, levels = current_samples), sample_num = as.numeric(sample), x_pos = sample_num + (match(graph, graph_types) - 2) * 0.15, end_adj = ifelse(start == end, start + 0.1, end))
    rect_data <- df_sub %>% distinct(sample, sample_num) %>% mutate(xmin = sample_num - 0.5, xmax = sample_num + 0.5)

    plot_list[[paste(current_method, current_type, sep = "_")]] <- ggplot() +
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = NA, color = "#ADB6B6FF", linewidth = 0.6) +
      geom_segment(data = range_sub, aes(x = x_pos, xend = x_pos, y = start, yend = end_adj, color = graph), linewidth = 0.7, arrow = arrow(type = "closed", length = unit(0.05, "inches"), ends = "both")) +
      geom_point(data = filter(df_sub, !is.na(Optimal_Cluster_by_ECS)), aes(x = x_pos, y = Optimal_Cluster_by_ECS), size = 2, color = "black") +
      scale_color_manual(values = color_palette, name = "Graph") +
      scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 2)) +
      scale_x_continuous(breaks = seq_along(current_samples), labels = current_samples, expand = expansion(mult = c(0, 0))) +
      labs(title = paste(current_method, "-", current_type), y = "Number of Clusters", x = NULL) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), legend.position = "right", plot.title = element_text(face = "bold"))
  }
}

wrap_plots(plot_list, ncol = 1, guides = "collect")
