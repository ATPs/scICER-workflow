#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

df <- read_tsv(
  "data/tables/benchmark/benchmark_summary_workflow_graph.tsv",
  col_types = cols(.default = col_character())
)

df_plot <- df %>%
  filter(method %in% c("RNA", "SCT", "sclens"), graph %in% c("snn", "knn", "umap")) %>%
  mutate(
    Consistent_Clustering_Ratio = as.numeric(Consistent_Clustering_Ratio),
    graph = factor(graph, levels = c("snn", "knn", "umap")),
    method = factor(method, levels = c("RNA", "SCT", "sclens"))
  )

graph_colors <- c("snn" = "#1F77B4", "knn" = "#FF7F0E", "umap" = "#2CA02C")
method_alpha <- c("RNA" = 1, "SCT" = 0.5, "sclens" = 0)
method_size <- c("RNA" = 0.8, "SCT" = 0.8, "sclens" = 1.2)

p1 <- ggplot(df_plot, aes(x = interaction(method, graph), y = Consistent_Clustering_Ratio)) +
  stat_summary(aes(color = graph, fill = graph, alpha = method, size = method), fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  geom_jitter(aes(color = graph), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.6, size = 0.6) +
  scale_color_manual(values = graph_colors, name = "Graph") +
  scale_fill_manual(values = graph_colors, name = "Graph") +
  scale_alpha_manual(values = method_alpha, name = "Workflow") +
  scale_size_manual(values = method_size, name = "Workflow") +
  labs(x = NULL, y = "Consistent Clustering Ratio") +
  theme_classic(base_size = 14) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "right")

p2 <- df_plot %>%
  mutate(sample_type = factor(sample_type, levels = c("Data with Hierarchical structure", "Data w/o Hierarchical structure"), labels = c("w hier.", "w/o hier."))) %>%
  ggplot(aes(x = sample_type, y = Consistent_Clustering_Ratio)) +
  stat_summary(aes(color = graph, fill = graph, alpha = method, size = method, group = interaction(method, graph)), fun = mean, geom = "bar", position = position_dodge(width = 0.9), width = 0.6) +
  stat_summary(aes(group = interaction(method, graph)), fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", position = position_dodge(width = 0.9), width = 0.15, color = "black") +
  geom_point(aes(color = graph, group = interaction(method, graph)), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), alpha = 0.6, size = 0.6) +
  scale_color_manual(values = graph_colors, name = "Graph") +
  scale_fill_manual(values = graph_colors, name = "Graph") +
  scale_alpha_manual(values = method_alpha, name = "Workflow") +
  scale_size_manual(values = method_size, name = "Workflow", guide = "none") +
  labs(x = NULL, y = "Consistent Clustering Ratio") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "right")

p1 / p2
