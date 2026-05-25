#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

df <- read_tsv(
  "data/tables/benchmark/benchmark_summary_primary.tsv",
  col_types = cols(.default = col_character())
)

df_plot <- df %>%
  filter(method %in% c("Julia", "umap_run_in_julia")) %>%
  mutate(
    Consistent_Clustering_Ratio = as.numeric(Consistent_Clustering_Ratio),
    method = recode(
      method,
      "Julia" = "scICER (Julia UMAP)",
      "umap_run_in_julia" = "scICE (reproduced)"
    ),
    method = factor(method, levels = c("scICE (reproduced)", "scICER (Julia UMAP)"))
  )

p1 <- ggplot(df_plot, aes(x = method, y = Consistent_Clustering_Ratio)) +
  stat_summary(aes(color = method), fun = mean, geom = "bar", fill = NA, width = 0.6, linewidth = 0.8) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black", linewidth = 0.8) +
  geom_point(color = "gray50", position = position_jitter(width = 0.15), alpha = 0.6, size = 0.6) +
  scale_color_manual(values = c("scICE (reproduced)" = "#4DBBD5FF", "scICER (Julia UMAP)" = "#925E9FFF"), name = NULL) +
  expand_limits(y = 0) +
  labs(x = NULL, y = "Consistent Clustering Ratio") +
  theme_classic(base_size = 14) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.length = unit(-0.15, "cm"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "top", legend.key = element_blank())

p2 <- df_plot %>%
  mutate(sample_type = factor(sample_type, levels = c("Data with Hierarchical structure", "Data w/o Hierarchical structure"), labels = c("w hier.", "w/o hier."))) %>%
  ggplot(aes(x = sample_type, y = Consistent_Clustering_Ratio)) +
  stat_summary(aes(color = method), fun = mean, geom = "bar", fill = NA, position = position_dodge(width = 0.8), width = 0.6, linewidth = 0.8) +
  stat_summary(aes(group = method), fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", color = "black", position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.8) +
  geom_point(color = "gray50", position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.6, size = 0.6) +
  scale_color_manual(values = c("scICE (reproduced)" = "#4DBBD5FF", "scICER (Julia UMAP)" = "#925E9FFF"), name = NULL) +
  labs(x = NULL, y = "Consistent Clustering Ratio") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks.length = unit(-0.15, "cm"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "top", legend.key = element_blank())

p1 / p2
