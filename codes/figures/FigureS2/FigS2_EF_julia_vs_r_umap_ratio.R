#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)

df <- read_tsv("data/tables/benchmark/benchmark_summary_primary.tsv", col_types = cols(.default = col_character()))

df_plot <- df %>%
  filter(method == "Julia" | (method == "sclens" & graph == "umap")) %>%
  mutate(
    Consistent_Clustering_Ratio = as.numeric(Consistent_Clustering_Ratio),
    method = case_when(method == "Julia" ~ "scICER (Julia UMAP)", method == "sclens" & graph == "umap" ~ "scICER (R UMAP)"),
    method = factor(method, levels = c("scICER (Julia UMAP)", "scICER (R UMAP)"))
  )

ggplot(df_plot, aes(x = method, y = Consistent_Clustering_Ratio)) +
  stat_summary(aes(color = method), fun = mean, geom = "bar", fill = NA, width = 0.6, linewidth = 0.8) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black", linewidth = 0.8) +
  geom_point(color = "gray50", position = position_jitter(width = 0.15), alpha = 0.6, size = 0.6) +
  scale_color_manual(values = c("scICER (Julia UMAP)" = "#925E9FFF", "scICER (R UMAP)" = "#00A087FF"), name = NULL) +
  expand_limits(y = 0) +
  labs(x = NULL, y = "Consistent Clustering Ratio") +
  theme_classic(base_size = 14) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.length = unit(-0.15, "cm"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "top", legend.key = element_blank())
