#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

full_usage_file <- "data/tables/runtime/pancreas/full_usage_memory.tsv"
full_window_file <- "data/tables/runtime/pancreas/full_window.tsv"
diet_usage_file <- "data/tables/runtime/pancreas/diet_usage_memory.tsv"
diet_window_file <- "data/tables/runtime/pancreas/diet_window.tsv"

summarize_memory <- function(usage_file, window_file, method) {
  usage <- read.delim(usage_file, check.names = FALSE)
  window_df <- read.delim(window_file, check.names = FALSE)
  usage_scice <- usage %>% filter(Timestamp_Seconds >= window_df$start[1], Timestamp_Seconds <= window_df$end[1])
  data.frame(Method = method, Peak_memory_GB = max(usage_scice$Memory_Usage_MB, na.rm = TRUE) / 1024)
}

df <- bind_rows(
  summarize_memory(full_usage_file, full_window_file, "Full"),
  summarize_memory(diet_usage_file, diet_window_file, "Diet")
)
df$Method <- factor(df$Method, levels = c("Full", "Diet"))

ggplot(df, aes(Method, Peak_memory_GB, fill = Method)) +
  geom_col(width = 0.55, color = "grey20", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", Peak_memory_GB)), vjust = -0.45, size = 4) +
  scale_fill_manual(values = c(Full = "#3B5B8A", Diet = "#C97C5D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Peak Memory Usage of scICE_clustering", x = NULL, y = "Peak memory (GB)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12), axis.text.x = element_text(face = "plain", size = 12), axis.text.y = element_text(size = 12), axis.line = element_line(color = "grey25", linewidth = 0.45), axis.ticks = element_line(color = "grey25", linewidth = 0.4), panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3), panel.grid.minor = element_blank(), legend.position = "none", plot.margin = margin(10, 15, 10, 15))
