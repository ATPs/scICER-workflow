#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(scales)
library(patchwork)

runtime_file <- "data/tables/runtime/scaling/runtime_summary_scicer.tsv"
reference_file <- "data/tables/runtime/scaling/runtime_summary_scice_reference.tsv"

df <- read_tsv(runtime_file, col_types = cols()) %>%
  mutate(Ncells = as.numeric(gsub(",", "", as.character(Ncells))), threads = as.integer(as.character(threads)))

df_reference <- read_tsv(reference_file, col_types = cols()) %>%
  mutate(Ncells = as.numeric(gsub(",", "", as.character(Ncells))), threads = as.integer(as.character(threads))) %>%
  select(sample, Ncells, threads, scice_total_time_sec)

df <- df %>%
  left_join(df_reference, by = c("sample_name" = "sample", "Ncells", "threads")) %>%
  mutate(threads = factor(threads, levels = sort(unique(threads)))) %>%
  arrange(threads, Ncells)

method_colors <- c("scICE" = "#4DBBD5FF", "scICER" = "#E377C2FF")
base_theme <- theme_classic(base_size = 14) + theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.4), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), strip.background = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 12), plot.title = element_text(size = 12), axis.line = element_line(colour = "black", linewidth = 0), axis.ticks = element_line(colour = "black"), legend.title = element_blank())
facet_threads <- facet_wrap(~ threads, nrow = 1, scales = "free_y")

p_time <- df %>%
  select(sample_name, Ncells, threads, scICE = scice_total_time_sec, scICER = scicer_total_time_sec) %>%
  pivot_longer(cols = c(scICE, scICER), names_to = "Method", values_to = "Time_sec") %>%
  ggplot(aes(x = Ncells, y = Time_sec, color = Method, group = Method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 3) +
  facet_threads +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = c(5000, 10000), labels = comma) +
  labs(x = "Number of Cells", y = "Total Time (seconds)", title = "Total Time vs Ncells") +
  base_theme

p_mem <- df %>%
  select(sample_name, Ncells, threads, scICE = scice_peak_memory_mb, scICER = scicer_peak_memory_mb) %>%
  pivot_longer(cols = c(scICE, scICER), names_to = "Method", values_to = "Memory_MB") %>%
  mutate(Memory_GB = Memory_MB / 1024) %>%
  ggplot(aes(x = Ncells, y = Memory_GB, color = Method, group = Method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 3) +
  facet_threads +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = c(5000, 10000), labels = comma) +
  labs(x = "Number of Cells", y = "Peak Memory (GB)", title = "Peak Memory vs Ncells") +
  base_theme

p_time / p_mem
