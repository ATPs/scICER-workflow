#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(scales)
library(patchwork)

trace_dir <- "data/tables/benchmark/resource_traces"
samples <- c("ZhengMix_785", "ZhengMix_2410", "ZhengMix8eq", "ZhengMix_7993", "ZhengMix_12073")
ncells <- c(785, 2410, 4000, 7993, 12073)
threads <- c(1, 5, 10, 20, 40)

read_trace_summary <- function(method, thread, sample) {
  file <- file.path(trace_dir, method, paste0("t", thread, "_", sample, ".txt"))
  if (!file.exists(file)) stop("Missing resource trace: ", file)

  lines <- readLines(file, warn = FALSE)
  summary_lines <- grep("^(Total_Time_Seconds|Peak_CPU_Percent|Average_CPU_Percent|Peak_PSS_MiB|Peak_PSS_MB|Return_Code|Interrupted|PSS_Complete_Samples|Complete_PSS_Samples|PSS_Incomplete_Samples|Incomplete_PSS_Samples)\\t", lines, value = TRUE)
  if (!length(summary_lines)) stop("Missing summary statistics: ", file)

  summary <- read_tsv(I(paste(summary_lines, collapse = "\n")), col_names = c("key", "value"), show_col_types = FALSE)
  get_value <- function(keys, required = TRUE) {
    value <- summary$value[match(keys, summary$key)]
    value <- value[!is.na(value)]
    if (!length(value)) {
      if (required) stop("Missing summary field ", paste(keys, collapse = "/"), ": ", file)
      return(NA_real_)
    }
    suppressWarnings(as.numeric(value[1]))
  }

  return_code <- get_value("Return_Code", required = FALSE)
  if (!is.na(return_code) && return_code != 0) stop("Benchmark failed: ", file)
  interrupted <- summary$value[match("Interrupted", summary$key)]
  if (!is.na(interrupted) && tolower(interrupted) == "true") stop("Benchmark was interrupted: ", file)

  peak_memory <- get_value(c("Peak_PSS_MiB", "Peak_PSS_MB"))
  complete_samples <- get_value(c("PSS_Complete_Samples", "Complete_PSS_Samples"), required = FALSE)
  if (!is.na(complete_samples) && complete_samples == 0) stop("No complete PSS samples: ", file)

  tibble(
    sample_name = sample,
    Ncells = ncells[match(sample, samples)],
    threads = thread,
    Method = method,
    Time_sec = get_value("Total_Time_Seconds"),
    Peak_CPU = get_value("Peak_CPU_Percent"),
    Average_CPU = get_value("Average_CPU_Percent"),
    Memory_MB = peak_memory
  )
}

df <- expand.grid(
  Method = c("scICE", "scICER"),
  threads = threads,
  sample_name = samples,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  do(read_trace_summary(.$Method, .$threads, .$sample_name)) %>%
  ungroup() %>%
  mutate(threads = as.integer(as.character(threads)), Memory_GB = Memory_MB / 1024) %>%
  arrange(threads, Ncells, Method)

method_colors <- c("scICE" = "#4DBBD5FF", "scICER" = "#E377C2FF")
base_theme <- theme_classic(base_size = 14) +
  theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        axis.line = element_line(colour = "black", linewidth = 0),
        axis.ticks = element_line(colour = "black"),
        legend.title = element_blank())
facet_threads <- facet_wrap(~ threads, nrow = 1, scales = "free_y")

make_plot <- function(value, y_label, title) {
  df %>%
    ggplot(aes(x = Ncells, y = .data[[value]], color = Method, group = Method)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 3) +
    facet_threads +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = c(5000, 10000), labels = comma) +
    labs(x = "Number of Cells", y = y_label, title = title) +
    base_theme
}

p_time <- make_plot("Time_sec", "Total Time (seconds)", "Total Time vs Ncells")
p_avg_cpu <- make_plot("Average_CPU", "Average Core CPU (%)", "Average Core CPU vs Ncells")
p_mem <- make_plot("Memory_GB", "Peak Memory (GB)", "Peak Memory vs Ncells")

p_time / p_avg_cpu / p_mem
