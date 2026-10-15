# Benchmark Resource Traces

This directory contains resource-usage traces for scICE (Julia) and scICER (R) benchmark runs. Each file records periodic CPU and memory measurements, followed by summary statistics for the run.

## File names

Files are organized by implementation and named `t<workers>_<dataset>.txt`. For example, `scICE/t10_ZhengMix_785.txt` contains measurements from an scICE run on the `ZhengMix_785` dataset using 10 workers. The worker counts in this collection are 1, 5, 10, 20, and 40.

## Measurements

The tab-separated table at the start of each file contains one row per sample. Memory measurements are in MiB.

| Column | Description |
| --- | --- |
| `Timestamp_Seconds` | Elapsed time since the monitored benchmark phase began, in seconds. |
| `CPU_Percent` | CPU usage as a percentage of one core. For example, 200% corresponds to two fully utilized cores. |
| `CPU_Cores` | CPU usage expressed as the equivalent number of fully utilized cores. |
| `Memory_Usage_MiB` | Primary memory measurement for the sample; these traces use PSS. |
| `Memory_Metric` | Memory metric used for `Memory_Usage_MiB`; `pss` denotes a complete PSS measurement and `pss_incomplete` an incomplete one. |
| `Memory_PSS_MiB` | Proportional set size (PSS) of the monitored process tree. Shared memory pages are apportioned among processes. |
| `PSS_Complete` | Indicates whether PSS was obtained for all monitored processes. |
| `PSS_Missing_Process_Count` | Number of monitored processes without an available PSS measurement. |
| `RSS_Sum_MiB` | Sum of resident set size (RSS) across monitored processes. Shared pages may be counted more than once. |
| `RSS_Complete` | Indicates whether RSS was obtained for all monitored processes. |
| `RSS_Missing_Process_Count` | Number of monitored processes without an available RSS measurement. |
| `Process_Count` | Number of processes observed in the monitored process tree. |
| `PIDs` | Comma-separated process IDs observed in the sample. |

`NA` indicates that a measurement is unavailable. PSS completeness is recorded in `PSS_Complete` and `PSS_Missing_Process_Count`.

## Run summary

The summary at the end of each file reports the run duration and overall resource usage.

| Field | Description |
| --- | --- |
| `Total_Time_Seconds` | Duration of the monitored benchmark phase, in seconds. |
| `Peak_CPU_Percent` / `Peak_CPU_Cores` | Maximum CPU usage as a percentage and in core equivalents. |
| `Average_CPU_Percent` / `Average_CPU_Cores` | Time-weighted mean CPU usage as a percentage and in core equivalents. |
| `Peak_PSS_MiB` / `Average_PSS_MiB` | Peak and time-weighted mean PSS in MiB. |
| `Measurements_Taken` | Number of samples collected. |
| `PSS_Complete_Samples` / `PSS_Incomplete_Samples` | Number of samples with complete and incomplete PSS measurements. |
| `Peak_RSS_Sum_MB` / `Average_RSS_Sum_MB` | Peak and time-weighted mean summed RSS, when available. |
| `Return_Code` | Exit code of the benchmark command; 0 indicates successful completion. |
| `Interrupted` | Indicates whether the run was interrupted, when reported. |

Some files use `MB` instead of `MiB` in summary field names, or use `Complete_PSS_Samples` and `Incomplete_PSS_Samples` for the PSS sample counts.
