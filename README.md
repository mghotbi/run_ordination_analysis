# run_ordination_analysis
A flexible ordination and PERMANOVA analysis pipeline for phyloseq objects, supporting custom grouping, interactions, and automatic sample filtering.

# Ordination and PERMANOVA Pipeline

This function provides a flexible ordination and PERMANOVA analysis pipeline for phyloseq objects, supporting custom grouping, interactions, and automatic sample filtering. The pipeline is fully adaptable to any gene marker or community dataset.

## How to Use

1. **Download or Source the Function**
   
You can source the R script directly from your local folder:

```r
source("run_ordination_analysis.R")

```

```r

source("https://raw.githubusercontent.com/mghotbi/run_ordination_analysis/Rhizosphere-nitrogen-fate/run_ordination_analysis.R")

```


```r
result <- run_ordination_analysis(
  ps_obj = ps_16s,
  marker_label = "Prokaryotic Community",
  group_var = "Treatment",
  interaction_var = "Genotype",
  replicate_var = "Replicate",
  treatment_levels = c("Control", "Flooding", "Herbivory", "DualStress"))

  result$plot
 print(result$permanova)

```
