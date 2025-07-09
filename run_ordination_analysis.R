# ordination_pipeline.R
# Author: Mitra Ghotbi
# Description:
# A flexible ordination and PERMANOVA pipeline for phyloseq objects.
# Supports dynamic grouping, interactions, and replicates with auto-filtering.
# Required Packages: phyloseq, vegan, ggplot2, dplyr, tibble

#' Run ordination analysis and PERMANOVA
run_ordination_analysis <- function(ps_obj,
                                    marker_label = "Community",
                                    group_var = "Group",
                                    interaction_var = "SubGroup",
                                    replicate_var = "Replicate",
                                    treatment_levels = NULL) {
  set.seed(23465)

  # Extract metadata
  meta_df <- as(sample_data(ps_obj), "data.frame")

  # Optional recoding: Normalize group_var
  meta_df[[group_var]] <- as.character(meta_df[[group_var]])
  meta_df[[group_var]][meta_df[[group_var]] == "FloodingHerbivory"] <- "DualStress"

  # Update phyloseq object
  sample_data(ps_obj)[[group_var]] <- factor(meta_df[[group_var]], levels = treatment_levels)
  sample_data(ps_obj)[[interaction_var]] <- as.factor(meta_df[[interaction_var]])
  sample_data(ps_obj)[[replicate_var]] <- as.factor(meta_df[[replicate_var]])

  # Filter samples with valid group_var
  keep_samples <- !is.na(sample_data(ps_obj)[[group_var]]) &
    sample_data(ps_obj)[[group_var]] != "" &
    sample_data(ps_obj)[[group_var]] != "NA"

  removed_samples <- sum(!keep_samples)
  cat("Removed", removed_samples, "samples due to missing or invalid group labels.\n")

  ps_obj <- prune_samples(keep_samples, ps_obj)

  if (nsamples(ps_obj) == 0) {
    stop("All samples removed after filtering. Please check your input and group_var.")
  }

  # Run ordination and collect results
  result <- compute_ordination_results(
    ps_obj = ps_obj,
    title_label = marker_label,
    group_var = group_var,
    interaction_var = interaction_var,
    replicate_var = replicate_var,
    treatment_levels = treatment_levels
  )

  print(result$plot)

  # Build PERMANOVA formula
  formula_str <- as.formula(paste("dist ~", group_var, "*", interaction_var))

  permanova_result <- vegan::adonis2(
    formula = formula_str,
    data = result$meta,
    strata = result$meta[[replicate_var]],
    permutations = 999,
    by = "terms"
  )

  print(permanova_result)

  return(list(
    plot = result$plot,
    permanova = permanova_result,
    meta = result$meta
  ))
}

#' Compute ordination plot and metadata
compute_ordination_results <- function(ps_obj,
                                       title_label = "Ordination Plot",
                                       group_var = "Group",
                                       interaction_var = "SubGroup",
                                       replicate_var = "Replicate",
                                       treatment_levels = NULL) {
  # Normalize counts (RLE)
  ps_norm <- normalization_set(ps_obj, method = "rle")$dat.normed

  # Ordination
  dist_bc <- phyloseq::distance(ps_norm, method = "bray")
  ord <- ordinate(ps_norm, method = "PCoA", distance = dist_bc)

  eig_vals <- ord$values$Relative_eig[1:2] * 100
  x_lab <- paste0("PCoA 1 (", round(eig_vals[1], 1), "%)")
  y_lab <- paste0("PCoA 2 (", round(eig_vals[2], 1), "%)")

  # Extract plotting data
  ord_plot_data <- plot_ordination(ps_obj, ord, color = group_var)$data

  ord_df <- ord_plot_data |>
    mutate(
      SampleID = rownames(ord_plot_data),
      GroupVar = factor(.data[[group_var]], levels = treatment_levels)
    ) |>
    filter(!is.na(GroupVar)) |>
    droplevels()

  # PERMANOVA metadata
  meta_df <- ord_df |>
    select(SampleID, GroupVar,
           all_of(group_var),
           all_of(interaction_var),
           all_of(replicate_var)) |>
    distinct()

  rownames(meta_df) <- meta_df$SampleID

  dist_matrix <- as.matrix(dist_bc)[rownames(meta_df), rownames(meta_df)]
  stopifnot(all(rownames(meta_df) == rownames(dist_matrix)))

  # Compute group centroids
  centroids <- ord_df |>
    group_by(GroupVar) |>
    summarise(
      PC1 = mean(Axis.1, na.rm = TRUE),
      PC2 = mean(Axis.2, na.rm = TRUE),
      .groups = "drop"
    )

  # Plot
  p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = GroupVar)) +
    geom_point(size = 2.2, alpha = 0.8) +
    stat_ellipse(type = "t", level = 0.95, linewidth = 0.5) +
    geom_point(data = centroids, aes(x = PC1, y = PC2, fill = GroupVar),
               shape = 21, size = 6, stroke = 1.3) +
    geom_text(data = centroids, aes(x = PC1, y = PC2, label = GroupVar),
              vjust = -1.3, fontface = "bold", size = 4, color = "black") +
    scale_color_manual(values = treatment_colors, na.translate = FALSE) +
    scale_fill_manual(values = treatment_colors, na.translate = FALSE) +
    labs(title = title_label, x = x_lab, y = y_lab) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  return(list(
    plot = p,
    dist = dist_matrix,
    meta = meta_df
  ))
}
