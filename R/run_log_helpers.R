start_console_log <- function(path, split = FALSE) {
  
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  
  sink(con, split = split)
  sink(con, type = "message")

  # return stopper
  function() {
    try(sink(type = "message"), silent = TRUE) 
    try(sink(), silent = TRUE) 
    try(close(con), silent = TRUE)
  }
}

# Powershell: Get-Content path -Wait

write_run_log_txt <- function(run_log, file) {
  
  fmt <- function(...) sprintf(...)
  
  # Helper function for formatting time
  fmt_time <- function(s) {
    if (is.null(s) || is.na(s)) {
      return("N/A")
    }
    if (s < 60) {
      sprintf("%.1f sec", s)
    } else if (s < 3600) {
      sprintf("%.1f min", s / 60)
    } else if (s < 86400) {
      sprintf("%.2f hrs", s / 3600)
    } else {
      sprintf("%.2f days", s / 86400)
    }
  }
  
  lines <- c(
    "============================================================",
    "MODEL RUN SUMMARY",
    "============================================================",
    "",
    fmt("Run type:        %s", run_log$run_type),
    fmt("Run ID:          %s", run_log$run_id),
    fmt("Output dir:      %s", run_log$out_dir),
    "",
    fmt("Start time:      %s", run_log$start_time),
    fmt("End time:        %s", run_log$end_time)
  )
  
  # Timing section
  if (!is.null(run_log$timing)) {
    lines <- c(
      lines,
      "",
      "Runtime breakdown:",
      fmt("  Data prep:    %s", fmt_time(run_log$timing$data_prep_s)),
      fmt("  Tuning:       %s", fmt_time(run_log$timing$tuning_s)),
      fmt("  CV fitting:   %s", fmt_time(run_log$timing$cv_s)),
      fmt("  Final fit:    %s", fmt_time(run_log$timing$final_fit_s)),
      fmt("  ─────────────────────────────"),
      fmt("  TOTAL:        %s (%.2f hrs)", 
          fmt_time(run_log$timing$total_s),
          run_log$timing$total_s / 3600)
    )
  } else if (!is.null(run_log$runtime_s)) {
    # Fallback for old format
    lines <- c(
      lines,
      fmt("Runtime (sec):   %.1f", run_log$runtime_s),
      fmt("Runtime (hrs):   %.2f", run_log$runtime_s / 3600)
    )
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "INPUT DATA",
    "------------------------------------------------------------",
    fmt("Burn severity:   %s", run_log$input_file),
    fmt("Predictor dir:   %s", run_log$metric_dir),
    fmt("ROI GPKG:        %s", run_log$roi_gpkg),
    fmt("ROI layer:       %s", run_log$roi_layer),
    "",
    fmt("Total samples:   %d", run_log$n_total),
    fmt("Train samples:   %d (%.1f%%)", run_log$n_train, 100 * run_log$n_train / run_log$n_total),
    fmt("Test samples:    %d (%.1f%%)", run_log$n_test, 100 * run_log$n_test / run_log$n_total),
    "",
    "------------------------------------------------------------",
    "CROSS-VALIDATION DESIGN",
    "------------------------------------------------------------",
    fmt("Initial test proportion: %.2f", run_log$test_prop),
    fmt("Spatial buffer (m):      %d", run_log$buffer_dist_m),
    fmt("V-folds:                 %d", run_log$v_folds),
    fmt("Blocks per fold:         %s", ifelse(!is.null(run_log$blocks_per_fold), run_log$blocks_per_fold, "N/A")),
    fmt("Tuning repeats:          %d", run_log$tune_reps),
    fmt("Fit repeats:             %d", run_log$fit_reps),
    fmt("Tuning resamples:        %d", run_log$n_tune_resamps),
    fmt("Fit resamples:           %d", run_log$n_fit_resamps),
    "",
    "------------------------------------------------------------",
    "MODEL SPECIFICATION",
    "------------------------------------------------------------",
    fmt("Model type:      %s", run_log$model$type),
    fmt("Engine:          %s", run_log$model$engine),
    fmt("Outcome:         %s", run_log$outcome),
    fmt("Predictors:      %d", run_log$n_predictors),
    fmt("Trees (tuning):  %d", run_log$model$ntree_tune),
    fmt("Trees (final):   %d", run_log$model$ntree_final),
    ""
  )
  
  # Parallelization info
  if (!is.null(run_log$parallelization)) {
    lines <- c(
      lines,
      "Parallelization:",
      fmt("  Outer threads (tuning/CV): %d", run_log$parallelization$outer_threads),
      fmt("  Inner threads (ranger):    %d", run_log$parallelization$inner_threads),
      fmt("  Total cores (tune/CV):     %d", 
          run_log$parallelization$outer_threads * run_log$parallelization$inner_threads),
      fmt("  Final fit threads:         %d", run_log$parallelization$final_threads),
      ""
    )
  } else if (!is.null(run_log$num_threads)) {
    # Fallback for old format
    lines <- c(lines, fmt("Threads used:    %d", run_log$num_threads), "")
  }
  
  # Engine args
  lines <- c(lines, "Engine arguments (final model):")
  for (nm in names(run_log$model$engine_args)) {
    val <- run_log$model$engine_args[[nm]]
    # Handle different types better
    if (is.logical(val)) {
      lines <- c(lines, fmt("  - %-20s %s", paste0(nm, ":"), toupper(as.character(val))))
    } else if (is.numeric(val)) {
      lines <- c(lines, fmt("  - %-20s %g", paste0(nm, ":"), val))
    } else {
      lines <- c(lines, fmt("  - %-20s %s", paste0(nm, ":"), val))
    }
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "TUNING",
    "------------------------------------------------------------",
    fmt("Method:           %s", run_log$tuning$method),
    fmt("Grid size:        %d", run_log$tuning$grid_size),
    fmt("Min_n range:      2 to %d", run_log$tuning$min_n_upper),
    fmt("Mtry range:       1 to %d", run_log$n_predictors),
    fmt("Selection rule:   %s", run_log$tuning$selection_rule),
    fmt("Metric optimized: %s", run_log$tuning$tune_metric),
    "",
    "Best hyperparameters:"
  )
  
  # Best params - handle numeric formatting
  for (nm in names(run_log$tuning$best_params)) {
    val <- run_log$tuning$best_params[[nm]]
    if (is.numeric(val)) {
      lines <- c(lines, fmt("  - %-20s %g", paste0(nm, ":"), val))
    } else {
      lines <- c(lines, fmt("  - %-20s %s", paste0(nm, ":"), val))
    }
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "PERFORMANCE",
    "------------------------------------------------------------"
  )
  
  # CV metrics
  if (!is.null(run_log$cv_metrics)) {
    lines <- c(lines, "Cross-validated performance:")
    cvm <- run_log$cv_metrics
    for (i in seq_len(nrow(cvm))) {
      # Add standard error if available
      if ("std_err" %in% names(cvm)) {
        lines <- c(
          lines,
          fmt(
            "  %-8s mean = %.4f, SE = %.4f",
            paste0(cvm$.metric[i], ":"),
            cvm$mean[i],
            cvm$std_err[i]
          )
        )
      } else {
        lines <- c(
          lines,
          fmt(
            "  %-8s %.4f",
            paste0(cvm$.metric[i], ":"),
            cvm$mean[i]
          )
        )
      }
    }
  }
  
  lines <- c(lines, "")
  
  # Test metrics
  if (!is.null(run_log$test_metrics)) {
    lines <- c(lines, "Test-set performance:")
    tm <- run_log$test_metrics
    for (i in seq_len(nrow(tm))) {
      lines <- c(
        lines,
        fmt(
          "  %-8s %.4f",
          paste0(tm$.metric[i], ":"),
          tm$.estimate[i]
        )
      )
    }
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "PREDICTOR VARIABLES",
    "------------------------------------------------------------"
  )
  
  # Add predictor list
  if (!is.null(run_log$predictor_names)) {
    lines <- c(lines, fmt("Total: %d predictors", length(run_log$predictor_names)))
    lines <- c(lines, "")
    
    # Print in columns for readability
    pred_list <- run_log$predictor_names
    n_cols <- 3
    n_rows <- ceiling(length(pred_list) / n_cols)
    
    for (i in seq_len(n_rows)) {
      row_items <- character(n_cols)
      for (j in seq_len(n_cols)) {
        idx <- (j - 1) * n_rows + i
        if (idx <= length(pred_list)) {
          row_items[j] <- sprintf("%-25s", pred_list[idx])
        } else {
          row_items[j] <- ""
        }
      }
      lines <- c(lines, paste0("  ", paste(row_items, collapse = "")))
    }
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "SEEDS (for reproducibility)",
    "------------------------------------------------------------"
  )
  
  for (nm in names(run_log$seeds)) {
    lines <- c(lines, fmt("  %-25s %d", paste0(nm, ":"), run_log$seeds[[nm]]))
  }
  
  lines <- c(
    lines,
    "",
    "============================================================",
    "END OF RUN LOG",
    "============================================================"
  )
  
  writeLines(lines, con = file)
  invisible(file)
}