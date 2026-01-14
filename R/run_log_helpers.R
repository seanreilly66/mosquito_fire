start_console_log <- function(path) {

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  
  sink(con)
  sink(con, type = "message")
  
  # return stopper
  function() {
    sink(type = "message")
    sink()
    close(con)
    
  }
}


write_run_log_txt <- function(run_log, file) {
  
  fmt <- function(...) sprintf(...)
  
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
    fmt("End time:        %s", run_log$end_time),
    fmt("Runtime (sec):   %.1f", run_log$runtime_s),
    "",
    "------------------------------------------------------------",
    "INPUT DATA",
    "------------------------------------------------------------",
    fmt("Training raster: %s", run_log$input_file),
    fmt("ROI GPKG:        %s", run_log$roi_gpkg),
    fmt("ROI layer:       %s", run_log$roi_layer),
    "",
    fmt("Total samples:   %d", run_log$n_total),
    fmt("Train samples:   %d", run_log$n_train),
    fmt("Test samples:    %d", run_log$n_test),
    "",
    "------------------------------------------------------------",
    "CROSS-VALIDATION DESIGN",
    "------------------------------------------------------------",
    fmt("Initial test proportion: %.2f", run_log$test_prop),
    fmt("Spatial buffer (m):      %d", run_log$buffer_dist_m),
    fmt("V-folds:                %d", run_log$v_folds),
    fmt("Tuning repeats:         %d", run_log$tune_reps),
    fmt("Fit repeats:            %d", run_log$fit_reps),
    fmt("Tuning resamples:       %d", run_log$n_tune_resamps),
    fmt("Fit resamples:          %d", run_log$n_fit_resamps),
    "",
    "------------------------------------------------------------",
    "MODEL SPECIFICATION",
    "------------------------------------------------------------",
    fmt("Model type:      %s", run_log$model$type),
    fmt("Engine:          %s", run_log$model$engine),
    fmt("Predictors:      %d", run_log$n_predictors),
    fmt("Trees (tuning):  %d", run_log$model$ntree_tune),
    fmt("Trees (final):   %d", run_log$model$ntree_final),
    fmt("Threads used:    %d", run_log$num_threads),
    "",
    "Engine arguments:"
  )
  
  # Engine args
  for (nm in names(run_log$model$engine_args)) {
    lines <- c(lines, fmt("  - %s: %s", nm, run_log$model$engine_args[[nm]]))
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "TUNING",
    "------------------------------------------------------------",
    fmt("Method:          %s", run_log$tuning$method),
    fmt("Grid size:       %d", run_log$tuning$grid_size),
    fmt("Selection rule:  %s", run_log$tuning$selection_rule),
    fmt("Metric optimized:%s", run_log$tuning$tune_metric),
    "",
    "Best hyperparameters:"
  )
  
  # Best params
  for (nm in names(run_log$tuning$best_params)) {
    lines <- c(lines, fmt("  - %s: %s", nm, run_log$tuning$best_params[[nm]]))
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
      lines <- c(
        lines,
        fmt(
          "  %s: mean = %.4f",
          cvm$.metric[i],
          cvm$mean[i]
        )
      )
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
          "  %s: %.4f",
          tm$.metric[i],
          tm$.estimate[i]
        )
      )
    }
  }
  
  lines <- c(
    lines,
    "",
    "------------------------------------------------------------",
    "SEEDS",
    "------------------------------------------------------------"
  )
  
  for (nm in names(run_log$seeds)) {
    lines <- c(lines, fmt("  %s: %d", nm, run_log$seeds[[nm]]))
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
