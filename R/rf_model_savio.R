home_root = file.path('/', 'global', 'home', 'users', 'seanreilly')
scratch_root = file.path('/', 'global', 'scratch', 'users', 'seanreilly', 'mosquito_fire')

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(tidymodels)
  library(spatialsample)
  library(sf)
  library(tictoc)
  # library(finetune)
  library(mirai)
})
source(file.path(home_root, 'R', 'run_log_helpers.R'))
source(file.path(home_root, 'R', 'plot_helpers.R'))
source(file.path(home_root, 'R', 'spatial_cv.R'))

# ==============================================================================
# Inputs
# ==============================================================================


burn_file <- file.path(scratch_root, "data", "burn_severity", "mosquito_cbi.tif") 
metric_dir <- file.path(scratch_root, "data", "fuelrasters") 

roi_gpkg = file.path(scratch_root, 'data', 'mosquito_fire.gpkg')
roi_layer = 'mosquito_study_area'

save_console = TRUE

# spatial block parameters
initial_prop = 0.2
v_folds = 5
blocks_per_fold = 3
tune_reps = 1
fit_reps = 10
buffer_dist = 2000
n_fold_plots = 10

# rf parameters
# ntrees defined in model setup below
grid_size = 125
min_n_upper = 200
mtry_frac = 0.5
view_tuning = TRUE
tune_metric = 'ccc'
recipe_info = TRUE

tune_trees = 200
cv_trees = 1000
final_trees = 5000


# session information
run_type <- 'rf_cbi_pt_fuels'
start_time <- Sys.time()
run_id <- format(start_time, "%Y%m%d_%H%M%S")

test_run = FALSE
sample_n = 1000


out_dir <- file.path(scratch_root, "outputs", "rf", paste0(run_type, "_", run_id))
dir.create(out_dir, recursive = TRUE)
prefix <- file.path(out_dir, paste0(run_type, '_', run_id))

# ==============================================================================
# Prep
# ==============================================================================

total_cores <- as.integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", 
    unset = Sys.getenv("SLURM_CPUS_ON_NODE", 
      unset = parallel::detectCores())))
outer_workers <- v_folds
ranger_threads <- max(1L, total_cores %/% outer_workers) 

# ------------------------------------------------------------------------------
# Initialize console log
# ------------------------------------------------------------------------------

if (save_console) {
  console_log <- start_console_log(paste0(prefix, '_console.txt'), split = FALSE)
}

tic('Full run')

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

message('Run: ', prefix)
message('1. Read data: ', Sys.time())

roi <- vect(roi_gpkg, roi_layer)

burn_rast <- rast(burn_file) %>%
  crop(roi) %>%
  mask(roi)
names(burn_rast) <- 'rdnbr'

pred_rast <- list.files(metric_dir, pattern = '\\.tif$', full.names = TRUE) %>%
  rast() %>%
  crop(roi) %>%
  mask(roi)

input_df <- c(burn_rast, pred_rast) %>%
  crop(roi) %>%
  mask(roi) %>%
  as.points() %>%
  st_as_sf() %>%
  mutate(st_coordinates(geometry) |> as_tibble())

# if (test_run) input_df <- input_df %>% sample_n(1000)

rm(burn_rast, pred_rast)
gc()

# ------------------------------------------------------------------------------
# Initial split
# ------------------------------------------------------------------------------

message('2. Spatial splits: ', Sys.time())

init_split <- spatial_init_split(
  roi,
  pts = input_df,
  test_prop = initial_prop,
  buffer_dist = buffer_dist,
  seed = 7298
)

train = training(init_split)
test = testing(init_split)

n_total <- nrow(input_df)
n_train <- nrow(train)
n_test <- nrow(test)

cat(
  sprintf(
    'Prop test: %.2f\n',
    n_test / (n_test + n_train)
  )
)

# --- Initial split plot ---
init_split_plot <- autoplot(init_split) + 
  ggtitle('Initial split') 

ggsave(plot = init_split_plot, 
       filename = paste0(prefix, '_initial_split.png'),
       width = 8, 
       height = 6, 
       dpi = 300)

# ------------------------------------------------------------------------------
# Spatial blocks
# ------------------------------------------------------------------------------

fold_plan <- spatial_block_folds(
  pts = train,
  v = v_folds,
  blocks_per_fold = blocks_per_fold,
  repeats = fit_reps,
  buffer_dist = buffer_dist,
  seed = 838
)

tune_folds <- spatial_block_rset(
  fold_plan,
  pts = train |> st_drop_geometry(),
  v = v_folds,
  subsample_repeats = tune_reps,
  subsample_seed = 4095
)

fit_folds <- spatial_block_rset(
  fold_plan,
  pts = train |> st_drop_geometry(),
  v = v_folds
)

n_tune_resamps <- nrow(tune_folds)
n_fit_resamps <- nrow(fit_folds)

# ==============================================================================
# Model recipe
# ==============================================================================

message('3. Model Specification: ', Sys.time())

rf_recipe <- recipe(rdnbr ~ ., train |> st_drop_geometry()) %>% 
  update_role(X, Y, new_role = 'id')

if (recipe_info) {
  
  prep(rf_recipe)
  rf_recipe$var_info
  
}

# ==============================================================================
# Model setup
# ==============================================================================

# Memory parallel:
num_threads_outer = outer_workers  # N data copies
num_threads_inner = ranger_threads  # Compute per copy
num_threads_final = total_cores


rf_model <- rand_forest(
  mtry = tune(),
  min_n = tune(),
  trees = !!tune_trees) %>%
  set_mode('regression') %>%
  set_engine('ranger', 
             num.threads = !!num_threads_inner,
             verbose = FALSE,
             splitrule = 'variance',
             replace = TRUE,
             sample.fraction = 0.95,
             min.bucket = 3,
             max.depth = tune())

rf_workflow <- workflow() %>%
  add_model(rf_model) %>%
  add_recipe(rf_recipe)

# Memory cleanup - large objects not needed anymore
if (exists('input_df')) rm(input_df)
if (exists('roi')) rm(roi)
gc()

# ==============================================================================
# Hyperparameter tuning
# ==============================================================================

message('4. Hyperparameter tuning: ', Sys.time())
tune_start_time <- Sys.time()

n_predictors <- sum(rf_recipe$var_info$role == 'predictor')

rf_grid <- expand.grid(
  mtry = seq(30, n_predictors, 5),
  min_n = c(5, 15, 25),
  max.depth = c(25, 35, 45)
)

tic('tune')
daemons(num_threads_outer)

set.seed(7129)
rf_tune <- rf_workflow %>%
  tune_grid(
    resamples = tune_folds,
    grid      = rf_grid,
    metrics   = yardstick::metric_set(rmse, rsq, ccc),
    control   = control_grid(
      parallel_over = "resamples",
      save_pred = FALSE,
      verbose = TRUE
    )
  )

toc()
daemons(0)


rf_tuned_param <- rf_tune %>%
  select_by_one_std_err(
    metric = "ccc", 
    max.depth,
    desc(min_n),
    mtry
  )

rf_tuned_model <- rf_model %>%
  set_args(trees = !!cv_trees)

rf_tuned_workflow <- workflow() %>%
  add_model(rf_tuned_model) %>%
  add_recipe(rf_recipe) %>%
  finalize_workflow(rf_tuned_param) 

# Tuning result visualization

if (view_tuning) {
  
  print('Tuning results')
  
  plot(rf_grid)

  rf_tune %>%
    show_best(metric = 'ccc')

}

# Write tuning outputs

grid_plot <- ggplot(
  data = rf_grid,
  mapping = aes(x = mtry, y = min_n, color = max.depth |> as.factor())) + 
  geom_jitter(width = 0.5, height = 0.5) +
  scale_color_discrete()

ggsave(plot = grid_plot, filename = paste0(prefix, '_grid.png'))

tune_plot <- rf_tune %>%
  collect_metrics() %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean , color = min_n)) +
  # geom_line(linewidth = 1.5, alpha = 0.6) +
  geom_point(size = 2) +
  facet_grid(.metric ~ max.depth, scales = 'free', axes = 'all') +
  scale_color_viridis_d(option = "plasma", begin = .9, end = 0)

ggsave(plot = tune_plot, filename = paste0(prefix, '_tune.png'))

saveRDS(rf_tune, paste0(prefix, '_tune.rds'))

write_csv(rf_tuned_param, paste0(prefix, '_best_params.csv'))

rf_tune %>% 
  collect_metrics() %>%
  write_csv(paste0(prefix, '_tune_metrics.csv'))

tune_end_time <- Sys.time()

# ==============================================================================
# CV model performance
# ==============================================================================

message('5. Fit CV model: ', Sys.time())
cv_start_time <- Sys.time()

tic('cv')
daemons(num_threads_outer)

set.seed(4285)

rf_cv <- rf_tuned_workflow %>%
  fit_resamples(
    resamples = fit_folds,
    metrics = yardstick::metric_set(ccc, rmse, rsq),
    control = control_resamples(
      save_pred = FALSE,
      pkgs = c('ranger', 'tidymodels', 'spatialsample'),
      save_workflow = FALSE
    ))

toc()
daemons(0)

saveRDS(rf_cv, paste0(prefix, '_tuned_cv_fit.rds'))

cv_results <- rf_cv %>%
  collect_metrics()

cv_results

write_csv(cv_results, paste0(prefix, '_cv_metrics.csv'))

cv_end_time <- Sys.time()

# ==============================================================================
# Final model fit
# ==============================================================================

message('6. Fit final model: ', Sys.time())
final_start_time <- Sys.time()

tic('final fit')

rf_final_model <- rand_forest(
  mtry = tune(),
  min_n = tune(),
  trees = !!final_trees) %>%
  set_mode('regression') %>%
  set_engine(
    'ranger',
    num.threads = !!num_threads_final,
    importance = 'permutation',
    write.forest = TRUE,
    keep.inbag = TRUE,
    node.stats = TRUE,
    seed = 8465,
    splitrule = 'variance',
    replace = TRUE,
    sample.fraction = 0.95,
    min.bucket = 3,
    max.depth = tune()
  )

rf_final_workflow <- rf_workflow %>%
  update_model(rf_final_model) %>%
  finalize_workflow(rf_tuned_param)

rf_fit <- rf_final_workflow %>%
  last_fit(
    split = init_split,
    metrics = yardstick::metric_set(rmse, rsq, ccc),
    control = control_last_fit(
      verbose = TRUE
    ))

toc() # Final fit

saveRDS(rf_fit, paste0(prefix, '_final_fit.rds'))
saveRDS(rf_final_workflow, paste0(prefix, '_final_workflow.rds'))

fit_results <- rf_fit %>%
  collect_metrics()

write_csv(fit_results, paste0(prefix, '_final_fit_metrics.csv'))

final_end_time <- Sys.time()
run_time = toc() # Full run

# ==============================================================================
# Performance evaluation
# ==============================================================================

final_preds <- rf_fit %>%
  collect_predictions()

fit_plot = final_preds %>%
  ggplot(aes(rdnbr, .pred)) +
  geom_abline(lty = 2,
              col = "gray",
              linewidth = 1.5) +
  geom_point(alpha = 0.5, color = 'firebrick') +
  labs(x = 'Observed', y = 'Predicted') +
  coord_obs_pred()

ggsave(plot = fit_plot, filename = paste0(prefix, '_predobs.png'))

# ==============================================================================
# Export logs
# ==============================================================================

if (save_console) {
  
  console_log()
  
}

message('Best param: ')
print(rf_tuned_param)

message('CV performance: ')
print(cv_results)

message('Final fit: ')
print(fit_results)

run_log <- list(
  
  # Metadata
  run_type = run_type,
  run_id = run_id,
  out_dir = out_dir,
  prefix = prefix,
  
  start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
  end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  
  # Timing breakdown (new)
  timing = list(
    data_prep_s = as.numeric(difftime(tune_start_time, start_time, units = 'secs')),
    tuning_s = as.numeric(difftime(tune_end_time, tune_start_time, units = 'secs')),
    cv_s = as.numeric(difftime(cv_end_time, cv_start_time, units = 'secs')),
    final_fit_s = as.numeric(difftime(final_end_time, final_start_time, units = 'secs')),
    total_s = as.numeric(difftime(Sys.time(), start_time, units = 'secs'))
  ),
  
  # Seeds
  seeds = list(
    init_split = 7298,
    fold_plan = 838, 
    tune_subsample = 4095,
    grid = 345,
    tuning = 7129,
    cv_fit = 4285,
    final_fit = 8465
  ),
  
  # Input data
  input_file = burn_file,
  metric_dir = metric_dir,
  roi_gpkg = roi_gpkg,
  roi_layer = roi_layer,
  
  # Data sizes
  n_total = n_total,
  n_train = n_train,
  n_test = n_test,
  
  # CV design
  test_prop = initial_prop,
  v_folds = v_folds,
  tune_reps = tune_reps,
  fit_reps = fit_reps,
  buffer_dist_m = buffer_dist,
  n_tune_resamps = n_tune_resamps, 
  n_fit_resamps  = n_fit_resamps,
  
  # Recipe / predictors
  outcome = "rdnbr",
  id_fields = c("geometry", "X", "Y"),
  n_predictors = n_predictors,
  predictor_names = rf_recipe$var_info %>%
    dplyr::filter(role == "predictor") %>%
    dplyr::pull(variable),
  
  # Model / tuning
  model = list(
    type = "Random Forest",
    engine = "ranger",
    ntree_tune = tune_trees,
    ntree_cv = cv_trees,
    ntree_final = final_trees,
    engine_args = list(
      num.threads = num_threads_final,
      importance = "permutation",
      write.forest = TRUE,
      keep.inbag = TRUE,
      node.stats = TRUE,
      seed = 8465
    )
  ),
  
  tuning = list(
    method = "tune_grid",
    grid_size = nrow(rf_grid),
    min_n_upper = min_n_upper,
    tune_metric = tune_metric,
    selection_rule = "best",
    best_params = as.list(rf_tuned_param[1, ])
  ),
  
  tuning_grid = rf_grid,
  
  # Performance
  cv_metrics = cv_results,
  test_metrics = fit_results,
  
  ## Parallelization settings
  parallelization = list(
    outer_threads = num_threads_outer,
    inner_threads = num_threads_inner,
    final_threads = num_threads_final
  )
)

saveRDS(run_log, paste0(prefix, "_run_log.rds"))

write_run_log_txt(
  run_log,
  file = paste0(prefix, "_run_log.txt")
)


# ==============================================================================
