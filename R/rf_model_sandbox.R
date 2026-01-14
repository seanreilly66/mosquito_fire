suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(tidymodels)
  library(spatialsample)
  library(sf)
  library(tictoc)
  library(finetune)
  library(mirai)
})
source('R/run_log_helpers.R')
source('R/plot_helpers.R')

# ==============================================================================
# Inputs
# ==============================================================================

input_file = 'data/temp/synth_train.tif'
# test_file = 'data/temp/synth_test.tif'

roi_gpkg = 'data/mosquito_fire.gpkg'
roi_layer = 'mosquito_study_area'

save_console = TRUE

# spatial block parameters
initial_prop = 0.2
v_folds = 10
tune_reps = 30
fit_reps = 100
buffer_dist = 2000
plot_folds = TRUE
plot_location = 'data/temp/folds'

# rf parameters
# ntrees defined in model setup below
grid_size = 150
min_n_upper = 80
view_tuning = TRUE
tune_metric = 'rmse'
recipe_info = TRUE

# session information
run_type <- 'rf_sandbox'
start_time <- Sys.time()
run_id <- format(start_time, "%Y%m%d_%H%M%S")
out_dir <- glue::glue('outputs/rf/{run_type}_{run_id}')
dir.create(out_dir, recursive = TRUE)
prefix <- file.path(out_dir, paste0(run_type, '_', run_id))

# ==============================================================================
# Prep
# ==============================================================================

# ------------------------------------------------------------------------------
# Initialize console log
# ------------------------------------------------------------------------------

if (save_console) {
  console_log <- start_console_log(paste0(prefix, '_console.txt'))
}

tic('Full run')

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

message('Run: ', prefix)
message('1. Read data: ', Sys.time())

roi <- vect(roi_gpkg, roi_layer)

input_df <- rast(input_file) %>%
  crop(roi) %>%
  mask(roi) %>%
  as.points() %>%
  st_as_sf() %>%
  mutate(st_coordinates(geometry) |> as_tibble())

# sandbox_test <- rast(test_file) %>%
#   crop(roi) %>%
#   mask(roi) %>%
#   as.points() %>%
#   st_as_sf() %>%
#   mutate(st_coordinates(geometry) |> as_tibble())

# ------------------------------------------------------------------------------
# Initial split
# ------------------------------------------------------------------------------

message('2. Spatial splits: ', Sys.time())

set.seed(7298)
init_split <- tidysdm::spatial_initial_split(
  data = input_df,
  prop = initial_prop,
  strategy = spatial_clustering_cv,
  buffer = buffer_dist,
  # n = 5
)

message('Initial split')
init_split

train = training(init_split)
test = testing(init_split)

init_split_plot <- autoplot(init_split) + 
  ggtitle('Initial split') 

ggsave(plot = init_split_plot, filename = paste0(prefix, '_initial_split.png'))


# ------------------------------------------------------------------------------
# Spatial blocks
# ------------------------------------------------------------------------------

set.seed(838)
tune_blocks <- spatial_block_cv(
  data = train,
  method = 'random',
  v = v_folds,
  repeats = tune_reps,
  buffer = buffer_dist
)

set.seed(4095)
fit_blocks <- spatial_block_cv(
  data = train,
  method = 'random',
  v = v_folds,
  repeats = fit_reps,
  buffer = buffer_dist
)

# Plot folds to file
if (plot_folds) {
  
  purrr::walk(fit_blocks$splits, 
              function(x) {
                p = autoplot(x)
                
                suppressMessages(
                  ggsave(
                    plot = p,
                    filename = tempfile(
                      fileext = '.png',
                      tmpdir = file.path(out_dir, 'folds')),
                    create.dir = TRUE))
              }
  )
}

# ==============================================================================
# Model recipe
# ==============================================================================

message('3. Model Specification: ', Sys.time())

rf_recipe <- recipe(rdnbr ~ ., train) %>% 
  update_role(rdnbr_true, new_role = 'reference') %>%
  update_role(geometry, X, Y, new_role = 'id')

if (recipe_info) {
  
  prep(rf_recipe)
  rf_recipe$var_info
  
}

# ==============================================================================
# Model setup
# ==============================================================================

rf_model <- rand_forest(
  mtry = tune(),
  min_n = tune(),
  trees = 1000) %>%
  set_mode('regression') %>%
  set_engine('ranger')

rf_workflow <- workflow() %>%
  add_model(rf_model) %>%
  add_recipe(rf_recipe)

num_threads = parallelly::availableCores() - 1

# ==============================================================================
# Hyperparameter tuning
# ==============================================================================

message('4. Hyperparameter tuning: ', Sys.time())

n_predictors <- sum(rf_recipe$var_info$role == 'predictor')

rf_params <- parameters(
  mtry(range = c(1, n_predictors)),
  min_n(range = c(2, min_n_upper))
)

set.seed(345)
rf_grid <- grid_space_filling(rf_params, size = grid_size)

tic('tune')
daemons(num_threads)

set.seed(7129)
rf_tune <- rf_workflow %>%
  tune_race_anova(
    resamples = tune_blocks,
    grid = rf_grid,
    metrics = yardstick::metric_set(rmse, rsq),
    control = finetune::control_race(
      verbose_elim = TRUE,
      pkgs = c('sf', 'ranger', 'finetune', 'tidymodels', 'spatialsample'),
      save_pred = TRUE,
      save_workflow = TRUE
    )
  )

toc()
daemons(0)

rf_tuned_param <- rf_tune %>%
  select_by_one_std_err(metric = tune_metric, desc(min_n), mtry) %>% 
  dplyr::slice(1)

rf_tuned_workflow <- rf_workflow %>%
  finalize_workflow(rf_tuned_param) 

# Tuning result visualization

if (view_tuning) {
  
  print('Tuning results')
  
  plot(rf_grid)

  rf_tune %>%
    show_best(metric = 'rmse')

}

# Write tuning outputs

grid_plot <- ggplot(rf_grid, mapping = aes(x = mtry, y = min_n)) + geom_point()

ggsave(plot = grid_plot, filename = paste0(prefix, '_grid.png'))

tune_plot <- rf_tune %>%
  collect_metrics() %>%
  mutate(mtry = factor(mtry)) %>%
  ggplot(aes(min_n, mean, color = mtry)) +
  geom_line(linewidth = 1.5, alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ .metric, scales = "free", ncol = 2) +
  scale_color_viridis_d(option = "plasma", begin = .9, end = 0)

ggsave(plot = tune_plot, filename = paste0(prefix, '_tune.png'))

saveRDS(rf_tune, paste0(prefix, '_tune.rds'))

write_csv(rf_tuned_param, paste0(prefix, '_best_params.csv'))

rf_tune %>% 
  collect_metrics() %>%
  write_csv(paste0(prefix, '_tune_metrics.csv'))

# ==============================================================================
# CV model performance
# ==============================================================================

message('5. Fit CV model: ', Sys.time())

tic('cv')
daemons(num_threads)

set.seed(4285)

rf_cv <- rf_tuned_workflow %>%
  fit_resamples(
    resamples = fit_blocks,
    metrics = yardstick::metric_set(rmse, rsq),
    control = control_resamples(
      save_pred = TRUE,
      pkgs = c('sf', 'ranger', 'finetune', 'tidymodels', 'spatialsample'),
      save_workflow = TRUE,
      parallel_over = 'everything'
    ))

toc()
daemons(0)

saveRDS(rf_cv, paste0(prefix, '_tuned_cv_fit.rds'))

cv_results <- rf_cv %>%
  collect_metrics()

cv_results

write_csv(cv_results, paste0(prefix, '_cv_metrics.csv'))

# ==============================================================================
# Final model fit
# ==============================================================================

message('6. Fit final model: ', Sys.time())

tic('final fit')

rf_final_model <- rand_forest(
  mtry = tune(),
  min_n = tune(),
  trees = 5000) %>%
  set_mode('regression') %>%
  set_engine(
    'ranger',
    num.threads = num_threads,
    importance = 'permutation',
    write.forest = TRUE,
    keep.inbag = TRUE,
    node.stats = TRUE,
    seed = 8465
  )

rf_final_workflow <- rf_workflow %>%
  update_model(rf_final_model) %>%
  finalize_workflow(rf_tuned_param)

rf_fit <- rf_final_workflow %>%
  last_fit(
    split = init_split,
    metrics = yardstick::metric_set(rmse, rsq),
    control = control_last_fit(
      verbose = TRUE
    ))

toc() # Final fit

saveRDS(rf_fit, paste0(prefix, '_final_fit.rds'))
saveRDS(rf_final_workflow, paste0(prefix, '_final_workflow.rds'))

fit_results <- rf_fit %>%
  collect_metrics()

write_csv(fit_results, paste0(prefix, '_final_fit_metrics.csv'))

run_time = toc() # Full run

# ==============================================================================
# Performance evaluation
# ==============================================================================

fit_results <- rf_fit %>%
  collect_metrics()

write_csv(fit_results, paste0(prefix, '_final_fit_metrics.csv'))

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
  runtime_s = as.numeric(difftime(Sys.time(), start_time, units = 'secs')),
  
  # Seeds
  seeds = list(
    init_split = 7298,
    tune_blocks = 838,
    fit_blocks = 4095,
    grid = 345,
    tuning = 7129,
    cv_fit = 4285,
    final_fit = 8465
  ),
  
  # Input data
  input_file = input_file,
  roi_gpkg = roi_gpkg,
  roi_layer = roi_layer,
  
  # Data sizes
  n_total = nrow(input_df),
  n_train = nrow(train),
  n_test = nrow(test),
  
  # CV design
  test_prop = initial_prop,
  v_folds = v_folds,
  tune_reps = tune_reps,
  fit_reps = fit_reps,
  buffer_dist_m = buffer_dist,
  n_tune_resamps = nrow(tune_blocks),
  n_fit_resamps  = nrow(fit_blocks),
  
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
    ntree_tune = 1000,
    ntree_final = 5000,
    engine_args = list(
      num.threads = num_threads,
      importance = "permutation",
      write.forest = TRUE,
      keep.inbag = TRUE,
      node.stats = TRUE,
      seed = 8465
    )
  ),
  
  tuning = list(
    method = "tune_race_anova",
    grid_size = grid_size,
    min_n_upper = min_n_upper,
    tune_metric = tune_metric,
    selection_rule = "one-standard-error",
    best_params = as.list(rf_tuned_param[1, ])
  ),
  
  # Performance
  cv_metrics = cv_results,
  test_metrics = fit_results,
  
  # Environment
  num_threads = num_threads
)

saveRDS(run_log, paste0(prefix, "_run_log.rds"))

write_run_log_txt(
  run_log,
  file = paste0(prefix, "_run_log.txt")
)

# ==============================================================================