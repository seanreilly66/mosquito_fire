# =============================================================================
# Exploratory continuous modeling: fuels -> rDNBR (random split; upper bound)
# Inputs:
#   - fuel rasters as a zip of .tif files
#   - burn severity as a standalone .tif
#   - study area polygon in a geopackage layer "mosquito_study_area"
# Output:
#   - model comparison table (null vs ridge vs RF)
#   - RF pred vs obs plot
#   - RF permutation importance table + VIP plot
#   - optional partial dependence for top predictors
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(tidyverse)
  library(tidymodels)
  library(ranger)
  library(vip)
})

# Sampling size for exploratory modeling (adjust to your RAM/time)
n_samp <- 150000

# Reproducibility
set.seed(1)


burn_file  <- 'data/burn_severity/mosquito_rdnbr_mtbs.tif'     # RdNBR response
metric_dir <- "data/fuelrasters"            # predictor rasters

roi_gpkg = 'data/mosquito_fire.gpkg'
roi_layer = 'mosquito_study_area'

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

roi <- vect(roi_gpkg, roi_layer)

burn_rast <- rast(burn_file) %>%
  crop(roi) %>%
  mask(roi)
names(burn_rast) <- 'rdnbr'

pred_rast <- rast(list.files(metric_dir, pattern="\\.tif$", full.names=TRUE)) %>%
  crop(roi) %>%
  mask(roi)

r_all <- c(burn_rast, pred_rast) %>%
  crop(roi) %>%
  mask(roi)

# ---------------------------- Sample pixels ----------------------------------

df <- spatSample(
  r_all,
  size = n_samp,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE,
  xy = TRUE
)

# Pred stats
cat("\nSampled rows:", nrow(df), "\n")
cat("Predictors:", ncol(df) - 3, "(excluding rdnbr, x, y)\n\n")
print(summary(df$rdnbr))

# Optional: remove extreme outliers if desired (comment out if not needed)
# df <- df |> filter(between(rdnbr, quantile(rdnbr, 0.001), quantile(rdnbr, 0.999)))

# -------------------------- Train/test split ---------------------------------
set.seed(2)
split <- initial_split(df, prop = 0.8)
train <- training(split)
test  <- testing(split)

# ------------------------------- Recipe --------------------------------------

rec <- recipe(rdnbr ~ ., data = train) |>
  update_role(x, y, new_role = "id") |>
  step_rm(x, y) |>
  step_zv(all_predictors())

# --------------------------- Model definitions -------------------------------

base_wf <- workflow() |>
  add_recipe(rec)

# Null model (predict mean)
null_mod <- null_model() |>
  set_mode("regression")

null_wf <- base_wf |> 
  add_model(null_mod)

null_fit <- fit(null_wf, train)

# Ridge regression baseline (handles multicollinearity; linear signal)

ridge_mod <- linear_reg(penalty = 0.01, mixture = 0) |> 
  set_engine("glmnet")

ridge_wf <- base_wf |>
  add_model(ridge_mod)

ridge_fit <- fit(ridge_wf, train)

# Random forest 
p <- ncol(train) - 3
rf_mod <- rand_forest(
  mode  = "regression",
  mtry  = max(1, floor(sqrt(p))),
  trees = 800,
  min_n = 5
) |>
  set_engine("ranger",
             importance = "permutation",
             num.threads = parallel::detectCores() - 1)

rf_wf <- base_wf %>% add_model(rf_mod)

rf_fit <- fit(rf_wf, train)

# ------------------------------ Evaluation -----------------------------------
reg_metrics <- metric_set(rmse, mae, rsq)

eval_model <- function(wf_fit, name){
  pred <- predict(wf_fit, test) |> bind_cols(test |> select(rdnbr))
  m <- reg_metrics(pred, truth = rdnbr, estimate = .pred) |> mutate(model = name)
  list(metrics = m, pred = pred)
}

res_null  <- eval_model(null_fit, "null")
res_ridge <- eval_model(ridge_fit, "ridge")
res_rf    <- eval_model(rf_fit, "rf")

metrics_tbl <- bind_rows(res_null$metrics, res_ridge$metrics, res_rf$metrics) |>
  select(model, .metric, .estimate) |>
  pivot_wider(names_from = .metric, values_from = .estimate) |>
  arrange(desc(rsq))

cat("\n================ Model comparison (random split) ================\n")
print(metrics_tbl)
cat("=================================================================\n\n")

# ----------------------- Predicted vs Observed plot --------------------------
ggplot(res_rf$pred, aes(x = rdnbr, y = .pred)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    title = "Random Forest: predicted vs observed (random split)",
    x = "Observed rDNBR",
    y = "Predicted rDNBR"
  ) +
  theme_minimal()

# ---------------------- RF variable importance -------------------------------
rf_fit_obj <- extract_fit_parsnip(rf_fit)$fit

imp <- tibble(
  feature = names(rf_fit_obj$variable.importance),
  importance = as.numeric(rf_fit_obj$variable.importance)
) |>
  arrange(desc(importance))

cat("\n================ RF permutation importance (top 25) =============\n")
print(imp |> slice_head(n = 25))
cat("=================================================================\n\n")

vip::vip(rf_fit_obj, num_features = 20)

# ---------------------- Optional: partial dependence -------------------------
# Partial dependence is easiest if we refit ranger directly on a data.frame.
# This section is optional; uncomment to run.

# suppressPackageStartupMessages({
#   library(pdp)
# })
#
# rf2 <- ranger(
#   rdnbr ~ .,
#   data = train |> select(-x, -y),
#   num.trees = 800,
#   min.node.size = 5,
#   mtry = max(1, floor(sqrt(ncol(train |> select(-x, -y)) - 1))),
#   importance = "permutation"
# )
#
# top_feats <- imp$feature[1:3]
# for (f in top_feats) {
#   pd <- partial(
#     object = rf2,
#     pred.var = f,
#     train = train |> select(-x, -y),
#     grid.resolution = 30
#   )
#   plot(pd, main = paste("PDP:", f))
# }

# ----------------------------- Quick readout ---------------------------------
# If rf rsq >> ridge rsq >> null, strong nonlinear relationship exists.
# If ridge improves over null but rf does not, relationship is largely linear.
# If rf ~ ridge ~ null, predictors have weak/no relationship at this scale.

invisible(metrics_tbl)
