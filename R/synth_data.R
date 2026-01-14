# Standalone synthetic raster dataset generator (variable-selection stress test)
# - Even grid raster (~1000 cells)
# - 10 continuous predictors (ALS proxies)
#   * 3 are truly irrelevant (pure noise / independent of y)
#   * 2 are strong predictors but highly collinear with other strong predictors (redundancy)
# - Outcome rdnbr has a known relationship to a subset of predictors
# - CRS in meters (UTM) so 2 km buffering makes sense
 
library(terra)


set.seed(12498)
output <- 'data/temp/synth_train.tif'

set.seed(98238)
output <- 'data/temp/synth_test.tif'

# ----------------------------
# 1) Even-grid raster (~1000 cells)
# ----------------------------
n_side <- 100
epsg <- 6339
crs_str <- paste0("EPSG:", epsg)

x_min <- 685976.3
x_max <- 714010.5
y_min <- 4308005
y_max <- 4340014

r0 <- rast(
  ncols = n_side, 
  nrows = n_side,
  xmin = x_min, xmax = x_max,
  ymin = y_min, ymax = y_max,
  crs  = crs_str
)

cell_ids <- 1:ncell(r0)
xy <- xyFromCell(r0, cell_ids)
x <- xy[, 1]; y <- xy[, 2]
n <- length(cell_ids)

# ----------------------------
# 2) Smooth spatial latent fields
# ----------------------------
dx <- max(1, x_max - x_min)
dy <- max(1, y_max - y_min)

s1 <- sin((x - x_min) / (0.14 * dx)) + cos((y - y_min) / (0.16 * dy))
s2 <- cos((x - x_min) / (0.22 * dx)) * sin((y - y_min) / (0.18 * dy))
s3 <- sin((x + y) / (0.30 * (dx + dy)))

s1 <- as.numeric(scale(s1 + rnorm(n, 0, 0.15)))
s2 <- as.numeric(scale(s2 + rnorm(n, 0, 0.15)))
s3 <- as.numeric(scale(s3 + rnorm(n, 0, 0.15)))

mk <- function(z) as.numeric(scale(z))

# ----------------------------
# 3) Predictors: designed for variable selection diagnostics
# ----------------------------
# Strong "true" predictors (signal-bearing)
als_01 <- mk( 1.0*s1 + rnorm(n, 0, 0.55))                 # strong, mostly s1
als_02 <- mk( 0.8*s2 + 0.2*s1 + rnorm(n, 0, 0.65))        # strong, mostly s2
als_03 <- mk(-0.7*s2 + 0.3*s3 + rnorm(n, 0, 0.70))        # strong, mix s2/s3
als_04 <- mk( 0.7*s3 + rnorm(n, 0, 0.65))                 # strong, mostly s3

# Collinear "decoy twins" (highly correlated with true predictors but redundant)
# These should compete with their twins in RF importance / selection.
als_05 <- mk( 0.92*als_01 + rnorm(n, 0, 0.20))            # ~ very high corr w/ als_01
als_06 <- mk( 0.90*als_02 + rnorm(n, 0, 0.25))            # ~ very high corr w/ als_02

# Moderate predictor (useful but weaker / nonlinear)
als_07 <- mk( 0.35*s1 + 0.25*s3 + rnorm(n, 0, 1.00))      # weaker, mixed

# Irrelevant predictors (should be dropped / low importance)
als_08 <- mk(rnorm(n, 0, 1.0))                            # pure IID noise
als_09 <- mk(rnorm(n, 0, 1.0))                            # pure IID noise
als_10 <- mk(0.05*s1 + rnorm(n, 0, 1.30))                 # effectively noise relative to signal

# ----------------------------
# 4) Known ground-truth relationship (only uses a subset)
# ----------------------------
hinge <- function(z) pmax(z - 0.3, 0)

# IMPORTANT:
# - Uses als_01, als_02, als_03, als_04, als_07
# - DOES NOT use als_05/als_06 (redundant twins)
# - DOES NOT use als_08/als_09/als_10 (irrelevant)
rdnbr_true <-
  150 +
  70*als_01 +                 # strong linear
  -55*als_03 +                # strong linear (negative)
  35*als_04 +                 # moderate linear
  28*als_01*als_02 +          # interaction
  -20*(als_02^2) +            # nonlinearity
  18*sin(pi*als_07) +         # weaker nonlinear signal
  20*hinge(als_04)            # threshold effect

# Residual spatial structure + noise (keeps spatial CV meaningful)
spatial_resid <- 16*s1 - 9*s2 + rnorm(n, 0, 7)
rdnbr <- rdnbr_true + spatial_resid + rnorm(n, 0, 25)

# ----------------------------
# 5) Write to raster stack
# ----------------------------
make_layer <- function(vals, template, name) {
  r <- template
  values(r) <- vals
  names(r) <- name
  r
}

r_stack <- c(
  make_layer(als_01, r0, "als_01_trueStrong"),
  make_layer(als_02, r0, "als_02_trueStrong"),
  make_layer(als_03, r0, "als_03_trueStrong"),
  make_layer(als_04, r0, "als_04_trueStrong"),
  make_layer(als_05, r0, "als_05_collinearTwin01"),
  make_layer(als_06, r0, "als_06_collinearTwin02"),
  make_layer(als_07, r0, "als_07_trueWeakNonlinear"),
  make_layer(als_08, r0, "als_08_irrelevant"),
  make_layer(als_09, r0, "als_09_irrelevant"),
  make_layer(als_10, r0, "als_10_irrelevantish"),
  make_layer(rdnbr_true, r0, "rdnbr_true"),
  make_layer(rdnbr, r0, "rdnbr")
)

# Optional: table view for modeling outside raster ops
dat <- data.frame(
  cell = cell_ids,
  x = x, y = y,
  epsg = epsg,
  als_01, als_02, als_03, als_04, als_05,
  als_06, als_07, als_08, als_09, als_10,
  rdnbr_true = rdnbr_true,
  rdnbr = rdnbr
)

# ----------------------------
# 6) Diagnostics you can run immediately
# ----------------------------
# Correlation check for the collinear twins (should be high)
cor_01_05 <- cor(dat$als_01, dat$als_05)
cor_02_06 <- cor(dat$als_02, dat$als_06)

cat("corr(als_01, als_05) =", round(cor_01_05, 3), "\n")
cat("corr(als_02, als_06) =", round(cor_02_06, 3), "\n")

# Sanity: irrelevant vars should be ~uncorrelated with rdnbr_true
cat("corr(als_08, rdnbr_true) =", round(cor(dat$als_08, dat$rdnbr_true), 3), "\n")
cat("corr(als_09, rdnbr_true) =", round(cor(dat$als_09, dat$rdnbr_true), 3), "\n")
cat("corr(als_10, rdnbr_true) =", round(cor(dat$als_10, dat$rdnbr_true), 3), "\n")

# ----------------------------
# 7) Save outputs (optional)
# ----------------------------
writeRaster(r_stack, output, overwrite = TRUE)

print(r_stack)
print(crs(r_stack))
summary(dat$rdnbr)
