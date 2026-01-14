# ==============================================================================
#
# Burn serverity raster variogram
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 
# Last commit: 
#
# Status: Needs documentation
#
# ==============================================================================
#
# Description:
#
# ==============================================================================
#
# User inputs:
#
# ==============================================================================
#
# Package dependencies:
#
# ==============================================================================
#
# Known problems:
#
# ==============================================================================

library(tidyverse)
library(sf)
library(gstat)
library(terra)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

burn_rast_file <- 'data/burn_severity/mosquito_rdnbr_mtbs.tif'

study_gpkg <- 'data/mosquito_fire.gpkg'
study_perim_layer <- 'mosquito_study_area'

burn_metric <- 'rdnbr'


# ==============================================================================
# ================================== Data prep =================================
# ==============================================================================

burn_rast <- rast(burn_rast_file)
study_perim <- st_read(study_gpkg, study_perim_layer) %>%
  st_zm

clip_rast <- crop(burn_rast, study_perim) %>%
  mask(study_perim)

# Aggregate cells to reduce to feasible processing size

target_n <- 20000L
fact <- max(1L, ceiling(sqrt(ncell(clip_rast) / target_n)))
agg_rast <- aggregate(clip_rast, fact = fact, fun = mean, na.rm = TRUE)



# Convert to sf df

df <- as.data.frame(clip_rast, xy = TRUE, na.rm = TRUE)
pts <- st_as_sf(df, coords = c("x","y"), crs = crs(agg_rast))
names(pts)[1] <- burn_metric

# Add numeric coords to remove landscape spatial trends

pts <- pts %>%
  mutate(.x = st_coordinates(geometry)[,1],
         .y = st_coordinates(geometry)[,2])

# ==============================================================================
# ================================== Variogram =================================
# ==============================================================================

bb   <- st_bbox(pts)
diag <- sqrt((bb$xmax-bb$xmin)^2 + (bb$ymax-bb$ymin)^2)
cutoff <- 0.5 * diag
width <- cutoff / 100

g <- gstat(
  id = 'rdnbr',
  formula = rdnbr ~ .x + .y,
  # formula = rdnbr ~ 1,
  data = st_drop_geometry(pts),
  locations = ~ .x + .y
  # ,
  # set = list(gls = 1)
)



vg_emp <- variogram(g, cutoff = cutoff, width = width, cressie = TRUE)

mdl = 'Exp'
mdl = 'Sph'
mdl = 'Mat'

v_init <- vgm(psill = max(vg_emp$gamma, na.rm=TRUE),
              model = mdl, 
              range = max(vg_emp$dist, na.rm=TRUE)/3,
              nugget = min(vg_emp$gamma, na.rm=TRUE))
v_fit <- fit.variogram(vg_emp, v_init)

plot(vg_emp, v_fit)

# Components
nugget       <- v_fit$psill[v_fit$model == "Nug"]
partial_sill <- sum(v_fit$psill[v_fit$model != "Nug"])
sill         <- nugget + partial_sill

# Effective range = distance where model reaches 95% of sill (model-agnostic)
dseq  <- data.frame(dist = seq(0, cutoff, length.out = 600))
vline <- variogramLine(v_fit, dist_vector = dseq$dist)
thr   <- 0.95 * sill
eff_range <- dseq$dist[min(which(vline$gamma >= thr))]

cat(sprintf("\nRdNBR variogram (%s): sill=%.3f (nugget=%.3f, partial=%.3f), eff_range=%.1f m\n",
            mdl, sill, nugget, partial_sill, eff_range))

#--- 4) Suggested spatial CV settings -----------------------------------------
# Blocks (square/hex): side/diameter ~ effective range
block_size_m   <- eff_range
buffer_radius_m<- eff_range  # for h-block/buffered CV

cat(sprintf("Suggested block side ≈ %.0f m; buffer radius ≈ %.0f m\n",
            block_size_m, buffer_radius_m))

















field_data <- read_csv('data/field/plot_field_metrics.csv')
field_loc <- read_sf('data/boundaries/ssu_3dforests.gdb', 'field_plots') %>%
  mutate(plot = as.numeric(plot)) %>%
  rename(x = Easting, y = Northing)

field_data <- left_join(field_data, field_loc)


var_names = c('biomass_sum', 'h_mean', 'cbh', 'densiometer_mean', 'cbd_mean', 'lai_mean')

field_data <- field_data %>%
  select(x, y, campaign, all_of(var_names)) 

var_labels = tibble(col_name = var_names,
                    label = c('Biomass', 'Mean height', 'CBH', 'CC', 'CBD', 'LAI'))

# Campaign filter

campaigns = unique(field_data$campaign)
field_data2 <- field_data

c_i = 2

field_data <- field_data2 %>%
  filter(campaign == campaigns[c_i])



# ============================== gstat variogram ===============================



plots_list <- list()
plot_label <- list()

i = 1
b_width = 100
mdl = 'Exp'





if (exists('gstat_fit')) rm('gstat_fit')

field_var <- var_names[i]
field_var

var_data <- field_data %>%
  filter(if_any(field_var, ~!is.na(.)))

gstat_var <- variogram(as.formula(glue('{field_var} ~ 1')),
                       locations = ~ x + y,
                       data = var_data,
                       cutoff = 5000,
                       width = b_width)

plot(gstat_var)

gstat_fit <- fit.variogram(gstat_var, vgm(500, mdl, 2000))
gstat_fit

gstat_fit_line <- variogramLine(gstat_fit, maxdist = 5000)

plot(gstat_var, model = gstat_fit)

plot_label[[field_var]] <- tibble(
  x = min(gstat_var$dist),
  y = max(gstat_var$gamma),
  label = glue('{var_labels %>%
  filter(col_name == field_var) %>%
  pull(label)}
Partial sill: {signif(gstat_fit$psill, 3)}
Range: {signif(gstat_fit$range, 3)}'))

plots_list[[field_var]] <- ggplot(mapping = aes(x = dist, y = gamma)) +
  geom_point(data = gstat_var) +
  geom_line(data = gstat_fit_line, color = 'firebrick') +
  labs(
    x = NULL,
    y = NULL
  ) +
  geom_vline(xintercept = gstat_fit$range, linetype = 'dashed') +
  geom_label(data = plot_label[[field_var]], aes(x = x, y = y, label = label), 
            hjust = 'inward', 
            vjust = 'inward', 
            family = 'serif', 
            fontface = 'plain',
            label.size = NA,
            size = 5)

i = i + 1



gstat_plot <- ggarrange(plotlist = plots_list, align = 'hv', ncol = 2, nrow = 3) %>%
  annotate_figure(left = text_grob('Semivariance', family = 'serif', size = 16, rot = 90),
                bottom = text_grob('Distance (m)', family = 'serif', size = 16),
                top = text_grob(glue('Model: {mdl}  Bin: {b_width}'), family = 'serif', size = 16))

gstat_plot


ggsave(gstat_plot,
       filename = glue('figures/gstat_variog_{mdl}_bin{b_width}.png'), width = 8.5, height = 10, units = 'in', dpi = 700)


# =============================== geoR variogram ===============================

# 
# field_geodata <- field_data %>%
#   select(x, y, all_of(field_var)) %>%
#   as.geodata(coords.col = 1:2)
# 
# plot(field_geodata)
# 
# br = seq(1, 5000, 100)
# 
# i_variog <- variog(
#   field_geodata,
#   breaks = br,
#   lamda = lmda,
#   # max.dist = 5000,
#   # pairs.min = 5,
#   # estimator.type = 'modulus'
# )
# 
# plot(i_variog)
# 
# i_ml <- variofit(
#   i_variog)
# 
# i_ml
# 
# plot(i_variog, main = 'geoR'); lines(i_ml, col='red')
# 
#   
