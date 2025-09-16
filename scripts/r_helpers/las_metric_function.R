library(tidyverse)
library(glue)

las_cld_metrics <- function(z, r) {
  
  # ----------------------------- Height metrics -------------------------------
  
  height_metrics = list(
    z_max = max(z, na.rm = TRUE),
    z_mean = mean(z, na.rm = TRUE))
  
  for (i in seq(0, 0.95, 0.05)) {
    height_metrics[glue('z_p{i*100}')] <- quantile(z,
                                                 probs = i,
                                                 names = FALSE,
                                                 na.rm = TRUE)
    
  }
  
  # ---------------------------- Variation metrics -----------------------------
  
  var_metrics = list(
      z_sd = sd(z, na.rm = TRUE),
      z_cv = sd(z, na.rm = TRUE) / mean(z, na.rm = TRUE),
      z_kurtosis = moments::kurtosis(z, na.rm = TRUE),
      z_skewness = moments::skewness(z, na.rm = TRUE),
      z_iqr = IQR(z, na.rm = TRUE),
      canopy_relief_ratio = (height_metrics$z_mean - height_metrics$z_p0) / 
        (height_metrics$z_max - height_metrics$z_p0)
    )
  
  # --------------------------- Vertical densities -----------------------------
  
  n_bins = 10
  h_thresh = 1
  
  bins = seq(h_thresh, 
             height_metrics$z_max, 
             (height_metrics$z_max - h_thresh) / n_bins)
  
  vertical_density = list()
  
  for (i in 1:(length(bins) - 1)) {
    val <- sum(z >= bins[i] & z < bins[i + 1]) / sum(z >= h_thresh)
    vertical_density[glue('d{str_pad(i, 2, pad = 0)}')] <- val 
  }

  # ------------------------------ Canopy cover --------------------------------
  
  canopy_cover = list(
    cc_1m = sum(z > 1) / length(z),
    cc_mean = sum(z > mean(z, na.rm = T)) / length(z)
  )
  
  # ------------------------------- Ladder fuels -------------------------------
  
  lf_min = 1
  lf_max = 8
  
  ladder_fuels <- list()
  
  for (i in seq(lf_min + 1, lf_max)) {
    
    val <- sum(z >= i - 1 & z < i) / sum(z < i)
    ladder_fuels[glue('lf_{i - 1}{i}')] <- val
    
  }

  # ------------------------------- First returns ------------------------------ 

  z1r <- z[r == 1]
  
  # four height bands: 4-8m, 8-16m, 16-32m, and 32+m 
  
  first_return_metrics = list(
    cc1r_48 = sum(z1r >= 4 & z1r < 8) / sum(z1r < 8),
    cc1r_816 = sum(z1r >= 8 & z1r < 16) / sum(z1r < 16),
    cc1r_1632 = sum(z1r >= 16 & z1r < 32) / sum(z1r < 32),
    cc1r_gt32 = sum(z1r >= 32) / length(z1r))
  
  # cc1r = list_c(first_return_metrics)
  # 
  # if (sum(cc1r == 0) == 4) {
  # 
  #   cc1r = rep(-9999, 4)
  # 
  # } else {
  # 
  #   cc1r = cc1r[1:max(which(cc1r > 0))]
  # 
  # }
  # 
  # val = ((sum((cc1r/sum(cc1r)) ^ 2) * length(cc1r)) ^ -1) * mean(cc1r)
  
  first_return_metrics['cc1r_3m'] = sum(z1r > 1) / length(z1r)
  
  # first_return_metrics['lddr_evn'] = val
  
  # --------------------------------- Return ---------------------------------- 
  
  return(
    c(
      height_metrics,
      var_metrics,
      vertical_density,
      canopy_cover,
      ladder_fuels,
      first_return_metrics
  ))

}