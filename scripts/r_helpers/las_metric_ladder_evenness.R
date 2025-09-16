library(tidyverse)
library(glue)

las_ladder_evenness <- function(z, r) {
  
  # ------------------------------- First returns ------------------------------ 
  
  z1r <- z[r == 1]
  
  # four height bands: 4-8m, 8-16m, 16-32m, and 32+m 
  
  first_return_cover = list(
    cc1r_48 = sum(z1r >= 4 & z1r < 8) / sum(z1r < 8),
    cc1r_816 = sum(z1r >= 8 & z1r < 16) / sum(z1r < 16),
    cc1r_1632 = sum(z1r >= 16 & z1r < 32) / sum(z1r < 32),
    cc1r_gt32 = sum(z1r >= 32) / length(z1r))
  
  cc1r = list_c(first_return_cover)

  if (sum(is.na(cc1r))) {
    
    cc1r = rep(NaN, 4)
    
  } else if (sum(cc1r == 0) == 4) {

    cc1r = rep(-9999, 4)

  } else {

    cc1r = cc1r[1:max(which(cc1r > 0))]

  }

  val = ((sum((cc1r/sum(cc1r)) ^ 2) * length(cc1r)) ^ -1) * mean(cc1r)
  
  # --------------------------------- Return ---------------------------------- 
  
  return(
    list(
      ladder_evenness = val
    ))
  
}