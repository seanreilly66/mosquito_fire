# ==============================================================================
#
# Plotting themes and functions
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 16 Dec 2024
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

# ==============================================================================
# ================================= Plot theme =================================
# ==============================================================================

theme_set(
  theme(
    text = element_text(family = 'serif', face = 'plain'),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    line = element_line(linewidth = 1),
    axis.line = element_line(),
    panel.background = element_rect(color = 'white'),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    legend.spacing = unit(0, "cm"),
    legend.margin = margin(0, 5, 0, 5),
    title = element_text(size = 12.8),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, 
                              hjust = 0, 
                              margin = margin(l = 0,
                                              b = 2))
  )
)


# ==============================================================================
# ============================ Plot helper function ============================
# ==============================================================================

bin_to_val = function(bin) {
  
  mat = bin %>%
    str_remove_all(pattern = '[^[:digit:],\\.]') %>%
    str_split(',', simplify = T)
  
  bin_df = tibble(
    bin_min = mat[, 1],
    bin_max = mat[, 2]) %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(
      bin_mid = mean(c_across(everything())))
  
  return(bin_df)
  
} 

# ==============================================================================
# =============================== Color palettes ===============================
# ==============================================================================

trtmnt_col_pal = c(
  'Control' = 'grey20',
  'Fire' = '#CC6677',
  'Mech' = '#537747',
  'Mech + Fire' = '#557EA2'
)

trtmnt_fill_pal = c(
  'Control' = 'white',
  'Fire' = '#CC6677',
  'Mech' = '#537747',
  'Mech + Fire' = '#557EA2'
)

# ==============================================================================