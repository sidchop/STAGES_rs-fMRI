remotes::install_github("LCBC-UiO/ggsegSchaefer")
library(ggsegSchaefer)
library(ggseg)

remotes::install_github("LCBC-UiO/ggseg3D")
library(ggseg3d)
x <- schaefer7_3d
x$ggseg_3d[[3]][[2]][2:32] <- "#9C755FFF"
x$ggseg_3d[[4]][[2]][2:32] <- "#9C755FFF"

x$ggseg_3d[[3]][[2]][33:69] <- "#F28E2BFF"
x$ggseg_3d[[4]][[2]][33:69] <- "#F28E2BFF"

x$ggseg_3d[[3]][[2]][70:92] <- "#E15759FF"
x$ggseg_3d[[4]][[2]][70:92] <- "#E15759FF"

x$ggseg_3d[[3]][[2]][93:114] <- "#76B7B2FF"
x$ggseg_3d[[4]][[2]][93:114] <- "#76B7B2FF"

x$ggseg_3d[[3]][[2]][115:127] <- "#59A14FFF"
x$ggseg_3d[[4]][[2]][115:127] <- "#59A14FFF"

x$ggseg_3d[[3]][[2]][128:149] <- "#EDC948FF"
x$ggseg_3d[[4]][[2]][128:149] <- "#EDC948FF"

x$ggseg_3d[[3]][[2]][150:201] <- "#B07AA1FF"
x$ggseg_3d[[4]][[2]][150:201] <- "#B07AA1FF"


ggseg3d(atlas = x, surface = "LCBC") %>% 
  pan_camera("right medial") %>%
  remove_axes()
