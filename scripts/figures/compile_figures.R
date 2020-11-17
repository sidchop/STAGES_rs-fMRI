#compile function, compiles plots for each effect


compile_figures <- function(brain, networkplot, circleplot, 
                            lh_m, lh_l, rh_m, rh_l, 
                            subcortl, subcortr, 
                            colorbar){
 surf <- cowplot::plot_grid(rh_l, lh_l,
                      rh_m, lh_m, 
                      nrow=2, ncol = 2)

  
#brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = -.01, y = 0.02 , scale = 0.5)
brainsl <- cowplot::ggdraw(surf) + cowplot::draw_plot(subcortl, x = -0.28, y = -0.35 , scale = 0.258)

brains <- cowplot::ggdraw(brainsl) + cowplot::draw_plot(subcortr, x = 0.28, y = -0.35 , scale = 0.258)

brains
#add colour bar
brains_cb <- cowplot::plot_grid(brains, cb, rel_widths = c(1, 0.1), rel_heights = c(1,0.3), 
                                ncol = 2)

ggsave(plot = brains_cb, filename = "temp", device = "png")

brains_x <- ggdraw() + draw_image("temp", scale = 1)
 # cp removed
all <- cowplot::plot_grid(brain, networkplot, brains_x, ncol = 3, nrow = 1, 
                          rel_widths = c(0.6,1,0.9), rel_heights = c(0.6,1,0.9))  
return(all)
  
}
