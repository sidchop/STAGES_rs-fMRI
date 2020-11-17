#make colour bar for surf and subcort renderings 

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
 # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1, cex.axis=3)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

library(RColorBrewer)
list(color = brewer.pal(9, "Reds"))
max.val <- c(63, 63, 58, 58, 21, 21, 6, 6, 21, 21, 12, 12, 5, 5)
ticks <- list(c(1, 15, 30, 45, 63), 
              c(1, 15, 30, 45, 63),
              c(1, 15, 30, 45, 58), 
              c(1, 15, 30, 45, 58), 
              c(1,10,20,30,41),
              c(1,10,20,30,41),
              c(1,3,6),
              c(1,3,6),
              c(1,5,10,15,21),             
              c(1,5,10,15,21),
              c(1,5,12),
              c(1,5,12),
              c(1,3,5),
              c(1,3,5))

setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/colourbars/")
for (m in 1:length(max.val)){
  tiff(filename = paste0("colourbar_",m,".tiff"), width = 150, height = 600)
  color.bar(lut = colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(100), 
            min = 1, 
            max = max.val[m], 
            ticks = ticks[[m]])
  dev.off()
}



