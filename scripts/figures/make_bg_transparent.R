#make image bg transparant 

#you will need magic:
#install.packages("magick")#
library(magick)

setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/surfaces/")
list <- list.files(pattern = '^analysis') #get list of file names you want

for (x in list) {
  
  pic  <- image_read(x)
  tpic <- image_transparent(pic, 'white')
  tpic_c <- image_crop(tpic, "810x700+250")
  tpic_c <- image_trim(tpic)
  image_write(tpic_c, path = paste0("t_",x), format = "png") # new file is t_$file_name
}

#flip subcortex images
setwd("~/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/")
list <- list.files(pattern = '^t_')
for (x in list) {
  
  pic  <- image_read(x)
  tpic <- image_flop(pic)
  #tpic_c <- image_crop(tpic, "810x700+250")
  image_write(tpic, path = paste0("f_",x), format = "png")
}
