#make_ggseg3d

#grey matter subcortical 3d
source("~/Dropbox/Sid/R_files/functions/sourceFolder.R")
sourceFolder("~/Dropbox/Sid/R_files/functions/ggseg3d_sid", recursive = T)

library(ggseg)
library(magrittr)
library(dplyr)
library(tidyr)
library(scales)


#for (hemi in c("l", "r")) {

#max.vect <- c(66, 115, 80, 58, 41, 10, 6, 5, 21, 15)
max.vect <- c(115, 115, 80, 80, 41, 41, 6, 6, 21, 21)

for (a in 1:10){
  remove(aseg_3d)
  load("~/Dropbox/Sid/R_files/functions/ggseg3d_sid/data/aseg_3d.rda")
  aseg_3d <- filter(aseg_3d, surf == "LCBC" & hemi == "subcort")
  aseg_3d <- tidyr::unnest(aseg_3d)
  degree <- as.matrix(read.table(paste0("~/Dropbox/Sid/R_files/STAGES_fmri/data/fsaverage/degree/analysis_subc_3daseg_",a,".txt")))
  #    ifelse(hemi=="l", degree[16:32] <- NA, degree[1:15] <- NA)
  data <- mutate(aseg_3d, p =  degree)
  #remove NA regions
  aseg_3d[which(is.na(data$p)),] <- NA
  aseg_3d <- drop_na(aseg_3d)
  
  data[which(is.na(data$p)),] <- NA
data <- drop_na(data)

#make 0 valus NA so wew can set them as grey in ggseg3d

data$p[data$p==0]<-NA
  #max(data$p, na.rm = T) #print max for color bar
  
x <- ggseg3d(.data = data, atlas = aseg_3d, colour = "p", text = "p",
             palette = c("light yellow",
                          "orange",
                          "red",
                          "dark red"),
             min.colour = 0, max.colour = max.vect[a]
             
)
x <-  remove_axes(x)

scene=list(camera = list(eye = list(x = 0, y = 1, z = -2.25)),
           aspectratio = list(x=1.6,y=1.6,z=1.6))

y <- plotly::layout(x, 
               scene = scene, 
               plot_bgcolor  = "rgba(0, 0, 0, 0)",
               paper_bgcolor = "rgba(0, 0, 0, 0)", 
               width = 810, height = 810)

y

z <- plotly::hide_colorbar(y)




plotly::orca(z, file =  paste0("analysis_sub_",a,".png"))

#  ggseg3d(.data = data, atlas = aseg_3d, colour = "p", text = "p",
#          na.alpha= 0,  palette = c("white")
#          , remove.axes = T, camera = list(x = 5, y = 90, z = 60), show.legend=F)
#  
  
}
#}

#some_data<- mutate(aseg_3d, p =  as.matrix(read.table("analysis_subc_3daseg_2.txt")))
#max(some_data$p, na.rm = T)
#ggseg3d(.data = some_data, atlas = aseg_3d, colour = "p", text = "p",
#        na.alpha= 0,  palette = c("light yellow", "yellow", "orange", "red", "dark red"), camera = "medial",
#        glassbrain = 0, glassbrain_hemisphere = "right", remove.axes = T)
#
#camera = "medial"

#Make legend
library(ggseg3d)
library(tidyr)
library(grDevices)
  remove(aseg_3d)
  aseg_3d <- filter(aseg_3d, surf == "LCBC" & hemi == "subcort")
  aseg_3d <- unnest(aseg_3d)
  degree <- c(NA, NA, NA, NA, 1, 2, 3, 4, NA, NA, NA, 5,6,7, NA, NA, NA, NA, NA, 1,2,3,4,5,6,7,NA, NA, NA, NA, NA, NA)
  
  #    ifelse(hemi=="l", degree[16:32] <- NA, degree[1:15] <- NA)
  data <- mutate(aseg_3d, p =  degree)
  #max(data$p, na.rm = T) #print max for color bar
  library(paletteer)
  x <- ggseg3d::ggseg3d(.data = data, atlas = aseg_3d,colour = "p", palette = c("#9C755FFF", "#F28E2BFF","#E15759FF",
                                                                       "#76B7B2FF", "#59A14FFF", "#EDC948FF",
                                                                       "#B07AA1FF", "#5575b3FF"),
               na.alpha= 0, na.colour = "grey")
               

  x
  
  scene=list(camera = list(eye = list(x = 0, y = 1, z = -2.25)),
             aspectratio = list(x=2.1,y=2.1,z=2.1))
  
  y <- plotly::layout(x, 
                      scene = scene, 
                      plot_bgcolor  = "rgba(0, 0, 0, 0)",
                      paper_bgcolor = "rgba(0, 0, 0, 0)", 
                      width = 810, height = 810)
  
  z <- plotly::hide_colorbar(y)
  
  
  
  
  plotly::orca(z, file =  paste0("legend_sub_.png"))
  