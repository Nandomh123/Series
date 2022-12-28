library(visdat) # visualización de datos
library(naniar)
library(readxl)
Humedad <- read.xlsx("BasesHumedad/Humedad.xlsx"
                     ,na.strings = "NA",rowNames = F)


png( "Imagenes/fig_a.png", width = 1000, height = 600, 
     units = "px", pointsize = 9, bg = "white", res = 140 )
vis_miss(Humedad,sort_miss = F )
dev.off()
