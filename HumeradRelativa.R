#######################################################################
############ Homogenization using the R Climatol package############
#######################################################################

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("blocklength","openxlsx","readr", "tidyverse", "climatol"
              , "writexl","readxl", "visdat","boot", "reshape","ggplot2")
ipak(packages)

# Lectura de la base de datos 
datos <- as.matrix(read.xlsx("BasesHumedad/Humedad.xlsx"
                             ,na.strings = "NA",rowNames = F))
estaciones <- read.xlsx("BasesHumedad/Coordenadas.xlsx")

# Ficheros de entrada
write(datos,'BasesHumedad/HR_2015-2017.dat')
write.table(estaciones, 'BasesHumedad/HR_2015-2017.est'
            , row.names=FALSE, col.names=FALSE)

# Homogenización de los datos Analisis Exploratorio diario
homogen('BasesHumedad/HR',2015,2017,expl=T)

# Convertir datos diarios a mensuales
dd2m(varcli="BasesHumedad/HR",2015,2017)

# Exploratorio de los datos mensules
homogen("BasesHumedad/HR-m",2015,2017,expl = TRUE)

# Exploratorio de los datos mensuales con valores obtenidos 
# de las anomalias normalizadas y de los 
# SNHT maximos para ventanas y para series completas.
homogen("BasesHumedad/HR-m",2015,2017,dz.min = -3.5, dz.max = 3
        , snht=0, snht2 = 8,std=2, vmin = 0)

# Ajustes de los datos diarios a partir de los mensuales
homogen("BasesHumedad/HR",2015,2017,dz.min = -6.5, dz.max = 3.5
        , vmin = 0, metad=TRUE)

# Análisis exploratorio de los datos
load("BasesHumedad/HR_2015-2017.rda")
View(est.c)

# Series homogéneas
dahstat('BasesHumedad/HR', 2015, 2017, stat='series')

#######################################################################
### Homo. by Climatol with the calculated SNHT calculated with MBB ###
#######################################################################

# Lectura de la base de datos con series homogenizadas
HR_2015_2017_series <- read_csv("C:/Users/Hp/Desktop/Tratamientos de R/R/Bootstrap estadístico/SerieHomogeneas/HR_2015-2017_series.csv")
HR1 <- rename(HR_2015_2017_series, c(E05="Alao",E09="Atillo",E06="Cumanda",E01="Espoch", E08="Matus",E07="Multitud",E04="Tixan",E03="Tunshi", E02="Urbina",`E02-2`="Urbina2"))

# Grafico de las series homogeneas
Numero <- c(1:1095)
HR <- select(HR1,-Date)
df <- cbind(Numero,HR)
df <- melt(df,id.vars = "Numero")
ggplot(df, aes(x = Numero, y = value, color = variable)) + geom_line()+labs(x="Observaciones",y="Humedad Relativa") + ggtitle("Series metorológicas de las 9 estaciones homogeneas") + theme(plot.title = element_text(hjust = 0.5))

# Creación del estadístico de prueba
Estadistico <- function(base){
  n <- length(base)
  A <- vector()
  B <- vector()
  T_a <- vector()
  for (i in 1:(n-1)){
    A[i] <- i*mean(base[1:i])*mean(base[1:i])
    B[i] <- (n-i)*mean(base[(i+1):n])*mean(base[(i+1):n])
    T_a[i] <- A[i]+B[i]
  }
  return(max(T_a))
}


#######################################################################
#################### Bootstrap por bloques moviles ####################
#######################################################################


########################### Estacion Alao ############################

Alao <- select(HR1,Alao)
Alao <- ts(Alao)

# Normalización de la serie
Alao <- (Alao-mean(Alao))/sd(Alao)

# Longuitud optima del bloque serie Alao
lbAlao <-  hhj(Alao, sub_sample = 10)

# Muestras MBB 
MBBAlao <- tsboot(Alao, Estadistico, R = 1000, l = 5, sim ="fixed" )
par(bg = "gray")
hist(MBBAlao$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Alao",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBAlao$t, c(0.95)), col= "blue",lwd=3)
dz <- density(MBBAlao$t)
lines(dz, col = "red", lwd = 2)
TAlao1 <- quantile(MBBAlao$t,c(0.95))
TAlao1
print(MBBAlao)

########################### Estacion Atillo ###########################

Atillo <- select(HR1,Atillo)
Atillo <- ts(Atillo)

# Normalización de la serie
Atillo <- (Atillo-mean(Atillo))/sd(Atillo)

# Longuitud optima del bloque serie Alao
lbAtillo <- hhj(Atillo,sub_sample = 10)

# Muestras MBB
MBBAtillo <- tsboot(Atillo, Estadistico, R = 1000, l = 10, sim ="fixed" )
par(bg = "gray")
hist(MBBAtillo$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie Atillo",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBAtillo$t, c(0.95)), col= "blue")
dz1 <- density(MBBAtillo$t)
lines(dz1, col = "red", lwd = 1)
TAtillo1 <- quantile(MBBAtillo$t,c(0.95))
TAtillo1
print(MBBAtillo)

########################## Estacion Cumanda ###########################

Cumanda <- select(HR1,Cumanda)
Cumanda <- ts(Cumanda)

# Normalización de la serie
Cumanda <- (Cumanda-mean(Cumanda))/sd(Cumanda)

# Longuitud optima del bloque serie 
lbCumanda <- hhj(Cumanda,sub_sample = 10)

# Muestras MBB
MBBCumanda<- tsboot(Cumanda, Estadistico, R = 1000, l = 5, sim ="fixed" )
hist(MBBCumanda$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Cumanda",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBCumanda$t, c(0.95)), col= "blue",lwd=3)
dz2 <- density(MBBCumanda$t)
lines(dz2, col = "red", lwd = 2)
TCumanda1 <- quantile(MBBCumanda$t,c(0.95))
TCumanda1
print(MBBCumanda)


########################### Estacion Espoch ###########################

Espoch <- select(HR1,Espoch)
Espoch <- ts(Espoch)

# Normalización de la serie
Espoch <- (Espoch-mean(Espoch))/sd(Espoch)


# Longuitud optima del bloque serie
lbEspoch <- hhj(Espoch, sub_sample = 10)

# Muestras MBB
MBBEspoch <- tsboot(Espoch, Estadistico, R = 1000, l = 5, sim ="fixed" )
hist(MBBEspoch$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Espoch",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBEspoch$t, c(0.95)), col= "blue",lwd=3)
dz3 <- density(MBBEspoch$t)
lines(dz3, col = "red", lwd = 2)
TEspoch1 <- quantile(MBBEspoch$t,c(0.95))
TEspoch1
print(MBBEspoch)

########################### Estacion Matus ############################

Matus <- select(HR1,Matus)
Matus <- ts(Matus)

# Normalización de la serie
Matus <- (Matus-mean(Matus))/sd(Matus)


# Longuitud optima del bloque serie Alao
lbMatus <- hhj(Matus,sub_sample = 10)

# Muestras MBB
MBBMatus <- tsboot(Matus, Estadistico, R = 1000, l = 5, sim ="fixed" )
hist(MBBMatus$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Matus",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBMatus$t, c(0.95)), col= "blue",lwd=3)
dz4 <- density(MBBMatus$t)
lines(dz4, col = "red", lwd = 2)
TMatus1 <- quantile(MBBMatus$t,c(0.95))
TMatus1
print(MBBMatus)

########################## Estacion Multitud ##########################

Multitud <- select(HR1,Multitud)
Multitud <- ts(Multitud)

# Normalización de la serie
Multitud <- (Multitud-mean(Multitud))/sd(Multitud)


# Longuitud optima del bloque serie Alao
lbMultitud<- hhj(Multitud,sub_sample = 10)

# Muestras MBB
MBBMultitud <- tsboot(Multitud, Estadistico, R = 1000, l = 10, sim ="fixed" )
hist(MBBMultitud$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Multitud",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBMultitud$t, c(0.95)), col= "blue",lwd=3)
dz5 <- density(MBBMultitud$t)
lines(dz5, col = "red", lwd = 1)
TMultitud1 <- quantile(MBBMultitud$t,c(0.95))
TMultitud1 
print(MBBMultitud)

############################ Estacion Tixan ###########################

Tixan <- select(HR1,Tixan)
Tixan <- ts(Tixan)

# Normalización de la serie
Tixan <- (Tixan-mean(Tixan))/sd(Tixan)

# Longuitud optima del bloque serie Alao
lbTixan <- hhj(Tixan,sub_sample = 10)

# Muestras MBB
MBBTixan <- tsboot(Tixan, Estadistico, R = 1000, l = 5, sim ="fixed" )
hist(MBBTixan$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Tixan",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBTixan$t, c(0.95)), col= "blue",lwd=3)
dz6 <- density(MBBTixan$t)
lines(dz6, col = "red", lwd = 2)
TTixan1 <- quantile(MBBTixan$t,c(0.95))
TTixan1
print(MBBTixan)

########################### Estacion Tunshi ###########################

Tunshi <- select(HR1,Tunshi)
Tunshi <- ts(Tunshi)

# Normalización de la serie
Tunshi <- (Tunshi-mean(Tunshi))/sd(Tunshi)

# Longuitud optima del bloque serie Tunshi
lbTunshi<- hhj(Tunshi, sub_sample = 10)

# Muestras MBB
MBBTunshi <- tsboot(Tunshi, Estadistico, R = 1000, l = 7, sim ="fixed" )
hist(MBBTunshi$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie Tunshi",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBTunshi$t, c(0.95)), col= "blue",lwd=3)
dz7 <- density(MBBTunshi$t)
lines(dz7, col = "red", lwd = 2)
TTunshi1 <- quantile(MBBTunshi$t,c(0.95))
TTunshi1
print(MBBTunshi)

########################### Estacion Urbina ###########################

Urbina2 <- select(HR1,Urbina2)
Urbina2 <- ts(Urbina2)

# Normalización de la serie
Urbina2 <- (Urbina2-mean(Urbina2))/sd(Urbina2)

# Longuitud optima del bloque serie Urbina
lbUrbina2<- hhj(Urbina2, sub_sample = 10)

# Muestras MBB
MBBUrbina2 <- tsboot(Urbina2, Estadistico, R = 1000, l = 5, sim ="fixed" )
hist(MBBUrbina2$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Urbina",main ="Distribución del Estadístico de las muestras MBB",ylab = "Densidad")
abline(v=quantile(MBBUrbina2$t, c(0.95)), col= "blue",lwd=3)
dz6 <- density(MBBUrbina2$t)
lines(dz7, col = "red", lwd = 2)
TUrbina1 <- quantile(MBBUrbina2$t,c(0.95))
TUrbina1
print(MBBUrbina2)


##### Hom. fijando el SNHT obtenido del promedio de las estaciones ####

# Promedio del umbral de correccion de inhomogeneidades de todas las estaciones
SNHT1 <- mean(c(TAlao1,TAtillo1,TCumanda1,TEspoch1,TMatus1,TMultitud1,TTixan1,TTunshi1,TUrbina1)) 
SNHT1

# Lectura y transformación de las Base de datos  
Datos <- read_xlsx("Propuesta1/Humedad.xlsx")
Datos$Espoch <- as.numeric(Datos$Espoch)
Datos$Cumanda <- as.numeric(Datos$Cumanda)
Datos$Multitud <- as.numeric(Datos$Multitud)
Datos$Tunshi <- as.numeric(Datos$Tunshi)
Base <- as.matrix(Datos)
estaciones <- read_xlsx("Propuesta1/Coordenadas.xlsx")

# Ficheros de entrada
write(Base,'Propuesta1/HRP1_2015-2017.dat')
write.table(estaciones, 'Propuesta1/HRP1_2015-2017.est', row.names=FALSE, col.names=FALSE)

# Aálisis exploratorio
homogen("Propuesta1/HRP1",2015,2017,expl = TRUE)

# Homogenizacion de las series fijando el umbral e correccion de inhomogeneidades SNHT
homogen("Propuesta1/HRP1",2015,2017, dz.min=-3.5, dz.max = 3.5, snht1 = SNHT1,vmin = 0,vmax = 100)

# Resumen estadístico
load('Propuesta1/HRP1_2015-2017.rda')
View(est.c)

# Series homogeneizadas 
dahstat('Propuesta1/HRP1', 2015, 2017, stat='series')


#######################################################################
####################### Bootstrap Estacionario ########################
#######################################################################

############################ Estacion Alao ############################

# Muestras SB 
SBAlao <- tsboot(Alao, Estadistico, R = 1000, l = 5, sim ="geom" )
par(bg = "gray")
hist(SBAlao$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Alao",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBAlao$t, c(0.95)), col= "blue",lwd=3)
dz <- density(SBAlao$t)
lines(dz, col = "red", lwd = 2)
TAlao <- quantile(SBAlao$t,c(0.95))
print(SBAlao)

########################### Estacion Atillo ###########################

# Muestras SB
SBAtillo <- tsboot(Atillo, Estadistico, R = 1000, l = 10, sim ="geom" )
par(bg = "gray")
hist(SBAtillo$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie Atillo",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBAtillo$t, c(0.95)), col= "blue",lwd=3)
dz1 <- density(SBAtillo$t)
lines(dz1, col = "red", lwd = 2)
TAtillo <- quantile(SBAtillo$t,c(0.95))
print(SBAtillo)

########################## Estacion Cumanda ###########################

# Muestras SB
SBCumanda<- tsboot(Cumanda, Estadistico, R = 1000, l = 5, sim ="geom" )
hist(SBCumanda$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Cumanda",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBCumanda$t, c(0.95)), col= "blue",lwd=3)
dz2 <- density(SBCumanda$t)
lines(dz2, col = "red", lwd = 2)
TCumanda <- quantile(SBCumanda$t,c(0.95))
print(SBCumanda)

########################### Estacion Espoch ###########################

# Muestras SB
SBEspoch <- tsboot(Espoch, Estadistico, R = 1000, l = 5, sim ="geom" )
hist(SBEspoch$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Espoch",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBEspoch$t, c(0.95)), col= "blue",lwd=3)
dz3 <- density(SBEspoch$t)
lines(dz3, col = "red", lwd = 2)
TEspoch <- quantile(SBEspoch$t,c(0.95))
print(SBEspoch)

########################### Estacion Matus ############################

# Muestras SB
SBMatus <- tsboot(Matus, Estadistico, R = 1000, l = 5, sim ="geom" )
hist(SBMatus$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Matus",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBMatus$t, c(0.95)), col= "blue",lwd=3)
dz4 <- density(SBMatus$t)
lines(dz4, col = "red", lwd = 2)
TMatus <- quantile(SBMatus$t,c(0.95))
print(SBMatus)

########################## Estacion Multitud ##########################

# Muestras SB
SBMultitud <- tsboot(Multitud, Estadistico, R = 1000, l = 10, sim ="geom" )
hist(SBMultitud$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Multitud",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBMultitud$t, c(0.95)), col= "blue",lwd=3)
dz5 <- density(SBMultitud$t)
lines(dz5, col = "red", lwd = 1)
TMultitud <- quantile(SBMultitud$t,c(0.95))
print(SBMultitud)

############################ Estacion Tixan ###########################

# Muestras SB
SBTixan <- tsboot(Tixan, Estadistico, R = 1000, l = 5, sim ="geom" )
hist(SBTixan$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Tixan",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBTixan$t, c(0.95)), col= "blue",lwd=3)
dz6 <- density(SBTixan$t)
lines(dz6, col = "red", lwd = 2)
TTixan <- quantile(SBTixan$t,c(0.95))
print(SBTixan)

########################### Estacion Tunshi ###########################

# Muestras SB
SBTunshi <- tsboot(Tunshi, Estadistico, R = 1000, l = 8, sim ="geom" )
hist(SBTunshi$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie Tunshi",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBTunshi$t, c(0.95)), col= "blue",lwd=3)
dz7 <- density(SBTunshi$t)
lines(dz7, col = "red", lwd = 2)
TTunshi <- quantile(SBTunshi$t,c(0.95))
print(SBTunshi)

########################### Estacion Urbina ###########################

# Muestras SB
SBUrbina2 <- tsboot(Urbina2, Estadistico, R = 1000, l = 8, sim ="geom" )
hist(SBUrbina2$t, prob=T,breaks = 50,col = "pink3",xlab="Estimación del estadístico para la serie de tiempo Urbina",main ="Distribución del Estadístico de las muestras SB",ylab = "Densidad")
abline(v=quantile(SBUrbina2$t, c(0.95)), col= "blue",lwd=3)
dz6 <- density(SBUrbina2$t)
lines(dz7, col = "red", lwd = 2)
TUrbina2 <- quantile(SBUrbina2$t,c(0.95))
print(SBUrbina2)

##### Hom. fijando el SNHT obtenido del promedio de las estaciones ####

# Promedio del umbral de corrección de inhomogeneidades de todas las estaciones
SNHT2 <- mean(c(TAlao,TAtillo,TCumanda,TMatus,TMultitud,TTixan,TTunshi,TUrbina2))
SNHT2

# Lectura y transformación de la base de datos
Datos1 <- read_xlsx("Propuesta1SB/Humedad.xlsx")
Datos1$Espoch <- as.numeric(Datos1$Espoch)
Datos1$Cumanda <- as.numeric(Datos1$Cumanda)
Datos1$Multitud <- as.numeric(Datos1$Multitud)
Datos1$Tunshi <- as.numeric(Datos1$Tunshi)
Base1 <- as.matrix(Datos1)
estaciones1 <- read_xlsx("Propuesta1SB/Coordenadas.xlsx")

# Ficheros de entrada
write(Base1,'Propuesta1SB/HRP1SB_2015-2017.dat')
write.table(estaciones1, 'Propuesta1SB/HRP1SB_2015-2017.est', row.names=FALSE, col.names=FALSE)

# Análisis exploratorio
homogen("Propuesta1SB/HRP1SB",2015,2017,expl = TRUE)

# Homogenización de las series fijando el umbral de correccion de inhomogeneidades SNHT
homogen("Propuesta1SB/HRP1SB",2015,2017, dz.min=-3.5, dz.max = 3.5, snht1 = SNHT2,vmin = 0,vmax = 100)

# Resumen estadistico
load('Propuesta1SB/HRP1SB_2015-2017.rda')
View(est.c)
# Series homogeneizadas y 
dahstat('Propuesta1SB/HRP1SB', 2015, 2017, stat='series')

# Gráficos de las series Homogeneas

Base <- read_csv("BasesHumedad/HR_2015-2017_series.csv")
BaseMMB <- read_csv("Propuesta1/HRP1_2015-2017_series.csv")
BaseSB <- read_csv("Propuesta1SB/HRP1SB_2015-2017_series.csv")
Fecha <- as.Date(Base$Date)
Urbina1 <- Base$`E02-2`
Urbina2 <- BaseMMB$`E02-2`
Urbina3 <- BaseSB$`E02-2`
Urbina <- data.frame(Fecha,Urbina1,Urbina3,Urbina2)
Urbina  
str(Urbina)
Urbina <- melt(Urbina,id.vars = "Fecha")
ggplot(Urbina, aes(x = Fecha, y = value, color = variable)) + geom_line()+labs(x="Observaciones",y="Humedad Relativa") + ggtitle("Series metorológicas de las 9 estaciones homogeneas") + theme(plot.title = element_text(hjust = 0.5))

Tunshi1 <- Base$E03
Tunshi2 <- BaseMMB$`E03-2`
Tunshi3 <- BaseSB$`E03-2`


Tunshi <- data.frame(Fecha,Tunshi1,Tunshi2,Tunshi3)
Tunshi <- melt(Tunshi,id.vars = "Fecha")
ggplot(Tunshi, aes(x = Fecha, y = value, color = variable)) + geom_line()+labs(x="Observaciones",y="Humedad Relativa") + ggtitle("Series metorológicas de las 9 estaciones homogeneas") + theme(plot.title = element_text(hjust = 0.5))
