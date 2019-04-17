#### Try to make spatial stuff happen. 
#the Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

#import our dataset (as date)
library(readxl)
Datasource <- read_excel("~/TestR/GenericRf/Datasource.xlsx", 
                         col_types = c("text", "text", "date", 
                                       "text", "numeric", "numeric", "text", 
                                       "numeric", "text"))
Dataset1<-Datasource #to keep this clean I'll assign it to something "working"


#_____4_1 looking at the Contours of max data. _____________________________
#
#https://stackoverflow.com/questions/11649580/how-to-convert-nad-83-coordinates-to-latitude-and-longitude-with-rgdal-package
#http://spatialreference.org/ref/esri/102311/
Maxvlaues<-Dataset1 %>% 
  unite(Well_Unit_Chem,Well,AquiferUnit,Chem,sep = "_") %>% 
  group_by(Well_Unit_Chem) %>% 
  filter(Value == max(Value)) %>% 
  filter(Date == max(Date)) %>% 
  separate(Well_Unit_Chem,c("Well","AquiferUnit","Chem"),sep = "_") %>% 
  arrange(Well)  

Maxvlaues<-Maxvlaues %>% subset(Chem=="TCE") %>% 
  subset(AquiferUnit=="Unit4") %>% #grabbing just unit4 to generate contours as a test
  subset(`Dup-Prim`=="Primary")
#View(Maxvlaues)

testPts <- transform(expand.grid(x=X,y=Y),z=Value)
v<-ggplot(Maxvlaues, aes(X,Y,z = Value))+
  geom_point(aes(colour=Value))+
  scale_color_gradient(low = "green4", high = "red",trans="log10")
V
#____________________________________________interpolating a grid
#note: limited, no extrapolation. linear interpolation only.
#(Spline cannot be fit to this type of data)
library(akima)
#  stat_density2d()

pts.grid <- interp(Maxvlaues$X, Maxvlaues$Y,Maxvlaues$Value)
minx<-min(Maxvlaues$X)-1000
maxx<-max(Maxvlaues$X)+1000
miny<-min(Maxvlaues$Y)-1000
maxy<-max(Maxvlaues$Y)+1000
pts.grid <- interp(Maxvlaues$X, Maxvlaues$Y,Maxvlaues$Value,
                   xo = seq(minx,maxx,length= 400),
                   yo = seq(miny,maxy,length= 400))

pts.grid2 <- expand.grid(x=pts.grid$x, y=pts.grid$y)
pts.grid2$z1 <- as.vector(pts.grid$z)
pts.grid2
ggplot(pts.grid2,aes(x,y))+
 # geom_point(aes(color=z1))+ #this adds points 
  stat_contour(data=na.omit(pts.grid2), binwidth=50, aes(x=x, y=y, z=z1,colour=..level..))+
  geom_point(data=Maxvlaues,aes(X,Y,colour=Value, size = Value))+# aes(colour=Value)+
  scale_color_gradient(low = "green4", high = "red",trans="log10")


#_____________________________________________________________________
#_____4_1 Thiessen Polygons and other interpolations _________________


#the Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

#import our dataset (as date)
library(readxl)
Datasource <- read_excel("~/TestR/GenericRf/Datasource.xlsx", 
                         col_types = c("text", "text", "date", 
                                       "text", "numeric", "numeric", "text", 
                                       "numeric", "text"))
Dataset1<-Datasource #to keep this clean I'll assign it to something "working"

Maxvlaues<-Dataset1 %>% 
  unite(Well_Unit_Chem,Well,AquiferUnit,Chem,sep = "_") %>% 
  group_by(Well_Unit_Chem) %>% 
  filter(Value == max(Value)) %>% 
  filter(Date == max(Date)) %>% 
  separate(Well_Unit_Chem,c("Well","AquiferUnit","Chem"),sep = "_") %>% 
  arrange(Well)  
#for now i'm just gonna grab one chem and one aquifer unit
Maxvlaues<-Maxvlaues %>% subset(Chem=="TCE") %>% 
  subset(AquiferUnit=="Unit4") %>% #grabbing just unit4 to generate contours as a test
  subset(`Dup-Prim`=="Primary")

Maxvlaues$logValue<-log10(Maxvlaues$Value)
Maxvlaues
#this may be not needed. 
###
library(sp)
library(sf)
library(rgdal)
library(tmap)
library(ggalt)
library(spatstat)  # Used for the dirichlet tessellation function
library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons


###OK THIS IS KEY!. now i have a real spatial set
####2872 is crs for California stplane in feet
p.sf<-sf::st_as_sf(Maxvlaues, coords = c("X", "Y"), crs = 2872)
#this is also very key. to use the examples i need stuff in spatial form, not as sf form. 
psf1<-as(p.sf,"Spatial")

#take a look a the data:
#tm_shape(p.sf)+tm_dots(col = "Value",palette="RdBu",size=0.7)

#need to make a ppp thing. this is a "point pattern" dataset, does not have spatial yet
#a  "window" is key
windowsX<-c(min(Maxvlaues$X)-1000,max(Maxvlaues$X)+1000)
WindowsY<-c(min(Maxvlaues$Y)-1000,max(Maxvlaues$Y)+1000)

pts.ppp<-ppp(Maxvlaues$X,Maxvlaues$Y,window=owin(windowsX,WindowsY))

#generate the polygons
th<-as(dirichlet(pts.ppp),"SpatialPolygons")
#and give it a spatial thing:
proj4string(th)<-proj4string(psf1)

## now to make a polygon to bound it
MV<-Maxvlaues
MV[,1:5]<-NULL
MV[,3:5]<-NULL
Convexhull1<-chull(MV)
coordTable<-as.matrix(MV[Convexhull1,])
library(mapview)
library(raster)
ch_bound<-coords2Polygons(coordTable,ID="A")
proj4string(ch_bound)<-proj4string(psf1)
chbound<-buffer(ch_bound, width = 500)
proj4string(chbound)<-proj4string(psf1)

#let's take a look:
tm_shape(th)+tm_polygons()+tm_shape(psf1)+tm_dots(col = "Value",palette="-RdBu",size=0.5)+
tm_shape(chbound)+tm_polygons(alpha=0.3)

#working with spatial data is hard! need to figure out how to join Z data to a spatialdata thing
th.z     <- over(th,psf1) #makes the thiesen polygons have the Z value.
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.clp   <- raster::intersect(chbound,th.spdf)

TP_out<-tm_shape(th.clp)+tm_polygons(col="logValue",palette="-Spectral",alpha=0.5)+
  tm_shape(psf1)+tm_dots(col = "logValue",palette="-Spectral",size=0.1)
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
tmap_save(tm=TP_out,filename = "Ouputfile.jpg",dpi=300)

#calculate mass using the area of the polygon

th.clp$SqFt<-area(th.clp)
th.clp$mass_kg<-th.clp$Value*th.clp$SqFt*0.25*.092903*25*1000/1000/1000/1000
th.clp$mass_kg<-signif(th.clp$mass_kg,digits = 2)
sum(th.clp$mass_kg)

TP_out<-tm_shape(th.clp)+
  tm_polygons(col="mass_kg",palette="-Spectral",alpha=0.5,title = "Total Mass (kg)")+
  tm_text("mass_kg",size = 0.7)+
  tm_shape(psf1)+tm_dots(col = "Value",palette="-Spectral",size=0.1,title = "Max sample (µg/L)")
TP_out
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
tmap_save(tm=TP_out,filename = "Ouputfile.jpg",dpi=300)




