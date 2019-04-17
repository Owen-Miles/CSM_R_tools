#_____________________________________________________________________
#_____4_1 Thiessen Polygons and other interpolations _________________

#the Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(sp)
library(sf)
library(rgdal)
library(tmap)
library(ggalt)
library(spatstat)  # Used for the dirichlet tessellation function
library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons


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

###OK THIS IS KEY!. now i have a real spatial set
####2872 is crs for California stplane in feet
p.sf<-sf::st_as_sf(Maxvlaues, coords = c("X", "Y"), crs = 2872)
psf1<-as(p.sf,"Spatial")#this is also very key. to use the examples i need stuff in spatial form, not as sf form. 


#take a look a the data:
#tm_shape(p.sf)+tm_dots(col = "Value",palette="RdBu",size=0.7)

#need to make a ppp thing. this is a "point pattern" dataset, does not have spatial yet
windowsX<-c(min(Maxvlaues$X)-1000,max(Maxvlaues$X)+1000)#a  "window" is key
WindowsY<-c(min(Maxvlaues$Y)-1000,max(Maxvlaues$Y)+1000)
pts.ppp<-ppp(Maxvlaues$X,Maxvlaues$Y,window=owin(windowsX,WindowsY))
th<-as(dirichlet(pts.ppp),"SpatialPolygons") #generate the polygons
proj4string(th)<-proj4string(psf1) #and give it a spatial designation

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
tm_shape(th)+tm_polygons()+tm_shape(psf1)+tm_dots(size=0.5)+
  tm_shape(chbound)+tm_polygons(alpha=0.1)

#working with spatial data is hard! need to figure out how to join Z data to a spatialdata thing
th.z     <- over(th,psf1) #makes the thiesen polygons have the Z value.
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.clp   <- raster::intersect(chbound,th.spdf)


#calculate mass using the area of the polygon

th.clp$SqFt<-area(th.clp)
th.clp$mass_kg<-th.clp$Value*th.clp$SqFt*0.25*.092903*25*1000/1000/1000/1000
th.clp$mass_kg<-signif(th.clp$mass_kg,digits = 2)
MassTP<-signif(sum(th.clp$mass_kg),2)

TP_out<-tm_shape(th.clp)+
  tm_polygons(col="logValue",palette="-Spectral",title = "Concentration (log) (µg/L)")+
  #tm_text("mass_kg",size = 0.7)+
  tm_shape(psf1)+tm_dots(size=0.1)+
  #tm_grid(alpha=0.2)+
  tm_layout(title = paste("TP=",MassTP,"kg"))+tm_legend(legend.outside=TRUE)
TP_out



#____________________IDW________________________
#dependent, for right now, on above.

library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function

grd <- as.data.frame(spsample(chbound, "regular", n=50000))

names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(chbound)
#note IDP is power value, maybe should be 10? what is value after tilda
P.idw <- gstat::idw(Value ~ 1, psf1, newdata=grd, idp=2.0)
# Convert to raster object then clip
r       <- raster(P.idw)
r.m     <- mask(r, chbound)
MassIDW<-signif(cellStats(r.m,stat = "mean")*area(chbound)*0.25*.092903*25*1000/1000/1000/1000,2)

P.idwl <- gstat::idw(logValue ~ 1, psf1, newdata=grd, idp=2.0)
# Convert to raster object then clip
rl       <- raster(P.idwl)
r.ml     <- mask(rl, chbound)

#linearIDW<-tm_shape(r.m) + 
 # tm_raster(palette = "-Spectral",  title="Concentration (linear) (µg/L)") + 
  #tm_shape(psf1) + tm_dots(size=0.2) +tm_legend(legend.outside=TRUE)+
  #tm_layout(title = bquote(paste(MassIDW,"kg")))
LogIDW<-tm_shape(r.ml) + 
  tm_raster(palette = "-Spectral",  title="Concentration (log) (µg/L)") + 
  tm_shape(psf1) + tm_dots(size=0.2) +tm_legend(legend.outside=TRUE)+
  tm_layout(title = paste("IDW=",MassIDW,"kg"))
LogIDW
?tm_legend

###_____________Kriging estimate______________
# Define the 1st order polynomial equation
f.1 <- as.formula(logValue ~ X + Y)
f.1

coordinates(psf1)
psf2<-psf1
df<-coordinates(psf1)
df[,2]
psf2$X<-df[,1]
psf2$Y<-df[,2]
# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
#var.smpl <- variogram(f.1, psf2, cloud = FALSE, cutoff=1000000, width=89900)
var.smpl <- variogram(f.1, psf2)
# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=14, model="Sph", nugget=0))
dat.fit  <- fit.variogram(var.smpl, vgm(psill=14, model="Sph", nugget=0))

# The following plot allows us to assess the fit
plot(var.smpl, dat.fit, xlim=c(0,10000))

# created in the earlier step)
dat.krg <- krige( f.1, psf2, grd, dat.fit)
r <- raster(dat.krg)
r.m <- mask(r, chbound)

krig_Mass<-signif(mean(10^values(r.m),na.rm=TRUE)*area(chbound)*0.25*.092903*25*1000/1000/1000/1000,2)
krig_Mass<-as.numeric(krig_Mass)
Krig_out<-tm_shape(r.m) + 
  tm_raster(palette = "-Spectral",title="Concentration (log) (µg/L)")+
  tm_shape(psf1) + tm_dots(size=0.2)+
  tm_legend(legend.outside=TRUE) +tm_layout(title = paste("Krig =",krig_Mass,"kg"))
Krig_out

Basemap<-tm_shape(chbound)+tm_polygons()+
  tm_shape(psf1)+tm_dots(col="logValue",size=0.7,palette="-Spectral",title = "Concentration (log) (µg/L)")+
  tm_legend(legend.outside=TRUE) +tm_layout(title = paste("Maximum TCE sample"))
Basemap
Compare<-tmap_arrange(Basemap,LogIDW, TP_out,Krig_out,ncol = 2,asp = 0.750,outer.margins = 0.01)
Compare
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
tmap_save(tm=Compare,filename = "Ouputfile.jpg",width = 20,height=20,units = "cm")



