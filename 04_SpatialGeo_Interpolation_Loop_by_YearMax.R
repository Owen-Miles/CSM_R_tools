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
library(gstat) # Use gstat's idw routine
library(mapview)
 # Used for the spsample function

Datasource <- read_excel("~/TestR/GenericRf/Datasource.xlsx", 
                         col_types = c("text", "text", "date", 
                                       "text", "numeric", "numeric", "text", 
                                       "numeric", "text"))
Dataset1<-Datasource #to keep this clean I'll assign it to something "working"
Dataset1<-subset(Datasource,`Dup-Prim`=="Primary") # having issues removing dups, so I'm doing it first thing

Dataset1bin<-Dataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>% summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)

Dataset1bin<-Dataset1bin %>% subset(AquiferUnit=="UnitA")
#Dataset1bin<-Dataset1bin %>% subset(yr>2005)
#Dataset1bin<-Dataset1bin %>% subset(Chem == "TCE")

Maxvlaues1<-Dataset1bin %>%unite("chem_yr",Chem,yr,sep = "_") 
Maxvlaues1<-Maxvlaues1[order(Maxvlaues1$chem_yr),]
Maxvlaues1
d= NULL 

for(i in unique(Maxvlaues1$chem_yr)){
  Maxvlaues<-subset(Maxvlaues1, chem_yr == i)  
  
Maxvlaues$logValue<-log10(Maxvlaues$Value)
p.sf<-sf::st_as_sf(Maxvlaues, coords = c("X", "Y"), crs = 2872)
psf1<-as(p.sf,"Spatial")#this is also very key. to use the examples i need stuff in spatial form, not as sf form. 
windowsX<-c(min(Maxvlaues$X)-1000,max(Maxvlaues$X)+1000)#a  "window" is key
WindowsY<-c(min(Maxvlaues$Y)-1000,max(Maxvlaues$Y)+1000)
pts.ppp<-ppp(Maxvlaues$X,Maxvlaues$Y,window=owin(windowsX,WindowsY))
th<-as(dirichlet(pts.ppp),"SpatialPolygons") #generate the polygons
proj4string(th)<-proj4string(psf1) #and give it a spatial designation

MV<-Maxvlaues
MV[,1:3]<-NULL
MV[,3:4]<-NULL
MV
Convexhull1<-chull(MV)
coordTable<-as.matrix(MV[Convexhull1,])

ch_bound<-coords2Polygons(coordTable,ID="A")
proj4string(ch_bound)<-proj4string(psf1)
chbound<-buffer(ch_bound, width = 500)
proj4string(chbound)<-proj4string(psf1)

tm_shape(th)+tm_polygons()+tm_shape(psf1)+tm_dots(size=0.5)+
  tm_shape(chbound)+tm_polygons(alpha=0.1)

th.z     <- over(th,psf1) #makes the thiesen polygons have the Z value.
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.clp   <- raster::intersect(chbound,th.spdf)


th.clp$SqFt<-area(th.clp)
th.clp$mass_kg<-th.clp$Value*th.clp$SqFt*0.25*.092903*25*1000/1000/1000/1000
th.clp$mass_kg<-signif(th.clp$mass_kg,digits = 2)
MassTP<-signif(sum(th.clp$mass_kg),2)

TP_out<-tm_shape(th.clp)+
  tm_polygons(col="logValue",palette="-Spectral",title = "Concentration (log) (µg/L)")+
  tm_shape(psf1)+tm_dots(size=0.1)+
  tm_layout(title = paste("TP=",MassTP,"kg"))+tm_legend(legend.outside=TRUE)
TP_out
#____________________IDW________________________
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
LogIDW<-tm_shape(r.ml) + 
  tm_raster(palette = "-Spectral",  title="Concentration (log) (µg/L)") + 
  tm_shape(psf1) + tm_dots(size=0.2) +tm_legend(legend.outside=TRUE)+
  tm_layout(title = paste("IDW=",MassIDW,"kg"))

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
var.smpl <- variogram(f.1, psf2)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=14, model="Sph", nugget=0))
dat.fit  <- fit.variogram(var.smpl, vgm(psill=14, model="Sph", nugget=0))

plot(var.smpl, dat.fit, xlim=c(0,10000))

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
  tm_legend(legend.outside=TRUE) +tm_layout(title = paste(i,"Sample Value"))
Compare<-tmap_arrange(Basemap,LogIDW, TP_out,Krig_out,ncol = 2,asp = 0.750,outer.margins = 0.01)
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/MassEstimates/")
tmap_save(tm=Compare,filename = paste0(i,"_MassMap_year_.jpg"),width = 20,height=20,units = "cm")
d=rbind(d, data.frame(i,MassIDW,MassTP,krig_Mass))
}
colnames(d)[colnames(d)=="MassIDW"] <- "IDW"
colnames(d)[colnames(d)=="MassTP"] <- "ThsnPoly"
colnames(d)[colnames(d)=="krig_Mass"] <- "Kriging"
d1<-d
d1<-separate(d1,i,c("Chem","yr"),sep = "_")
library(xlsx)
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/MassEstimates/")
write.xlsx(d1, "MassEstimates.xlsx")
d2<-unite(d1,"chem_yr",Chem,yr,sep = "_")
d3<-gather(d2,"Interpolation","Mass_Est",2:4)
d3<-separate(d3,chem_yr,c("Chem","yr"))
d3$yr<-as.numeric(d3$yr)
massestplot<-ggplot(d3,aes(yr,Mass_Est))+geom_point()+geom_line()+aes(colour=Interpolation)+
  scale_y_log10()+facet_wrap(~Chem)+
  theme_minimal()+ylab("Mass Kg") +xlab("Date")
ggsave(massestplot,file="MassEstimateSummaryPlot.png",width=20,height=14,units="cm")

massestplot<-ggplot(d3,aes(yr,Mass_Est))+geom_point()+geom_line()+aes(colour=Interpolation)+
  scale_y_log10()+facet_wrap(~Chem)+
  theme_minimal()+ylab("Mass Kg") +xlab("Date")+coord_polar()
massestplot

