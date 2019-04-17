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
MMDatabase <- read_excel("~/TestR/GenericRf/MolarMassDatabase.xlsx")
Dataset1<-Datasource #to keep this clean I'll assign it to something "working"
MDataset1<-left_join(Dataset1,MMDatabase,by = "Chem")
MDataset1$ValueM<-MDataset1$Value/MDataset1$MolarMass
#adds µmol/L column
MDataset1$Value<-MDataset1$ValueM
MDataset1$MolarMass<-NULL
MDataset1$ValueM<-NULL

#________________________________________________________________________
#________06.1_Violin plots of data by year (max of year) __________________

#first part just "bins" the data. using for other plots like this, a little redundant.
Dataset1bin<-MDataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>% summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)
Dataset1bin

library(gganimate)
#note LEvel order is key for other graphics
level_order <- c('PCE', 'TCE', 'cDCE',"VC")
Animate1<-ggplot(Dataset1bin)+aes(x=factor(Chem, level = level_order),Value)+aes(color=Chem)+
  geom_violin(scale="width")+geom_jitter(width=0.01,height = 0.01)+
  scale_y_log10()+transition_time(yr)+ease_aes('cubic-in-out')+
  labs(title = 'Date: {frame_time}', x = 'Chemical', y = 'Concentration (µmol/L)')
animate(Animate1, nframes = 350, fps = 13, height = 450, width =750)
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
anim_save("ViolinPlotsMolar.gif")

#________08_RadialPlots  of data in awell by year ___________________
#Note this doesnt work really well right now...
Dataset1bin<-MDataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>%  summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)

Dataset1bin<-subset(Dataset1bin,Well=="MW-32")
Dataset1bin<-subset(Dataset1bin,AquiferUnit=="Unit2")
Radial1<-ggplot(Dataset1bin,aes(Chem,Value))+geom_bar(stat = "identity",aes(fill=Chem))+
  scale_y_log10()+coord_polar()+
  labs(title = 'MW-31', x = 'Chemical', y = 'Concentration (µmol/L')+
  facet_wrap(~ yr)
Radial1
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(Radial1,file=paste0("RadialPlot_onewellMolar.png"),width=25,height=25,units="cm")



#_______010 converting to mol from ug and generating barplots

Dataset1bin<-MDataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>%  summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)
Dataset1bin<-subset(Dataset1bin,Well=="MW-32")
Dataset1bin<-subset(Dataset1bin,AquiferUnit=="Unit2")
Dataset1bin

level_order <- c('PCE', 'TCE', 'cDCE',"VC")
barplots<-ggplot(Dataset1bin,aes(yr,Value))+geom_bar(stat = "identity",position="fill")+
  aes(color=factor(Chem, level = level_order))+labs(color="Compound")+
  aes(fill=factor(Chem, level = level_order))+labs(fill="Compound")+
  labs(title = paste("Well: MW-32"), x = 'year', y = 'Molar Ratio')+
  theme_minimal()
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(barplots,file=paste0("Barplot_aWell.png"),width=25,height=10,units="cm")
barplots

#________________________________________________________________________
#________06.1_Violin plots of data by year (max of year) __________________

#first part just "bins" the data. using for other plots like this, a little redundant.
Dataset1bin<-MDataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>% summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)
Dataset1bin

#note LEvel order is key for other graphics
level_order <- c('PCE', 'TCE', 'cDCE',"VC")
Dataset1bin<-subset(Dataset1bin,yr=="2010")
viol1<-ggplot(Dataset1bin)+aes(x=factor(Chem, level = level_order),Value)+aes(color=Chem)+
  geom_violin(scale="width")+geom_jitter(width=0.01,height = 0.01)+
  scale_y_log10()+
  labs(title = 'Date: 2010', x = 'Chemical', y = 'Concentration (µmol/L)')
ggsave(viol1,file=paste0("viol1.png"),width=25,height=10,units="cm")


Dataset2bin<-Dataset1
Dataset2bin$bin <- cut(Dataset2bin$Date, breaks = "years")
Dataset2bin<-Dataset2bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>% summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)
Dataset2bin

#note LEvel order is key for other graphics
level_order <- c('PCE', 'TCE', 'cDCE',"VC")
Dataset2bin<-subset(Dataset2bin,yr=="2010")
viol2<-ggplot(Dataset2bin)+aes(x=factor(Chem, level = level_order),Value)+aes(color=Chem)+
  geom_violin(scale="width",alpha=0.2)+geom_jitter(width=0.01,height = 0.01)+
  scale_y_log10()+
  labs(title = 'Date: 2010', x = 'Chemical', y = 'Concentration (µg/L)')
viol2
ggsave(viol2,file=paste0("viol2.png"),width=25,height=10,units="cm")


