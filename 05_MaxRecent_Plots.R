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

#_____________________________________________________________________
#____________01_Max value in each well Table of max________________

# saves an excel file with the max value for each well and chemical 
Maxvlaues<-Dataset1 %>% 
  unite(Well_Unit_Chem,Well,AquiferUnit,Chem,sep = "_") %>% 
  group_by(Well_Unit_Chem) %>% 
  filter(Value == max(Value)) %>% 
  filter(Date == max(Date)) %>% 
  separate(Well_Unit_Chem,c("Well","AquiferUnit","Chem"),sep = "_") %>% 
  arrange(Well)  
Maxvlaues$MaxValue<-Maxvlaues$Value
Maxvlaues$Value<-NULL
Maxvlaues<-Maxvlaues %>% unite(key1,Well,Chem,sep = "_")

Mostrecent<-Dataset1
Mostrecent$'Dup-Prim'<-NULL
Mostrecent
Mostrecent<-Mostrecent %>% 
  unite(Well_Unit_Chem,Well,AquiferUnit,Chem,sep = "_") %>% 
  group_by(Well_Unit_Chem) %>% 
  filter(Date == max(Date)) %>% 
  filter(Value == max(Value)) %>% 
  separate(Well_Unit_Chem,c("Well","AquiferUnit","Chem"),sep = "_") %>% 
  arrange(Well) %>% unique() 
Mostrecent$RecentValue<-Mostrecent$Value
Mostrecent$Value<-NULL
Mostrecent<-Mostrecent %>% unite(key1,Well,Chem,sep = "_")
joined<-left_join(Mostrecent,Maxvlaues,by="key1")
joined$AquiferUnit.x<-NULL
joined$Date.x<-NULL
joined$X.x<-NULL
joined$Y.x<-NULL
joined$Value<-NULL
joined$Flag.x<-NULL
joined$AquiferUnit.y<-NULL
joined$Date.y<-NULL
joined$`Dup-Prim`<-NULL
joined$X.y<-NULL
joined$Y.y<-NULL
joined$Flag.y<-NULL
joined<-joined %>% unique() %>% separate(key1,c("Well","Chem"),sep="_")
joined
plotjoin<-joined
#plotjoin<-subset(joined,Chem=="TCE")
plotjoin$Reduced<-plotjoin$MaxValue/plotjoin$RecentValue
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
library(ggrepel)

for(i in unique(plotjoin$Chem)){
  subtoplot<-subset(plotjoin, Chem == i)
  p<-ggplot(subtoplot,aes(MaxValue,RecentValue,label=Well))+geom_point()+
    scale_x_log10(breaks =c(0.1,1,10,100,1000,10000,100000),minor_breaks=NULL,labels=c("0.1","1","10","100","1,000","10,000","100,000"),limits=c(0.1,100000))+
    scale_y_log10(breaks =c(0.1,1,10,100,1000,10000,100000),minor_breaks=NULL,labels=c("0.1","1","10","100","1,000","10,000","100,000"),limits=c(0.1,100000))+
    aes(color=Reduced)+
    scale_color_gradient(low = "red", high = "green4",trans="log10")+
    theme_minimal()+
    #geom_label_repel(aes(label=Well))+
    geom_text_repel(data = subset(subtoplot,MaxValue>1000),size = 3,box.padding = 0.7)+
    coord_fixed(ratio = 1)+labs(color="Reduction")+
    ylab("Most recent sample (µg/L)") + 
    xlab("Historical max sample (µg/L)")+ggtitle(i)+
  ggsave(p,file=paste0(i,"_HistMax_Current.png"),width=20,height=14,units="cm")
}
  

