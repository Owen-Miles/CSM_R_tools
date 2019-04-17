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
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
library(xlsx)
write.xlsx(Maxvlaues, "MaxChemSample.xlsx")
Maxvlaues
#_____________________________________________________________________
#____________02_Max value in each well maps of max by compound________
Maxvlaues<-Dataset1 %>% 
  unite(Well_Unit_Chem,Well,AquiferUnit,Chem,sep = "_") %>% 
  group_by(Well_Unit_Chem) %>% 
  filter(Value == max(Value)) %>% 
  filter(Date == max(Date)) %>% 
  separate(Well_Unit_Chem,c("Well","AquiferUnit","Chem"),sep = "_") %>% 
  arrange(Well) 
map_max<-ggplot(Maxvlaues)+geom_point(aes(X,Y),alpha=0.3)+
  aes(size=Value)+scale_radius(range=c(0.5,10),trans="log10")+aes(colour=Value)+
  scale_color_gradient(low = "green4", high = "red",trans="log10")+
  facet_wrap(~Chem)+theme_bw()+xlab("X (ft)")+ylab("Y (ft)")+
  coord_fixed(ratio = 1)+
  ggtitle("Spatial Distribution: Max of each well")
map_max
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(map_max,file="Map_of_Maxs.png",width=30,height=30,units="cm")

#________________________________________________________________________
#________03_Total VOCs, max of each well, mapped, by aquifer unit________

###This section just calculates a max VOCs y spreading data (could be useful elsewhere)
Tvocs<-Dataset1
Tvocs$`Dup-Prim`<-NULL
Tvocs$Flag<-NULL
Tvocs<-Tvocs %>% unite(key1,Well,AquiferUnit,Date,X,Y,sep = "_") %>% 
  arrange(key1) %>% unite(key2,key1,Chem, sep = "_o_") %>% 
  group_by(key2) %>% summarise(Value=max(Value)) %>% 
  separate(key2,c("key1","Chem"),sep = "_o_") %>%
  spread(Chem,Value) %>% 
  separate(key1,c("Well","AquiferUnit","Date","X","Y"),convert=TRUE,sep = "_")
Tvocs$TotalVocs<-Tvocs$cDCE+Tvocs$PCE+Tvocs$TCE+Tvocs$VC

#this section makes the max value
Maxvlaues<-Tvocs %>% unite(key1,Well,AquiferUnit,X,Y,sep = "_") %>% 
  group_by(key1) %>%  summarize(TotalVocs = max(TotalVocs)) %>% 
  separate(key1,c("Well","AquiferUnit","X","Y"),convert=TRUE,sep = "_")
#this was a thought, perhaps plot the wells in the background, kind of hard to get working
Wellunique<-Maxvlaues
Wellunique$AquiferUnit<-NULL
Wellunique$TotalVocs<-NULL

Maxvlaues
map_max<-ggplot(Maxvlaues)+
  # this line adds a reference point of all wlels, may be confusing, comment out?
  geom_point(data = Wellunique,aes(X,Y), size=0.1,alpha=0.1)+
  geom_point(aes(X,Y,size=TotalVocs,colour=TotalVocs),alpha=0.4)+
  scale_radius(range=c(0.75,10),trans="log10")+
  scale_color_gradient(low = "green4", high = "red",trans="log10")+
  facet_wrap(~AquiferUnit)+theme_bw()+xlab("X (ft)")+ylab("Y (ft)")+
  coord_fixed(ratio = 1)+
  ggtitle("Spatial Distribution: Max of each well")
map_max
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(map_max,file="Map_of_Maxs_AquiferUnit.png",width=40,height=40,units="cm",dpi=400)



#________________________________________________________________________
#________04_Histogram of all samples_____________________________________

#adjust the secbinwidth to change the number of days
secbinwidth <- 1*30*24*60*60  #note, binsize 30 days in seconds, because R thinks in seconds,
histgram <- ggplot(Dataset1) + geom_histogram(aes(Date),binwidth = secbinwidth)+
  ggtitle(paste("Histogram of Samples: bin width", secbinwidth/60/60/24,"days",sep=" "))+
  theme_minimal()+ylab("Count")
histgram
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(histgram,file="Histogram_of_Samples.png",width=16,height=14,units="cm")

#________________________________________________________________________
#________05_Loop that generates log-TC plots of all compounds____________
Dataset2<-Dataset1
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/TC-Plots")
#the limits were determined by the histogram, visually. can also make that automatic
xlims<-as.POSIXct(c("1995-01-01","2019-01-01"))
level_order <- c('PCE', 'TCE', 'cDCE',"VC")
for(i in unique(Dataset2$Well)){
  subtoplot<-subset(Dataset2, Well == i)
  #sampnum<-length(subtoplot$Value)  #if we want to only do wells with some min number of samples
  printplot<- ggplot(subtoplot)+
    geom_line(aes(Date,Value))+
    #for loess smooth lines:
    # geom_smooth(aes(Date,Value),weight = 0.2, alpha = 0.2)+
    geom_point(aes(Date,Value,shape=Flag))+
    aes(color=factor(Chem, level = level_order))+
    scale_y_log10(breaks =c(0.1,1,10,100,1000,10000,100000),minor_breaks=NULL,labels=c("0.1","1","10","100","1,000","10,000","100,000"),limits=c(0.1,100000))+
    scale_x_datetime(limits=xlims)+
    ylab("Concentration (µg/L)") + 
    xlab("Date")+ggtitle(paste("Well:",i,"Aquifer unit:",subtoplot$AquiferUnit))+labs(color="Compound")+ theme_minimal()
  printplot
  ggsave(printplot,file=paste0(i,".png"),width=20,height=14,units="cm")
}


#________________________________________________________________________
#________06_Violin plots of data by year (max of year) __________________

#first part just "bins" the data. using for other plots like this, a little redundant.
Dataset1bin<-Dataset1
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
  labs(title = 'Date: {frame_time}', x = 'Chemical', y = 'Concentration (µg/L')
animate(Animate1, nframes = 350, fps = 13, height = 450, width =750)
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
anim_save("ViolinPlots.gif")



#________________________________________________________________________
#________07_RadialPlots  of data by year (max of year) __________________
#first part just "bins" the data. using for other plots like this, a little redundant.
library(gganimate)
Dataset1bin<-Dataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>%  summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)

#need some sort of filter for this. otherwise we get too many plots.
#for now, say there are four wells i'm particularly interested in:
radialplotlist<-c("MW-21","MW-20","MW-23")
                  #"MW-6","MW-20", "MW-38","MW-88")
# could select by number of them?
# - wide range in concentration (at least 3 orders of magnitude)
# - could select the first few and save them, iterate, whatever

##
# SOMETHING IS WRONG?! it's not plotting correctly. 
# Coercing certain wells for some reason?

Dataset1bin<-subset(Dataset1bin,Well %in% radialplotlist)
#View(Dataset1bin)
Animate1<-ggplot(Dataset1bin,aes(Chem,Value))+geom_bar(stat = "identity",aes(fill=Chem))+
  scale_y_log10()+coord_polar()+
  transition_time(yr)+ease_aes("cubic-in-out")+
  labs(title = 'Date: {frame_time}', x = 'Chemical', y = 'Concentration (µg/L')+
  facet_wrap(~ Well)
animate(Animate1, nframes = 250, fps = 12, height = 400, width =800, start_pause = 5,end_pause = 5)
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
anim_save("RadialPlots.gif")

##
#________08_RadialPlots  of data in a year __________________
#first part just "bins" the data. using for other plots like this, a little redundant.
Dataset1bin<-Dataset1
Dataset1bin$bin <- cut(Dataset1bin$Date, breaks = "years")
Dataset1bin<-Dataset1bin %>% separate(bin,c("yr","m","d")) %>% 
  unite(Key1,Well,AquiferUnit,Chem,yr,X,Y,sep = "_") %>% 
  group_by(Key1) %>%  summarize(Value=max(Value)) %>% 
  separate(Key1,c("Well","AquiferUnit","Chem","yr","X","Y"),sep = "_",convert = TRUE) %>% 
  arrange(Well)

radialplotlist<-c("MW-21","MW-20","MW-23","MW-6","MW-20", "MW-38","MW-88",
                  "MW-22","MW-23","MW-28","MW-25","MW-26","MW-31","MW-32",
                  "MW-33","MW-34","MW-35","MW-37")

Dataset1bin<-subset(Dataset1bin,Well %in% radialplotlist)
Dataset1bin<-subset(Dataset1bin,yr=2010)
Radial1<-ggplot(Dataset1bin,aes(Chem,Value))+geom_bar(stat = "identity",aes(fill=Chem))+
  scale_y_log10()+coord_polar()+
  labs(title = 'Date: 2010', x = 'Chemical', y = 'Concentration (µg/L)')+
  facet_wrap(~ Well)
Radial1
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(Radial1,file=paste0("RadialPlot_aYear.png"),width=25,height=25,units="cm")


#________08_RadialPlots  of data in awell by year __________________
#first part just "bins" the data. using for other plots like this, a little redundant.
Dataset1bin<-Dataset1
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
  labs(title = 'MW-31', x = 'Chemical', y = 'Concentration (µg/L)')+
  facet_wrap(~ yr)
Radial1
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(Radial1,file=paste0("RadialPlot_onewell.png"),width=25,height=25,units="cm")


