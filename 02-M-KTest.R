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

##----- Mann Kendall - Single well that I select
library(Kendall)
KendallData<-Dataset1
KendallData<-KendallData %>% subset(Well =="MW-1") %>% subset(Chem=="PCE")
outputk<-MannKendall(KendallData$Value)
Kendalloutput<-unlist(outputk, recursive = TRUE, use.names = TRUE)
InitialDF<-data.frame(Kendalloutput)
CountDF<-data.frame(Kendalloutput)
KendallDF<-rbind(t(InitialDF),t(CountDF))
row.names(KendallDF)<-NULL
KendallDF


# ----- Mann Kendall on a wells with at least 10 samples, and with 
KendallData<-Dataset1
library(Kendall)
attach(KendallData)
KendallData<- KendallData %>% unite("wellchem",Well,Chem,sep = "_")
KendallData$X<-NULL
KendallData$Y<-NULL
KendallData
#set up the data frame that will accept the MK data. 
#getting an error in MK, need to scrub out the data that is less than 3. S, unique count, remove the 
j<-length(unique(KendallData$wellchem))
mkout<-data.frame(name=rep("nam",j),
                  Samp_Num=rep(0,j),ND_num=rep(0,j),
                  tau=rep(0.1,j), sl = rep(0.1,j), 
                  S = rep(0.1,j), D = rep(0.1,j), 
                  varS =rep(0.1,j), COV=rep(0.1,j),
                  CF=rep(0.1,j),Descrip=rep("character",j))
mkout$name<-as.character(mkout$name)
mkout$Descrip<-as.character(mkout$Descrip)
mkout
k=1
#iterates telling you M-K statistics and the number of samples involved. 
for(i in unique(KendallData$wellchem)){
 # i<-c("MW-1_PCE")
  submk<-subset(KendallData,wellchem == i)
  sampnum<-length(submk$Value)
  NDnum<-sum(submk$Flag == "ND")
  
  if (sampnum>5){
    kout<-MannKendall(submk$Value)
    kout<-unlist(kout,recursive = TRUE, use.names = TRUE)
  } else {
    kout<-c(0,0,0,0,0)
  }
  mkout[k,4:8]<-kout #MK statistics
  mkout[k,3]<-NDnum 
  mkout[k,2]<-sampnum #number of total samples
  mkout[k,1]<-i #name
  covar<-mkout[k,5]
  covar<-as.numeric(covar)
  covar<-(1-(covar/2))
  if (sampnum>5){
    mkout[k,10]<-covar
    mkout[k,9]<- sd(submk$Value)/mean(submk$Value) #cov for MK trend determination
    } else{
    mkout[k,9:10]<-c(0,0)}#CF for MK trend determination
  S<-as.numeric(mkout[k,6])
  Cov1<-as.numeric(mkout[k,9])
  covar<-(mkout[k,10]) #THIS Line F-s it up? why?
  if ((S<0)&(covar>0.95)) {res1<-c("Decreasing")
  }else if ((S<0)&(covar>0.90)) {res1<-c("Probably Decreasing")
  }else if ((S>0)&(covar<=0.90)) {res1<-c("No Trend")
  }else if ((S>0)&(covar>0.95)){res1<-c("Increasing")
  }else if ((S>0)&(covar>0.90)){res1<-c("Probably Increasing")
  }else if ((S<=0)&((covar<0.9)&(Cov1>=1))) {res1<-c("No Trend")
  }else if ((S<=0)&((covar<0.9)&(Cov1<1))) {res1<-c("Stable")
  }else {res1<-c("NA")}
  mkout[k,11]<-res1
  k=k+1
}
mkout<-mkout %>% arrange(name) %>% separate(name,c("Name","Chem"),sep = "_")
mkout
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
write.csv(mkout,file = "Mann-KendallSummary_UNVALIDATED.csv")

#___________________________________
## Subsection: make a plot of this!
DatasetXY<-Dataset1
DatasetXY$Date<-NULL
DatasetXY$`Dup-Prim` <-NULL
DatasetXY$Chem<-NULL
DatasetXY$Value<-NULL
DatasetXY$Flag<-NULL
colnames(DatasetXY)[1] <- "Name"
DatasetXY
mkoutXY<-left_join(mkout,DatasetXY,by = "Name")
level_order <- c('Increasing', 'Probably Increasing', 'Stable','Probably Decreasing','Decreasing','No Trend')
map_MK<-ggplot(mkoutXY)+geom_point(aes(X,Y),alpha=0.5)+
  aes(colour=factor(Descrip,level = level_order))+
  scale_colour_discrete(name = "Result")+
 # aes(shape=factor(Descrip,level = level_order))+
  facet_wrap(~Chem)+theme_bw()+xlab("X (ft)")+ylab("Y (ft)")+ 
  coord_fixed(ratio = 1)+
  ggtitle("Mann-Kendall Test Result")
map_MK
setwd("C:/Users/omiles/Documents/TestR/GenericRf/Output/")
ggsave(map_MK,file="MK_resultAll.png",width=30,height=30,units="cm")







#### Eventually get this wroking: export the plots with MK statistics
# xample of plotting Mann-kendall
for(i in unique(SYS_LOC_CODE)){
  #  i<- "TW-12"
  subtoplot<-subset(All_data1, SYS_LOC_CODE == i)
  printplot<- ggplot(subtoplot)+
    geom_line(aes(SAMPLE_DATE,REPORT_RESULT_VALUE))+
    geom_point(aes(SAMPLE_DATE,REPORT_RESULT_VALUE,shape=DETECT_FLAG))+
    aes(color=CHEMICAL_NAME)+
    scale_y_log10(breaks =c(0.1,1,10,100,1000,10000,100000),minor_breaks=NULL,labels=c("0.1","1","10","100","1,000","10,000","100,000"),limits=c(0.1,100000))+
    scale_x_datetime(limits=xlims)+
    ylab("Concentration (µg/L)") + 
    xlab("Date")+ggtitle(i)+labs(color="Compound")+ theme_minimal()
  printplot
  ggsave(printplot,file=paste0(i,".png"),width=20,height=14,units="cm")
}


