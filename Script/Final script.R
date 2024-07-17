library(survival)
library(ggplot2)
library(DataCombine)
library(stringr)
library("survminer")
library(car)
library(dplyr)
library(emmeans)
library(glmmTMB)
library(readxl)
library(tidyverse)
library(lme4)
library(writexl)
library(plyr)
library(vegan)
library(pairwiseAdonis)
library(corrplot)
library(FactoMineR)
library(rstatix)
library(hms)
library(fmsb)
library(mvabund)
library(reshape2)


cytcolors<-c(Linareva="#BC4B07",Mahana="#D3BC38",Manava="#5BA7A3",Nursery="#424242")
sem <- function(x) sd(x)/sqrt(length(x))

##Temperature graph
ManaTemp<-read.csv("Data/manahobo1.csv",header=F,sep=",")
ManaTemp$timestamp<-as.POSIXct(ManaTemp$V2,format="%m/%d/%Y %I:%M:%OS %p")
ManaTemp$Manava<-ManaTemp$V3
MahaTemp<-read.csv("Data/mahahobo1.csv",header=F,sep=",")
MahaTemp$timestamp<-as.POSIXct(MahaTemp$V2,format="%m/%d/%Y %I:%M:%OS %p")
MahaTemp$Mahana<-MahaTemp$V3
LinaTemp<-read.csv("Data/linhobo4.csv",header=F,sep=",")
LinaTemp<-LinaTemp[,1:5]
LinaTemp$timestamp<-as.POSIXct(LinaTemp$V2,format="%m/%d/%Y %I:%M:%OS %p")
LinaTemp$Linareva<-LinaTemp$V3
NursTemp<-read.csv("Data/nurseryhobo.csv",header=F,sep=",")
NursTemp$timestamp<-as.POSIXct(NursTemp$V2,format="%m/%d/%Y %I:%M:%OS %p")
NursTemp$Nursery<-NursTemp$V3

total<-inner_join(ManaTemp[,6:7],MahaTemp[,6:7],by="timestamp")
total<-inner_join(total,LinaTemp[,6:7],by="timestamp")
total<-inner_join(total,NursTemp[,6:7],by="timestamp")

total<-total %>% mutate(timestamp = str_replace_all(timestamp, "0021", "2021"))
total$timestamp<-as.POSIXct(total$timestamp)
temp<-ggplot(total,aes(x=timestamp))+ylab("Water temperature (°C)")+xlab("")+
  geom_line(aes(y = Linareva), color = "#BC4B07")+
  geom_line(aes(y = Mahana), color = "#D3BC38")+
  geom_line(aes(y = Manava), color = "#5BA7A3")+
  geom_line(aes(y = Nursery), color = "#424242")+
  scale_x_datetime(date_labels = "%Y-%m",date_breaks="30 days")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
temp
ggsave("Output/Temp.pdf",temp,width = 9, height = 4)
ggsave("Output/Temp.svg",temp,width = 9, height = 4)

##Temperature stats
library(tseries)
adf.test(total$Linareva)
adf.test(total$Mahana)
adf.test(total$Manava)
adf.test(total$Nursery)
##all show constant variance over the 3 measured months

M=cor(total[,2:5],use="pairwise.complete.obs")
res1 <- cor.mtest(total[,2:5], conf.level = 0.95)
pvalues<-as.data.frame(res1$p)
corrplot(M, method = 'color',tl.col = 'black',p.mat = res1$p,sig.level = 0.05,insig = "blank")
cov(total[,2:5],use="pairwise.complete.obs")

total2<-gather(total,"Site","Temp",Manava,Mahana,Linareva,Nursery)
total2$Day<-as.Date(total2$timestamp,format="%Y-%m-%d")
totalmean<-ddply(total2,c("Day","Site"),summarise,Mean=mean(Temp),Max=max(Temp),Min=min(Temp))
totalmean$Range<-totalmean$Max-totalmean$Min
shapiro.test(totalmean$Mean)
bartlett.test(Mean ~ Site, data = totalmean)
anova_test(totalmean,Mean~Site)
shapiro.test(totalmean$Max)
bartlett.test(Max ~ Site, data = totalmean)
pairwise.t.test(totalmean$Max, totalmean$Site)
shapiro.test(totalmean$Min)
kruskal.test(totalmean$Min,totalmean$Site)
pairwise.wilcox.test(totalmean$Min,totalmean$Site)
shapiro.test(totalmean$Range)
kruskal.test(totalmean$Range,totalmean$Site)
pairwise.wilcox.test(totalmean$Range,totalmean$Site)
aggregate(Range~Site,totalmean,mean)
aggregate(Range~Site,totalmean,sem)

##long term temperature
chenal<-read_excel("Data/Chenal_2018.xlsx")
chenal$time<-paste(chenal$date,substr(chenal$Heure,12,20))                  
colnames(chenal)[2]<-'Temp_chenal'
chenal$time<-as.POSIXct(chenal$time,format="%Y-%m-%d %H:%M:%OS",tz="Pacific/Tahiti")
chenal$Temp_chenal<-as.numeric(chenal$Temp_chenal)
chenal20<-subset(chenal, grepl('^2020', chenal$date))
F1<-read_excel("Data/F1_2018-2021.xlsx")
F1$timestamp<-paste(F1$date,substr(F1$heure,12,20),F1$pm)                  
F1$time<-as.POSIXct(F1$timestamp,format="%Y-%m-%d %I:%M:%OS %p",tz="Pacific/Tahiti")
F1$Temp_F<-as.numeric(F1$Temp_F)
total<-left_join(F1,chenal,by="time")

temp<-ggplot(total,aes(x=time))+ylab("Water temperature (°C)")+xlab("")+
  geom_line(aes(y = Temp_F), color = "#5BD9D9")+
  geom_line(aes(y = Temp_chenal), color = "#f0aa7f")+
  scale_x_datetime(date_labels = "%Y-%m",date_breaks="1 month")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  geom_hline(yintercept=31,color="#bf4d58")
temp
ggsave("Output/longtemp.jpg",temp,width = 9, height = 4)


##Nutrients
nuts<-read_xlsx("Data/Nutrients.xlsx")
sand<-nuts %>% filter(Sample=="Sand")
water<-nuts %>% filter(Sample=="Water")

pca1<-PCA(sand[,5:9])
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
sand$pc1 <- pca1$ind$coord[, 1]
sand$pc2 <- pca1$ind$coord[, 2]
Pca1<-ggplot(data = sand, aes(x = pc1, y = pc2,color=Site)) +
  geom_point()+
  scale_color_manual(values=cytcolors)+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1*2, y = 0, 
                                                    yend = Dim.2*2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*2.5, y =  Dim.2*2.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (48.45% variance)",y="PCA2 (35.77% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="none")+
  ggtitle("Sand")
Pca1

pca2<-PCA(water[,5:9])
pca.vars <- pca2$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
water$pc1 <- pca2$ind$coord[, 1]
water$pc2 <- pca2$ind$coord[, 2]
Pca2<-ggplot(data = water, aes(x = pc1, y = pc2,color=Site)) +
  geom_point()+
  scale_color_manual(values=cytcolors)+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1*2, y = 0, 
                                                    yend = Dim.2*2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*2.5, y =  Dim.2*2.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (41.58% variance)",y="PCA2 (28.32% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right")+
  ggtitle("Water")
Pca2

nutrients<-ggarrange(Pca1,Pca2,widths=c(0.8,1))
nutrients
ggsave("Output/nutrients.jpg",nutrients,width = 7, height = 3)

kruskal.test(sand$'SI(OH)4',sand$Site)
kruskal.test(water$PO4,water$Site)
library(dunn.test)
dunn.test(sand$'SI(OH)4',sand$Site,method="holm")
dunn.test(water$'SI(OH)4',water$Site,method="holm")


##Algal nitrogen
dt1 <-read.csv("Data/mcr_lter_macroalgaeN.csv",header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Year",     
                 "Site",     
                 "Habitat",     
                 "Genus",     
                 "Dry_Weight",     
                 "C",     
                 "H",     
                 "N",     
                 "CN_ratio",     
                 "Comment"    ), check.names=TRUE)
dt2<-dt1 %>% filter(Genus=="Turbinaria" & Year==2020 & Habitat=="Back Reef")
dt3<-dt2 %>% filter(Site %in% c("LTER 1","LTER 2","LTER 5","LTER 6"))
replace<-as.data.frame(c("LTER 1","LTER 2","LTER 5","LTER 6"))
colnames(replace)[1]<-"Site"
replace$New<-c("Nursery","Manava","Linareva","Mahana")
dt3<-FindReplace(dt3,"Site",replace,"Site","New",exact=TRUE, vector=FALSE)
shapiro.test(dt3$N)
bartlett.test(N ~ Site, data = dt3)
anova_test(dt3,N~Site)
pairwise.t.test(dt3$N,dt3$Site)
aggregate(N~Site,dt3,mean)
aggregate(N~Site,dt3,sem)

#Transects stats Table S1
transects<-read_xlsx("Data/transects.xlsx")
transects2<-ddply(transects,c('Site','transect','Category'),summarise,sum=sum(Percent))
transects3<-spread(transects2,Category,sum)
transects3<-transects3 %>% replace(is.na(.), 0)
adonis2(transects3[,3:7]~transects3$Site,data=transects3)
pairwise.adonis(transects3[,3:7], factors = transects3$Site, p.adjust.m = "BH")
shapiro.test(transects3$'Live coral')
kruskal.test(transects3$'Live coral',transects3$Site)
pairwise.wilcox.test(transects3$'Live coral',transects3$Site)
aggregate(transects3$'Live coral'~Site,transects3,mean)
aggregate(transects3$'Live coral'~Site,transects3,sem)
shapiro.test(transects3$Macroalgae)
bartlett.test(Macroalgae ~ Site, data = transects3)
anova_test(transects3,Macroalgae~Site)
pairwise.t.test(transects3$Macroalgae,transects3$Site)
aggregate(Macroalgae~Site,transects3,mean)
aggregate(Macroalgae~Site,transects3,sem)

transects4<-ddply(transects[transects$Genus!='NA',],c('Site','transect','Genus'),summarise,sum=sum(Percent))
transects5<-spread(transects4,Genus,sum)
transects5<-transects5 %>% replace(is.na(.), 0)
transects5$sum<-rowSums(transects5[,3:9])
transects6<-transects5
transects6[,3:9]<-transects5[,3:9]*100/transects5[,10]
transects6$sum<-rowSums(transects6[,3:9])
adonis2(transects6[,3:9]~transects6$Site,data=transects6)
pairwise.adonis(transects6[,3:9], factors = transects6$Site, p.adjust.m = "BH")
pairwise.wilcox.test(transects6$Pavona,transects6$Site)
aggregate(Pavona~Site,transects6,mean)
aggregate(Pavona~Site,transects6,sem)
pairwise.wilcox.test(transects6$Porites,transects6$Site)
aggregate(Porites~Site,transects6,mean)
aggregate(Porites~Site,transects6,sem)

genera<-ddply(transects,c('Site','transect','Genus'),summarise,n=length(Genus))
genera2<-ddply(genera,c('Site','transect'),summarise,genus=length(Genus)-1)
shapiro.test(genera2$genus)
bartlett.test(genus ~ Site, data = genera2)
anova_test(genera2,genus~Site)
aggregate(genus~Site,genera2,mean)
aggregate(genus~Site,genera2,sem)
mean(genera2$genus)
sem(genera2$genus)

cyth<-transects %>% filter(ID2=="CYTHEREA")
cyth2<-ddply(cyth,c('Site','transect'),summarise,sum=sum(Percent))
cyth2<-rbind(cyth2,c("Linareva",3,0))
cyth2<-rbind(cyth2,c("Linareva",6,0))
mean(as.numeric(cyth2$sum))
sem(as.numeric(cyth2$sum))

##Transects Graphs
bentcolors<-c("Dead coral"="#E4E8F5","Fauna"="#f09ebc", "Live coral"="#F55839CC", "Macroalgae"="#82ad45", "Sediment"="#EDC285")
bentcolors2<-c("#FAB908","#F2D679","#a3ad32","#c286db","#e68540","#d450a6","#6FCAF7")

transects2<-ddply(transects,c('Site','Category'),summarise,sum=sum(Percent)/6)
bent1<-ggplot(data=transects2, aes(fill=Category, y=sum, x=fct_rev(Site)))+ 
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values=bentcolors)+
  labs(fill="Benthic categories", x=NULL,y="Cover percentage",title = "A. Benthic composition")+
  theme(panel.background = element_rect(fill='white'),
        plot.title = element_text(hjust = 0.5))

transects4<-ddply(transects,c('Site','Genus'),summarise,sum=sum(Percent)/6)
transects4<-transects4 %>% filter(Genus!="NA")
bent2<-ggplot(data=transects4, aes(fill=Genus, y=sum, x=fct_rev(Site)))+ 
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values=bentcolors2)+
  labs(fill="Live coral genus", x=NULL,y="Coral percentage",title = "B. Coral composition")+
  theme(panel.background = element_rect(fill='white'),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.33))
bent<-ggarrange(bent1,bent2)
bent
ggsave("Output/benthos.pdf",bent,width = 9, height = 4)
ggsave("Output/benthos.svg",bent,width = 9, height = 4)

##Spiderplot
df<-ddply(totalmean,("Site"),summarise,Mean_Temp=mean(Mean), Daily_Range=mean(Range))
df$Coral_Cover<-c(mean(transects3$'Live coral'[transects3$Site=="Linareva"]),
                  mean(transects3$'Live coral'[transects3$Site=="Mahana"]),
                  mean(transects3$'Live coral'[transects3$Site=="Manava"]),0)
df$Algal_Cover<-c(mean(transects3$Macroalgae[transects3$Site=="Linareva"]),
                  mean(transects3$Macroalgae[transects3$Site=="Mahana"]),
                  mean(transects3$Macroalgae[transects3$Site=="Manava"]),0)
df$Generic_Richness<-c(mean(genera2$genus[genera2$Site=="Linareva"]),
                       mean(genera2$genus[genera2$Site=="Mahana"]),
                       mean(genera2$genus[genera2$Site=="Manava"]),0)
df$Nitrogen<-c(mean(dt3$N[dt3$Site=="Linareva"]),
               mean(dt3$N[dt3$Site=="Mahana"]),
               mean(dt3$N[dt3$Site=="Manava"]),
               mean(dt3$N[dt3$Site=="Nursery"]))
max_min <- data.frame(
  Mean_Temp=c(28,26),
  Daily_Range=c(2,0),
  Coral_Cover=c(100,0),
  Algal_Cover=c(100,0),
  Generic_Richness=c(10,0),
  Nitrogen=c(0.8,0.6))
rownames(max_min) <- c("Max", "Min")

df2<-df[,-1]
row.names(df2)<-df$Site
df2<-rbind(max_min,df2)
colnames(df2)

create_beautiful_radarchart <- function(data, color = c("#BC4B07","#D3BC38","#5BA7A3","#424242")  , 
                                        vlabels = c( "Mean Temperature \n (26-28°C)",
                                                     "Daily Range \n (0-2°C)",
                                                     "Coral cover \n (0-100%)",
                                                     "Algal cover \n (0-100%)",
                                                     "Coral diversity \n (0-10 genera)",
                                                     "Nitrogen \n (0.6-0.8%)"), vlcex = 0.7,
                                        caxislabels = c(0, 25, 50, 75, '100'), title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.1), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


op <- par(mar = c(1, 1, 1, 1))
create_beautiful_radarchart(df2)
legend(
  x = "right", legend = rownames(df2[-c(1,2),]), horiz = F,
  bty = "n", pch = 15 , col = c("#BC4B07","#D3BC38","#5BA7A3","#424242"),
  text.col = "black", cex = 0.8, pt.cex = 1,y.intersp=2
)

## save as pdf width 9 height 5
## save as svg width 570 height 360

##spiderplot with added nutrients
df<-ddply(totalmean,("Site"),summarise,Mean_Temp=mean(Mean), Daily_Range=mean(Range))
df$Coral_Cover<-c(mean(transects3$'Live coral'[transects3$Site=="Linareva"]),
                  mean(transects3$'Live coral'[transects3$Site=="Mahana"]),
                  mean(transects3$'Live coral'[transects3$Site=="Manava"]),0)
df$Algal_Cover<-c(mean(transects3$Macroalgae[transects3$Site=="Linareva"]),
                  mean(transects3$Macroalgae[transects3$Site=="Mahana"]),
                  mean(transects3$Macroalgae[transects3$Site=="Manava"]),0)
df$Generic_Richness<-c(mean(genera2$genus[genera2$Site=="Linareva"]),
                       mean(genera2$genus[genera2$Site=="Mahana"]),
                       mean(genera2$genus[genera2$Site=="Manava"]),0)
df$Nitrogen<-c(mean(dt3$N[dt3$Site=="Linareva"]),
               mean(dt3$N[dt3$Site=="Mahana"]),
               mean(dt3$N[dt3$Site=="Manava"]),
               mean(dt3$N[dt3$Site=="Nursery"]))
df$NH4<-c(mean(sand$NH4[sand$Site=="Linareva"]),
               mean(sand$NH4[sand$Site=="Mahana"]),
               mean(sand$NH4[sand$Site=="Manava"]),
               mean(sand$NH4[sand$Site=="Nursery"]))
df$SIO<-c(mean(sand$'SI(OH)4'[sand$Site=="Linareva"]),
          mean(sand$'SI(OH)4'[sand$Site=="Mahana"]),
          mean(sand$'SI(OH)4'[sand$Site=="Manava"]),
          mean(sand$'SI(OH)4'[sand$Site=="Nursery"]))
df$NO3<-c(mean(water$NO3[water$Site=="Linareva"]),
          mean(water$NO3[water$Site=="Mahana"]),
          mean(water$NO3[water$Site=="Manava"]),
          mean(water$NO3[water$Site=="Nursery"]))
df$PO4<-c(mean(water$PO4[water$Site=="Linareva"]),
          mean(water$PO4[water$Site=="Mahana"]),
          mean(water$PO4[water$Site=="Manava"]),
          mean(water$PO4[water$Site=="Nursery"]))

max_min <- data.frame(
  Mean_Temp=c(28,26),
  Daily_Range=c(2,0),
  Coral_Cover=c(100,0),
  Algal_Cover=c(100,0),
  Generic_Richness=c(10,0),
  Nitrogen=c(0.8,0.6),
  NH4=c(13,0),
  SIO=c(10,0),
  NO3=c(0.5,0),
  PO4=c(0.5,0))
rownames(max_min) <- c("Max", "Min")

df2<-df[,-1]
row.names(df2)<-df$Site
df2<-rbind(max_min,df2)
colnames(df2)

create_beautiful_radarchart <- function(data, color = c("#BC4B07","#D3BC38","#5BA7A3","#424242")  , 
                                        vlabels = c( "Mean Temperature \n (26-28°C)",
                                                     "Daily Range \n (0-2°C)",
                                                     "Coral cover \n (0-100%)",
                                                     "Algal cover \n (0-100%)",
                                                     "Coral diversity \n (0-10 genera)",
                                                     "Nitrogen \n (0.6-0.8%)",
                                                     "NH4 in sand \n (0-13 µmol/L)",
                                                     "Si(OH)4 in sand \n (0-10 µmol/L)",
                                                     "N03 in water \n (0-0.5 µmol/L)",
                                                     "PO4 in water \n (0-0.5 µmol/L)"),
                                        vlcex = 0.7,
                                        caxislabels = c(0, 25, 50, 75, '100'), title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.1), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


op <- par(mar = c(1, 1, 1, 1))
create_beautiful_radarchart(df2)
legend(
  x = "right", legend = rownames(df2[-c(1,2),]), horiz = F,
  bty = "n", pch = 15 , col = c("#BC4B07","#D3BC38","#5BA7A3","#424242"),
  text.col = "black", cex = 0.8, pt.cex = 1,y.intersp=2
)

## save as pdf width 9 height 5
## save as svg width 570 height 360



##Survival
monitoring<-read_xlsx("Data/monitoring.xlsx")
surv<-monitoring %>% filter(dead==1)
surv1<-monitoring %>% filter (month==27 & live==1)
surv2<-rbind(surv,surv1)
surv2$title<-"Survival"
surv2$genotype<-as.factor(surv2$genotype)

cox <- coxph(Surv(month, dead) ~ site*genotype*treatment, data = surv2)
Anova(cox)

surv3<-surv2 %>% filter(treatment=="Transplant")
surv3$time<-as.factor(surv3$month)
cox <- coxph(Surv(month, dead) ~ site*genotype*time, data = surv3)
Anova(cox)
cox <- coxph(Surv(month, dead) ~ time, data = surv3)
pwpm(emmeans(cox,"time"))
cox <- coxph(Surv(month, dead) ~ genotype, data = surv3)
pwpm(emmeans(cox,"genotype"))

fit <- survfit(Surv(month, dead) ~ site, data = surv2)
survbysite<-ggsurvplot(fit, data = surv2,
                       conf.int = TRUE,
                       conf.int.style="ribbon",
                       conf.int.alpha=0.1,
                       xlab="Months",
                       legend.title="Site",
                       legend.labs=c("Linareva","Mahana","Manava","Nursery"),
                       palette=cytcolors,
                       legend= c(0.2, 0.2))+
  guides(color = guide_legend(override.aes = list(shape = NA)))
survbysite
ggsave("Output/survival.pdf",survbysite$plot,width = 6, height = 5)
ggsave("Output/survival.svg",survbysite$plot,width = 6, height = 5)


fit <- survfit(Surv(month, dead) ~ site+genotype, data = surv2)
survplot<-ggsurvplot(fit, data = surv2,
                     conf.int = TRUE,
                     conf.int.style="ribbon",
                     conf.int.alpha=0.1,
                     xlab="Months",
                     legend.title="Site",
                     facet.by="genotype",
                     legend.labs=c("Linareva","Mahana","Manava","Nursery"),
                     palette=cytcolors,
                     legend= c(0.85, 0.13),
                     alpha=0.9)+
  guides(color = guide_legend(override.aes = list(shape = NA)))
survplot
survplot2 <- survplot + 
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 14))
survplot2  

ggsave("Output/survgen.pdf",survplot2,width = 7, height = 7)
ggsave("Output/survgen.svg",survplot2,width = 7, height = 7)

#can we group genotypes with similar survival curves in transplantation sites?
surv4<-surv3[surv3$genotype=="10",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
surv4<-surv3[surv3$genotype=="2",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
pwpm(emmeans(cox,"site"))
surv4<-surv3[surv3$genotype=="3",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
pwpm(emmeans(cox,"site"))
surv4<-surv3[surv3$genotype=="4",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
surv4<-surv3[surv3$genotype=="5",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
surv4<-surv3[surv3$genotype=="7",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
surv4<-surv3[surv3$genotype=="8",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)
surv4<-surv3[surv3$genotype=="9",]
cox <- coxph(Surv(month, dead) ~ site, data = surv4)
Anova(cox)

surv4<-surv3[surv3$site=="Linareva",]
cox <- coxph(Surv(month, dead) ~ genotype, data = surv4)
Anova(cox)
surv4<-surv3[surv3$site=="Mahana",]
cox <- coxph(Surv(month, dead) ~ genotype, data = surv4)
Anova(cox)
surv4<-surv3[surv3$site=="Manava",]
cox <- coxph(Surv(month, dead) ~ genotype, data = surv4)
Anova(cox)
pwpm(emmeans(cox,"genotype"))

growth<-read_xlsx("Data/GROWTH.xlsx")
growth$Genotype<-as.factor(growth$Genotype)
growth$colony<-as.factor(growth$colony)
growth$days<-as.numeric(growth$date2-growth$date1)
growth$gf<-(growth$final_size-growth$Initial_size)*100/(growth$Initial_size*growth$days)
growth1<-growth %>% filter(Time!="total" & Treatment=="Transplant")
sumg<-ddply(growth1,c("Site","State"),summarize, total=length(State))
bysite<-aggregate(Time~Site,growth1,length)
sumg<-left_join(sumg,bysite,by="Site")
sumg$percent<-sumg$total*100/sumg$Time

mod1<-glmer(survival~final_size+Genotype*Site+Time+State+gf+(1|colony),growth1,family=binomial)
Anova(mod1,type=3)
emmeans(mod1,pairwise~Time)
aggregate(survival~Time, growth1, mean)

a<-growth1 %>% mutate(prob=survival) %>% ggplot(aes(final_size,prob))+
  geom_point(alpha=0.2)+
  geom_smooth(method="glm",method.args=list(family="binomial"),color="#465F3D")+
  xlab(expression(paste("Surface in cm"^{2})))+
  ylab('Probability to survive')+
  theme_bw()+
  xlim(0,800)+
  scale_y_continuous(breaks = seq(0,1,0.2))
a
growth1<-growth %>% filter(Time!="total" & Treatment=="Control")
mod1<-glmer(survival~final_size+(1|colony),growth1,family=binomial)
Anova(mod1)

aggregate(Initial_size~Time, growth1, mean)
aggregate(Initial_size~Time, growth1, sem)
growth2<-growth1 %>% filter(Time=="March 2021")

##Growth
growth<-read_xlsx("Data/GROWTH.xlsx")
growth$Genotype<-as.factor(growth$Genotype)
growth$Year<-as.factor(growth$Year)
growth$days<-as.numeric(growth$date2-growth$date1)
growth$gf<-(growth$final_size-growth$Initial_size)*100/(growth$Initial_size*growth$days)

## models
finalgrowth<-growth %>% filter(Time=="total")
hist(finalgrowth$gf)
final<-finalgrowth %>% filter(Site!="Nursery")
mod1<-lm(gf~Genotype*Site*State+Initial_size,final)
Anova(mod1)
plot(mod1)
final<-finalgrowth %>% filter(site!="Nursery" & State=="OK")
mod1<-lm(gf~Genotype*Site+Initial_size,final)
Anova(mod1)
plot(mod1)

growth1<-growth %>% filter(Time!="total")
mod1<-lmer(gf~State+Genotype*Time*Treatment+Initial_size+(1|colony),growth1)
Anova(mod1)
aov<-Anova(mod1)
aov$predictor<-row.names(aov)
aov$Chisq<-round(aov$Chisq,2)
aov$'Pr(>Chisq)'<-round(aov$'Pr(>Chisq)',3)
aov$'Pr(>Chisq)'[aov$'Pr(>Chisq)'==0.000]<-"<0.001"
aov$predictor<-str_replace_all(aov$predictor,":","*")
aov$predictor<-str_replace_all(aov$predictor,"_"," ")
write_xlsx(aov,"Output/Growthmodel.xlsx")

growth2<-growth1 %>% filter(State=="OK")
mod1<-lmer(gf~Initial_size+Genotype*Time*Treatment+(1|colony),growth2)
Anova(mod1)
aov<-Anova(mod1)
aov$predictor<-row.names(aov)
aov$Chisq<-round(aov$Chisq,2)
aov$'Pr(>Chisq)'<-round(aov$'Pr(>Chisq)',3)
aov$'Pr(>Chisq)'[aov$'Pr(>Chisq)'==0.000]<-"<0.001"
aov$predictor<-str_replace_all(aov$predictor,":","*")
write_xlsx(aov,"Output/Growthmodel.xlsx")

growth2<-growth1 %>% filter(Treatment=="Transplant" & State=="OK")
mod1<-lmer(gf~Initial_size+Genotype*Site*Time+(1|colony),growth2)
Anova(mod1)
aov<-Anova(mod1)
aov$predictor<-row.names(aov)
aov$Chisq<-round(aov$Chisq,2)
aov$'Pr(>Chisq)'<-round(aov$'Pr(>Chisq)',3)
aov$'Pr(>Chisq)'[aov$'Pr(>Chisq)'==0.000]<-"<0.001"
aov$predictor<-str_replace_all(aov$predictor,":","*")
write_xlsx(aov,"Output/Growthmodel.xlsx")
pairwise.wilcox.test(growth2$gf,growth2$Site)

growth2<-growth1 %>% filter(Treatment=="Transplant" & State!="OK")
mod1<-lmer(gf~Initial_size+Genotype*Site*Time+(1|colony),growth2)
Anova(mod1)
aov<-Anova(mod1)
aov$predictor<-row.names(aov)
aov$Chisq<-round(aov$Chisq,2)
aov$'Pr(>Chisq)'<-round(aov$'Pr(>Chisq)',3)
aov$'Pr(>Chisq)'[aov$'Pr(>Chisq)'==0.000]<-"<0.001"
aov$predictor<-str_replace_all(aov$predictor,":","*")
write_xlsx(aov,"Output/Growthmodel.xlsx")

##interaction site*period
growth2<-growth1 %>% filter(Site=="Linareva" & State=="OK")
pairwise.wilcox.test(growth2$gf, growth2$Time)
growth2<-growth1 %>% filter(Site=="Mahana" & State=="OK")
pairwise.wilcox.test(growth2$gf, growth2$Time)
ggplot(growth2, aes(x = Time,y=gf))+
  geom_boxplot()
growth2<-growth1 %>% filter(Site=="Manava" & State=="OK")
pairwise.wilcox.test(growth2$gf, growth2$Time)

growth1<-growth %>% filter(Time!="total" & Treatment=="Transplant")
growth2<-growth1 %>% filter(Time=="March 2021")
pairwise.wilcox.test(growth2$gf, growth2$Site)
growth2<-growth1 %>% filter(Time=="June 2021")
growth2<-growth1 %>% filter(Time=="October 2021")
growth2<-growth1 %>% filter(Time=="February 2022")
growth2<-growth1 %>% filter(Time=="June 2022")
growth2<-growth1 %>% filter(Time=="October 2022")
growth2<-growth1 %>% filter(Time=="February 2023")
growth2<-growth1 %>% filter(Time=="April 2023")
aggregate(gf~Year,growth1,mean)
aggregate(gf~Year,growth1,sem)

growth1<-growth %>% filter(Time!="total" & Treatment=="Transplant" & State=="OK")
growth2<-growth1 %>% filter(Time=="October 2021")
pairwise.wilcox.test(growth2$gf, growth2$Site)
aggregate(gf~Year,growth1,mean)
aggregate(gf~Year,growth1,sem)


growth1<-growth %>% filter(Time!="total")
growth2<-growth1 %>% filter(Site=="Nursery")
cor.test(growth2$gf,growth2$Initial_size)
growth2<-growth1 %>% filter(Site!="Nursery")
cor.test(growth2$gf,growth2$Initial_size)

growth1<-growth %>% filter(Time!="total")
growth2<-ddply(growth1,c("Site","Time"),summarise,mean=mean(gf),sem=sem(gf))
growth2$mean<-as.numeric(growth2$mean)
growth2$sem<-as.numeric(growth2$sem)
Months<-data.frame(c("March 2021","June 2021","October 2021",
  "February 2022","June 2022","October 2022",
  "February 2023","April 2023"),c(1,4,8,12,16,20,24,26))
colnames(Months)[1]<-"Time"
colnames(Months)[2]<-"month"
growth2$month<-growth2$Time
growth2<-FindReplace(growth2,"month",Months,from="Time",to="month",exact=TRUE,vector=FALSE)
initial<-data.frame("Site"=c("Linareva","Mahana","Manava","Nursery"),
                    "Time"=c(0,0,0,0),
                    "mean"=c(0.203,0.203,0.203,0.203),
                    "sem"=c(0.048,0.048,0.048,0.048),
                    "month"=c(0,0,0,0))
growth3<-rbind(growth2,initial)
growth3$month<-factor(growth3$month,c("0","1","4","8","12","16","20","24","26"))

fgrowth<-ggplot(growth3, aes(x = month,y=mean,group=Site))+
  geom_line(aes(col=Site))+
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem,fill=Site),alpha=0.4) +
  ylab("Growth (% per day)")+
  xlab("Months after transplantation")+
  theme_bw()+
  theme(legend.position = c(0.7, 0.2))+
  scale_color_manual(values=cytcolors)+
  scale_fill_manual(values=cytcolors)
fgrowth

fgrowth<-ggplot(growth3, aes(x = month,y=mean,group=Site))+
  geom_line(aes(col=Site))+
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem,fill=Site),alpha=0.4) +
  ylab("Growth (% per day)")+
  xlab("Months after transplantation")+
  theme_bw()+
  theme(legend.position = c(0.7, 0.2))+
  scale_color_manual(values=cytcolors)+
  scale_fill_manual(values=cytcolors)
fgrowth


yearly<-read_xlsx("Data/yearlygrowth.xlsx")
yearly$genotype<-as.factor(yearly$genotype)
yearly$year<-paste("Year",yearly$year,sep=" ")
yearly$days<-as.numeric(yearly$date2-yearly$date1)
yearly$gf<-(yearly$final_size-yearly$initial_size)*100/(yearly$initial_size*yearly$days)
yearly2<-yearly %>% filter(site!="Nursery")
aggregate(gf~year, yearly2, mean)
aggregate(gf~year, yearly2, sem)
pairwise.wilcox.test(yearly2$gf,yearly2$site)

yearplot<-ggplot(yearly, aes(x=year,y=gf,fill=site)) +
  geom_boxplot()+
  theme_bw()+
  ylab("Growth (% per day)")+
  xlab("")+
  scale_fill_manual(values=cytcolors)+
  theme(legend.position = "none")
yearplot

gplot<-ggarrange(fgrowth,yearplot, widths = c(1, 0.5))
gplot
ggsave("Output/Fig4B.pdf",gplot,width = 11, height = 4)
ggsave("Output/Fig4B.svg",gplot,width = 11, height = 4)

##Maturity stats
maturity<-read_xlsx("Data/maturity.xlsx")
maturity$maturity<-0
maturity$maturity[maturity$mature=="MATURE"]<-1
maturity$genotype<-as.factor(maturity$genotype)
maturity$colony<-as.factor(maturity$colony)
maturity$year<-as.factor(maturity$year)
mod1<-glm(maturity~genotype+Treatment+year+size,maturity,family=binomial)
Anova(mod1)

maturity2<-maturity %>% filter(Treatment=="Transplant")
mod1<-glmer(maturity~genotype+site+size+year+(1|colony),maturity2,family=binomial)
Anova(mod1,type=3)
anova<-Anova(mod1)
aov<-as.data.frame(anova)
aov$comparison<-row.names(aov)
aov$`Pr(>Chisq)`<-ifelse(aov$`Pr(>Chisq)`<0.001,"<0.001",round(aov$`Pr(>Chisq)`,digits=3))
aov$'LR Chisq'<-round(aov$'LR Chisq',digits=2)
aov %>% write_xlsx("Output/maturityanova.xlsx")
mod1<-glmer(maturity~genotype+(1|colony),maturity2,family=binomial)
emmeans(mod1,pairwise~genotype)

maturity2<-maturity %>% filter(year=="2021" & Treatment=="Transplant")
kruskal.test(maturity2$maturity,maturity2$site)
mod1<-glm(maturity~size+genotype+site,maturity2,family=binomial)
Anova(mod1)
mod1<-glm(maturity~genotype,maturity2,family=binomial)
emmeans(mod1,pairwise~genotype)

maturity2<-maturity %>% filter(year=="2022" & Treatment=="Transplant")
mod1<-glm(maturity~size+genotype+site,maturity2,family=binomial)
Anova(mod1)
kruskal.test(maturity2$maturity,maturity2$site)
mod1<-glm(maturity~genotype,maturity2,family=binomial)
emmeans(mod1,pairwise~genotype)

maturity2<-maturity %>% filter(Treatment=="Control")
mod1<-glm(maturity~size+year,maturity2,family=binomial)
Anova(mod1)
kruskal.test(maturity2$mature,maturity2$site)
aggregate(size~mature+year,maturity2,mean)
aggregate(size~mature+year,maturity2,sem)
maturity2<-maturity %>% filter(Treatment=="Control" & year=="2022")
kruskal.test(maturity2$size,maturity2$mature)
maturity2<-maturity %>% filter(Treatment=="Control")


##Maturity graph
maturity2<-ddply(maturity,c("site","year"),summarise,percent=sum(maturity)/length(maturity)*100,se=sem(maturity)*100)
mat1<-ggplot(maturity2,aes(y=percent,x=site,fill=site)) +
  geom_bar(colour="black",stat="identity")+
  geom_linerange(aes(x=site, ymin=percent-se, ymax=percent+se))+
  theme_bw()+
  ylab("Percent maturity")+
  xlab("Site")+
  facet_grid(~year)+
  theme(legend.position = "none",
        axis.text.x = element_blank())+
  xlab("")+
  scale_fill_manual(values=cytcolors)
mat1

maturity$mature<-str_to_sentence(maturity$mature)
maturity3<-ddply(maturity[maturity$size>0,],c("site","year","mature"),summarise,meansize=mean(size,na.rm=TRUE),
                 se=sem(size))
dodge <- position_dodge2(width = 0.9,preserve = "single")
mat2<-ggplot(maturity3[-15,], aes(x=mature,y=meansize, fill=site)) +
  geom_bar(colour="black",stat="identity",position =dodge)+
  geom_linerange(aes(x=mature, ymin=meansize-se, ymax=meansize+se),
                 position =dodge)+
  ylab(expression(paste("Surface in cm"^{2})))+
  theme_bw()+
  theme(legend.position = c(0.2,0.65),axis.title.x = element_blank())+
  scale_fill_manual(values=cytcolors)+
  facet_grid(~year)
mat2 
mat<-ggarrange(mat1,mat2,widths=c(0.5,1))
mat
ggsave("Output/maturity.pdf",mat,width = 10, height = 3)
ggsave("Output/maturity.svg",mat,width = 10, height = 3)

##binary logistic regression to predict size for 95% maturity
maturity2<-maturity %>% filter(Treatment=="Transplant")
b<-maturity2 %>% mutate(prob=maturity) %>% ggplot(aes(size,prob))+
  geom_point(alpha=0.2)+
  geom_smooth(method="glm",method.args=list(family="binomial"),color="#465F3D")+
  xlab(expression(paste("Surface in cm"^{2})))+
  ylab('Probability to be mature')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,1,0.2))
b
binom<-ggarrange(a,b)
binom
ggsave("Output/sizeregression.pdf",binom,width = 7, height = 3)
ggsave("Output/sizeregression.svg",binom,width = 7, height = 3)

##Symbiont diversity
sites<-ddply(monitoring,c("site","colony","genotype"),summarize,n=length(colony))
Abundance <- read.table(file="Data/Seqs_relative.txt", sep = "\t", fill=T,header=T,dec=",")
Abundance<-Abundance[,c(2, 40:240)]
Abundance$colony<-map(strsplit(Abundance$sample_name, split = "_"), 1)
Abundance$time<-map(strsplit(Abundance$sample_name, split = "_"),2)
Abundance$sample_name[duplicated(Abundance$sample_name)]
Abundance<-Abundance[-c(190,194,195),]
Abundance$genotype<-substr(Abundance$colony,1,nchar(Abundance$colony)-2)
Abundance$site<-as.character(Abundance$colony)
Abundance<-FindReplace(Abundance,"site",sites,"colony","site",exact=TRUE, vector=FALSE)
Abundance[nchar(Abundance$site)<5,206]<-"Nursery"

Abund2<-Abundance
Abund2[,2:202] <- mutate_all(Abund2[,2:202], function(x) as.numeric(as.character(x)))
Abund2[,2:202]<-ceiling(Abund2[,2:202])
Abund2$sum<-rowSums(Abund2[,2:202])
mean(Abund2$sum)
sem(Abund2$sum)



##20 most abundant classes
Abund<-gather(Abundance,strain,percent,c(2:202))
ITS2 <- mutate_all(Abundance[,2:202], function(x) as.numeric(as.character(x)))
ITS2<-ITS2 %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
max(ITS2[193,1:201])
sum<-read.table(file="Data/Seqs_relative.txt", sep = "\t", fill=T,header=T,dec=",")
sum<-sum[1,40:240]
sum[1,]<-ITS2[193,]
sum <- mutate_all(sum, function(x) as.numeric(as.character(x)))
sum2<-gather(sum,strain,percent,c(1:201))
sum2$clade<-"A"
sum2$clade[79:148]<-"C"
sum2$clade[149:201]<-"D"
sum(sum2$percent)
sum3<-ddply(sum2,c("clade"),summarize,percent=sum(percent)/192)
sum(sum3$percent)
sum2$new<-"Other A"
sum2$new[79:148]<-"Other C"
sum2$new[149:201]<-"Other D"
unique(sum2$new)
sum2$new[sum2$percent>1.1]<-sum2$strain[sum2$percent>1.1]
unique(sum2$new)
unique(sum2$strain)
sum2<-as.data.frame(sum2)

ITS2_2<-Abund %>% filter(percent>0)
ITS2_2<-as.data.frame(ITS2_2)
ITS2_2<-FindReplace(ITS2_2,"strain",sum2,"strain","new",exact=TRUE, vector=FALSE)
unique(ITS2_2$strain)
ITS2_2$percent<-as.numeric(ITS2_2$percent)
ITS2_2$sample_name<-as.character(ITS2_2$sample_name)
ITS2_2$colony<-as.character(ITS2_2$colony)
ITS2_2$time<-as.character(ITS2_2$time)
ITS2_3<-ddply(ITS2_2,c("sample_name","colony","time","site","genotype","strain"),summarise,percent=sum(percent))
ITS2_4<-spread(ITS2_3,strain,percent)
ITS2_4[is.na(ITS2_4)]<-0
ITS2_4<-ITS2_4[-32,]
ITS2_4<-ITS2_4 %>% filter(colony!="421")
colnames(ITS2_4)[24]<-"C175"
colnames(ITS2_4)[25]<-"C53"

colors2<-c("#465F3D","#BC4B07" , "#5B1D91","#8BF0A4" ,"#D3BC38", "#6E4465", "#328F37","#E3423F","#F2D679" ,
           "#444E6E", "#B88BF7", "#7DFAF6","#E6A640", "#FAB908", "#5BA7A3","#E03E97", "" , "#1B4217",
           "#732D3D", "#6FCAF7")
colors3<-c("#BC4B07","#6E4465","#D3BC38",  "#328F37","#a3ad32","#46bf3d","#e3d12b","#8bf571",
            "#0a6e24","#E6A640","#234f20", "#FAB908","#33914b","#F2D679","#1e73b0","#3ca3bd", "#3453A8" ,
           "#732D3D","#1B4217","#6FCAF7")


ITS2_4$name<-substr(ITS2_4$site, 0, 3)
graph<-gather(ITS2_4,strain,percent,c(6:25))
graph$genotype<-factor(graph$genotype,c("2","3","4","5","7","8","9","10"))
fig<-ggplot(graph, aes(x=colony, y=percent*100, fill=strain)) +
  geom_bar(stat="identity", width=0.9)+
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),legend.text.align = 0,
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines"))+
  facet_grid(time~genotype+fct_rev(name),scale='free')+
  ylab("Symbiont strain relative abundance %")+
  scale_fill_manual(values=colors3)+
  xlab("Colony")
fig

ggsave("Output/symbionts.jpg",fig,width = 10, height = 5.5)
ggsave("Output/symbionts.pdf",fig,width = 10, height = 5.5)
ggsave("Output/symbionts.svg",fig,width = 10, height = 5.5)

#ITS2 profiles
ITS2Profil <- read.table("Data/ITS2Profiles.txt", sep = "\t", fill=T,header=T,dec=",")
ITS2Profil <- ITS2Profil[-c(1:5,13,196,200:202),c(2:length(ITS2Profil))]
names(ITS2Profil) <- as.character(unlist(ITS2Profil[1,]))
ITS2Profil<-ITS2Profil[-1,]
colnames(ITS2Profil)[1]<-"sample_name"

ITS2Profil$colony<-map(strsplit(ITS2Profil$sample_name, split = "_"), 1)
ITS2Profil$time<-map(strsplit(ITS2Profil$sample_name, split = "_"),2)
ITS2Profil$genotype<-substr(ITS2Profil$colony,1,nchar(ITS2Profil$colony)-2)
ITS2Profil$site<-as.character(ITS2Profil$colony)
ITS2Profil<-FindReplace(ITS2Profil,"site",sites,"colony","site",exact=TRUE, vector=FALSE)
ITS2Profil[nchar(ITS2Profil$site)<5,35]<-"Nursery"
unique(ITS2Profil$site)
##remove A16 profile which is 0
ITS2Profil<-ITS2Profil[,-11]
ITS2Profil$name<-substr(ITS2Profil$site, 0, 3)

ITS2Profil[,2:30] <- mutate_all(ITS2Profil[,2:30], function(x) as.numeric(as.character(x)))
ITS2Profil$sum<-rowSums(ITS2Profil[,2:30])
ITS2Profil[,2:30]<-ITS2Profil[,2:30]/ITS2Profil[,36]
ITS2Profil$sum2<-rowSums(ITS2Profil[,2:30])

Profilesum <- mutate_all(ITS2Profil[,2:30], function(x) as.numeric(as.character(x)))
Profilesum<- Profilesum%>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
which(Profilesum[192,]<1)+1

ITS2Profil2<-ITS2Profil[,-c(which(Profilesum[192,]<1)+1)]
ITS2Profil2$sum<-rowSums(ITS2Profil2[,2:19])
ITS2Profil2[,2:19]<-ITS2Profil2[,2:19]/ITS2Profil2[,25]
ITS2Profil2$sum2<-rowSums(ITS2Profil2[,2:19])

graph<-gather(ITS2Profil2,profile,percent,c(2:19))
graph$percent<-as.numeric(graph$percent)
graph$colony<-as.character(graph$colony)
graph$time<-as.character(graph$time)
graph$genotype<-factor(graph$genotype,c("2","3","4","5","7","8","9","10"))
colors18<-c("#BC4B07","#6E4465","#db6c04","#c70223","#732D3D","#E3423F","#f09ed6",
           "#D3BC38",  "#328F37","#a3ad32","#46bf3d","#e3d12b","#234f20",
           "#1e73b0","#3ca3bd", "#3453A8" ,"#1623cc","#6FCAF7")


fig<-ggplot(graph, aes(x=colony, y=percent*100, fill=profile)) +
  geom_bar(stat="identity", width=0.9)+
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),legend.text.align = 0,
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines"))+
  facet_grid(time~genotype+fct_rev(name),scale='free')+
  ylab("Symbiont profile relative abundance %")+
  scale_fill_manual(values=colors18)+
  xlab("Colony")+
  guides(fill=guide_legend(ncol =1))
fig

ggsave("Output/symbiontprofiles.jpg",fig,width = 12, height = 5)
ggsave("Output/symbiontprofiles.svg",fig,width = 12, height = 5)

Abund3<-ITS2Profil
Abund3[,2:30]<-ceiling(Abund3[,2:30])
Abund3$sum<-rowSums(Abund3[,2:30])
mean(Abund3$sum)
sem(Abund3$sum)

#stats with full diversity
Abundance<-as.data.frame(Abundance)
Abundance$time<-as.character(Abundance$time)
Abundance$genotype<-as.character(Abundance$genotype)
Abundance$site2<-Abundance$site
Abundance$site2[Abundance$time=="T0"]<-"Nursery"
Abundance[,2:202]<-mutate_all(Abundance[,2:202], function(x) as.numeric(as.character(x)))
Abundance$treatment<-"Transplant"
Abundance$treatment[Abundance$site2=="Nursery"]<-"Control"
Abundance$Origin<-Abundance$genotype
Abundance<-FindReplace(Abundance,"Origin",genotypes,from="Gen",to="Origin",exact=TRUE,vector=FALSE)
dist <- vegdist(Abundance[,2:202], method="bray") 

permanova<-adonis2(dist ~ genotype*treatment*time, data=Abundance, permutations=9999)
permanova
permanova<-adonis2(dist ~ Origin+genotype*site2*time, data=Abundance, permutations=9999)
permanova<-adonis2(dist ~ genotype*site2*time, data=Abundance, permutations=9999)
permanova
model<-Abundance %>% filter(site2=="Nursery")
dist <- vegdist(model[,2:202], method="bray") 
permanova <- adonis2(dist ~ genotype*time, data=model, permutations=9999)
permanova

model<-Abundance %>% filter(time=="T0")
aggregate(time~genotype,model,length)
dist <- vegdist(model[,2:202], method="bray")
adonis2(dist ~ genotype*site, data=model, permutations=9999)
pairwise.adonis(dist, factors = model$genotype, p.adjust.m = "BH")
pairwise.adonis(dist, factors = model$site, p.adjust.m = "BH")

model<-Abundance %>% filter(time!="T0")
dist <- vegdist(model[,2:202], method="bray")
permanova <- adonis2(dist ~ genotype*site*time, data=model, permutations=9999)
permanova
model<-Abundance %>% filter(time=="T12")
dist <- vegdist(model[,2:202], method="bray")
permanova <- adonis2(dist ~ genotype*site, data=model, permutations=9999)
permanova
model<-Abundance %>% filter(time=="T27")
dist <- vegdist(model[,2:202], method="bray")
permanova <- adonis2(dist ~ genotype*site, data=model, permutations=9999)
permanova <- adonis2(dist ~ genotype*treatment, data=model, permutations=9999)
permanova

model<-Abundance %>% filter(time!="T27" & genotype=="7")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$site, p.adjust.m = "BH")
pairwise.adonis(dist, factors = model$treatment, p.adjust.m = "BH")
kruskal.test(model$A1ee,model$treatment)
kruskal.test(model$C3ae,model$time)

model<-Abundance %>% filter(genotype=="10")
dist <- vegdist(model[,2:202], method="bray")
adonis2(dist ~ site2*time, data=model, permutations=9999)
pairwise.adonis(dist, factors = model$time, p.adjust.m = "BH")
pairwise.adonis(dist, factors = model$site2, p.adjust.m = "BH")
pairwise.wilcox.test(model$A1ee,model$time)
pairwise.wilcox.test(model$D1,model$time)
pairwise.wilcox.test(model$C3ae,model$time)

model<-Abundance %>% filter(genotype=="10" & site=="Nursery")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$time, p.adjust.m = "BH")

model<-Abundance %>% filter(genotype=="10" & time=="T12")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$site, p.adjust.m = "BH")

model<-Abundance %>% filter(genotype=="3")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$time, p.adjust.m = "BH")

symsp<-mvabund(Abundance[,2:202])
mod1 <- manyglm(symsp ~ Abundance$genotype*Abundance$site2*Abundance$time, family = "negative_binomial")
anova(mod1)


symsp<-mvabund(model[,2:202])
mod1 <- manyglm(symsp ~ model$treatment, family = "negative_binomial")
anova(mod1)
mod1 <- manyglm(symsp ~ model$site, family = "negative_binomial")
anova(mod1)
anova(mod1, pairwise.comp=model$site)

model<-Abundance %>% filter(time=="T12" & site=="Nursery")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$genotype, p.adjust.m = "BH")
symsp<-mvabund(model[,2:202])
mod1 <- manyglm(symsp ~ model$genotype, family = "negative_binomial")
anova(mod1, pairwise.comp=model$genotype)

model<-Abundance %>% filter(time=="T12" & site=="Linareva")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$genotype, p.adjust.m = "BH")
symsp<-mvabund(model[,2:202])
mod1 <- manyglm(symsp ~ model$genotype, family = "negative_binomial")
anova(mod1, pairwise.comp=model$genotype)

model<-Abundance %>% filter(time=="T12" & site=="Mahana")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$genotype, p.adjust.m = "BH")
symsp<-mvabund(model[,2:202])
mod1 <- manyglm(symsp ~ model$genotype, family = "negative_binomial")
anova(mod1, pairwise.comp=model$genotype)

model<-Abundance %>% filter(time=="T12" & site=="Manava")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$genotype, p.adjust.m = "BH")
symsp<-mvabund(model[,2:202])
mod1 <- manyglm(symsp ~ model$genotype, family = "negative_binomial")
anova(mod1, pairwise.comp=model$genotype)

model<-Abundance %>% filter(site=="Manava")
dist <- vegdist(model[,2:202], method="bray")
adonis2(dist ~ genotype*time, data=model, permutations=9999)
pairwise.adonis(dist, factors = model$time, p.adjust.m = "BH")

model<-Abundance %>% filter(time=="T12" & genotype=="10")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$treatment, p.adjust.m = "BH")
kruskal.test(model$D1,model$treatment)

model<-Abundance %>% filter(time=="T27")
dist <- vegdist(model[,2:202], method="bray")
pairwise.adonis(dist, factors = model$site, p.adjust.m = "BH")

#link symbiont diversity with time of death and growth
ITS2_4$Lifespan<-ITS2_4$colony
ITS2_4<-FindReplace(ITS2_4,"Lifespan",surv,"colony","month",exact=TRUE, vector=FALSE)
ITS2_4$Lifespan[nchar(ITS2_4$Lifespan)>2]<-30
ITS2_4$Growth<-ITS2_4$colony
ITS2_4$Size<-ITS2_4$colony
ITS2_4$Maturity<-ITS2_4$colony
ITS2_4[ITS2_4$time=="T0",]<-FindReplace(ITS2_4[ITS2_4$time=="T0",],"Growth",yearly[yearly$year=="Year 1",],"colony","gf",exact=TRUE, vector=FALSE)
ITS2_4[ITS2_4$time=="T0",]<-FindReplace(ITS2_4[ITS2_4$time=="T0",],"Size",yearly[yearly$year=="Year 1",],"colony","initial_size",exact=TRUE, vector=FALSE)
ITS2_4[ITS2_4$time!="T0",]<-FindReplace(ITS2_4[ITS2_4$time!="T0",],"Growth",yearly[yearly$year=="Year 2",],"colony","gf",exact=TRUE, vector=FALSE)
ITS2_4[ITS2_4$time=="T12",]<-FindReplace(ITS2_4[ITS2_4$time=="T12",],"Size",yearly[yearly$year=="Year 2",],"colony","initial_size",exact=TRUE, vector=FALSE)
ITS2_4[ITS2_4$time=="T27",]<-FindReplace(ITS2_4[ITS2_4$time=="T27",],"Size",yearly[yearly$year=="Year 2",],"colony","final_size",exact=TRUE, vector=FALSE)
ITS2_4[nchar(ITS2_4$Growth)<5,]<-FindReplace(ITS2_4[nchar(ITS2_4$Growth)<5,],"Size",yearly[yearly$year=="Year 1",],"colony","final_size",exact=TRUE, vector=FALSE)
ITS2_4[nchar(ITS2_4$Growth)<5,]<-FindReplace(ITS2_4[nchar(ITS2_4$Growth)<5,],"Growth",yearly[yearly$year=="Year 1",],"colony","gf",exact=TRUE, vector=FALSE)
ITS2_4$Size[nchar(ITS2_4$Growth)<5]<-NA
ITS2_4$Growth[nchar(ITS2_4$Growth)<5]<-NA
maturity<-as.data.frame(maturity)
ITS2_4[ITS2_4$time!="T0",]<-FindReplace(ITS2_4[ITS2_4$time!="T0",],"Maturity",maturity[maturity$year=="2022",],"colony","maturity",exact=TRUE, vector=FALSE)
ITS2_4$Maturity[nchar(ITS2_4$Maturity)>1]<-NA

ITS2_4[,c(6:25,27:30)]<-mutate_all(ITS2_4[,c(6:25,27:30)], function(x) as.numeric(as.character(x)))
M=cor(ITS2_4[,c(6:25,27:30)],use="pairwise.complete.obs")
res1 <- cor.mtest(ITS2_4[,c(6:25,27:30)], conf.level = 0.95)
pvalues<-as.data.frame(res1$p)
corrplot(M, method = 'color',tl.col = 'black',p.mat = res1$p,sig.level = 0.05,insig = "blank")
#save svg 500x500
##simplify dataframe by grouping clades
ITS2_4$CladeA<-rowSums(ITS2_4[,c(6,7,21)])
ITS2_4$CladeC<-rowSums(ITS2_4[,c(8:17,22,24,25)])
ITS2_4$CladeD<-rowSums(ITS2_4[,c(18:20,23)])
data<-ITS2_4[,-c(6:26)]
data$Death<-1
data$Death[data$Lifespan==30]<-0


graph<-gather(data,strain,percent,c(9:11))
graph$genotype<-factor(graph$genotype,c("2","3","4","5","7","8","9","10"))
fig<-ggplot(graph, aes(x=colony, y=percent*100, fill=strain)) +
  geom_bar(stat="identity", width=0.9)+
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),legend.text.align = 0,
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines"))+
  facet_grid(time~genotype+fct_rev(site),scale='free')+
  ylab("Symbiont strain relative abundance %")+
  xlab("Colony")
fig

M=cor(data[,6:13],use="pairwise.complete.obs")
res1 <- cor.mtest(data[,6:13], conf.level = 0.95)
pvalues<-as.data.frame(res1$p)
corrplot(M, method = 'color',type = 'lower',tl.col = 'black',p.mat = res1$p,sig.level = 0.05,insig = "blank")
#save svg 600x520
data2<-data %>% filter(site!="Nursery")
M=cor(data2[,6:13],use="pairwise.complete.obs")
res1 <- cor.mtest(data2[,6:13], conf.level = 0.95)
pvalues<-as.data.frame(res1$p)
corrplot(M, method = 'color',type = 'lower',tl.col = 'black',p.mat = res1$p,sig.level = 0.05,insig = "blank")

##do initial clades influence survival pattern
data0<-data[data$time=="T0",c(2,9:12)]
data0$colony<-as.factor(data0$colony)
surv2$colony<-as.factor(surv2$colony)
surv3<-inner_join(surv2,data0,by="colony")
cox <- coxph(Surv(month, dead) ~ site*genotype+CladeA+CladeC+CladeD, data = surv3)
Anova(cox)

##PCA

pca1<-PCA(data[,6:13])
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
data$pc1 <- pca1$ind$coord[, 1]
data$pc2 <- pca1$ind$coord[, 2]
data$Genotype<-data$genotype
data$Genotype<-factor(data$Genotype, c("2","3","4","5","7","8","9","10"))
data$Site<-data$site
data$Site[data$time=="T0"]<-"Nursery"
data$Time<-data$time

Pca1<-ggplot(data = data, aes(x = pc1, y = pc2,color=Genotype,shape=Site,group=Genotype)) +
  geom_point()+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1*3, y = 0, 
                                                    yend = Dim.2*3),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  stat_ellipse( aes(color = Genotype), 
                show.legend = FALSE, 
                level = 0.95)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*3.5, y =  Dim.2*3.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (34.82% variance)",y="PCA2 (19.58% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right")
Pca1

colors8<-c("#BC4B07","#6E4465","#c70223","#f09ed6",
            "#D3BC38",  "#328F37",
           "#6FCAF7", "#1623cc")


colors18<-c("#BC4B07","#6E4465","#db6c04","#c70223","#732D3D","#E3423F","#f09ed6",
            "#D3BC38",  "#328F37","#a3ad32","#46bf3d","#e3d12b","#234f20",
            "#1e73b0","#3ca3bd", "#3453A8" ,"#1623cc","#6FCAF7")
Pca2<-ggplot(data = data, aes(x = pc1, y = pc2,color=Genotype,shape=Time,group=Time)) +
  geom_point()+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1*3, y = 0, 
                                                    yend = Dim.2*3),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  scale_color_manual(values=colors8)+
  stat_ellipse( aes(group = Time), 
                show.legend = FALSE, 
                level = 0.95)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*3.5, y =  Dim.2*3.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (34.82% variance)",y="PCA2 (19.58% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right")
Pca2
Pca<-ggarrange(Pca1,Pca2)
Pca
ggsave("Output/PCA.pdf",Pca,width=10,height=4)
ggsave("Output/PCA.svg",Pca,width=10,height=4)
ggsave("Output/PCA.jpg",Pca2,width=5,height=4)

##PCA by time
model<-data %>% filter(time=="T0")
pcaT0<-PCA(model[,c(6:8,10:13)])
pca.varsT0 <- pcaT0$var$coord %>% data.frame
pca.varsT0$vars <- rownames(pca.varsT0)
pca.vars.mT0 <- melt(pca.varsT0, id.vars = "vars")
model$pc1 <- pcaT0$ind$coord[, 1]
model$pc2 <- pcaT0$ind$coord[, 2]
PcaT0<-ggplot(data = model, aes(x = pc1, y = pc2,color=Genotype)) +
  geom_point()+
  geom_segment(data = pca.varsT0, inherit.aes =F, aes(x = 0, xend = Dim.1*3, y = 0, 
                                                    yend = Dim.2*3),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  scale_color_manual(values=colors8)+
  geom_text(data = pca.varsT0, inherit.aes =F,
            aes(x = Dim.1*3.5, y =  Dim.2*3.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (32.89% variance)",y="PCA2 (27.38% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5))+ggtitle("T0")+ylim(-4.5,4.5)+xlim(-4.5,4.5)
PcaT0

modelT12<-data %>% filter(time=="T12")
pcaT12<-PCA(modelT12[,c(6:13)])
pca.varsT12 <- pcaT12$var$coord %>% data.frame
pca.varsT12$vars <- rownames(pca.varsT12)
pca.vars.mT12 <- melt(pca.varsT12, id.vars = "vars")
modelT12$pc1 <- pcaT12$ind$coord[, 1]
modelT12$pc2 <- pcaT12$ind$coord[, 2]
PcaT12<-ggplot(data = modelT12, aes(x = pc1, y = pc2,color=Genotype)) +
  geom_point()+
  geom_segment(data = pca.varsT12, inherit.aes =F, aes(x = 0, xend = Dim.1*3, y = 0, 
                                                    yend = Dim.2*3),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  scale_color_manual(values=colors8)+
  geom_text(data = pca.varsT12, inherit.aes =F,
            aes(x = Dim.1*3.5, y =  Dim.2*3.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (35.76% variance)",y="PCA2 (28.41% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5))+ggtitle("T12")+ylim(-4.5,4.5)+xlim(-4.5,4.5)
PcaT12

modelT27<-data %>% filter(time=="T27")
pcaT27<-PCA(modelT27[,c(7:12)])
pca.varsT27 <- pcaT27$var$coord %>% data.frame
pca.varsT27$vars <- rownames(pca.varsT27)
pca.vars.mT27 <- melt(pca.varsT27, id.vars = "vars")
modelT27$pc1 <- pcaT27$ind$coord[, 1]
modelT27$pc2 <- pcaT27$ind$coord[, 2]
PcaT27<-ggplot(data = modelT27, aes(x = pc1, y = pc2,color=Genotype)) +
  geom_point()+
  geom_segment(data = pca.varsT27, inherit.aes =F, aes(x = 0, xend = Dim.1*3, y = 0, 
                                                    yend = Dim.2*3),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  scale_color_manual(values=colors8)+
  geom_text(data = pca.varsT27, inherit.aes =F,
            aes(x = Dim.1*3.5, y =  Dim.2*3.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (44.61% variance)",y="PCA2 (28.93% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right",
                        plot.title = element_text(hjust = 0.5))+ggtitle("T27")+ylim(-4.5,4.5)+xlim(-4.5,4.5)
PcaT27
Pca<-ggarrange(PcaT0,PcaT12,PcaT27)
Pca
ggsave("Output/PCAbytime.jpg",Pca,width=7,height=7)

datafin<-data %>% filter(colony %in% data$colony[data$time=="T27"])
pca1<-PCA(datafin[,7:12])
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
datafin$pc1 <- pca1$ind$coord[, 1]
datafin$pc2 <- pca1$ind$coord[, 2]

Pca<-ggplot(data = datafin, aes(x = pc1, y = pc2,color=genotype,shape=site)) +
  geom_point()+
  facet_grid(~time)+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1, y = 0, 
                                                    yend = Dim.2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*1.5, y =  Dim.2*1.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (38.95% variance)",y="PCA2 (22.17% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right")
Pca
ggsave("Output/PCA2.jpg",Pca,width=8,height=4)
ggsave("Output/PCA2.svg",Pca,width=8,height=4)

#evolution of symbiont genera
genera<-ITS2_4[,c(3,4,5,31,32,33)]
colnames(genera)[4]<-"Symbiodinium"
colnames(genera)[5]<-"Cladocopium"
colnames(genera)[6]<-"Durusdinium"
genera$treatment<-"Transplant"
genera$treatment[genera$site=="Nursery"]<-"Control"
model<-genera %>% filter(time=="T12" & genotype=="10")
kruskal.test(model$Cladocopium,model$treatment)
kruskal.test(model$Durusdinium,model$treatment)
model<-genera %>% filter(time=="T27")
kruskal.test(model$Symbiodinium,model$site)
model<-genera %>% filter(time!="T0")
kruskal.test(model$Cladocopium,model$site)
model<-genera %>% filter(genotype=="3")
pairwise.wilcox.test(model$Durusdinium,model$time)
pairwise.wilcox.test(model$Symbiodinium,model$time)
pairwise.wilcox.test(model$Cladocopium,model$time)

genera2<-gather(genera,Genus,percent,c(4:6))
meangenus<-ddply(genera2,c("genotype","time","Genus"),summarize,perc=mean(percent))

genera3<-ddply(genera2,c("time","Genus"),summarize,perc=mean(percent),x=sem(percent))
genera3$Genus<-factor(genera3$Genus,c("Symbiodinium", "Cladocopium" , "Durusdinium"))
colors3<-c("#BC4B07","#a3ad32","#1e73b0")                 
genus<-ggplot(genera3, aes(x = time, y = perc, col=Genus))+
  geom_point(size = 4)+
  geom_line(lwd=0.9)+
  ylab("Relative abundance %")+
  geom_errorbar(aes(ymin=perc-x, ymax=perc+x),width = 0.1)+
  scale_color_manual(values = colors3)+
  theme_bw()+
  theme(legend.text = element_text(face = "italic"))

ggsave("Output/genera_evolution.pdf",genus,width=5,height=4)
ggsave("Output/genera_evolution.svg",genus,width=5,height=4)

