# Clonality 2017 samples ########

#ask Ford:
#why is barplot not adding to 1 but boxplot is?
            # how to check for shuffling/ shifting? Create new column, if changed from CD TO D, or CD to C then shuffling, 
# if changed from D to C  vs then shifting

#Josh: Kriging prop D, bleaching


################################ EDGER ##################################################################


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()


BiocManager::install(c("edgeR", "limma"))

absolute_counttype<-read_excel("142_20210309_DBV_20210309T054204.seqs.absolute.abund_and_meta.xlsx") %>%
  filter (post_med_absolute>0) %>%
  select(c(2, 38:320))%>%
  rename(ID=sample_name)

#need package data.table
trans<-dcast(melt(absolute_counttype, id.vars = "ID"), variable ~ ID)

row.names(trans) <- trans$variable  #to make first column as name of rows
trans1<-trans%>%
  select(-c(1))


#adding group:

  groupC<-c("C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C",	"C")

groupD<-c("D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D",	"D")
group<-c("B",groupC,groupD,"F")

###

y <- DGEList(counts=trans1)

y <- DGEList(counts=trans1, group=group)

z <- calcNormFactors(y)

write.table(z[["counts"]], file="testofnormalizedfiles.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(z[["samples"]], file="testofnormalizedsamples.txt", row.names = TRUE, col.names = TRUE, sep = "\t" )

#open files in text and convert to xlsx and fix first line

#now need to calculate the read number of reads per sample/symbiont

numbereadssymbiont<-read_excel("testofnormalizedfiles.xlsx")

transposed<-dcast(melt(numbereadssymbiont, id.vars = "symbiont"), variable ~ symbiont)%>%
  rename(ID="variable") ##not working

numbereadfactor<-read_excel("testofnormalizedsamples.xlsx")


numbereadfactortotalnumber<-numbereadfactor%>%
  mutate(realreadnumber=lib.size*norm.factors)#to know read number

colnames(numbereadfactor) <- c("ID","group","lib.size", "norm.factors") #to change the column names
  
numbereadfactor$ID <- factor(numbereadfactor$ID)

#merging both in a single table

mergedtable<-left_join(numbereadfactor,transposed, by= "ID")
  
####now trying to multiply each value for the corresponding factor

row.names(transposed) <- transposed$ID  #to make first column as name of rows
transposed1<-transposed%>%
  select(-c(1))

f<-numbereadfactor%>%
  select(1,4)

row.names(f) <- f$ID  #to make first column as name of rows
f2<-f %>%
  select(c(2))


new_sweep<-sweep(transposed1, MARGIN=1,f2$norm.factors, "*")
ac<-16*2.098093e-01
bb<-2019*2.265828e-01  ##just checking that the values are actually correct, it seems they are!

#now this is the correct table with read number for the types

library(dplyr)
absolute_type_normalized <- tibble::rownames_to_column(new_sweep, "ID") #using this as the absolute type #


#checking how edge R changes the figures of read number per symb types
trial<-absolute_type_normalized%>%
mutate(Totalrow=rowSums (.[2:283] )) 

trial$sumtypes <- rowSums(trial > 0) 


library("ggpubr")
ggscatter(trial, y="Totalrow", x="sumtypes", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "number of symb types", ylab = "read number") #still some correlation!





########################################## Loading packages####################################################################################

dependencies <- c("tidyverse",
                  "plotly",
                  "here",
                  "readxl",
                  "janitor",
                  "vegan", 
                  "ggplot2", 
                  "dplyr"
)

#check if packages are installed - load if so, install+load if not)
for (i in dependencies) {
  if (i %in% row.names(installed.packages())){
    eval(bquote(library(.(i))))
    message(paste("loaded package",i))
  } else {
    install.packages(i)
    eval(bquote(library(.(i))))
    message(paste("installed and loaded package",i))
  }
}

setwd("~/Desktop/Clonality")

#Processeddata<-readRDS("Processeddata") #it already had 2019, but bleaching scores changed, no dont use it for bleaching score




########################################################### CALCULATION READ NUMBERS####################
setwd("~/Desktop/Clonality")


absol_sum2019<-read_excel("142_20210309_DBV_20210309T054204.seqs.absolute.abund_and_meta.xlsx") %>%
  filter (post_med_absolute>0) %>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  select(-c(1,6:9, 11:18,21:40))%>%
  filter(year=="2019")%>%
  mutate(Totalrow=rowSums (.[8:290] )) %>%
  mutate(meantotal=mean(Totalrow))%>%
  mutate(alllog=log10(Totalrow))%>%
  mutate(meanlog=mean(alllog))%>%
  mutate(stdlog=sd(alllog))%>%
  mutate("2std"=meanlog+(2*stdlog))%>%
  mutate("-2std"=meanlog-(2*stdlog))%>%
  filter(Totalrow<=26290.09)%>%
  filter(Totalrow>=58.64595)%>%
  select(sample_name)

#496 samples
backtransform2=10^4.419792; backtransform2 #26290.09
backtransformnegative2=10^1.768238; backtransformnegative2 # 58.64595


#for 2018, based on samples from chapter 1
#making the full table tidy, excluding samples 2018 and 2019 then merging 

absol_sum2018<-read_excel("142_20210309_DBV_20210309T054204.seqs.absolute.abund_and_meta.xlsx") %>%
  filter (post_med_absolute>0) %>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  select(-c(1,6:9, 11:18,21:40))%>%
  filter(year=="2017")%>%
  mutate(Totalrow=rowSums (.[8:290] )) %>%
  filter(Totalrow>=250)%>% #standardizing number of reads to 2 std above and below mean
  filter(Totalrow<=16560) %>%
  select(sample_name)
  
mergedabsolute<-bind_rows(absol_sum2018, absol_sum2019)

##not used__
onlyID<-absol_sum2019%>%
  select(ID) #now 472 samples

onlysymbiont<-absol_sum2019%>%
  select(1,5:287)
onlysymbiont$sumtypes <- rowSums(onlysymbiont > 0)
totalreadandtypes<- onlysymbiont%>%
  mutate(Totalrow=rowSums (.[2:284] ))

#visualization read number vs types of symbionts


library("ggpubr")
ggscatter(totalreadandtypes, y="Totalrow", x="sumtypes", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "number of symb types", ylab = "read number")
##_____________________


#################################### uploading dataset ####################################
#relative
Symportalrel<-read_excel("142_20210309_DBV_20210309T054204.seqs.relative.abund_and_meta.xlsx") %>%
  filter (post_med_absolute>0)%>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  select(-c(1,6:40))  %>%
  left_join(mergedabsolute,., by="sample_name")%>%
  select(-c(1,3))

#absolute
Symportalab<-read_excel("142_20210309_DBV_20210309T054204.seqs.absolute.abund_and_meta.xlsx") %>%
 filter (post_med_absolute>0) %>%
 separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
 separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  select(-c(1,6:40)) %>%
  left_join(mergedabsolute,., by="sample_name")%>%
  select(-c(1,3))


#Symportalab<-absolute_type_normalized%>%
  #rename(sample_name=ID)%>%
  #separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  #separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  #select(-c(1,3)) #using absolute_type_normalized instead


#profile relative
Symportalpro<-read_excel("142_20210309_DBV_20210309T054204.profiles.relative.abund_and_meta.xlsx",skip=6)%>%
  dplyr:: rename(sample_name = 2)%>%
  select( -c(1))%>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  left_join(mergedabsolute,., by="sample_name")%>%
  select(-c(1,3))

Symportalproab<-read_excel("142_20210309_DBV_20210309T054204.profiles.absolute.abund_and_meta.xlsx",skip=6)%>%
  dplyr::rename(sample_name= 2)%>%
  select( -c(1))%>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  left_join(mergedabsolute,., by="sample_name")%>%
  select(-c(1,3))


#Metadata

Meta<- read_excel("Coral_Sample_Metadata.xlsx")%>% #for 2017
  select(-c(1,4: 9, 13, 16:21))%>%
  mutate(ID= as.character(ID))

depthmeta<- read_excel("Survey_Sites.xlsx")   %>%
  inner_join(., Meta, by = "site")%>%
  mutate(ID= as.character(ID), block=as.character(block))

meta19<- read_excel("MontiporaBleachingScoring.xlsx", sheet= "Aggregated Bleaching Scores" )%>%
  select(c(1,7)) %>%
  rename(bleachscore=Mean) %>%
  rename(ID=COLONY)%>%
  mutate(year=2019)%>%
  mutate(ID= as.character(ID), year=as.character(year))


###making tables tidy


Relative<-Symportalrel %>%
  pivot_longer(c("B1": "F3.1"), names_to = "symbtypes", values_to = "relative")%>%
  mutate(relative= as.numeric(relative)) %>%
  filter(relative != 0) 


Absolute<-Symportalab %>%
  pivot_longer(c("B1": "F3.1"), names_to = "symbtypes", values_to = "absolute")%>%
  mutate(absolute= as.numeric(absolute)) %>%
  filter(absolute != 0)  

Profile<-  Symportalpro %>%
  pivot_longer(c(3:32), names_to = "profiles", values_to = "relativepro")%>%
  mutate(relativepro= as.numeric(relativepro)) %>%
  filter(relativepro != 0)  

ProfileAb<-Symportalproab%>%
  select(c(1:32))%>%
  pivot_longer(c(3:32), names_to = "Profiles", values_to = "absoluteepro")%>%
  mutate(absoluteepro= as.numeric(absoluteepro)) %>%
  filter(absoluteepro != 0)  


#Merging all tables into a single table, with all information and saving it as RDS

Processeddata<-inner_join(Relative,Absolute,by=c("ID","year","symbtypes"))%>%
  left_join(.,Profile,by=c("ID","year"))%>%
  left_join(.,depthmeta,by=c("ID"))%>%
  left_join(.,meta19,by=c("ID","year"))%>%
  left_join(.,ProfileAb, by=c("ID", "year"))%>%
  mutate(block=as.factor(block),ID=as.factor(ID),site=as.factor(site))%>%
  mutate(clade = case_when(
    str_detect(symbtypes, "A") ~ "A",
    str_detect(symbtypes, "B") ~ "B",
    str_detect(symbtypes, "C") ~ "C",
    str_detect(symbtypes, "D") ~ "D",
    str_detect(symbtypes, "E") ~ "E",
    str_detect(symbtypes, "F") ~ "F", 
    TRUE ~ NA_character_))%>%
  filter(clade=="C"| clade=="D") 


str(Processeddata)
saveRDS(Processeddata,"Processeddata")


###############################################checking numbers-not used####################################

Symportalab<-read_excel("142_20210309_DBV_20210309T054204.seqs.absolute.abund_and_meta.xlsx") %>%
  filter (post_med_absolute>0) %>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  select(-c(1,2,4,6:40))%>%
  filter(year==2019)



Symportalproab<-read_excel("142_20210309_DBV_20210309T054204.profiles.absolute.abund_and_meta.xlsx",skip=6)%>%
  rename(sample_name= 2)%>%
  select( -c(1))%>%
  separate(sample_name,into=c("year","ID"),sep="_ID",remove=FALSE)%>%
  separate(year, into=c("Clo", "year"),sep="o",remove=FALSE) %>%
  select(-c(1,3))    %>%
  filter(year==2019)  
      
m2017<-right_join(Meta, Symportalab, by="ID")%>%
  select(B1:F3.1)

m2019<-right_join(Meta, Symportalab, by="ID")%>%
  select(B1:F3.1)


 
m2019



sum(m$B1)
b<-m2019[, colSums(m2019 != 0) > 0]

f<-Symportalproab[, colSums(Symportalproab != 0) > 0]




################################## proportion of symbionts ##################################################################################

library(tidyverse)
Processeddata<-readRDS("Processeddata")

Dataall<-Processeddata %>%
  group_by(ID,year)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, year, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads) %>% 
  spread(clade,proportion)%>%
  select( ID, year,C,D, site,symbtypes, relative, absolute, bleachscore)%>%
  distinct()%>%
  group_by(ID,year)%>%
  fill(c(C,D),.direction="up")%>%
  drop_na(C)%>%
  mutate(D=case_when(C==1~0,TRUE~as.numeric(D)))%>%
  dplyr::rename (prop_c=C,prop_d=D)

table(Dataall$year)


#############################################josh test #############################################################
#Dataall%>%filter(year==2019)%>%select(ID)%>%distinct()

# JOSH TEST

tempmeta<-readRDS("tempmeta")

temp <- aggregate(absolute ~ ID + year + clade, Processeddata, sum)
tempC <- temp[temp$clade=="C",]
names(tempC) <- c("ID", "year", "C", "absC")
tempD <- temp[temp$clade=="D",]
names(tempD) <- c("ID", "year", "D", "absD")

temp2 <- merge(tempC, tempD)
pp <- Processeddata[c("ID", "year", "bleachscore", "block", "site")]
pp <- pp[!duplicated(pp),]
pp$site <- as.character(pp$site)
temp3 <- merge(temp2, pp, all.x=TRUE)
temp3 <- temp3[temp3$year==2019,]


temp3 <- merge(temp3, tempmerged)
pp <- tempmeta[c("site", "depth_m")]
pp <- pp[!duplicated(pp),]
temp3 <- merge(temp3, pp, all.x=TRUE)
write.csv(temp3, "josh_data.csv", row.names = FALSE)

temp3$prop_D <- temp3$absD / (temp3$absD + temp3$absC)
temp3$bleachscore_norm <- temp3$bleachscore / 3
temp3$block <- factor(temp3$block)

quartz()
plot(jitter(prop_D, 400) ~ jitter(bleachscore, 4), temp3, col="grey")


mod1 <- glm(cbind(absD, absC) ~ bleachscore + block, temp3, family="binomial")
mod2 <- glm(cbind(absD, absC) ~ T_DHW + block, temp3, family="binomial")
summary(mod2)



panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
pairs(temp3[c("T_DHW", "T_mean", "S_mean", "depth_m")], upper.panel=panel.cor)

pairs(temp3[c("T_DHW", "T_mean", "T_min", "T_max", "T_range", "T_drange", "T_dsd", "S_mean", "S_min", "S_max", "S_range", "Std", "depth_m")], upper.panel=panel.cor)

library(hier.part)
hier.part(temp3$prop_D, temp3[c("T_DHW", "T_mean", "T_max", "T_range", "T_drange", "S_mean", "S_min", "Std", "depth_m")])
mod2 <- glm(cbind(absD, absC) ~ T_mean + T_range + S_mean + depth_m + block, temp3, family="binomial")
drop1(mod2, test="Chisq")
pairs(temp3[c("T_mean", "T_range", "S_mean", "depth_m")], upper.panel=panel.cor)
hier.part(temp3$prop_D, temp3[c("T_mean", "T_range", "S_mean", "depth_m", "bleachscore")])
hier.part(temp3$bleachscore, temp3[c("T_mean", "T_range", "S_mean", "depth_m", "prop_D")])

#depth explains only 2 % variation
bs <- seq(0, 3, 0.1)
pred1 <- predict(mod1, list(bleachscore=bs, block=rep("1", length(bs))), type="response", se.fit=TRUE)
lines(bs, pred1$fit, col="#4B0055")
pred5 <- predict(mod1, list(bleachscore=bs, block=rep("5", length(bs))), type="response", se.fit=TRUE)
lines(bs, pred5$fit,col="#00588B")
pred2 <- predict(mod1, list(bleachscore=bs, block=rep("2", length(bs))), type="response", se.fit=TRUE)
lines(bs, pred2$fit, col="#009B95")
pred3 <- predict(mod1, list(bleachscore=bs, block=rep("3", length(bs))), type="response", se.fit=TRUE)
lines(bs, pred3$fit, col="#53CC67")
pred4 <- predict(mod1, list(bleachscore=bs, block=rep("4", length(bs))), type="response", se.fit=TRUE)
lines(bs, pred4$fit, col="#FDE333")

#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"


# other way around

quartz()
b<-plot(jitter(bleachscore_norm, 4) ~ jitter(prop_D, 500), temp3, col="grey");b

mod <- glm(bleachscore_norm ~ prop_D + block, temp3, family="binomial")
summary(mod)

pd <- seq(0, 1, 0.05)
pred1 <- predict(mod, list(prop_D=pd, block=rep("1", length(pd))), type="response", se.fit=TRUE)
lines(pd, pred1$fit, col="black")
pred5 <- predict(mod, list(prop_D=pd, block=rep("5", length(pd))), type="response", se.fit=TRUE)
lines(pd, pred5$fit,col="red")
pred2 <- predict(mod, list(prop_D=pd, block=rep("2", length(pd))), type="response", se.fit=TRUE)
lines(pd, pred2$fit, col="blue")
pred3 <- predict(mod, list(prop_D=pd, block=rep("3", length(pd))), type="response", se.fit=TRUE)
lines(pd, pred3$fit, col="purple")
pred4 <- predict(mod, list(prop_D=pd, block=rep("4", length(pd))), type="response", se.fit=TRUE)
lines(pd, pred4$fit, col="green")





#polygon(c(bs, rev(bs)), c(predA$fit + predA$se.fit, rev(predA$fit - predA$se.fit)), border=rgb(1,0,0,0.3), col=rgb(1,0,0,0.3))

#############################################boxplot proportion of D per site ########################################

Processeddata<-readRDS("Processeddata")

Dataall<-Processeddata %>%
  group_by(ID,year)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, year, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  group_by(ID, year, clade) %>% 
  mutate(proportion=cladereads/totalreads) %>% 
  group_by(ID, year) %>% 
  spread(clade,proportion)%>%
  select( ID, year,C,D, site,symbtypes, relative, absolute, bleachscore)%>%
  distinct()%>%
  group_by(ID,year)%>%
  fill(c(C,D),.direction="up")%>%
  drop_na(C)%>%
  mutate(D=case_when(C==1~0,TRUE~as.numeric(D)))%>%
  rename (prop_c=C,prop_d=D)



quartz()
PropD = ggplot(Dataall) + 
  theme_classic() +
  geom_boxplot(aes(x = site, y = prop_d)) + 
  facet_wrap (~year, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 7), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 10));PropD


checkD<-Processeddata%>%
  filter(ID=="256")


a<- Dataall%>%
  spread(year, prop_d)%>%
  clean_names() %>%
  select(id,site,x2017,x2019)%>%
  distinct()%>%
  group_by(id)%>%
  fill(c(x2017,x2019),.direction="up")%>%
  drop_na()%>%
  separate(site,into=c("block","trash"),sep="_",remove=FALSE)%>%select(-trash)%>%
  mutate(diff= x2019-x2017)
  

str(a)
saveRDS(a,"a")

a<-readRDS("a")
  
#scatterplot propd D 2017, propD 2019
quartz()
k<-ggplot(a, aes(x=x2017, y=x2019, color=block))+
  scale_fill_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
  scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
  geom_point(alpha=0.7)+
  theme_classic()+
  geom_abline(slope=1)+
  xlab("Proportion of D in 2018")+
  ylab("Proportion of D in 2019")+
  #scale_fill_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270"))+
  #scale_color_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270"))+
  theme(axis.text.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 9), 
        axis.title.x = element_text(size = 9),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)); k




#boxplot difference in prop per block
quartz(w=5,h=4)
m<-ggplot(a,aes(site,diff))+
  geom_boxplot(alpha=0.3)+
  #geom_point(aes(site,diff))+
  theme_classic()+
  facet_wrap (~block, ncol=5, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, size = 7.5, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = "Difference % D")+
    theme(legend.position="none") ;  m +theme(strip.background = element_rect(fill="#4B0055"))

q4<- sequential_hcl(5, palette="Viridis"); q4


m<-ggplot(a,aes(site,diff, fill=block))+
  geom_boxplot(alpha=0.3)+
  #geom_point(aes(site,diff))+
  theme_classic()+
  facet_wrap (~block, ncol=5, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, size = 7.5, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = "Difference % D")+
  scale_fill_viridis_d(name= "block") +
  theme(legend.position="none") ;  m

  #scale_fill_manual(values = colours1)


  n<-ggplot(a, aes(block,diff, fill=block))+
  geom_boxplot(alpha=0.3)+
  #geom_point(aes(site,diff))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = "Difference % D")+
    scale_fill_viridis_d(name="Block")+
  theme(legend.position="none") ; n
#scale_fill_manual(values = colours1)


#for permanova


permanovaD<-na.omit(a)%>%
  ungroup()%>%
  select(site, block, diff)%>%
  distinct()%>%
  mutate(diffpositive=diff+1)


numericaD<-permanovaD%>%
  ungroup()%>%
  select(diffpositive)


siteblockD<-permanovaD%>%
  ungroup()%>%
  select(1:2)%>%#p.env contains only the environmental (independant) variables 
  mutate(site=as.factor(site))%>%
  mutate(block=as.factor(block)) 

#running the permanova test

#for all sites in block (general)
permanovaDall <- adonis(numericaD ~ block, data=siteblockD, permutations=999, method="bray"); permanovaDall

#for each pairwises of blocks
airPermanovaD<- pairwise.adonis2(numericaD ~ block,data=siteblockD, permutations=999, method="bray"); Permanova 

#not signigificant

#testing if proportion of D is different among both years

b<-a%>%
  ungroup()%>%
  select("x2017", "x2019")

write.csv(b, "b.csv")

# t test one tailed t test -2.92593. The p-value is .001773. The result is significant at p < .05.
#t test one tailed The t-value is 2.92593. The p-value is .003546. The result is significant at p < .05.
  


#not used
#quartz()
#one_row <- plot_grid(m,NULL, k, labels = c('A',NULL,'B'),align="h",axis="tb",nrow = 1, label_size = 12,rel_widths = c(1,0.05, 0.6)); one_row
#quartz()
#one_row <- plot_grid(m,NULL, n, labels = c('A',NULL,'B'),align="h",axis="tb",nrow = 1, label_size = 12,rel_widths = c(1,0.05, 0.6)); one_row

#trying another cowplot option

cow_n_bleaching <- plot_grid(NULL, n, NULL, bleaching,
                    # A negative rel_height shrinks space between elements
                    rel_heights = c(0.1, 1, 0.15, 1),
                    ncol = 1,
                    label_y = 1.07,
                    labels = c("","B","","C"),
                    label_size = 12);cow_n_bleaching


cow_m <- plot_grid(NULL, m,
                   rel_heights = c(0.1, 1.85),
                   ncol = 1,
                   label_y = 1 + 0.07 /1.85,
                   labels = c("","A"),
                   label_size = 12)

quartz()
cow_mn <- plot_grid(cow_m,cow_n_bleaching,
                    rel_widths = c(3,2),
                    nrow = 1); cow_mn

#adding coral figure

icon_coral <- "coral.png"
icon_bleached<-"bleached.png"

quartz()
cow_final <- ggdraw() +
  draw_plot(cow_mn) +
  draw_image(image = icon_coral, x = 0.155, y = -0.1, scale = 0.080)+
  draw_image(image = icon_bleached, x = 0.166, y = -0.34, scale = 0.080); cow_final



### % of D in 2018 and 2019

p<-ggplot(a, aes(block,x2017, fill=block))+
  geom_boxplot(alpha=0.3)+
  #geom_point(aes(site,diff))+
  theme_classic()+
  theme(axis.title.x=element_blank(), 
              axis.text.x = element_blank(),
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = " %D 2018")+
  scale_fill_viridis_d(name="Block")+
  theme(legend.position="none") ; p


p9<-ggplot(a, aes(block,x2019, fill=block))+
  geom_boxplot(alpha=0.3)+
  #geom_point(aes(site,diff))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = "% D 2019")+
  scale_fill_viridis_d(name="Block")+
  theme(legend.position="none") ; p9


quartz()
cow_yw <- plot_grid(NULL, p, p9,
                    # A negative rel_height shrinks space between elements
                    rel_heights = c(0.1, 1, 1.2),
                    ncol = 1,
                    label_y = 1.07,
                    labels = c("","A","B"),
                    label_size = 12);cow_yw

cow_p <- plot_grid(NULL, bl,
                   rel_heights = c(0.1, 1.85),
                   ncol = 1,
                   label_y = 1 + 0.07 /1.85,
                   labels = c("","C"),
                   label_size = 12)

quartz()
cow_ppl <- plot_grid(cow_yw,cow_p,
                    rel_widths = c(2,3),
                    nrow = 1); cow_ppl


########################################## Bleaching per block #######################################


bleachblock<-Dataall%>%
  filter(year==2019)%>%
  separate(site,into=c("block","trash"),sep="_",remove=FALSE)%>%select(-trash)

bleaching<-ggplot(bleachblock, aes(block,bleachscore, fill=block))+
  geom_boxplot(alpha=0.2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = "Bleach score")+
  annotate(geom="text", x= 1, y=3.5, label="*2",color="gray70", size=3.2)+   
  annotate(geom="text", x= 2, y=3.5, label="*1,4,5",color="gray70", size=3,2)+
  annotate(geom="text", x= 3, y=3.5, label="*4,5",color="gray70", size=3.2)+
  annotate(geom="text", x= 4, y=3.5, label="*2",color="gray70", size=3.2)+
annotate(geom="text", x= 5, y=3.5, label="*2,3",color="gray70", size=3.2)+
scale_fill_viridis_d(name="block")+
  theme(legend.position="none");bleaching

quartz()
bl<-ggplot(bleachblock,aes(site,bleachscore, fill=block))+
  geom_boxplot(alpha=0.3)+
  #geom_point(aes(site,diff))+
  theme_classic()+
  facet_wrap (~block, ncol=5, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, size = 7.5, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 9)) +
  labs(x = "", y = "Bleaching score")+
  scale_fill_viridis_d(name= "block") +
  theme(legend.position="none") ;  bl




#doing a permanoa to test significance difference

nona<-na.omit(bleachblock)%>%
  ungroup()%>%
  select(ID, site, block, bleachscore)%>%
  distinct()

numericalble<-nona%>%
  ungroup()%>%
  select(bleachscore)
  

siteblockble<-nona%>%
  ungroup()%>%
  select(1:3)%>%#p.env contains only the environmental (independant) variables 
  mutate(site=as.factor(site))%>%
  mutate(block=as.factor(block)) 

#running the permanova test

#for all sites in block (general)
Permanovable <- adonis(numericalble ~ block, data=siteblockble, permutations=999, method="bray"); Permanovable

#for each pairwises of blocks
Permanovable <- pairwise.adonis2(numericalble ~ block,data=siteblockble, permutations=999, method="bray"); Permanovable 


###trying heat map bleaching per block


bleachgr<-read_excel("bleachingR2.xlsx")%>%
  select(c(2:6))

heatmap(as.matrix(bleachgr[, -1], Colv = NA, scale="column", scale="row"))









########################################## Major Symbiont type (relative) ############################

types<- Dataall%>%
  dplyr::select(symbtypes) %>% 
  distinct() 


Majorsymbtypesall <-Processeddata%>% 
  group_by(year)%>%
  mutate(typessummary = case_when(
    str_detect(symbtypes, "D1") ~ "D1",
    str_detect(symbtypes, "D4") ~ "D4",
    str_detect(symbtypes, "D6") ~ "D6",
    str_detect(symbtypes, "C17") ~ "C17",
    str_detect(symbtypes, "D3") ~ "D3", 
    str_detect(symbtypes, "_D") ~ "D", 
    str_detect(symbtypes, "C31") ~ "C31", 
    str_detect(symbtypes, "C21") ~ "C21", 
    str_detect(symbtypes, "C15") ~ "C15", 
    str_detect(symbtypes, "_C") ~ "C", 
    str_detect(symbtypes, "C3") ~ "C3", 
    str_detect(symbtypes, "C1") ~ "C1", 
    str_detect(symbtypes, "D2") ~ "D2" ))%>% 
  group_by(year, site)%>% 
  mutate(totalreadstype=sum(absolute)) %>% 
  group_by(year, site, typessummary) %>% 
  mutate(summarytypes=sum(absolute)) %>% 
  mutate(proportionsymb=summarytypes/totalreadstype) %>% 
  spread(typessummary,proportionsymb)%>%
  mutate(D1 = replace_na(D1, 0)) %>%
  mutate(D4 = replace_na(D4, 0))   %>%  
  mutate(D6 = replace_na(D6, 0))   %>%
  mutate(C17 = replace_na(C17, 0))   %>%
  mutate(D3 = replace_na(D3, 0))   %>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C31 = replace_na(C31, 0))%>%
  mutate(C21 = replace_na(C21, 0))%>%
  mutate(C15 = replace_na(C15, 0))%>%
  mutate(C = replace_na(C, 0))%>%
  mutate(C3 = replace_na(C3, 0))%>%
  mutate(C1 = replace_na(C1, 0))%>%
  mutate(D2 = replace_na(D2, 0))%>%
  dplyr::select(D1,D4,D6,C17,D3,D,C31,C21,C15,C,C3,C1,D2, site, block, year) %>% 
  distinct() 



#making it tidy
#real tidysym
tidysym<-Majorsymbtypesall %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0)  

tidysym2017<- tidysym%>%
  filter(year== "2017")

tidysym2019<- tidysym%>%
  filter(year== "2019")


colours1=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461", "#E8956D", "#F0B384", "#F5D1A8")


quartz(height = 4)
mx4 = ggplot(tidysym2017, aes(x = site, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  #facet_wrap (~block, ncol=8, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm"))  + 
  scale_y_continuous(expand = c(0,0)) + 
  #ggtitle(label = "A") +
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)); mx4


quartz(height = 4)
mx5 = ggplot(tidysym2019, aes(x = site, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  #facet_wrap (~block, ncol=8, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm"))  + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 2)); mx5


#panel 2017, 2019

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('2017','', '2019'), align = "v" ,label_size = 12, ncol=1, rel_heights = c(1.1,0,1))

#____________________________trying side by side

str(tidysym)
saveRDS(tidysym,"tidysym")

test<- tidysym%>%
  filter(site=="1_1"| site=="1_2"| site=="1_3" | site=="2_1"| site=="2_2"| site=="2_3")


colours1=c("#26185F" ,"#005D9E" ,"#0095AF", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461", "#E8956D", "#F0B384", "#F5D1A8")


label<-c("1_10", "", "1_2", "","1_3","", "1_4","", "1_6","", "1_9","", "2_1","", "2_2","", "2_3","",  "2_4","", "2_7","", "2_8", "","3_2", "", "3_3", "", "3_4", "", "3_5", "", "3_6", "", "3_7", "" ,
         "4_11", "", "4_14", "", "4_4", "", "4_5", "", "4_6", "", "4_8","","5_1", "", "5_2", "", "5_3","",  "5_6", "", "5_7", "","5_8", "")

alpha=year

mx4 = ggplot(tidysym, aes(x = interaction(year, site),
                       fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity", width = 0.9, position= position_stack()) +
  geom_vline(xintercept = c(0.5,2.5,4.5,6.5,8.5,10.5, 12.5, 14.5, 16.5, 18.5,20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5, 36.5, 38.5, 40.5, 42.5, 44.5, 46.6, 48.5, 50.5, 52.5, 54.5, 56.5,58.5),color= "white", size=3)+
  #geom_vline(xintercept = c(0.5,2.5,4.5,6.5,8.5,10.5, 12.5, 14.5, 16.5, 18.5,20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5, 36.5, 38.5, 40.5, 42.5, 44.5, 46.6, 48.5, 50.5, 52.5, 54.5, 56.5,58.5),color= "white", size=5)+
  facet_wrap (~block, ncol=5, scales="free") +
  theme(axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm"))  + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(labels= label) +
  #ggtitle(label = "A") +
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
 # scale_alpha_manual(values=c(0.5,1))+
  #scale_pattern_manual(values = c(2019 = "stripe", 2017 = "none"))+
  guides(fill = guide_legend(nrow = 1)); mx4

###


######################################### trying pattern#############

pattern1=c("stripe", "none")

library(ggplot2)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
ggplot(data = tidysym, mapping=aes(x =interaction(year, site), fill = majortypes, y =proportion, pattern = year)) +
  geom_bar_pattern(position = position_dodge(preserve = "single"),
                   stat='identity'), 
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = colours1) +
  scale_pattern_manual(values = pattern1) +
  labs(x = "", y="Relative proportion", fill = "ITS2 types", pattern = "year") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))


remotes::install_github("coolbutuseless/ggpattern")



#####new trial

remotes::install_github("coolbutuseless/ggpattern")
quartz()

ggplot(data=mydata)+
  geom_col(aes(x=letters, y=amount, fill=letters)) +
  geom_col_pattern(
    aes(letters, amount, pattern_fill = split,fill=letters), 
    pattern = 'stripe',
    colour  = 'black'
  )+ 
  scale_fill_manual(values=colfill)

geom_bar(stat = "identity", width = 0.9, position= position_stack()) +


################### from stackoverflow

library(ggplot2)
library(ggpattern)
install.packages('ggplot2')
install.packages('scales')
install.packages('grid')
install.packages('glue')
install.packages('rlang')
install.packages('sf')
install.packages('png')

# install.packages("remotes")
remotes::install_github("trevorld/gridpattern")
remotes::install_github("coolbutuseless/ggpattern")

colors2<- c("#004533", "#006F4B", "#32985A" ,"#62BB63" ,"#A2D489", "#CFE8AE", "#EEF7CF", "#802A07", "#B03F00" ,"#DE5600", "#F27E00", "#FBA453", "#FFC693")

quartz()
ggplot(tidysym, aes(x = year, y = proportion, fill = majortypes, pattern = year)) +
  #geom_col_pattern(pattern_size = .25) +
  #scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors2) +
  facet_wrap(~site, nrow = 1, strip.position =  "bottom") +
  theme(
    strip.placement = "outside",
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 3),
    panel.spacing.x = unit(0, "pt"),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 3),
    axis.title.y = element_text(size = 7),
    legend.title = element_text(size = 7),
    axis.text.y = element_text(colour = "black", size = 7),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"), 
    legend.margin = margin(-0.4, 0, 0, 0, unit = "cm")) +
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(pattern = "none")),
         pattern = guide_legend(nrow = 1, override.aes = list(fill = NA)))


########################################## Major Symbiont type (profile) ####################################

Processeddata<-readRDS("Processeddata")

DataallPro<-Processeddata %>% 
  drop_na(Profiles)%>%
  group_by(year, site)%>% 
  mutate(totalreadstype=sum(absoluteepro))%>%
  group_by(year, site, Profiles) %>% 
  mutate(summarytypes=sum(absoluteepro)) %>% 
  mutate(proportionsymb=summarytypes/totalreadstype)%>%
  select(site,Profiles,proportionsymb, block, year)%>%
  distinct()


colours1<-c( "#14505C", "#175762" ,"#1B5F67", "#1F666C", "#236E71", "#287576", "#2D7D7B", "#328580", "#388C84", "#3E9488", "#459B8C", "#4CA38F", "#54AA92", "#5CB295", "#65B998", "#6EC09A", "#7AC79D", "#88CCA1", "#96D2A6", "#A3D8AB",
             "#AFDDB1", "#BBE2B7", "#C7E5BE", "#D6765D" ,"#F8DFC1", "#772C4B")


DataallPro2017<- DataallPro%>%
  filter(year== "2017")

DataallPro2019<- DataallPro%>%
  filter(year== "2019")

quartz(height = 6)
mx4 = ggplot(DataallPro2017, aes(x = site, fill = Profiles, y = proportionsymb)) + 
  geom_bar(stat = "identity") + 
  facet_wrap (~block, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7), legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm")) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "Profiles") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 10)); mx4




quartz(height = 6)
mx5 = ggplot(DataallPro2019, aes(x = site, fill = Profiles, y = proportionsymb)) + 
  geom_bar(stat = "identity") + 
  facet_wrap (~block, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7), legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm")) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "Profiles") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 10)); mx5

#side by side


colours2<-c("#26185F","#0B4151", "#064957" ,"#01505D" ,"#005862", "#006068", "#00686D", "#007171", "#007975", "#018179",
"#11897D", "#1F9180", "#2C9983", "#39A185" ,"#46A988", "#53B18A" ,"#60B98C", "#6EC08E" ,"#7CC890",
"#8ACF92" ,"#99D794" ,"#A7DE97", "#B6E59A", "#C5EC9E", "#D4F3A3", "#E8956D", "#F0B384", "#F5D1A8", "#9C3E5D", "#67223F")

label<-c("1_10", "", "1_2", "","1_3","", "1_4","", "1_6","", "1_9","", "2_1","", "2_2","", "2_3","",  "2_4","", "2_7","", "2_8", "","3_2", "", "3_3", "", "3_4", "", "3_5", "", "3_6", "", "3_7", "" ,
         "4_11", "", "4_14", "", "4_4", "", "4_5", "", "4_6", "", "4_8","","5_1", "", "5_2", "", "5_3","",  "5_6", "", "5_7", "","5_8", "")



######################################### panel #####

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,pro, labels = c('A','', 'B') ,align = "v" ,label_size = 12, ncol=1, rel_heights = c(1,0,1.3))



##### side by side## use this one!

label<-c("1_10", "", "1_2", "","1_3","", "1_4","", "1_6","", "1_9","", "2_1","", "2_2","", "2_3","",  "2_4","", "2_7","", "2_8", "","3_2", "", "3_3", "", "3_4", "", "3_5", "", "3_6", "", "3_7", "" ,
         "4_11", "", "4_14", "", "4_4", "", "4_5", "", "4_6", "", "4_8","","5_1", "", "5_2", "", "5_3","",  "5_6", "", "5_7", "","5_8", "")
quartz(height = 4)
pro = ggplot(DataallPro, aes(x = interaction(year, site), fill = Profiles, y =proportionsymb)) + 
  geom_bar(stat = "identity", width = 0.9, position= position_stack()) +
  geom_vline(xintercept = c(0.5,2.5,4.5,6.5,8.5,10.5, 12.5, 14.5, 16.5, 18.5,20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5, 36.5, 38.5, 40.5, 42.5, 44.5, 46.6, 48.5, 50.5, 52.5, 54.5, 56.5,58.5),color= "white", size=3)+
  #facet_wrap (~block, ncol=8, scales="free") +
  theme(axis.text.x = element_text(size = 7, colour = "black", vjust = 0.5, hjust = -0.0), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm"))  + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(labels= label) +
  #ggtitle(label = "A") +
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours2) +
  scale_alpha_manual(values=c(0.5,1))+
  #scale_pattern_manual(values = c(2019 = "stripe", 2017 = "none"))+
  guides(fill = guide_legend(nrow = 10)); pro



######################################## bar plot relative site ###############

types<- Dataall%>%
  group_by(ID)%>%
  ungroup()%>%
  dplyr::select(symbtypes) %>% 
  distinct() 


Majorsymbtypesall <-Dataall%>% 
  mutate(typessummary = case_when(
    str_detect(symbtypes, "D1") ~ "D1",
    str_detect(symbtypes, "D4") ~ "D4",
    str_detect(symbtypes, "D6") ~ "D6",
    str_detect(symbtypes, "C17") ~ "C17",
    str_detect(symbtypes, "D3") ~ "D3", 
    str_detect(symbtypes, "_D") ~ "D", 
    str_detect(symbtypes, "C31") ~ "C31", 
    str_detect(symbtypes, "C21") ~ "C21", 
    str_detect(symbtypes, "C15") ~ "C15", 
    str_detect(symbtypes, "_C") ~ "C", 
    str_detect(symbtypes, "C3") ~ "C3", 
    str_detect(symbtypes, "C1") ~ "C1", 
    str_detect(symbtypes, "D2") ~ "D2" ))%>% 
  group_by(year, ID)%>% 
  mutate(totalreadstype=sum(absolute)) %>% 
  group_by( year, ID, typessummary) %>% 
  mutate(summarytypes=sum(absolute)) %>% 
  mutate(proportionsymb=summarytypes/totalreadstype) %>% 
  spread(typessummary,proportionsymb)%>%
  mutate(D1 = replace_na(D1, 0)) %>%
  mutate(D4 = replace_na(D4, 0))   %>%  
  mutate(D6 = replace_na(D6, 0))   %>%
  mutate(C17 = replace_na(C17, 0))   %>%
  mutate(D3 = replace_na(D3, 0))   %>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C31 = replace_na(C31, 0))%>%
  mutate(C21 = replace_na(C21, 0))%>%
  mutate(C15 = replace_na(C15, 0))%>%
  mutate(C = replace_na(C, 0))%>%
  mutate(C3 = replace_na(C3, 0))%>%
  mutate(C1 = replace_na(C1, 0))%>%
  mutate(D2 = replace_na(D2, 0))%>%
  dplyr::select(D1,D4,D6,C17,D3,D,C31,C21,C15,C,C3,C1,D2, site, ID, block, year) %>% 
  distinct()  

#making it tidy

tidysym<-Majorsymbtypesall %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==1)

tidysym2017<- tidysym%>%
  filter(year== "2017")

tidysym2019<- tidysym%>%
  filter(year== "2019")

colours1=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")

colours2=c("#14505C" ,"#226B70" , "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")


"#14505C" C
"#226B70" C1
"#348781" C15
"#4CA38F"C17
"#6BBE99"C21
"#9AD4A8"C3
"#C7E5BE"C31
"#B13F63"D
"#CC5762"D1
"#DE7461"D2
"#E8956D"D3
"#F0B384"D4
"#F5D1A8"D6



tidysym2017$site_f = factor(tidysym2017$site, levels=c('1_2','1_3','1_4','1_6', '1_9', '1_10'))
tidysym2019$site_f = factor(tidysym2019$site, levels=c('1_2','1_3','1_4','1_6', '1_9', '1_10'))


quartz(height = 4)
mx4 = ggplot(tidysym2017, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site_f, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  guides(fill = guide_legend(nrow = 7)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)); mx4

quartz(height = 4)
mx5 = ggplot(tidysym2019, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site_f, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  guides(fill = guide_legend(nrow = 7)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours2) +
  guides(fill = guide_legend(nrow = 1)); mx5



#panel 2017, 2019

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('2017','', '2019'), align = "v" ,label_size = 12, ncol=1, rel_heights = c(1.2,0,1))




#########block 2

"#14505C" C
"#226B70" C1
"#348781" C15
"#4CA38F"C17
"#6BBE99"C21
"#9AD4A8"C3
"#C7E5BE"C31
"#B13F63"D
"#CC5762"D1
"#DE7461"D2
"#E8956D"D3
"#F0B384"D4
"#F5D1A8"D6

tidysym<-Majorsymbtypesall %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==2)

tidysym2017<- tidysym%>%
  filter(year== "2017")

tidysym2019<- tidysym%>%
  filter(year== "2019")

colours1=c("#14505C" ,"#226B70" , "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")


colours2=c("#14505C" ,"#226B70", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")

quartz(height = 4)
mx4 = ggplot(tidysym2017, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)) ; mx4

quartz(height = 4)
mx5 = ggplot(tidysym2019, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours2) +
  guides(fill = guide_legend(nrow = 1)) ; mx5



#panel 2017, 2019

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('2017','', '2019'), align = "v" ,label_size = 12, ncol=1, rel_heights = c(1.2,0,1))



#block 3

"#14505C" C
"#226B70" C1
"#348781" C15
"#4CA38F"C17
"#6BBE99"C21
"#9AD4A8"C3
"#C7E5BE"C31
"#B13F63"D
"#CC5762"D1
"#DE7461"D2
"#E8956D"D3
"#F0B384"D4
"#F5D1A8"D6

tidysym<-Majorsymbtypesall %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==3)

tidysym2017<- tidysym%>%
  filter(year== "2017")

tidysym2019<- tidysym%>%
  filter(year== "2019")



colours1=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")


quartz(height = 4)
mx4 = ggplot(tidysym2017, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)) ; mx4

quartz(height = 4)
mx5 = ggplot(tidysym2019, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)) ; mx5


#panel 2017, 2019

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('2017','', '2019'), align = "v" ,label_size = 12, ncol=1, rel_heights = c(1.2,0,1))



####block 4

"#14505C" C
"#226B70" C1
"#348781" C15
"#4CA38F"C17
"#6BBE99"C21
"#9AD4A8"C3
"#C7E5BE"C31
"#B13F63"D
"#CC5762"D1
"#DE7461"D2
"#E8956D"D3
"#F0B384"D4
"#F5D1A8"D6

tidysym<-Majorsymbtypesall %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==4)

tidysym2017<- tidysym%>%
  filter(year== "2017")

tidysym2019<- tidysym%>%
  filter(year== "2019")

colours1=c("#14505C" ,"#226B70" , "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762","#E8956D", "#F0B384", "#F5D1A8")


quartz(height = 4)
mx4 = ggplot(tidysym2017, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)) ; mx4

quartz(height = 4)
mx5 = ggplot(tidysym2019, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)) ; mx5


#panel 2017, 2019

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('2017','', '2019'), align = "v" ,label_size = 12, ncol=1, rel_heights = c(1.2,0,1))


####block 5

"#14505C" C
"#226B70" C1
"#348781" C15
"#4CA38F"C17
"#6BBE99"C21
"#9AD4A8"C3
"#C7E5BE"C31
"#B13F63"D
"#CC5762"D1
"#DE7461"D2
"#E8956D"D3
"#F0B384"D4
"#F5D1A8"D6

tidysym<-Majorsymbtypesall %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==5)

tidysym2017<- tidysym%>%
  filter(year== "2017")

tidysym2019<- tidysym%>%
  filter(year== "2019")


colours1=c("#14505C" ,"#226B70" , "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")

colours2=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")

quartz(height = 4)
mx4 = ggplot(tidysym2017, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)) ; mx4



quartz(height = 4)
mx5 = ggplot(tidysym2019, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=6, scales="free") +
  theme(axis.title.y = element_text(size = 7),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="top",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours2) +
  guides(fill = guide_legend(nrow = 1)) ; mx5



#panel 2017, 2019

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('2017','', '2019'), align = "v" ,label_size = 12, ncol=1, rel_heights = c(1,0,1.2))



########################################nMDS relative ######################################

wide<- Symportalrel%>%
  right_join(depthmeta,., by="ID")


numerical<-wide%>%
  select(11:293)%>%
  select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data


meta<-wide%>%
  select(1:10) %>%
  rename(Block=block)

out<-bind_cols(meta,numerical)%>%
  group_by(ID)

nmds_results <- metaMDS(comm = numerical,  # Define the community data 
                        distance = "bray",       # Specify a bray-curtis distance
                        try = 100) #Run 100 stress 0.06557  excelent fit



set.seed(216)
NMDS=metaMDS(numerical, distance='bray') #0.06024865 
nmdsall<- bind_cols(meta, NMDS$points[,1],NMDS$points[,2]) %>%
  rename(NMDS1=11, NMDS2=12)


q4<- sequential_hcl(7, palette="Blues")




a<- nmdsall%>%
  group_by(year, Block)

quartz()
g<-ggplot(nmdsall, aes(x=NMDS1, y=NMDS2, color=year)) + 
  stat_ellipse(aes(fill = year, color=year), geom = "polygon", alpha = 0.2) +
  scale_fill_manual(labels=c("2018", "2019"), values=c("#F8B83C","#5087C1"))+
  scale_color_manual(values=c("#F8B83C","#5087C1"), labels=c("2018", "2019"))+
  geom_point(alpha = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5), 
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        axis.text.x = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        # axis.line = element_line(colour = "black", size=0.5),
        legend.position ="right",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'))+
  geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
  geom_vline(xintercept = 0, linetype="dotted", color="grey") +
  #labs(color = "Year", labels=c("2018", "2019"), fill = "Year")  +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x= -0.5, y=-0.60, label="Stress=0.061",
           color="grey", size=3); g




  #nmds 2019 only
  
  
  wide<- Symportalrel%>%
    right_join(depthmeta,., by="ID")%>%
    filter(year=="2019")
  
  
  numerical<-wide%>%
    select(13:295)%>%
    select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data
  
  
  meta<-wide%>%
    select(1:12) %>%
    rename(Block=block)
  
  out<-bind_cols(meta,numerical)%>%
    group_by(ID)
  
  nmds_results <- metaMDS(comm = numerical,  # Define the community data 
                          distance = "bray",       # Specify a bray-curtis distance
                          try = 100) # stress=0.061
  
  
  
  set.seed(216)
  NMDS=metaMDS(numerical, distance='bray') 
  nmdsall9<- bind_cols(meta, NMDS$points[,1],NMDS$points[,2]) %>%
    rename(NMDS1=13, NMDS2=14)
  
  
  
  q4<- sequential_hcl(5, palette="Viridis"); q4
  
  
  quartz()
  year9<-ggplot(nmdsall9, aes(x=NMDS1, y=NMDS2, color=Block)) + 
    stat_ellipse(aes(fill = Block,color=Block), geom = "polygon", alpha = 0.2) +
    scale_fill_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
    scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
    geom_point(alpha = 0.5) +
    #geom_point(aes(shape=year))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5), 
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          legend.title = element_text(size = 7), 
          axis.text.y = element_text(colour = "black", size = 7),
          axis.text.x = element_text(colour = "black", size = 7),
          panel.border = element_blank(),
         # axis.line = element_line(colour = "black", size=0.5),
          legend.position ="NULL",
          legend.text = element_text(size = 7, colour = "black"),
          legend.key.size = unit(0.3, 'cm'))+
    geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
    geom_vline(xintercept = 0, linetype="dotted", color="grey") +
    labs(color = "Blocks", fill = "Blocks")  +
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x= -0.5, y=-0.60, label="Stress=0.061",
             color="grey", size=3)+
    annotate(geom="text", x=-0.8 , y=0.7, label="2019",
             color="black", size=3.5); year9
  
  
  
  #####2018 only
  

  
  wide<- Symportalrel%>%
    right_join(depthmeta,., by="ID")%>%
    filter(year=="2017")
  
  
  numerical<-wide%>%
    select(13:295)%>%
    select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data
  
  
  meta<-wide%>%
    select(1:12) %>%
    rename(Block=block)
  
  out<-bind_cols(meta,numerical)%>%
    group_by(ID)
  
  nmds_results <- metaMDS(comm = numerical,  # Define the community data 
                          distance = "bray",       # Specify a bray-curtis distance
                          try = 100) # stress=0.061
  
  
  
  set.seed(216)
  NMDS=metaMDS(numerical, distance='bray') 
  nmdsall8<- bind_cols(meta, NMDS$points[,1],NMDS$points[,2]) %>%
    rename(NMDS1=13, NMDS2=14)
  

  quartz()
  year8<-ggplot(nmdsall8, aes(x=NMDS1, y=NMDS2, color=Block)) + 
    stat_ellipse(aes(fill = Block,color=Block), geom = "polygon", alpha = 0.2) +
    scale_fill_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
    scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
    geom_point(alpha = 0.5) +
    #geom_point(aes(shape=year))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5), 
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          legend.title = element_text(size = 7), 
          axis.text.y = element_text(colour = "black", size = 7),
          axis.text.x = element_text(colour = "black", size = 7),
          panel.border = element_blank(),
          # axis.line = element_line(colour = "black", size=0.5),
          legend.position ="right",
          legend.text = element_text(size = 7, colour = "black"),
          legend.key.size = unit(0.3, 'cm'))+
    geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
    geom_vline(xintercept = 0, linetype="dotted", color="grey") +
    labs(color = "Blocks", fill = "Blocks")  +
    theme(plot.title = element_text(hjust = 0.5))+
    annotate(geom="text", x= -0.5, y=-0.50, label="Stress=0.062",
             color="grey23", size=3) 
    #annotate(geom="text", x=-0.8 , y=0.7, label="2018",
             #color="black", size=3.5); 
  year8
  
  
  
  
  
  
  
  
  #panel 
  quartz()
  top_row <- plot_grid(year8, year9, g, NULL, labels = c('A', 'B', 'C', ''),ncol = 2, label_size = 12, align = "h",rel_widths = c(1.15, 1, 1.15)); top_row
  
  one_row <- plot_grid(year8, year9, g, labels = c('A', 'B', 'C'),ncol = 3, label_size = 12, align = "h",rel_widths = c(1.15, 1, 1.15)); one_row
  quartz()
  plot_grid(top_row, g, labels = c('','C' ), label_size = 12, ncol = 1, rel_widths = c(2, 1))
  
  plot_grid(b, a, g, NULL,  labels = c('A', 'B', 'C','' ), label_size = 12, align = "h", ncol = 2)

##########################################permanovas #####################################
install.packages("devtools")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


numerical<-wide%>%
  select(13:295)%>%
  select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data


siteblock<-wide%>%
  select(1:3, 12)%>%#p.env contains only the environmental (independant) variables 
  mutate(site=as.factor(site))%>%
  mutate(block=as.factor(block)) %>%
  mutate(year=as.factor(year))

#running the permanova test

#for all sites in block (general)
Permanova1 <- adonis(numerical ~ year, data=siteblock, permutations=999, method="bray"); Permanova1 # 0.049 #it doesnt seem significna t from the nmds

#for each pairwises of blocks
Permanova <- pairwise.adonis2(numerical ~ year: block,data=siteblock, permutations=999, method="bray"); Permanova #gives only 1 result, of year

#for each pairwises of blocks
Permanova <- pairwise.adonis2(numerical ~ block:year,data=siteblock, permutations=999, method="bray"); Permanova 


##2019 only 

wide9<-wide%>%
  filter(year=="2019")

numerical9<-wide9%>%
  select(13:295)%>%
  select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data


siteblock9<-wide9%>%
  select(1:3, 12)%>%#p.env contains only the environmental (independant) variables 
  mutate(site=as.factor(site))%>%
  mutate(block=as.factor(block))
  

#for all sites in block (general)
Permanova1 <- adonis(numerical9 ~ site:block+block,data=siteblock9, permutations=999, method="bray"); Permanova1 

#for each pairwises of blocks
Permanova <- pairwise.adonis2(numerical9 ~ block,data=siteblock9, permutations=999, method="bray"); Permanova 


########

bleachh<- Symportalrel%>%
  right_join(depthmeta,., by="ID")%>%
  filter(year=="2019")%>%
  right_join(meta19, by= "ID")%>%
  dplyr::select(1:12, 296)

numericalbleach<- bleachh%>%
  select(13)

Permanova1 <- adonis(numericalbleach ~ site:block+block,data=bleachh, permutations=999, method="bray"); Permanova1 

#for each pairwises of blocks
Permanova <- pairwise.adonis(numericalbleach ~ block, data=bleachh, permutations=999, method="bray"); Permanova 


########################################### hierarhical clustering analysis ####

numerical<-wide%>%
  filter(year==2019)%>%
  select(13:295)%>%
  select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data


siteblock<-wide%>%#p.env contains only the environmental (independant) variables
  select(1:3, 12)%>%
  filter(year== "2019")%>%
  mutate(site=as.factor(site))%>%
  mutate(block=as.factor(block)) %>%
  mutate(year=as.factor(year))

d <- dist(numerical,  method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)
#trying to color it

#______
#using another script

numericallbel<-wide%>%
  filter(year==2019)%>%
  select(1,13:295)

str(numericallbel)

library(dendextend)
numericallbel2 <- numericallbel[,-1]
row.names(numericallbel2) <- numericallbel$block




dend <- numericallbel2 %>% dist %>% hclust %>% as.dendrogram
dend <- dend %>% set("labels_colors", as.numeric(iris[,1]), order_value = TRUE) %>%
  set("labels_cex", .1)
par(mar = c(4,1,0,8))
plot(dend, horiz = T)







########################################## Sediments ##################################

setwd("~/Desktop/Clonality") 

#will work on the mastersedimentation file


mastersedimentation<-read_excel("mastersedimentation.xlsx")


#mean_sediment<-mastersedimentation%>%
#gather(., key= "sedsite", value= "value", "g/day1":"g/day7", factor_key=TRUE) %>%
# na.omit() #making it long

cleansed<- mastersedimentation %>% dplyr::select(c("site","S_mean","S_min", "S_max", "S_range", "Std")) 




Sed = ggplot(cleansed) + 
  theme_classic() +
  geom_point(aes(x = site, y = S_range)) + 
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 7), 
        legend.text = element_text(size = 9, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 10));Sed

########################################## Temperature ###################################
library(tidyverse)
library(janitor)
library(lubridate)
library(scales)

setwd("~/Desktop/Clonality")

raw<-read_csv("KBaySiteTemps.csv") %>%
  drop_na()%>%clean_names()%>% mutate(date_time=ymd_hms(date_time))%>%
  mutate(date=date(date_time))%>%
  dplyr::group_by(date,site)%>% dplyr::summarise(mean_temp=mean(temp)) %>%
  mutate(group=case_when(site=="mean_temp" ~ "mean",
                         site!="mean_temp" ~ "sites"))

testdaily<-read_csv("KBaySiteTemps.csv") %>%
  drop_na()%>%
  clean_names()%>% 
  mutate(date_time=ymd_hms(date_time))%>%
  mutate(date=date(date_time))%>%
  filter(date>=as_datetime("2019-01-27")) %>%
  select(date, temp, site) 
  #filter(site=="3_2")  #min temp is not 27!


kaneoherawt <-read_csv("KBaySiteTemps.csv") %>%
  mutate(date_time=ymd_hms(date_time))%>% 
  pivot_wider(names_from= site, values_from= temp)


dailytest19<- kaneoherawt %>%
  mutate(mean=rowMeans(.[,2:31],na.rm=TRUE))%>%
  select(., -nerr) %>%
  gather(site,temp,-date_time)%>%filter(temp!="NA")%>%
  group_by(site,date_time=floor_date(date_time,"day"))%>%summarise(temp=mean(temp))%>%
  mutate(group=case_when(site=="mean" ~ "Mean",
                         site!="mean" ~ "KBaySites")) %>%
  rename(.,"Site"="site", "Time"="date_time", "Temp"="temp") %>%
  filter(Time>=as_datetime("2019-01-20"))
 

install.packages(scales)


library(imputeTS)


pcaprep<- dailytest19%>%
  ungroup() %>%
  filter(Site!= "mean") %>%
  select(-group) %>% spread (Time, Temp)%>%
  filter(Site!="2_2")

meta<- pcaprep%>% select(Site)
meta$Site2 = meta$Site

meta2<-meta %>%
  separate(Site2,into=c("block","no"),sep="_",remove=FALSE)%>%
  select(Site, block) %>%
  rename(Block=block)

library(imputeTS)
pcadata<- pcaprep%>% select(-c(Site)) %>% mutate_if (., is.factor, as.character)%>%
  na_mean(., option="mean") #replaces missing values with mean 



library(factoextra)
pca<- prcomp(pcadata)
summary(pca)
axes<- fviz_pca_ind(pca, axes=c(1,2))
plotdata<- cbind(meta2, axes$data) %>%
  rename(Dim1= "x") %>%
  rename(Dim2= "y")

##variables
#var$coord: coordinates of variables to create a scatter plot
#var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
#var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).


################# test nmds following another script, adds elipsis into figure


set.seed(216)
NMDS=metaMDS(pcadata, distance='bray') 
nmdstemp<- bind_cols(meta2, NMDS$points[,1],NMDS$points[,2]) %>%
  rename(NMDS1=3, NMDS2=4)

a2019<-ggplot(nmdstemp) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none")+
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_point(aes(x=NMDS1, y=NMDS2, color=Block), size=2) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=Block), type='t')+
  scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333")); a2019 #type=t does not assume a normal distribution or a known standard deviation





nmds_results <- metaMDS(comm = numerical,  # Define the community data 
                        distance = "bray",       # Specify a bray-curtis distance
                        try = 100) 

####permanova temperature

permatemp<- cbind(meta2, pcadata)


#for all sites in block (general)
Permanova1temp <- adonis(pcadata ~ Site:Block+Block,data=permatemp, permutations=999, method="bray"); Permanova1temp 


#for each pairwises of blocks
Permanovatemp <- pairwise.adonis2(pcadata ~ Block,data=permatemp, permutations=999, method="bray"); Permanovatemp 

######2018 elypsis 


library(tidyverse)
library(janitor)
library(lubridate)
library(scales)

setwd("~/Desktop/Clonality")

raw<-read_csv("KBaySiteTemps.csv") %>%
  drop_na()%>%clean_names()%>% mutate(date_time=ymd_hms(date_time))%>%
  mutate(date=date(date_time))%>%
  dplyr::group_by(date,site)%>% dplyr::summarise(mean_temp=mean(temp)) %>%
  mutate(group=case_when(site=="mean_temp" ~ "mean",
                         site!="mean_temp" ~ "sites"))

testdaily<-read_csv("KBaySiteTemps.csv") %>%
  drop_na()%>%
  clean_names()%>% 
  mutate(date_time=ymd_hms(date_time))%>%
  mutate(date=date(date_time))%>%
  filter(date<=as_datetime("2019-01-27")) %>%
  select(date, temp, site) 
#filter(site=="3_2")  #min temp is not 27!


kaneoherawt <-read_csv("KBaySiteTemps.csv") %>%
  mutate(date_time=ymd_hms(date_time))%>% 
  pivot_wider(names_from= site, values_from= temp)


dailytest18<- kaneoherawt %>%
  mutate(mean=rowMeans(.[,2:31],na.rm=TRUE))%>%
  select(., -nerr) %>%
  gather(site,temp,-date_time)%>%filter(temp!="NA")%>%
  group_by(site,date_time=floor_date(date_time,"day"))%>%summarise(temp=mean(temp))%>%
  mutate(group=case_when(site=="mean" ~ "Mean",
                         site!="mean" ~ "KBaySites")) %>%
  rename(.,"Site"="site", "Time"="date_time", "Temp"="temp") %>%
  filter(Time<=as_datetime("2019-01-20"))


install.packages(scales)


library(imputeTS)


pcaprep<- dailytest18%>%
  ungroup() %>%
  filter(Site!= "mean") %>%
  select(-group) %>% spread (Time, Temp)%>%
  filter(Site!="2_2")

meta<- pcaprep%>% select(Site)
meta$Site2 = meta$Site

meta2<-meta %>%
  separate(Site2,into=c("block","no"),sep="_",remove=FALSE)%>%
  select(Site, block) %>%
  rename(Block=block)

library(imputeTS)
pcadata<- pcaprep%>% select(-c(Site)) %>% mutate_if (., is.factor, as.character)%>%
  na_mean(., option="mean") #replaces missing values with mean 

library(factoextra)
pca<- prcomp(pcadata)
summary(pca)
axes<- fviz_pca_ind(pca, axes=c(1,2))
plotdata<- cbind(meta2, axes$data) %>%
  rename(Dim1= "x") %>%
  rename(Dim2= "y")


################# test nmds following another script, adds elipsis into figure


set.seed(216)
NMDS=metaMDS(pcadata, distance='bray') 
nmdstemp<- bind_cols(meta2, NMDS$points[,1],NMDS$points[,2]) %>%
  rename(NMDS1=3, NMDS2=4)

a2018<-ggplot(nmdstemp) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="right")+
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_point(aes(x=NMDS1, y=NMDS2, color=Block), size=2) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=Block), type='t')+
  scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333")); a2018 #type=t does not assume a normal distribution or a known standard deviation




#type=t does not assume a normal distribution or a known standard deviation



install.packages(cowplot)
quartz()
plot_grid(a2018, a2019, labels = c('2018', '2019'), label_size = 10)




########################################## Distance-based redundancy analysis (db-RDA) ####
###***need to check why there are differences in results (2019 only?), (when in 2019?), 
#(half 2018-on?), (the whole time plus 2019?), (excluding samples low read numbers 2019?), 
#(doing mean of some sites?, or not doing that at all?)

###############calculating temp mean

#2019 only
setwd("~/Desktop/Clonality")
raw<-read_csv("KBaySiteTemps.csv") %>%
filter(site=="1_10")

library(lubridate)

#2019 only
x<-raw %>%
  filter(date_time>=as_datetime("2019-01-20"))  %>%
  filter(site!= "nerr") %>%
  drop_na() %>% clean_names() %>% 
  mutate(date_time=ymd_hms(date_time))%>%
  mutate(date=date(date_time))%>%
  dplyr::group_by(date,site)%>% 
  dplyr::mutate(mean_temp=mean(temp),dailyrange= max(temp)-min(temp), sd=sd(temp)) %>%
  group_by(site)%>%
  mutate(T_mean=mean(mean_temp), T_max=max(mean_temp), T_min=min(mean_temp), T_range= T_max-T_min,T_drange=mean(dailyrange),T_dsd=mean(sd, na.rm=TRUE))%>%
 # mutate(T_mean=case_when(site=="1_3" ~ mean,
                        # site=="1_9" ~ mean,
                       #  site== "1_10" ~ mean, 
                        #  site== "2_2" ~ mean,
                       #  site=="3_3" ~ mean, 
                        #site== "_3"~ mean, TRUE~as.numeric(T_mean))) %>%
  #mutate(T_min=case_when(site=="3_2"~ min,
                      # site=="4_11"~ min,
                       # site=="4_14" ~ min, 
                        # site=="2_2" ~ min, 
                       #site=="2_3"~ min, TRUE~as.numeric(T_min))) %>%
 #mutate(T_max=case_when(site== "4_8"~ max, TRUE~as.numeric(T_max)))%>%
  select(site, T_mean, T_min, T_max, T_range, T_drange, T_dsd)%>%
  distinct()

#quartz()
#ggplot(x)+
#geom_line(aes(x=date_time, y=temp))+geom_vline(aes(xintercept='2018-07-20')) +
#facet_wrap(~ site)



#mean max, min, range are always per site


#ggplot(x%>%sample_n(1000))+geom_point(aes(T_drange,T_dsd))



DHW<- raw %>%
  filter(date_time>=as_datetime("2019-01-20"))%>%
  #filter(date_time>=as_datetime("2018-07-20")) %>%
  mutate(dif=as.numeric(lag(date_time, 1L)-date_time)) %>% 
  filter(dif=="-10") %>%
  filter(temp>28.5) %>%
  group_by(site)%>%
  tally() %>%
  mutate(T_DHW=((n/6)/24)/7) %>%
  filter(site!= "nerr") %>%
  select(-c(2))


tempmerged<-left_join(DHW,x, by="site")
write.csv(tempmerged, "tempmerged.csv", row.names = FALSE)

tempmerged<-left_join(DHW,x, by="site")%>%
  left_join(., cleansed, by="site")


tempmeta<-left_join(wide, tempmerged, by= "site")%>%
  select(-year)


str(tempmeta)
saveRDS(tempmeta,"tempmeta")


################################## dbRA visualization ####

Processeddata<-readRDS("Processeddata")
tempmeta<-readRDS("tempmeta")


library(tidyverse) #need to unload plyr before running this
Data<-Processeddata %>% 
  filter(year==2019)%>%
  dplyr::group_by(ID)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads)%>%
  select(ID,clade,proportion, year, site)%>% 
  distinct()%>%
  spread(clade,proportion)%>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C= replace_na(C, 0))%>%
  mutate(Majorclade=case_when (
    C>0.8~"C",
    D>0.8~"D",
    TRUE~"CD"))%>%
  select(ID, year, Majorclade, site) %>%
  distinct() 

library(tidyverse)
testando<-Processeddata %>% 
  filter(year==2019)%>%
  dplyr::group_by(ID)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads)%>%
  select(ID,clade,proportion, year, site)%>% 
  distinct()%>%
  spread(clade,proportion)%>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C=replace_na(C, 0))%>%
  mutate(Majorclade=case_when (
    C>=0.99~"C",
    D>=0.99~"D",
    TRUE~"CD"))%>%
  select(ID, year, Majorclade, site, C,D) %>%
  distinct()%>% 
  filter(Majorclade=="CD")%>%
  mutate(Majorclade2=case_when (
    C>0.8~"C",
    D>0.8~"D",
    TRUE~"CD"))
  
  


library(plyr) 
count(testando, c("Majorclade2")) #C=226, CD 140, D 106 
count(Data, c("Majorclade")) #C=254, CD=63,D=155


tempmetaCD<-left_join(tempmeta, Data, by="ID")%>%
  dplyr::filter(year=="2019")





#changing order
tempmetaCD. <- tempmetaCD[c(1,3:4,2,293:307, 10:292)]%>%
  select(-c(19))
  


species<-tempmetaCD.[,19:301] 
environment= tempmetaCD.[,4:16] 
scaled_clim = scale(environment)

species001= (species + 0.001)

str(species001)
str(environment)
rankindex(environment, species001, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman") #Bra is 0.06920659

quartz()
dbRDA = capscale(species001 ~ depth_m + T_DHW + T_mean  + T_max + T_min  +  T_range + T_drange + T_dsd + S_mean + S_min + S_max + S_range + Std ,environment, dist="bray")

#bra is 0.1061

plot(dbRDA) 
anova(dbRDA) # is the model significant? #yes 0.001

#Permutation tests to access significant of constraints

set.seed(999)
anova(dbRDA) # overall test of the significant of the analysis 0.001
anova(dbRDA, by="axis", perm.max=200) # test axes for significance #it bugs
anova(dbRDA, by="terms", permu=999) # test for sign. environ. variables 




#making the graph in ggplot
x <- as.data.frame(scores(dbRDA, display = "sites"))
tempmetaCD.$CAP1 <- x$CAP1
tempmetaCD.$CAP2 <- x$CAP2 
d <- data.frame(Variables = rownames(dbRDA$CCA$biplot), dbRDA$CCA$biplot)%>%
  filter(Variables=="depth_m"|Variables=="T_DHW" | Variables== "T_range"|Variables=="T_drange"| Variables=="T_dsd"| Variables=="S_min")


quartz()
i<-ggplot(tempmetaCD., aes(x= CAP1, y= CAP2, color=Majorclade)) + 
  stat_ellipse(aes(fill = Majorclade,color=Majorclade), geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values=c("#14505C", "orange1", "#772C4B"))+
  scale_color_manual(values=c("#14505C", "orange1", "#772C4B"))+
  geom_point(alpha = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  coord_cartesian(xlim=c(-4, 5), ylim=c(-7, 8))  +
  geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
  geom_vline(xintercept = 0, linetype="dotted", color="grey") +
  labs(color = "Major Clades", fill = "Major Clades")  +
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_segment(data = d, aes(x = 0, y = 0, xend = (CAP1*5),yend = (CAP2*6)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate("text", x = (d$CAP1*5), y = (d$CAP2*5),label = d$Variables);i


##with enviro from Ford

setwd("~/Desktop/MRDS_map/")
enviro<-readRDS("environmental_summary")%>%select(-sample, -block)%>%distinct()

tempmetaCDF<-left_join(Data,enviro, by="site")%>%
  dplyr::filter(year=="2019")

#changing order
tempmetaCDF. <- tempmetaCD[c(1,2,295:306, 12:294, 308)]


species<-tempmetaCD.[,15:297] 
environment= tempmetaCD.[,2:14] 
scaled_clim = scale(environment)

species001= (species + 0.001)

str(species001)
str(environment)
rankindex(environment, species001, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman") #Bra is 0.06920659

quartz()
dbRDA = capscale(species001 ~ depth_m + T_DHW + T_mean + T_max + T_min + T_range + T_drange + T_dsd +  S_mean+ S_min +S_max  + Std ,environment, dist="bray")

#bra is 0.1061

plot(dbRDA) 
anova(dbRDA) # is the model significant? #yes 0.001

#Permutation tests to access significant of constraints

set.seed(999)
anova(dbRDA) # overall test of the significant of the analysis 0.001
anova(dbRDA, by="axis", perm.max=200) # test axes for significance #it bugs
anova(dbRDA, by="terms", permu=999) # test for sign. environ. variables 




#making the graph in ggplot
x <- as.data.frame(scores(dbRDA, display = "sites"))
tempmetaCD.$CAP1 <- x$CAP1
tempmetaCD.$CAP2 <- x$CAP2 
d <- data.frame(Variables = rownames(dbRDA$CCA$biplot), dbRDA$CCA$biplot)%>%
  filter(Variables=="depth_m"| Variables=="T_mean"| Variables=="T_max"| Variables=="T_min"|Variables=="T_DHW" | Variables== "S_min"| Variables=="T_dsd"| Variables=="Std")


quartz()
i<-ggplot(tempmetaCD., aes(x= CAP1, y= CAP2, color=Majorclade)) + 
  stat_ellipse(aes(fill = Majorclade,color=Majorclade), geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values=c("#14505C", "orange1", "#772C4B"))+
  scale_color_manual(values=c("#14505C", "orange1", "#772C4B"))+
  geom_point(alpha = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  coord_cartesian(xlim=c(-4, 5), ylim=c(-7, 8))  +
  geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
  geom_vline(xintercept = 0, linetype="dotted", color="grey") +
  labs(color = "Major Clades", fill = "Major Clades")  +
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_segment(data = d, aes(x = 0, y = 0, xend = (CAP1*5),yend = (CAP2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate("text", x = (d$CAP1*5), y = (d$CAP2*5),label = d$Variables);i



################################# shuffling and shifting ##################

library(tidyverse)
Data<-Processeddata %>% 
  group_by(ID, year)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, year, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads)%>%
  select(ID,clade,proportion, year, totalreads)%>%distinct()%>%
  group_by(ID, year)%>%
  spread(clade,proportion)%>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C=replace_na(C, 0))%>%
  mutate(Majorclade=case_when(
    C>0.51~"C",
    D>0.51~"D",
    TRUE~"CD"))%>%
  select(ID, year, Majorclade,totalreads )%>%
  distinct() 

tempmetaCD<-left_join(tempmeta, Data, by="ID")%>%
  dplyr::filter(year=="2019")


a<-Data %>% group_by(ID, year) %>% mutate(new = + (n_distinct(Majorclade) > 1)) #create new colum that says if c d hav changed in 2017 and 2019. 0=the same, 1 changed

b<-table(a$new); b



#it is saying no change!
###trying again



############################## Bleaching scores ###################################################

Dataall_<-Processeddata %>% 
  group_by(ID)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads) %>% 
  spread(clade,proportion)%>%
  dplyr::rename(prop_c=C,prop_d=D)%>%
  mutate(prop_d = replace_na(prop_d, 0)) %>%
  mutate(prop_c = 1-prop_d) %>%
  group_by(ID)%>%
  mutate(sumclade=prop_c + prop_d)%>%
  select("ID", "year", "prop_c", "prop_d", "block","bleachscore" )%>%
  distinct() 

bleachscore<-Dataall_%>%
  filter(year=="2019")%>%
  rename(bleach= bleachscore)%>%
  rename(Block= block)
  
  
bleach = bleachscore$bleach
prop_d = bleachscore$prop_d

bleachscore$bleach <- as.numeric(bleachscore$bleach)


res <- cor.test(bleach, prop_d, 
                method = "pearson"); res


#0.40 correlation

############################# scatterplot bleaching vs prop_d #################


g1 <- glm(prop_d~bleach*Block,family='binomial',  data=bleachscore)
summary(g1)

step(g1, test="LRT") #suggests variables that actually should be included

#why Block 1 does not show?
#odds of a coral bleaching according to block, temp, prop d, sedimentation, depth
g1 <- glm(I(bleach <=1 ) ~ Block,
          family = binomial(), 
          data = bleachscore)






quartz()
ggplot(bleachscore, aes(x=bleach, y=prop_d, color = Block, shape=Block)) + 
  #geom_point() +
geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  theme_classic() +
  scale_shape_manual(values=c(16,16,16, 16, 16)) + 
  scale_color_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270" )) +
  labs(x = "Bleach score", y = "Proportion of Durusdinium", fill = "Blocks")

#prop of D in block 1 and 5 remains almost constant, no higher prop of d in colonies who did not bleach
#corals in those locations are resistnt to bleaching not because of symbiont 


ggplot(data=bleachscore, aes(x=bleach, y=prop_d, ymin=sckT.lo, ymax=sckT.up, fill=type, linetype=type)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5) + 
  scale_x_log10() + 
  scale_y_log10() + 
  xlab(as.expression(expression( paste("Radius (", R[500], ")") ))) + 
  ylab("Scaled Temperature")

quartz()
p3 <- ggplot(bleachscore, aes(x=bleach, y=prop_d, color=Block)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum(); p3



#quartz()
#p <- ggplot(bleachscore, aes(bleach, prop_d, fill=Block)) +
  #geom_point(shape=21, size=2.5) + 
  #stat_smooth(method=lm, size=1.2, span=.6,linetype = "dashed", lwd=0.5, alpha=0.1, aes (color=Block)) +
  #theme_classic() +
  #scale_color_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270" )) +
  #labs(x = "Bleach score", y = "Proportion of Durusdinium")+
  #scale_fill_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270")); p




quartz()
p <- ggplot(bleachscore, aes(bleach, prop_d, fill=Block)) +
  #geom_point(shape=21, size=2.5) + 
  stat_smooth(method=lm, alpha=0.1, lwd=0.5, aes (color=Block)) +
  theme_classic() +
  scale_color_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270" )) +
  labs(x = "Health score", y = "Proportion of Durusdinium")+
  scale_fill_manual(values=c("#F8B83C","#5087C1","#8273B5","#398E88","#C53270")); p

######## another trial
library(ggplot2)












  
#So now the fill and colour aesthetics are mapped to different colour schemes. Either remove the offending scale_fill_brewer or add another in for colour

p <- p + scale_colour_brewer(palette = "Set1")


  
theme_classic()


############################ heatmap #############


bleachtrial<-bleachscore%>%
  select ( bleach, prop_d, Block)%>%
  distinct()

bleachtrial$prop_d.cat <- cut(bleachtrial$prop_d , breaks = 5, labels = c("20", "40", "60", "80","100%"))

bleachtrial<-bleachscore%>%
  select( bleach, prop_d, Block)%>%
  distinct()

bleachtrial$bleach <- as.numeric(bleachtrial$bleach)
bleachtrial$site <- as.factor(bleachtrial$site)
bleachtrial$prop_d<-as.numeric(bleachtrial$prop_d)
bleachtrial$Block<-as.factor(bleachtrial$Block)

quartz()
p <- ggplot(bleachtrial, aes(x=Block, y=bleach)) +
  geom_point(aes(color = prop_d)) +
  scale_fill_continuous(high="#26185F", low="#FCFFDD") +
  labs(title = "Bleach score per site",
       y = "bleach score");p



q3<- sequential_hcl(5, palette="YlGnBu"); q3 


####
bleachtrial<-bleachscore%>%
  select ( prop_d,prop_c, Block)%>%
  distinct() %>%
  gather(., prop_c, measurement, control:cond2, factor_key=TRUE)

  
  

##############################looking more closely at 1 and 5

block15<- bleachscore%>%
  filter(Block==1 | Block==5)





############################ Kaneohe Map ##################################
load(ggmap)
register_google(key = "AIzaSyAVD2dghJco5L0fYjv3v-um33JnQmMrQko",write = TRUE)


lat <- c(21.42, 21.50)
lon <- c(-157.85, -157.77)
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 13, maptype = "satellite", source = "google")
plot(map)
ggmap(map)+
  labs(x = 'Longitude', y = 'Latitude')

setwd("~/Desktop/Clonality")
mapdata<- read_excel("kaneoheclopoints.xlsx")
mapdata$Block <- as.character(mapdata$Block)


#adding points
quartz()
ggmap(map)+
  geom_point(data = mapdata,
             mapping = aes(x = Longitude, 
                           y = Latitude,
                           fill = Block), 
             size=3,
             pch = 23,
             colour = "black", 
             alpha = 0.9 ) +
  scale_color_manual(values=c("#F8B83C", "#5087C1","#8273B5","#398E88", "#C53270")) +
  scale_fill_manual(values=c("#F8B83C", "#5087C1","#8273B5","#398E88", "#C53270")) +
  theme_bw(base_size=9) +
  theme(legend.key.size = unit(0.9, "cm"))

quartz()
ggmap(map)+
  geom_point(data = mapdata,
             mapping = aes(x = Longitude, 
                           y = Latitude,
                           fill = Block), 
             size=3,
             pch = 23,
             colour = "black", 
             alpha = 0.9 ) +
  scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333")) +
  scale_fill_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333")) +
  theme_bw(base_size=9) +
  theme(legend.key.size = unit(0.9, "cm"))


#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"


theme_set(theme_bw() + 
            theme(legend.key=element_blank()))+
  theme(axis.text.x = element_text(size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7), legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 7))



#how to change size of letters in lon?
#how to not show lon and lat? only the values?
#how to show less numbers of lat and long?


quartz()

install.packages("usmap")
library(usmap)

usmap::plot_usmap(include = c("HI"))

## --------------------- black and white

library(sf)
library(rgdal)
library(ggplot2)

# Load shapefile of Oahu
oahu <- readOGR("~/Desktop/Oahu-coastline-shp/Coastline.shp")

# ggplot map:
oahu.spt <- spTransform(oahu, CRS("+proj=longlat +datum=WGS84"))
oahu.for <- fortify(oahu.spt)

setwd("~/Desktop/Clonality")
mapdata<- read_excel("kaneoheclopoints.xlsx")
mapdata$Block <- as.character(mapdata$Block)
mapdata$Site <- as.character(mapdata$Site)


quartz()
ggplot() +
  geom_polygon(data=oahu.for,aes(long, lat, group = group)) +
  theme_bw() +
  coord_sf(xlim = c(-157.86,-157.75), ylim = c(21.413, 21.5),
           expand = FALSE) +
  geom_point(data = mapdata,
             mapping = aes(x = Longitude, 
                           y = Latitude,
                           fill = Block),
             size=3,
             pch = 21,
             colour = "black", 
             alpha = 0.9 ) +
  scale_color_manual(values=c("#8EBC50","#F8DC5D" ,"#0E376F","#3A6BA5","#F99F00")) +
  scale_fill_manual(values=c("#8EBC50","#F8DC5D" ,"#0E376F","#3A6BA5","#F98800")) +
  theme_bw(base_size=9) +
  theme(legend.key.size = unit(0.9, "cm"))
    

#"#F39300","#BF004D" ,"#88F0E6","#520B4D" , "#294A49"
#"#2F327D" purple dark

library(colorspace)
library(dplyr)
library(ggplot2)
library(readxl)

q3<- sequential_hcl(9, palette="GnBu"); q3 
q3<- sequential_hcl(9, palette="YlOrRd"); q3 

#"#7D0025" dark red
#"#F39300" orange
#"#F6B732" "#F9D67E" yellows


#adding points
quartz()
ggmap(map)+
  geom_point(data = mapdata,
             mapping = aes(x = Longitude, 
                           y = Latitude,
                           fill = Block), 
             size=3,
             pch = 23,
             colour = "black", 
             alpha = 0.9 ) +
  scale_color_manual(values=c("#F8B83C", "#5087C1","#8273B5","#398E88", "#C53270")) +
  scale_fill_manual(values=c("#F8B83C", "#5087C1","#8273B5","#398E88", "#C53270")) +
  theme_bw(base_size=9) +
  theme(legend.key.size = unit(0.9, "cm"))


############################## Tryng interpolation ################################################

load(ggmap)
register_google(key = "AIzaSyAVD2dghJco5L0fYjv3v-um33JnQmMrQko",write = TRUE)


lat <- c(21.42, 21.50)
lon <- c(-157.85, -157.77)
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 13, maptype = "satellite", source = "google")
plot(map)
ggmap(map)+
  labs(x = 'Longitude', y = 'Latitude')

setwd("~/Desktop/Clonality")
#mapdata<- read_excel("kaneoheclopoints.xlsx") %>%

mapdata<- Processeddata19%>%
  select("ID", "year", "depth_m", "block", "lat", "lon", "bleachscore")%>%
  filter(year==2019) %>%
  distinct() %>%
  group_by(block)%>%
  mutate(meanblock=mean(bleachscore))
  

mapdata$bleaching <- as.factor(mapdata$bleaching)


#adding points
quartz()
ggmap(map)+
  geom_point(data = mapdata,
             mapping = aes(x = lon, y = lat, fill = bleaching), size=3,
             pch = 23,
             colour = bleaching,size=bleaching,
             alpha = 0.9 )  +
  scale_colour_brewer(palette ="Set1")+
  theme_bw(base_size=9) +
  theme(legend.key.size = unit(0.9, "cm"))

#it doesnt work!!

############################## kriging using script from rpubs #########


  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
  library(sp)
  library(gstat)


mapdata<- Processeddata19 %>%
  select("ID", "year", "depth_m", "block", "lat", "lon", "bleachscore", "site") %>%
  filter(year==2019) %>%
  distinct() 

mapdata <- mapdata[!is.na(mapdata$bleachscore),]
boxplot(bleachscore ~ site, mapdata, las=2)

q<- mapdata%>%
  group_by(site)%>%
  mutate(meansite=mean(site))

site <- aggregate(bleachscore ~ site, mapdata, mean)
site$sd <- aggregate(bleachscore ~ site, mapdata, sd)$bleachscore
site$lat <- aggregate(lat ~ site, mapdata, mean)$lat
site$lon <- aggregate(lon ~ site, mapdata, mean)$lon



quartz()
c<-mapdata %>% 
  as.data.frame %>% 
  ggplot(aes(lat, lon)) + geom_point(aes(size=bleachscore), position=position_jitter(h=0.01, w=0.01), color="blue", shape= 24, alpha=0.5) + 
  ggtitle("Bleach score") + coord_equal() + theme_bw(); c



coordinates(site) <- ~ lat + lon
class(site)
plot(site)
axis(1)
axis(2)
#fitting the variogram

lzn.vgm <- variogram(log(bleachscore) ~ 1, site) # calculates sample variogram values, does not work, probably cause of missing data 
lzn.vgm <- variogram(sqrt(bleachscore) ~ 1, site) # calculates sample variogram values, does not work, probably cause of missing data 
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1))

plot(lzn.vgm, lzn.fit) # plot the sample values, along with the fit model

 # Make grid for predictions
res <- 100
lon <- seq(min(site$lon), max(site$lon), length.out=res)
lat <- seq(min(site$lat), max(site$lat), length.out=res)
gg <- expand.grid(lon=lon, lat=lat)

coordinates(gg) <- ~ lat + lon # step 3 above
lzn.kriged <- krige(sqrt(bleachscore) ~ 1, site, gg, model=lzn.fit)

lzn.kriged %>% as.data.frame %>%
  ggplot(aes(x=lat, y=lon)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw()

plot(site)

plot1 <- site %>% as.data.frame %>%
  ggplot(aes(lat, lon)) + geom_point(size=1) + coord_equal() + 
  ggtitle("Points with measurements")



pred_bleachscore <- krige.conv(data=site$bleachscore, coords=ls1[c("lat", "lon")], locations=gg, krige=krige.control(cov.pars=c(ml_max$par["sigmasq",]$values, ml_max$par["phi",]$values), nugget=0))

# this is clearly gridded over the region of interest
plot2 <- meuse.grid %>% as.data.frame %>%
  ggplot(aes(x, y)) + geom_point(size=1) + coord_equal() + 
  ggtitle("Points at which to estimate")

library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)



########################## from Callum #####

# KRIGING
library(geoR)
library(raster)
library(rgdal)


# Load data
mapdata<- Processeddata19 %>%
  select("ID", "year", "depth_m", "block", "lat", "lon", "bleachscore", "site") %>%
  filter(year==2019) %>%
  distinct() 

mapdata <- mapdata[!is.na(mapdata$bleachscore),]


#calculating mean bleachscore per site
site <- aggregate(bleachscore ~ site, mapdata, mean)
site$sd <- aggregate(bleachscore ~ site, mapdata, sd)$bleachscore
site$lat <- aggregate(lat ~ site, mapdata, mean)$lat
site$lon <- aggregate(lon ~ site, mapdata, mean)$lon



# sites <- na.omit(sites)

map
site$bleachscore <- site$mmm
site$bleachscore <- log10(site$bleachscore + 0.01)
plot(lat ~ lon, site, cex=sqrt(site$bleachscore))

# Calculate variogram
ml_dhw <- likfit(data=site$bleachscore, coords=site[c("lon", "lat")], ini=c(0.5, 0.005), fix.nug = TRUE)
summary(ml_dhw)

# Make grid for predictions
res <- 100
lon <- seq(min(site$lon)-0.005, max(site$lon)+0.005, length.out=res)
lat <- seq(min(site$lat)-0.005, max(site$lat)+0.005, length.out=res)
gg <- expand.grid(lon=lon, lat=lat)

# Krige
pred_dhw <- krige.conv(data=site$bleachscore, coords=site[c("lon", "lat")], locations=gg, krige=krige.control(cov.pars=c(ml_dhw$par["sigmasq",]$values, ml_dhw$par["phi",]$values), nugget=0))

library(rgdal)
gis1 <- readOGR("~/Desktop/Clonality/Oahu", "Oahu2", verbose=FALSE)
gis2 <- readOGR("~/Desktop/Clonality/Oahu", "Oahu3", verbose=FALSE)

mat <- matrix(pred_dhw$predict, res, res)
mat <- apply(t(mat), 2, rev)
mat <- raster(mat, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=crs(gis1))

mat <- mask(mat, gis1, inverse=TRUE)

png("output/krige_mmm.png", width=6, height=6.45, units = 'in', res = 300)

quartz()
plot(mat, col=heat.colors(100, rev=TRUE, alpha=0.5), main="Bleaching score")
plot(gis1, add=TRUE)
plot(gis2, add=TRUE, border=rgb(0,0,0,0.1))
text(sites$lon, sites$lat, sites$site, cex=0.7)

dev.off()


################

setwd("~/Desktop")
Processeddatameso<-readRDS("Alldata_mesocosm")

m<- Processeddatameso%>%
  select(sample_name, fastq_fwd_file_name.x, Spp, Site)%>%
  separate(sample_name,into=c("dataset","sampleID"),sep="ID",remove=FALSE)

sort(m$sampleID, decreasing=T)

ndx = order(m$sampleID, decreasing = T)



#######Not used box plot############




quartz()
D<- ggplot(a, aes(x = site, y = prop_d)) +
  geom_bar(stat = "identity") + 
  facet_wrap (~year, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7), legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), 
        legend.margin=margin(-0.4,0,0,0,unit="cm")) +
  scale_y_continuous(expand = c(0,0)); D 


geom_point(alpha=0.3, color="tomato", position = "jitter") +
  geom_boxplot(alpha=0) + coord_flip() + facet_wrap( ~ sex)

#it works! year side by side
quartz()
trial<-ggplot(Dataall, aes(factor(site), prop_d)) +
  geom_boxplot(aes(colour = year), varwidth = TRUE)+
  # geom_point(alpha=0.2, color=year, position = "jitter")+
  theme(legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(0.1, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm")) +
  scale_y_continuous(expand = c(0,0)); trial



boxD + theme_classic()

colours=c("#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47")

a<-ggplot(Dataall, aes(x=block, y=prop_d, fill=block)) + 
  geom_boxplot(alpha=0.3) +
  facet_wrap (~year, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.5, hjust = 1),     
        axis.title.y = element_text(size = 13), legend.title = element_text(size = 13), 
        legend.text = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Proportion of D")+
  scale_fill_manual(values=colours) +
  theme_classic(); a




boxD + theme_classic()

colours=c("#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47","#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47", "#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47", "#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47", "#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47", "#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47", "#00798c","#8D96a3", "#66a182", "#2e4057", "#DAFF47" )

a<-ggplot(Dataall, aes(x=site, y=prop_d, fill=block)) + 
  geom_boxplot(alpha=0.3) +
  facet_wrap (~year, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.5, hjust = 1),      axis.title.y = element_text(size = 13), legend.title = element_text(size = 13), 
        legend.text = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Proportion of D")+
  scale_fill_manual(values=colours) +theme_classic(); a

#barplot prop D




#does not add to 1
Dataallyear<-Processeddata %>% 
  group_by(year, block)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(year,block, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads) %>% 
  spread(clade, proportion)%>%
  dplyr::rename(prop_c=C,prop_d=D)%>%
  mutate(prop_d = replace_na(prop_d, 0)) %>%
  mutate(prop_c = 1-prop_d) 

#prop of d does not add to 1, want to do it by year. Only 2 bars, with std 

quartz()
a<-ggplot(Dataallyear, aes(x=year, y=prop_d)) + 
  geom_boxplot(alpha=0.3) +
  facet_wrap (~year, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.5, hjust = 1),  
        axis.title.y = element_text(size = 13), legend.title = element_text(size = 13), 
        legend.text = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Proportion of D")+
  scale_fill_manual(values=colours) +theme_classic(); a


quartz()
mx4 = ggplot(Dataallyear, aes(x = year, fill = year, y = prop_d)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 7),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm"))  + 
  scale_y_continuous(expand = c(0,0)) + 
  #ggtitle(label = "A") +
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  #scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 1)); mx4
