
rm(list=ls())
## open the libraries
library(tidyverse);library(scales);library(ggthemes);library(ggplot2);library(cowplot);library(ggpubr)
library(betapart)
library(nlme); library(writexl); library(xtable); library(emmeans);library(piecewiseSEM)

## color needed 
colorblind_pal()(8)
show_col(colorblind_pal()(8))
## environment needed
pd<-position_dodge(width=0.5)
mt<-theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())

## set up the work directory 
dir.data<-"H:/resistance and recovery/R codes/raw data/"
dir.graphs<-"H:/resistance and recovery/R codes/graphs/"
setwd(dir.graphs)
## add data for biomass and cover
load("data for control and NPK.rdata")
d8<-read.csv("biomass and cover data for 55 nutnet sites.csv")
#################################################################################################
##################################### growing season data #######################################
#################################################################################################
## find the start and end of the experiment (year_trt 1 to ...)
## note that we have data from look.us and lagoas.br from year_trt 2!!!!
gs<-read.csv("information for 55 nutnet sites.csv")
gs4 <- d8 %>% group_by(site_code) %>% summarise_at("year", list(max, min)) %>%
  mutate(min.y=ifelse((site_code == "look.us" | site_code == "lagoas.br"), fn2 - 1, fn2)) %>% mutate(fn2 = NULL) %>%
  dplyr::rename(site_lyear = fn1, site_fyear = min.y, fn1 = NULL, fn2 = NULL) %>%
  merge(y = gs, by = c("site_code"))
length(unique(gs4$site_code))

#################################################################################################
###################################### water balance data #######################################
#################################################################################################
## get the data of precipitation and potential evaporation
prec <- read.csv(paste0(dir.data, "CRU-monthly-pre-1901-2021.csv"), header = T) %>%dplyr::rename(prec = value, year=years, month=months)
table(prec$year)
pet<-read.csv(paste0(dir.data, "CRU-monthly-pet-1901-2021.csv"), header = T, sep=",") %>%dplyr::rename(pet = value, year=years, month=months)
  
### note that potential evaporation was calculated as average per day, need to transform to monthly data. 
month.df <- data.frame(month = c("Jan","Feb", "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                       month1 = seq(1, 12, 1),
                       days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

dd0 <- prec %>% merge(y = pet, by=c("site_code", "year", "month", "date")) %>%
  merge(y = month.df, by=c("month"))  %>%
  mutate(pet1 = pet * days) %>%
  arrange(site_code, year, month1)

dd <- dd0 %>%
  bind_rows(dd0 %>% filter(site_code == "msla.us") %>% mutate(site_code = "msla_2.us")) %>%
  bind_rows(dd0 %>% filter(site_code == "msla.us") %>% mutate(site_code = "msla_3.us"))

## add the growing season data 
dd1a <- dd %>% merge(y = gs4[ , c("site_code", "gs_start_month", "harvest_month")], by = "site_code") %>% arrange(site_code, year, month1)

# at most sites, growing seasons are within a year
# at some sites, growing season start in Oct. or Nov. and peak biomass season was in the next year
# do it separately for these two situations
dd1 <- dd1a %>% mutate(dif = harvest_month - gs_start_month, one.v.two = ifelse(dif < 0, "two", "one"))
 
## summarize growing season precipitation and evapotranspiration for sites with growing season within one year
colnames(dd1)
s_dd2.within <- dd1 %>%
  filter(one.v.two == "one") %>%filter(month1>=gs_start_month & month1 <= harvest_month)%>%
  group_by(site_code, year) %>% 
  summarise_at(c("prec", "pet1"), list(sum))
unique(s_dd2.within$site_code)
range(s_dd2.within$year)

## summarize growing season precipitation and evapotranspiration for sites with growing season across two years
## months in the late season of the previous year should be counted to the next year in which peak biomass is recorded!!!! 
## in this case, data in 1901 would be incomplete for the peak biomass in that year, because the starting 
## growing months in the previous year (1900) were not available. Thus, delete data in 1901.
# delete year 2022, because we don't have climate data in that year.
## summarize growing season precipitation and evapotranspiration during growing season 
s_dd2.across <- dd1 %>%
  filter(one.v.two == "two")%>%filter(month1>=gs_start_month | month1 <= harvest_month)%>%
  mutate(year1 = ifelse(month1 >= gs_start_month, year+1, year)) %>%
  filter(year1 != 1901) %>% filter(year1 != 2022) %>% 
  group_by(site_code, year1) %>% 
  summarise_at(c("prec", "pet1"), list(sum))
unique(s_dd2.across$site_code)
range(s_dd2.across$year1)

## for sites where growing season is within a year, whole growing season was recorded in 1901,  
## to be consistent, delete data in 1901 
## merge all the sites 
# Compute water balance and standardize it 
s_dd3 <- s_dd2.within %>% filter(year!=1901) %>%  dplyr::rename(year1=year) %>% 
  bind_rows(s_dd2.across) %>%
  mutate(dif.re0 = prec - pet1) %>%
  group_by(site_code) %>%mutate(spei = scale(dif.re0))
range(s_dd3$year1)
length(unique(s_dd3$site_code))

## check the distribution of the data 
ggplot(s_dd3, aes(year1, spei)) + theme_bw() + geom_point() +
  facet_wrap(~site_code, scale="free_y")+
  geom_hline(yintercept = 0, linetype="dashed", color="red")
che <- s_dd3 %>% group_by(site_code) %>% summarise_at("spei", list(mean, sd))## 

## add the start and end of the experiments 
s_dd4 <- merge(s_dd3, gs4[,c("site_code", "site_fyear", "site_lyear")], by=c("site_code")) 

#################################################################################################
##################################compare water balance at all sites#############################
#################################################################################################
## calculate water balance from 2007 to 2021 for all sites
## add the growing season data 
avg.bal1 <- s_dd4 %>% filter(year1 > 2006) %>% group_by(site_code) %>%
  summarise_at("dif.re0", list(mean)) %>%
  dplyr::rename(water.balance.last.15 = dif.re0) %>% 
  merge(y = gs4[ , c("site_code", "gs_start_month",  "harvest_month")], by = "site_code")

#################################################################################################
############################## categorize dry and wet climate events  ###########################
#################################################################################################
## include the pre-treatment year (year_trt 0)
s_dd5 <- s_dd4 %>% filter(year1 >= site_fyear - 1 & year1 <= site_lyear) %>% mutate(year = year1) 

s_dd5$climate<-"Normal"
s_dd5$climate[s_dd5$spei>=1.28]<-"Extreme wet"
s_dd5$climate[s_dd5$spei<=-1.28]<-"Extreme dry"
s_dd5$climate[s_dd5$spei>=0.67 & s_dd5$spei<1.28]<-"Moderate wet"
s_dd5$climate[s_dd5$spei<=-0.67 & s_dd5$spei>-1.28]<-"Moderate dry"
## check whether all sites have extreme and normal years 
s_dd5_n<-subset(s_dd5, climate=="Normal" & year>=site_fyear)
s_dd5_e<-subset(s_dd5, climate!="Normal" & year>=site_fyear)
## find sites without normal years or extreme years 
ss1<-unique(s_dd5_n$site_code)
ss2<-unique(s_dd5_e$site_code)
ss3<-levels(factor(dd1$site_code))
(sss1 <- ss3[which(!ss3 %in% ss1)])## 
(sss2 <- ss3[which(!ss3 %in% ss2)])## 
## bind these two data sets
(sss3<-c(sss1, sss2))
## delete sites without extreme or normal events, also add year_trt
s_dd6<-s_dd5%>%filter(!site_code%in%sss3)%>%mutate(year_trt=year-site_fyear+1)
colnames(s_dd6)
# focusing on sites with at least four-year data 
# indicate missing data 
# get sites and years with biomass data
# add number of years 
num.year.at.least.4<-d8%>%select(site_code, year, year_trt)%>%distinct()%>%group_by(site_code)%>%
  summarise(number.of.years=length(year_trt))%>%filter(number.of.years>3)
intact.bio <- d8 %>% select(site_code, year_trt) %>% distinct() %>% mutate(bio.data = "Intact")
s_dd7<-s_dd6%>%filter(site_code %in%num.year.at.least.4$site_code)%>%
  merge(y=intact.bio, by=c("site_code", "year_trt"), all.x=T)%>%
  mutate(bio.data1=ifelse(is.na(bio.data), "Missing", bio.data))%>%
  mutate(bio.data2=ifelse(year_trt==0, "Intact", bio.data1))%>%
  filter(!site_code %in%c("kibber.in"))
 
# relevel climate events
s_dd7$climate<-factor(s_dd7$climate, levels=c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
(pp<-ggplot(s_dd7, aes(year_trt, spei, color=climate, shape=bio.data2))+geom_point(size=5)+mt+theme_cowplot(font_size = 25)+
    facet_wrap(~site_code, scale="free_x", ncol=5)+
    scale_x_continuous(expand=c(0,0), limits = c(-0.5,15), breaks = seq(0,16,2))+
    scale_color_manual(values = c("#E69F00","#F0E442",  "#000000", "#56B4E9",  "#0072B2"))+
    scale_shape_manual(values=c(19, 17))+
    geom_vline(xintercept = 0, linetype="dashed")+
    labs(x=NULL, y="SPEI", shape="Biomass data", color="Climate extremes"))
# delete site kibber.in, this site has only one normal year, but not biomass data 
# check climate events
table(s_dd7$climate)
# ggsave(pp, file="climate extremes.pdf", width = 21, height = 29.7, dpi=600)

# write.csv(s_dd7, file="raw data of moderate and extreme growing seasons at 55 sites.csv")

# also save the trend of SPEI
(pp.spei.trend<-s_dd3%>%filter(site_code %in%num.year.at.least.4$site_code)%>%
    filter(!site_code %in%c("kibber.in"))%>%
    ggplot(aes(year1, spei)) +mt+theme_cowplot(font_size = 20)+
    geom_point()+
    facet_wrap(~site_code, scale="free_y", ncol=5)+
    geom_hline(yintercept = 0, linetype="dashed", color="red")+
    labs(x="Year", y="SPEI"))
# ggsave(pp.spei.trend, file="SPEI over time.pdf", width = 21, height = 29.7, dpi=600)

###########################################################################################
################ calculate resistance and recovery for all non-extreme years ##############
###########################################################################################
data.bio.rich<-d8%>%mutate(biomass=live_mass)%>%
  select(site_code, block, plot, trt, year_trt, year, biomass, richness)%>%distinct()%>%
  pivot_longer(cols = c("biomass", "richness"), names_to ="community.property", values_to = "property.value" )
# merge biomass and richness with climate extremes 
d.select<-s_dd7%>%filter(year_trt!=0)%>%select("site_code", "year", "year_trt", "spei", "climate")%>% 
  merge(data.bio.rich, by=c("site_code", "year", "year_trt"))%>%
  mutate(variable.id=paste(site_code, block, trt, community.property, sep="_"))%>%arrange(variable.id)

## calculate mean values in normal years 
d.select_norm<-d.select%>%filter(climate=="Normal")%>% group_by(variable.id) %>% dplyr::summarise(avg.property=mean(property.value))

d.select_r1<-d.select%>%merge(d.select_norm, by=c("variable.id"))%>% mutate(deviation.property=abs(property.value - avg.property), resistance=avg.property/abs(property.value - avg.property))
 
# get resistance data by deleting normal climate events 
resis1<-d.select_r1%>%filter(climate!="Normal")%>%mutate(stability.facets="resistance")%>%dplyr::rename(values=resistance)%>%
  select(variable.id, year, year_trt, climate, spei, stability.facets, values)

# make sure that recovery is always calculated as first year compared with the second year following the first
recov<-c()
for(pl in unique(d.select_r1$variable.id)){
  # pl<-"lancaster.uk_2_Control_all_biomass"
  t.data2<-d.select_r1%>%filter(variable.id==pl)%>%arrange(year_trt)

  t.recovery<-c()
  for(i in 1:nrow(t.data2)){
    avg.property<-t.data2$avg.property[1]
    t.recovery[i]<-abs(t.data2$property.value[i] - avg.property)/abs(t.data2$property.value[i+1]-avg.property)
  }
  # make sure that recovery is always calculated as first year compared with the second year
  t.data3<-t.data2%>%mutate(consecutive=lead(year_trt)- year_trt, recovery=ifelse(consecutive==1, t.recovery, NA))
  recov<-rbind(recov, t.data3)
}  

recov1<-recov%>%filter(climate!="Normal")%>%mutate(stability.facets="recovery")%>%dplyr::rename(values=recovery)%>%
  select(variable.id, year, year_trt, climate, spei,  stability.facets, values)
 
# add resistance and recovery together 
resis.recov<-resis1%>%bind_rows(recov1)%>%merge(d.select%>%select(site_code, block, trt, community.property, variable.id)%>%distinct(), by=c("variable.id"))

###################### resistance and recovery for composition ############################# 

d.select.com<-d8%>%
  merge(s_dd7[,c("site_code", "year", "spei", "climate")], by=c("site_code", "year"))%>%
  mutate(community.property="composition", variable.id=paste(site_code, block, trt,community.property, sep="_"))%>%arrange(variable.id)

## calculate mean abundance for each species in normal years 
s_n<-d.select.com%>%filter(climate=="Normal")%>%group_by(variable.id, standard_taxon)%>%summarise(cover.avg=mean(max_cover))

## calculate resistance and recovery
rr.com<-c()
for(pl in unique(d.select.com$variable.id)){
  # pl<-"rook.uk_1_NPK_all_composition"; yr<-1
  ref.data<-s_n%>%filter(variable.id==pl)%>%dplyr::rename(max_cover=cover.avg)%>%mutate(reference.com="yes")
  data<-d.select.com%>%filter(variable.id==pl)
  sim.t<-c()
  for(yr in unique(data$year_trt)){
    data1<-data%>%filter(year_trt==yr)%>%dplyr::select(variable.id, standard_taxon, max_cover)%>%mutate(reference.com="no")
    ## add two communities together 
    t.com<-data1%>%bind_rows(ref.data)%>%pivot_wider(names_from = standard_taxon, values_from = max_cover)%>%
      replace(is.na(.), 0)%>%mutate(variable.id=NULL, reference.com=NULL)%>%dplyr::select(which(!colSums(.)==0))
    if(nrow(t.com)<2) next
    dis<-beta.pair.abund(t.com, index.family = "bray")
    sim<-as.numeric(1-dis$beta.bray)
    sim.temp<-data.frame(variable.id=pl, year_trt=yr, sim=sim)
    sim.t<-rbind(sim.t, sim.temp)
  }
  rr.com<-rbind(rr.com, sim.t)
}
rr.com1<-d.select.com%>%dplyr::select("variable.id", "site_code", "trt", "block", "year_trt", "climate", "spei")%>%distinct()%>%
  merge(y=rr.com, by=c("variable.id", "year_trt"))%>%mutate(resistance=sim) # higher values indicate higher resistance
# check resistance 
check.resis<-rr.com1%>%filter(climate!="Normal")%>%dplyr::select(site_code, year_trt)%>%unique()
# 
resis_com<-rr.com1%>%filter(climate!="Normal")%>%mutate(stability.facets="resistance")%>%dplyr::rename(values=sim)%>%
  dplyr::select(variable.id, year_trt, climate, spei, stability.facets, values)

## calculate recovery. a community has higher recovery when a community has higher similarity
## one year after the climate event than that under the climate event (both compared with the average of normal years)
recov_com<-rr.com1%>%arrange(variable.id, year_trt)%>%group_by(variable.id)%>%
  # make sure that recovery is always calculated as first year compared with the second year
  mutate(recovery=lead(sim)/sim, consecutive=lead(year_trt)- year_trt, recovery1=ifelse(consecutive==1, recovery, NA))%>%
  mutate(stability.facets="recovery")%>%dplyr::rename(values=recovery1)%>%
  filter(climate!="Normal")%>%dplyr::select(variable.id, year_trt, climate, spei, stability.facets, values)
# check recovery 
# check.recov<-recov_com%>%filter(climate!="Normal")%>%mutate(dif=recovery-recovery1)%>%filter(dif!=0)

# add resistance and recovery for three community aspects 
colnames(resis_com)
colnames(recov_com)
resis.recov.com<-resis_com%>%bind_rows(recov_com)%>%
   merge(d.select.com%>%select(site_code, block, trt, variable.id)%>%distinct(), by=c("variable.id"))%>%
  mutate(community.property="composition")

# add this to richness and biomass data
colnames(resis.recov)
colnames(resis.recov.com)
resis.recov1<-resis.recov%>%select(variable.id, year_trt, climate, spei, stability.facets, values, site_code,  block,  trt, community.property)%>%distinct()%>%
  bind_rows(resis.recov.com)
colnames(resis.recov1)
# select years for resistance and recovery
# for resistance and recovery, if the previous year is not the same extreme events, this year should not be included due to confounding effects
# for recovery, if the next year is a different extreme event, this year should not be included due to confounding effects
# the pre-treatment year and years with missing biomass should be included for selection 
num.events<-as.data.frame(table(s_dd7$climate))%>%mutate(data.type="original")

select.years.extremes<-s_dd7%>%select(site_code, year_trt, climate, spei)%>%filter(climate!="Normal")%>%
  mutate(spei1=as.numeric(spei),climate.num.0.67=case_when((spei1>=0.67)~1, (spei1<=- 0.67)~ -1, TRUE~99))%>%
  group_by(site_code)%>%
  mutate(consecutive.later=lead(year_trt)-year_trt, climate.type.later=lead(climate.num.0.67)+climate.num.0.67,
         consecutive.previous=lag(year_trt)-year_trt, climate.type.previous=lag(climate.num.0.67)+climate.num.0.67)%>%
  mutate(across(c("consecutive.later", "climate.type.later", "consecutive.previous", "climate.type.previous"), ~replace_na(.,999))) %>%
  mutate(year.resis=ifelse((consecutive.previous==-1 & climate.type.previous==0), NA, year_trt))%>%
  mutate(year.recov=ifelse((consecutive.later==1 & climate.type.later==0), NA, year.resis))%>%
  ungroup()%>%select(site_code, climate, year.resis, year.recov)%>%
  pivot_longer(cols = c("year.resis", "year.recov"), names_to ="year.used" , values_to ="year_trt")%>%
  mutate(stability.facets=ifelse(grepl("resis", year.used), "resistance", "recovery"), year.used=NULL)%>%filter(!is.na(year_trt))
num.events.exc.consecutive<-select.years.extremes%>%filter(year_trt!=0)%>%group_by(stability.facets, climate)%>%summarise(N=length(year_trt))
  
s_dd8_year.used.or.not<-select.years.extremes%>%select(site_code, year_trt)%>%distinct()%>%mutate(year.used="Yes")%>%
  merge(y=s_dd7[,c("site_code", "year_trt", "climate", "bio.data2", "spei")], by=c("site_code", "year_trt"), all.y=T)%>%
  mutate(year.used.1=ifelse(is.na(year.used), "No", year.used))%>%
  mutate(year.used.2=ifelse((year_trt==0), "No", year.used.1), year.used.3=ifelse((climate=="Normal"), "Yes", year.used.2))

cut<-c("0.67 and 1.28")
(pp1<-s_dd8_year.used.or.not%>%
    ggplot(aes(year_trt, spei, color=climate, shape=bio.data2, alpha=year.used.3))+geom_point(size=7)+
    theme_cowplot(font_size = 25)+panel_border()+
    facet_wrap(~site_code,  ncol=5, scale="free_x")+
    scale_x_continuous(expand=c(0,0), limits = c(-0.5,15), breaks = seq(0,16,2))+
    scale_color_manual(values = c("#E69F00","#F0E442",  "#000000", "#56B4E9",  "#0072B2"))+
    scale_shape_manual(values=c(19, 17))+
    scale_alpha_manual(values=c(0.1, 1))+
    geom_vline(xintercept = 0, linetype="dashed")+
    theme(legend.title=element_blank())+
    labs(x=NULL, shape=NULL, y=paste0("SPEI during growing season based on cutoff of ", cut)))
# ggsave(pp1, file="SPEI based on cutoff 0.67 and 1.28 indicating years used.pdf", width = 21, height = 29.7, dpi=600)

# when two same climate extremes happen consecutively, recovery only calculate for the later year 
# the later year should be followed by either a less extreme year or normal year
select.resis.recov<-c()
for(cut in c(0.67, 1, 1.28, 1.5, 2)){
  # cut<-1.28
     resis.recov2<-resis.recov1%>%
        mutate(spei1=as.numeric(spei), climate1=case_when((spei1>=cut)~"Wet",
                                (spei1<=- cut)~"Dry",
                                TRUE~"between.normal.and.extreme"))%>%
     arrange(stability.facets, variable.id, year_trt)%>%group_by(stability.facets, variable.id)%>%
     mutate(climate.num=case_when((spei1>=cut)~1,
                                 (spei1<=- cut)~ -1,
                                 TRUE~0), consecutive=lead(year_trt)-year_trt, climate.type=lead(climate.num)-climate.num, consecutive1=ifelse(is.na(consecutive), 999, consecutive))%>%
    mutate(values3=ifelse((stability.facets=="recovery" & consecutive1==1 & climate.type==0), NA, values))%>%
    merge(select.years.extremes, by=c("stability.facets", "site_code", "year_trt", "climate"), all.y=T)%>%
    filter(climate1%in% c("Wet", "Dry"))%>%
    select(variable.id, year_trt, climate1, spei1, stability.facets, values3, site_code,  block,  trt, community.property)%>%distinct()%>%
     filter(!is.na(values3))%>%mutate(cutoff=cut)
  
  select.resis.recov<-rbind(select.resis.recov, resis.recov2)
}

(check.year.0.67<-select.resis.recov%>%ungroup()%>%
    select(cutoff, site_code, year_trt, stability.facets, climate1, spei1)%>%filter(cutoff==0.67)%>%distinct()%>%
    ggplot(aes(year_trt, spei1, color=climate1, shape=stability.facets))+theme_cowplot(font_size = 25)+panel_border()+
    geom_point(size=5, alpha=0.7, position = position_dodge2(w = 0.5))+
    facet_wrap(~site_code,  ncol=5, scale="free_x")+
    scale_x_continuous(expand=c(0,0), limits = c(0.5,15), breaks = seq(0,16,2))+
    scale_color_manual(values = c("#E69F00",  "#0072B2"))+
      geom_vline(xintercept = 0, linetype="dashed")+
    theme(legend.title=element_blank())+
    labs(x=NULL, shape=NULL, y="SPEI"))
# ggsave(check.year.0.67, file="years used for resistance and recovery based on cutoff 0.67.pdf", width = 21, height = 29.7, dpi=600)

# average resistance during dry and wet and recovery from dry and wet for each site 
colnames(select.resis.recov)
s_resis.recov<-select.resis.recov%>%mutate(stability.facets1=paste(stability.facets, climate1, sep="_"))%>%
    filter(!values3%in%c("Inf", "-Inf"))%>%filter(!is.na(values3))%>%
    group_by(cutoff, site_code, trt, block, community.property, stability.facets1, variable.id)%>%
    summarise(values=mean(values3))

############################################################################################
######################calculate invariability during experimental years ####################
############################################################################################
# delete years at some sites that were not included for resistance and recovery 
stb<- d.select%>%group_by(variable.id) %>% 
  summarise_at(c("property.value"), list(N=~length(.), avg=~mean(.), sd=~sd(.), invariability=~mean(.)/sd(.)))%>%filter(N>=3)
range(stb$N)
####calculate detrended stability
## a linear term may be enough for most sites, a quadratic term may be ok as well 
sd.d<-c()
for(pl in unique(d.select$variable.id)){
  # pl<-"msla_2.us_2_Control"
  data<-d.select %>%filter(variable.id==pl)%>%filter(!is.na(property.value))
  if(nrow(data)<3) next
  mod<-lm(property.value~year_trt, data=data, na.action = na.omit)
  sd.lm<-sd(residuals(mod))
  mod1<-lm(property.value~year_trt+I(year_trt^2), data=data, na.action = na.omit)
  sd.q<-sd(residuals(mod1))
  temp<-data.frame(variable.id=pl, sd.lm=sd.lm, sd.q=sd.q)
  sd.d<-rbind(sd.d, temp)
}
## add detrended sd to stb data
stb2<-stb%>%merge(y=sd.d, by=c("variable.id"))%>% mutate(invariability.d=avg/sd.lm)%>%
  merge(d.select%>%select(site_code, block, trt, community.property, variable.id)%>%distinct(), by=c("variable.id"))%>%
   dplyr::select(site_code, trt, block, community.property, invariability, invariability.d)%>%
   pivot_longer(cols = c("invariability", "invariability.d"), names_to ="stability.facets1" , values_to ="values")

## calculate community dissimilarity for each plot over the experimental years 
inv.com<-c()
for(pl in unique(d.select.com$variable.id)){
  # pl<-"rook.uk_2_NPK_all_composition"
  
  ## delete never occurred species to speed up calculation
  data<-d.select.com%>%filter(variable.id==pl)%>%dplyr::select(variable.id, year_trt, standard_taxon, max_cover)%>%pivot_wider(names_from = standard_taxon, values_from = max_cover)%>%
    replace(is.na(.), 0)%>%mutate(variable.id=NULL, year_trt=NULL)%>%dplyr::select(which(!colSums(.)==0))
 
  dis<-beta.multi.abund(data, index.family = "bray")
  sim<-1-dis$beta.BRAY # higher values indicate higher invariability 
  sim.temp<-data.frame(variable.id=pl, sim=sim)
  inv.com<-rbind(inv.com, sim.temp)
}
inv.com1<-d.select.com%>%dplyr::select("variable.id", "site_code", "trt", "block", "community.property")%>%distinct()%>%
  merge(y=inv.com, by=c("variable.id"))%>%mutate(values=sim, stability.facets1="invariability")%>%
  dplyr::select(site_code, trt, block, community.property, stability.facets1, values)

# add all temporal invariability from three community aspects 
stb3<-stb2%>%bind_rows(inv.com1)
table(stb3$community.property)

############################################################################################
############################ add all stability facets together  ############################
############################################################################################
# using the sames sites for resistance and recovery (which differ using different cutoffs) for invariability
colnames(s_resis.recov)
colnames(stb3)

data.stability.facets<-c()
for (cut in  c(0.67, 1, 1.28, 1.5, 2)){
  # cut<-0.67
  temp.data.cut<- s_resis.recov%>%select(site_code, trt, block,community.property, stability.facets1,  values, cutoff)%>%
   filter(cutoff==cut)
 n.sites<-length(unique(subset(temp.data.cut, cutoff==cut)$site_code))## check how many sites were included for analyses
 all.stability<-stb3%>%filter(site_code%in%temp.data.cut$site_code)%>%mutate(cutoff=cut)%>%bind_rows(temp.data.cut)%>%mutate(n.sites=n.sites)         
 data.stability.facets<-rbind(data.stability.facets, all.stability)
  }
unique(data.stability.facets$n.sites)


#############################################################################################
############ nutrient addition effects on stability for all aspects and  facets  ############
#############################################################################################
colnames(data.stability.facets)
data.stability.facets.sub<-data.stability.facets%>%filter(!is.na(values))%>%filter(!values%in% c("Inf", "-Inf"))
all.data.l_1<-data.stability.facets.sub%>%mutate(variable.id=paste(cutoff, community.property, stability.facets1, sep="_"))%>% 
    mutate(values1=ifelse(community.property=="composition", values, log(values)))

# plot raw data for stability facets for different sites
(pp.trt.raw<-all.data.l_1%>%filter(cutoff %in% c(0.67))%>%filter(!(stability.facets1=="invariability" & community.property=="biomass"))%>%
    filter(!(stability.facets1=="invariability" & community.property=="richness"))%>%
    mutate(stability.facets2=ifelse(stability.facets1 %in% c("invariability", "invariability.d"), "invariability", stability.facets1))%>%
    ggplot(aes(stability.facets2, values1, shape=community.property, color=trt ))+ theme_cowplot(font_size = 20)+panel_border()+
    stat_summary(fun=mean, geom="point", size=2, position=position_dodge(width=0.7), alpha=0.6)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=position_dodge(width=0.7))+
    facet_wrap(~site_code, ncol=5)+
    scale_shape_manual(values = c(16, 18, 17))+
    guides(color = guide_legend(title = "Treatment", override.aes = list(shape = NA)),
           shape = guide_legend(title = "Commmunity aspects"))+
    theme(legend.position ="bottom", axis.text.x = element_text(angle = 90, vjust=0.7, hjust=1))+
    labs(x=NULL, y=NULL, color=NULL, shape=NULL))

# ggsave(pp.trt.raw, file="raw data for each stability facets based on cutoff 0.67.pdf", width = 21, height = 29.7, dpi=600)


# look at the data distribution 
# all.data.l_1%>%filter(cutoff==1.28)%>%ggplot()+geom_density(aes(x=values1, color=trt))+facet_wrap(~variable.id, scales = "free")

eff.trt.facets<-c()
for(i in unique(all.data.l_1$variable.id)){
  # i<-"1.28_richness_recovery_Dry"   
  data1<-all.data.l_1%>%filter(variable.id==i)%>%filter(!values1%in%c("Inf", "-Inf"))
  # record how many sites have these two variables 
  n.sites<-length(unique(data1$site_code))
  check.values<-unique(data1$values1)
  # ggplot(data1)+geom_histogram(aes(x=values1))+facet_wrap(~trt)
  mod<-lme(values1 ~ trt, random=~1|site_code/block,  data=data1)
  # plot(mod)
  t<-xtable(summary(mod)$tTable)
  rs<-rsquared(mod)
  get.sd.sites<-as.numeric(VarCorr(mod)[,"StdDev"][2])
  get.sd.block<-as.numeric(VarCorr(mod)[,"StdDev"][4])
  t1<-t%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, get.sd.block=get.sd.block,  number.sites=n.sites, terms=rownames(.),  variable.id=i)
  
  eff.trt.facets<-rbind(eff.trt.facets, t1) 
}

# save the full table 
eff.trt.facets.full<-all.data.l_1%>%ungroup()%>%dplyr::select(cutoff, community.property, stability.facets1, variable.id)%>%unique()%>%
  merge(eff.trt.facets, by=c("variable.id"))%>%
  dplyr::rename(p='p-value')%>%dplyr::select(cutoff, community.property, stability.facets1, terms, Value, Std.Error,  DF,  "t-value", "p", r2m, r2c, get.sd.block, get.sd.sites)%>%
  arrange(cutoff, community.property, stability.facets1)
table(eff.trt.facets.full$community.property)
eff.trt.facets.full[,5:13]<-round(eff.trt.facets.full[,5:13], 2)

eff.trt.facets.sig<-eff.trt.facets.full%>%filter(p<=0.05  & terms!="(Intercept)")
# write_xlsx(eff.trt.facets.sig, path="significant treatment effects on stability facets.xlsx", col_names = TRUE)
# write_xlsx(eff.trt.facets.full, path="full table of treatment effects on stability facets.xlsx", col_names = TRUE)

######### plot treatment effects on five stability facets of three community aspects ########
# delete intercept
eff.trt1<-eff.trt.facets.full%>%filter(terms!="(Intercept)")
## calculate 95% confidence intervals 
eff.trt1$lower<-eff.trt1$Value-1.96*eff.trt1$Std.Error
eff.trt1$upper<-eff.trt1$Value+1.96*eff.trt1$Std.Error
eff.trt1$trt1<-eff.trt1$terms
eff.trt1$trt<-gsub("trt", "", eff.trt1$trt1)
unique(eff.trt1$cutoff)
eff.trt4<-eff.trt1
eff.trt4$stability.facets1<-factor(eff.trt4$stability.facets1, levels = c("invariability.d", "invariability", "resistance_Dry", "resistance_Wet", "recovery_Dry", "recovery_Wet"))

##  plot all stability.facets together in one figure 
(pp.trt<-eff.trt4%>%ggplot(aes(stability.facets1, Value, color=community.property))+theme_bw(base_size = 20)+
    facet_wrap(~cutoff , nrow=4, scales = "free_y")+
    geom_point(position=pd, pch=19, size=5, alpha=0.6)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    geom_hline(yintercept = 0, linetype ="dotted") +
    guides(color=guide_legend(nrow=1))+
    theme(legend.position = c(0.8, 0.2), axis.text.x = element_text(angle = 15))+
    labs(x=NULL, y="Effect size of nutrient addition", color=NULL))

# to get the legend for pentagons (the effect size is wrong in this figure)
colour.crosswalk <- c("positive" = "black", "negative" = "red")
(pp.trt.all<-eff.trt4%>%filter((cutoff %in% c(1.28)))%>%mutate(Value.scaled=abs(Value)*10)%>%
    mutate(direction=ifelse(Value>=0, "positive", "negative"))%>%
    ggplot(aes(stability.facets1, Value.scaled,  colour = as.character(direction), size = (Value.scaled/10)))+theme_bw(base_size = 16)+
    facet_wrap(~community.property)+
    geom_point(position=pd, pch=15,alpha=0.6)+
    scale_colour_manual(values = colour.crosswalk) +
    labs(x=NULL, y="Effect size of nutrient addition", size="Effect size", color="Direction"))
legend<-get_legend(pp.trt.all)

# plotting using pentagons
source("functions to draw pentagen.R")
facet_coordinates<-read.csv("facet_coordinates.csv")
colnames(facet_coordinates)[1]<-"facet"
facet_coordinates$Facet<-c("invariability", "resistance_Dry", "resistance_Wet", "recovery_Dry", "recovery_Wet", "Nutrient.addition")

cut<-c(1.28, 0.67)
facet_data <- eff.trt4 %>% filter(cutoff %in% cut)%>% filter(!(stability.facets1=="invariability" & community.property=="biomass"))%>%
  filter(!(stability.facets1=="invariability" & community.property=="richness"))%>%
  mutate(stability.facets2=gsub(".d", "", stability.facets1))%>%
  mutate(stability.facets1= stability.facets2)%>%
  mutate(From_Facet="Nutrient.addition", To_Facet=stability.facets1, value=Value)%>%
  merge(y = facet_coordinates, by.x = "From_Facet", by.y = "Facet") %>%
  merge(y = facet_coordinates, by.x = "To_Facet", by.y = "Facet",
        suffixes = c("", "end")) %>%
  mutate(direction = ifelse(value >= 0, "positive", "negative"), 
         width = as.integer((abs(value) * 10)))%>%
  mutate(id=paste(cutoff, community.property, sep="_"), id1=paste(cutoff, sep="_"))

list_plots <- vector('list', length(unique(facet_data$id)))
for (i in unique(facet_data$id)){
    # i<-"1.28_all_richness"
    facet_data_temp<-facet_data%>%filter(id==i)%>% mutate(sig=case_when(p<=0.05 ~ "significant",  TRUE~"non-significant"))
    facet_data_temp$sig<-factor(facet_data_temp$sig, levels = c("non-significant", "significant"))
    fontsize <- 16 / .pt; size.npk<-20 / .pt
    
     list_plots[[i]]<-gg.pentagon(data = facet_data_temp)+ 
      geom_label(aes( x = 0, y = 0, label = "NPK"), fill="white", size = size.npk)
  }

for(i in c(1.28, 0.67)){ 
  # i<-1.28
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  # rename each element 
  names(list_plots_temp)<-gsub( "(.*)_(.*)", "\\2", names(list_plots_temp))
(pp.trt.effects<-plot_grid(list_plots_temp$biomass,
                           list_plots_temp$composition,
                           list_plots_temp$richness,
                           legend,
                           rel_widths = c(4.5,4.5,4.5,1.5),
          nrow=1, vjust = 5, hjust=-0.1, 
           labels = c("A (Biomass)", "B (Composition)", "C (Richness)"), label_fontface = "bold", label_size = 18))

 ggsave(pp.trt.effects, height=6, width=13.3, dpi=600, file=paste0("effects on stability facets for ", i , ".png"))
}

#############################################################################################
###calculate correlation between stability facets within biomass, composition, and richness within sites ##
#############################################################################################
# focus on the whole community and cutoff of 0.67 and 1.28 sd
data.stability.facets.relation<-data.stability.facets.sub%>%filter(cutoff %in% c(0.67, 1.28))%>%
  filter(!(stability.facets1=="invariability" & community.property=="biomass"))%>%
  filter(!(stability.facets1=="invariability" & community.property=="richness"))%>%
  mutate(stability.facets2=gsub(".d", "", stability.facets1))%>%
  mutate(stability.facets1= stability.facets2, stability.facets2=NULL, n.sites=NULL)%>%distinct()

all.data.l_1<-data.stability.facets.relation%>%
  mutate(variable.id=paste(cutoff, community.property, site_code,  trt, sep="_"))%>%
  pivot_wider(names_from=stability.facets1, values_from = values)

corr.stab.facet<-c()
for(id in unique(all.data.l_1$variable.id)){ 
  # id<-"1.28_biomass_all_bnch.us_NPK"
  data<-all.data.l_1%>%filter(variable.id==id)
  # calculate correlation between stability facets 
  inv_resis.d<-cor(data$invariability, data$resistance_Dry)
  inv_resis.w<-cor(data$invariability, data$resistance_Wet)
  inv_recov.d<-cor(data$invariability, data$recovery_Dry)
  inv_recov.w<-cor(data$invariability, data$recovery_Wet)
  resis.d_resis.w<-cor(data$resistance_Dry, data$resistance_Wet)
  resis.d_recov.d<-cor(data$resistance_Dry, data$recovery_Dry)
  resis.d_recov.w<-cor(data$resistance_Dry, data$recovery_Wet)
  resis.w_recov.d<-cor(data$resistance_Wet, data$recovery_Dry)
  resis.w_recov.w<-cor(data$resistance_Wet, data$recovery_Wet)
  recov.d_recov.w<-cor(data$recovery_Dry, data$recovery_Wet)
  # make a data frame 
  corr.temp<-data.frame(variable.id=id,
                        inv_resis.d=inv_resis.d, inv_resis.w=inv_resis.w, inv_recov.d=inv_recov.d, inv_recov.w=inv_recov.w, resis.d_resis.w=resis.d_resis.w,
                        resis.d_recov.d=resis.d_recov.d, resis.d_recov.w=resis.d_recov.w, resis.w_recov.d=resis.w_recov.d, resis.w_recov.w, recov.d_recov.w)
  corr.stab.facet<-rbind(corr.stab.facet, corr.temp)
}

# change data to long version 
corr.stab.facet.l<-all.data.l_1%>%ungroup()%>%select(cutoff, community.property, site_code,  trt, variable.id)%>%distinct()%>%
  merge(corr.stab.facet, by=c("variable.id"))%>%mutate(variable.id=NULL)%>% 
  pivot_longer(cols =inv_resis.d:recov.d_recov.w, names_to = "correlation.type")%>%filter(!is.na(value))%>%
  mutate(variable.id=paste(cutoff, community.property, correlation.type,  sep="_"))
hist(corr.stab.facet.l$value)
# check the data (correlation coefficients should not be |1| which is perfect fit or 0 which show no trends at all)
che.range<-corr.stab.facet.l%>%filter(value %in% c(0, -1, 1))
# look at the data distribution 
# corr.stab.facet.l%>%ggplot()+geom_density(aes(x=value, color=trt))+facet_wrap(~interaction(correlation.type, community.property), scales = "free")

out.corr<-c(); estimated.ci<-c()
for(id in unique(corr.stab.facet.l$variable.id)){ 
  # id<-"1.28_richness_recov.d_recov.w"
  data2<-corr.stab.facet.l%>%filter(variable.id==id)
  # use sites that have both control and NPK treatments
  check.trt<-data2%>%group_by(site_code)%>%summarise(N=length(value))%>%filter(N==2)
  data3<-data2%>%filter(site_code%in%check.trt$site_code)
  # record how many sites have these two variables 
  n.sites<-length(unique(data3$site_code))
  # record how many treatments 
  n.trt<-length(unique(data3$trt))
  if(n.sites<3|n.trt<2) next
  mod<-lme(value ~ trt,random=~1|site_code, data=data3)
  t<-xtable(summary(mod)$tTable)
  rs<-rsquared(mod)
  get.sd.sites<-as.numeric(VarCorr(mod)[,"StdDev"][1])
  t1<-t%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, number.sites=n.sites, terms=rownames(.),  variable.id=id)
  out.corr<-rbind(out.corr, t1)
  
  # estimated mean and 95% confidence ci
  emm<-emmeans(mod, specs = "trt")
  ci<-confint(emm)
  ci1<-ci%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, number.sites=n.sites,  variable.id=id)
  estimated.ci<-rbind(estimated.ci, ci1)
}

out.corr1<-corr.stab.facet.l%>%select(cutoff, community.property, correlation.type, variable.id)%>%distinct()%>%
  merge(out.corr, by=c("variable.id"))%>%mutate(variable.id=NULL)%>%mutate(across(4:11, \(x) round(x, 2)))%>%dplyr::rename(p='p-value')%>%
  arrange(cutoff,community.property,  terms, correlation.type)
table(out.corr1$community.property)
## sort the significant effects 
range(out.corr1$number.sites)
lme.relationship.facet.sig<-out.corr1%>%filter(p<=0.05)

# save the tables
# write_xlsx(lme.relationship.facet.sig, path="significant treatment effects on relationships among stability facets.xlsx", col_names = TRUE)
# write_xlsx(out.corr1, path="full table of treatment effects on relationships among stability facets.xlsx", col_names = TRUE)

################ plot treatment effects on correlation among stability facets  ##############
colnames(estimated.ci)
estimated.ci1<-corr.stab.facet.l%>%select(cutoff, community.property, correlation.type, variable.id)%>%distinct()%>%
  merge(estimated.ci, by=c("variable.id"))%>%mutate(across(9:13,\(x) round(x, 2)))%>%mutate(variable.id=NULL, direction=ifelse(emmean>=0, "positive", "negative"))
# write_xlsx(estimated.ci1%>%mutate(number.sites=NULL, direction=NULL, life.form=NULL), path="estimated relationships among stability facets using emmeans.xlsx", col_names = TRUE)

# get the legend for pentagons (the effect size is wrong in this figure)
(pp.cor<-estimated.ci1%>%filter(cutoff==1.28)%>%mutate(Value.scaled=abs(emmean)*10)%>%
    ggplot(aes(correlation.type, Value.scaled,  colour = factor(direction), size = (Value.scaled/10)))+theme_bw(base_size = 16)+
    facet_grid(trt~community.property, scales = "free_y")+
    geom_point(position=pd, pch=15, alpha=0.6)+
    #geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=pd) +
    geom_hline(yintercept = 0, linetype ="dotted") +
    scale_colour_manual(values = colour.crosswalk) +
    theme(axis.text.x = element_text(angle=20, hjust=1))+
    #theme(legend.position = c(0.8, 0.2))+
    labs(x=NULL, y="Correlation coeffecient", size="Correlation\n coefficient", color="Direction"))
legend.cor<-get_legend(pp.cor)

# plotting using pentagons
facet_data_cor_stab <- estimated.ci1%>%filter(cutoff %in% c(0.67, 1.28))%>%
  mutate(From_Facet=case_when(grepl("inv_", correlation.type)~"invariability", 
                              grepl("resis.d_", correlation.type)~"resistance_Dry", 
                              grepl("resis.w_", correlation.type)~"resistance_Wet",
                              grepl("recov.d_", correlation.type)~"recovery_Dry",
                              grepl("recov.w_", correlation.type)~"recovery_Wet"),
         To_Facet=case_when(grepl("_inv", correlation.type)~"invariability", 
                             grepl("_resis.d", correlation.type)~"resistance_Dry", 
                             grepl("_resis.w", correlation.type)~"resistance_Wet",
                             grepl("_recov.d", correlation.type)~"recovery_Dry",
                             grepl("_recov.w", correlation.type)~"recovery_Wet"),
         value=emmean)%>%
  merge(y = facet_coordinates, by.x = "From_Facet", by.y = "Facet") %>%
  merge(y = facet_coordinates, by.x = "To_Facet", by.y = "Facet",
        suffixes = c("", "end")) %>%
  mutate(direction = ifelse(value >= 0, "positive", "negative"),
         width = as.integer((abs(value) * 10) ), 
        sig=ifelse(((upper.CL/lower.CL)>0|upper.CL==0|lower.CL==0), "significant", "non-significant"))%>%
  mutate(id=paste(cutoff,  community.property, trt, sep="_"), id1=paste(cutoff, sep="_"))

list_plots <- vector('list', length(unique(facet_data_cor_stab$id)))

for (i in unique(facet_data_cor_stab$id)){
  # i<-"1.28_richness_NPK"
  facet_data_cor_stab_temp<-facet_data_cor_stab%>%filter(id==i)
  list_plots[[i]]<-gg.pentagon(data = facet_data_cor_stab_temp)
}

for(i in unique(facet_data_cor_stab$id1)){ 
  # i<-"1.28"
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  # rename each element 
  names(list_plots_temp)<-gsub( "(.*)_(.*)_(.*)", "\\2\\3", names(list_plots_temp))
  (pp.cor<-plot_grid(list_plots_temp$biomassControl,
                             list_plots_temp$compositionControl,
                             list_plots_temp$richnessControl,
                             
                             list_plots_temp$biomassNPK,
                              list_plots_temp$compositionNPK,
                             list_plots_temp$richnessNPK,
                              nrow=2))
  (pp.cor1<-plot_grid(pp.cor, legend.cor, rel_widths = c(8.6, 1.4)))
  (pp.cor2<-ggpubr::annotate_figure(pp.cor1, top=text_grob("Biomass                                        Composition                                    Richness                          ", size=18, face = "bold") ,
                                            left=text_grob("NPK                               Control", rot=90, face = "bold", size=18)))
  
  ggsave(pp.cor2, height=6.64, width=13.3, dpi=600, file=paste0("relationships between stability facets ", i , ".png"))
}

#############################################################################################
##correlation between biomass and composition, biomass and richness, composition and richness within sites#
#############################################################################################
all.data.l_1<-data.stability.facets.relation%>%
  mutate(variable.id=paste(cutoff, stability.facets1, site_code,  trt, sep="_"))%>%
  pivot_wider(names_from=community.property, values_from = values)

corr.asp<-c()
for(id in unique(all.data.l_1$variable.id)){ 
  # id<-"1.5_invariability_all_saline.us_NPK"
  data1<-all.data.l_1%>%filter(variable.id==id)
  # calculate correlation between community aspects
  bio_com<-cor(data1$biomass, data1$composition)
  bio_div<-cor(data1$biomass, data1$richness)
  com_div<-cor(data1$composition, data1$richness)
  # make a data frame 
  corr.temp<-data.frame(variable.id=id, bio_com=bio_com, bio_div=bio_div, com_div=com_div)
  corr.asp<-rbind(corr.asp, corr.temp)
}
# change data to long version 
corr.asp.l<-all.data.l_1%>%ungroup()%>%select(cutoff, stability.facets1, site_code,  trt, variable.id)%>%distinct()%>%
  merge(corr.asp, by=c("variable.id"))%>%mutate(variable.id=NULL)%>% 
  pivot_longer(cols =c("bio_com", "bio_div", "com_div"), names_to = "correlation.type")%>%filter(!is.na(value))%>%
  mutate(variable.id=paste(cutoff, stability.facets1, correlation.type,  sep="_"))

## look at treatment effects
trt.corr.asp<-c(); estimated.ci.asp<-c()

for(id in unique(corr.asp.l$variable.id)){ 
  # id<-"0.67_recovery_Dry_bio_div_all"
  data2<-corr.asp.l%>%filter(variable.id==id)
  # use sites that have both control and NPK treatments
  check.trt<-data2%>%group_by(site_code)%>%summarise(N=length(value))%>%filter(N==2)
  data3<-data2%>%filter(site_code%in%check.trt$site_code)
  # record how many sites have these two variables 
  n.sites<-length(unique(data3$site_code))
  # record how many treatments 
  n.trt<-length(unique(data3$trt))
  if(n.sites<3|n.trt<2) next
  mod<-lme(value ~ trt,random=~1|site_code, data=data3)
  t<-xtable(summary(mod)$tTable)
  rs<-rsquared(mod)
  get.sd.sites<-as.numeric(VarCorr(mod)[,"StdDev"][1])
  t1<-t%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, number.sites=n.sites, terms=rownames(.),  variable.id=id)
  trt.corr.asp<-rbind(trt.corr.asp, t1)
  
  # estimated mean and 95% confidence ci
  emm<-emmeans(mod, specs = "trt")
  ci<-confint(emm)
  ci1<-ci%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, get.sd.sites=get.sd.sites, number.sites=n.sites,  variable.id=id)
  estimated.ci.asp<-rbind(estimated.ci.asp, ci1)
}
trt.corr.asp1<-corr.asp.l%>%select(cutoff, stability.facets1, correlation.type, variable.id)%>%distinct()%>%
  merge(trt.corr.asp, by=c("variable.id"))%>%mutate(variable.id=NULL)%>%mutate(across(4:11,\(x) round(x, 2)))%>%dplyr::rename(p='p-value')%>%
  arrange(cutoff, stability.facets1,  terms, correlation.type)
table(trt.corr.asp1$stability.facets1)

# sort the significant effects 
range(trt.corr.asp1$number.sites)
lme.relationship.aspect.sig<-trt.corr.asp1%>%filter(p<=0.05)
# save the full tables
# write_xlsx(lme.relationship.aspect.sig, path="significant treatment effects on relationships among community aspects.xlsx", col_names = TRUE)
# write_xlsx(trt.corr.asp1, path="full table of treatment effects on relationships among community aspects.xlsx", col_names = TRUE)

################ plot treatment effects on correlation among community aspects ##############
colnames(estimated.ci.asp)
unique(estimated.ci.asp$correlation.type)
estimated.ci.asp1<-corr.asp.l%>%select(cutoff, stability.facets1 , correlation.type, variable.id)%>%distinct()%>%
  merge(estimated.ci.asp, by=c("variable.id"))%>%mutate(across(6:13, \(x) round(x, 2)))%>%mutate(variable.id=NULL)
# write_xlsx(estimated.ci.asp1%>%mutate(number.sites=NULL), path="estimated relationships among community aspects using emmeans.xlsx", col_names = TRUE)

(pp.cor.asp<-estimated.ci.asp1%>%ggplot(aes(correlation.type, emmean, color=trt))+theme_bw(base_size = 20)+
    facet_grid(cutoff~stability.facets1, scales = "free_y")+
    geom_point(position=pd, size=5, pch=19, alpha=0.6)+
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=pd) +
     geom_hline(yintercept = 0, linetype ="dotted") +
     theme(axis.text.x = element_text(angle=20, hjust=1))+
    labs(x=NULL, y="Correlation coeffecient", color=NULL))

# plotting using triangle
aspect_coordinates<-facet_coordinates%>%filter(facet %in% c("Facet 1", "Facet 4", "Facet 5"))%>%
  mutate(Aspect=c("biomass",  "richness", "composition"), facet=NULL, Facet=NULL,
         x=c(0, -0.866, 0.866), y=c(1, -0.5, -0.5))

Aspect_data_cor <- estimated.ci.asp1%>%
  mutate(From_Aspect=case_when(grepl("bio_", correlation.type)~"biomass", 
                              grepl("div_", correlation.type)~"richness", 
                              grepl("com_", correlation.type)~"composition"),
                            
         To_Aspect=case_when(grepl(".bio", correlation.type)~"biomass", 
                             grepl("_div", correlation.type)~"richness", 
                             grepl("_com", correlation.type)~"composition"),
         value=emmean)%>%
  merge(y = aspect_coordinates, by.x = "From_Aspect", by.y = "Aspect") %>%
  merge(y = aspect_coordinates, by.x = "To_Aspect", by.y = "Aspect",
        suffixes = c("", "end")) %>%
  mutate(direction = ifelse(value >= 0, "positive", "negative"),
         width = as.integer((abs(value) * 10)),
         sig=ifelse(((upper.CL/lower.CL)>0|upper.CL==0|lower.CL==0), "significant", "non-significant"))%>%
  mutate(id=paste(cutoff,  stability.facets1, trt, sep="_"), id1=paste(cutoff, sep="_"))

list_plots <- vector('list', length(unique(Aspect_data_cor$id)))
for (i in unique(Aspect_data_cor$id)){
  # i<-"0.67_all_recovery_Dry_Control"
  Aspect_data_cor_temp<-Aspect_data_cor%>%filter(id==i)
  list_plots[[i]]<-gg.triangle(data = Aspect_data_cor_temp)
}

for(i in unique(Aspect_data_cor$id1)){ 
  # i<-"0.67"
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  # rename each element 
  names(list_plots_temp)<-gsub("^[^_]*_", "", names(list_plots_temp))
  (pp.asp<-plot_grid(list_plots_temp$invariability_Control,
                             list_plots_temp$resistance_Dry_Control,
                             list_plots_temp$resistance_Wet_Control,
                             list_plots_temp$recovery_Dry_Control,
                             list_plots_temp$recovery_Wet_Control,
                             
                             list_plots_temp$invariability_NPK,
                             list_plots_temp$resistance_Dry_NPK,
                             list_plots_temp$resistance_Wet_NPK,
                             list_plots_temp$recovery_Dry_NPK,
                             list_plots_temp$recovery_Wet_NPK,
                             NA, NA, NA,NA,NA,
                     rel_heights = c(4,4,2),
                             nrow=3))
  (pp.asp1<-plot_grid(pp.asp, legend.cor, rel_widths = c(8.6, 1.4)))
  (pp.asp2<-ggpubr::annotate_figure(pp.asp1, 
                                    top=text_grob("  Invariability            Resistance_Dry            Resistance_Wet            Recovery_Dry            Recovery_Wet                           ", size=16, face = "bold") ,
                                            left=text_grob("                     NPK                              Control", rot=90, face = "bold", size=16)))
  
  ggsave(pp.asp2, height=6.64, width=13.3, dpi=600, file=paste0("relationships between community Aspects ", i , ".png"))
}

###########################################################################################
##### show biomass at each site during different growing seasons using cutoff of 0.67sd #######
###########################################################################################

climate.extreme.cut0.67<-s_dd7%>%merge(d8%>%select(site_code, block, trt, year_trt, live_mass, richness)%>%distinct(), by=c("site_code", "year_trt"))

unique(climate.extreme.cut0.67$climate)
colour.crosswalk <- c("Extreme dry"="#E69F00", "Moderate dry"="#F0E442",  "Normal"="#000000",  "Moderate wet"="#56B4E9",  "Extreme wet"="#0072B2")
colour.crosswalk1 <- c("Dry"="#E69F00", "Normal"="#000000",   "Wet"="#0072B2")

climate.extreme.cut0.67_each.site_avg<-climate.extreme.cut0.67%>%filter(climate=="Normal")%>%group_by(site_code, trt)%>%
  dplyr::summarise(avg_bio=mean(live_mass), avg_rich=mean(richness))
(pp.each.site.bio<-ggplot(climate.extreme.cut0.67, aes(year_trt, live_mass, color=climate, shape=trt))+theme_cowplot(font_size = 25)+panel_border()+
    facet_wrap(~site_code, scales="free_y", ncol=5)+
    stat_summary(fun=mean, geom="point", size=5, position=pd, alpha=0.6)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(climate.extreme.cut0.67_each.site_avg, mapping=aes(yintercept =avg_bio, linetype=trt), alpha=0.6)+
    scale_colour_manual(values = colour.crosswalk) +
    scale_x_continuous(expand=c(0,0), limits = c(0.5,15), breaks = seq(0,15,2))+
    theme(axis.text.y = element_text(angle=30))+
    labs(x="Years after nutrient addition", y="Aboveground biomass", color="Climate\n extremes", shape="Treatments", linetype="Treatments"))
# ggsave("pp.each.site.bio.pdf", width = 21, height = 29.7, dpi=600)

# focus on site look.us 
look_each.block_avg<-climate.extreme.cut0.67%>%filter(site_code=="look.us" & climate=="Normal")%>%
  group_by(site_code, block, trt)%>% dplyr::summarise(avg_bio=mean(live_mass), avg_rich=mean(richness))
 
(pp.look.bio<-climate.extreme.cut0.67 %>%filter(site_code=="look.us")%>%
    #mutate(climate0=ifelse(grepl("wet", climate), "Wet", climate), climate1=ifelse(grepl("dry", climate0), "Dry", climate0), climate2=ifelse(grepl("dry", climate1), "Dry", climate1))%>%
    mutate(climate1=case_when(grepl("dry", climate)~"Dry", grepl("wet", climate)~"Wet", TRUE~"Normal"))%>%
    select(block, trt, year_trt, live_mass, climate1, spei)%>%
    pivot_longer(cols = c("live_mass", "spei"))%>%mutate(block1=ifelse(name=="live_mass",  block, "SPEI"), trt1=ifelse(name=="live_mass",  trt, "SPEI"))%>%
    merge(y=look_each.block_avg, by=c("block", "trt"))%>%mutate(avg_bio1=ifelse(name=="live_mass", avg_bio, NA))%>%
    merge(select.resis.recov%>%ungroup()%>%
            select(cutoff, site_code, year_trt, stability.facets, spei1)%>%filter(site_code=="look.us" & cutoff==0.67)%>%distinct(), by=c("site_code", "year_trt"), all.x=T)%>%
    mutate(stability.facets1=ifelse(trt1=="SPEI", "SPEI", stability.facets), stability.facets2=replace_na(stability.facets1, "Not used"), stability.facets3=ifelse(climate1=="Normal", "Normal level", stability.facets2))%>%
    ggplot(aes(year_trt, value))+theme_cowplot(font_size = 20)+panel_border()+
    facet_wrap(~block1, nrow=4, scale="free_y")+
    geom_jitter(size=5, position=pd, aes(color=climate1, shape=trt1, alpha=stability.facets3))+geom_line(aes(group=trt1), position=pd, alpha=0.2)+
    geom_hline(aes(yintercept =avg_bio1, linetype=trt))+
    scale_colour_manual(values = colour.crosswalk1) +
    scale_alpha_manual(values = c(1, 0.1, 0.3, 0.6, 1))+
    scale_x_continuous(expand=c(0,0), limits = c(1.5,14.5), breaks = seq(0,15,1))+
    theme(axis.text.y = element_text(angle=30))+
    labs(x="Years after nutrient addition", y="SPEI                            Aboveground biomass   ", color="Climate\n extremes", shape=NULL, linetype="Treatments", alpha="Years used"))
# ggsave(pp.look.bio, dpi=600,width=13.3, height=6.64,  file="above ground biomass at site look.us.png")

# show an example how resistance and recovery is calculated using block 1 at site look.us 
(look.year.resis<-select.years.extremes%>%filter( site_code=="look.us" & year_trt >0 )%>%filter(stability.facets=="resistance"))

cut<-0.67; max.year.look<-d8%>%select(site_code, year_trt)%>%group_by(site_code)%>%summarise(n.max=max(year_trt))%>%filter(site_code=="look.us") 

look.year.recov<-s_dd7%>%filter(site_code=="look.us" & year_trt >0)%>%select(site_code, year_trt, climate, spei)%>%distinct()%>%filter(climate!="Normal")%>%
  mutate(spei1=as.numeric(spei))%>%
  arrange(site_code, year_trt)%>%group_by(site_code)%>%
  mutate(climate.num=case_when((spei1>=cut)~1,
                               (spei1<=- cut)~ -1,
                               TRUE~0), consecutive=lead(year_trt)-year_trt, climate.type=lead(climate.num)-climate.num, consecutive1=ifelse(is.na(consecutive), 999, consecutive))%>%
  mutate(year.recov=ifelse((consecutive1==1 & climate.type==0), NA, year_trt))%>%
  merge(select.years.extremes%>%filter(stability.facets=="recovery" & year_trt!=0), by=c("site_code", "year_trt", "climate"))%>%filter(year_trt!=max.year.look$n.max)%>%filter(!is.na(year.recov))

look.rr<-recov%>%filter( site_code=="look.us" & block==1)%>%mutate(climate1=case_when(climate %in% c("Extreme dry", "Moderate dry")~"Dry",  climate %in% c("Extreme wet", "Moderate wet")~"Wet", TRUE ~ "Normal"))%>%
  mutate(resistance1=ifelse(year_trt %in% look.year.resis$year_trt, resistance, NA), recovery1=ifelse(year_trt %in% look.year.recov$year.recov, recovery, NA))%>%
  group_by(community.property, trt, climate1)%>%mutate(resistance.avg=mean(resistance1, na.rm=T), recovery.avg=mean(recovery1, na.rm=T))%>%
  select(community.property, year_trt, climate1, spei, trt, property.value, avg.property, deviation.property, resistance, recovery, resistance1, recovery1, resistance.avg, recovery.avg)%>%
  arrange(community.property, trt)

colnames(recov_com)
look.rr.com<-rr.com1%>% filter(site_code=="look.us" & block==1)%>%
  mutate(climate1=case_when(climate %in% c("Extreme dry", "Moderate dry")~"Dry",  climate %in% c("Extreme wet", "Moderate wet")~"Wet", TRUE ~ "Normal"))%>%
  mutate(recovery=lead(sim)/sim, consecutive=lead(year_trt)- year_trt, recovery1=ifelse(consecutive==1, recovery, NA))%>%
   mutate(community.property="composition")%>%ungroup()%>%
   mutate(resistance1=ifelse(year_trt %in% look.year.resis$year_trt, resistance, NA), recovery1=ifelse(year_trt %in% look.year.recov$year.recov, recovery, NA))%>%
   group_by(community.property, trt, climate1)%>%mutate(resistance.avg=mean(resistance1, na.rm=T), recovery.avg=mean(recovery1, na.rm=T))%>%
   select(community.property, year_trt, climate1, spei, trt, sim, resistance, recovery,resistance1, recovery1, resistance.avg, recovery.avg)%>% 
   arrange(community.property, year_trt, trt)

# write.csv(look.rr, file="an example for calculating resistance and recovery for biomass and richness using site look.us.csv")
# write.csv(look.rr.com, file="an example for calculating resistance and recovery for community composition using site look.us.csv")

###########################################################################################
# average biomass, richness, composition during and one year after climate events using cutoff of 0.67sd #
###########################################################################################

bio.rich<- data.bio.rich%>%mutate(community.property=gsub("live_mass", "biomass", community.property))%>%
  select("site_code", "trt", "block", "year_trt", "community.property", "property.value")%>%distinct()

climate.extreme.cut0.67<-select.years.extremes%>%filter(stability.facets=="resistance")%>%
  merge(bio.rich, by=c("site_code", "year_trt"))%>%mutate(stability.facets=NULL)

## indicate the mean biomass in normal years for both treatments 
with(climate.extreme.cut0.67, tapply(property.value, list(community.property, trt, climate), mean))
## select biomass During 
climate.extreme.cut0.67.u1<-climate.extreme.cut0.67%>%mutate(events="During")

## summarize biomass and richness one year after 
## add biomass During 
## calculate biomass difference from the normal level for each site each treatment 
## add average biomass under different treatments different blocks 
climate.extreme.cut0.67_block_avg<-d.select%>%filter(climate!="Normal")%>%filter(community.property %in%c("biomass", "richness"))%>%
  select("site_code", "trt", "block", "year_trt", "climate", "community.property", "property.value")%>%distinct()%>%
  group_by(site_code, community.property, block, trt)%>%
  dplyr::summarise(avg_property=mean(property.value))
# calculate average across sites
climate.extreme.cut0.67_t<-climate.extreme.cut0.67_block_avg%>%group_by(community.property, trt)%>%
  # summarise(avg=mean(avg_property), N=length(avg_property), sd=sd(avg_property), se=sd/sqrt(N))%>%
  # mutate(upper=avg+se, lower=avg-se)%>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(.$avg_property))))
t.c.bio<-climate.extreme.cut0.67_t%>%filter(community.property=="biomass" & trt=="Control")%>%ungroup()%>%select(Mean)%>%distinct()
t.n.bio<-climate.extreme.cut0.67_t%>%filter(community.property=="biomass" & trt=="NPK")%>%ungroup()%>%select(Mean)%>%distinct()
(t.n.bio -t.c.bio)/t.c.bio
t.c.rich<-climate.extreme.cut0.67_t%>%filter(community.property=="richness" & trt=="Control")%>%ungroup()%>%select(Mean)%>%distinct()
t.n.rich<-climate.extreme.cut0.67_t%>%filter(community.property=="richness" & trt=="NPK")%>%ungroup()%>%select(Mean)%>%distinct()
(t.n.rich -t.c.rich)/t.c.rich

# when two climate extremes happen consecutively, recovery only calculate for the later year
# climate events happen in the last year of a site cannot calculate recovery because no data the year after
# sites with only one normal growing season cannot calculate recovery, because |y(e+1)-y(n)| cannot be calculated 
max.year<-d8%>%select(site_code, year_trt)%>%group_by(site_code)%>%summarise(n.max=max(year_trt))
cut<-0.67
years.recovery<-s_dd7%>%select(site_code, year_trt, climate, spei)%>%distinct()%>%filter(climate!="Normal")%>%
  mutate(spei1=as.numeric(spei), climate1=case_when((spei1>=cut)~"Wet",
                                                    (spei1<=- cut)~"Dry",
                                                    TRUE~"between.normal.and.extreme"))%>%
  arrange(site_code, year_trt)%>%group_by(site_code)%>%
  mutate(climate.num=case_when((spei1>=cut)~1,
                               (spei1<=- cut)~ -1,
                               TRUE~0), consecutive=lead(year_trt)-year_trt, climate.type=lead(climate.num)-climate.num, consecutive1=ifelse(is.na(consecutive), 999, consecutive))%>%
  mutate(year.recov=ifelse((consecutive1==1 & climate.type==0), NA, year_trt))%>%
  merge(select.years.extremes%>%filter(stability.facets=="recovery" & year_trt!=0), by=c("site_code", "year_trt", "climate"), all.y=T)%>%
  filter(!is.na(year.recov))%>%merge(max.year, by=c("site_code"))%>%filter(year_trt!=n.max)
length(unique(years.recovery$site_code))

# During last extreme before a normal season
last.extreme.for.recovery<-years.recovery%>%select(site_code, year_trt, climate)%>%distinct()%>%
  arrange(site_code, year_trt)%>%
  merge(y=climate.extreme.cut0.67[,c("site_code", "trt", "year_trt", "block", "community.property", "property.value")], by=c("site_code",  "year_trt"))%>%distinct()%>%
  mutate(events="During last extreme")

climate.extreme.cut0.67.au<-years.recovery%>%select(site_code, year_trt, climate)%>%distinct()%>%
  arrange(site_code, year_trt)%>%
  mutate(year_trt1=year_trt+1)%>%
  mutate(year_trt=year_trt1, year_trt1=NULL)%>%
  merge(y=bio.rich[,c("site_code", "trt", "year_trt", "block","community.property", "property.value")], by=c("site_code",  "year_trt"))%>%distinct()%>%
  mutate(events="One year after")%>%bind_rows(climate.extreme.cut0.67.u1)%>%
  bind_rows(last.extreme.for.recovery[,c("site_code", "year_trt", "climate", "trt", "block", "community.property", "property.value", "events")])%>%
  merge(climate.extreme.cut0.67_block_avg, by=c("community.property", "site_code", "block", "trt"))%>%mutate(dif_property=(property.value-avg_property), dev_property=abs(property.value-avg_property))%>%
  mutate(property.value=NULL, avg_property=NULL)%>%
  pivot_longer(cols=c("dif_property", "dev_property"))%>%arrange(name, trt, events)%>%
  mutate(reference.line=ifelse(name=="dev_property", NA, 0), climate1=case_when(grepl("dry", climate)~"Dry", grepl("wet", climate)~"Wet", TRUE~"Normal"))

# relevel the variables 
climate.extreme.cut0.67.au$events<-factor(climate.extreme.cut0.67.au$events, levels = c("During", "During last extreme",  "One year after"))
## biomass in each climate event average across sites During and one year after
(pp.u.o<-climate.extreme.cut0.67.au%>%filter(community.property=="biomass")%>%
    mutate(name0=ifelse(name=="dif_property", "Change in biomass", "Deviation in biomass"), name1=as.factor(name0), name1=relevel(name1, ref="Change in biomass"))%>%
    ggplot(aes(events, value, color=trt))+theme_cowplot(font_size = 20)+panel_border()+
    facet_grid(name1~climate1, scale="free_y", space="free_x")+
    stat_summary(fun=mean, geom="point", size=5, position=pd)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(aes(yintercept=reference.line), linetype="dashed", linewidth=1, alpha=0.6)+
    theme(axis.text.x = element_text(angle=35, hjust=1))+
    labs(x=NULL, y=NULL, color=NULL))
# ggsave(pp.u.o, dpi=600, width=13.3, height=6.64, file="aboveground biomass During and one year after.png")

(pp.rich.u.o<-climate.extreme.cut0.67.au%>%filter(community.property=="richness")%>%
    mutate(name0=ifelse(name=="dif_property", "Change in richness", "Deviation in richness"), name1=as.factor(name0), name1=relevel(name1, ref="Change in richness"))%>%
    ggplot( aes(events, value, color=trt))+theme_cowplot(font_size = 20)+panel_border()+
    facet_grid(name1~climate1, scale="free_y", space="free_x")+
    stat_summary(fun=mean, geom="point", size=5, position=pd)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(aes(yintercept=reference.line), linetype="dashed", linewidth=1, alpha=0.6)+
    theme(axis.text.x = element_text(angle=35, hjust=1))+
    labs(x=NULL, y=NULL, color=NULL))
# ggsave(pp.rich.u.o, dpi=600, width=13.3, height=6.64, file="Species richness During and one year after.png")

###########################################################################################
#######compare community dissimilarity during and one year after climate events ###########
###########################################################################################
composition<-d8%>%
  merge(s_dd7[,c("site_code", "year_trt", "climate")], by=c("site_code", "year_trt"))%>%
  mutate(plot.id=paste(site_code, block, trt, sep="_"))

## calculate mean abundance for each species in normal years 
s_n<-composition%>%filter(climate=="Normal")%>%group_by(plot.id, standard_taxon)%>%summarise(cover.avg=mean(max_cover))

## calculate community dissimilarity for each plot over normal growing seasons 
inv.com.normal<-c()
for(pl in unique(composition$plot.id)){
  data<-composition%>%filter(plot.id==pl & climate=="Normal")%>%select(plot.id, year_trt, standard_taxon, max_cover)%>%pivot_wider(names_from = standard_taxon, values_from = max_cover)%>%
    replace(is.na(.), 0)%>%mutate(plot.id=NULL, year_trt=NULL)%>%select(which(!colSums(.)==0))## delete never occurred species to speed up calculation
  if(nrow(data)<2) next
  dis<-beta.multi.abund(data, index.family = "bray")
  sim<-1-dis$beta.BRAY # higher values indicate higher invariability 
  sim.temp<-data.frame(plot.id=pl, sim=sim)
  inv.com.normal<-rbind(inv.com.normal, sim.temp)
}
inv.com.normal1<-composition%>%select("plot.id", "site_code", "trt", "block")%>%distinct()%>%
  merge(y=inv.com.normal, by=c("plot.id"))%>%mutate(climate="Normal", events="During", plot.id=NULL) # higher values indicate higher resistance

## calculate resistance and recovery
rr.com<-c()
for(pl in unique(composition$plot.id)){
  # pl<-"valm.ch_2_Control"; yr<-2
  ref.data<-s_n%>%filter(plot.id==pl)%>%dplyr::rename(max_cover=cover.avg)%>%mutate(reference.com="yes")
  data<-composition%>%filter(plot.id==pl)
  sim.t<-c()
  for(yr in unique(data$year_trt)){
    data1<-data%>%filter(year_trt==yr)%>%select(plot.id, standard_taxon, max_cover)%>%mutate(reference.com="no")
    ## add two communities together 
    t.com<-data1%>%bind_rows(ref.data)%>%pivot_wider(names_from = standard_taxon, values_from = max_cover)%>%
      replace(is.na(.), 0)%>%mutate(plot.id=NULL, reference.com=NULL)%>%select(which(!colSums(.)==0))
    if(nrow(t.com)<2) next
    dis<-beta.pair.abund(t.com, index.family = "bray")
    sim<-as.numeric(1-dis$beta.bray)
    sim.temp<-data.frame(plot.id=pl, year_trt=yr, sim=sim)
    sim.t<-rbind(sim.t, sim.temp)
  }
  rr.com<-rbind(rr.com, sim.t)
}
rr.com1<-composition%>%select("plot.id", "site_code", "trt", "block", "year_trt", "climate")%>%distinct()%>%
  merge(y=rr.com, by=c("plot.id", "year_trt")) # higher values indicate higher resistance
composition.u1<-rr.com1%>%filter(climate!="Normal")%>%mutate(events="During", plot.id=NULL) # similarity during climate events 

# During last extreme before a normal season
last.extreme.for.recovery.composition<-years.recovery%>%select(site_code, year_trt, climate)%>%distinct()%>%
  arrange(site_code, year_trt)%>%
  merge(y=rr.com1[,c("site_code", "trt", "year_trt", "block", "sim")], by=c("site_code",  "year_trt"))%>%distinct()%>%
  mutate(events="During last extreme")

composition.au<-years.recovery%>%select(site_code, year_trt, climate)%>%distinct()%>%
  arrange(site_code, year_trt)%>%
  mutate(year_trt1=year_trt+1)%>%
  mutate(year_trt=year_trt1, year_trt1=NULL)%>%
  merge(y=rr.com1[,c("site_code", "trt", "year_trt", "block", "sim")], by=c("site_code",  "year_trt"))%>%distinct()%>%
  mutate(events="One year after")%>%bind_rows(composition.u1)%>%bind_rows(last.extreme.for.recovery.composition)%>%
  bind_rows(inv.com.normal1)%>%  mutate(climate1=ifelse(grepl("wet", climate), "Wet", ifelse(grepl("dry", climate), "Dry", "Normal")))
composition.au$climate1<-factor(composition.au$climate1, levels = c("Normal", "Dry", "Wet"))
composition.au$events<-factor(composition.au$events, levels = c("During", "During last extreme",  "One year after"))

(pp.composition.u.o<-ggplot(composition.au, aes(events, sim, color=trt))+
    theme_cowplot(font_size = 20)+panel_border()+
    theme(axis.text.x = element_text(angle=35, hjust=1))+
    facet_wrap(~climate1, scale="free_x")+
    stat_summary(fun=mean, geom="point", size=5, position=pd)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x=NULL, y="Community similarity", color=NULL))
# ggsave(pp.composition.u.o, dpi=600, width=13.3, height=6.64, file="composition similarity During and one year after.png")

###########################################################################################
########################### sort the information for sites used ###########################
###########################################################################################
exp.yr<-s_dd7%>%filter(year_trt>0)%>%group_by(site_code)%>%summarise_at("year_trt", list(min, max))%>%
  filter(!fn2<4)%>%mutate(years.used=paste(fn1, fn2, sep="-"))%>%
  mutate(years.used1=ifelse(site_code=="bnch.us", "1,3-11", ifelse(site_code=="mcla.us", "1,3-12", years.used)))
# read the data of water balance during experimental years used 
wat.bal1<-avg.bal1%>%mutate(growing.season=paste(gs_start_month, harvest_month, sep="-"))

## indicate sites that only have dry or wet events or both during experimental years
## here excluding sites or years not used because of missing biomass
dw<-c()
for(s in unique(s_dd7$site_code)){
  # s<-"cbgb.us"
  data<-s_dd7%>%filter(climate!="Normal")%>% mutate(climate1=ifelse(grepl("wet", climate), "Wet", "Dry"))%>%
    filter(site_code==s)%>%select(site_code, climate1, year_trt)%>%distinct()
  if(nrow(data)==0){
    next
  }else{
    N<-length(unique(data$climate1))
    data$group.site<-ifelse(N==2, "dry.and.wet", "dry.or.wet")
    dw<-rbind(dw, data)
  }
}

dw$climate2<-"Dry and Wet growing seasons"
dw$climate2[dw$climate1=="Wet" & dw$group.site=="dry.or.wet"]<-"Wet growing seasons only"
dw$climate2[dw$climate1=="Dry" & dw$group.site=="dry.or.wet"]<-"Dry growing seasons only"
# relevel climate1 extremes 
dw$climate2<-factor(dw$climate2, levels=c("Dry growing seasons only", "Wet growing seasons only", "Dry and Wet growing seasons"))
# check sample size 
table(dw$climate2)

## merge with the data for water balance and experimental years used 
site.inf<-gs%>%select(site_code, habitat, continent, country, latitude, longitude, first_nutrient_year)%>%distinct()%>%
  merge(y=wat.bal1[,c("site_code", "growing.season", "water.balance.last.15")], by="site_code")%>%
  merge(y=exp.yr[,c("site_code", "years.used")], by=c("site_code"))%>%
  merge(y=dw%>%select("site_code", "climate2")%>%distinct(), by=c("site_code"))

# seems that first nutrient year is not always the same to first fence year 
site.inf1<-site.inf%>%mutate(first_fenced_year=NULL)%>%dplyr::rename(first.experimental.year=first_nutrient_year)
site.inf1$water.balance.last.15<-round(site.inf1$water.balance.last.15, 2)
site.inf1$latitude<-round(site.inf1$latitude, 2)
site.inf1$longitude<-round(site.inf1$longitude, 2)
# write_xlsx(site.inf1, path="information for sites used.xlsx")

###########################################################################################
####################################  world map ###########################################
###########################################################################################

# Create world map for all sites and experimental years
WorldData <- map_data('world') %>%
  filter(region != "Antarctica") %>%
  fortify()

(pp.map<-ggplot(data = WorldData, aes(x = long, y = lat)) +
    geom_map(map = WorldData, aes(map_id = region), size = 0.5, fill = "gray99", color = "light grey") +
    geom_point(data = site.inf1, aes(x = longitude, y = latitude, fill = water.balance.last.15, shape = climate2), alpha=0.5, size=5, position=position_jitter(h=0.5, w=0.5)) +
    scale_fill_gradient2(low = "#E69F00", high = "#0072B2") +
    scale_shape_manual(values = c(21, 22, 25)) +
    labs(x = NULL, y = NULL, fill = "Water balance", shape = "Climate extremes") +
    scale_y_continuous(breaks=c(), expand = c(0, 0)) +
    scale_x_continuous(breaks=c(), expand = c(0, 0)) +
    theme_cowplot(font_size = 20)+panel_border()+
    #guides(fill = guide_colourbar(order = 1), shape = guide_legend(order = 2)) +
    theme(legend.position = c(0.02, 0.3), legend.background = element_blank()))
# ggsave(pp.map, file="map1.png", height = 6.64, width = 13.3, dpi = 600)

# the end 
