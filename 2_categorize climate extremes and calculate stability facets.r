## note that we have data from year_trt 2 for look.us and lagoas.br !!!!
## correction from site pi for growing season (start to end) for the following sites
## kilp.fi  6-9(originally 6-8)
## saana.fi 6-9(originally 6-8)
# stability facets in richness and biomass (but not for composition) is on log scale 

rm(list=ls())
## open the library
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

#################################################################################################
##################################### growing season data #######################################
#################################################################################################
# growing season for each site using data from Sid 
# the harvest month refer to the end of growing season in Sid's data 
gs0<-read.csv(paste0(dir.data, "Weather_site_inventory_20200921.csv"), header=T, sep=",") 
# add growing season for site msla_2.us and msla_3.us, they are very close to site msla.us
gs <- gs0 %>%
  bind_rows(gs0 %>% filter(site_code == "msla.us") %>% mutate(site_code = "msla_2.us")) %>%
  bind_rows(gs0 %>% filter(site_code == "msla.us") %>% mutate(site_code = "msla_3.us")) %>%
  mutate(X = NULL) 
## check sites with and without growing season data 
gs1<- gs %>% filter(site_code %in% d8$site_code) %>% dplyr::select(site_code, harvest_month, gs_start_month, gs_len)
(gs.unknown <- unique(d8$site_code[ ! d8$site_code %in% gs$site_code])) %>% sort()

## check with site PI except for the site "azi.cn"    "barta.us" "lake.us" 
pi<-read.csv(paste0(dir.data, "pi-contact-list-8-March-2021.csv"), sep=",", head=T)
gs.unknown.pi <- pi %>% filter(pi$site_code %in% gs.unknown)
## also check site lagoas.br
#lagoas.br<-pi%>%filter(site_code =="lagoas.br")
## I communicated with site PIs from those sites where growing season was not known, some PIs responded. 
#  sites     start, peak, end,   notes
# gilb.za,   10,    4 ,    4,      peak may be in April, but they usually harvest in Feb or Mar    
# hnvr.us    4 ,    7,     8, 
# lagoas.br  10,    3,     4,
# potrok.ar  10,    4,     4,
# summ.za    9,     4,     5,      peak may be in April, but they usually harvest in Feb or Mar
# neba.jp    4,     8,     9,      generally start at end of April, although there is some fluctuation due to the amount of snowfall, and peaks at July-August (we collect samples in July)
# lake.us    5,            9,       # harvest in August 
# ahth.is    6,     8,     9,
# amlr.is    6,     8,     9,
# also included two sites where site pi corrected growing season 

add.site<-data.frame(t(data.frame(c("gilb.za", 10, 4), c("hnvr.us", 4, 8),
                                  c("potrok.ar", 10, 4),c("lagoas.br", 10, 4),
                                  c("summ.za", 9, 5), c("neba.jp", 4, 9), 
                                  c("lake.us", 6, 9),
                                  c("ahth.is", 6, 9), c("amlr.is", 6, 9),
                                  c("kilp.fi", 6, 9), c("saana.fi", 6, 9))))
names(add.site)<-c("site_code", "gs_start_month", "harvest_month")
rownames(add.site)<-NULL
add.site[,2:3]<-lapply(add.site[,2:3], as.numeric)

colnames(gs1)
gs3 <- gs1 %>% filter(!site_code %in% c("kilp.fi", "saana.fi")) %>% 
  dplyr::select(site_code, gs_start_month, harvest_month) %>% bind_rows(add.site)
## find the start and end of the experiment (year_trt 1 to ...)
## note that we have data from look.us and lagoas.br from year_trt 2!!!!
gs4 <- d8 %>% group_by(site_code) %>% summarise_at("year", list(max, min)) %>%
  mutate(min.y=ifelse((site_code == "look.us" | site_code == "lagoas.br"), fn2 - 1, fn2)) %>% mutate(fn2 = NULL) %>%
  dplyr::rename(site_lyear = fn1, site_fyear = min.y, fn1 = NULL, fn2 = NULL) %>%
  merge(y = gs3, by = c("site_code"))
length(unique(gs4$site_code))

## check sites that still do not have growing season data 
che.sites <- d8 %>% select(site_code, year_trt) %>% filter(! site_code %in% gs4$site_code) %>% distinct()
unique(che.sites$site_code)# sites do not have growing season data, these sites do not have long-term time series 

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
# Compute water balance (BAL) and standardize it 
s_dd3 <- s_dd2.within %>% filter(year!=1901) %>%  dplyr::rename(year1=year) %>% 
  bind_rows(s_dd2.across) %>%
  mutate(dif.re0 = prec - pet1) %>%
  group_by(site_code) %>%mutate(dif.re = scale(dif.re0))
range(s_dd3$year1)
length(unique(s_dd3$site_code))

## check the means 
ggplot(s_dd3, aes(year1, dif.re)) + theme_bw() + geom_point() +
  facet_wrap(~site_code, scale="free_y")+
  geom_hline(yintercept = 0, linetype="dashed", color="red")
che <- s_dd3 %>% group_by(site_code) %>% summarise_at("dif.re", list(mean, sd))## 

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
s_dd5$climate[s_dd5$dif.re>=1.28]<-"Extreme wet"
s_dd5$climate[s_dd5$dif.re<=-1.28]<-"Extreme dry"
s_dd5$climate[s_dd5$dif.re>=0.67 & s_dd5$dif.re<1.28]<-"Moderate wet"
s_dd5$climate[s_dd5$dif.re<=-0.67 & s_dd5$dif.re>-1.28]<-"Moderate dry"
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
s_dd6<-s_dd5%>%filter(!site_code%in%sss3)%>%mutate(year_trt=year-site_fyear+1)%>%dplyr::rename(spei=dif.re)
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
# ggsave(pp, width=13.3, height=6.64, dpi=600, file="moderate and extreme growing seasons.png")
ggsave(pp, file="climate extremes.pdf", width = 21, height = 29.7, dpi=600)

# write.csv(s_dd7, file="raw data of moderate and extreme growing seasons at 55 sites.csv")

###########################################################################################
###calculate resistance and recovery for all non-extreme years for different functional groups##############
###########################################################################################
# load biomass and richness for each functional group
load("richness and biomass for each functional group.rdata")
str(rich.bio.each.functional.group)

# add number of years 
num.year<-d8%>%select(site_code, year, year_trt)%>%distinct()%>%group_by(site_code)%>%
  summarise(number.of.years=length(year_trt))
# check whether biomass were sorted to different functional groups during experiments
che.stru<-rich.bio.each.functional.group%>%group_by(site_code, trt, block,  life.form, community.property)%>%
   summarise(N=length(property.value))%>%merge(num.year, by=c("site_code"))%>%
  mutate(dif=N-number.of.years)%>%filter(dif!=0)%>%filter(dif< -1)
unique(che.stru$site_code)
che.stru1<-rich.bio.each.functional.group%>%filter(site_code %in% c("kbs.us", "lancaster.uk", "shps.us"))%>%arrange(site_code, block, trt, community.property, life.form, year_trt)
# kbs.us only has functional group data in year 1, lancaster.uk only has functional group data in year 9 
# shps.us data in year 1 and 2 are missing in block 3. 

# merge biomass and richness in each functional groups with climate extremes (focus on total biomass and richness)
d.select<-s_dd7%>%filter(year_trt!=0)%>%select("site_code", "year", "year_trt", "spei", "climate")%>% 
  full_join(rich.bio.each.functional.group, by=c("site_code", "year", "year_trt"), multiple = "all")%>%
  filter(life.form=="all")%>%
  mutate(variable.id=paste(site_code, block, trt, life.form, community.property, sep="_"))%>%arrange(variable.id)

## calculate mean values in normal years 
d.select_norm<-d.select%>%filter(climate=="Normal")%>% group_by(variable.id) %>% dplyr::summarise(avg.property=mean(property.value))

d.select_r1<-d.select%>%merge(d.select_norm, by=c("variable.id"))%>% mutate(resistance=avg.property/abs(property.value - avg.property))
 
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
resis.recov<-resis1%>%bind_rows(recov1)%>%merge(d.select%>%select(site_code, block, trt, life.form, community.property, variable.id)%>%distinct(), by=c("variable.id"))

###################### resistance and recovery for composition ############################# 

d.select.com<-d8%>%filter(life.form!="OTHER")%>%mutate(life.form="all")%>%
  merge(s_dd7[,c("site_code", "year", "spei", "climate")], by=c("site_code", "year"))%>%
  mutate(community.property="composition", variable.id=paste(site_code, block, trt, life.form, community.property, sep="_"))%>%arrange(variable.id)

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
rr.com1<-d.select.com%>%dplyr::select("variable.id", "site_code", "trt", "block", "year_trt", "climate", "spei", "life.form")%>%distinct()%>%
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
   merge(d.select.com%>%select(site_code, block, trt, life.form, variable.id)%>%distinct(), by=c("variable.id"))%>%
  mutate(community.property="composition")

# add this to richness and biomass data
colnames(resis.recov)
colnames(resis.recov.com)
resis.recov1<-resis.recov%>%select(variable.id, year_trt, climate, spei, stability.facets, values, site_code,  block,  trt, life.form, community.property)%>%distinct()%>%
  bind_rows(resis.recov.com)

# test whether direction (dry and wet) and intensity impact resistance and recovery (log transform richness and biomass)
resis.recov.dir.int<-resis.recov1%>%filter(life.form=="all")%>%
  mutate(intensity=abs(spei), direction=ifelse(grepl("wet", climate), "wet", "dry"))%>%
  mutate(log.values=ifelse(community.property=="composition", values, log(values)))%>%filter(!is.na(log.values))%>%filter(!log.values %in% c("Inf", "-Inf"))

test.continuous<-c()
for(pro in unique(resis.recov.dir.int$community.property)){ 
 for (i in unique(resis.recov.dir.int$stability.facets)){
  # i<-"resistance"
  data.temp<-resis.recov.dir.int%>%filter(stability.facets==i & community.property==pro)
  n.sites<-length(unique(data.temp$site_code))
  
  mod<-lme(log.values ~ trt*direction*intensity, random=~1|site_code/year_trt, data=data.temp)
  plot(mod)
  t<-xtable(anova(mod))
  rs<-rsquared(mod)
  t1<-t%>%mutate(r2m=rs$Marginal, r2c=rs$Conditional, n.sites=n.sites, stability.facets=i, community.property=pro, terms=rownames(.))
  test.continuous<-test.continuous%>%bind_rows(t1)
  }
}
test.continuous[,3:6]<-lapply(test.continuous[,3:6], function(x) round(x, 2))
colnames(test.continuous)[4]<-"p"
test.continuous1<-test.continuous%>%select(stability.facets, community.property, terms, 'F-value', numDF, denDF, p, r2m, r2c, n.sites)
test.continuous.sig<-test.continuous%>%filter(p<=0.1)
# write_xlsx(test.continuous1, path="test interaction between treatment, direction, and intensity of climate extremes on resistance and recovery.xlsx")

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
# ggsave(pp1, dpi=600, width=13.3, height=6.64, file=paste0("SPEI based on cutoff ", cut, " indicating years used.png"))
ggsave(pp1, file="SPEI based on cutoff 0.67 and 1.28 indicating years used.pdf", width = 21, height = 29.7, dpi=600)

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
    select(variable.id, year_trt, climate1, spei1, stability.facets, values3, site_code,  block,  trt, life.form, community.property)%>%distinct()%>%
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
# ggsave(check.year.0.67, dpi=600, width=13.3, height=6.64, file="years used for resistance and recovery based on cutoff 0.67.png")
ggsave(check.year.0.67, file="years used for resistance and recovery based on cutoff 0.67.pdf", width = 21, height = 29.7, dpi=600)

# average resistance during dry and wet and recovery from dry and wet for each site 
colnames(select.resis.recov)
s_resis.recov<-select.resis.recov%>%mutate(stability.facets1=paste(stability.facets, climate1, sep="_"))%>%
    filter(!values3%in%c("Inf", "-Inf"))%>%filter(!is.na(values3))%>%
    group_by(cutoff, site_code, trt, block, life.form, community.property, stability.facets1, variable.id)%>%
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
  merge(d.select%>%select(site_code, block, trt, life.form, community.property, variable.id)%>%distinct(), by=c("variable.id"))%>%
   dplyr::select(site_code, trt, block, life.form, community.property, invariability, invariability.d)%>%
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
inv.com1<-d.select.com%>%dplyr::select("variable.id", "site_code", "trt", "block", "life.form", "community.property")%>%distinct()%>%
  merge(y=inv.com, by=c("variable.id"))%>%mutate(values=sim, stability.facets1="invariability")%>%
  dplyr::select(site_code, trt, block, life.form, community.property, stability.facets1, values)

# add all temporal invariability from three community aspects 
stb3<-stb2%>%bind_rows(inv.com1)
table(stb3$community.property)
check<-stb2%>%filter(grepl("richness", community.property))%>%filter(stability.facets1=="invariability")%>%
    merge(inv.com1[,c("site_code", "trt", "block", "life.form", "values")], by=c("site_code", "trt", "block", "life.form"), all.x=T)%>%
  filter(is.na(values.y))%>%   filter(!is.na(values.x))
unique(check$life.form)
# when species occur only in one year, but not in other years, it is not possible to calculate invariability in composition

############################################################################################
############################ add all stability facets together  ############################
############################################################################################
# using the sames sites for resistance and recovery (which differ using different cutoffs) for invariability
colnames(s_resis.recov)
colnames(stb3)

data.stability.facets<-c()
for (cut in  c(0.67, 1, 1.28, 1.5, 2)){
  # cut<-0.67
  temp.data.cut<- s_resis.recov%>%select(site_code, trt, block,life.form, community.property, stability.facets1,  values, cutoff)%>%
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
# delete cutoff of 2 sd, because model does not converge for invariability graminoids in composition 
# for cutoff of 1.5 and 2, maybe just focus on the whole community 
all.data.l_1<-data.stability.facets.sub%>%mutate(variable.id=paste(cutoff, community.property, stability.facets1, life.form, sep="_"))%>% 
  mutate(comb.cut.form=paste(cutoff, life.form))%>%filter(!comb.cut.form %in% c("1.5 FORB", "1.5 GRAMINOID", "1.5 WOODY", "1.5 LEGUME", "2 FORB", "2 GRAMINOID", "2 LEGUME", "2 WOODY"))%>%
   mutate(values1=ifelse(community.property=="composition", values, log(values)))
# look at the data distribution 
# all.data.l_1%>%filter(cutoff==1.28)%>%ggplot()+geom_density(aes(x=values1, color=trt))+facet_wrap(~variable.id, scales = "free")

eff.trt.facets<-c()
for(i in unique(all.data.l_1$variable.id)){
  # i<-"1.28_richness_recovery_Dry_all"   
  data1<-all.data.l_1%>%filter(variable.id==i)%>%filter(!values1%in%c("Inf", "-Inf"))
  # record how many sites have these two variables 
  n.sites<-length(unique(data1$site_code))
  check.values<-unique(data1$values1)
  if(length(check.values)<3|n.sites<3) next
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
eff.trt.facets.full<-all.data.l_1%>%ungroup()%>%dplyr::select(cutoff, community.property, stability.facets1, life.form, variable.id)%>%unique()%>%
  merge(eff.trt.facets, by=c("variable.id"))%>%
  dplyr::rename(p='p-value')%>%dplyr::select(cutoff, community.property, stability.facets1, life.form, terms, Value, Std.Error,  DF,  "t-value", "p", r2m, r2c, get.sd.block, get.sd.sites)%>%
  arrange(cutoff, community.property, stability.facets1, life.form)
table(eff.trt.facets.full$community.property)
eff.trt.facets.full[,6:14]<-round(eff.trt.facets.full[,6:14], 2)

eff.trt.facets.sig<-eff.trt.facets.full%>%filter(p<=0.05  & terms!="(Intercept)")
 write_xlsx(eff.trt.facets.sig, path="significant treatment effects on stability facets.xlsx", col_names = TRUE)
 write_xlsx(eff.trt.facets.full, path="full table of treatment effects on stability facets.xlsx", col_names = TRUE)

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
(pp.trt<-eff.trt4%>%filter(cutoff %in% c(1.28))%>%ggplot(aes(stability.facets1, Value, color=community.property))+theme_bw(base_size = 20)+
    facet_wrap(~life.form , nrow=4, scales = "free_y")+
    geom_point(position=pd, pch=19, size=5, alpha=0.6)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    geom_hline(yintercept = 0, linetype ="dotted") +
    guides(color=guide_legend(nrow=1))+
    theme(legend.position = c(0.8, 0.2), axis.text.x = element_text(angle = 15))+
    labs(x=NULL, y="Effect size of nutrient addition", color=NULL))
# ggsave(pp.trt, width=13.3, height=6.64, dpi=600, file="effects on stability facets in biomass, richness and composition.png")

# to get the legend for pentagons (the effect size is wrong in this figure)
colour.crosswalk <- c("positive" = "black", "negative" = "red")
(pp.trt.all<-eff.trt4%>%filter((cutoff %in% c(1.28) & life.form=="all"))%>%mutate(Value.scaled=abs(Value)*10)%>%
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
  mutate(id=paste(cutoff, life.form, community.property, sep="_"), id1=paste(cutoff, life.form, sep="_"))

list_plots <- vector('list', length(unique(facet_data$id)))
for (i in unique(facet_data$id)){
    # i<-"1.28_all_richness"
    facet_data_temp<-facet_data%>%filter(id==i)%>% mutate(sig=case_when(p<=0.05 ~ "significant",  TRUE~"non-significant"))
    facet_data_temp$sig<-factor(facet_data_temp$sig, levels = c("non-significant", "significant"))
   if(grepl("all", i)){fontsize <- 16 / .pt; size.npk<-20 / .pt}else{fontsize <- 11 / .pt; size.npk<-16 / .pt}
    list_plots[[i]]<-gg.pentagon(data = facet_data_temp)+ 
      geom_label(aes( x = 0, y = 0, label = "NPK"), fill="white", size = size.npk)
  }

for(i in c("0.67_all" , "1.28_all")){ 
  # i<-"1.28_all"
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  # rename each element 
  names(list_plots_temp)<-gsub( "(.*)_(.*)_(.*)", "\\3", names(list_plots_temp))
(pp.trt.effects<-plot_grid(list_plots_temp$biomass,
                           list_plots_temp$composition,
                           list_plots_temp$richness,
                           legend,
                           rel_widths = c(4.5,4.5,4.5,1.5),
          nrow=1, vjust = 5, hjust=-0.1, 
           labels = c("A (Biomass)", "B (Composition)", "C (Richness)"), label_fontface = "bold", label_size = 18))

 ggsave(pp.trt.effects, height=6, width=13.3, dpi=600, file=paste0("effects on stability facets for ", i , ".png"))
}

for(i in c(0.67, 1.28)){ 
  # i<-"1.28"
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  list_plots_temp1<-list_plots_temp[-grep("all",  names(list_plots_temp))]
  # rename each element 
  names(list_plots_temp1)<-gsub( "(.*)_(.*)_(.*)", "\\2\\3", names(list_plots_temp1))
  (pp.trt.effects<-plot_grid(list_plots_temp1$GRAMINOIDbiomass,
                             list_plots_temp1$FORBbiomass,
                             list_plots_temp1$LEGUMEbiomass,
                             list_plots_temp1$WOODYbiomass,
                             list_plots_temp1$GRAMINOIDcomposition,
                             list_plots_temp1$FORBcomposition,
                             list_plots_temp1$LEGUMEcomposition,
                             list_plots_temp1$WOODYcomposition,
                             list_plots_temp1$GRAMINOIDrichness,
                             list_plots_temp1$FORBrichness,
                             list_plots_temp1$LEGUMErichness,
                             list_plots_temp1$WOODYrichness,
                             nrow=3, ncol=4))
  (pp.trt.effects1<-plot_grid(pp.trt.effects,
                              legend, rel_widths = c(8.8, 1.2)))
                              
  (pp.trt.effects2<-ggpubr::annotate_figure(pp.trt.effects1,
                                            top=text_grob("Graminoids                                   Forbs                                   Legumes                                   Woody                            ", size=15, face = "bold") ,
                                            left=text_grob("Richness                      Composition                      Biomass", rot=90, face = "bold", size=15)))
 
     ggsave(pp.trt.effects2, height=6.64, width=13.3, dpi=600,  file=paste0("effects on stability facets for functional groups ", i , ".png"))
}

#############################################################################################
###calculate correlation between stability facets within biomass, composition, and richness within sites ##
#############################################################################################
# focus on the whole community and cutoff of 0.67 and 1.28 sd
data.stability.facets.relation<-data.stability.facets.sub%>%filter(life.form=="all" & cutoff %in% c(0.67, 1.28))%>%
  filter(!(stability.facets1=="invariability" & community.property=="biomass"))%>%
  filter(!(stability.facets1=="invariability" & community.property=="richness"))%>%
  mutate(stability.facets2=gsub(".d", "", stability.facets1))%>%
  mutate(stability.facets1= stability.facets2, stability.facets2=NULL)

all.data.l_1<-data.stability.facets.relation%>%
  mutate(variable.id=paste(cutoff, community.property, life.form, site_code,  trt, sep="_"))%>%
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
corr.stab.facet.l<-all.data.l_1%>%ungroup()%>%select(cutoff, community.property, life.form, site_code,  trt, variable.id)%>%distinct()%>%
  merge(corr.stab.facet, by=c("variable.id"))%>%mutate(variable.id=NULL)%>% 
  pivot_longer(cols =inv_resis.d:recov.d_recov.w, names_to = "correlation.type")%>%filter(!is.na(value))%>%
  mutate(variable.id=paste(cutoff, community.property, correlation.type, life.form,  sep="_"))
hist(corr.stab.facet.l$value)
# check the data (correlation coefficients should not be |1| which is perfect fit or 0 which show no trends at all)
che.range<-corr.stab.facet.l%>%filter(value %in% c(0, -1, 1))

# look at the data distribution 
colnames(corr.stab.facet.l)
corr.stab.facet.l%>%ggplot()+geom_density(aes(x=value, color=trt))+facet_wrap(~interaction(correlation.type, community.property), scales = "free")

out.corr<-c(); estimated.ci<-c()
for(id in unique(corr.stab.facet.l$variable.id)){ 
  # id<-"1.28_richness_recov.d_recov.w_all"
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

out.corr1<-corr.stab.facet.l%>%select(cutoff, community.property, correlation.type, life.form, variable.id)%>%distinct()%>%
  merge(out.corr, by=c("variable.id"))%>%mutate(variable.id=NULL)%>%mutate(across(5:12, \(x) round(x, 2)))%>%dplyr::rename(p='p-value')%>%
  arrange(cutoff,community.property,  life.form, terms, correlation.type)
table(out.corr1$community.property)
## sort the significant effects 
range(out.corr1$number.sites)
lme.relationship.facet.sig<-out.corr1%>%filter(p<=0.05)

# save the tables
 write_xlsx(lme.relationship.facet.sig, path="significant treatment effects on relationships among stability facets.xlsx", col_names = TRUE)
 write_xlsx(out.corr1, path="full table of treatment effects on relationships among stability facets.xlsx", col_names = TRUE)

################ plot treatment effects on correlation among stability facets  ##############
colnames(estimated.ci)
unique(estimated.ci$correlation.type)
estimated.ci1<-corr.stab.facet.l%>%select(cutoff, community.property, correlation.type, life.form, variable.id)%>%distinct()%>%
  merge(estimated.ci, by=c("variable.id"))%>%mutate(across(7:14, round, 2))%>%mutate(variable.id=NULL, direction=ifelse(emmean>=0, "positive", "negative"))
# also save a copy
estimated.ci2<-estimated.ci1%>%mutate(number.sites=NULL, direction=NULL, life.form=NULL)
write_xlsx(estimated.ci2, path="estimated relationships among stability facets using emmeans.xlsx", col_names = TRUE)

(pp.cor<-estimated.ci1%>%filter(cutoff==1.28 &life.form=="all")%>%mutate(Value.scaled=abs(emmean)*10)%>%
    ggplot(aes(correlation.type, Value.scaled,  colour = factor(direction), size = (Value.scaled/10)))+theme_bw(base_size = 16)+
    facet_grid(trt~community.property, scales = "free_y")+
    geom_point(position=pd, pch=15, alpha=0.6)+
    #geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=pd) +
    geom_hline(yintercept = 0, linetype ="dotted") +
    scale_colour_manual(values = colour.crosswalk) +
    theme(axis.text.x = element_text(angle=20, hjust=1))+
    #theme(legend.position = c(0.8, 0.2))+
    labs(x=NULL, y="Correlation coeffecient", size="Correlation\n coefficient", color="Direction"))
# ggsave(pp.cor, width=13.3, height=6.64, dpi=600, file="effects on correlation among stability facets in biomass, richness and composition.png")
legend.cor<-get_legend(pp.cor)

# plotting using pentagons
facet_data_cor_stab <- estimated.ci1%>%filter(cutoff %in% c(0.67, 1.28) &life.form=="all")%>%
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
  mutate(id=paste(cutoff,  life.form, community.property, trt, sep="_"), id1=paste(cutoff, life.form, sep="_"))

list_plots <- vector('list', length(unique(facet_data_cor_stab$id)))

for (i in unique(facet_data_cor_stab$id)){
  # i<-"1.28_biomass"
  facet_data_cor_stab_temp<-facet_data_cor_stab%>%filter(id==i)
  list_plots[[i]]<-gg.pentagon(data = facet_data_cor_stab_temp)
}

for(i in unique(facet_data_cor_stab$id1)){ 
  # i<-"1.28_all"
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  # rename each element 
  names(list_plots_temp)<-gsub( "(.*)_(.*)_(.*)_(.*)", "\\3\\4", names(list_plots_temp))
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
  mutate(variable.id=paste(cutoff, stability.facets1, life.form, site_code,  trt, sep="_"))%>%
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
corr.asp.l<-all.data.l_1%>%ungroup()%>%select(cutoff, stability.facets1, life.form, site_code,  trt, variable.id)%>%distinct()%>%
  merge(corr.asp, by=c("variable.id"))%>%mutate(variable.id=NULL)%>% 
  pivot_longer(cols =c("bio_com", "bio_div", "com_div"), names_to = "correlation.type")%>%filter(!is.na(value))%>%
  mutate(variable.id=paste(cutoff, stability.facets1, correlation.type, life.form,  sep="_"))

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
trt.corr.asp1<-corr.asp.l%>%select(cutoff, stability.facets1, correlation.type, life.form, variable.id)%>%distinct()%>%
  merge(trt.corr.asp, by=c("variable.id"))%>%mutate(variable.id=NULL)%>%mutate(across(5:12, round, 2))%>%dplyr::rename(p='p-value')%>%
  arrange(cutoff, stability.facets1,  life.form, terms, correlation.type)
table(trt.corr.asp1$stability.facets1)

# sort the significant effects 
range(trt.corr.asp1$number.sites)
lme.relationship.aspect.sig<-trt.corr.asp1%>%filter(p<=0.05)%>% filter(life.form=="all")
# also save the full tables
 write_xlsx(lme.relationship.aspect.sig, path="significant treatment effects on relationships among community aspects.xlsx", col_names = TRUE)
 write_xlsx(trt.corr.asp1, path="full table of treatment effects on relationships among community aspects.xlsx", col_names = TRUE)

################ plot treatment effects on correlation among community aspects ##############
colnames(estimated.ci.asp)
unique(estimated.ci.asp$correlation.type)
estimated.ci.asp1<-corr.asp.l%>%select(cutoff, stability.facets1 , correlation.type, life.form, variable.id)%>%distinct()%>%
  merge(estimated.ci.asp, by=c("variable.id"))%>%mutate(across(7:14, round, 2))%>%mutate(variable.id=NULL)
# also save a copy
estimated.ci.asp2<-estimated.ci.asp1%>%mutate(number.sites=NULL, life.form=NULL)
write_xlsx(estimated.ci.asp2, path="estimated relationships among community aspects using emmeans.xlsx", col_names = TRUE)

(pp.cor.asp<-estimated.ci.asp1%>%filter(cutoff==1.28 &life.form=="all")%>%ggplot(aes(correlation.type, emmean, color=trt))+theme_bw(base_size = 20)+
    facet_wrap(~stability.facets1, nrow=4, scales = "free_y")+
    geom_point(position=pd, size=5, pch=19, alpha=0.6)+
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=pd) +
    guides(color=guide_legend(nrow=1))+
    geom_hline(yintercept = 0, linetype ="dotted") +
     theme(axis.text.x = element_text(angle=20, hjust=1))+
    #theme(legend.position = c(0.8, 0.2))+
    labs(x=NULL, y="Correlation coeffecient", color=NULL))
# ggsave(pp.cor.asp, width=13.3, height=6.64, dpi=600, file="effects on correlation among community aspects for each stability facet.png")

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
  mutate(id=paste(cutoff,  life.form, stability.facets1, trt, sep="_"), id1=paste(cutoff, life.form, sep="_"))

list_plots <- vector('list', length(unique(Aspect_data_cor$id)))
for (i in unique(Aspect_data_cor$id)){
  # i<-"0.67_all_recovery_Dry_Control"
  Aspect_data_cor_temp<-Aspect_data_cor%>%filter(id==i)
  list_plots[[i]]<-gg.triangle(data = Aspect_data_cor_temp)
}

for(i in unique(Aspect_data_cor$id1)){ 
  # i<-"0.67_all"
  list_plots_temp<-list_plots[grep(i, names(list_plots))]
  # rename each element 
  names(list_plots_temp)<-gsub( "(.*)_(.*)_(.*)_(.*)", "\\2\\3\\4", names(list_plots_temp))
  (pp.asp<-plot_grid(list_plots_temp$allinvariabilityControl,
                             list_plots_temp$resistanceDryControl,
                             list_plots_temp$resistanceWetControl,
                             list_plots_temp$recoveryDryControl,
                             list_plots_temp$recoveryWetControl,
                             
                             list_plots_temp$allinvariabilityNPK,
                             list_plots_temp$resistanceDryNPK,
                             list_plots_temp$resistanceWetNPK,
                             list_plots_temp$recoveryDryNPK,
                             list_plots_temp$recoveryWetNPK,
                             NA, NA, NA,NA,NA,
                     rel_heights = c(4,4,2),
                             nrow=3))
  (pp.asp1<-plot_grid(pp.asp, legend.cor, rel_widths = c(8.6, 1.4)))
  (pp.asp2<-ggpubr::annotate_figure(pp.asp1, 
                                    top=text_grob("  Invariability            Resistance_Dry            Resistance_Wet            Recovery_Dry            Recovery_Wet                           ", size=16, face = "bold") ,
                                            left=text_grob("                     NPK                              Control", rot=90, face = "bold", size=16)))
  
  ggsave(pp.asp2, height=6.64, width=13.3, dpi=600, file=paste0("relationships between community Aspects ", i , ".png"))
}


#################################################################################################
################################# add results to a word document  ################################
#################################################################################################
library(officer)
Table.S1<-readxl::read_xlsx(path="information for sites used.xlsx", 1)
Table.S2<-readxl::read_xlsx(path="combination of extreme events and years used for resistance and recovery.xlsx")
Table.S2<-Table.S2%>%filter(!One.year.after%in%c("Moderate dry", "Moderate wet"))%>%
  filter(!One.year.before%in%c("Moderate dry", "Moderate wet"))%>%
  mutate(One.year.before1=gsub("Extreme d", "D", One.year.before), One.year.before2=gsub("Extreme w", "W", One.year.before1))%>%
  mutate(One.year.after1=gsub("Extreme d", "D", One.year.after), One.year.after2=gsub("Extreme w", "W", One.year.after1))%>%
  mutate(an.extreme.season1=gsub("Extreme d", "D", an.extreme.season), an.extreme.season2=gsub("Extreme w", "W", an.extreme.season1))%>%
  select(One.year.before2, an.extreme.season2, One.year.after2, year.resis, year.recov)

Table.S3<-readxl::read_xlsx(path="full table of treatment effects on stability facets.xlsx", 1)
Table.S3<-Table.S3%>%filter(cutoff==0.67 & life.form=="all" & stability.facets1!="invariability")%>%
  mutate(cutoff=NULL, life.form=NULL, n.sites=NULL, stability.facets1=ifelse(stability.facets1=="invariability.d", "invariability", stability.facets1))

Table.S4<-readxl::read_xlsx(path="estimated relationships among stability facets using emmeans.xlsx", 1)
Table.S4<-Table.S4%>%filter(cutoff==0.67)%>% mutate(cutoff=NULL)%>%
  arrange(community.property, correlation.type)

Table.S5<-readxl::read_xlsx(path="estimated relationships among community aspects using emmeans.xlsx", 1)
Table.S5<-Table.S5%>%filter(cutoff==0.67)%>%mutate(cutoff=NULL)%>%
  arrange(stability.facets1, correlation.type)

Table.S6<-readxl::read_xlsx(path="site PIs for sites used.xlsx")%>%mutate(email=NULL)

# adjust the names for the tables 
colnames(Table.S1)<-c("site_code",  "habitat", "continent", "latitude", "longitude", "first experimental year", "growing season", "water balance during last 15 years",   "years used")
colnames(Table.S2)<-c("One year before",   "an dry or wet growing season",  "One year after",    "year for resistance",         "year for recovery" )
colnames(Table.S3)<-c("community property",  "stability facet",  "terms",  "Value",  "Std Error",  "DF",
                      "t-value", "p", "r2m",   "r2c", "sd (block)",   "SD (site)" )
colnames(Table.S4)<-c("community property",  "correlation type", "trt",  "emmean", "SE",  "df",  
                      "lower.CL",  "upper.CL",  "r2m", "r2c",  "SD (site)" )
colnames(Table.S5)<-c("stability facet",  "correlation type", "trt",  "emmean", "SE",  "df",  
                      "lower.CL",  "upper.CL",  "r2m", "r2c",  "SD (site)" )

colnames(Table.S6)

fig.2<-list.files(pattern = c("effects on stability facets for 0.67_all"))
fig.3<-list.files(pattern = c("relationships between stability facets 0.67_all"))
fig.4<-list.files(pattern = c("relationships between community Aspects 0.67_all"))

fig.7<-list.files(pattern = c("moderate and extreme growing seasons"))
fig.8<-list.files(pattern = c("years used for resistance and recovery based on cutoff 0.67"))
fig.9<-list.files(pattern = c("above ground biomass at each site using cutoff of 1.28"))
fig.10<-list.files(pattern = c("richness at each site using cutoff of 1.28"))
fig.11<-list.files(pattern = c("above ground biomass at site look.us"))

fig.12<-list.files(pattern = c("aboveground biomass During and one year after"))
fig.13<-list.files(pattern = c("Species richness During and one year after"))
fig.14<-list.files(pattern = c("composition similarity During and one year after"))

fig.15<-list.files(pattern = c("effects on stability facets for functional groups 0.67"))
fig.16<-list.files(pattern = c("effects on stability facets for hill numbers 0.67"))

fig.17<-list.files(pattern = c("effects on stability facets for 1.28_all"))
fig.18<-list.files(pattern = c("relationships between stability facets 1.28_all"))
fig.19<-list.files(pattern = c("relationships between community Aspects 1.28_all"))

all.figs<-c(fig.2, fig.3, fig.4,  fig.7, fig.8, fig.9, fig.10, fig.11, fig.12, fig.13, fig.14, fig.15, fig.16, fig.17,fig.18, fig.19)

my.doc<-read_docx()
my.doc<-body_remove(my.doc) %>% cursor_end()
(word_size <- docx_dim(my.doc))

my.doc%>%body_add_par("Relationships between functional traits and resistance and recovery", style="heading 1")%>%
  body_add_par("  ")

for(i in 1:length(all.figs)){
  # i<-2
  my.doc%>%body_add_img(src=all.figs[i], width = 6.2, height=3.1, style="centered", pos="after")
  my.doc%>%body_add_par(gsub(".png", "", basename(all.figs[i])), style='graphic title', pos="after")%>%
    body_add_par("See Table S for test statistics.", style='Normal', pos="after")%>%body_add_par("  ")
  
}

my.doc%>% body_add_par("Table S1. geolocation of sites used",
                       style = "table title", pos="after")%>%body_add_par("  ")
my.doc%>%body_add_table(Table.S1, style='Table Professional', pos="after")%>%body_add_par("  ")


my.doc%>% body_add_par("Table S2. Years used for calculating resistance and recovery",
                       style = "table title", pos="after")
my.doc%>%body_add_table(Table.S2, style='Table Professional', pos="after")%>%body_add_par("  ")


my.doc%>% body_add_par("Table S3. An anova table summarizing fixed effects of nutrient addition on stability facets in community aspects",
                       style = "table title", pos="after")%>%
  body_add_par("Models were specified as lme (y~trt, random=~1|sites/blocks). 
                      r2m (marginal R2): proportion of variance explained by the fixed effects in the model;
                       r2c (conditional R2): proportion of variance explained by the fixed and random effects" ,
               style='Normal', pos="after")%>%body_add_par("  ")
my.doc%>%body_add_table(Table.S3, style='Table Professional', pos="after")%>%body_add_par("  ")


my.doc%>% body_add_par("Table S4. An anova table summarizing fixed effects of nutrient addition on relationships between stability facets in each community aspect",
                       style = "table title", pos="after")%>%
  body_add_par("Models were specified as lme (y~trt, random=~1|sites)
                      . r2m (marginal R2): proportion of variance explained by the fixed effects in the model;
                       r2c (conditional R2): proportion of variance explained by the fixed and random effects; 
                       SD (sites): standard deviation for random effects of sites.", style='Normal', pos="after")%>%body_add_par("  ")
my.doc%>%body_add_table(Table.S4, style='Table Professional', pos="after")%>%body_add_par("  ")


my.doc%>% body_add_par("Table S5.  An anova table summarizing fixed effects of nutrient addition on relationships between community aspects in each stability facet",
                       style = "table title", pos="after")%>%
  body_add_par("Models were specified as lme (y~trt, random=~1|sites). 
                      r2m (marginal R2): proportion of variance explained by the fixed effects in the model;
                       r2c (conditional R2): proportion of variance explained by the fixed and random effects; 
                       SD (sites): standard deviation for random effects of sites.", style='Normal', pos="after")%>%body_add_par("  ")

my.doc%>%body_add_table(Table.S5, style='Table Professional', pos="after")%>%body_add_par("  ")


my.doc%>% body_add_par("Table S6. Principal investigators contributing data but who are not authors",
                       style = "table title", pos="after")%>%body_add_par("  ")

my.doc%>%body_add_table(Table.S6, style='Table Professional', pos="after")%>%body_add_par("  ")


print(my.doc, target="stability facets and their correlations.docx")


# also save all tables into excel 
names_of_dfs<-c("Table.S1", "Table.S2", "Table.S3", "Table.S4", "Table.S5", "Table.S6")
# writexl:: write_xlsx(setNames(lapply(names_of_dfs,get),names_of_dfs), path = "all tables.xlsx")

# the end 
