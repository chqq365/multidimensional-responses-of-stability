
###########################################################################################
################ calculate resistance and recovery for all non-extreme years ##############
###########################################################################################
table(num.year.at.least.4$number.of.years)
# select sites with more than 10 years of treatments
num.year.at.least.10<-num.year.at.least.4%>%filter(number.of.years >=10)

# check whether effects differ using 5-year data and 10-year data 
summary.eff.trt.facets<-c(); summary.estimated.relation.stability.facets<-c(); summary.estimated.relation.community.aspects<-c()
for(dura in c(5, 10)){ 
# dura<-10
yr.i<-num.year.at.least.10%>%filter(number.of.years>= dura)
# delete sites that do not have one extreme or normal year 
s_dd7_check.extreme<-s_dd7 %>%filter(site_code %in%yr.i$site_code)%>%filter(year_trt>0 & year_trt< dura +1)%>%
  mutate(extreme.events=ifelse(climate=="Normal", "normal", "extreme"))%>%
  group_by(site_code, extreme.events)%>%summarise(N=length(year_trt))%>%filter(N>0)%>%
  group_by(site_code)%>%summarise(N.unique=length(extreme.events))%>%filter(N.unique>1)

s_dd7_select<-s_dd7 %>%filter(site_code %in%yr.i$site_code)%>%filter(year_trt>0 & year_trt< dura +1) %>%filter(site_code %in%s_dd7_check.extreme$site_code)

data.bio.rich<-d8%>%filter(site_code %in%yr.i$site_code)%>%filter(year_trt< dura +1)%>%
 mutate(biomass=live_mass)%>% select(site_code, block, plot, trt, year_trt, year, biomass, richness)%>%distinct()%>%
  pivot_longer(cols = c("biomass", "richness"), names_to ="community.property", values_to = "property.value" )

# merge biomass and richness with climate extremes 
d.select<-s_dd7_select%>%select("site_code", "year", "year_trt", "spei", "climate")%>% 
  merge(data.bio.rich, by=c("site_code", "year", "year_trt"))%>%
  mutate(variable.id=paste(site_code, block, trt, community.property, sep="_"))%>%arrange(variable.id)

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
resis.recov<-resis1%>%bind_rows(recov1)%>%merge(d.select%>%select(site_code, block, trt, community.property, variable.id)%>%distinct(), by=c("variable.id"))

###################### resistance and recovery for composition ############################# 

d.select.com<-d8%>%filter(site_code %in%yr.i$site_code)%>%filter(year_trt< dura +1)%>%
  merge(s_dd7_select[,c("site_code", "year", "spei", "climate")], by=c("site_code", "year"))%>%
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

select.years.extremes<-s_dd7_select%>%select(site_code, year_trt, climate, spei)%>%filter(climate!="Normal")%>%
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
  merge(y=s_dd7_select[,c("site_code", "year_trt", "climate", "bio.data2", "spei")], by=c("site_code", "year_trt"), all.y=T)%>%
  mutate(year.used.1=ifelse(is.na(year.used), "No", year.used))%>%
  mutate(year.used.2=ifelse((year_trt==0), "No", year.used.1), year.used.3=ifelse((climate=="Normal"), "Yes", year.used.2))

# when two same climate extremes happen consecutively, recovery only calculate for the later year 
# the later year should be followed by either a less extreme year or normal year
select.resis.recov<-c()
for(cut in c(0.67, 1.28)){
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
    ggplot(aes(year_trt, spei1, color=climate1, shape=stability.facets))+theme_cowplot()+panel_border()+
    geom_point(size=5, alpha=0.7, position = position_dodge2(w = 0.5))+
    facet_wrap(~site_code,  scale="free_x")+
    #scale_x_continuous(expand=c(0,0), limits = c(0.5,15), breaks = seq(0,16,2))+
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
for (cut in  c(0.67, 1.28)){
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
  dplyr::rename(p='p-value')%>%dplyr::select(cutoff, community.property, stability.facets1, terms, Value, Std.Error,  DF,  "t-value", "p", r2m, r2c, number.sites, get.sd.block, get.sd.sites)%>%
  arrange(cutoff, community.property, stability.facets1)
table(eff.trt.facets.full$community.property)
eff.trt.facets.full[,5:13]<-round(eff.trt.facets.full[,5:13], 2)


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

estimated.ci1<-corr.stab.facet.l%>%select(cutoff, community.property, correlation.type, variable.id)%>%distinct()%>%
  merge(estimated.ci, by=c("variable.id"))%>%mutate(across(9:13,\(x) round(x, 2)))%>%mutate(variable.id=NULL, direction=ifelse(emmean>=0, "positive", "negative"))
# write_xlsx(estimated.ci1%>%mutate(number.sites=NULL, direction=NULL, life.form=NULL), path="estimated relationships among stability facets using emmeans.xlsx", col_names = TRUE)

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
estimated.ci.asp1<-corr.asp.l%>%select(cutoff, stability.facets1 , correlation.type, variable.id)%>%distinct()%>%
  merge(estimated.ci.asp, by=c("variable.id"))%>%mutate(across(6:13, \(x) round(x, 2)))%>%mutate(variable.id=NULL)
# write_xlsx(estimated.ci.asp1%>%mutate(number.sites=NULL), path="estimated relationships among community aspects using emmeans.xlsx", col_names = TRUE)

summary.eff.trt.facets<-summary.eff.trt.facets%>%bind_rows(eff.trt.facets.full%>%mutate(duration=dura))
summary.estimated.relation.stability.facets<-summary.estimated.relation.stability.facets%>%bind_rows(estimated.ci1%>%mutate(duration=dura))
summary.estimated.relation.community.aspects<-summary.estimated.relation.community.aspects%>%bind_rows(estimated.ci.asp1%>%mutate(duration=dura))
}

summary.eff.trt.facets.sig<-summary.eff.trt.facets%>%filter(p<=0.05  & terms!="(Intercept)" & stability.facets1!="invariability" & cutoff==0.67)

summary.estimated.relation.stability.facets.sig<-summary.estimated.relation.stability.facets%>%
  mutate(sig=ifelse(((upper.CL/lower.CL)>0|upper.CL==0|lower.CL==0), "significant", "non-significant"))%>%filter(sig=="significant" & cutoff==0.67)

summary.estimated.relation.community.aspects.sig<-summary.estimated.relation.community.aspects%>%
  mutate(sig=ifelse(((upper.CL/lower.CL)>0|upper.CL==0|lower.CL==0), "significant", "non-significant"))%>%filter(sig=="significant" & cutoff==0.67)

# write_xlsx(summary.eff.trt.facets, path="effects of nutrient addition on stability facets in the short and long term.xlsx", col_names = TRUE)
# write_xlsx(summary.estimated.relation.stability.facets, path="estimated relationships among stability facets using emmeans in the short and long term.xlsx", col_names = TRUE)
# write_xlsx(summary.estimated.relation.community.aspects, path="estimated relationships among community aspects using emmeans in the short and long term.xlsx", col_names = TRUE)

