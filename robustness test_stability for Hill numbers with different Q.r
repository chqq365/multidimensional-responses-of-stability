
rm(list=ls())
## open the libraries
library(tidyverse);library(ggplot2);library(cowplot);library(ggpubr)
library(chemodiv)
library(nlme); library(writexl); library(xtable); library(emmeans);library(piecewiseSEM)
## environment needed
pd<-position_dodge(width=0.5)

## set up the work directory 
dir.data<-"H:/resistance and recovery/R codes/raw data/"
dir.graphs<-"H:/resistance and recovery/R codes/graphs/"
setwd(dir.graphs)
## add data for biomass and cover
d8<-read.csv("biomass and cover data for 55 nutnet sites.csv")

d8.w<-d8%>%mutate(X=NULL, Family=NULL)%>%filter(max_cover>0)%>%
  pivot_wider(names_from = "standard_taxon", values_from = "max_cover")%>%distinct()
cover<-d8.w[,9:ncol(d8.w)]
cover[is.na(cover)]<-0

diversity.q<-c()
for (i in 0:2){
  # i<-0
  div_temp <-calcDiv(cover, type="HillDiv", q=i)%>%bind_cols(d8.w%>%select(1:8)%>%distinct())%>%
    mutate(q=i)
  diversity.q<-rbind(diversity.q, div_temp)
}

###########################################################################################
###calculate resistance and recovery for all non-extreme years for different diversity ####
###########################################################################################
climate.extremes<-read.csv("raw data of moderate and extreme growing seasons at 55 sites.csv")

colnames(diversity.q)[1]<-"property.value"
str(diversity.q)
d.select<-climate.extremes%>%filter(year_trt!=0)%>%select("site_code", "year", "year_trt", "spei", "climate")%>% 
  full_join(diversity.q, by=c("site_code", "year", "year_trt"), multiple = "all")%>%
  mutate(variable.id=paste(site_code, block, trt, q, sep="_"))%>%arrange(variable.id)

## calculate mean values in normal years 
d.select_norm<-d.select%>%filter(climate=="Normal")%>% group_by(variable.id) %>% dplyr::summarise(avg.property=mean(property.value))

d.select_r1<-d.select%>%merge(d.select_norm, by=c("variable.id"))%>% mutate(resistance=avg.property/abs(property.value - avg.property))

# get resistance data by deleting normal climate events 
resis1<-d.select_r1%>%filter(climate!="Normal")%>%mutate(stability.facets="resistance")%>%dplyr::rename(values=resistance)%>%
  select(variable.id, year, year_trt, climate, spei, stability.facets, values)

# make sure that recovery is always calculated as first year compared with the second year
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
resis.recov1<-resis1%>%bind_rows(recov1)%>%merge(d.select%>%select(site_code, block, trt, q, variable.id)%>%distinct(), by=c("variable.id"))

colnames(resis.recov1)
# for resistance and recovery, if the previous year is not the same extreme events, this year should not be included due to confounding effects
# for recovery, if the next year is a different extreme event, this year should not be included due to confounding effects
# the pre-treatment year and years with missing biomass should be included for selection 
select.years.extremes<-climate.extremes%>%select(site_code, year_trt, climate, spei)%>%filter(climate!="Normal")%>%
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
    select(variable.id, year_trt, climate1, spei1, stability.facets, values3, site_code,  block,  trt, q)%>%distinct()%>%
    filter(!is.na(values3))%>%mutate(cutoff=cut)
  
  select.resis.recov<-rbind(select.resis.recov, resis.recov2)
}

# average resistance during dry and wet and recovery from dry and wet for each site 
colnames(select.resis.recov)
s_resis.recov<-select.resis.recov%>%mutate(stability.facets1=paste(stability.facets, climate1, sep="_"))%>%
  filter(!values3%in%c("Inf", "-Inf"))%>%filter(!is.na(values3))%>%
  group_by(cutoff, site_code, trt, block, q, stability.facets1, variable.id)%>%
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
  merge(d.select%>%select(site_code, block, trt, q, variable.id)%>%distinct(), by=c("variable.id"))%>%
  dplyr::select(site_code, trt, block, q, invariability, invariability.d)%>%
  pivot_longer(cols = c("invariability", "invariability.d"), names_to ="stability.facets1" , values_to ="values")

############################################################################################
############################ add all stability facets together  ############################
############################################################################################
# using the sames sites for resistance and recovery (which differ using different cutoffs) for invariability
colnames(s_resis.recov)
colnames(stb2)

data.stability.facets<-c()
for (cut in  c(0.67, 1.28)){
  # cut<-0.67
  temp.data.cut<- s_resis.recov%>%select(site_code, trt, block,q, stability.facets1,  values, cutoff)%>%
    filter(cutoff==cut)
  n.sites<-length(unique(subset(temp.data.cut, cutoff==cut)$site_code))## check how many sites were included for analyses
  all.stability<-stb2%>%filter(site_code%in%temp.data.cut$site_code)%>%mutate(cutoff=cut)%>%bind_rows(temp.data.cut)%>%mutate(n.sites=n.sites)         
  data.stability.facets<-rbind(data.stability.facets, all.stability)
}
unique(data.stability.facets$n.sites)

#############################################################################################
############ nutrient addition effects on stability for all aspects and  facets  ############
#############################################################################################
# run analyses 
colnames(data.stability.facets)
data.stability.facets.sub<-data.stability.facets%>%filter(!is.na(values))%>%filter(!values%in% c("Inf", "-Inf"))

all.data.l_1<-data.stability.facets.sub%>%mutate(variable.id=paste(cutoff, stability.facets1, q, sep="_"))%>% 
  mutate(values1=log(values))
# look at the data distribution 
all.data.l_1%>%filter(cutoff==1.28)%>%ggplot()+geom_density(aes(x=values1, color=trt))+facet_wrap(~variable.id, scales = "free")

eff.trt.facets<-c()
for(i in unique(all.data.l_1$variable.id)){
  # i<-"1.28_com_recovery_Dry_all" 
  data1<-all.data.l_1%>%filter(variable.id==i)%>%filter(!values1%in%c("Inf", "-Inf"))
  # record how many sites have these two variables 
  n.sites<-length(unique(data1$site_code))
  check.values<-unique(data1$values1)
  if(length(check.values)<3|n.sites<3) next
  # ggplot(data1)+geom_histogram(aes(x=values1))+facet_wrap(~trt)
  mod<-lme(values1 ~ trt, random=~1|site_code/block,  data=data1)
  # plot(mod)
  t1<-xtable(summary(mod)$tTable)
  rs<-rsquared(mod)
  t1$r2m<-rs$Marginal
  t1$r2c<-rs$Conditional
  t1$terms<-rownames(t1)
  t1$variable.id<-i
  t1$n.sites<-n.sites
  eff.trt.facets<-rbind(eff.trt.facets, t1) 
}

# save the full table 
eff.trt.facets.full<-all.data.l_1%>%ungroup()%>%dplyr::select(cutoff,  stability.facets1, q, variable.id)%>%unique()%>%
  merge(eff.trt.facets, by=c("variable.id"))%>%
  dplyr::rename(p='p-value')%>%dplyr::select(cutoff, q, stability.facets1, terms, Value, Std.Error,  DF,  "t-value", "p", r2m, r2c, n.sites)%>%
  arrange(cutoff, q, stability.facets1)
table(eff.trt.facets.full$q)
eff.trt.facets.full[,6:12]<-round(eff.trt.facets.full[,6:12], 2)

eff.trt.facets.sig<-eff.trt.facets.full%>%filter(p<=0.05  & terms!="(Intercept)" & cutoff==1.28)
# write_xlsx(eff.trt.facets.sig, path="significant treatment effects on stability in hill numbers.xlsx", col_names = TRUE)
# write_xlsx(eff.trt.facets.full, path="full table of treatment effects on stability in hill numbers.xlsx", col_names = TRUE)

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
(pp.trt<-eff.trt4%>%filter(cutoff %in% c(1.28))%>%ggplot(aes(stability.facets1, Value, color=factor(q)))+theme_bw(base_size = 20)+
    geom_point(position=pd, size=5, pch=19, alpha=0.6)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    geom_hline(yintercept = 0, linetype ="dotted") +
    guides(color=guide_legend(nrow=1))+
    theme(legend.position = c(0.8, 0.2), axis.text.x = element_text(angle = 15))+
    labs(x=NULL, y="Effect size of nutrient addition", color="Q"))
# ggsave(pp.trt, width=13.3, height=6.64, dpi=600, file="effects on stability facets in diversity with different Q.png")

# to get the legend for pentagons (the effect size is wrong in this figure)
colour.crosswalk <- c("positive" = "black", "negative" = "red")
(pp.trt.all<-eff.trt4%>%filter((cutoff %in% c(1.28) ))%>%mutate(Value.scaled=abs(Value)*10)%>%
    mutate(direction=ifelse(Value>=0, "positive", "negative"))%>%
    ggplot(aes(stability.facets1, Value.scaled,  colour = as.character(direction), size = (Value.scaled/10)))+theme_bw(base_size = 16)+    geom_point(position=pd, pch=15,alpha=0.6)+
    facet_wrap(~q , nrow=3, scales = "free_y")+
    scale_colour_manual(values = colour.crosswalk) +
    labs(x=NULL, y="Effect size of nutrient addition", size="Effect size", color="Direction"))
legend<-get_legend(pp.trt.all)

# plotting using pentagons
source("functions to draw pentagen.R")
facet_coordinates<-read.csv("facet_coordinates.csv")
colnames(facet_coordinates)[1]<-"facet"
facet_coordinates$Facet<-c("invariability", "resistance_Dry", "resistance_Wet", "recovery_Dry", "recovery_Wet", "Nutrient.addition")

cut<-c(1.28, 0.67)
facet_data <- eff.trt4 %>% filter(cutoff %in% cut)%>% filter(!(stability.facets1=="invariability"))%>%
  mutate(stability.facets2=gsub(".d", "", stability.facets1))%>%
  mutate(stability.facets1= stability.facets2)%>%
  mutate(From_Facet="Nutrient.addition", To_Facet=stability.facets1, value=Value)%>%
  merge(y = facet_coordinates, by.x = "From_Facet", by.y = "Facet") %>%
  merge(y = facet_coordinates, by.x = "To_Facet", by.y = "Facet",
        suffixes = c("", "end")) %>%
  mutate(direction = ifelse(value >= 0, "positive", "negative"), 
         width = as.integer((abs(value) * 10) ))%>%
  mutate(id=paste(cutoff, q, sep="_"))


list_plots <- vector('list', length(unique(facet_data$id)))
for (i in unique(facet_data$id)){
  # i<-"0.67_0"
  facet_data_temp<-facet_data%>%filter(id==i)%>% mutate(sig=case_when(p<=0.05 ~ "significant",  TRUE~"non-significant"))
  facet_data_temp$sig<-factor(facet_data_temp$sig, levels = c("non-significant", "significant"))
  fontsize <- 12 / .pt; size.npk<-16 / .pt
  
  list_plots[[i]]<-gg.pentagon(data = facet_data_temp)+ 
    geom_label(aes( x = 0, y = 0, label = "NPK"), fill="white", size = size.npk)
}

for(i in c(0.67, 1.28)){ 
  # i<-"1.28"
  list_plots_temp1<-list_plots[grep(i, names(list_plots))]
   # rename each element 
  names(list_plots_temp1)<-gsub( "(.*)_(.*)", "\\2", names(list_plots_temp1))
  (pp.trt.effects<-plot_grid(list_plots_temp1$`0`,
                             list_plots_temp1$`1`,
                             list_plots_temp1$`2`,
                              nrow=1,  vjust=7,
                             labels = c("Q = 0", "Q = 1", "Q = 2"), label_fontface = "bold", label_size = 18))
  (pp.trt.effects1<-plot_grid(pp.trt.effects,
                              legend, rel_widths = c(8.8, 1.2)))
 
  ggsave(pp.trt.effects1, height=6.64, width=13.3, file=paste0("effects on stability facets for hill numbers ", i , ".pdf"))
}

# the end
