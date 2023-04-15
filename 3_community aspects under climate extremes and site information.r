###########################################################################################
###continue from line 224 of 2_categorize climate events and information for sites used.R ##############
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
# ggsave(pp.each.site.bio, dpi=600,width=13.3, height=6.64,  file="above ground biomass at each site using cutoff of 1.28.png")
 ggsave("pp.each.site.bio.pdf", width = 21, height = 29.7, dpi=600)


(pp.each.site.rich<-ggplot(climate.extreme.cut0.67, aes(year_trt, richness, color=climate, shape=trt))+theme_cowplot(font_size = 20)+
    facet_wrap(~site_code, scales="free_y", ncol=7)+
    stat_summary(fun=mean, geom="point", size=2, position=pd, alpha=0.6)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(climate.extreme.cut0.67_each.site_avg, mapping=aes(yintercept =avg_rich, linetype=trt), alpha=0.6)+
    scale_colour_manual(values = colour.crosswalk) +
    scale_x_continuous(expand=c(0,0), limits = c(0,15), breaks = seq(0,15,2))+
    theme(axis.text.y = element_text(angle=30))+
     labs(x="Years after nutrient addition", y="Species richness", color="Climate\n extremes", shape="Treatments", linetype="Treatments"))
# ggsave(pp.each.site.rich, dpi=600,width=13.3, height=6.64,  file="richness at each site using cutoff of 1.28.png")
 ggsave("pp.each.site.rich.pdf", width = 21, height = 29.7, dpi=600)
 
# focus on site look.us 
look_each.block_avg<-climate.extreme.cut0.67%>%filter(site_code=="look.us" & climate=="Normal")%>%
  group_by(site_code, block, trt)%>%
   dplyr::summarise(avg_bio=mean(live_mass), avg_rich=mean(richness))
 
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
 ggsave(pp.look.bio, dpi=600,width=13.3, height=6.64,  file="above ground biomass at site look.us.png")
 
(pp.look.rich<-climate.extreme.cut0.67 %>%filter(site_code=="look.us")%>%
   select(block, trt, year_trt, richness, climate, spei)%>%
    pivot_longer(cols = c("richness", "spei"))%>%mutate(block1=ifelse(name=="richness",  block, "SPEI"), trt1=ifelse(name=="richness",  trt, "SPEI"))%>%
    merge(y=look_each.block_avg, by=c("block", "trt"))%>%mutate(avg_rich1=ifelse(name=="richness", avg_rich, NA))%>%
    ggplot(aes(year_trt, value))+theme_bw(base_size = 18)+
    facet_wrap(~block1, nrow=4, scales="free_y")+
    geom_point(size=3, position=pd, alpha=0.6, aes(color=climate, shape=trt1))+geom_line(aes(group=trt1), position=pd, alpha=0.2)+
    geom_hline(aes(yintercept =avg_rich1, linetype=trt), alpha=0.6)+
    scale_colour_manual(values = colour.crosswalk) +
    scale_x_continuous(expand=c(0,0), limits = c(0,15), breaks = seq(0,15,1))+
    theme(axis.text.y = element_text(angle=30))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x="Years after nutrient addition", y="Species richness", color="Climate\n extremes", shape="Treatments", linetype="Treatments"))
# ggsave(pp.look.rich, dpi=600,width=13.3, height=6.64,  file="Species richness at site look.us.png")

 
###########################################################################################
# biomass under different climate events under and one year after climate events using cutoff of 0.67#
###########################################################################################
 
bio.rich<- rich.bio.each.functional.group%>%filter(life.form=="all")%>%filter(community.property %in%c("biomass", "richness"))%>%
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
## add average biomass During different treatments different blocks 
climate.extreme.cut0.67_block_avg<-d.select%>%filter(climate!="Normal")%>%filter(life.form=="all")%>%filter(community.property %in%c("biomass", "richness"))%>%
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
 ggsave(pp.u.o, dpi=600, width=13.3, height=6.64, file="aboveground biomass During and one year after.png")

(pp.rich.u.o<-climate.extreme.cut0.67.au%>%filter(community.property=="richness")%>%
    mutate(name0=ifelse(name=="dif_property", "Change in richness", "Deviation in richness"), name1=as.factor(name0), name1=relevel(name1, ref="Change in richness"))%>%
    ggplot( aes(events, value, color=trt))+theme_cowplot(font_size = 20)+panel_border()+
    facet_grid(name1~climate1, scale="free_y", space="free_x")+
    stat_summary(fun=mean, geom="point", size=5, position=pd)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    geom_hline(aes(yintercept=reference.line), linetype="dashed", linewidth=1, alpha=0.6)+
    theme(axis.text.x = element_text(angle=35, hjust=1))+
     labs(x=NULL, y=NULL, color=NULL))
 ggsave(pp.rich.u.o, dpi=600, width=13.3, height=6.64, file="Species richness During and one year after.png")

###########################################################################################
#######compare community dissimilarity during and one year after climate events ###########
###########################################################################################
composition<-d8%>%filter(life.form!="OTHER")%>%mutate(life.form="all")%>%
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

# Fig. S9
(pp.composition.u.o<-ggplot(composition.au, aes(events, sim, color=trt))+
    theme_cowplot(font_size = 20)+panel_border()+
    theme(axis.text.x = element_text(angle=35, hjust=1))+
    facet_wrap(~climate1, scale="free_x")+
    stat_summary(fun=mean, geom="point", size=5, position=pd)+
    stat_summary(fun.data =mean_cl_boot, geom="errorbar", width=0.1, position=pd)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x=NULL, y="Community similarity", color=NULL))
 ggsave(pp.composition.u.o, dpi=600, width=13.3, height=6.64, file="composition similarity During and one year after.png")

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

##### Read the data for geolocation 
all.s <- read.csv(paste0(dir.data, 'comb-by-plot-clim-soil-2022-11-15.csv'), header = T)
colnames(all.s)
## select variables needed 
## merge with the data for water balance and experimental years used 
site.inf<-all.s%>%select(site_code, habitat, continent, country, latitude, longitude, first_nutrient_year)%>%distinct()%>%
  merge(y=wat.bal1[,c("site_code", "growing.season", "water.balance.last.15")], by="site_code")%>%
  merge(y=exp.yr[,c("site_code", "years.used")], by=c("site_code"))%>%
  merge(dw[, c("site_code", "climate2")], group_by=c("site_code"))
# seems that first nutrient year is not always the same to first fence year 
site.inf1<-site.inf%>%mutate(first_fenced_year=NULL)%>%dplyr::rename(first.experimental.year=first_nutrient_year)
site.inf1$water.balance.last.15<-round(site.inf1$water.balance.last.15, 2)
site.inf1$latitude<-round(site.inf1$latitude, 2)
site.inf1$longitude<-round(site.inf1$longitude, 2)
# table S1
 write_xlsx(site.inf1, path="information for sites used.xlsx")

###########################################################################################
#duration and climate1 events sites deleting sites without biomass but keep sites and years that dry and wet happen consecutively#
###########################################################################################
# doing this is to keep a balanced view about the dry and wet climate1 extremes 
# using cutoff of 0.67
dd3<-s_dd7%>%group_by(site_code)%>%summarise(num_yr=length(year_trt))%>%distinct()%>%
  merge(dw[, c("site_code", "climate2")], group_by=c("site_code"))%>%distinct()%>%
  merge(y=exp.yr[,c("site_code", "years.used")], by=c("site_code"))
range(dd3$num_yr)
(pp<-ggplot(dd3)+geom_histogram(aes(x=num_yr))+facet_wrap(~climate2)+theme_bw(base_size=18)+
    #scale_x_continuous(expand=c(0,0), limits = c(3,15), breaks = seq(1,15,3))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x="Experimental duration", y="Number of sites"))

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
 ggsave(pp.map, file="map1.png", height = 6.64, width = 13.3, dpi = 600)
# Fig. S1
(pp.map.du<-plot_grid(pp.map, pp, nrow=2, rel_heights = c(4, 2)))
 ggsave(pp.map.du, height = 6.64, width = 13.3, file="duration and experimental sites.png")

###########################################################################################
############################ check site PIs for acknowledgement### ########################
###########################################################################################
site.pis<-read.csv(paste0(dir.data, 'pi-contact-list-8-March-2021.csv'), header = T)
authors<-readxl::read_xlsx("author contribution.xlsx")%>%select(1,2,3, 4,5, 7)
colnames(authors)<-c("Order", "First.Name",  "Last.Name",   "Full.name", "Affiliations",  "Email")

authors1<-authors%>%mutate(pi.names=toupper(paste(First.Name, Last.Name)))%>%
  mutate(author.ordered=paste(Full.name, Order), affiliation.ordered=paste(Order, Affiliations))%>%arrange(Order)

site.pis1<-site.pis%>%filter(site_code%in%dd3$site_code)%>%
  mutate(PI.name=paste(firstname, lastname, sep=" "), PI.name1=toupper(PI.name))%>%
  filter(!PI.name1 %in% c(authors1$pi.names))%>%
  filter(!PI.name %in%c("Eric Seabloom", "W. Harpole", "Anu Eskelinen", "Pedro Daleo"))%>%
  select(site_code, PI.name, institution, email)%>%arrange(site_code)

# table. S6
 write_xlsx(site.pis1, path="site PIs for sites used.xlsx")  

# check for those whose sites were included due to use cutoff of 0.67 sd in the main text
site.cut.1.28<-all.data.l_1%>%filter(cutoff==1.28)%>%select(site_code)%>%distinct()
site.pis.remind<-site.pis%>%filter(site_code%in%dd3$site_code)%>%
  filter(!site_code %in% site.cut.1.28$site_code)%>%
  select(site_code, firstname, lastname, institution, email)%>%arrange(site_code)
email.list<-(site.pis.remind$email)
email.list1<-paste(email.list, collapse =";")
