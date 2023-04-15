
###########################################################################################
############# construct all possible combinations for extreme growing seasons ##############
###########################################################################################

all.combi<-expand.grid(One.year.before=c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"), 
                       an.extreme.season=c("Extreme dry", "Extreme wet"), 
                       One.year.after=c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

all.combi.num<-expand.grid(One.year.before.num=c(-10, -5, 0, 5, 10), 
                           an.extreme.season.num=c(-10, 10), 
                           One.year.after.num=c(-10, -5, 0, 5, 10))

all.combi.1<-all.combi%>%bind_cols(all.combi.num)%>%
  arrange(an.extreme.season, One.year.before, One.year.after)%>%
  mutate(climate.type.later=an.extreme.season.num+One.year.after.num,
         climate.type.previous= One.year.before.num+an.extreme.season.num)%>%
  mutate(year.resis=ifelse(climate.type.previous%in%c(0, -5, 5), "N", "Y"))%>%
  mutate(year.recov=ifelse(climate.type.later%in%c(0, -5, 5), "N", year.resis))%>%
  mutate(year.recov1=ifelse(abs(One.year.after.num)<abs(an.extreme.season.num),  year.recov, "N"))%>%
  select(One.year.before, an.extreme.season, One.year.after, year.resis, year.recov1)%>%dplyr::rename(year.recov=year.recov1)

# writexl::write_xlsx(all.combi.1, path="combination of extreme events and years used for resistance and recovery.xlsx")


#############################################################################################
##### draw examples of treatments effects on stability facets and their relationships ########
#############################################################################################
stability.facets.list<-c("invariability", "resistance_Dry", "resistance_Wet", "recovery_Dry", "recovery_Wet")

# plotting using pentagons
source("functions to draw pentagen.R")
fontsize <- 18 / .pt

facet_coordinates<-read.csv("facet_coordinates.csv")
colnames(facet_coordinates)[1]<-"facet"
facet_coordinates$Facet<-c(stability.facets.list, "Nutrient.addition")
# treatment effects 
facet_values<-data.frame(c(	rep("Nutrient.addition", 5)),
                         stability.facets.list,
                         c(rep(0.2, 5)),
                         c(rep("A", 5)))
colnames(facet_values)<-c("From_Facet","To_Facet","value", "example")

(fig.a<-facet_values %>%
    merge(y = facet_coordinates, by.x = "From_Facet", by.y = "Facet") %>%
    merge(y = facet_coordinates, by.x = "To_Facet", by.y = "Facet",
          suffixes = c("", "end")) %>%
    mutate(direction = ifelse(value >= 0, "positive", "negative"), 
           width = as.integer((abs(value) * 5) + 1), sig="marginal-significant")%>%
    gg.pentagon()+
    geom_label(aes( x = 0, y = 0, label = "NPK"), fill="white", size = 20 / .pt))

# relationships 
facet_values<-as.data.frame(t(combn(stability.facets.list, 2)))
facet_values$value<-0.2
colnames(facet_values)<-c("From_Facet","To_Facet","value")

(fig.b<-facet_values %>%
    merge(y = facet_coordinates, by.x = "From_Facet", by.y = "Facet") %>%
    merge(y = facet_coordinates, by.x = "To_Facet", by.y = "Facet",
          suffixes = c("", "end")) %>%
    mutate(direction = ifelse(value >= 0, "positive", "negative"), 
           width = as.integer((abs(value) * 5) + 1), sig="marginal-significant")%>%
    gg.pentagon())

# relationships between aspects 
aspect_coordinates<-facet_coordinates%>%filter(facet %in% c("Facet 1", "Facet 4", "Facet 5"))%>%
  mutate(Aspect=c("biomass",  "richness", "composition"), facet=NULL, Facet=NULL,
         x=c(0, -0.866, 0.866), y=c(1, -0.5, -0.5))

aspect_values<-as.data.frame(t(combn( c("biomass", "richness", "composition"), 2)))
aspect_values$value<-0.2
colnames(aspect_values)<-c("From_aspect","To_aspect","value")

(fig.c<-aspect_values %>%
    merge(y = aspect_coordinates, by.x = "From_aspect", by.y = "Aspect") %>%
    merge(y = aspect_coordinates, by.x = "To_aspect", by.y = "Aspect",
          suffixes = c("", "end")) %>%
    mutate(direction = ifelse(value >= 0, "positive", "negative"), 
           width = as.integer((abs(value) * 5) + 1), sig="marginal-significant")%>%
    gg.triangle())

(all.plots<-plot_grid(fig.a, fig.b, fig.c, labels = "AUTO", 
                      nrow=1, vjust = 8.5,  label_fontface = "bold", label_size = 18 ))
ggsave(all.plots, height=6.64, width=13.3, dpi=600, file="components for conceptual figure.png")


