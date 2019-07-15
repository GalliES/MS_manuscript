
#####################################
##### Make all the subsets here #####
#####################################
#install.packages("openxlsx")
library(openxlsx)
wb = createWorkbook()


files_names <- list.files(getwd(), pattern='.fcs$', full=FALSE)
files_names

ff <- read.FCS("discovery_cohort.fcs", transformation = F, truncate_max_range = F)
data_all <- exprs(ff)
desc <- description(ff)
colnames(data_all)

dim(data_all)
head(data_all)

# read in panel
panel <- read.xls("FILES/panel_discovery.xlsx", sheet = 1, header = TRUE, verbose = FALSE)

(cytokines <- make.names(as.vector(panel$Antigen[ panel$Category == "cytokines"])))
(cytokines <- cytokines[-6])

# read in metadata
md <- read.xls("FILES/meta_data_discovery_cohort.xlsx", sheet = 1, header = T, verbose = F)
colnames(md)
data_all <- merge(data_all, md[,c(2,3,6,9)], by= 'gate_source')

###boxplot
comb <- ggplot(md, aes(x = diagnosis.general, y = age, ymin = 0, ymax = max(md$age))) +
  #geom_line(linetype = "solid",  size= 0.5, aes(group = patient_num))+
  geom_point(size = 4, aes(color = diagnosis.general)) +
  geom_boxplot(aes(fill = diagnosis.general), lwd = 0.5, color = "white", outlier.shape = NA) +
  stat_boxplot(geom = 'errorbar', aes(color = diagnosis.general)) +
  #facet_wrap(~ diagnosis.general, scales = "fixed", nrow=1) +
  #geom_text(data = wp, aes(y = value*1.3, x = x, label = p.value, group = NULL), size = 2) + 
  #geom_blank(data=data_dummy) +
  #geom_text(data = wp, aes(y = value*1.2, x = x, label = p.value.nindc, group = NULL), size = 3) +  
  #geom_text(data = wp, aes(y = value*1.1, x = x, label = p.value.hc, group = NULL), size = 3) +      
  theme_bar2 + 
  theme(strip.text = element_text(size = rel(1.5))) +
  theme(panel.spacing = unit(2.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, color = "black"))
comb



# save plots
ggsave(filename = "Fig.S3A_age.pdf", plot = comb, 
       useDingbats = F, scale = 1, width = 3*1.65, height = 3.5, units = c("in"))


t <- pairwise.wilcox.test(data = md , md$age, g = md$diagnosis.general, p.adjust.method = "BH", paired = F, correct = T, exact = F)
t$p.value


# select data & calculate frequencies
data_df <- data_all
data_melt <- melt(data_df, id.vars = "collaborator_id", measure.vars = cytokines)
cyto_freq <- ddply(data_melt, .(collaborator_id, variable), function(x){
  Freq <- length(x$value[x$value >= 0.325])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})


colnames(cyto_freq)[3] <- 'value'
df_melt <- cyto_freq

head(df_melt)
max(subset(df_melt, variable == 'IFN.g')$value)
write.xlsx(df_melt, 'extfig3_cyo.xlsx')



