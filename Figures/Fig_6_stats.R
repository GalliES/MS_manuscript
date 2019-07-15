######################################
##### load all neccessary things #####
######################################
library(openxlsx)
wb = createWorkbook()


# load all functions and themes
setwd("/Users/Edoardo Galli/Documents/MS_daclizumab/FILE/Combina/")

source("functions_library_colorcode.R")

files_names <- list.files(getwd(), pattern='.fcs$', full=FALSE)
files_names


##### ##### ##### ##### ##### 
##### Stats playing ##### 
##### ##### ##### ##### ##### 
freq <- read.xlsx('FILES/CNS_cohort.xlsx')
freq
# read in metadata
md <- read.xlsx("FILES/meta_data_CNS_cohort.xlsx")
md
tail(md)
fm <- merge(freq, md , by = 'patient')

#write.xlsx(fm, 'CNS_data_clean.xlsx')

# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)
tail(df_plot)



colnames(df_plot)
cor_plot <- droplevels(subset(df_plot, diagnosis %in% c("ms")))
#write.xlsx(cor_plot, 'corplot.xlsx')
cor_plot
table(cor_plot$pt_id)

heat <- vector( length = 9)
for(i in 2:10) {
  cor_melt <- melt(cor_plot, 'week', measure.vars = colnames(df_plot[,c(2:10)]))
  g <- subset(cor_melt, variable == colnames(df_plot)[i])
  p <- subset(g,week == 'PBMC')
  c <- subset(g, week == 'CSF')
  print(cor(p[,'value'],c[,'value']))
  heat_mat[i,] <- as.vector(cor(p[,'value'],c[,'value']))
  }


# melt for ggplot
df_melt <- melt(cor_plot, measure.vars = colnames(cor_plot)[2:18])
head(df_melt)
tail(df_melt)
str(df_melt)


# ordering and piecharts
data_melt_drop <- droplevels(subset(df_melt, diagnosis %in% c("ms") & patient != 21))


wres_all <- ddply(data_melt_drop, .(variable), function(x) {
  before <- subset(x, week == 'PBMC')
  before <- before[c(as.vector(order(before$patient, decreasing = F))),]
  after <-subset(x, week == 'CSF')
  after <- after[c(as.vector(order(after$patient, decreasing = F))),]
  res <- wilcox.test(before[,'value'], after[,'value'],p.adjust.method = "BH", paired = TRUE)
  res$p.value
  with(res, data.frame(method, t(as.vector(p.value))))
})
wres_all


###adjust p.value
wres_all$BH.p.value <- p.adjust(wres_all$t.as.vector.p.value.., method = "BH")
wres_all



# add stars
t <- vector()
for(i in 1:nrow(wres_all)) {
  if (is.na(wres_all$BH.p.value[i])) (t[i] <- "")
  else if (wres_all$BH.p.value[i] <= 0.01)  (t[i] <- "***")
  else if (wres_all$BH.p.value[i] <= 0.05)  (t[i] <- "**")
  else if (wres_all$BH.p.value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres_all$stars <- t
wres_all$p.value <- round(wres_all$p.value, 3)
wres_all

write.xlsx(wres_all, 'STATS/Fig_6_stat.xlsx')



# make p-value data.frame
wp <- wres_all
wp$p.value <- paste("P = ", round(wres_all$BH.p.value, 4), t,  sep = "")
wp$x <- 1.5  
wp$value <- ddply(data_melt_drop, .(variable), function(x) {max(x$value, na.rm = T)})$V1
wp


# make an overview plot
comb <- ggplot(data_melt_drop, aes(x = week, y = value, ymin = 0, ymax = value*1.2)) +
  #geom_line(linetype = "solid",  size= 0.5, aes(group = patient_num))+
  geom_point(size = 4, aes(color = variable)) +
  geom_boxplot(aes(fill = variable), lwd = 0.5, color = "white", outlier.shape = NA) +
  stat_boxplot(geom = 'errorbar', aes(color = variable)) +
  facet_wrap(~ variable, scales = "free", nrow=2) +
  #geom_text(data = wp, aes(y = value*1.3, x = x, label = p.value, group = NULL), size = 3) + 
  #geom_blank(data=data_dummy) +
  #geom_text(data = wp, aes(y = value*1.2, x = x, label = p.value.nindc, group = NULL), size = 3) +  
  #geom_text(data = wp, aes(y = value*1.1, x = x, label = p.value.hc, group = NULL), size = 3) +      
  theme_bar2 + 
  theme(strip.text = element_text(size = rel(1.5))) +
  theme(panel.spacing = unit(2.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, color = "black"))
comb


# save plots
#ggsave(filename = "OUTPUTnat/temp_cd4_cytokinepositivity_all_boxplot.pdf", plot = comb, 
#       useDingbats = F, scale = 1.75, width = 7*1.65, height = 2*3.5, units = c("in"))



# rename and round
names(wres_all)[3]  <- c("p.value")
wres_all[,3] <- round(wres_all[,3], 4)
wres_all


### Summary
sum <- ddply(data_melt_drop, .(variable, diagnosis, week) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum


write.xlsx(sum, 'STATS/Fig_6_sum.xlsx')


###adjust p.value
wres_all$BH.p.value <- p.adjust(wres_all$p.value, method = "BH")
wres_all



# add stars
t <- vector()
for(i in 1:nrow(wres_all)) {
  if (is.na(wres_all$BH.p.value[i])) (t[i] <- "")
  else if (wres_all$BH.p.value[i] <= 0.01)  (t[i] <- "***")
  else if (wres_all$BH.p.value[i] <= 0.05)  (t[i] <- "**")
  else if (wres_all$BH.p.value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres_all$stars <- t
wres_all$p.value <- round(wres_all$p.value, 3)
wres_all

# make p-value data.frame
wp <- wres_all
wp$p.value <- paste("P = ", round(wres_all$BH.p.value, 4), t,  sep = "")
wp$x <- 1.5  
wp$value <- ddply(data_melt_drop, .(variable), function(x) {max(x$value, na.rm = T)})$V1
wp


str(data_melt_drop)
data_melt_drop$week <- as.factor(data_melt_drop$week)
effect_size <- ddply(data_melt_drop, .(variable) , function(x) {
  t2 <- wilcox_test(value ~ week, data = x, distribution = "exact", conf.int = T)
  p <- pvalue(t2)
  z <- t2@statistic@teststatistic
  r <- as.numeric(z)/sqrt(nrow(x))
  d <- 2*r/sqrt(1-r^2)
  y <- x[!is.na(x$value),]
  cles <- cles.fnc(variable = "value", group = "week", baseline = "PBMC", 
                   data = y, print = F)
  data.frame(p.value = p, z = -z, r = -r, r2 = r^2, d = -d, cles)
})
effect_size

write.xlsx(effect_size, 'STATS/Fig_6_effect_size.xlsx')



fig <- 'Fig_6'
addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres_all)
addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =effect_size)


##############################################################################



##Save excel workbook
saveWorkbook(wb, "STATS/Figure_6.xlsx")


