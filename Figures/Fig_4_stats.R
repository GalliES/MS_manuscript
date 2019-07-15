
#####################################
##### Make all the subsets here #####
#####################################
library(openxlsx)
wb = createWorkbook()

files_names <- list.files(getwd(), pattern='.fcs$', full=FALSE)
files_names

ff <- read.FCS("validation_cohort.fcs", transformation = F, truncate_max_range = F)
data_all <- exprs(ff)
desc <- description(ff)
colnames(data_all)

head(data_all)
dim(data_all)

# read in panel
panel <- read.xls("FILES/panel_validation.xlsx", sheet = 1, header = TRUE, verbose = FALSE)

(cytokines <- make.names(as.vector(panel$Antigen[ panel$Category == "cytokines"])))

# read in metadata
md <- read.xls("FILES/meta_data_validation_cohort.xlsx", sheet = 1, header = T, verbose = F)
colnames(md)
head(data_all)
data_all <- data_all[data_all[,'manual_labels'] == 1,]


######################################
##### CellCNN_gating_CIS_cohort ######
######################################


fm <- ddply(data.frame(data_all), .(patient), function(x)
{df <- data.frame(x) 
(nrow(df[df[,"intersection_3_runs"] ==1 ,])/nrow(df[df[,"patient"],]))*100 }
)
colnames(fm)[2] <- "freq"
fm


# check and rename
colnames(fm) <- c("patient", c("signature_pop"))
head(fm)


# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2])
head(df_melt)



# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm



# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
#df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[c(4,1,6,3,5,2,7)])
df_melt_2 <- merge(df_melt, md , by = 'patient')

# plot boxplots
df_melt_2$Diagnosis <- factor(df_melt_2$Diagnosis, levels = levels(df_melt_2$Diagnosis)[c(2,1,4,5,3)])

b2 <- ggplot(data = df_melt_2, aes(x = Diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  

# save
ggsave(filename = "OUTPUT/Fig_4B_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary
fig <- 'Fig_4B'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum

write.xlsx(sum, 'STATS/Fig_4B_pop_boxplots_sum.xlsx')


#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$Diagnosis, p.adjust.method = "BH", paired = F, correct = T, exact = F)
  with(t, data.frame(method, t(as.vector(p.value)))) })
wres

colnames(wres) [c(3:6,8:10,13:14,18)] <- c('HCvsCIS','HCvsMS','HCvsNINDC',
                                                 'HCvsINDC','CISvsMS','CISvsNINDC',
                                           'CISvsINDC','MSvsNINDC','MSvsINDC','NINDCvsINDC')
wres <- wres[,c(1,2,3:6,8:10,13:14,18)]

wres
wres <- melt(wres, id.vars = c("variable",'method' ), measure.vars = colnames(wres)[-c(1,2)])

###decide the * pValue system
wres[,4] <- round(wres[,4], 6)
t <- vector()
for(i in 1:nrow(wres)) {
  if (wres$value[i] <= 0.01)  (t[i] <- "***")
  else if (wres$value[i] <= 0.05)  (t[i] <- "**")
  else if (wres$value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres$stars <- t
wres

write.xlsx(wres, 'STATS/Fig_4B_pop_boxplots_stat.xlsx')

addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)



#####################################
##### Manual_gating_CIS_cohort ######
#####################################

fm <- ddply(data.frame(data_all), .(patient), function(x)
{df <- data.frame(x) 
(nrow(df[df[,"CXCR4"] > 0.4 & df[,"IFN.g"] > 0.25 & df[,"GM.CSF"] > 0.25  & df[,"TNF.a"] > 0.25& df[,"IL.2"] > 0.25,])/nrow(df[df[,"patient"],]))*100 }
)
colnames(fm)[2] <- "freq"
fm


# check and rename
colnames(fm) <- c("patient", c("signature_pop"))
head(fm)


# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2])
head(df_melt)



# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm



# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
#df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[c(4,1,6,3,5,2,7)])
df_melt_2 <- merge(df_melt, md , by = 'patient')

# plot boxplots
df_melt_2$Diagnosis <- factor(df_melt_2$Diagnosis, levels = levels(df_melt_2$Diagnosis)[c(2,1,4,5,3)])

b2 <- ggplot(data = df_melt_2, aes(x = Diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  

# save
ggsave(filename = "OUTPUT/Fig_4E_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary
fig <- 'Fig_4E'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum

write.xlsx(sum, 'STATS/Fig_4E_pop_boxplots_sum.xlsx')


#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$Diagnosis, p.adjust.method = "BH", paired = F, correct = T, exact = F)
  with(t, data.frame(method, t(as.vector(p.value)))) })
wres

colnames(wres) [c(3:6,8:10,13:14,18)] <- c('HCvsCIS','HCvsMS','HCvsNINDC',
                                           'HCvsINDC','CISvsMS','CISvsNINDC',
                                           'CISvsINDC','MSvsNINDC','MSvsINDC','NINDCvsINDC')
wres <- wres[,c(1,2,3:6,8:10,13:14,18)]

wres
wres <- melt(wres, id.vars = c("variable",'method' ), measure.vars = colnames(wres)[-c(1,2)])

###decide the * pValue system
wres[,4] <- round(wres[,4], 6)
t <- vector()
for(i in 1:nrow(wres)) {
  if (wres$value[i] <= 0.01)  (t[i] <- "***")
  else if (wres$value[i] <= 0.05)  (t[i] <- "**")
  else if (wres$value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres$stars <- t
wres
write.xlsx(wres, 'STATS/Fig_4E_pop_boxplots_stat.xlsx')


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)

##############################################################################



##Save excel workbook
saveWorkbook(wb, "STATS/Figure_4.xlsx")
