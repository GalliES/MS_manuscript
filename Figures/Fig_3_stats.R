
#####################################
##### Make all the subsets here #####
#####################################
library(openxlsx)
wb = createWorkbook()


files_names <- list.files(getwd(), pattern='.fcs$', full=FALSE)
files_names

ff <- read.FCS("discovery_cohort.fcs", transformation = F, truncate_max_range = F)
data_all <- exprs(ff)
desc <- description(ff)
colnames(data_all)


# read in panel
panel <- read.xls("FILES/panel_discovery.xlsx", sheet = 1, header = TRUE, verbose = FALSE)

(cytokines <- make.names(as.vector(panel$Antigen[ panel$Category == "cytokines"])))
(cytokines <- cytokines[-6])

# read in metadata
md <- read.xls("FILES/meta_data_discovery_cohort.xlsx", sheet = 1, header = T, verbose = F)
colnames(md)
head(data_all)

####adjust th_labels with reduced data
data <- data_all[data_all[,'manual_labels']==4, ]
dim(data)
th_labels <- readRDS(file = "FILES/th_labels.rds")
length(th_labels)
# select data
data_cd4 <- cbind(data, th_labels)
data_cd4 <- merge(data_cd4, md , by = 'gate_source')
data_cd4 <- droplevels(subset(data_cd4, diagnosis == 'RRMS'| diagnosis == 'NINDC'))



###################################################
##### Demonstrate makeup of selected Th cells #####
###################################################



# select only selected cells
data <- data_cd4[data_cd4$intersection_3_runs == 1,]

# calculate frequencies
fm <- ddply(data, .(gate_source), function(x) {
  Freq <- rep(0, 4)
  for(i in 1:4) {
    div <- nrow(x)
    num <- nrow(x[x[,"cd4_labels"] == i,])
    if(i %in% x[,"cd4_labels"]) {
      Freq[i] <- num/div*100
    } else {Freq[i] <-  0}               
  }
  Freq
})
fm
colnames(fm) <- c("gate_source", "Tnaive", "Teff", "Tem", "Tcm")
fm



# exclude samples with less than 50 cells
N <- ddply(data, .(gate_source), nrow)$V1
fm[N < 10, 2:ncol(fm)] <- NA
fm



# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:5])
head(df_melt)


# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
cluster_names <- make.names(c("Tnaive", "Teff", "Tem", "Tcm"))
fmm$name <- as.factor(cluster_names)
fmm



# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm$name <- factor(fmm$name, levels = fmm$name)
fmm

# make summary pie chart
pie(fmm$median, col = db1[c(4,3,1,6)])


# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = variable, y = value, fill = variable, ymax = 100)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  #facet_wrap('variable', scales = 'fix', ncol = 4) +
  scale_fill_manual(values = db_mem) +
  scale_color_manual(values = db_mem) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
#NINDC=19
#RRMS=27

# save
ggsave(filename = "OUTPUT/Fig_3A_gm_per_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig = 'Fig_3A'

sum <- ddply(df_melt_2, .(variable) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_3A_pop_boxplots_sum.xlsx')


#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$diagnosis, p.adjust.method = "BH", paired = F, correct = T, exact = F)
  with(t, data.frame(method, t(as.vector(p.value)))) })

wres
wres <- melt(wres, id.vars = c("variable",'method' ), measure.vars = colnames(wres)[-c(1,2)])

###decide the * pValue system
wres[,4] <- round(wres[,4], 4)
t <- vector()
for(i in 1:nrow(wres)) {
  if (wres$value[i] <= 0.01)  (t[i] <- "***")
  else if (wres$value[i] <= 0.05)  (t[i] <- "**")
  else if (wres$value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres$stars <- t
wres

write.xlsx(wres, 'STATS/Fig_3A_pop_boxplots_stat.xlsx')

addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)

wres <- ddply(df_melt_2, .(variable), function(x) {
  t2 <- wilcox_test(value ~ diagnosis, data = x, distribution = "exact", conf.int = T)
  p <- pvalue(t2)
  z <- t2@statistic@teststatistic
  r <- as.numeric(z)/sqrt(nrow(x))
  d <- 2*r/sqrt(1-r^2)
  y <- x[!is.na(x$value),]
  cles <- cles.fnc(variable = "value", group = "diagnosis", baseline = "RRMS", 
                   data = y, print = F)
  data.frame(p.value = p, z = -z, r = -r, r2 = r^2, d = -d, cles)
})
wres

write.xlsx(wres, 'STATS/Fig_3A_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################



###### Fig.3B
data_df <- droplevels(subset(data_cd4, cd4_labels == 3))
fm <- ddply(data_df, .(gate_source), function(x)
{df <- data.frame(x) 
(nrow(df[df[,"intersection_3_runs"] ==1 ,])/nrow(df[df[,"gate_source"],]))*100 }
)
colnames(fm)[2] <- "freq"
fm



# check and rename
colnames(fm) <- c("gate_source", c("signature_pop"))
head(fm)



# exclude samples with less than 100 cells
fm$N <- ddply(data_df, .(gate_source), nrow)$V1
fm[fm$N <= 199, 2] <- NA
fm



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
df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  



nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
#NINDC=30
#RRMS=31


# save
ggsave(filename = "OUTPUT/Fig_3B_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig = 'Fig_3B'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_3B_pop_boxplots_sum.xlsx')


#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$diagnosis, p.adjust.method = "BH", paired = F, correct = T, exact = F)
  with(t, data.frame(method, t(as.vector(p.value)))) })

wres
wres <- melt(wres, id.vars = c("variable",'method' ), measure.vars = colnames(wres)[-c(1,2)])

###decide the * pValue system
wres[,4] <- round(wres[,4], 4)
t <- vector()
for(i in 1:nrow(wres)) {
  if (wres$value[i] <= 0.01)  (t[i] <- "***")
  else if (wres$value[i] <= 0.05)  (t[i] <- "**")
  else if (wres$value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres$stars <- t
wres

write.xlsx(wres, 'STATS/Fig_3B_pop_boxplots_stat.xlsx')


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)

wres <- ddply(df_melt_2, .(variable), function(x) {
  t2 <- wilcox_test(value ~ diagnosis, data = x, distribution = "exact", conf.int = T)
  p <- pvalue(t2)
  z <- t2@statistic@teststatistic
  r <- as.numeric(z)/sqrt(nrow(x))
  d <- 2*r/sqrt(1-r^2)
  y <- x[!is.na(x$value),]
  cles <- cles.fnc(variable = "value", group = "diagnosis", baseline = "RRMS", 
                   data = y, print = F)
  data.frame(p.value = p, z = -z, r = -r, r2 = r^2, d = -d, cles)
})
wres

write.xlsx(wres, 'STATS/Fig_3B_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)
##############################################################################


###Fig. 3C
# calculate frequencies for all donors
# select only selected cells
data_cd4$th_labels <- as.factor(data_cd4$th_labels)
data_df <- data_cd4[data_cd4$intersection_3_runs == 1,]
t <- table(data_df[,"gate_source"], data_df[,"th_labels"])
fm <- prop.table(t, 1)*100
apply(fm, 2, median)

# cd4 exclude stuff
t2 <- table(data[,"gate_source"])
excl <- names(which(t2 < 10))
fm[rownames(fm) %in% excl,] <- NA
fm <- data.frame(fm)
fm


# cast and rename
fm <- dcast(fm, Var1 ~ Var2, value.var = "Freq") 
colnames(fm) <-  c("gate_source", 'Tn', "Th1", "Th2", "Th17", "Th22", "ThGM", "Tfh","Treg", "Other")
fm


# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)


# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:10])
head(df_melt)
tail(df_melt)

# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
#cluster_names <- make.names(c("Tnaive", "Teff", "Tem", "Tcm"))
#fmm$name <- as.factor(cluster_names)
fmm



# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
#fmm$name <- factor(fmm$name, levels = fmm$name)
fmm



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = variable, y = value, fill = variable, ymin = 0)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  #facet_wrap('variable', scales = 'fix', ncol = 11) +
  scale_fill_manual(values = db_th) +
  scale_color_manual(values = db_th) +
  #ylim(0,70) +
  theme_bar2
b2  



nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
#NINDC=19
#RRMS=27

# save
ggsave(filename = "OUTPUT/Fig_3C_gm_per_cyt_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig = 'Fig_3C'

sum <- ddply(df_melt_2, .(variable) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_3C_pop_boxplots_sum.xlsx')


#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$diagnosis, p.adjust.method = "BH", paired = F, correct = T, exact = F)
  with(t, data.frame(method, t(as.vector(p.value)))) })

wres
wres <- melt(wres, id.vars = c("variable",'method' ), measure.vars = colnames(wres)[-c(1,2)])

###decide the * pValue system
wres[,4] <- round(wres[,4], 4)
t <- vector()
for(i in 1:nrow(wres)) {
  if (wres$value[i] <= 0.01)  (t[i] <- "***")
  else if (wres$value[i] <= 0.05)  (t[i] <- "**")
  else if (wres$value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres$stars <- t
wres

write.xlsx(wres, 'STATS/Fig_3C_pop_boxplots_stat.xlsx')


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)

wres <- ddply(df_melt_2, .(variable), function(x) {
  t2 <- wilcox_test(value ~ diagnosis, data = x, distribution = "exact", conf.int = T)
  p <- pvalue(t2)
  z <- t2@statistic@teststatistic
  r <- as.numeric(z)/sqrt(nrow(x))
  d <- 2*r/sqrt(1-r^2)
  y <- x[!is.na(x$value),]
  cles <- cles.fnc(variable = "value", group = "diagnosis", baseline = "RRMS", 
                   data = y, print = F)
  data.frame(p.value = p, z = -z, r = -r, r2 = r^2, d = -d, cles)
})
wres

write.xlsx(wres, 'STATS/Fig_3C_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)
##############################################################################



##Save excel workbook
saveWorkbook(wb, "STATS/Figure_3.xlsx")



