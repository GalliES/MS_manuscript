
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


######################## select data

data_df <- data.frame(data_all)

# calculate frequencies
fm <- ddply(data_df, .(patient), function(x) {
  Freq <- rep(0, 6)
  for(i in 1:6) {
    div <- nrow(x)
    num <- nrow(x[x[,"manual_labels"] == i,])
    if(i %in% x[,"manual_labels"]) {
      Freq[i] <- num/div*100
    } else {Freq[i] <-  0}               
  }
  Freq
})
fm



# check and rename
colnames(fm) <- c("patient", c("CD4", "CD8", "gdT.cells", "B.cell", "NK", "Myeloid"))
head(fm)



# exclude samples with less than 100 cells
fm$N <- ddply(data_df, .(patient), nrow)$V1
fm[fm$N <= 50, 2:8] <- NA
fm



# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:7])
head(df_melt)


# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
cluster_names <- make.names(c("CD4", "CD8", "gdT.cells", "B.cell", "NK", "Myeloid"))
fmm$name <- as.factor(cluster_names)
fmm



# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm$name <- factor(fmm$name, levels = fmm$name)
fmm



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])
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


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'patient')
table(t$Diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig.S2D_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_S8D'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum


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


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)


##############################################################################


####Fig. S8F-I
######################## select data

data_df <- data.frame(data_all)
data_df <- subset(data_df, manual_labels == 1)
data_cd4 <- merge(data_df, md, by = 'patient')
data.us <- data.frame(data_all[data_all[,'patient']==73,])
data_cd4_us <- subset(data.us, manual_labels == 1)


# make cytokine vector
(cytokines <- make.names(as.vector(panel$Antigen[panel$CD4 == "yes" & panel$Category == "cytokines"])))

# select data & calculate frequencies
per.vector <- apply(data_cd4_us[,cytokines], 2, function(x) quantile(x, 0.995, names = F))
per.vector[c(1,2,3,6,9,10,12,13)] <- c(0.4,0.25,0.4,0.4,0.25,0.25,0.25,0.4)

data_df <- data_cd4
#data_df <- data.frame(merge(data_cd4,md,by='patient'))
data_melt <- melt(data_df, id.vars = c("patient"), 
                  measure.vars = cytokines)

cyto_freq <- ddply(data_melt, .(patient, variable), function(x){
  Freq <- length(x$value[x$value >= as.numeric(per.vector[x$variable])])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})

### Merge md
fm <- dcast(cyto_freq, patient+N ~ variable, value.var = "freq")
#write.csv2(fm, "cytokines_pos_0.4_1.csv")

# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
# melt for ggplot
df_melt <- melt(df_plot, measure.vars = cytokines)

# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm$name <- as.factor(cytokines)
fmm

# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm$name <- factor(fmm$name, levels = fmm$name)
fmm

# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])
df_melt_2 <- merge(df_melt, md , by = 'patient')

# plot boxplots
df_melt_2$Diagnosis <- factor(df_melt_2$Diagnosis, levels = levels(df_melt_2$Diagnosis)[c(2,1,4,5,3)])

b2 <- ggplot(data = df_melt_2, aes(x = Diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'patient')
table(t$Diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig.S8F_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_S8F'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum


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


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)


##############################################################################

######################## select data

data_df <- data.frame(data_all)
data_df <- subset(data_df, manual_labels == 2)
data_cd4 <- merge(data_df, md, by = 'patient')
data.us <- data.frame(data_all[data_all[,'patient']==73,])
data_cd4_us <- subset(data.us, manual_labels == 2)


# make cytokine vector
(cytokines <- make.names(as.vector(panel$Antigen[panel$CD4 == "yes" & panel$Category == "cytokines"])))

# select data & calculate frequencies
per.vector <- apply(data_cd4_us[,cytokines], 2, function(x) quantile(x, 0.995, names = F))
per.vector[c(1,2,3,6,9,10,12,13)] <- c(0.4,0.25,0.4,0.4,0.25,0.25,0.25,0.4)

data_df <- data_cd4
#data_df <- data.frame(merge(data_cd4,md,by='patient'))
data_melt <- melt(data_df, id.vars = c("patient"), 
                  measure.vars = cytokines)

cyto_freq <- ddply(data_melt, .(patient, variable), function(x){
  Freq <- length(x$value[x$value >= as.numeric(per.vector[x$variable])])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})

### Merge md
fm <- dcast(cyto_freq, patient+N ~ variable, value.var = "freq")
#write.csv2(fm, "cytokines_pos_0.4_1.csv")

# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
# melt for ggplot
df_melt <- melt(df_plot, measure.vars = cytokines)

# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm$name <- as.factor(cytokines)
fmm

# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm$name <- factor(fmm$name, levels = fmm$name)
fmm

# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])
df_melt_2 <- merge(df_melt, md , by = 'patient')

# plot boxplots
df_melt_2$Diagnosis <- factor(df_melt_2$Diagnosis, levels = levels(df_melt_2$Diagnosis)[c(2,1,4,5,3)])

b2 <- ggplot(data = df_melt_2, aes(x = Diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'patient')
table(t$Diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig.S8G_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_S8G'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum


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


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)


##############################################################################



######################## select data

data_df <- data.frame(data_all)
data_df <- subset(data_df, manual_labels == 5)
data_cd4 <- merge(data_df, md, by = 'patient')
data.us <- data.frame(data_all[data_all[,'patient']==73,])
data_cd4_us <- subset(data.us, manual_labels == 5)


# make cytokine vector
(cytokines <- make.names(as.vector(panel$Antigen[panel$CD4 == "yes" & panel$Category == "cytokines"])))

# select data & calculate frequencies
per.vector <- apply(data_cd4_us[,cytokines], 2, function(x) quantile(x, 0.995, names = F))
per.vector[c(1,2,3,6,9,10,12,13)] <- c(0.4,0.25,0.4,0.4,0.25,0.25,0.25,0.4)

data_df <- data_cd4
#data_df <- data.frame(merge(data_cd4,md,by='patient'))
data_melt <- melt(data_df, id.vars = c("patient"), 
                  measure.vars = cytokines)

cyto_freq <- ddply(data_melt, .(patient, variable), function(x){
  Freq <- length(x$value[x$value >= as.numeric(per.vector[x$variable])])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})

### Merge md
fm <- dcast(cyto_freq, patient+N ~ variable, value.var = "freq")
#write.csv2(fm, "cytokines_pos_0.4_1.csv")

# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
# melt for ggplot
df_melt <- melt(df_plot, measure.vars = cytokines)

# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm$name <- as.factor(cytokines)
fmm

# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm$name <- factor(fmm$name, levels = fmm$name)
fmm

# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])
df_melt_2 <- merge(df_melt, md , by = 'patient')

# plot boxplots
df_melt_2$Diagnosis <- factor(df_melt_2$Diagnosis, levels = levels(df_melt_2$Diagnosis)[c(2,1,4,5,3)])

b2 <- ggplot(data = df_melt_2, aes(x = Diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'patient')
table(t$Diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig.S8H_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_S8H'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum


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


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)


##############################################################################


######################## select data

data_df <- data.frame(data_all)
data_df <- subset(data_df, manual_labels == 4)
data_cd4 <- merge(data_df, md, by = 'patient')
data.us <- data.frame(data_all[data_all[,'patient']==73,])
data_cd4_us <- subset(data.us, manual_labels == 4)


# make cytokine vector
(cytokines <- make.names(as.vector(panel$Antigen[panel$CD4 == "yes" & panel$Category == "cytokines"])))

# select data & calculate frequencies
per.vector <- apply(data_cd4_us[,cytokines], 2, function(x) quantile(x, 0.995, names = F))
per.vector[c(1,2,3,6,9,10,12,13)] <- c(0.4,0.25,0.4,0.4,0.25,0.25,0.25,0.4)

data_df <- data_cd4
#data_df <- data.frame(merge(data_cd4,md,by='patient'))
data_melt <- melt(data_df, id.vars = c("patient"), 
                  measure.vars = cytokines)

cyto_freq <- ddply(data_melt, .(patient, variable), function(x){
  Freq <- length(x$value[x$value >= as.numeric(per.vector[x$variable])])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})

### Merge md
fm <- dcast(cyto_freq, patient+N ~ variable, value.var = "freq")
#write.csv2(fm, "cytokines_pos_0.4_1.csv")

# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
# melt for ggplot
df_melt <- melt(df_plot, measure.vars = cytokines)

# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm$name <- as.factor(cytokines)
fmm

# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm$name <- factor(fmm$name, levels = fmm$name)
fmm

# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])
df_melt_2 <- merge(df_melt, md , by = 'patient')

# plot boxplots
df_melt_2$Diagnosis <- factor(df_melt_2$Diagnosis, levels = levels(df_melt_2$Diagnosis)[c(2,1,4,5,3)])

b2 <- ggplot(data = df_melt_2, aes(x = Diagnosis, y = value, fill = variable, ymax = 1.5)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'patient')
table(t$Diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig.S8I_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_S8I'

sum <- ddply(df_melt_2, .(variable, Diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum


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


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)


##############################################################################




##Save excel workbook
saveWorkbook(wb, "STATS/Ext_Figure_8.xlsx")

