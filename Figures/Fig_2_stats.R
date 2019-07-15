
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
data_all <- merge(data_all, md[,c(2,7)], by= 'gate_source')


data_df <- droplevels(subset(data_all, diagnosis == 'RRMS' | diagnosis == 'NINDC'))
data_df <- subset(data_df, GM.CSF >= 0.25) 


# calculate frequencies
fm <- ddply(data_df, .(gate_source), function(x) {
  Freq <- rep(0, 7)
  for(i in 1:7) {
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
colnames(fm) <- c("gate_source", c("CD8", "NKT", "NK", "CD4", "gdT.cells",
                                   "B.cell", "Myeloid"))
head(fm)



# exclude samples with less than 100 cells
fm$N <- ddply(data_df, .(gate_source), nrow)$V1
fm[fm$N <= 50, 2:8] <- NA
fm



# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:8])
head(df_melt)


# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
cluster_names <- make.names(c("CD8 T cells", "NKT cells", "NK cells", "CD4 T cells", "gdT cells",
                              "B cells", "Monocytes"))
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

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable, ymax = 100)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 7) +
  scale_fill_manual(values = db_lin) +
  scale_color_manual(values = db_lin) +
  #ylim(0,70) +
  theme_bar2
b2  



nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig_2B_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_2B'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2B_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2B_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2B_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################



data_df <- droplevels(subset(data_all, diagnosis == 'RRMS' | diagnosis == 'NINDC'))

fm <- ddply(data_df, .(gate_source, manual_labels), function(x) {
  cell_n <- nrow(x)
  gm <- nrow(x[x[,"GM.CSF"] >= 0.25,])
  Freq = gm/cell_n*100
})



# exclude samples with less than 50 cells
N <- ddply(data_df, .(gate_source), nrow)$V1
fm[N < 50, 3] <- NA
fm



# cast and rename
fm <- dcast(fm, gate_source ~ manual_labels, value.var = "V1") 
colnames(fm)[-1] <- c("CD8", "NKT", "NK", "CD4", "gdT.cells",
                      "B.cell", "Myeloid")
fm


# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:8])
head(df_melt)


# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
cluster_names <- make.names(c("CD8 T cells", "NKT cells", "NK cells", "CD4 T cells", "gdT cells",
                              "B cells", "Monocytes"))
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

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db_lin) +
  scale_color_manual(values = db_lin) +
  #ylim(0,70) +
  theme_bar2
b2  

nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=31
##RRMS=31


# save
ggsave(filename = "OUTPUT/Fig_2C_gm_per_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary
fig <- 'Fig_2C'
sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2C_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2C_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2C_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################




# select cd4 only
data_cd4 <- data_df[data_df[,"manual_labels"] == 4,]
table(data_cd4[,"cd4_labels"])

# select HD and GM+
#data_cd4_hd <- data_cd4[data_cd4[,"gate_source"] %in% hd_id,]
dim(data_cd4_gm)


# calculate frequency of positive cells for all donors and all clusters
fm <- ddply(data_cd4, .(gate_source), function(x) {
  cell_n <- nrow(x)
  gm <- as.numeric(table(x$cd4_labels))
  Freq = gm/cell_n*100
  #if(cell_n < 50) {Freq <- NA} else {Freq <- Freq}
})
colnames(fm) <- c("gate_source", "Tnaive", "Teff", "Tem", "Tcm")
fm



# exclude samples with less than 50 cells
N <- ddply(data_cd4, .(gate_source), nrow)$V1
fm[N < 50, 2:ncol(fm)] <- NA
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



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
#df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[])

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable, ymax = 100)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db_mem) +
  scale_color_manual(values = db_mem) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=31
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig_2D_gm_per_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_2D'
sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2D_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2D_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2D_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################



data_cd4_gm <- data_cd4[data_cd4[,"GM.CSF"] >= 0.25,] # For CD4 T cells this is appropriate
dim(data_cd4_gm)
data_cd4_gm$cd4_labels <- as.factor(data_cd4_gm$cd4_labels)
# calculate frequency of positive cells for all donors and all clusters
fm <- ddply(data.frame(data_cd4_gm), .(gate_source), function(x) {
  cell_n <- nrow(x)
  gm <- as.numeric(table(x$cd4_labels))
  Freq = gm/cell_n*100
  #if(cell_n < 50) {Freq <- NA} else {Freq <- Freq}
})
colnames(fm) <- c("gate_source", "Tnaive", "Teff", "Tem", "Tcm")
fm



# exclude samples with less than 50 cells
N <- ddply(data_cd4, .(gate_source), nrow)$V1
fm[N < 1000, 2:ncol(fm)] <- NA
fm
N


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



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
#df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[])

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable, ymin = 0, ymax = 100)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 7) +
  scale_fill_manual(values = db_mem) +
  scale_color_manual(values = db_mem) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=29
##RRMS=30

# save
ggsave(filename = "OUTPUT/Fig_2E_gm_per_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_2E'
sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2E_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2E_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2E_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################



# select data & calculate frequencies
data_cd4_gm <- data_cd4[data_cd4[,"GM.CSF"] >= 0.25,] # For CD4 T cells this is appropriate
dim(data_cd4_gm)
data_cd4_gm$cd4_labels <- as.factor(data_cd4_gm$cd4_labels)

data_df <- data_cd4_gm
data_melt <- melt(data_df, id.vars = "gate_source", measure.vars = cytokines[-7])
cyto_freq <- ddply(data_melt, .(gate_source, variable), function(x){
  Freq <- length(x$value[x$value >= 0.325])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})


colnames(cyto_freq)[3] <- 'value'
df_melt <- cyto_freq
# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
#cluster_names <- make.names(c("CD8 T cells", "NKT cells", "NK cells", "CD4 T cells", "gdT cells",
#                              "B cells", "Monocytes"))
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

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  geom_point(size = 4, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 12) +
  scale_fill_manual(values = db_cy) +
  scale_color_manual(values = db_cy) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=29
##RRMS=30

# save
ggsave(filename = "OUTPUT/Fig_2G_pop_boxplots.pdf", plot = b2, width = 1.5*6, height = 2*3, 
       scale = 2.5, useDingbats = F)


########################
###### STAT block ######
### Summary

fig <- 'Fig_2G'
sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2G_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2G_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2G_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################



# select the data and the cytokine
data <- data.frame(data_cd4)
(labels <- cytokines[cytokines != "GM.CSF"])
cutoff <- 0.25 #for GM-CSf
min_cells <- 50
i <- 1


# select temp data 
(i <- i+1) #from second time onwards
labels[i]
data_temp <- data[data[,labels[i]] >= cutoff,]
dim(data_temp)



# loop through all the gate sources
pos <- ddply(data_temp, .(gate_source), function(x) {
  c(freq = as.numeric(prop.table(table(x$GM.CSF >= cutoff))["TRUE"]*100),
    N = nrow(x))
})
pos



# for few cells set NA
pos[is.na(pos$freq), "freq"] <- 0
pos[pos$N < min_cells, "freq"] <- NA
pos



# combine with the rest (first time)
#gm_pos <- pos
#gm_pos <- gm_pos[,-3]
#colnames(gm_pos)[2] <- labels[i]



# combine with the rest (from second on)
gm_pos <- merge(gm_pos, pos[,1:2], by = "gate_source")
colnames(gm_pos)[i+1] <- labels[i]
head(gm_pos)

#########
# save or load
#saveRDS(gm_pos, file = "FILES/gm_pos_total_th.rds")
fm <- readRDS(file = "FILES/gm_pos_total_th.rds")


# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)



# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:12])
head(df_melt)


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
#df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[])

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable, ymin = 0, ymax = 100)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 11) +
  scale_fill_manual(values = db_cy) +
  scale_color_manual(values = db_cy) +
  #ylim(0,70) +
  theme_bar2
b2  

nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=13
##RRMS=20

# save
ggsave(filename = "OUTPUT/Fig_2H_gm_per_cyt_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_2H'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2H_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2H_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2H_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################

####adjust th_labels with reduced data
data_all <- exprs(ff)
data <- data_all[data_all[,'manual_labels']==4, ]
dim(data)
th_labels <- readRDS(file = "FILES/th_labels.rds")
length(th_labels)
# select data
data_cd4 <- cbind(data, th_labels)
data_cd4 <- merge(data_cd4, md , by = 'gate_source')
#th_labels <- data_cd4$th_labels
#saveRDS(th_labels, 'FILES/th_labels_post_md_merge.rds')

#### select data and merge with th_labels 
### @ careful with merging and in case run the code above
#th_labels <- readRDS(file = "FILES/th_labels_post_md_merge.rds")
#data_cd4 <- droplevels(subset(data_all, manual_labels == 4 ))
#data_cd4 <- cbind(data_cd4, th_labels)
data_cd4 <- droplevels(subset(data_cd4, diagnosis == 'NINDC' | diagnosis == 'RRMS'))


data <- data_cd4[data_cd4[, 'GM.CSF'] >= 0.325,]
data_df <- data.frame(data)
data_df <- droplevels(subset(data_df, th_labels != 1))
dim(data_df)
head(data_df)
# calculate frequencies

# calculate frequencies for all donors
t <- table(data_df[,"gate_source"], data_df[,"th_labels"])
fm <- prop.table(t, 1)*100
apply(fm, 2, median)



# cd4 exclude stuff
t2 <- table(data[,"gate_source"])
excl <- names(which(t2 < 50))
fm[rownames(fm) %in% excl,] <- NA
fm <- data.frame(fm)
fm


# cast and rename
fm <- dcast(fm, Var1 ~ Var2, value.var = "Freq") 
colnames(fm) <-  c("gate_source",  "Th1", "Th2", "Th17", "Th22", "ThGM", "Tfh","Treg", "Other")
fm


# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
head(df_plot)


# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:9])
head(df_melt)


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

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable, ymin = 0, ymax = 80)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  facet_wrap('variable', scales = 'fix', ncol = 11) +
  scale_fill_manual(values = db_th) +
  scale_color_manual(values = db_th) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=22
##RRMS=29

# save
ggsave(filename = "OUTPUT/Fig_2J_gm_per_cyt_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

fig <- 'Fig_2J'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_2J_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2J_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_2J_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################


##Save excel workbook
saveWorkbook(wb, "STATS/Figure_2.xlsx")

