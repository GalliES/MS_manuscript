
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



# read in panel
panel <- read.xls("FILES/panel_discovery.xlsx", sheet = 1, header = TRUE, verbose = FALSE)

(cytokines <- make.names(as.vector(panel$Antigen[ panel$Category == "cytokines"])))
(cytokines <- cytokines[-6])

# read in metadata
md <- read.xls("FILES/meta_data_discovery_cohort.xlsx", sheet = 1, header = T, verbose = F)
colnames(md)
head(data_all)
data_all <- merge(data_all, md[,c(2,7,9)], by= 'gate_source')

data_df <- droplevels(subset(data_all, diagnosis == 'RRMS' | diagnosis == 'NINDC'))

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
fm[fm$N <= 199, 2:8] <- NA
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
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[c(4,1,6,3,5,2,7)])
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

dim(table(df_melt_2$gate_source))



nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=31
##RRMS=31


# save
ggsave(filename = "OUTPUT/Fig_1C_pop_boxplots.pdf", plot = b2, width = 1.5*7, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary
fig <- 'Fig_1C'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_1C_pop_boxplots_sum.xlsx')

#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$diagnosis, p.adjust.method = "none", paired = F, correct = T, exact = F)
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

write.xlsx(wres, 'STATS/Fig_1C_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_1C_pop_boxplots_effect_size.xlsx')


addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)


##############################################################################



# select data & calculate frequencies
data_df <- droplevels(subset(data_all, diagnosis == 'RRMS' | diagnosis == 'NINDC'))
data_melt <- melt(data_df, id.vars = "gate_source", measure.vars = cytokines)
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
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[c(10,12,8,7,11,5,1,3,4,6,9,2)])

df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])

b2 <- ggplot(data = df_melt_2, aes(x = diagnosis, y = value, fill = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  geom_point(size = 4, aes(colour = variable)) +
  facet_wrap('variable', scales = 'free', ncol = 12) +
  scale_fill_manual(values = db_cy) +
  scale_color_manual(values = db_cy) +
  #ylim(0,70) +
  theme_bar2
b2  


t <- subset(df_melt_2, variable == 'GM.CSF')
table(t$diagnosis)
##NINDC=31
##RRMS=31

# save
ggsave(filename = "OUTPUT/Fig_1E_pop_boxplots.pdf", plot = b2, width = 1.5*6, height = 2*3, 
       scale = 2.5, useDingbats = F)


########################
###### STAT block ######
### Summary
fig <- 'Fig_1E'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_1E_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_1E_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_1E_pop_boxplots_effect_size.xlsx')


addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)



##############################################################################


####Fig.1F
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

dim(table(df_melt_2$gate_source))


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=31
##RRMS=31


# save
ggsave(filename = "OUTPUT/Fig_1F_pop_boxplots.pdf", plot = b2, width = 1.5, height = 3, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary
fig <- 'Fig_1F'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_1F_pop_boxplots_sum.xlsx')


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

write.xlsx(wres, 'STATS/Fig_1F_pop_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_1F_pop_boxplots_effect_size.xlsx')



addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################



# select data
data_df <- data_all
data_df <- droplevels(subset(data_all, diagnosis == 'RRMS' | diagnosis == 'NINDC'))
data_df_s <- subset(data_df, intersection_3_runs == 1)


# calculate frequencies
fm <- ddply(data_df_s, .(gate_source), function(x) {
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
fm$N <- ddply(data_df_s, .(gate_source), nrow)$V1
fm[fm$N <= 50, 2:8] <- NA
fm

nrow(fm[complete.cases(fm),])



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

b2 <- ggplot(data = df_melt_2, aes(x = variable, y = value)) +
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  #facet_wrap('variable', scales = 'fix', ncol = 7) +
  scale_fill_manual(values = db_lin) +
  scale_color_manual(values = db_lin) +
  #ylim(0,70) +
  theme_bar2
b2  


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=9
##RRMS=19


# save
ggsave(filename = "OUTPUT/Fig_1H_sig_pop_boxplots.pdf", plot = b2, width = 2*3.2, height = 2*2.1, 
       scale = 0.75, useDingbats = F)


########################
###### STAT block ######
### Summary
fig <- 'Fig_1H'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_1H_sig_pop_boxplots_sum.xlsx')


addWorksheet(wb, paste(fig,'_sum',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)


##############################################################################


# select data
data_df <- droplevels(subset(data_all, diagnosis == 'RRMS' | diagnosis == 'NINDC'))
data_df_s <- subset(data_df, intersection_3_runs == 1  & manual_labels == 4  )

# select data & calculate frequencies
data_melt <- melt(data_df_s, id.vars = "gate_source", measure.vars = cytokines)
cyto_freq <- ddply(data_melt, .(gate_source, variable), function(x){
  Freq <- length(x$value[x$value >= 0.325])/length(x$value)*100
  if(length(Freq) == 0) Freq <- 0
  c(freq = Freq,
    N = length(x$value))
})


colnames(cyto_freq)[3] <- 'value'


# exclude samples with less than 50 cells
df_melt <- droplevels(subset(cyto_freq, N >= 50))
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

b2 <- ggplot(data = df_melt_2, aes(x = variable, y = value, fill = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  geom_point(size = 4, aes(colour = variable)) +
  #facet_wrap('variable', scales = 'free', ncol = 12) +
  scale_fill_manual(values = db_cy) +
  scale_color_manual(values = db_cy) +
  #ylim(0,70) +
  theme_bar2
b2  

nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=9
##RRMS=19

# save
ggsave(filename = "OUTPUT/Fig_1I_sign_cyto_boxplots.pdf", plot = b2, width = 1.5*6, height = 2*3, 
       scale = 2.5, useDingbats = F)


########################
###### STAT block ######
### Summary
fig <- 'Fig_1I'

sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum
write.xlsx(sum, 'STATS/Fig_1I_sign_cyto_boxplots_sum.xlsx')


addWorksheet(wb, paste(fig,'_sum',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)


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

write.xlsx(wres, 'STATS/Fig_1I_sign_cyto_boxplots_stat.xlsx')


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

write.xlsx(wres, 'STATS/Fig_1I_sign_cyto_boxplots_effect_size.xlsx')

##############################################################################


##Save excel workbook
saveWorkbook(wb, "STATS/Figure_1.xlsx")



