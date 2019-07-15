##### ##### ##### ##### ##### 
##### Stats dmf cohort ##### 
##### ##### ##### ##### ##### 
library(openxlsx)
wb = createWorkbook()

fm <- read.xlsx('FILES/dmf_cohort.xlsx')
str(fm)

# show as bargraph plot
df_plot <- data.frame(fm)
df_plot$samples <- 1:nrow(fm)
colnames(df_plot)[23] <- 'patient'
df_plot <- droplevels(subset(df_plot, sample_id != 'US'))
head(df_plot)
tail(df_plot)


# melt for ggplot
data_melt <- melt(df_plot, measure.vars = colnames(df_plot)[2:19])

# calculate summary frequencies
fmm <- ddply(data_melt, .(variable), function(x) {c(median = median(x$value, na.rm = T),
                                                    sem = b.median(x$value, 1000),
                                                    n = as.numeric(table(is.na(x$value))["FALSE"]))}) 
fmm
order <- as.numeric(rownames(fmm[order(fmm$median, decreasing = T),]))

# reoder for plotting
data_melt$variable <- factor(data_melt$variable, levels = levels(data_melt$variable)[order])
data_melt_drop <- data_melt
data_melt_drop$treated <- as.factor(data_melt_drop$treated)
data_melt_drop$treated <- factor(data_melt_drop$treated, levels = levels(data_melt_drop$treated)[c(2,1)])

wres_all <- ddply(data_melt_drop, .(variable), function(x) {
  before <- subset(x, treated == 'ut')
  before <- before[c(as.vector(order(before$patient, decreasing = F))),]
  after <-subset(x, treated == 't')
  after <- after[c(as.vector(order(after$patient, decreasing = F))),]
  res <- wilcox.test(before[,'value'], after[,'value'], paired = TRUE)
  res$p.value
  with(res, data.frame(method, t(as.vector(p.value))))
})
wres_all


# make an overview plot
comb <- ggplot(droplevels(subset(data_melt_drop,  patient != 21)), aes(x = treated, y = value, ymin = 0, ymax = value*1.2)) +
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
#ggsave(filename = "OUTPUTnat/temp_cd4_cytokinepositivity_all_boxplot.pdf", plot = comb, useDingbats = F, scale = 1.75, width = 7*1.65, height = 2*3.5, units = c("in"))



# rename and round
names(wres_all)[3]  <- c("p.value")
wres_all[,3] <- round(wres_all[,3], 4)
wres_all


### Summary
sum <- ddply(data_melt_drop, .(variable, treated) , function(x) {
  c(median = median(x$value, na.rm = T), 
    sem = b.median(x$value, 1000), 
    sd = sd(x$value, na.rm = T), 
    n = length(na.omit(x$value)),
    max = max(x$value, na.rm = T),
    min = min(x$value, na.rm = T)
  )})
sum

write.csv(sum, 'STATS/Fig_5_sign_cyto_boxplots_sum.xlsx')


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


write.xlsx(wres_all, 'STATS/Fig_5_stat.xlsx')


# make p-value data.frame
wp <- wres_all
wp$p.value <- paste("P = ", round(wres_all$BH.p.value, 4), t,  sep = "")
wp$x <- 1.5  
wp$value <- ddply(data_melt_drop, .(variable), function(x) {max(x$value, na.rm = T)})$V1
wp



#####change here for different graphs
colnames(fm)
temp <- droplevels(subset(data_melt_drop, variable %in% colnames(fm)[2:7]))
temp <- droplevels(subset(data_melt_drop, variable %in% colnames(fm)[8:9]))
temp$variable <- factor(temp$variable, levels = levels(temp$variable)[c(2,1)])

temp <- droplevels(subset(data_melt_drop, variable %in% colnames(fm)[c(10,15)]))
temp <- droplevels(subset(data_melt_drop, variable %in% colnames(fm)[2:7]))

# plot boxplots  geom_line(linetype = "dashed", color = "black", size= 1)+
b5 <- ggplot(temp, aes(x = treated, y = value, ymin = 0, ymax = value*1.2)) +
  geom_line(linetype = "solid",  size= 0.5, aes(group = patient_num))+
  geom_point(size = 4, aes(color = variable)) +
  geom_boxplot(aes(fill = variable), lwd = 0.5, color = "white", outlier.shape = NA) +
  stat_boxplot(geom = 'errorbar', aes(color = variable)) +
  facet_wrap(~ variable, scales = "fix", nrow= 1) +
  geom_text(data = droplevels(subset(wp, variable %in% make.names(levels(temp$variable)))), aes(y = value*1.3, x = x, label = p.value, group = NULL), size = 3) + 
  #geom_blank(data=data_dummy) +
  #geom_text(data = wp, aes(y = value*1.1, x = x, label = p.value, group = NULL), size = 3) + 
  #geom_text(data = wp, aes(y = value*1.2, x = x, label = p.value.nindc, group = NULL), size = 3) +  
  #geom_text(data = wp, aes(y = value*1.1, x = x, label = p.value.hc, group = NULL), size = 3) +      
  theme_bar2 + 
  theme(strip.text = element_text(size = rel(1.5))) +
  theme(panel.spacing = unit(2.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, color = "black"))

b5  

# save
ggsave(filename = "signature_cytokines_boxplot_exclude_responder.pdf", plot = b5, 
       useDingbats = F, scale = 1.75, width = 2*1.65, height = 1*3.5, units = c("in"))

str(data_melt_drop)
effect_size <- ddply( data_melt_drop , .(variable) , function(x) {
  t2 <- wilcox_test(value ~ treated, data = x, distribution = "exact", conf.int = T)
  p <- pvalue(t2)
  z <- t2@statistic@teststatistic
  r <- as.numeric(z)/sqrt(nrow(x))
  d <- 2*r/sqrt(1-r^2)
  y <- x[!is.na(x$value),]
  cles <- cles.fnc(variable = "value", group = "treated", baseline = "ut", 
                   data = y, print = F)
  data.frame(p.value = p, z = -z, r = -r, r2 = r^2, d = -d, cles)
})
effect_size

write.xlsx(wres, 'STATS/Fig_5_effect_size.xlsx')


fig <- 'Fig_5'
addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres_all)
addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =effect_size)


##############################################################################



##Save excel workbook
saveWorkbook(wb, "STATS/Figure_5.xlsx")
