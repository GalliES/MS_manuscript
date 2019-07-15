######################################
##### load all neccessary things #####
######################################

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
data_all <- merge(data_all, md[,c(2,3,6,9)], by= 'gate_source')





##############################
##### calculating bh-SNE #####
##############################



# select the data
data <- data_all



# clustering channels
clustering_cols <- make.names(as.vector(panel$Antigen[panel$Category != "cytokines"]))
clustering_cols



# create subsample vector
data_naive <- data[data[,"diagnosis"] == 'HC',]
n_sub <- 10000
n <- nrow(data_naive)
set.seed(123)
ix <- sample(1:n, n_sub)

data_dac <- data[data[,"diagnosis"] == 'NINDC',]
n_sub <- 10000
n <- nrow(data_dac)
set.seed(123)
ixx <- sample(1:n, n_sub)

data_hd <- data_all[data_all[,"diagnosis"] == 'RRMS',]
n_sub <- 10000
n <- nrow(data_hd)
set.seed(123)
ixxx <- sample(1:n, n_sub)

t <- rbind(data_dac[ixx,], data_naive[ix,], data_hd[ixxx,])


# prepare data for Rtsne
data_rtsne <- t[ , clustering_cols]

# run bh SNE
set.seed(123)
out_rtsne <- Rtsne.multicore(data_rtsne, dims = 2, perplexity = 50, theta = 0.1, 
                   max_iter = 10000, verbose = T, pca = F, check_duplicates = T)

#saveRDS(out_rtsne, 'rtsne_long.RDS')
out_rtsne <- readRDS('FILES/rtsne_long.RDS')


# prepare the tSNE data
tsne <- as.data.frame(out_rtsne$Y)
colnames(tsne) <- c("tSNE1", "tSNE2")



# plot tSNE black
t1 <- ggplot(tsne, aes(x = tSNE1, y = tSNE2)) +
        geom_point(size = 0.5) +
        coord_fixed(ratio = 1) +
        theme_tsne + 
        theme_facs
t1




# save the plots as png
ggsave(filename = "OUTPUT/tsne_final.png", plot = t1, scale = 1.5)


# prepare the expression data
data_plot <- cbind(t[,c(2:36)], tsne)
#colnames(data_plot) [11] <- "CD7"

data_melt <- melt(data_plot, id.vars = c("tSNE1", "tSNE2"))



# plot tSNEs with expression overlayed
t2 <- ggplot(data_melt, aes(x = tSNE1, y = tSNE2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = jet.colors(100), limits = c(0,1)) +
  facet_wrap(~ variable, ncol = 7, scales = "fixed") +
  theme(plot.background = element_rect(color="white"),
        title = element_text(size = rel(1.1)),
        legend.position = "right",
        legend.key.size = unit(5,"point"),
        legend.key.width = unit(5,"point"),
        plot.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme_tsne
t2

write.xlsx(data_melt, 'temp_ext_fig_1A.xlsx')

# save the plots as png
ggsave(filename = "OUTPUT/tsne_expression_final.png", t2, scale = 1, 
       width = 7*2, height = 5*2.3, units = c("in"))


# prepare the expression data
data_plot <- cbind(t[,c(2:37,46)], tsne) 
data_melt <- melt(data_plot, id.vars = c("tSNE1", "tSNE2",'diagnosis','manual_labels'))

data_dff <- data_plot
data_dff[data_dff[,"manual_labels"] == 1,"manual_labels"] <- "CD8" 
data_dff[data_dff[,"manual_labels"] == 2,"manual_labels"] <- "NKT" 
data_dff[data_dff[,"manual_labels"] == 3,"manual_labels"] <- "NK" 
data_dff[data_dff[,"manual_labels"] == 4,"manual_labels"] <- "CD4" 
data_dff[data_dff[,"manual_labels"] == 5,"manual_labels"] <- "gdT" 
data_dff[data_dff[,"manual_labels"] == 6,"manual_labels"] <- "Bcell" 
data_dff[data_dff[,"manual_labels"] == 7,"manual_labels"] <- "Myeloid" 
data_dff$manual_labels <- factor(data_dff$manual_labels, levels = c("CD8",'NKT','NK','CD4', "gdT", "Bcell",   "Myeloid")) 
data_dff$diagnosis <- factor(data_dff$diagnosis, levels = c('HC','OND','RRMS')) 

temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))
t3 <- ggplot(subset(data_dff, diagnosis != 'no'), aes(x = tSNE1, y = tSNE2, color = as.factor(manual_labels))) +
  geom_point(size = 0.001) +
  coord_fixed(ratio = 1) +
  #scale_x_discrete(limits = c("CD4","CD8", "TCRgd", "B cell", "NK", "Myeloid"))+
  scale_fill_manual(values = db1[c(1,5,4,3,6,7)]) +
  #scale_fill_manual(values = db1) +
  #scale_color_manual(values = db1[c(4,   6,2,5,3,1,4)]) +
  #scale_colour_manual(name = NULL, values = db1) +
  #scale_fill_manual(values = db1[c(5,4,3,1,7,6)]) +
  #scale_colour_gradientn(colours = jet.colors(100), limits = c(0,6)) +
  #facet_wrap(~ diagnosis, ncol = 3, scales = "free") +
  theme_facs+
  #facet_wrap("weeks") +
  theme_tsne +
  theme(plot.background = element_rect(color="white"),
        title = element_text(size = rel(1.1)),
        legend.position = "right",
        legend.key.size = unit(5,"point"),
        legend.key.width = unit(5,"point"),
        plot.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))
t3


write.xlsx(data_dff, 'temp_ext_fig_1B.xlsx')

#dev.off()
t3 <- ggplot(subset(data_dff, diagnosis != 'no'), aes(x = tSNE1, y = tSNE2, color =  as.factor(manual_labels))) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_x_discrete(limits = c("CD4","CD8", "TCRgd", "B cell", "NK", "Myeloid"))+
  scale_fill_manual(values = db1[c(1,5,4,3,6,7)]) +
  #scale_fill_manual(values = db1[c(5,4,3,1,7,6)]) +
  #scale_colour_gradientn(colours = jet.colors(100), limits = c(0,6)) +
  facet_wrap(~ diagnosis, ncol = 3, scales = "free") +
  theme_facs+
  #facet_wrap("weeks") +
  theme_tsne +
  theme(plot.background = element_rect(color="white"),
        title = element_text(size = rel(1.1)),
        legend.position = "right",
        legend.key.size = unit(5,"point"),
        legend.key.width = unit(5,"point"),
        plot.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))
t3



# save the plots as png
ggsave(filename = "clusters_tsne_final.png", plot = t3, scale = 1.5, 
       width = 3*2, height = 1*2.3, units = c("in"))
ggsave(filename = "OUTPUT/clusters_tsne_final.pdf", plot = t3, scale = 1.5)



library(openxlsx)
wb = createWorkbook()

addWorksheet(wb, 'Ext_Fig_1A')
addWorksheet(wb, 'Ext_Fig_1B')
writeData(wb, sheet = 'Ext_Fig_1A', x = data_melt)
writeData(wb, sheet = 'Ext_Fig_1B', x = data_dff)

##Save excel workbook
saveWorkbook(wb, "Ext_Figure_1.xlsx")
getwd()
