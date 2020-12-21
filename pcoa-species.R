library(vegan)
library(ggplot2)
library(plyr)
spe1#species data
spe1data.t<-t(spe1data)
group <- read.delim('group.txt',  sep = '\t', stringsAsFactors = F)
spe1data.t<-t(otuFW[,-1])
spe1data.t<-t(spe1[,-1])
#distance <- vegdist(spe1data.t, method = 'euclidean')
distance <- vegdist(spe1data.t, method = 'bray')
distance <- vegdist(spe1data.t, method = 'canberra')
#distance <- vegdist(spe1data.t, method = 'jaccard')
pcoa <- cmdscale(distance, k = 5, eig = T)
point <- data.frame(pcoa$point)
#point.plsda<-data.frame(ma.plsda$variates$X)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

#pcoa <- cmdscale(as.dist(dis), k = (nrow(dis) - 1), eig = T)
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$name <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

sample_site <- merge(sample_site, group, by = 'name', all.x = T)

sample_site$group <-c(rep("Control",5),rep("CR",5),rep("HF",5))
group_border <- ddply(sample_site, 'group', function(df) df[chull(df[[2]], df[[3]]), ])

pcoa_plot <- 
  
  ggplot(sample_site, aes(PCoA1, PCoA2, group = group))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  #geom_polygon(data = group_border, aes(fill = group)) +
  
  geom_point( size = 1.5) + 
  stat_ellipse(aes(x = PCoA1, y = PCoA2, fill = group), geom = "polygon", alpha = 1/2,levles = 0.95) +
  scale_fill_manual(values = c('blue2',  'yellow2','red')) +
  #guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) +
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%')) 

ggsave('PCoA-b.png', pcoa_plot)

###plsda
ma.plsda<-mixOmics::plsda(t(spe1data),group[,2],mode = "classic",ncomp = 3)
point.plsda<-data.frame(ma.plsda$variates$X)
sample_site.plsda <- data.frame(point.plsda)[1:2]
sample_site.plsda$name <- rownames(data.ma)
names(sample_site.plsda)[1:2] <- c('PC1', 'PC2')

sample_site.plsda <- merge(sample_site.plsda, group, by = 'name', all.x = T)

sample_site.plsda$group <- factor(sample_site.plsda$group, levels = c('Control', 'CR','HF'))
group_border <- ddply(sample_site.plsda, 'group', function(df) df[chull(df[[2]], df[[3]]), ])

plsda_plot <- ggplot(sample_site.plsda, aes(PC1, PC2, group = group))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  #geom_polygon(data = group_border, aes(fill = group)) +
  geom_point( size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = c('blue2', 'yellow2', 'red')) +
  stat_ellipse(aes(x = PC1, y = PC2, fill = group), geom = "polygon", alpha = 1/2, levles = 0.95) +
  labs(x = paste('PC1 axis1: ', round(100 * ma.plsda$explained_variance$X[1], 2), '%'), y = paste('PC2 axis2: ', round(100 * ma.plsda$explained_variance$X[2], 2), '%')) 
ggsave('PLSDA-3group.png', plsda_plot, width = 6, height = 5)





##HF-C
spe1data.thc<-t(spe1data[,c(1:10)]) 
distance.hc <- vegdist(spe1data.thc, method = 'bray')
distance.hc <- vegdist(spe1data.thc, method = 'canberra')
pcoa.hc <- cmdscale(distance.hc, k = 5, eig = T)
point.hc <- data.frame(pcoa.hc$point)
pcoa_eig.hc <- (pcoa.hc$eig)[1:2] / sum(pcoa.hc$eig)

#pcoa <- cmdscale(as.dist(dis), k = (nrow(dis) - 1), eig = T)
sample_site.hc <- data.frame({pcoa.hc$point})[1:2]
sample_site.hc$name <- rownames(sample_site.hc)
names(sample_site.hc)[1:2] <- c('PCoA1', 'PCoA2')
#gro<-c(rep("blue2",5),rep("green3",5),rep("red",5))

sample_site.hc <- merge(sample_site.hc, group[c(1:10),], by = 'name', all.x = T)

sample_site.hc$group <- factor(sample_site.hc$group, levels = c('Control','HF'))
group_border.hc <- ddply(sample_site.hc, 'group', function(df) df[chull(df[[2]], df[[3]]), ])

pcoa_plot <- ggplot(sample_site.hc, aes(PCoA1, PCoA2, group = group))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
#  geom_polygon(data = group_border.hc, aes(fill = group)) +
  geom_point( size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = c('blue2',  'red')) +
  stat_ellipse(aes(x = PCoA1, y = PCoA2, fill = group), geom = "polygon", alpha = 1/2, levles = 0.95) +
  
 # guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) +
  
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%')) 
ggsave('PCoA-HFvsC-bray.png', pcoa_plot, width = 6, height = 5)


##CR-C
spe1data.thc<-t(spe1data[,c(6:15)]) 
distance.hc <- vegdist(spe1data.thc, method = 'bray')
distance.hc <- vegdist(spe1data.thc, method = 'canberra')
pcoa.hc <- cmdscale(distance.hc, k = 5, eig = T)
point.hc <- data.frame(pcoa.hc$point)
pcoa_eig.hc <- (pcoa.hc$eig)[1:2] / sum(pcoa.hc$eig)

#pcoa <- cmdscale(as.dist(dis), k = (nrow(dis) - 1), eig = T)
sample_site.hc <- data.frame({pcoa.hc$point})[1:2]
sample_site.hc$name <- rownames(sample_site.hc)
names(sample_site.hc)[1:2] <- c('PCoA1', 'PCoA2')
#gro<-c(rep("blue2",5),rep("green3",5),rep("red",5))

sample_site.hc <- merge(sample_site.hc, group[c(6:15),], by = 'name', all.x = T)

sample_site.hc$group <- factor(sample_site.hc$group, levels = c('Control','CR'))
group_border.hc <- ddply(sample_site.hc, 'group', function(df) df[chull(df[[2]], df[[3]]), ])

pcoa_plot <- ggplot(sample_site.hc, aes(PCoA1, PCoA2, group = group))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
 # geom_polygon(data = group_border.hc, aes(fill = group)) +
  geom_point( size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = c('blue2',  'green3')) +
  stat_ellipse(aes(x = PCoA1, y = PCoA2, fill = group), geom = "polygon", alpha = 1/2, levles = 0.95) +
  
  #guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) +
  
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%')) 
ggsave('PCoA-CRvsC-canberra.png', pcoa_plot, width = 6, height = 5)


##CR-H
spe1data.thc<-t(spe1data[,-c(6:10)]) 
distance.hc <- vegdist(spe1data.thc, method = 'bray')
distance.hc <- vegdist(spe1data.thc, method = 'canberra')
pcoa.hc <- cmdscale(distance.hc, k = 5, eig = T)
point.hc <- data.frame(pcoa.hc$point)
pcoa_eig.hc <- (pcoa.hc$eig)[1:2] / sum(pcoa.hc$eig)

#pcoa <- cmdscale(as.dist(dis), k = (nrow(dis) - 1), eig = T)
sample_site.hc <- data.frame({pcoa.hc$point})[1:2]
sample_site.hc$name <- rownames(sample_site.hc)
names(sample_site.hc)[1:2] <- c('PCoA1', 'PCoA2')
#gro<-c(rep("blue2",5),rep("green3",5),rep("red",5))

sample_site.hc <- merge(sample_site.hc, group[-c(6:10),], by = 'name', all.x = T)

sample_site.hc$group <- factor(sample_site.hc$group, levels = c('HF','CR'))
group_border.hc <- ddply(sample_site.hc, 'group', function(df) df[chull(df[[2]], df[[3]]), ])

pcoa_plot <- ggplot(sample_site.hc, aes(PCoA1, PCoA2, group = group))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
 # geom_polygon(data = group_border.hc, aes(fill = group)) +
  geom_point( size = 1.5, alpha = 0.8) + 
  stat_ellipse(aes(x = PCoA1, y = PCoA2, fill = group), geom = "polygon", alpha = 1/2, levles = 0.95) +
  scale_fill_manual(values = c('red',  'green3')) +
#  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) +
  
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%')) 
ggsave('PCoA-CRvsHF-canberra.png', pcoa_plot, width = 6, height = 5)





