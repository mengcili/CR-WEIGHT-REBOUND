


changBAcon1<-changBA[c(1:18),c(2:49)]
changBAper1<-changBA[c(1:18),c(50:97)]
changBAcon2<-changBA[c(19:36),c(2:49)]
changBAper2<-changBA[c(19:36),c(50:97)]
changBAcon3<-changBA[c(37:54),c(2:49)]
changBAper3<-changBA[c(37:54),c(50:97)]
changBAcon4<-changBA[c(55:72),c(2:49)]
changBAper4<-changBA[c(55:72),c(50:97)]
changBAcon5<-changBA[c(73:90),c(2:49)]
changBAper5<-changBA[c(73:90),c(50:97)]



datBAtest<-changBAper3
datBAtest<-changBAcon3
cecplot1<-cecplot[,-4]
datBAtest<-conBACec[,-1]
data.ma2<-scale(datBAtest,center = F, scale = TRUE)
Y<-c(rep(c("KPD","vehicle","PD"),c(5,5,5)))
rownames(data.ma2)<-c("H1","H2","H3","H4","H5","H6","C1","C2","C3","C4","C5","C6","CR1","CR2","CR3","CR4","CR5","CR6")

y.group<-kjl_group[,2]
data.ma2<-scale(dkjl[,-1],center = F, scale = TRUE)
ma.plsda<-mixOmics::plsda(data.ma2,y.group,mode = "classic",ncomp = 3)
ma.pca<-mixOmics::pca(t(otuFW[,-1]),ncomp = 4)
vips4<-vip(ma.plsda)
vipd1<-vips4[,1]
sort(vipd1,decreasing = T)

ma.opls<-opls(data.ma2,y.group)

opls(data.ma2,y.group,fig.pdfC="plsda")


plotIndiv(ma.pca, comp = c(1,2),group=group[,2],ind.names = F,legend=TRUE,
          ellipse =TRUE,point.lwd = 1,style = "ggplot2",col=c("red","pink","blue","green"),pch=16,title = "PCA",legend.title = "Group")
plotIndiv(ma.plsdal2, comp = c(1,2),ind.names = F,legend=TRUE,ellipse =TRUE,
          point.lwd = 1,style = "ggplot2",col=c("red","green","blue"),pch=16,title = "PLS-DA",legend.title = "Group")

#point.plsda<-data.frame(ma.plsda$variates$X)
point.opls<-data.frame(ma.opls@scoreMN)
sample_site.plsda <- data.frame(point.plsda)[1:2]
sample_site.plsda$name <- rownames(data.ma2)
names(sample_site.plsda)[1:2] <- c('PC1', 'PC2')
group2<-rbind(group,c("CR6","CR"),c("H6","HF"),c("C6","Control"))
sample_site.plsda <- merge(sample_site.plsda, group2, by = 'name', all.x = T)

sample_site.plsda$group <- factor(sample_site.plsda$group, levels = c('Control', 'CR','HF'))
group_border <- ddply(sample_site.plsda, 'group', function(df) df[chull(df[[2]], df[[3]]), ])

plsda_plot <- ggplot(sample_site.plsda, aes(PC1, PC2, group = group))+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  
  geom_point(col=c("blue2")) +
  geom_point( size = 1.5, alpha = 0.8) + 
  scale_fill_manual(values = c('blue2', 'green3', 'red')) +
  stat_ellipse(aes(x = PC1, y = PC2, fill = group), geom = "polygon", alpha = 1/2,level=0.95) +
  labs(x = paste('PC1 axis1: ', round(100 * ma.plsda$explained_variance$X[1], 2), '%'), y = paste('PC2 axis2: ', round(100 * ma.plsda$explained_variance$X[2], 2), '%')) 



ggsave('oPLSDA-con5.png', plsda_plot, width = 6, height = 5)









