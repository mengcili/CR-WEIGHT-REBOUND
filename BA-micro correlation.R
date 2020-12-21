cormic

corBA

cormic1<-t(cor2data[c(1:9),-1])
corBA1<-t(cor2data[c(12:38),-1])
resu.r<-matrix(NA,nrow=25,ncol=9)
resu.p<-matrix(NA,nrow=25,ncol=9)
resu.rmic<-matrix(NA,nrow=25,ncol=9)
resu.pmic<-matrix(NA,nrow=25,ncol=9)
rownames(resu.p)<-rownames(resu.r)<-colnames(corBA)
colnames(resu.p)<-colnames(resu.r)<-colnames(cormic)
rownames(resu.p)<-rownames(resu.r)<-c("CA"   ,     "DCA"    ,   "X3ketoCA" , "X7ketoDCA" ,"X23norDCA", "CDCA"  ,   
"aMCA"    ,  "bMCA"  ,    "LCA"   ,    "UDCA"     , "wMCA",      "bUDCA"   , 
 "GCA"   ,    "GDCA"   ,   "TCA"   ,    "TDCA"   ,   "TCDCA"  ,   "TaMCA" ,   
 "TbMCA"  ,   "TLCA"    ,  "TUDCA"  ,   "THCA"  ,    "THDCA"   ,  "TDHCA",    
 "TwMCA"  )
for(i in 1:9){
  for(j in 1:25){
   c1<-as.data.frame(cormic[,i])
   c2<-as.data.frame(corBA[,j])
   corre<-cor.test(as.numeric(as.matrix(c1)),as.numeric(as.matrix(c2)),alternative = "two.sided","pearson")
   resu.r[j,i]<-corre$estimate
   resu.p[j,i]<-corre$p.value
   
  cc2<- as.numeric(as.matrix(c2))
  cc1<- as.numeric(as.matrix(c1))
   micxy<-minerva::mine(cc2,cc1)
   
   
   micr<-micxy$MIC
   micp<-0
   
   for(ii in seq_len(101))
   {
      #bootstrap
      options(warn = -1)
      bootx<-matrix(sample(as.numeric(as.matrix(c1)),replace = TRUE))
      booty<-matrix(sample(as.numeric(as.matrix(c2)),replace = TRUE))
      if(sd(bootx)==0){
         bootx<-matrix(sample(as.numeric(as.matrix(c1)),replace = FALSE))
      }
      if(sd(booty)==0){
         booty<-matrix(sample(as.numeric(as.matrix(c2)),replace = FALSE))
      }
      tmp<-minerva::mine(booty,bootx)
      MICtp<-tmp$MIC
      tempM<-ifelse(micr<=MICtp,1,0)
      micp<-tempM+micp
   }
   micp<-micp/101
   
   resu.rmic[j,i]<-micr
   resu.pmic[j,i]<-micp
   
   
   
   
   
   
   }
}



rmic<-cbind(resu.rmic,resu.pmic)

resu.r
resu.padj<-p.adjust(resu.p,"bonferroni")
resu.padjmic<-p.adjust(resu.pmic,"fdr")
dim(resu.padj)<-dim(resu.p)
matp<-matrix(ifelse(resu.padj<0.05,"*",""))
matp<-matrix(ifelse(resu.p<0.05,"*",""))
dim(resu.padjmic)<-dim(matp)<-dim(resu.p)
rownames(matp)<-rownames(resu.padj)<-rownames(resu.r)
colnames(matp)<-colnames(resu.padj)<-colnames(resu.r)
color = colorRampPalette(rev(brewer.pal(n = 3.4, name =
                                           "RdGy")))(100)
pheatmap(t(resu.r),color = color,scale = "none", display_numbers = t(matp),fontsize_row=8,fontsize_col=6,border_color=NA,annotation_legend = T)



corr2<-read.csv(file.choose())
simpleNetwork(corr2, Source = 1, Target = 2, height = NULL, width = NULL,
              linkDistance = 50, charge = -30, fontSize = 7, fontFamily = "serif",
              linkColour = "#666", nodeColour = "#3182bd", opacity = 0.6, zoom = F)



# serum vs 9species
cor2data<-read.csv(file.choose())
cormic1<-t(cor2data[c(1:9),-1])
corBA1<-t(cor2data[c(12:38),-1])
resu.r1<-matrix(NA,nrow=27,ncol=9)
resu.p1<-matrix(NA,nrow=27,ncol=9)
for(i in 1:9){
   for(j in 1:27){
      c1<-as.data.frame(cormic1[,i])
      c2<-as.data.frame(corBA1[,j])
      corre<-cor.test(as.numeric(as.matrix(c1)),as.numeric(as.matrix(c2)),alternative = "two.sided","spearman")
      resu.r1[j,i]<-corre$estimate
      resu.p1[j,i]<-corre$p.value
      
   }
}
resu.padj_b<-p.adjust(resu.p1,"bonferroni")
resu.padj_f<-p.adjust(resu.p1,"fdr")
dim(resu.padj_b)<-dim(resu.p1)
dim(resu.padj_f)<-dim(resu.p1)

allr2<-cbind(resu.r1,resu.p1,resu.padj_b,resu.padj_f)

write.csv(allr2,"re2.csv")


x<-cbind(colnames(corBA),t(corBA))
y<-cbind(colnames(cormic),t(cormic))
gramm2<-naiveGramm(y,x,F)
alregramm2<-cbind(gramm2$r,gramm2$p,gramm2$type)

set.seed(100)
xdata<-as.data.frame(x[,-1])
ydata<-as.data.frame(y[,-1])
xname<-x[,1]
yname<-y[,1]
options(warn = -1)
rresult<-presult<-matrix(0,
                         nrow = nrow(xdata),ncol = nrow(ydata))
ltresult<-lnresult<-matrix(0,
                           nrow = nrow(xdata),ncol = nrow(ydata))
rownames(presult) <-rownames(rresult)<-rownames(lnresult)<-xname
colnames(presult)<-colnames(rresult)<-colnames(lnresult)<-yname

for(i in seq_len(nrow(xdata))){
   for(j in seq_len(nrow(ydata))){
      x1<-as.numeric(t(xdata[i,]))
      y1<-as.numeric(t(ydata[j,]))
      lmx<-lm(x1~y1)
      summ<-summary(lmx)
      pvalue<-summ$coefficients[2,4]
      rvalue<-lmx$coefficients[2]*sd(y1)/sd(x1)
      if(pvalue<alpha&rvalue>r){
         presult[i,j]<-pvalue
         rresult[i,j]<-rvalue
         r2value<-rvalue^2
         ltresult[i,j]<-r2value
         lnresult[i,j]<-"linear"
      }
      else{
         #MIC
         warnings('off')
         micxy<-minerva::mine(x1,y1)
         micr<-micxy$MIC
         micp<-0
         for(ii in seq_len(101))
         {
            bootx<-matrix(sample(x1,replace = FALSE))
            booty<-matrix(sample(y1,replace = FALSE))
            tmp<-minerva::mine(bootx,booty)
            MICtp<-tmp$MIC
            tempM<-ifelse(micr<=MICtp,1,0)
            micp<-tempM+micp
         }
         micp<-micp/101
         prsmic<-cor.test(y1,x1,method="pearson")
         prsr<-prsmic$estimate
         ltmic<-1-micr+prsr^2
         presult[i,j]<-micp
         rresult[i,j]<-micr
         lnresult[i,j]<-"nonlinear"
      }
   }
}

results2<-list()
results2[["r"]]<-rresult
results2[["p"]]<-presult
results2[["type"]]<-lnresult

pheatmap(rresult,color = color,fontsize_row=8,border_color=NA,annotation_legend = T,width = 3,height = 5)


