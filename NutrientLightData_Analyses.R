##Nutrient and light Data Analyses for paper titled 'Shifting limitation of primary production: experimental support for a 
##new model in lake ecosystems"

rm(list=ls())
library(plyr)
library(xts)

##read in data files
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Rscripts/GitHub Materials/')
final=read.csv('compiledNutrientLightData.csv', header = T, stringsAsFactors = F, sep = ',')

##calculating the mean of all the variables and storing in new data frame called 'avgs' 
for (i in 5:ncol(final)) {
  if (i==5) {
    avgs=aggregate(final[,i], by=list(final$Site, final$Date.Sample), FUN = mean, na.rm = T)
    colnames(avgs)=c('Site', 'Date.Sample', 'DOC')
  } else {
    temp=aggregate(final[,i], by=list(final$Site, final$Date.Sample), FUN = mean, na.rm = T)
    avgs[,i-2]=temp[,3]
  }
}
avgs=avgs[order(avgs$Date.Sample),]
colnames(avgs)[4:10]=c('TP', 'SRP', 'nitrate', 'TN', 'chlorophyll', 'kd', 'Ibar')

##calculating stnd deviation of all the variables and storing in new data frame 'stdevs' 
for (i in 5:ncol(final)) {
  if (i==5) {
    stdevs=aggregate(final[,i], by=list(final$Site, final$Date.Sample), FUN = sd, na.rm = T)
    colnames(stdevs)=c('Site', 'Date.Sample', 'DOC')
  } else {
    temp=aggregate(final[,i], by=list(final$Site, final$Date.Sample), FUN = sd, na.rm = T)
    stdevs[,i-2]=temp[,3]
  }
}
stdevs=stdevs[order(stdevs$Date.Sample),]
colnames(stdevs)[4:10]=c('TP', 'SRP', 'nitrate', 'TN', 'chlorophyll', 'kd', 'Ibar')

##calculating difference between the average DOC concentration of high and low chromophoric treatments to find
##average effect of aquashade on DOC concentration
AS=avgs[which(avgs$Site == 'A' | avgs$Site == 'B'),]
meanAS=aggregate(AS$DOC, by=list(AS$Date.Sample), FUN=mean)
non=avgs[which(avgs$Site == 'C' | avgs$Site == 'D'),]
meanNon=aggregate(non$DOC, by=list(non$Date.Sample), FUN=mean)
diff=meanAS$x-meanNon$x
avgDiff=mean(diff)

##getting table prepped to run repeated measures anova and mixed-effect ancovas
dataAnov=final
dataAnov=dataAnov[order(dataAnov$Date.Sample, decreasing = FALSE), ]
dataAnov$time=rep(1:10, each=12)
dataAnov$time=as.factor(dataAnov$time)
dataAnov$as=rep(NA)
dataAnov$nuts=rep(NA)
for (i in 1:nrow(dataAnov)) {
  if (dataAnov[i,3] == 'A') {
    dataAnov[i,14]=1
    dataAnov[i,15]=0
  } 
  if (dataAnov[i,3] == 'B') {
    dataAnov[i,14]=1
    dataAnov[i,15]=1
  } 
  if (dataAnov[i,3] == 'C') {
    dataAnov[i,14]=0
    dataAnov[i,15]=0
  } 
  if (dataAnov[i,3] == 'D') {
    dataAnov[i,14]=0
    dataAnov[i,15]=1
  }
}
dataAnov$as=as.factor(dataAnov$as)
dataAnov$nuts=as.factor(dataAnov$nuts)
dataAnov$Lake.ID=as.factor(dataAnov$Lake.ID)

##repeated measures anova for DOC with random mesocosm effect
library(nlme)
docAnova=lme(DOC ~ nuts*as*time,random=~1|Lake.ID, data=dataAnov)

### ANCOVA with DOC for all other nutrient variables
nutrients=c('TP', 'SRP', 'Nitrate', 'TN', 'Chlorophyll', 'kd', 'Ibar')
pars=array(NA,c(12,8,7))

TPlm=lme(TP~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$TP),])
pars[,,1]=as.matrix(coef(TPlm))

SRPlm=lme(SRP~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$SRP),])
pars[,,2]=as.matrix(coef(SRPlm))

NITlm=lme(nitrate~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$nitrate),])
pars[,,3]=as.matrix(coef(NITlm))

TNlm=lme(TN~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$TN),])
pars[,,4]=as.matrix(coef(TNlm))

CHLlm=lme(chlorophyll~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$chlorophyll),])
pars[,,5]=as.matrix(coef(CHLlm))

KDlm=lme(kd~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$kd),])
pars[,,6]=as.matrix(coef(KDlm))

LITlm=lme(Ibar~nuts*as*DOC,random=~1|Lake.ID,data=dataAnov[!is.na(dataAnov$Ibar),])
pars[,,7]=as.matrix(coef(LITlm))

names=c('TP (ug/L)', 'SRP (ug/L)', 'Nitrate (ug/L)', 'TN (ug/L)', 'Chlorophyll (ug/L)', 'kd (m-1)', 'Ibar')
maxY=c(48, 25, 1500, 2200, 40, 9, 1)
treats=c('A', 'B', 'C', 'D')

# for generating predicted values
predDOC=seq(5,45,0.1)

design=matrix(c(
  int=rep(1,12),
  nuts=c(0,1,0,0,1,1,0,1,1,1,0,0),
  as=c(1,1,1,1,0,1,0,1,0,0,0,0),
  doc=rep(1,12),
  nutsAs=c(0,1,0,0,0,1,0,1,0,0,0,0),
  nutsDoc=c(0,1,0,0,1,1,0,1,1,1,0,0),
  asDoc=c(1,1,1,1,0,1,0,1,0,0,0,0),
  nutsAsDoc=c(0,1,0,0,0,1,0,1,0,0,0,0)),nrow=12,byrow=FALSE)

predTP=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,1]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predTP[,i]=X%*%c(1,pars[1,2:8,1])
}

predSRP=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,2]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predSRP[,i]=X%*%c(1,pars[1,2:8,2])
}

predNIT=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,3]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predNIT[,i]=X%*%c(1,pars[1,2:8,3])
}

predTN=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,4]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predTN[,i]=X%*%c(1,pars[1,2:8,4])
}

predCHL=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,5]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predCHL[,i]=X%*%c(1,pars[1,2:8,5])
}

predKD=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,6]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predKD[,i]=X%*%c(1,pars[1,2:8,6])
}

predLIT=matrix(NA,12,length(predDOC))
for(i in 1:length(predDOC)){
  X=design
  X[,1]=pars[,1,7]
  X[,c(4,6:8)]=X[,c(4,6:8)]*predDOC[i]
  predLIT[,i]=X%*%c(1,pars[1,2:8,7])
}

##########b##############################################################################
##fgenerating figure 2 for the paper
colors=rep('goldenrod1',nrow(dataAnov))
colors[dataAnov$Site=="B"]='slateblue3'
colors[dataAnov$Site=="C"]='violetred2'
colors[dataAnov$Site=="D"]='darkorange1'

pchs=rep(15,nrow(dataAnov))
pchs[dataAnov$Site=="B"]=15
pchs[dataAnov$Site=="C"]=15
pchs[dataAnov$Site=="D"]=15

lwds=rep(3,nrow(dataAnov))
lwds[dataAnov$Site=="C"]=3
lwds[dataAnov$Site=="D"]=3

##putting Figure 2A-C on one pdf
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('figure2.pdf', width = 6, height = 21)
par(mfrow=c(3,1))
##figure 3a
par(mar=c(5.1,7, 3.1, 2.1) + 0.1)
plot(as.Date(avgs[which(avgs$Site == 'A'),2]), avgs[which(avgs$Site == 'A'),3], 
     col = 'goldenrod1', ylim=c(5,49), pch=15, cex=1.75, type='b', lwd=3, 
     xlab='time (month-day)', ylab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'), cex.lab=1.75, cex.axis=1.75, xaxt = 'n')
arrows(as.Date(avgs[which(avgs$Site == 'A'),2]), avgs[which(avgs$Site == 'A'),3]-stdevs[which(stdevs$Site == 'A'),3], as.Date(avgs[which(avgs$Site == 'A'),2]), avgs[which(avgs$Site == 'A'),3]+stdevs[which(stdevs$Site == 'A'),3], length=0.05, angle=90, code=3, lwd=2)
points(as.Date(avgs[which(avgs$Site == 'B'),2]), avgs[which(avgs$Site == 'B'),3], col = 'slateblue3', pch=15, cex=1.75, type='b', lwd=3)
arrows(as.Date(avgs[which(avgs$Site == 'B'),2]), avgs[which(avgs$Site == 'B'),3]-stdevs[which(stdevs$Site == 'B'),3], as.Date(avgs[which(avgs$Site == 'B'),2]), avgs[which(avgs$Site == 'B'),3]+stdevs[which(stdevs$Site == 'B'),3], length=0.05, angle=90, code=3, lwd=2)
points(as.Date(avgs[which(avgs$Site == 'C'),2]), avgs[which(avgs$Site == 'C'),3], col = 'violetred2', pch=15, cex=1.75, type='b', lwd=3)
arrows(as.Date(avgs[which(avgs$Site == 'C'),2]), avgs[which(avgs$Site == 'C'),3]-stdevs[which(stdevs$Site == 'C'),3], as.Date(avgs[which(avgs$Site == 'C'),2]), avgs[which(avgs$Site == 'C'),3]+stdevs[which(stdevs$Site == 'C'),3], length=0.05, angle=90, code=3, lwd=2)
points(as.Date(avgs[which(avgs$Site == 'D'),2]), avgs[which(avgs$Site == 'D'),3], col = 'darkorange1', pch=15, cex=1.75, type='b', lwd=3)
arrows(as.Date(avgs[which(avgs$Site == 'D'),2]), avgs[which(avgs$Site == 'D'),3]-stdevs[which(stdevs$Site == 'D'),3], as.Date(avgs[which(avgs$Site == 'D'),2]), avgs[which(avgs$Site == 'D'),3]+stdevs[which(stdevs$Site == 'D'),3], length=0.05, angle=90, code=3, lwd=2)
legend(as.Date('2018-05-30'), 50, legend=c("low chromo., low P:C stoich.", "low chromo., high P:C stoich.", 'high chromo., low P:C stoich.', 'high chromo., high P:C stoich.'), 
       col=c("violetred2", "darkorange1", 'goldenrod1', 'slateblue3'), pch = 15, lwd=3, pt.cex=2, cex=1.25, box.lty = 0)
axis.Date(side = 1, x = as.Date(avgs[which(avgs$Site == 'A'),2]),
          format= "%m-%d", cex.axis=1.5)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0.25, di[1]), from="in", to="user")
y <- grconvertY(c(0, 20.5), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "A"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
##figure 3b
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAnov$DOC,dataAnov$TP,pch=pchs,col=colors,lwd=lwds,cex=1.75,
     xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'), ylab=expression('total phosphorus'*' (ug'*' L'^-1*')'), cex.axis=1.75, cex.lab=1.75)
lines(predDOC,colMeans(predTP[c(1,3,4),]),lwd=3,col='goldenrod1')
lines(predDOC,colMeans(predTP[c(2,6,8),]),lwd=3,col='slateblue3')
lines(predDOC,colMeans(predTP[c(7,11,12),]),lwd=3,col='violetred2')
lines(predDOC,colMeans(predTP[c(5,9,10),]),lwd=3,col='darkorange1')
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0.25, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "B"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)

##figure 3c
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAnov$DOC,dataAnov$kd,pch=pchs,col=colors,lwd=lwds,cex=1.75,
     xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'), ylab=expression('vertical light attenuation'*' (m'^-1*')'), cex.axis=1.75, cex.lab=1.75)
lines(predDOC,colMeans(predKD[c(1,3,4),]),lwd=3,col='goldenrod1')
lines(predDOC,colMeans(predKD[c(2,6,8),]),lwd=3,col='slateblue3')
lines(predDOC,colMeans(predKD[c(7,11,12),]),lwd=3,col='violetred2')
lines(predDOC,colMeans(predKD[c(5,9,10),]),lwd=3,col='darkorange1')
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0.25, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "C"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)

dev.off()

##supplementary figure S2 (TN w. ANCOVA fit)
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('supp2.pdf', width = 6, height = 7)
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1)
plot(dataAnov$DOC,dataAnov$TN,pch=pchs,col=colors,lwd=lwds,cex=1.75,
     xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'), ylab=expression('total nitrogen'*' (ug'*' L'^-1*')'), cex.axis=1.5, cex.lab=1.5)
lines(predDOC,colMeans(predTN[c(1,3,4),]),lwd=3,col='goldenrod1')
lines(predDOC,colMeans(predTN[c(2,6,8),]),lwd=3,col='slateblue3')
lines(predDOC,colMeans(predTN[c(7,11,12),]),lwd=3,col='violetred2')
lines(predDOC,colMeans(predTN[c(5,9,10),]),lwd=3,col='darkorange1')
legend(3, 2100, legend=c("low chromo., low P:C stoich.", "low chromo., high P:C stoich.", 'high chromo., low P:C stoich.', 'high chromo., high P:C stoich.'), 
       col=c("violetred2", "darkorange1", 'goldenrod1', 'slateblue3'), pch = 15, lwd=3, pt.cex=2, cex=1.25, box.lty = 0)
dev.off()
