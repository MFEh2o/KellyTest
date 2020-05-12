##GPP data bootstrapping and multiple model testing for the GPP-DOC productivity relationship
##for paper titled 'Shifting limitation of primary production: experimental support for a 
##new model in lake ecosystems"

rm(list=ls())

setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Rscripts/GitHub Materials/')
dataAR=read.csv('gppData.csv', header = T, stringsAsFactors = F, sep = ',')

# fixed-effect hump-shaped model - 12 different ARs and random mesocosm effect
library(nlme)

##multiple model testing
# fixed-effect  - 12 different ARS (quadratic w/ intxns), Random effect of mesosocsm to control for repeated measures
full_ARfixed=lme(GPP~Nuts*AS+DOC+DOC:Nuts+DOC:AS+DOC:Nuts:AS+DOC2+DOC2:Nuts+DOC2:AS+DOC2:Nuts:AS+meso:GPPpre,random=~1|meso,data=dataAR)
full_ARfixed_ml=update(full_ARfixed, . ~ ., method = 'ML')

##no treatment affects but has AR as fixed effect (quadrtic w/out intxns)
doconlyQuadratic_ARfixed=lme(GPP~DOC+DOC2+meso:GPPpre,random=~1|meso,data=dataAR)
doconlyQuadratic_ARfixed_ml=update(doconlyQuadratic_ARfixed, . ~ ., method = 'ML')

#fixed effect linear w/ intxns
fullLinear_ARfixed=lme(GPP~Nuts*AS+DOC+DOC:Nuts+DOC:AS+DOC:Nuts:AS+meso:GPPpre,random=~1|meso,data=dataAR)
fullLinear_ARfixed_ml=update(fullLinear_ARfixed, . ~ ., method = 'ML')

#fixed effect linear w/out intxns
doconlyLinear_ARfixed=lme(GPP~DOC+meso:GPPpre, random=~1|meso, data=dataAR)
doconlyLinear_ARfixed_ml=update(doconlyLinear_ARfixed, . ~ ., method = 'ML')

#fixed effect intercept w/intxns
fullIntercept_ARfixed=lme(GPP~Nuts*AS+meso:GPPpre, random=~1|meso, data=dataAR)
fullIntercept_ARfixed_ml=update(fullIntercept_ARfixed, . ~ ., method = 'ML')

#fixed effect intercept w/out intxns
doconlyIntercept_ARfixed=lme(GPP~meso:GPPpre, random=~1|meso, data=dataAR)
doconlyIntercept_ARfixed_ml=update(doconlyIntercept_ARfixed, . ~ ., method = 'ML')

##run 5 ANOVAs for likelihood ratio test between fullQuadratic AR fixed and 5 other models
anova(doconlyQuadratic_ARfixed_ml, full_ARfixed_ml)
anova(fullLinear_ARfixed_ml, full_ARfixed_ml)
anova(doconlyLinear_ARfixed_ml, full_ARfixed_ml)
anova(fullIntercept_ARfixed_ml, full_ARfixed_ml)
anova(doconlyIntercept_ARfixed_ml, full_ARfixed_ml)

# predicted vs. observed plot
predicted=predict(full_ARfixed)
poFit=lm(predicted~dataAR$GPP) # R2=0.61, p<2e-16; slope = 0.62, intercept= 73
summary(poFit)

##get residuals 
dataAR$resids=residuals(full_ARfixed)

### Bootstrapping globally
# generating bootstrapped simulations and model parameter estimates
BS=array(NA,c(nrow(coefficients(full_ARfixed)),ncol(coefficients(full_ARfixed)),1000))
BSsimGPP=matrix(NA, nrow(dataAR),1000) ##save all new generated GPP estimates for plotting later
for(i in 1:1000){
  BSdata=dataAR
  BSdata$GPP=predicted+sample(dataAR$resids,length(dataAR$resids),replace=FALSE)
  BSsimGPP[,i]=BSdata$GPP
  BS[,,i]=as.matrix(coefficients(lme(GPP~Nuts*AS+DOC+DOC:Nuts+DOC:AS+DOC:Nuts:AS+DOC2+DOC2:Nuts+DOC2:AS+DOC2:Nuts:AS+meso:GPPpre,random=~1|meso,data=BSdata)))
}

BSpars=BS[,c(1:10,23:24),]
BSparsAR=BS[,11:22,]


# generate predicted GPP for each set of BS parameters
DOCpred=seq(5,45,0.1)

design=matrix(c(rep(1,12),
                rep(c(0,1,0,1),each=3),
                rep(c(1,1,0,0),each=3),
                rep(1,12),
                rep(1,12),
                rep(c(0,1,0,1),each=3)*rep(c(1,1,0,0),each=3),
                rep(c(0,1,0,1),each=3),
                rep(c(1,1,0,0),each=3),
                rep(c(0,1,0,1),each=3),
                rep(c(1,1,0,0),each=3),
                rep(c(0,1,0,1),each=3)*rep(c(1,1,0,0),each=3),
                rep(c(0,1,0,1),each=3)*rep(c(1,1,0,0),each=3)),12,12,byrow=FALSE)

colnames(design)=rownames(BSpars)

maxgppBSar=matrix(NA,12,1000)
cdocBSar=maxgppBSar
GPPpredBSar=array(NA,c(12,length(DOCpred),1000))
for(j in 1:1000){
  pars=BSpars[,,j]
  ARpars=BSparsAR[1,,j]
  GPP0=BSsimGPP[dataAR$Day==2,j]
  
  X=design
  X[,1]=pars[,1]
  X[,c(4,7,8,11)]=X[,c(4,7,8,11)]*DOCpred[1]
  X[,c(5,9,10,12)]=X[,c(5,9,10,12)]*DOCpred[1]^2
  
  GPPpredBSar[,1,j]=X%*%c(1,pars[1,2:ncol(pars)])+ARpars*GPP0
  for(i in 2:length(DOCpred)){
    X=design
    X[,1]=pars[,1]
    X[,c(4,7,8,11)]=X[,c(4,7,8,11)]*DOCpred[i]
    X[,c(5,9,10,12)]=X[,c(5,9,10,12)]*DOCpred[i]^2
    GPPpredBSar[,i,j]=X%*%c(1,pars[1,2:ncol(pars)])+ARpars*GPPpredBSar[,(i-1),j]
  }
  
  # skipping first 20 observations because sometimes get higher GPP at super low DOC  and this isn't the vertex of the hupm
  maxgppBSar[,j]=apply(GPPpredBSar[,21:length(DOCpred),j],1,max)
  for(i in 1:12){
    cdocBSar[i,j]=DOCpred[which(GPPpredBSar[i,,j]==maxgppBSar[i,j])]  
  }
}

# generating 95% intervals for bootstrapped predictions of GPP 
asPolyCoordsAR=cbind(c(DOCpred,rev(DOCpred)),rep(NA,length(DOCpred)*2))
bothPolyCoordsAR=asPolyCoordsAR
controlPolyCoordsAR=asPolyCoordsAR
nutsPolyCoordsAR=asPolyCoordsAR

for(i in 1:length(DOCpred)){
  asPolyCoordsAR[c(i,(2*length(DOCpred)+1-i)),2]=quantile(GPPpredBSar[1:3,i,],probs=c(0.025,0.975))
  bothPolyCoordsAR[c(i,(2*length(DOCpred)+1-i)),2]=quantile(GPPpredBSar[4:6,i,],probs=c(0.025,0.975))
  controlPolyCoordsAR[c(i,(2*length(DOCpred)+1-i)),2]=quantile(GPPpredBSar[7:9,i,],probs=c(0.025,0.975))
  nutsPolyCoordsAR[c(i,(2*length(DOCpred)+1-i)),2]=quantile(GPPpredBSar[10:12,i,],probs=c(0.025,0.975))
}

# generated predicted values of model
predDOC=seq(5,45,0.1)
pred_arGPP=matrix(0,length(predDOC),12)
mesos=unique(dataAR$meso)
pars=coef(full_ARfixed)[,c(1:10,23:24)]
parsAR=coef(full_ARfixed)[,11:22]

for(i in 1:length(mesos)){
  cur=dataAR[dataAR$meso==mesos[i],]
  
  nuts=unique(cur$Nuts)
  as=unique(cur$AS)
  
  pred_arGPP[1,i]=pars[i,1]+pars[i,2]*nuts+pars[i,3]*as+pars[i,4]*predDOC[1]+pars[i,5]*predDOC[1]^2+pars[i,6]*nuts*as+pars[i,7]*nuts*predDOC[1]+pars[i,8]*as*predDOC[1]+pars[i,9]*nuts*predDOC[1]^2+pars[i,10]*as*predDOC[1]^2+pars[i,11]*nuts*as*predDOC[1]+pars[i,12]*nuts*as*predDOC[1]^2+parsAR[i,i]*cur$GPPpre[1]
  for(j in 2:nrow(pred_arGPP)){
    pred_arGPP[j,i]=pars[i,1]+pars[i,2]*nuts+pars[i,3]*as+pars[i,4]*predDOC[j]+pars[i,5]*predDOC[j]^2+pars[i,6]*nuts*as+pars[i,7]*nuts*predDOC[j]+pars[i,8]*as*predDOC[j]+pars[i,9]*nuts*predDOC[j]^2+pars[i,10]*as*predDOC[j]^2+pars[i,11]*nuts*as*predDOC[j]+pars[i,12]*nuts*as*predDOC[j]^2+parsAR[i,i]*pred_arGPP[(j-1),i]
  }
}




# boxplots of vertices w/ observed values
obs_maxGPPar=apply(pred_arGPP[21:length(DOCpred),],2,max)
obs_cDOCar=numeric(12)
for(i in 1:12){
  obs_cDOCar[i]=predDOC[which(pred_arGPP[,i]==obs_maxGPPar[i])]
}
names.box=c(expression('low chromo. \nlow P:C stoich.'), expression('low chromo. \nhigh P:C stoich.'), expression('high chromo. \nlow P:C stoich.'), expression('high chromo. \nhigh P:C stoich.'))

##code for regenerating figure 3 in paper
##figure 3a
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('figure3.pdf', width = 12, height = 14)
par(mfrow=c(2,2))
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1)
plot(predDOC,rowMeans(pred_arGPP[,1:3]),type='l',lwd=6,col='goldenrod1',ylim=c(0,650),xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5)
lines(predDOC,rowMeans(pred_arGPP[,4:6]),lwd=6,col='slateblue3')
lines(predDOC,rowMeans(pred_arGPP[,7:9]),lwd=6,col='violetred2')
lines(predDOC,rowMeans(pred_arGPP[,10:12]),lwd=6,col='darkorange1')

polygon(nutsPolyCoordsAR[,1],nutsPolyCoordsAR[,2],col='darkorange1',density=30, border=NA)
polygon(bothPolyCoordsAR[,1],bothPolyCoordsAR[,2],col='slateblue3',density=30, border=NA)
polygon(controlPolyCoordsAR[,1],controlPolyCoordsAR[,2],col='violetred2',density=60, border=NA)
polygon(asPolyCoordsAR[,1],asPolyCoordsAR[,2],col='goldenrod1',density=60, border=NA)
points(21.1,150, pch=19, lwd=10, col = 'goldenrod2', bg='black')
points(17.1,227, pch=19, lwd=10, col = 'slateblue3', bg='black')
points(26.3,278, pch=19, lwd=10, col = 'violetred2', bg='black')
points(25,357, pch=19, lwd=10, col = 'darkorange1', bg='black')
legend(4,650, legend=c("low chromo., low P:C stoich.", "low chromo., high P:C stoich.", 'high chromo., low P:C stoich.', 'high chromo., high P:C stoich.'), 
       col=c("violetred2", "darkorange1", 'goldenrod2', 'slateblue3'), pch=19, lwd=4, pt.cex=1.25, cex=1, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0.25, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "A"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
#figure 3b
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR$GPP,predicted,xlim=c(0,900),ylim=c(0,900),
     xlab=expression('observed GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'),ylab=expression('predicted GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5)
abline(a=0,b=1,lwd=3)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, 13.5), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "B"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
##figure 3c
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1, xpd=FALSE)
boxplot(as.vector(cdocBSar)~rep(rep(c('3','4','1','2'),each=3),1000),outline=FALSE,ylab=expression(' critical DOC'*' (mg'*' L'^-1*')'),xlab="treatment", cex.axis=1.5, cex.lab=1.5, 
        col=c('violetred2', 'darkorange1', 'goldenrod1', 'slateblue3'), xaxt='n')
axis(side = 1, at = seq_along(names.box), labels = names.box, tick = FALSE)
points(rep(c(3,4,1,2),each=3),obs_cDOCar,lwd=2)
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
##figure 3d
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1, xpd=FALSE)
boxplot(as.vector(maxgppBSar)~rep(rep(c('3','4','1','2'),each=3),1000),outline=FALSE,ylab=expression('maximum GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), xlab="treatment", cex.axis=1.5, cex.lab=1.5, 
        col=c('violetred2', 'darkorange1', 'goldenrod1', 'slateblue3'), xaxt='n')
axis(side = 1, at = seq_along(names.box), labels = names.box, tick = FALSE)
points(rep(c(3,4,1,2),each=3),obs_maxGPPar,lwd=2)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "D"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)

dev.off()

##code for generating supplemental figure 3 with observed GPP data with fit from statistical model 
##fits for each treatment
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('figure.s3.pdf', width = 12, height = 14)
par(mfrow=c(2,2))
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1)
plot(dataAR[which(dataAR$Treatment == "C"),8], dataAR[which(dataAR$Treatment == "C"),4],pch=16,col='violetred2',ylim=c(0,max(dataAR$GPP)),xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, xlim=c(min(dataAR$DOC), max(DOCpred)), cex=1.5)
lines(predDOC,rowMeans(pred_arGPP[,7:9]),lwd=6,col='violetred2')
points(26.3,278, pch=19, lwd=10, col = 'black', bg='black')
legend(2,800, legend='low chromo., low P:C stoich.', 
       cex=1.5, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0.25, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "A"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
#figure s3b
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR[which(dataAR$Treatment == "D"),8], dataAR[which(dataAR$Treatment == "D"),4],pch=16,col='darkorange1', ylim=c(0,max(dataAR$GPP)), xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, xlim=c(min(dataAR$DOC), max(DOCpred)), cex=1.5)
lines(predDOC,rowMeans(pred_arGPP[,10:12]),lwd=6,col='darkorange1')
points(25,357, pch=19, lwd=10, col = 'black', bg='black')
legend(2,800, legend='low chromo., high P:C stoich.', 
       cex=1.5, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, 13.5), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "B"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
##figure s3c
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR[which(dataAR$Treatment == "A"),8], dataAR[which(dataAR$Treatment == "A"),4],pch=16,col='goldenrod1',ylim=c(0,max(dataAR$GPP)),xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, xlim=c(min(dataAR$DOC), max(DOCpred)), cex=1.5)
lines(predDOC, rowMeans(pred_arGPP[,1:3]), lwd=6,col='goldenrod1')
points(21.1,150, pch=19, lwd=10, col = 'black', bg='black')
legend(2,750, legend='high chromo., low P:C stoich.', 
       cex=1.5, box.lty = 0)
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
##figure s3d
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR[which(dataAR$Treatment == "B"),8], dataAR[which(dataAR$Treatment == "B"),4],pch=16,col='slateblue3',ylim=c(0,max(dataAR$GPP)),xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, xlim=c(min(dataAR$DOC), max(DOCpred)), cex=1.5)
lines(predDOC,rowMeans(pred_arGPP[,4:6]),lwd=6,col='slateblue3')
points(17.1,227, pch=19, lwd=10, col = 'black', bg='black')
legend(2,750, legend='high chromo., high P:C stoich.', 
      cex=1.5, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "D"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)

dev.off()


##plotting GPP vs. time for each treatment individually
##black line is the day of the 'hump' of HB lake daily GPP
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('figure.gppTime.pdf', width = 12, height = 14)
par(mfrow=c(2,2))
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1)
plot(dataAR[which(dataAR$Treatment == "C"),3], dataAR[which(dataAR$Treatment == "C"),4],pch=16,col='violetred2',ylim=c(0,max(dataAR$GPP)),xlab=expression('day of experiment'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, cex=1.5)
abline(v=49, col='violetred2', lwd=3)
abline(v=32, col='black', lwd=3)
legend(2,800, legend='low chromo., low P:C stoich.', 
       cex=1.5, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0.25, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "A"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
#figure s3b
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR[which(dataAR$Treatment == "D"),3], dataAR[which(dataAR$Treatment == "D"),4],pch=16,col='darkorange1', ylim=c(0,max(dataAR$GPP)), xlab=expression('day of experiment'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, cex=1.5)
abline(v=46, col='darkorange1', lwd=3)
abline(v=32, col='black', lwd=3)
legend(2,800, legend='low chromo., high P:C stoich.', 
       cex=1.5, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, 13.5), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "B"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)
##figure s3c
par(mar=c(5.1, 7, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR[which(dataAR$Treatment == "A"),3], dataAR[which(dataAR$Treatment == "A"),4],pch=16,col='goldenrod1',ylim=c(0,max(dataAR$GPP)),xlab=expression('day of experiment'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, cex=1.5)
abline(v=38, col='goldenrod1', lwd=3)
abline(v=32, col='black', lwd=3)
legend(2,750, legend='high chromo., low P:C stoich.', 
       cex=1.5, box.lty = 0)
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
##figure s3d
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1, xpd=FALSE)
plot(dataAR[which(dataAR$Treatment == "B"),3], dataAR[which(dataAR$Treatment == "B"),4],pch=16,col='slateblue3',ylim=c(0,max(dataAR$GPP)),xlab=expression('day of experiment'),
     ylab=expression('GPP'*' (mg'*' C'*' m'^-2*' day'^-1*')'), cex.axis=1.5, cex.lab=1.5, cex=1.5)
abline(v=34, col='slateblue3', lwd=3)
abline(v=32, col='black', lwd=3)
legend(2,750, legend='high chromo., high P:C stoich.', 
       cex=1.5, box.lty = 0)
par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "D"
x <- x[1] + strwidth(txt, cex=2.8) 
y <- y[2] - strheight(txt, cex=2.8) 
text(x, y, txt, cex=2.8)

dev.off()
