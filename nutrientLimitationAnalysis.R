##Script for mesocosm SRP and nitrate model testing for nturient limitation
##for paper titled 'Shifting limitation of primary production: experimental support for a 
##new model in lake ecosystems"

rm(list=ls())

setwd("C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Rscripts/GitHub Materials")
#setwd("~/Documents/Research/People/Students/current/Olson_Carly/mesocosms/manuscript/EcoLetters_Submitted/Revisions/nutrient_limitation_analyses/")
final=read.csv('compiledNutrientLightData.csv', header = T, stringsAsFactors = F)

##data frame in form to use with ANOVAs
dataAnov=final
dataAnov=dataAnov[order(dataAnov$Date.Sample, decreasing = FALSE), ]
dataAnov$time=rep(1:10, each=12)
dataAnov$time=as.factor(dataAnov$time)
dataAnov$AS=rep(NA)
dataAnov$Nuts=rep(NA)
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
dataAnov$DOC2=dataAnov$DOC*dataAnov$DOC
dataAnov$DOY=as.numeric(strftime(dataAnov$Date.Sample, format = '%j'))
dataAnov$DOY2=dataAnov$DOY*dataAnov$DOY

# can't have 0's with gamma distribution...
dataAnov[dataAnov[,7]==0,7]=0.1

#########piecewise regression on srp data vs. doc (acute)
###A=+aquashade, b=+aquashade & +nutrients, c=control (no treatments administered), d=+nutrients
##a=acute (model 3), o=obtuse (model 4), l=linear (model 2), i=intercept (model 1)
##p=srp, n=nitrate


##model testing: srp vs. doc

# intercept (model 1)
intFit<-function(p,x,y){
  sum(-dgamma(x=y,scale=p[1]/exp(p[2]),shape=exp(p[2]),log=TRUE))
}
guess=c(4,1)
fitA.p.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),7])
fitB.p.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),7])
fitC.p.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),7])
fitD.p.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),7])

#linear (model 2)
linReg<-function(p,x,y){
  yhat=p[1]*x+p[2]
  shape=yhat/exp(p[3])
  sum(-dgamma(x=y,shape=shape,scale=exp(p[3]), log=TRUE))
}
guess=c(1.2,0.5,1)
fitA.p.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),7])
fitB.p.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),7])
fitC.p.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),7])
fitD.p.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),7])

##model 3
bpReg<-function(p,x,y){
  yhat=ifelse(x<p[1],p[2],p[3]*x+(p[2]-p[3]*p[1]))
  shape=yhat/exp(p[4])
  sum(-dgamma(x=y,shape=shape,scale=exp(p[4]),log=TRUE))
}
guess=c(15,2,0.5,1)
fitA.p.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),7])
fitB.p.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),7])
fitC.p.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),7])
fitD.p.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),7])

#model 4
bpReg2<-function(p,x,y){
  yhat=ifelse(x<p[1], p[2]*x+p[3], p[4]*x+(-p[1]*p[4]+(p[3]+p[1]*p[2])))
  shape=yhat/exp(p[5])
  sum(-dgamma(x=y,shape=shape,scale=exp(p[5]),log=TRUE))
}
guess=c(15,0.3,1,2,1)
fitA.p.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),7])
fitB.p.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),7])
fitC.p.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),7])
fitD.p.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),7])

##model comparison for srp vs. doc MOdels 1-4

##+aquashade treatment (A)
#AICs of each model
fitA.p.i$value*2+2*length(fitA.p.i$par)
fitA.p.l$value*2+2*length(fitA.p.l$par)
fitA.p.a$value*2+2*length(fitA.p.a$par)
fitA.p.o$value*2+2*length(fitA.p.o$par)

##plot fits from the four models
plot(dataAnov[dataAnov$Site=="A",5],dataAnov[dataAnov$Site=="A",7])
lines(4:44,rep(fitA.p.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitA.p.l$par[2]+fitA.p.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitA.p.a$par[1],fitA.p.a$par[2],fitA.p.a$par[3]*(4:44)+(fitA.p.a$par[2]-fitA.p.a$par[3]*fitA.p.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitA.p.o$par[1],fitA.p.o$par[2]*(4:44)+fitA.p.o$par[3],fitA.p.o$par[4]*(4:44)+(-fitA.p.o$par[1]*fitA.p.o$par[4]+(fitA.p.o$par[3]+fitA.p.o$par[1]*fitA.p.o$par[2]))),lwd=2,col='purple')

##+aquashade & +nutrients treatment (B)
#AICs of each model
fitB.p.i$value*2+2*length(fitB.p.i$par)
fitB.p.l$value*2+2*length(fitB.p.l$par)
fitB.p.a$value*2+2*length(fitB.p.a$par)
fitB.p.o$value*2+2*length(fitB.p.o$par)

#plot fits from the four models
plot(dataAnov[dataAnov$Site=="B",5],dataAnov[dataAnov$Site=="B",7])
lines(4:44,rep(fitB.p.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitB.p.l$par[2]+fitB.p.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitB.p.a$par[1],fitB.p.a$par[2],fitB.p.a$par[3]*(4:44)+(fitB.p.a$par[2]-fitB.p.a$par[3]*fitB.p.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitB.p.o$par[1],fitB.p.o$par[2]*(4:44)+fitB.p.o$par[3],fitB.p.o$par[4]*(4:44)+(-fitB.p.o$par[1]*fitB.p.o$par[4]+(fitB.p.o$par[3]+fitB.p.o$par[1]*fitB.p.o$par[2]))),lwd=2,col='purple')

#control (C)
# AICs of each model
fitC.p.i$value*2+2*length(fitC.p.i$par)
fitC.p.l$value*2+2*length(fitC.p.l$par)
fitC.p.a$value*2+2*length(fitC.p.a$par)
fitC.p.o$value*2+2*length(fitC.p.o$par)

#plot fits from four models
plot(dataAnov[dataAnov$Site=="C",5],dataAnov[dataAnov$Site=="C",7])
lines(4:44,rep(fitC.p.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitC.p.l$par[2]+fitC.p.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitC.p.a$par[1],fitC.p.a$par[2],fitC.p.a$par[3]*(4:44)+(fitC.p.a$par[2]-fitC.p.a$par[3]*fitC.p.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitC.p.o$par[1],fitC.p.o$par[2]*(4:44)+fitC.p.o$par[3],fitC.p.o$par[4]*(4:44)+(-fitC.p.o$par[1]*fitC.p.o$par[4]+(fitC.p.o$par[3]+fitC.p.o$par[1]*fitC.p.o$par[2]))),lwd=2,col='purple')

#+nutrients treatment (D)
# AICs of each model
fitD.p.i$value*2+2*length(fitD.p.i$par)
fitD.p.l$value*2+2*length(fitD.p.l$par)
fitD.p.a$value*2+2*length(fitD.p.a$par)
fitD.p.o$value*2+2*length(fitD.p.o$par)

##plots fits from four models
plot(dataAnov[dataAnov$Site=="D",5],dataAnov[dataAnov$Site=="D",7])
lines(4:44,rep(fitD.p.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitD.p.l$par[2]+fitD.p.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitD.p.a$par[1],fitD.p.a$par[2],fitD.p.a$par[3]*(4:44)+(fitD.p.a$par[2]-fitD.p.a$par[3]*fitD.p.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitD.p.o$par[1],fitD.p.o$par[2]*(4:44)+fitD.p.o$par[3],fitD.p.o$par[4]*(4:44)+(-fitD.p.o$par[1]*fitD.p.o$par[4]+(fitD.p.o$par[3]+fitD.p.o$par[1]*fitD.p.o$par[2]))),lwd=2,col='purple')

# put all AICs for SRP with gamma likelihood in matrix for comparing to normal likelihood
AICs_gammaP=cbind(
  c(fitA.p.i$value*2+2*length(fitA.p.i$par),
    fitA.p.l$value*2+2*length(fitA.p.l$par),
    fitA.p.a$value*2+2*length(fitA.p.a$par),
    fitA.p.o$value*2+2*length(fitA.p.o$par)),
  c(fitB.p.i$value*2+2*length(fitB.p.i$par),
    fitB.p.l$value*2+2*length(fitB.p.l$par),
    fitB.p.a$value*2+2*length(fitB.p.a$par),
    fitB.p.o$value*2+2*length(fitB.p.o$par)),
  c(fitC.p.i$value*2+2*length(fitC.p.i$par),
    fitC.p.l$value*2+2*length(fitC.p.l$par),
    fitC.p.a$value*2+2*length(fitC.p.a$par),
    fitC.p.o$value*2+2*length(fitC.p.o$par)),
  c(fitD.p.i$value*2+2*length(fitD.p.i$par),
    fitD.p.l$value*2+2*length(fitD.p.l$par),
    fitD.p.a$value*2+2*length(fitD.p.a$par),
    fitD.p.o$value*2+2*length(fitD.p.o$par))
)

colnames(AICs_gammaP)=c("A","B","C","D")
rownames(AICs_gammaP)=c("intercept","linear","acuteBP","obtuseBP")



#####nitrate vs. DOC

#model 1
guess=c(300,1)
fitA.n.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),8])
fitB.n.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),8])
fitC.n.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),8])
fitD.n.i=optim(guess,intFit,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),8])

#model 2
guess=c(30,0,1)
fitA.n.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),8])
fitB.n.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),8])
fitC.n.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),8])
fitD.n.l=optim(guess,linReg,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),8])

##model 3
guess=c(25,50,30,1)
fitA.n.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),8])
fitB.n.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),8])
fitC.n.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),8])
fitD.n.a=optim(guess,bpReg,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),8])

##model 4
guess=c(25,10,-20,30,1)
fitA.n.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='A'),5],y=dataAnov[which(dataAnov$Site=='A'),8])
fitB.n.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='B'),5],y=dataAnov[which(dataAnov$Site=='B'),8])
fitC.n.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='C'),5],y=dataAnov[which(dataAnov$Site=='C'),8])
fitD.n.o=optim(guess,bpReg2,x=dataAnov[which(dataAnov$Site=='D'),5],y=dataAnov[which(dataAnov$Site=='D'),8])


##model comparison for nitrate vs. doc models 1-4

##+aquashade treatment (A)
# AICs for each model
fitA.n.i$value*2+2*length(fitA.n.i$par)
fitA.n.l$value*2+2*length(fitA.n.l$par)
fitA.n.a$value*2+2*length(fitA.n.a$par)
fitA.n.o$value*2+2*length(fitA.n.o$par)

##plot model fits
plot(dataAnov[dataAnov$Site=="A",5],dataAnov[dataAnov$Site=="A",8])
lines(4:44,rep(fitA.n.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitA.n.l$par[2]+fitA.n.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitA.n.a$par[1],fitA.n.a$par[2],fitA.n.a$par[3]*(4:44)+(fitA.n.a$par[2]-fitA.n.a$par[3]*fitA.n.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitA.n.o$par[1],fitA.n.o$par[2]*(4:44)+fitA.n.o$par[3],fitA.n.o$par[4]*(4:44)+(-fitA.n.o$par[1]*fitA.n.o$par[4]+(fitA.n.o$par[3]+fitA.n.o$par[1]*fitA.n.o$par[2]))),lwd=2,col='purple')

#+aquashade & +nutrient treatment (B)
# AICs for each model
fitB.n.i$value*2+2*length(fitB.n.i$par)
fitB.n.l$value*2+2*length(fitB.n.l$par)
fitB.n.a$value*2+2*length(fitB.n.a$par)
fitB.n.o$value*2+2*length(fitB.n.o$par)

##plot model fits
plot(dataAnov[dataAnov$Site=="B",5],dataAnov[dataAnov$Site=="B",8])
lines(4:44,rep(fitB.n.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitB.n.l$par[2]+fitB.n.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitB.n.a$par[1],fitB.n.a$par[2],fitB.n.a$par[3]*(4:44)+(fitB.n.a$par[2]-fitB.n.a$par[3]*fitB.n.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitB.n.o$par[1],fitB.n.o$par[2]*(4:44)+fitB.n.o$par[3],fitB.n.o$par[4]*(4:44)+(-fitB.n.o$par[1]*fitB.n.o$par[4]+(fitB.n.o$par[3]+fitB.n.o$par[1]*fitB.n.o$par[2]))),lwd=2,col='purple')

##control (C)
# AICs for each model
fitC.n.i$value*2+2*length(fitC.n.i$par)
fitC.n.l$value*2+2*length(fitC.n.l$par)
fitC.n.a$value*2+2*length(fitC.n.a$par)
fitC.n.o$value*2+2*length(fitC.n.o$par)

##plot model fits
plot(dataAnov[dataAnov$Site=="C",5],dataAnov[dataAnov$Site=="C",8])
lines(4:44,rep(fitC.n.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitC.n.l$par[2]+fitC.n.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitC.n.a$par[1],fitC.n.a$par[2],fitC.n.a$par[3]*(4:44)+(fitC.n.a$par[2]-fitC.n.a$par[3]*fitC.n.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitC.n.o$par[1],fitC.n.o$par[2]*(4:44)+fitC.n.o$par[3],fitC.n.o$par[4]*(4:44)+(-fitC.n.o$par[1]*fitC.n.o$par[4]+(fitC.n.o$par[3]+fitC.n.o$par[1]*fitC.n.o$par[2]))),lwd=2,col='purple')

##+nutrients treatment (D)
# AICs for each model
fitD.n.i$value*2+2*length(fitD.n.i$par)
fitD.n.l$value*2+2*length(fitD.n.l$par)
fitD.n.a$value*2+2*length(fitD.n.a$par)
fitD.n.o$value*2+2*length(fitD.n.o$par)

##plot model fits
plot(dataAnov[dataAnov$Site=="D",5],dataAnov[dataAnov$Site=="D",8])
lines(4:44,rep(fitD.n.i$par[1],length(4:44)),lwd=2)
lines(4:44,fitD.n.l$par[2]+fitD.n.l$par[1]*(4:44),lwd=2,col='red')
lines(4:44,ifelse((4:44)<fitD.n.a$par[1],fitD.n.a$par[2],fitD.n.a$par[3]*(4:44)+(fitD.n.a$par[2]-fitD.n.a$par[3]*fitD.n.a$par[1])),lwd=2,col='green')
lines(4:44,ifelse((4:44)<fitD.n.o$par[1],fitD.n.o$par[2]*(4:44)+fitD.n.o$par[3],fitD.n.o$par[4]*(4:44)+(-fitD.n.o$par[1]*fitD.n.o$par[4]+(fitD.n.o$par[3]+fitD.n.o$par[1]*fitD.n.o$par[2]))),lwd=2,col='purple')

##nitrate AIC values
AICs_gammaN=cbind(
  c(fitA.n.i$value*2+2*length(fitA.n.i$par),
    fitA.n.l$value*2+2*length(fitA.n.l$par),
    fitA.n.a$value*2+2*length(fitA.n.a$par),
    fitA.n.o$value*2+2*length(fitA.n.o$par)),
  c(fitB.n.i$value*2+2*length(fitB.n.i$par),
    fitB.n.l$value*2+2*length(fitB.n.l$par),
    fitB.n.a$value*2+2*length(fitB.n.a$par),
    fitB.n.o$value*2+2*length(fitB.n.o$par)),
  c(fitC.n.i$value*2+2*length(fitC.n.i$par),
    fitC.n.l$value*2+2*length(fitC.n.l$par),
    fitC.n.a$value*2+2*length(fitC.n.a$par),
    fitC.n.o$value*2+2*length(fitC.n.o$par)),
  c(fitD.n.i$value*2+2*length(fitD.n.i$par),
    fitD.n.l$value*2+2*length(fitD.n.l$par),
    fitD.n.a$value*2+2*length(fitD.n.a$par),
    fitD.n.o$value*2+2*length(fitD.n.o$par))
)

colnames(AICs_gammaN)=c("A","B","C","D")
rownames(AICs_gammaN)=c("intercept","linear","acuteBP","obtuseBP")

###########figure S5
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


##srp 
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('supp5a.pdf', width = 6, height = 7)
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1)
plot(dataAnov$DOC,dataAnov$SRP,pch=pchs,col=colors,lwd=lwds,cex=1.75,
     xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'), ylab=expression('soluble reactive phosphorus'*' (ug'*' L'^-1*')'), cex.axis=1.5, cex.lab=1.5, ylim=c(0, 25))
lines(4:48,ifelse((4:48)<fitA.p.a$par[1],fitA.p.a$par[2],fitA.p.a$par[3]*(4:48)+(fitA.p.a$par[2]-fitA.p.a$par[3]*fitA.p.a$par[1])),lwd=3,col='goldenrod1')
lines(4:48,ifelse((4:48)<fitB.p.a$par[1],fitB.p.a$par[2],fitB.p.a$par[3]*(4:48)+(fitB.p.a$par[2]-fitB.p.a$par[3]*fitB.p.a$par[1])),lwd=3,col='slateblue3')
lines(4:48,rep(fitC.p.i$par[1],length(4:48)),lwd=3, col='violetred2')
lines(4:48,ifelse((4:48)<fitD.p.a$par[1],fitD.p.a$par[2],fitD.p.a$par[3]*(4:48)+(fitD.p.a$par[2]-fitD.p.a$par[3]*fitD.p.a$par[1])),lwd=3,col='darkorange1')
legend(3, 25, legend=c("low chromo., low P:C stoich.", "low chromo., high P:C stoich.", 'high chromo., low P:C stoich.', 'high chromo., high P:C stoich.'), 
       col=c("violetred2", "darkorange1", 'goldenrod1', 'slateblue3'), pch = 15, lwd=3, pt.cex=2, cex=1.25, box.lty = 0)
dev.off()

#nitrate
setwd('C:/Users/Carly/Documents/MFE/UNDERC 2018/Mesocosm 2018/Paper/Revisions/Figures/')
pdf('supp5b.pdf', width = 6, height = 7)
par(mar=c(5.1, 5, 3.1, 2.1) + 0.1)
plot(dataAnov$DOC,dataAnov$nitrate,pch=pchs,col=colors,lwd=lwds,cex=1.75,
     xlab=expression('dissolved organic carbon'*' (mg'*' L'^-1*')'), ylab=expression('nitrate'*' (ug'*' L'^-1*')'), cex.axis=1.5, cex.lab=1.5, ylim=c(0, 2000))
lines(4:48,ifelse((4:48)<fitA.n.a$par[1],fitA.n.a$par[2],fitA.n.a$par[3]*(4:48)+(fitA.n.a$par[2]-fitA.n.a$par[3]*fitA.n.a$par[1])),lwd=3,col='goldenrod1')
lines(4:48,ifelse((4:48)<fitB.n.a$par[1],fitB.n.a$par[2],fitB.n.a$par[3]*(4:48)+(fitB.n.a$par[2]-fitB.n.a$par[3]*fitB.n.a$par[1])),lwd=3,col='slateblue3')
lines(4:48,ifelse((4:48)<fitC.n.o$par[1],fitC.n.o$par[2]*(4:48)+fitC.n.o$par[3],fitC.n.o$par[4]*(4:48)+(-fitC.n.o$par[1]*fitC.n.o$par[4]+(fitC.n.o$par[3]+fitC.n.o$par[1]*fitC.n.o$par[2]))),lwd=3,col='violetred2')
lines(4:48,ifelse((4:48)<fitD.n.o$par[1],fitD.n.o$par[2]*(4:48)+fitD.n.o$par[3],fitD.n.o$par[4]*(4:48)+(-fitD.n.o$par[1]*fitD.n.o$par[4]+(fitD.n.o$par[3]+fitD.n.o$par[1]*fitD.n.o$par[2]))),lwd=3,col='darkorange1')
legend(3, 2000, legend=c("low chromo., low P:C stoich.", "low chromo., high P:C stoich.", 'high chromo., low P:C stoich.', 'high chromo., high P:C stoich.'), 
       col=c("violetred2", "darkorange1", 'goldenrod1', 'slateblue3'), pch = 15, lwd=3, pt.cex=2, cex=1.25, box.lty = 0)
dev.off()



