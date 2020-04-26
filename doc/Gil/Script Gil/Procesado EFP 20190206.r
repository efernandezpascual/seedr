source("funciones.R")

classlist<-function(...,class=class)
{
  ls<-list(...)
  class(ls)<-class
  return(ls)
}

# data: data.frame with columns Time, G, PG
# Returns CDFs and interpolated CDFs for all cuts
cdf<-function(data)
{
   data<-rbind(data[1,],data)
   data[1,c("Time","G")]<-c(0,0)    # Force a row with Time=0
   f<-factor(data$Time)    # For arbitrary times

   data$Time<-as.numeric(f)
   nTime<-nlevels(f)
   cfd<-data.frame(Time=as.numeric(levels(f)), Prob=rep(0,nTime))

   sm<-get.groups(data,data$Dish)
   PG<-0
   for (i in 1:length(sm))
   {
     if (nrow(sm[[i]]))
     {
       PG<-PG+max(sm[[i]]$PG)
       for (j in 1:nTime)
       {
         if (sum(sm[[i]]$Time<=j)>0)
cfd$Prob[j]<-cfd$Prob[j]+max(sm[[i]]$G[sm[[i]]$Time<=j])
       }
     }
   }
   cfd$Prob<-cfd$Prob/PG
   n<-nrow(cfd)
   cfd$ni<-round((cfd$Prob-c(0,cfd$Prob[1:(n-1)]))*PG)

   pg<-seq(0,sum(cfd$ni)/PG,by=1/PG)


   return(list(cdf=cfd,cdfInt=data.frame(Time=quantileCFD(cfd,pg),
Prob=pg)))

}


# cfd columns Time and prob
quantileCFD<-function(cfd,prob=0.5)
{
  pos<-c()
  pos0<-c()
  for (i in 1:length(prob))
  {
    posA<-match(FALSE,cfd[["Prob"]]<prob[i], nomatch=NA)
    pos0A<-posA-1
    if (is.na(posA))	# Above maximum, use last and point of previous jump (to avoid flat areas)
    {
      posA<-which.max(cfd[["Prob"]])
      pos0A<-match(FALSE,cfd[["Prob"]]<cfd[posA,"Prob"], nomatch=NA)-1
    }
    if (pos0A==0)
    {
      posA<-match(FALSE,cfd[["Prob"]]<=prob[i], nomatch=NA)
      pos0A<-posA-1
    }
    pos<-c(pos,posA)
    pos0<-c(pos0,pos0A)
  }
  p<-cfd[pos,"Prob"]
  q<-cfd[pos0,"Prob"]
  y<-cfd[pos,"Time"]
  x<-cfd[pos0,"Time"]

  m<-x+(prob-q)*(y-x)/(p-q)
  return(m)
}

InterpolateFunction<-function(f,val=0, ymin=NA, ymax=NA)
{
  ret<-c()
  for (vl in val)
  {
    pos<-match(FALSE,f[,1]<vl, nomatch=NA)
    pos0<-pos-1
    if (is.na(pos))	# Above maximum, use last and point of previous jump (to avoid flat areas)
    {
      pos<-which.max(f[,1])
      pos0<-match(FALSE,f[,1]<f[pos,1], nomatch=NA)-1
    }
    if (pos==1)
    {
      pos0<-1
      pos<-match(FALSE,f[,1]>f[1,1], nomatch=NA)+1
    }
    q<-f[pos,2]
    y<-f[pos,1]
    p<-f[pos0,2]
    x<-f[pos0,1]

    m<-p+(vl-x)*(q-p)/(y-x)
    ret<-c(ret,m)
  }
  if (!is.na(ymin)) ret[ret<ymin]<-ymin
  if (!is.na(ymax)) ret[ret>ymax]<-ymax
  return(ret)
}

# Germination times and treatments: FTx and FXb computation
GerminationTT<-function(ddt)
{
  dg<-get.groups(ddt,ddt$Treatment)

  info<-data.frame(Treatment=NA,gamma=NA)
  cdf.list<-list()
  cdfInt.list<-list()
  for (i in 1:length(dg)) 
  {
    cc<-cdf(dg[[i]])
    cdf.list[[i]]<-cc$cdf
    cdfInt.list[[i]]<-cc$cdfInt
    info[i,]<-c(dg[[i]]$Treatment[1],max(cdf.list[[i]]$Prob))
  }
  names(cdf.list)<-names(dg)
  names(cdfInt.list)<-names(dg)
  return(list(Base=info,GerminationList=cdf.list,GerminationIntList=cdfInt.list))
}

# GmP: Parametric Gamma estimates or empirical ones (data.frame with Treatment - Gamma values to be tested)
# GerminationList: List of Germination Time distribution functions for the different treatments
ModelBasedTheta<-function(GmP,GerminationList)
{
  Theta<-c()
  ID<-c()
  k<-nrow(GmP)
  l<-length(GerminationList)
  for (i in 1:l)
  {
    for (j in 1:k)
    {
      if (max(GerminationList[[i]]$Prob)>=GmP$gamma[j])
      {
        ID<-c(ID,i)
        Theta<-c(Theta,((as.numeric(names(GerminationList)[i])-GmP$Treatment[j]))*quantileCFD(GerminationList[[i]],prob=GmP$gamma[j]))
      }
    }
  }
  ret<-data.frame(ID=ID,Theta=Theta)
  ret[ret$Theta>0,]
}

ModelBasedGerminationTime<-function(x,Theta,Time,FXb)
{
  ret<-data.frame(XbTime=x-(Theta/Time))
  if (is.data.frame(FXb)) ret$F<-InterpolateFunction(FXb,ret$XbTime,ymin=0, ymax=1)
  else ret$F<-predict(FXb,data.frame(Treatment=ret$XbTime))
  ret
}

# Transforms FTx function into FXb for a fixed Theta according to the model
FTxToFXb<-function(germ, Theta)
{
  tx<-c()
  ty<-c()
  id<-c()
  k<-nrow(germ$Base)
  for (i in 1:k)
  {
    nonull<-germ$GerminationList[[i]]$Time>0
    tx<-c(tx,germ$Base[i,1]-(Theta/germ$GerminationList[[i]][nonull,1]))
    ty<-c(ty,germ$GerminationList[[i]][nonull,2])
    id<-c(id,rep(i,sum(nonull)))
  }
  data.frame(ID=id,X=tx,F=ty)
}

MinCVTheta<-function(germ)
{
  k<-nrow(germ$Base)
  R<-c()
  X<-c()
  for (i in 1:k)
  {
    X<-c(X,rep(germ$Base[i,1],sum(germ$GerminationList[[i]]$ni)))
    OK<-germ$GerminationList[[i]]$ni>0
    R<-c(R,rep(1/germ$GerminationList[[i]]$Time[OK], times=round(germ$GerminationList[[i]]$ni[OK])))
  }
  n<-length(X)
  mR<- mean(R)*n/(n-1)
  mX<- mean(X)*n/(n-1)
  
  b<- 2*var(R)*(mR)^2 -2*(mR^2)*var(R)          
  c<- -2*cov(R,X)*mR^2-4*var(R)*mX*mR +2*mX*var(R)*mR+4*(mR^2)*cov(R,X)                         
  d<- 2*var(R)*(mX)^2+4*cov(R,X)*mX*mR-2*var(X)*mR^2-4*cov(R,X)*mX*mR                    
  e<- -2*cov(R,X)*(mX)^2+2*var(X)*mX*mR
  roots<-polyroot(c(e,d,c,b))
  roots<-Re(roots[(Im(roots)<1E-10)&(Re(roots)>0)])
  roots
}

MinVarQuantile<-function(germ,grid=NA,removeQ=NA)
{
  if (is.na(grid))
  {
    grid<-sort(unique(do.call(rbind, germ$GerminationList)$Prob))
    grid<-grid[-c(1,length(grid))]
  }
  
  n<-length(grid)
  k<-nrow(germ$Base)
  R<-matrix(nrow=k,ncol=n)
  X<-matrix(nrow=k,ncol=n)
  for (i in 1:k)
  {
    R[i,]<-1/quantileCFD(germ$GerminationList[[i]],grid)
    X[i,]<-germ[["Base"]][["Treatment"]][i]
    NOK<-(grid<min(germ$GerminationList[[i]]$Prob[-1]))|(grid>max(germ$GerminationList[[i]]$Prob))
    R[i,NOK]<-NA
  }

  # Remove elements of grid with no data
  nj<-colSums(!is.na(R))
  grid<-grid[nj>0]
  R<-R[,nj>0]
  X<-X[,nj>0]

  m<-colMeans(R,na.rm=TRUE)
  m2<-colMeans(R^2,na.rm=TRUE)
  MeanVarR<-mean(m2-m^2)
  MR<-rowMeans(R,na.rm=TRUE)
  CXR<-cov(germ[["Base"]][["Treatment"]],MR)*(k-1)/k
  Theta<-CXR/MeanVarR
  vr<-apply(X-Theta*R,2,var)

  if (!is.na(removeQ))
  {
    trm<-quantile(vr,removeQ)
    OK<-vr<trm

    m<-colMeans(R[,OK],na.rm=TRUE)
    m2<-colMeans(R[,OK]^2,na.rm=TRUE)
    MeanVarR<-mean(m2-m^2)
    MR<-rowMeans(R[,OK],na.rm=TRUE)
    CXR<-cov(germ[["Base"]][["Treatment"]],MR)*(k-1)/k
    Theta<-CXR/MeanVarR
    vr<-apply(X-Theta*R,2,var)
  }
  df<-as.data.frame(t(X-Theta*R))
  names(df)<-germ[["Base"]][["Treatment"]]
  return(list(Theta=Theta,info=cbind(data.frame(X=sort(rowMeans(df,na.rm=TRUE)),F=grid,Xorg=rowMeans(df,na.rm=TRUE),var=vr),df)))
}


MinCVTheta<-function(germ)
{
  k<-nrow(germ$Base)
  R<-c()
  X<-c()
  for (i in 1:k)
  {
    OK<-germ$GerminationIntList[[i]]$Prob>0
    X<-c(X,rep(germ$Base[i,1],length(germ$GerminationIntList[[i]]$Time[OK])))
    R<-c(R,1/germ$GerminationIntList[[i]]$Time[OK])
  }
  n<-length(X)
  mR<- mean(R)*n/(n-1)
  mX<- mean(X)*n/(n-1)
  
  b<- 2*var(R)*(mR)^2 -2*(mR^2)*var(R)          
  c<- -2*cov(R,X)*mR^2-4*var(R)*mX*mR +2*mX*var(R)*mR+4*(mR^2)*cov(R,X)                         
  d<- 2*var(R)*(mX)^2+4*cov(R,X)*mX*mR-2*var(X)*mR^2-4*cov(R,X)*mX*mR                    
  e<- -2*cov(R,X)*(mX)^2+2*var(X)*mX*mR
  roots<-polyroot(c(e,d,c,b))
  roots<-Re(roots[(Im(roots)<1E-10)&(Re(roots)>0)])
  roots
}


MinVarQuantile<-function(germ,grid=NA,removeQ=NA)
{
  if (is.na(grid))
  {
    grid<-sort(unique(do.call(rbind, germ$GerminationIntList)$Prob))
    grid<-grid[-c(1,length(grid))]
  }
  
  n<-length(grid)
  k<-nrow(germ$Base)
  R<-matrix(nrow=k,ncol=n)
  X<-matrix(nrow=k,ncol=n)
  for (i in 1:k)
  {
    R[i,]<-1/quantileCFD(germ$GerminationIntList[[i]],grid)
    X[i,]<-germ[["Base"]][["Treatment"]][i]
    NOK<-(grid<min(germ$GerminationIntList[[i]]$Prob[-1]))|(grid>max(germ$GerminationIntList[[i]]$Prob))
    R[i,NOK]<-NA
  }

  # Remove elements of grid with no data
  nj<-colSums(!is.na(R))
  grid<-grid[nj>0]
  R<-R[,nj>0]
  X<-X[,nj>0]

  m<-colMeans(R,na.rm=TRUE)
  m2<-colMeans(R^2,na.rm=TRUE)
  MeanVarR<-mean(m2-m^2)
  MR<-rowMeans(R,na.rm=TRUE)
  CXR<-cov(germ[["Base"]][["Treatment"]],MR)*(k-1)/k
  Theta<-CXR/MeanVarR
  vr<-apply(X-Theta*R,2,var)

  if (!is.na(removeQ))
  {
    repeat
    {
      trm<-quantile(vr,removeQ)
      OK<-vr<trm

      m<-colMeans(R[,OK],na.rm=TRUE)
      m2<-colMeans(R[,OK]^2,na.rm=TRUE)
      MeanVarR<-mean(m2-m^2)
      MR<-rowMeans(R[,OK],na.rm=TRUE)
      CXR<-cov(germ[["Base"]][["Treatment"]],MR)*(k-1)/k
      Theta<-CXR/MeanVarR
      vr<-apply(X-Theta*R,2,var)
      if (sum(OK!=(vr<quantile(vr,removeQ)),na.rm=TRUE)<1) break
    } 
  }
  df<-as.data.frame(t(X-Theta*R))
  names(df)<-germ[["Base"]][["Treatment"]]
  return(list(Theta=Theta,info=cbind(data.frame(X=sort(rowMeans(df,na.rm=TRUE)),F=grid,Xorg=rowMeans(df,na.rm=TRUE),var=vr),df)))
}




#dat <- read.table("Hídrico.csv", header = T, sep=";")
#head(dat)
#names(dat)<-c("Grouping","Population","Treatment","Dish","Time","G","PG") 
#dat$Gr<-factor(paste0(dat$Grouping,dat$Population))
#dat<-dat[,c(8,3:7)]
#names(dat)[1]<-"Grouping"
#summary(dat)
#dt<-get.groups(dat,dat$Grouping)


#dat <- read.table("Stephanie.csv", header = T, sep=";")
#head(dat)
#dt<-get.groups(dat,dat$Grouping)
#dt[[7]]<-NULL
#head(dt[[7]])


subset(read.csv("Water.csv"), 
       Experiment == "Dormancy") -> dat
dat[] <- lapply(dat, function(x) if(is.factor(x)) factor(x) else x) 
dat[, 1:6] -> dat
dt<-get.groups(dat,dat$Grouping)

#dat <- read.table("Combined no DMSO.csv", header = T, sep=",")
#head(dat)
#names(dat)<-c("Grouping","Population","Treatment","Dish","Time","G","PG") 
#dat$Gr<-factor(paste0(dat$Grouping,dat$Population))
#dat<-dat[,c(8,3:7)]
#names(dat)[1]<-"Grouping"
#summary(dat)
#dt<-get.groups(dat,dat$Grouping)


gr<-6
germ<-GerminationTT(dt[[gr]])

Th<-MinVarQuantile(germ,removeQ=0.9)
Theta<-Th$Theta
Theta<-MinCVTheta(germ)

plot(Th$info[,1:2])

k<-ncol(Th$info)-4
plot(unlist(Th$info[,-(1:4)]),rep(Th$info$F,k),col=rep(rainbow(k),each=nrow(Th$info)),main=paste0("JoinConversion -- ",as.character(dt[[gr]][1,1])))
legend("topleft",legend=germ$Base[,1],fill=rainbow(max(FXb2$ID)))
EG<-Th$info[,1:2]
lines(EG$X,EG$F,col="black",lwd=3)


EG<-Th$info[,1:2]

library(robustbase)
library(MASS)
library(Hmisc)  #For weighted measures: wtd.mean, wtd.var ...

estadisticos.media(t(Th$info[2,-(1:3)]))
huber(t(Th$info[2,-(1:3)]))
Mwgt(t(Th$info[2,-(1:3)]),cc=0.5,psi="Huber")
weighted.mean(t(Th$info[2,-(1:3)]),Mwgt(t(Th$info[2,-(1:3)]),cc=0.5,psi="huber")) 

apply(Th$info[,-(1:3)],1,huber)

pdf("HídricoThetaCVJoin.pdf",width=10,height=6)
for (gr in 1:length(dt))
{

# Computation of empirical distribution functions (GerminationTimes and BaseGermination)
germ<-GerminationTT(dt[[gr]])

plot(c(0,max(dt[[gr]]$Time)),c(0,1),type="n",main="Experimental Distribution Functions")
k<-length(germ$GerminationIntList)
cl<-rainbow(k)
for (pos in 1:k)
{
  lines(germ$GerminationIntList[[pos]],col=cl[pos],lwd=3)
}
legend("bottomright",fill=cl,legend=germ$Base[,1], pt.cex=1,cex=0.8)

# Parametric Gamma estimation
library(drc)
Gm <- drm(gamma ~ Treatment, data = germ$Base, fct = G.3())
GmP<-data.frame(Treatment=seq(min(germ$Base$Treatment),max(germ$Base$Treatment),length.out=20))
GmP$gamma<-predict(Gm,GmP)

GmP<-germ$Base

Th<-ModelBasedTheta(GmP,germ$GerminationList)
estadisticos.mediana(Th)
estadisticos.media(Th)

Theta<-mean(Th$Theta)

# Theta<-MinCVTheta(germ)

# Estimation of the CFDs by using the model - direct approach (use of estimated FXb directly)
# FXb<-data.frame(Treatment=germ$Base$Treatment)
# FXb$gamma<-predict(Gm,FXb)
FXb<-germ$Base
par.old<-par(mfrow=c(2,4),oma=c(0,0,4,0))
for (pos in 1:length(germ$GerminationList))
{
  tval<-seq(0,max(germ$GerminationList[[pos]][,1]),by=0.1)  # Grid for time points
  FTxE<-ModelBasedGerminationTime(FXb[pos,1],Theta,tval,FXb)  #Empirical Gamma
  FTxP<-ModelBasedGerminationTime(FXb[pos,1],Theta,tval,Gm) #Parametrised Gamma
  Empirical<-InterpolateFunction(germ$GerminationList[[pos]],tval)

  plot(c(min(tval),max(tval)),c(0,1),type="n",main=GmP[pos,1])
  lines(tval,Empirical,lty=2)
  points(germ$GerminationList[[pos]])
  OK<-FTxE$XbTime>=min(germ$Base$Treatment)
  lines(tval[OK],FTxE$F[OK],col="green")
  lines(tval[!OK],FTxE$F[!OK],col="red")
  OK<-FTxP$XbTime>=min(germ$Base$Treatment)
  lines(tval[OK],FTxP$F[OK],col="blue")
  lines(tval[!OK],FTxP$F[!OK],col="red")
  legend("bottomright",fill=c("green","blue","red"),legend=c("Empirical","Parametric","Extrapolation"), pt.cex=1,cex=0.8)
}
par(par.old)
title(paste0("Direct approach -- ",as.character(dt[[gr]][1,1])),  line = -2, outer = TRUE)


## Alternative approach 
# Transform FTx empirical data into FXb space and use combined information to estimate FXb
# A previous estimation of Theta is needed (option - look for the best one)
Theta<-MinCVTheta(germ)
FXb2<-FTxToFXb(germ, Theta)

# Transformation and computation of Theta optimaly
Th<-MinVarQuantile(germ,removeQ=0.9)
Theta2<-Th$Theta

# pdf(paste0(c("TBase_",as.character(dg[[1]][1,1]),"_",as.character(dg[[1]][1,2]),".pdf"),collapse=""),width=10,height=6)
library(drc)
OK<-FXb2$F>0
plot(FXb2$X[OK],FXb2$F[OK],col=rainbow(max(FXb2$ID))[FXb2$ID[OK]],main=paste0("FTxToFXb (MinVarQuantile) -- ",as.character(dt[[gr]][1,1])))
lines(germ$Base)
legend("topleft",legend=germ$Base[,1],fill=rainbow(max(FXb2$ID)))

fm <- drm(F ~ X, data = FXb2, fct = G.3())
grid<-data.frame(X=seq(min(FXb2$X),max(FXb2$X),by=0.1))
lines(grid$X,predict(fm,grid))

k<-ncol(Th$info)-4
plot(unlist(Th$info[,-(1:4)]),rep(Th$info$F,k),col=rep(rainbow(k),each=nrow(Th$info)),main=paste0("JoinConversion -- ",as.character(dt[[gr]][1,1])))
legend("topleft",legend=germ$Base[,1],fill=rainbow(max(FXb2$ID)))
EG<-Th$info[,1:2]
lines(EG$X,EG$F,col="black",lwd=3)

FXb<-germ$Base

par.old<-par(mfrow=c(2,4),oma=c(0,0,4,0))
for (pos in 1:length(germ$GerminationList))
{
  tval<-seq(0,max(germ$GerminationList[[pos]][,1]),by=0.1)  # Grid for time points
  FTxE<-ModelBasedGerminationTime(FXb[pos,1],Theta2,tval,EG)  #Empirical Gamma (join conversion)
  FTxP<-ModelBasedGerminationTime(FXb[pos,1],Theta,tval,fm) #Parametrised Gamma
  Empirical<-InterpolateFunction(germ$GerminationList[[pos]],tval,ymax=1)

  plot(c(min(tval),max(tval)),c(0,1),type="n",main=GmP[pos,1])
  lines(tval,Empirical,lty=2)
  points(germ$GerminationList[[pos]])
  lines(tval,FTxE$F,col="green")

  
  lines(tval,FTxP$F,col="blue")
  legend("bottomright",fill=c("green","blue"),legend=c("Empirical","Parametric"), pt.cex=1,cex=0.8)
}
par(par.old)
title(paste0("Alternative approaches (Th Est + Parametric vs join estimation) -- ",as.character(dt[[gr]][1,1])),  line = -2, outer = TRUE)

## Alternative approach - fitting theta


}
dev.off()

### Two population case
dr<-dt[[gr]]
dr[["Population"]]<-factor(dr[["Population"]])
dgg<-get.groups(dr,dr$Population)
dg<-get.groups(dgg[[1]],dgg[[1]]$Treatment)
dg<-get.groups(dgg[[2]],dgg[[2]]$Treatment)

### Single population case
dg<-get.groups(dt[[gr]],dt[[gr]]$Treatment)

info<-data.frame(Treatment=NA,gamma=NA)
cdf.list<-list()
for (i in 1:length(dg)) 
{
  cdf.list[[i]]<-cdf(dg[[i]])
  info[i,]<-c(dg[[i]]$Treatment[1],max(cdf.list[[i]]$Prob))
}

## Tiempo de crecimiento
# t.crec<--1.3
# for (i in 1:length(cdf.list)) 
# {
#   cdf.list[[i]]<-cdf.list[[i]][cdf.list[[i]][,1]>t.crec,]
#   cdf.list[[i]][,1]<-cdf.list[[i]][,1]-t.crec
# }

fm <- drm(gamma ~ Treatment, data = info, fct = G.3())
grid<-seq(min(info$Treatment),max(info$Treatment),length.out=20)
plot(info)
lines(grid,predict(fm,data.frame(Treatment=grid)))

or<-order(info$gamma)	### Reorder, has sense if no big differences as an improvement of the estimation (check Dette ordering) 

or<-1:(length(info$gamma))

n.info<-info
info
n.info$gamma<-n.info$gamma[or]
n.cdf.list<-cdf.list[or]

n.info$gamma<-sort(n.info$gamma)

# n.info<-n.info[-c(2,3,4),]
# n.cdf.list[[3]]<-NULL
# n.cdf.list[[4]]<-NULL
# n.cdf.list[[2]]<-NULL
## Assume Gaussian, estimate n.info$gamma (quantiles of Gaussian) and addapt maximums of distributions
# cfd<-data.frame(Time=n.info[,1],Prob=n.info[,2])
# s1<-quantileCFD(cfd,prob=pnorm(-1))
# m<-quantileCFD(cfd,prob=0.5)
# s2<-quantileCFD(cfd,prob=pnorm(1))
# n.info$gamma<-pnorm(n.info[,1],m,(s2-s1)/2)
# for (i in 1:k)
# {
#   n.cdf.list[[i]]$Prob<-n.info$gamma[i]*n.cdf.list[[i]]$Prob/max(n.cdf.list[[i]]$Prob)
# }

# Estimation of inverse of cdf at each gamma
p<-1
Theta<-c()
ID<-c()
k<-nrow(n.info)
for (i in 2:k)
{
  for (j in 1:(i-1))
  {
    ID<-c(ID,i)
    Theta<-c(Theta,((n.info$Treatment[i]-n.info$Treatment[j])^(p))*quantileCFD(n.cdf.list[[i]],prob=n.info$gamma[j]))
  }
}
data.frame(ID,Theta)

estadisticos.media(Theta)
Th2<-by.groups(Theta,ID,estadisticos.media)
Th<-mean(Th2[,1])

#Pintar para cada uno lo estimado con el modelo (Theta y n.info) frente a lo obtenido 
Th<-mean(Theta)
Th<-median(Theta)

pdf(paste0(c("Modelo_",as.character(dg[[1]][1,1]),"_",as.character(dg[[1]][1,2]),".pdf"),collapse=""),width=10,height=6)
par.old<-par(mfrow=c(2,4),oma=c(0,0,4,0))
for (pos in 1:length(n.cdf.list))
{
tval<-seq(0,max(n.cdf.list[[pos]][,1]),by=0.1)
FTx<-pmax(0,InterpolateFunction(n.info,n.info[pos,1]-(Th/tval)^(1/p)))
Empirical<-InterpolateFunction(n.cdf.list[[pos]],tval)

plot(c(min(tval),max(tval)),c(0,1),type="n",main=n.info[pos,1])
lines(tval,Empirical,lty=2)
points(n.cdf.list[[pos]])
lines(tval,FTx)
}
par(par.old)
title(paste0("Direct approach -- ",as.character(dg[[1]][1,1])), side = 3, line = -2, outer = TRUE)



# Th<-0.7
# Th<-c(1,1,1,1,1,1,1,1)
Th<-mean(Theta)
Th<-median(Theta)
Th<-rep(Th,k)
# Th<-rep(1,k)
#Th<-Th2
# Th[3:5]<-Th[3:5]*1.8
# Th[5]<-0.2
tx<-c()
ty<-c()
id<-c()
for (i in 1:k)
{
  for (j in 1:nrow(n.cdf.list[[i]]))
  {
    if (n.cdf.list[[i]][j,1]>0)
    {
      tx<-c(tx,n.info[i,1]-(Th[i]/(n.cdf.list[[i]][j,1]))^(1/p))
      ty<-c(ty,n.cdf.list[[i]][j,2])
      id<-c(id,i)
    }
  }
}

# pdf(paste0(c("TBase_",as.character(dg[[1]][1,1]),"_",as.character(dg[[1]][1,2]),".pdf"),collapse=""),width=10,height=6)
plot(tx,ty,col=rainbow(k)[id],main=paste0("Joint graph -- ",as.character(dg[[1]][1,1])))
lines(n.info)
#lines(ksmooth(tx, ty, bandwidth = 0.3, x.points=tx))
legend("topleft",legend=n.info[,1],fill=rainbow(k))
library(drc)
or<-order(tx)
tx<-tx[or]
ty<-ty[or]
id<-id[or]
fm <- drm(y ~ x, data = data.frame(x=tx,y=ty), fct = G.3())
# fm <- drm(y ~ x, data = data.frame(x=tx,y=ty), fct = W1.3())
lines(tx,predict(fm))
# LL.2LL.3LL.3uLL.4 LL.5 W1.2 W1.3 W1.4 ?W2.2 W2.3 W2.4 BC.4 BC.5 LL2.2 LL2.3 LL2.3u LL2.4 LL2.5 AR.2 AR.3 MM.2 MM.3 G.2, G.3, G.3u G.4
dev.off()

pdf(paste0(c("Modelo2_",as.character(dg[[1]][1,1]),"_",as.character(dg[[1]][1,2]),".pdf"),collapse=""),width=10,height=6)
par.old<-par(mfrow=c(2,4))
for (pos in 1:length(n.cdf.list))
{
tval<-seq(0,max(n.cdf.list[[pos]][,1]),by=0.1)
FTx<-pmax(0,InterpolateFunction(n.info,n.info[pos,1]-Th[pos]/tval))
FTx<-predict(fm,data.frame(x=n.info[pos,1]-Th[pos]/tval))
Empirical<-InterpolateFunction(n.cdf.list[[pos]],tval)
plot(tval,pmax(FTx,Empirical),type="n",main=n.info[pos,1],ylim=c(0,1))
lines(tval,Empirical,lty=2)
points(n.cdf.list[[pos]])
lines(tval,FTx)
}
par(par.old)
dev.off()



ksmooth(tx, ty, bandwidth = 0.5, x.points=tx)









