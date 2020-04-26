setwd ("C:/EFP/Trabayu/Activos/2018 FICYT3/GT2 - R tool/Gil/Simulación")
source("/home/gil/Documentos/Actual/R/funciones.R")

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
  data[1,c("Time","G")]<-c(0,0)	# Force a row with Time=0
  f<-factor(data$Time)	# For arbitrary times
  
  data$Time<-as.numeric(f)
  nTime<-nlevels(f)
  cfd<-data.frame(Time=as.numeric(levels(f)), Prob=rep(0,nTime))
  
  sm<-get.groups(data,data$Dish)
  PG<-0
  for (i in 1:length(sm))
  {
    PG<-PG+max(sm[[i]]$PG)
    for (j in 1:nTime)
    {
      if (sum(sm[[i]]$Time<=j)>0) cfd$Prob[j]<-cfd$Prob[j]+max(sm[[i]]$G[sm[[i]]$Time<=j])
    }
  }
  cfd$Prob<-cfd$Prob/PG
  n<-nrow(cfd)
  cfd$ni<-round((cfd$Prob-c(0,cfd$Prob[1:(n-1)]))*PG)
  
  pg<-seq(0,sum(cfd$ni)/PG,by=1/PG)
  
  
  return(list(cdf=cfd,cdfInt=data.frame(Time=quantileCFD(cfd,pmax(0,pg-(0.5/PG))), Prob=pg)))
}

# cfd columns Time and prob
quantileCFD<-function(cfd,prob=0.5,extrapolation=TRUE)
{
  pos<-c()
  pos0<-c()
  nas<-rep(FALSE,length(prob))
  for (i in 1:length(prob))
  {
    posA<-match(FALSE,cfd[["Prob"]]<prob[i], nomatch=NA)
    pos0A<-posA-1
    if (is.na(posA))	# Above maximum, use last and point of previous jump (to avoid flat areas)
    {
      nas[i]<-!extrapolation
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
  m[nas]<-NA
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



MinVarQuantile<-function(germ,grid=NA,removeQ=NA,max.iter=4)
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

  # Remove elements of grid with no data or with only one data point
  nj<-colSums(!is.na(R))
  Ugrid<-grid[nj>1]
  UR<-R[,nj>1]
  UX<-X[,nj>1]
  nj<-nj[nj>1]
  
  # Check if rows with no data
  ni<-rowSums(!is.na(UR))
  UR<-UR[ni>0,] 
  UX<-UX[ni>0,]
  k<-sum(ni>0)
  
  m<-colMeans(UR,na.rm=TRUE)
  m2<-colMeans(UR^2,na.rm=TRUE)
#   MeanVarR<-mean(m2-m^2)
  MeanVarR<-weighted.mean(m2-m^2,nj)
  MR<-rowMeans(UR,na.rm=TRUE)
  CXR<-cov(germ[["Base"]][["Treatment"]][ni>0],MR)*(k-1)/k
  Theta<-CXR/MeanVarR
  vr<-apply(UX-Theta*UR,2,var)

  if (!is.na(removeQ))
  {
    iter<-0
    repeat
    {
      trm<-quantile(vr,removeQ)
      OK<-vr<trm
      RR<-UR[,OK]
      RX<-UX[,OK]

      # Check if rows with no data
      nj<-colSums(!is.na(RR))
      ni<-rowSums(!is.na(RR))
      RR<-RR[ni>0,] 
      RX<-RX[ni>0,]
      k<-sum(ni>0)
      
      m<-colMeans(RR,na.rm=TRUE)
      m2<-colMeans(RR^2,na.rm=TRUE)
#       MeanVarR<-mean(m2-m^2)
      MeanVarR<-weighted.mean(m2-m^2,nj)     
      MR<-rowMeans(RR,na.rm=TRUE)
      CXR<-cov(germ[["Base"]][["Treatment"]][ni>0],MR)*(k-1)/k
      Theta<-CXR/MeanVarR
      vr<-apply(UX-Theta*UR,2,var)
      iter<-iter+1
      if (sum(OK!=(vr<quantile(vr,removeQ)),na.rm=TRUE)<1) break
      if (iter>=max.iter) break
    } 
  }
  df<-as.data.frame(t(X-Theta*R))
  vr<-apply(X-Theta*R,2,var)
  names(df)<-germ[["Base"]][["Treatment"]]
  return(list(Theta=Theta,info=cbind(data.frame(X=sort(rowMeans(df,na.rm=TRUE)),F=grid,Xorg=rowMeans(df,na.rm=TRUE),var=vr),df)))
}


# 
# 
# dat <- read.table("Hídrico.csv", header = T, sep=";")
# names(dat)<-c("Grouping","Population","Treatment","Dish","Time","G","PG") 
# dat$Gr<-factor(paste0(dat$Grouping,dat$Population))
# dat<-dat[,c(8,3:7)]
# names(dat)[1]<-"Grouping"
# 
# Time<-sort(unique(round(dat[dat$Grouping=="AETR029AT","Time"])))  # Observation times
# 
# # pdf("SimulationThetaCVJoin.pdf",width=10,height=6)
# # for (kkk in 1:10)
# # {

SimulateGerminatio<-function(Treatments=c(-1.6,-1,-0.8,-0.6,-0.4,-0.2,0), Ptheta=33, Time=c(7,14,18,21,28,33,37,41,47,52,57,62,76,89,101,111,126,135,151,197,341,532))
{
  x<-Treatments
  k<-length(x)

  shape<-2
  scale<-1
  Xb<- -rweibull(k*100,shape,scale) # +0.4 para inducir dormición, -0.4 para que varios tratamientos lleguen al 100%


  data<-data.frame(Treatment=rep(x,each=100),Dish=NA,Xb=Xb,G=NA,PG=25,Tx=NA)        
  data$Dish<-factor(LETTERS[rep(rep(1:4,each=25),times=k)])
  data$Tx<-Ptheta/(data$Treatment-Xb)
  data$Tx[data$Tx<0]<-NA

  #Censored times - common censoring for all treatments
  data$cTx<-Time[as.numeric(cut(data$Tx,c(0,Time)))]  
  data$cXb<-data$Treatment-Ptheta/data$cTx
  
  # Frequencies computation - Transform data to the usual input
  pos<-1
  gr<-get.groups(data,data$Treatment)
  dat<-data.frame() 
  for (Tr in gr)
  {
    gr2<-get.groups(Tr,Tr$Dish)
    for (Ds in gr2)
    {
      fr<-freq(cut(Ds$Tx,c(0,Time)),type="ni",with.na=FALSE)
      nr<-length(fr)
      dat<-rbind(dat,data.frame(Treatment=Ds$Treatment[1],Dish=Ds$Dish[1],Time=Time,G=cumsum(as.numeric(fr)),PG=25))   
    }
  }
  sim=dat
  
  # From input (dat) get interpolation of censored times!
  germ<-GerminationTT(dat)
  
  # Compute the interpolated Tx for each treatment (just for comparative purpouses) and its transformation
  gr<-get.groups(data,data$Treatment)
  data<-data.frame()
  for (i in 1:length(gr))
  {
    nm<-names(germ$GerminationIntList)[[i]]
    dgerm<-germ$GerminationIntList[[nm]][-1,]
    if (nrow(dgerm)< nrow(gr[[nm]])) dgerm[ (nrow(dgerm)+1):nrow(gr[[nm]]),]<-NA
    data<-rbind(data,cbind(gr[[nm]][order(gr[[nm]]$Tx),],data.frame(iTx=dgerm$Time,iXb=gr[[nm]]$Treatment-Ptheta/dgerm$Time,Prob=dgerm$Prob)) )
  }
  
  grid<-seq(minimo(data$Xb),maximo(data$Xb),length.out=1000) 
  XbTheoretical<-data.frame(Xb=seq(minimo(data$Xb),maximo(data$Xb),length.out=1000), Prob=1-pweibull(-grid,2,1))
  md<-by.groups(data$iXb,data$Prob,media) # Mean of interpolated Xb by levels
  md[,1]<-isoreg(md[,1])$yf
  md$Prob<-as.numeric(row.names(md))
  names(md)[1]<-"Xb"

  list(sim=sim, dtorig=data, XbTheoretical=XbTheoretical, XbEstimated=md)
}


# Simulación concreta con gráficos
sm<-SimulateGerminatio(Ptheta=33,Time=seq(0.01,1000,by=1))  
dat<-sm$sim
data<-sm$dtorig

iTr<-as.numeric(factor(data$Treatment))  
plot(data$iXb,data$Prob,pch="*",col=rainbow(7)[iTr], main="Censored data: interpolated and transformed back with true theta")  # Plot of different true transformed interpolated distribution functions
legend("topleft",fill=rainbow(7),legend=levels(factor(data$Treatment)))
lines(ecdf(data$Xb),pch=".",lwd=2,col="brown")  # Empirical based on true data (best possible estimation)
lines(sm$XbTheoretical,lwd=2) # Theoretical distribution function
lines(sm$XbEstimated,col="blue",lwd=2)  # Empirical estimated by levelwise mean of transformed (using true theta) interpolated data
legend("top",lwd=2,col=c("black","brown","blue"),legend=c("Theoretical","Empirical - Original","Empirical - Estimated"))
  
  
  
# SIMULACION  
# Start of estimation procedure  
# From input (dat) get interpolation of censored times!
aa<-replicate(100,{
 sm<-SimulateGerminatio(Ptheta=1,Time=seq(0.01,20,by=0.05))  
 dat<-sm$sim
 data<-sm$dtorig
 germ<-GerminationTT(dat)
 
 GmP<-germ$Base
 c(media=media(ModelBasedTheta(GmP,germ$GerminationIntList))$Theta, minCV=MinCVTheta(germ), minVQ=MinVarQuantile(germ,removeQ=NA)$Theta)
})

  
estadisticos.media(as.data.frame(t(aa)))
estadisticos.mediana(as.data.frame(t(aa)))
cuantiles(as.data.frame(t(aa)),prob=seq(0,1,by=0.1))










MinCVTheta(germ)

Theta<-media(ModelBasedTheta(GmP,germ$GerminationIntList))$Theta
FXb2<-FTxToFXb(germ, Theta)

# Transformation and computation of Theta optimaly
#Th<-MinVarQuantile(germ,removeQ=0.9)
Th<-
Theta2<-Th$Theta

estadisticos.media(Th)
Th<-ModelBasedTheta(GmP,germ$GerminationIntList)
estadisticos.media(Th)

print(estadisticos.media(Th))

head(data)


mecdf<-function(x)
{
  n <- length(x)
  x <- sort(x)
  if (length(x) < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/n, 
      method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}



plot(mecdf(data$cTx[(data$Treatment==0)&(data$Dish=="D")]),col="blue")
plot(mecdf(data$Tx[(data$Treatment==0)&(data$Dish=="D")]),add=T)
plot(mecdf(data$cTx[(data$Treatment==0)&(data$Dish=="D")]),col="blue",add=T)

plot(mecdf(data$iTx[(data$Treatment==0)&(data$Dish=="D")]),col="blue")
plot(mecdf(data$Tx[(data$Treatment==0)&(data$Dish=="D")]),add=T)
plot(mecdf(data$iTx[(data$Treatment==0)&(data$Dish=="D")]),col="blue",add=T)


plot(mecdf(data$cXb[(data$Treatment==0)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment==0)]),add=T)
plot(mecdf(data$cXb[(data$Treatment==0)]),col="blue",add=T)

plot(mecdf(data$iXb[(data$Treatment==0)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment==0)]),add=T)
plot(mecdf(data$iXb[(data$Treatment==0)]),col="blue",add=T)

plot(mecdf(data$cXb[(data$Treatment== -0.2)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -0.2)]),add=T)
plot(mecdf(data$cXb[(data$Treatment== -0.2)]),col="blue",add=T)

plot(mecdf(data$iXb[(data$Treatment== -0.2)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -0.2)]),add=T)
plot(mecdf(data$iXb[(data$Treatment== -0.2)]),col="blue",add=T)

plot(mecdf(data$cXb[(data$Treatment== -0.4)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -0.4)]),add=T)
plot(mecdf(data$cXb[(data$Treatment== -0.4)]),col="blue",add=T)

plot(mecdf(data$iXb[(data$Treatment== -0.4)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -0.4)]),add=T)
plot(mecdf(data$iXb[(data$Treatment== -0.4)]),col="blue",add=T)


plot(mecdf(data$cTx[(data$Treatment== -0.4)]),col="blue")
plot(mecdf(data$Tx[(data$Treatment== -0.4)]),add=T)
plot(mecdf(data$cTx[(data$Treatment== -0.4)]),col="blue",add=T)

plot(mecdf(data$cXb[(data$Treatment== -0.6)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -0.6)]),add=T)
plot(mecdf(data$cXb[(data$Treatment== -0.6)]),col="blue",add=T)

plot(mecdf(data$iXb[(data$Treatment== -0.6)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -0.6)]),add=T)
plot(mecdf(data$iXb[(data$Treatment== -0.6)]),col="blue",add=T)


plot(mecdf(data$cXb[(data$Treatment== -1.6)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -1.6)]),add=T)
plot(mecdf(data$cXb[(data$Treatment== -1.6)]),col="blue",add=T)

plot(mecdf(data$iXb[(data$Treatment== -1.6)]),col="blue")
plot(mecdf(data$Xb[(data$Treatment== -1.6)]),add=T)
plot(mecdf(data$iXb[(data$Treatment== -1.6)]),col="blue",add=T)


plot(mecdf(data$iXb),col="blue")
plot(mecdf(data$Xb),add=T)
plot(mecdf(data$iXb),col="blue",add=T)

data[data$iXb< -2,]
estadisticos.media(data$iXb-data$Xb)
estadisticos.media(data$cXb-data$Xb)

data[(data$Treatment== -0.4),]


plot(Ds$Xb,Ds$cXb)

germ<-GerminationTT(dat)


plot(ecdf(data$cTx[(data$Treatment==0)]),col="blue")
plot(ecdf(data$Tx[(data$Treatment==0)]),add=T)
plot(ecdf(data$cTx[(data$Treatment==0)]),col="blue",add=T)

plot(ecdf(data$Tx))
plot(ecdf(data$cTx),add=T,col="blue")


lines(germ$GerminationList[['0']][1:2])


head(dat)
summary(dat)
dt<-list()
dt[[1]]<-dat


# pdf("SimulationThetaCVJoin.pdf",width=10,height=6)
for (gr in 1:1)
{

histograma(Xb,main=paste0("Xb distribution - Simulated from -Weibull(",as.character(shape),",",as.character(scale),")"))

plot(ecdf(Xb),main=paste0("Empirical Xb distribution - Simulated from -Weibull(",as.character(shape),",",as.character(scale),")"))

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
Th<-ModelBasedTheta(GmP,germ$GerminationIntList)
estadisticos.mediana(Th)
print(estadisticos.media(Th))

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
title(paste0("Direct approach -- Theta=",as.character(Theta)),  line = -2, outer = TRUE)


## Alternative approach 
# Transform FTx empirical data into FXb space and use combined information to estimate FXb
# A previous estimation of Theta is needed (option - look for the best one)
Theta<-MinCVTheta(germ)
FXb2<-FTxToFXb(germ, Theta)

# Transformation and computation of Theta optimaly
#Th<-MinVarQuantile(germ,removeQ=0.9)
Th<-MinVarQuantile(germ,removeQ=NA)
Theta2<-Th$Theta

# pdf(paste0(c("TBase_",as.character(dg[[1]][1,1]),"_",as.character(dg[[1]][1,2]),".pdf"),collapse=""),width=10,height=6)
library(drc)
OK<-FXb2$F>0
plot(FXb2$X[OK],FXb2$F[OK],col=rainbow(max(FXb2$ID))[FXb2$ID[OK]],main=paste0("FTxToFXb (MinCVTheta) -- Theta=",as.character(Theta)))
lines(germ$Base)
legend("topleft",legend=germ$Base[,1],fill=rainbow(max(FXb2$ID)))

fm <- drm(F ~ X, data = FXb2, fct = G.3())
grid<-data.frame(X=seq(min(FXb2$X),max(FXb2$X),by=0.1))
lines(grid$X,predict(fm,grid))

k<-ncol(Th$info)-4
plot(unlist(Th$info[,-(1:4)]),rep(Th$info$F,k),col=rep(rainbow(k),each=nrow(Th$info)),main=paste0("JoinConversion -- Theta=",as.character(Theta2)))
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












seed.germinatio<-function(dt,min.points=2)
{
  dg<-get.groups(dt,dt$Treatment)
  med<-data.frame(Treatment=NA,cdf=NA)
  for (i in dt) 
  {
    m<-quantileIntervalCensoredData(i,prob=prob) 
    med<-rbind(med,cbind(i$Treatment[1],1/m$q))
  }
  names(med)<-c("Treatment","Rate")
  keyTreatments(med,min.points=min.points)
}



sg<-list()
for(i  in seq_along(dt))
{
  sg[[names(dt)[[i]]]]<-seed.germinatio(dt[[i]],min.points=nrow,prob=0.30)
}
names(sg)
plot(sg[["MOMO"]])





# data: data.frame with columns Time, G, PG
quantileIntervalCensoredData<-function(data,prob=0.5)
{
  f<-factor(data$Time)	# For arbitrary times 
  data$Time<-as.numeric(f)-1
  Time<-sort(unique(c(0,data$Time)))
  Time
  cfd<-data.frame(Time=Time, Prob=rep(0,length(Time)), Tot=rep(0,length(Time)))
  for (i in 1:nrow(data))
  {
    cfd[data[i,"Time"]+1,"Prob"]<-cfd[data[i,"Time"]+1,"Prob"]+data[i,"G"]
    cfd[data[i,"Time"]+1,"Tot"]<-cfd[data[i,"Time"]+1,"Tot"]+data[i,"PG"]
  }
  cfd[["Prob"]]<-cfd[["Prob"]]/cfd[["Tot"]]
  cfd[1,"Prob"]<-0
  cfd$Time<-as.numeric(levels(f)) # Restore original time values  

  if (max(data$G)==0) return(list(q=Inf,cfd=cfd))
  
  pos<-match(FALSE,cfd[["Prob"]]<prob, nomatch=NA)
  pos0<-pos-1
  if (is.na(pos))	# Above maximum, use last and point of previous jump (to avoid flat areas)
  {
    pos<-which.max(cfd[["Prob"]])
    pos0<-match(FALSE,cfd[["Prob"]]<cfd[pos,"Prob"], nomatch=NA)-1
  }
  p<-cfd[pos,"Prob"]
  q<-cfd[pos0,"Prob"]
  y<-cfd[pos,"Time"]
  x<-cfd[pos0,"Time"]

  m<-x+(prob-q)*(y-x)/(p-q)
  return(list(q=m,cfd=cfd))
}


keyTreatments<-function(med,min.points=2)
{
  if (is.function(min.points)) min.points<-min.points(med)
  if (nrow(med)==0) # Nothing to do
  {
    return(classlist(data=med,keyTr=data.frame(TB=NA,TP=NA,TE=NA,row.names = NULL),R1=NA,R2=NA,class="germinatio"))
  }

  # Sort in ascending treatments order
  ord<-order(med$Treatment)
  med<-med[ord,]
  md<-med
  
  # Remove extra 0's on the tails (preserve 1 cero at most in each tail)
  p0<-which.max(med$Rate>0)-1
  if (p0==0) p0<-1
  n<-nrow(med)
  p1<-which.max(med$Rate[n:1]>0)-1
  if (p1==0) p1<-1
  p1<-(n-p1)+1
  med<-med[p0:p1,]

  # Determine the optimum peak (minimum total sum of squared residual)
  if (max(med$Rate)<1E-10) # No germination for any treatment
  {
    return(classlist(data=md,keyTr=data.frame(TB=max(med$Treatment),TP=NA,TE=NA,row.names = NULL),R1=NULL,R2=NULL,class="germinatio"))
  }
  
  if (nrow(med)==1) # For this element Rate is greater than 0	
  {
    return(classlist(data=md,keyTr=data.frame(TB=NA,TP=med$Treatment[1],TE=NA,row.names = NULL),R1=NULL,R2=NULL,class="germinatio"))
  }

  
  n<-nrow(med)
  LT<-lm(Rate~Treatment,data=med)
  csum<-sum(LT$residuals^2)
  if (LT$coeff[2]>0) {peak<-n} else {peak<-1}
  
  LD1<-NULL
  LD2<-NULL
  if (min.points<n)
  {
    for (i in min.points:(n-min.points+1))
    {
      L1<-lm(Rate~Treatment,data=med[1:i,])
      L2<-lm(Rate~Treatment,data=med[i:n,])
      
      if (L1$coeff[2]<=0) L1$residuals<-rep(Inf,length(L1$residuals))
      if (L2$coeff[2]>=0) L2$residuals<-rep(Inf,length(L2$residuals))
      val<-sum(c(L1$residuals,L2$residuals)^2)
      if (val<csum)
      {
	csum<-val
	peak<-i
	LD1<-L1
	LD2<-L2
      }
    }
  }

  if (is.null(LD1))
  {
    if (peak==1) 
    {
      return(classlist(data=md,keyTr=data.frame(TB=NA,TP=med$Treatment[peak],TE=-LT$coeff[1]/LT$coeff[2],row.names = NULL),R1=NULL,R2=LT,class="germinatio"))
    }
    else
    {
      return(classlist(data=md,keyTr=data.frame(TB=-LT$coeff[1]/LT$coeff[2],TP=med$Treatment[peak],TE=NA,row.names = NULL),R1=LT,R2=NULL,class="germinatio"))
    }
  }
  else
  {
    return(classlist(data=md,keyTr=data.frame(TB=-LD1$coeff[1]/LD1$coeff[2],TP=(LD1$coeff[1]-LD2$coeff[1])/(LD2$coeff[2]-LD1$coeff[2]),TE=-LD2$coeff[1]/LD2$coeff[2],row.names = NULL),R1=LD1,R2=LD2,class="germinatio"))
  }
}

seed.germinatio<-function(datA,prob=0.5,min.points=2)
{
  dt<-get.groups(datA,datA$Treatment)
  med<-data.frame()
  for (i in dt) 
  {
    m<-quantileIntervalCensoredData(i,prob=prob) 
    med<-rbind(med,cbind(i$Treatment[1],1/m$q))
  }
  names(med)<-c("Treatment","Rate")
  keyTreatments(med,min.points=min.points)
}

plot.germinatio<-function(kTr,...)
{
  ymax=max(kTr$data["Rate"])+0.001
  if (!is.null(kTr$R1)) ymax=max(predict(kTr$R1,list(Treatment=kTr$keyTr[1,2])),kTr$data$Rate)
  if (!is.null(kTr$R2)) ymax=max(predict(kTr$R2,list(Treatment=kTr$keyTr[1,2])),kTr$data$Rate)

  xmin<-min(kTr$data["Treatment"],kTr$keyTr,na.rm=TRUE)
  xmax<-max(kTr$data["Treatment"],kTr$keyTr,na.rm=TRUE)
  marg<-(xmax-xmin)*0.1
  plot(kTr$data[["Treatment"]],kTr$data[["Rate"]],xlim=c(xmin-marg,xmax+marg),ylim=c(0,ymax*1.1),xlab="Treatment",ylab="Rate")
  if (!is.null(kTr$R1)) 
  {
    ymax<-predict(kTr$R1,list(Treatment=kTr$keyTr[1,2]))
    points(kTr$R1$model[,2],kTr$R1$model[,1],col="blue",pch=3)
    lines(c(kTr$keyTr[1,1],kTr$keyTr[1,2]),c(0,ymax))
    points(c(kTr$keyTr[1,1],kTr$keyTr[1,2]),c(0,ymax),pch=20)
  }
  if (!is.null(kTr$R2))
  {
    ymax<-predict(kTr$R2,list(Treatment=kTr$keyTr[1,2]))
    points(kTr$R2$model[,2],kTr$R2$model[,1],col="red",pch=4)
    lines(c(kTr$keyTr[1,2],kTr$keyTr[1,3]),c(ymax,0))  
    points(c(kTr$keyTr[1,2],kTr$keyTr[1,3]),c(ymax,0),pch=20)
  }
}


dat <- read.table("S1 Example Data.txt", header = T)
head(dat)
summary(dat)

dt<-get.groups(dat,dat$Grouping)
sg<-list()
for(i  in seq_along(dt))
{
  sg[[names(dt)[[i]]]]<-seed.germinatio(dt[[i]])
}

plot(sg[["A"]])

sum(sg[["A"]]$R1$residuals^2)

x<-dt[["A"]]
m<-quantileIntervalCensoredData(x[x$Treatment==16.25,])
m$q
plot(m$cfd$Time,m$cfd$Prob)


dat <- read.table("Malva.csv", header = T, sep=";")
m<-quantileIntervalCensoredData(dat)
m$q
plot(m$cfd$Time,m$cfd$Prob)

dat <- read.table("Arroxos.csv", header = T, sep=";")
m<-quantileIntervalCensoredData(dat)
m$q
plot(m$cfd$Time,m$cfd$Prob)

dat <- read.table("Stephanie.csv", header = T, sep=";")
head(dat)

dt<-get.groups(dat,dat$Grouping)
dt[[7]]<-NULL

dt[[7]][1:200,]

sg<-list()
for(i  in seq_along(dt))
{
  sg[[names(dt)[[i]]]]<-seed.germinatio(dt[[i]],min.points=nrow,prob=0.30)
}
names(sg)
plot(sg[["MOMO"]])


x<-dt[["MOMO"]]
x[x$Treatment==-0.4,]

m<-quantileIntervalCensoredData(x[x$Treatment==-0.6,])
plot(m$cfd$Time,m$cfd$Prob)


x<-sg[["MOMO"]]$data

class(sg[["ANCO"]])

data<-dt[[8]]
seed.germinatio(data,min.points=nrow)
head(data)

datA<-data

points(k$R2$model[,1],k$R2$model[,2],pch=3,col="red")

points(c(25,26),0.03,pch=3,col="red")

k$keyTr[1,3]
predict(k$R1,list(Treatment=k$keyTr[1,2]))
k[["keyTr"]]

plot(med[["Treatment"]],med[["Rate"]])


m<-quantileIntervalCensoredData(dat[dat$Grouping=="A",])
m$q

plot(m$cfd[["Time"]],m$cfd[["Prob"]])


Dish<-get.groups(dat,dat$Dish)

frequency.table(Dish[[1]]$G)

dt<-Dish[[1]]

dat


toSurvival<-function(data,time="Time",G="G",PG="PG")
{
  tm<-c(0,(data[,"Time"]),Inf)
  nm<-c(0,data[,"G"],Dish[[1]][1,"PG"])

  n<-length(tm)
  intervals<-data.frame(Init=tm[1:(n-1)],End=tm[2:n],Freq=nm[2:n]-nm[1:(n-1)])
  intervals<-intervals[intervals[,3]>0,]
  n<-nrow(intervals)
  intervals<-intervals[rep.int(1:n,intervals[,3]),1:2]
  rownames(intervals)<-1:nrow(intervals)
  intervals
}

censoreddata<-rbind(toSurvival(Dish[[1]]),toSurvival(Dish[[2]]),toSurvival(Dish[[3]]),toSurvival(Dish[[4]]))

quantile(rowMeans(censoreddata),probs=0.5,type=4)


tabla.frecuencias((censoreddata[,2]))

tt<-table(censoreddata[,2],useNA="no")
t2<-cumsum(tt)
t2<-t2/max(t2)
val<-0.99
pos<-match(FALSE,t2<val, nomatch=length(t2))
val-t2[


cens<-Surv(censoreddata[,1],censoreddata[,2],type="interval2")


fit<-survfit(cens~1)
plot(fit)

median(fit)


summary(fit)

cbind(fit$time,1-fit$surv)

quantile(fit,probs=0.8,confint=FALSE)

tm<-c(0,(Dish[[1]][,"Time"]),Inf)
nm<-c(0,Dish[[1]][,"G"],Dish[[1]][1,"PG"])

n<-length(tm)
intervals<-cbind(tm[1:(n-1)],tm[2:n],nm[2:n]-nm[1:(n-1)])
intervals<-intervals[intervals[,3]>0,]
n<-nrow(intervals)
censoreddata<-intervals[rep.int(1:n,intervals[,3]),1:2]
n<-nrow(censoreddata)

cens<-Surv(censoreddata[,1],censoreddata[,2],type="interval2")

fit<-survfit(cens~1)

summary(fit)

quantile(fit,probs=0.5,confint=FALSE)
plot(fit)
axis(side = 1, at = 0:80)

names(fit)
tabla.frecuencias(dat[["Time"]])

library(survival) 
# Mirar ajuste de densidad vs distribución
survfit # Kapplan-Meyer???
