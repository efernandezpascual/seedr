setwd ("/home/gil/Documentos/Actual/Investigación Actual/INDUROT/Eduardo Fernandez Pascual")
source("/home/gil/Documentos/Actual/R/funciones.R")

library(data.table)
library(seedr)

# Temp<-read.table("Temperatura.csv",sep=",",header=TRUE)
# head(Temp)
# head(peg)

head(centaury)

germinatio.data<-function(d, t, g, pg, x, reps = NULL, groups = NULL)
{
    dd <- data.table(d)
#     dd<-dd[-6]  ### Ejemplo para quitar un tiempo en un grupo
    dd<-dd[, .(germinated=sum(germinated),germinable=sum(germinable)), by = c(reps, groups,x,t)] # Se agrupa tb por t para unificar no grupos

#     dd[temperature==20] ### Observar el 6
    dd[,germinable:=max(germinable),by = c(x,reps, groups)] # Unificar los germinables
#     dd[temperature==20] ### Observar el 6

    setorderv(dd, c(reps, groups,x,t))
#     dd[temperature==20] ### Observar el 6


    dd[, cumulative := cumsum(get(g)), by = c(x, reps, groups)]
    dd[, germination := cumulative / get(pg)]
#     dd[temperature==20]
#     dd[temperature==0]

    return(list(data=dd,t=t,x=x,reps=reps,groups=groups))
}

# cdf columns Time and prob
quantileCDF<-function(cdf,prob=0.5,fx="times",fy="germination",extrapolate.prange=1)
{
  pos<-c()
  pos0<-c()
  for (i in 1:length(prob))
  {
    posA<-match(FALSE,cdf[[fy]]<prob[i], nomatch=NA)
    pos0A<-posA-1
    if (is.na(posA))	# Above maximum, use last and point of previous jump (to avoid flat areas)
    {
#       posA<-which.max(cdf[[fy]])
      posA<-length(cdf[[fy]])
      pos0A<-match(FALSE,cdf[[fy]]<cdf[[fy]][posA], nomatch=NA)-1
    }
    if (pos0A==0)
    {
      posA<-match(FALSE,cdf[[fy]]<=prob[i], nomatch=NA)
      pos0A<-posA-1
    }
    pos<-c(pos,posA)
    pos0<-c(pos0,pos0A)
  }
  p<-cdf[[fy]][pos]
  q<-cdf[[fy]][pos0]
  y<-cdf[[fx]][pos]
  x<-cdf[[fx]][pos0]

  m<-x+(prob-q)*(y-x)/(p-q)
  m[m>(extrapolate.prange*max(cdf[[fx]]))]<-NA
  return(m)
}

gd<-germinatio.data(d = centaury, t = "times", g = "germinated", pg = "germinable", x = "temperature", reps = "dish", groups = c("species","population"))

gd<-germinatio.data(d = centaury, t = "times", g = "germinated", pg = "germinable", x = "temperature", groups = c("species","population"))



quantileCDF(gd$data[temperature==20],prob=(1:10)/10)


as.data.frame(dt1)

qt<-split(dt.q,dt.q[,c(gd$reps, gd$groups, "q"),with=FALSE])
qt<-qt[[9]]
qt<-qt[qt$q==0.4,]

optimum.Ts<-function(x,r,pos,min.ptos=3)
{
  n<-length(r)

  a1<-NULL
  b1<-NULL
  if (pos>=min.ptos)
  {
    yy<-x[1:pos]
    xx<-r[1:pos]
    a1<-cov(xx,yy)/var(xx)
    b1<-mean(yy)-a1*mean(xx)
  }

  a2<-NULL
  b2<-NULL
  if ((n-pos+1)>=min.ptos)
  {
    yy<-x[pos:n]
    xx<-r[pos:n]
    a2<-cov(xx,yy)/var(xx)
    b2<-mean(yy)-a2*mean(xx)
  }

  na<-as.double(NA)
  if (is.null(b1))
  {
    if (is.null(b2)) return(list(Tb=na,Tc=na,To=na,Ro=na,Tmin=min(x,na.rm=FALSE),Tmax=max(x,na.rm=FALSE),Nb=pos,Nc=n-pos+1))
    return(list(Tb=na,Tc=b2,To=a2*r[pos]+b2,Ro=r[pos],Tmin=min(x,na.rm=FALSE),Tmax=max(x,na.rm=FALSE),Nb=pos,Nc=n-pos+1))
  } else
  {
    if (is.null(b2)) return(list(Tb=b1,Tc=na,To=a1*r[pos]+b1,Ro=r[pos],Tmin=min(x,na.rm=FALSE),Tmax=max(x,na.rm=FALSE),Nb=pos,Nc=n-pos+1))
    Ro=(b2-b1)/(a1-a2)
    return(list(Tb=b1,Tc=b2,To=a1*Ro+b1,Ro=Ro,Tmin=min(x,na.rm=FALSE),Tmax=max(x,na.rm=FALSE),Nb=pos,Nc=n-pos+1))
  }
}

max.R2<-function(x,r,min.ptos=3)
{
  n<-length(r)
  if (n<min.ptos) return(NA)
  opt<-as.integer(0)
  opt.val<-0
  for (pos in min.ptos:(n-min.ptos+1))
  {
    R2<-0

    yy<-x[1:pos]
    xx<-r[1:pos]
    cv<-cov(xx,yy)
    if (cv<=0) next
    R2<-(cv^2)/(var(xx)*var(yy))

    yy<-x[pos:n]
    xx<-r[pos:n]
    cv<-cov(xx,yy)
    if (cv>=0) next
    R2<-R2+(cv^2)/(var(xx)*var(yy))

    if (R2>opt.val)
    {
      opt<-pos
      opt.val<-R2
    }
  }

  cv<-cov(x,r)
  R2<-(cv^2)/(var(x)*var(r))
  if (cv>0) if (R2>opt.val) return(n)
  if (cv<0) if (R2>opt.val) return(1L)
  return(opt)
}

Gauss.fit<-function(x,germination)
{
  filter<-x>0
  x<-x[filter]
  germination<-germination[filter]

  probit = qnorm(germination, 0, 1)
  filter<-is.finite(probit)
  x<-x[filter]
  germination<-germination[filter]
  probit<-probit[filter]

  hlm <- lm.fit(matrix(c(x,rep(1,length(x))),ncol=2),probit)
  m<-hlm$coefficients[1]
  b<-hlm$coefficients[2]

  theta50 <- -b/m
  sigma <- 1/m

  data.frame(theta50, sigma)
}



estimate.Tb<-function(gd,grid.x=(1:10)/10,min.ptos=3,method=c("Max R2","Max value"),extrapolate.prange=1.2)
{
  dt.q<-gd$data[,.(q=grid.x, r=1/quantileCDF(cdf=.SD,prob=grid.x,fx="times",fy="germination",extrapolate.prange=extrapolate.prange)),by=c(gd$reps, gd$groups,gd$x)]
  dt.q<-dt.q[!is.na(dt.q$r)]

  if ("Max value" %in% method)  # Tb<=Tmin, Tc>=Tmax for all quantiles
  {
    Ts<-dt.q[,optimum.Ts(.SD[[gd$x]],r,which.max(r),min.ptos=min.ptos),by=c(gd$reps, gd$groups, "q")]
    Topt<-Ts[,.(Tb=min(c(mean(Tb,na.rm=TRUE),min(Tmin,na.rm=TRUE))),Tc=max(c(mean(Tc,na.rm=TRUE),max(Tmax),na.rm=TRUE)),To=mean(To,na.rm=TRUE)),by=c(gd$reps, gd$groups)]
  }

  if ("Max R2" %in% method)  # Tb<=Tmin, Tc>=Tmax for all quantiles
  {
    Ts<-dt.q[,optimum.Ts(.SD[[gd$x]],r,max.R2(.SD[[gd$x]],r,min.ptos=min.ptos),min.ptos=min.ptos),by=c(gd$reps, gd$groups, "q")]
    Topt<-Ts[,.(Tb=min(c(mean(Tb,na.rm=TRUE),min(Tmin,na.rm=TRUE))),Tc=max(c(mean(Tc,na.rm=TRUE),max(Tmax),na.rm=TRUE)),To=mean(To,na.rm=TRUE)),by=c(gd$reps, gd$groups)]
  }

  gd$data<-Topt[gd$data,on=c(gd$groups)]  # Left join
  gd$data[,suboptimal:=(get(gd$x)<=To)*1] # Suboptimal variable
  gd$data[,thetag:=get(gd$t)*((get(gd$x)-Tb)*suboptimal+(Tc-get(gd$x))*(1-suboptimal)),by=c(gd$groups)] # thetag computation


  ddd<-split(gd$data,gd$data[,c(gd$groups,"suboptimal"),with=FALSE])[[3]]
  x<-ddd[["thetag"]]
  germination<-ddd[["germination"]]
  germination


}

gd$data[200:260,]

plot(range(dt1$temperature),range(dt1$r),type="n")
ddd<-split(dt1,dt1$q)

pl<-rainbow(length(ddd))
for (i in 1:length(ddd)) points(ddd[[i]][,-2],col=pl[i])


dd[order(get(t)),]
dd[dd$temperature==20][order(get(t))]


?physiotime
brad<-physiotime(d = peg, t = "times", g = "germinated", pg = "germinable", x = "psi", reps = "dish",groups = c("species", "temperature"),method = "bradford")

brad<-physiotime(d = grasses, t = "times", g = "germinated", pg = "germinable", x = "psi", reps = "dish",groups = c("species", "temperature"),method = "bradford")

str(brad)
head(peg)

bb<-unclass(brad[[1]])

plot(bb[["data"]][["psi"]],bb[["data"]][["psibg"]])




brad[1]

class(brad[[1]])<-"bradford"

summary(brad)

plot(brad)

bradford(listx[[1]], t, g, pg, x, reps, groups)
listx[[1]]->d
d
grasses
ls()

y<-x
dd[order(get(y)),]
