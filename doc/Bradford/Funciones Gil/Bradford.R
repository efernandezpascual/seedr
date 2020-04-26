setwd ("/home/gil/Documentos/Actual/Investigación Actual/INDUROT/Eduardo Fernandez Pascual")
source("/home/gil/Documentos/Actual/R/funciones.R")

X<-read.table("Bradford.csv",header=TRUE,sep=",",dec=".")
X<-X[,c(1:5)]
names(X)<-c("Temp","psi","T","Germination","F")
X$probit<-qnorm(X$F,0,1)

Bradford<-function(X,Y,Z)
{
  Sxy<-cov(X,Y)
  Sxz<-cov(X,Z)
  Sxx<-var(X)
  Syz<-cov(Y,Z)
  Syy<-var(Y)
  Szz<-var(Z)

  a <- -2*(Sxz^2)*Syz+2*Sxy*Sxz*Szz
  b<- 2*(Sxz^2)*Syy +4*Sxy*Syz*Sxz-2*(Sxy^2)*Szz-4*Syz*Sxy*Sxz
  c<-2*Syz*Sxy^2 - 2*Sxy*Sxz*Syy
    
#   a <-  -2*(cov(Z,X)^2)*cov(Y,Z)+2*cov(Y,X)*cov(Z,X)*var(Z)
#   b<- 2*(cov(Z,X)^2)*var(Y) +4*(cov(Y,X)*cov (Z,Y))*cov (Z,X)- 2*(cov(X,Y)^2)*var(Z)-4*cov (Z,Y)*cov (Y,X)*cov (Z,X)
#   c<-2*(cov(Z,Y))*cov(Y,X)^2 - 2*cov(Y,X)*cov(Z,X)*var(Y)
 
  max(c((-b+sqrt(b^2-4*a*c))/(2*a),(-b-sqrt(b^2-4*a*c))/(2*a)))
}

Bradford(X$probit,X$psi,1/X$T)







  
R2<-function(X,Y)
{
  (cov(X,Y)^2)/(var(X)*var(Y))
}

R2(X=X$probit,Y=X$psi-45.18/X$T)

# R2(45.18,Y=X$probit,H=X$psi,T=X$T)

grid<-seq(0,50,0.1)
r2<-vector("numeric",length(grid))
for (i in 1:length(grid))
{
  r2[i]<-R2(X=X$probit,Y=X$psi-grid[i]/X$T)
}
plot(grid,r2)

grid[which.max(r2)]
Bradford(X$probit,X$psi,1/X$T)

