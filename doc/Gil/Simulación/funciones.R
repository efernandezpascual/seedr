#  Useful functions for Estadística (grado) and Metodos_Estadisticos

# ------------------------------------------------------------------
# Note: Using R commander if the type of a variable is modified manually then R commander 
# menus are not updated automatically. In order to update the menus according to the current 
# kind of variables information execute: activateMenus() 
# ------------------------------------------------------------------

# Force enough digits to avoid problems
options(digits=10)

# ------------------------------------------------------------------
# FACTORS
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Transformation of a "numeric" factor to a vector of numbers
# ------------------------------------------------------------------
factor.to.numeric <- function(x)
{
  lv<-as.numeric(levels(x))
  if (sum(is.na(lv))>0) stop("Factor contains non-numeric levels")
  return(lv[x])
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Reorder the levels of a factor according to the frequencies
# Should be used only for nominal variables (encoded as a factor)
# Input: x is the factor to be reordered
#	 you can specify decreasing=TRUE to use a decreasing order
# Output: the reordered factor
# ------------------------------------------------------------------
factor.by.frequency <- function(x, ...)
{
  y<-table(factor(x))	# To get the frequency table of x
  z<-factor(x,levels=names(y)[order(y,...)])
  return(z)
}
sort.factor.by.frequency<-factor.by.frequency
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Introduce data by frequency
# ------------------------------------------------------------------
# x = data frame with a variable that contains the frequencies for each case
# freq = name of the variable that contains the frequencies - by default last column
data.by.frequency <- function(x,freq=names(x)[length(x)])
{
  if (!is.data.frame(x)) stop("need data.frame")
  as.data.frame(lapply(x,rep,x[,freq]))
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Returns the frequency tables (as explained in theory part)
# ------------------------------------------------------------------
# x = variable to be analyzed (if numeric then it will be converted to factor)
# ------------------------------------------------------------------
frequency.table<- function(x)
{
	if (!is.factor(x)) x<-factor(x,ordered=TRUE)
	tb<-table(x,useNA="always")
	nm<-names(tb)

	na.pos<-(1:length(tb))[is.na(names(tb))]
	nm[na.pos]<-"<NA>"
	names(tb)<-nm
	tbn<-tb[-na.pos]

	tbn.df<-data.frame(ni=as.vector(tbn))
	rownames(tbn.df)<-nm[-na.pos]
	tb.df<-data.frame(ni=as.vector(tb))
	rownames(tb.df)<-nm

	# Relative frequencies
	tbn.df$fi<-prop.table(tbn)
	tb.df$fi<-prop.table(tb)	

	# Cumulative - only for ordered factors
	if (is.ordered(x))
	{
		tbn.df$Ni<-cumsum(tbn.df$ni)
		tbn.df$Fi<-cumsum(tbn.df$fi)
		tb.df$Ni<-cumsum(tb.df$ni)
		tb.df$Fi<-cumsum(tb.df$fi)
	}

	rw<-tb.df[1,]
	rw[]<-NA
	rw[2]<-1
	rwn<-rw
	rw[1]<-sum(tb.df[,1])
	rwn[1]<-sum(tbn.df[,1])
	tb.df<-rbind(tb.df,"Total"=rw)
	tbn.df<-rbind(tbn.df,"Total"=rwn)

	# Structure each column as a table
	for (i in 1:length(tb.df))
	{
		tb.df[,i]<-as.table(array(tb.df[,i]))
		names(tb.df[,i])<-rownames(tb.df)
		tbn.df[,i]<-as.table(array(tbn.df[,i]))
		names(tbn.df[,i])<-rownames(tbn.df)
	}

	list(CompleteTable=tb.df,WithoutNATable=tbn.df)
}
tabla.frecuencias<-frequency.table
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Separate access to each column in a frequency table (without total)
# with.na=TRUE or FALSE
# type="ni" or "fi" or "Ni" or "Fi"
# The column is calculated even if it has no sense (for instance Fi for non-ordered factors, or Ni with NA values...)
# ------------------------------------------------------------------ 
freq<-function(x,with.na=TRUE,type="ni")
{	
  if (is.data.frame(x))
  {
    y<-list()
    nc<-c()
    for (i in names(x)) 
    {
      y[[i]]<-freq(x[[i]],with.na=with.na,type=type)
      nc<-c(nc,ncol(y[[i]]))
    }
    if (max(nc)-min(nc)>0) return(y)
    return(join.groups(y))
  }
  if (with.na) useNA<-"always" else useNA<-"no"
  if (!is.factor(x)) x<-as.factor(x)
  ret<-as.data.frame.matrix(t(as.matrix(table(x,useNA=useNA))))
  row.names(ret)<-type
  names(ret)[is.na(names(ret))]<-"<NA>"
  if (type=="fi") ret<-ret/sum(ret)
  if (type=="Ni") ret[1,]<-cumsum(as.numeric(ret[1,]))
  if (type=="Fi") ret[1,]<-cumsum(as.numeric(ret[1,]))/sum(ret)
  return(ret)
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# It tries to return the by.groups generated data in compact form
# ------------------------------------------------------------------
join.groups<-function(y,title=NA)
{
  nr<-nrow(y[[1]])
  nc<-ncol(y[[1]])
  if ((nr>1)&(nc>1)) return(y)

  ret<-c()
  for (i in 1:length(y))
  {
    if (is.null(names(y[[i]]))) names(y[[i]])<-title
    if (nr>1) y[[i]]<-t(y[[i]])
    ret<-rbind(ret,y[[i]])
  } 
  row.names(ret)<-names(y)
  return(ret)	
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Splits the data.frame 'data' by the levels of factor 'group'  
# ------------------------------------------------------------------
get.groups<-function(data,groups,...)
{
  split(data,groups)
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Applies function 'fun' to 'data' for each level of the variable 'groups'
# It tries to return the values in a compact (data.frame) form - 
# otherwise it returns a list indexed by the levels of the factor groups
# ------------------------------------------------------------------
by.groups<-function(data,groups,fun,...)
{
  nm<-deparse(substitute(fun))
  dtname<-deparse(substitute(data))
  x<-split(data,groups)
  y<-list()
	flag<-FALSE
  for (i in names(x)) 
  {
    y[[i]]<-fun(x[[i]],...)
		if (class(y[[i]])=="reglist")
		{
			flag<-TRUE	#Cannot be converted to data.frame and thus join.groups will not work
      names(y[[i]])<-dtname
		}
		else if (!is.data.frame(y[[i]]))
    {
      y[[i]]<-as.data.frame(y[[i]])
      names(y[[i]])<-dtname
    }
  }
	if (flag) return(y)
  return(join.groups(y,nm))
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Sample size and lost values summary
# ------------------------------------------------------------------
sample.size<-function(x)
{
  n<-length(x)
  n.na<-sum(is.na(x))
  n<-n-n.na
  return(data.frame(N=n,n.NA=n.na,row.names=""))
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Basic descriptive statistics - by default NA values are omitted
# ------------------------------------------------------------------
formals(mean.default)$na.rm = TRUE
formals(var)$na.rm=TRUE
quantile.default<-stats:::quantile.default
formals(quantile.default)$na.rm = TRUE

# Uncomment next line for using mid point by default in case of non-unique quantile
# formals(quantile.default)$type = 2

# ------------------------------------------------------------------
# mean function - redefined for factors (numeric), data.frames and matrix (per column)
# ------------------------------------------------------------------
mean<-function(x,trim=0,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...) 
{
  mn<-function(x,trim,na.rm,row.name,...) # Force interpretation of initial values
  {
    UseMethod("mean")
  }
  mn(x,trim,na.rm,row.name,...)
}
mean.factor<-function(x,...)
{
  mean(factor.to.numeric(x),...)
}
mean.data.frame<-function(x,trim,na.rm,row.name,...)
{
  y<-x[1,,drop=FALSE]
  row.names(y)<-row.name
  for (i in names(x)) y[[i]]<-mean(x[[i]],trim,na.rm,row.name,...)
  return(y)
}
mean.matrix<-function(x,trim,na.rm,row.name,...) mean.data.frame(as.data.frame(x),trim,na.rm,row.name,...)
media <- mean
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Clasical variance function: By default var computes the sample variance
# ------------------------------------------------------------------
variance <- function (x, row.name=deparse(match.call()[[1]]), ...) 
{
  if (is.factor(x)) x<-factor.to.numeric(x)
  mean(x^2,row.name=row.name, ...) - mean(x, row.name=row.name, ...)^2
}
varianza<-variance
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Standard Deviation
# ------------------------------------------------------------------
standard.deviation <- function (x, row.name=deparse(match.call()[[1]]), ...)
{
  sqrt(variance(x, row.name=row.name, ...))
}
desviacion.tipica<-standard.deviation
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Skewness/Asymmetry
# ------------------------------------------------------------------
skewness <- function (x, row.name=deparse(match.call()[[1]]), ...)
{
  if (is.factor(x)) x<-factor.to.numeric(x)
  y<-mean(x,row.name=row.name)
  if (is.data.frame(y)) y<-y[rep(1,nrow(x)),]
  
  mean((x-y)^3, row.name=row.name,...) / standard.deviation(x, row.name=row.name, ...)^3
}
asymmetry <- skewness
asimetria <- skewness
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Kurtosis
# ------------------------------------------------------------------
kurtosis <- function (x, row.name=deparse(match.call()[[1]]), ...)
{
  if (is.factor(x)) x<-factor.to.numeric(x)
  y<-mean(x,row.name=row.name)
  if (is.data.frame(y)) y<-y[rep(1,nrow(x)),]

  mean((x-y)^4, row.name=row.name, ...) / variance(x, row.name=row.name, ...)^2 - 3
}
curtosis <- kurtosis
apuntamiento <- kurtosis
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Coefficient of Variation (CV)
# ------------------------------------------------------------------
CV <- function (x, row.name=deparse(match.call()[[1]]), ...)
{
  standard.deviation(x, row.name=row.name, ...) / abs(mean(x, row.name=row.name, ...))
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Full set of statistics related with the mean
# Factors are coerced to numeric vectors
# ------------------------------------------------------------------
mean.stats <- function(data,...)
{
  if (is.data.frame(data)) 
  {
    if (ncol(data)==1) return(mean.stats(data[[1]]))
    else 
    {
      ret<-c()
      for (i in names(data)) 
      {
	y<-mean.stats(data[[i]])
	row.names(y)<-i
	ret<-rbind(ret,y)
      }
      return(as.data.frame(ret))
    }
  }
  if (is.factor(data)) data<-factor.to.numeric(data)
  return(data.frame(mean=mean(data),standard.deviation=standard.deviation(data),variance=variance(data),CV=CV(data),skewness=skewness(data),kurtosis=kurtosis(data),sample.size(data),row.names=""))
}
estadisticos.media<-function(data)
{
  ret<-mean.stats(data)
  names(ret)[1:6]<-c("media","desviacion.tipica","varianza","CV","asimetria","cursosis")
  return(ret)
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Ordered statistics: 
# For factors min and max values are returned for each probability 
# ------------------------------------------------------------------

#formals(median.default)$na.rm=TRUE
#formals(quantile.default)$na.rm=TRUE

minimo<-function(data,row.name="Mínimo")
{
	if (is.factor(data)) return(min(factor.to.numeric(data),na.rm=TRUE))
	if (is.data.frame(data))
	{
    y<-data[1,]
		row.names(y)<-row.name
    for (i in names(data)) 
    {
			y[[i]][1]<-min(data[[i]],na.rm=TRUE)
		}		
		return(y)
	}
	min(data,na.rm=TRUE)
}

maximo<-function(data,row.name="Máximo")
{
	if (is.factor(data)) return(max(factor.to.numeric(data),na.rm=TRUE))
	if (is.data.frame(data))
	{
    y<-data[1,]
		row.names(y)<-row.name
    for (i in names(data)) 
    {
			y[[i]][1]<-max(data[[i]],na.rm=TRUE)
		}		
		return(y)
	}
	max(data,na.rm=TRUE)
}

recorrido<-function(data,row.name="Recorrido")
{
	maximo(data,row.name)-minimo(data,row.name)
}

amplitud<-recorrido

minimum<-minimo
maximum<-maximo
variation.range<-recorrido

formals(minimum)$row.name<-"Minimum"
formals(maximum)$row.name<-"Maximum"
formals(variation.range)$row.name<-"Range"



# ------------------------------------------------------------------
# Quantile function: use type=2 for mid point with numerical variables
# ------------------------------------------------------------------
quantile<-function(x, probs = seq(0, 1, 0.25), na.rm=TRUE, ...)
{
  if (is.factor(x))
  {
    if (na.rm) x<-x[!is.na(x)]
    n<-length(x)
    if (n==0L) return(x[FALSE][NA])
    if (!is.ordered(x)) warning("factor is coerced to ordered factor")
    num<-as.numeric(x)
    x<-x[order(num)]
    index <- 1 + (n - 1) * probs
    lo <- floor(index)
    hi <- ceiling(index)
    r<-data.frame(min_value=x[lo],max_value=x[hi])
    row.names(r)<-paste(100*probs,"%",sep="")
    return(r)
  }
  if (is.data.frame(x)) 
  {
    y<-c()
    for (i in names(x)) 
    {
      if (is.factor(x[[i]])) 
      {
	stop("multiple factors are not allowed yet")
	return
      }
      y<-rbind(y,quantile(x[[i]],probs=probs,na.rm=na.rm,...))
    }
    y<-as.data.frame(t(y))
    names(y)<-names(x)
    return(y)
  }
  UseMethod("quantile")
}
cuantiles<-quantile
quantiles<-quantile
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Quartiles 
# ------------------------------------------------------------------
quartiles <- function(x, ...)
{
  quantile(x,probs=c(0.25,0.5,0.75))
}
cuartiles <- quartiles
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Deciles
# ------------------------------------------------------------------
deciles <- function(x, ...)
{
  quantile(x,probs=(1:9)/10)
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Median
# ------------------------------------------------------------------
median<-function(x,na.rm=TRUE,row.name=deparse(match.call()[[1]]), ...)
{
  if (is.factor(x)|is.data.frame(x))
  {
    ret<-quantile(x,probs=0.5,...)
    row.names(ret)<-row.name
    return(ret)
  }
  mn<-function(x,na.rm,...) # Force interpretation of initial values
  {
    UseMethod("median")
  }
  mn(x,na.rm,...)
}
mediana<-median
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Inter Quartile Range: IQR
# ------------------------------------------------------------------
IQR <- function (x,na.rm=TRUE,row.name=deparse(match.call()[[1]]), ...)
{
  if (is.factor(x))
  {
    y<-quantile(x,probs=c(0.25, 0.75),na.rm)
    zlow<-as.numeric(y$low)
    if (prod(as.character(factor(zlow))==as.character(y$low)))	#Factor is numerical - compute IQR
    {
      zhigh<-as.numeric(y$high)
      ret<-data.frame(min_value=zlow[2]-zhigh[1], max_value=zhigh[2]-zlow[1])
      row.names(ret)<-row.name
      return(ret)
    }
    else return(x[FALSE][NA])
  }
  ret<-quantile(x, probs=c(0.25, 0.75), na.rm, ...)
  if (is.data.frame(x))
  {
    ret<-diff(as.matrix(ret))
    row.names(ret)<-row.name
    return(as.data.frame(ret))
  }
  names(ret)<-NULL
  diff(ret)
}
iqr<-IQR
RIC <- IQR
recorrido.intercuartilico<-IQR
interquartile.range<-IQR
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# MEDA
# ------------------------------------------------------------------
MEDA <- function(x,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...)
{
  if (is.factor(x)) stop("need numeric data")
  m<-median(x,na.rm=TRUE,row.name=row.name,...)
  if (is.data.frame(x)) 
  {
    m<-m[rep(1,nrow(x)),]
  }
  return(median(abs(x-m),na.rm=TRUE,row.name=row.name,...))
}
meda<-MEDA
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Variation Index based on median and meda
# ------------------------------------------------------------------
median.VI <- function(x,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...)
{
  if (is.factor(x)) stop("need numeric data")
  return(MEDA(x,na.rm=na.rm, row.name=row.name,...)/abs(median(x,na.rm=na.rm, row.name=row.name,...)))
}
median.variationindex <- median.VI
indice.variacion.mediana <- median.VI
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Full set of statistics related with the median
# Factors are coerced to numeric vectors
# ------------------------------------------------------------------
median.stats <- function(x)
{
  if (is.factor(x)) x<-factor.to.numeric(x)
  if (is.data.frame(x)) 
  {
    if (ncol(x)==1) return(median.stats(x[[1]]))
    else 
    {
      ret<-c()
      for (i in names(x)) 
      {
	y<-median.stats(x[[i]])
	row.names(y)<-i
	ret<-rbind(ret,y)
      }
      return(as.data.frame(ret))
    }
  }
  
  m<-median(x)
  md<-median(abs(x-m))
  return(data.frame(median=m,MEDA=md,median.VI=md/abs(m),IQR=IQR(x),sample.size(x),row.names=""))
}
estadisticos.mediana<-function(x)
{
  ret<-median.stats(x)
  names(ret)[1:4]<-c("mediana","MEDA","indice.variacion.mediana","RIC")
  return(ret)
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Returns the trimmed data set - put NA for each discarded data (value.low=NA,value.upp=NA)
# trim=Overall minimum proportion of data to be trimmed (trim/2 at each side)
# ------------------------------------------------------------------
trim.data<-function(x,trim=0.05,value.low=NA,value.upp=NA)
{
  if (is.data.frame(x))
  {
    for (i in names(x)) x[[i]]<-trim.data(x[[i]],trim=trim,value.low=value.low,value.upp=value.upp)
    return(x)
  }

  per<-order(x)
  n<-length(x)-sum(is.na(x))
  trim<-trim/2
  m<-ceiling(n*trim)
  M<-floor(n*(1-trim))+1
  x[per[1:m]]<-value.low
  x[per[M:n]]<-value.upp
  return(x)
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Trimmed mean function
trim.mean <- function(x,trim=0.05,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...)
{
  mean(trim.data(x,trim=trim),na.rm=na.rm, row.name=row.name,...)
}
media.recortada<-trim.mean
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Trimmed Clasical variance function
trim.variance <- function(x,trim=0.05,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...)
{
  variance(trim.data(x,trim=trim),na.rm=na.rm, row.name=row.name,...)
}
varianza.recortada<-trim.variance
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Trimmed Standard Deviation
trim.standard.deviation <- function (x,trim=0.05,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...)
{
  standard.deviation(trim.data(x,trim=trim),na.rm=na.rm, row.name=row.name,...)
}
desviacion.tipica.recortada<-trim.standard.deviation
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Trimmed Coefficient of Variation (CV)
trim.CV <- function (x,trim=0.05,na.rm=TRUE,row.name=deparse(match.call()[[1]]),...)
{
  CV(trim.data(x,trim=trim),na.rm=na.rm, row.name=row.name,...)
}
CV.recortado<-trim.CV
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Computes the mean, standard deviation, variance and CV for the trimmed data
# Factors are coerced to numeric vectors
# ------------------------------------------------------------------
trim.stats<-function(x,trim=0.05)
{
  mean.stats(trim.data(x,trim=trim))
}
estadisticos.recortados<-function(x,trim=0.05)
{
  estadisticos.media(trim.data(x,trim=trim))
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Moda - Only for factors
# ------------------------------------------------------------------
moda <- function(x,row.name=deparse(match.call()[[1]]))
{
  if (is.data.frame(x))
  {
    y<-x[1,,drop=FALSE]
    row.names(y)<-row.name
    for (i in names(x)) y[[i]]<-moda(x[[i]],row.name)
    return(y)
  }
  if (!is.factor(x)) stop("need factor data")
  y<-table(x)
  names(y[y==maximo(y)])
}

mode.statistic<-moda





# ------------------------------------------------------------------
# Extends barplot functionality to data.frames
# ------------------------------------------------------------------
formals(barplot.default)$beside = TRUE
barplot.data.frame<-function(height,legend.text=(nrow(height)>1),...)
{
  barplot.default(as.matrix(height),legend.text=legend.text,...)
}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Histogram with density plot
# dens=    
#       "kernel": Kernel density estimation with BW by cross-validation
#       "normal": normal density fitted by sample mean and sample standard deviation
#       function: specified function to be plotted. Ex: dens=function(x){dnorm(x,mean=3,sd=1)}
#       NULL=No density plot
# dens.col: color for the density function
# dens.lwd: width of the line for plotting the density function
# xlab: label for the X axis
# main: main title of the graphic
histogram <- function(x, dens="kernel", dens.col="red", dens.lwd=2, 
main=paste("Histogram of",xlab,eval(if (!is.function(dens)) {if (dens=="kernel") "and density estimation (Bw=CV)" else if (dens=="normal") "and Gaussian estimation"})), xlab=xname,...)
{
  xname<-paste(deparse(substitute(x), 500), collapse = "\n")
  if (is.data.frame(x)|is.factor(data)) stop("need numeric data")
  x<-x[!is.na(x)]
  if (!is.null(dens)) hist(x,freq=FALSE,xlab=xlab,main=main,...) else hist(x,xlab=xlab,main=main,...)
  
  if (is.function(dens))
  {
    s<-sd(x)
    grid<-seq(min(x)-s,max(x)+s,length.out=400)   
    lines(grid,dens(grid),col=dens.col,lwd=dens.lwd)
  }
  else
  {
    if (dens=="kernel")
    {
      d<-density(x,bw=bw.bcv(x, nb = 1000))
      lines(d$x,d$y,col=dens.col,lwd=dens.lwd)
    }
    if (dens=="normal")
    {
      m<-mean(x)
      s<-sd(x)
      grid<-seq(min(x)-s,max(x)+s,length.out=400)
      d<-dnorm(grid,m,s)
      lines(grid,d,col=dens.col,lwd=dens.lwd)
    }  
  }
}
histograma<-histogram
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Graphical and numerical diagnostic for linear regression
# ------------------------------------------------------------------
regresion.diagnostico<-function(modelo)
{
  rs<-scale(modelo[["residuals"]])
  fs<-scale(modelo[["fitted.values"]])
  
  # Identify those Zres outside [-3,3]
  pos<-abs(rs)>3
  x<-fs[pos]
  y<-rs[pos]
  lab<-row.names(rs)[pos]
  ps<-2*(x<0)+2

  oldpar <- par(oma=c(0,0,0,0), mfrow=c(2,2))
  plot(modelo$model[,c(2,1)],main="Nube de puntos y modelo estimado")
  abline(modelo[["coefficients"]],lwd=2, col="red")
  points(modelo$model[pos,c(2,1)],pch=19,col="red")  
  
  plot(fs,rs,xlab="Zpred",ylab="Zres",main="Residuos vs. valores predichos")
  abline(h=c(-3,0,3), lty = 3, lwd=2, col = "gray50")
  if (length(x))
  {
    points(x,y,pch=19,col="red")  
    text(x,y,lab,pos=ps)  
  }
  
  ylim <- range(rs, na.rm = TRUE)
  ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
  qq<-qqnorm(rs,ylab = "Zres", ylim = ylim, xlab="Cuantiles de la distribución normal")
  qqline(rs, lty = 3, lwd=2, col = "gray50")
  points(qq$x[pos],qq$y[pos],pch=19,col="red") 
  
  histograma(rs, main="Residuos (st.) vs. curva normal",xlab="Zres", ylab="Densidad")
  grid<-seq(ylim[1],ylim[2],length.out=100)
  d<-dnorm(grid)
  lines(grid,d,col="blue",lwd=3,lty=4)
  legend("topleft",legend=c("Estimación kernel","Curva normal"),fill=c("red","blue"),cex=0.8)
  par(oldpar)
  
  library(lmtest)
  dw<-dwtest(modelo$model[,1] ~ modelo$model[,2], alternative="two.sided")
  dw$data.name<-paste(eval(names(modelo$model)[1]),"~",eval(names(modelo$model)[2]))
  return(dw)
}
regression.diagnostics<-regresion.diagnostico
# ------------------------------------------------------------------




# ------------------------------------------------------------------
# Usual Regression Models
# Linear model is always included for comparative purpouses
# Input: x and y variables
#				 data: data.frame containing x and y variables
#				 models: models to be analyzed - currently "Logaritmico","Inverso","Cuadratico","Cubico","Potencial" and/or "Exponencial"
# ------------------------------------------------------------------
Regression.Models<-function(data,x,y,models=c("Logaritmico","Inverso","Cuadratico","Cubico","Potencial","Exponencial"),transform=TRUE)
{
	envir=parent.frame()
  models<-c("Lineal",models)
  
  xname<-deparse(substitute(x))
  yname<-deparse(substitute(y))
  dataname<-deparse(substitute(data))

  ret<-list()
  if ('Lineal' %in% models) {ret[["Lineal"]]=eval(parse(text=paste("lm(formula=",yname,"~",xname,", data=",dataname,")",sep="")),envir=envir); ret[["Lineal"]]$Desc<-"Lineal: y(x)=b0+b1*x"}
  if ('Logaritmico' %in% models) {ret[["Logaritmico"]]=eval(parse(text=paste("lm(formula=",yname,"~I(log(",xname,")), data=",dataname,")",sep="")),envir=envir);ret[["Logaritmico"]]$Desc<-"Logarítmico: y(x)=b0+b1*ln(x)"} 
  if ('Inverso' %in% models) {ret[["Inverso"]]=eval(parse(text=paste("lm(formula=",yname,"~I(1/",xname,"), data=",dataname,")",sep="")),envir=envir); ret[["Inverso"]]$Desc<-"Inverso: y(x)=b0+b1/x"}
  if ('Inverso0' %in% models) {ret[["Inverso0"]]=eval(parse(text=paste("lm(formula=",yname,"~I(1/",xname,")-1, data=",dataname,")",sep="")),envir=envir); ret[["Inverso0"]]$Desc<-"Inverso0: y(x)=b0/x"}
  if ('Cuadratico' %in% models) {ret[["Cuadratico"]]=eval(parse(text=paste("lm(formula=",yname,"~I(",xname,") + I(",xname,"^2), data=",dataname,")",sep="")),envir=envir); ret[["Cuadratico"]]$Desc<-"Cuadrático: y(x)=b0+b1*x+b2*x^2"}
  if ('Cubico' %in% models) {ret[["Cubico"]]=eval(parse(text=paste("lm(formula=",yname,"~",xname,"+ I(",xname,"^2) + I(",xname,"^3), data=",dataname,")",sep="")),envir=envir); ret[["Cubico"]]$Desc<-"Cúbico: y(x)=b0+b1*x+b2*x^2+b3*x^3"}
  if ('Potencial' %in% models) {ret[["Potencial"]]=eval(parse(text=paste("lm(formula=log(",yname,")~I(log(",xname,")), data=",dataname,")",sep="")),envir=envir); ret[["Potencial"]]$Desc<-"Potencial: y(x)=b0*x^b1"}
  if ('Exponencial' %in% models) {ret[["Exponencial"]]=eval(parse(text=paste("lm(formula=log(",yname,")~",xname,", data=",dataname,")",sep="")),envir=envir); ret[["Exponencial"]]$Desc<-"Exponencial: y(x)=b0*exp(b1*x)"} 

  if (!transform) #Non-log transformation - nls estimation instead
  {
    if ('Potencial' %in% models) {start<-list(b0=exp(ret[['Potencial']]$coefficients[1]), b1=ret[['Potencial']]$coefficients[2]);   
    aux=eval(parse(text=paste("nls(formula=",yname,"~b0*",xname,"^b1, data=data,start=start)",sep=""))); 
    ret[["Potencial"]]$coefficients[1]<-log(coef(aux)[1]); ret[["Potencial"]]$coefficients[2]<-coef(aux)[2];} 
    if ('Exponencial' %in% models) {start<-list(b0=exp(ret[['Exponencial']]$coefficients[1]), b1=ret[['Exponencial']]$coefficients[2]);   
    aux=eval(parse(text=paste("nls(formula=",yname,"~b0*exp(b1*",xname,"), data=data,start=start)",sep=""))); 
    ret[["Exponencial"]]$coefficients[1]<-log(coef(aux)[1]); ret[["Exponencial"]]$coefficients[2]<-coef(aux)[2];} 
  }
  
  
  class(ret)<-"reglist"
  
  return(ret)
}

Modelos.Regresion<-Regression.Models

# ------------------------------------------------------------------
# Specialized methods for reglist class
# ------------------------------------------------------------------
coefficients.reglist<-function(x)
{
  n<-length(x)
  nm<-names(x)

  pos<-1  
  res<-data.frame(b0=rep(NA,n),b1=rep(NA,n),b2=rep(NA,n),b3=rep(NA,n),"R.squared"=rep(NA,n))
  for (i in x)
  {
    s<-summary(i)
    row.names(res)[pos]<-i$Desc
    res[pos,1:4]<-i$coefficients[1:4]
    res[pos,5]<-s$r.squared  
    if (grepl("Potencial", i$Desc)||grepl("Exponencial", i$Desc))
    {
      res[pos,1]<-exp(res[pos,1])
      #Recompute also the R.squared coefficient
      ydata<-subset(x[[1]]$model,select=1)
      xdata<-subset(x[[1]]$model,select=2)
      yest<-exp(predict(i,xdata))
      res[pos,5]<-1-variance((ydata-yest))/variance(ydata)
    }
    pos<-pos+1
  }
  if (!("Cubico" %in% nm)) res$b3<-NULL
  if (!("Cuadratico" %in% nm)) res$b2<-NULL
  return(res)
}

coef.reglist<-coefficients.reglist

print.reglist<-function(x)
{
  res<-coefficients(x)
  aa<-format.data.frame(res, digits = NULL, na.encode = FALSE, justify="right")
  aa[is.na(res)]<-gsub("NA","  ",aa[is.na(res)])
  print(aa)
}

summary.reglist<-function(x)
{
  lapply(x,summary)
}

my.predict<-function(x,newdata)
{
  if (grepl("Potencial", x$Desc)||grepl("Exponencial", x$Desc)) return(exp(predict(x,newdata)))
  return(predict(x,newdata))
}

predict.reglist<-function(object, newdata)
{
  if (missing(newdata) || is.null(newdata)) newdata<-subset(object[[1]]$model,select=2)
  x<-sapply(object,my.predict,newdata)
  if (is.vector(x))
  {
    x<-matrix(x,ncol=length(x),dimnames=list(NA,names(object)))
  }
  row.names(x)<-newdata[[1]]
  x
}

plot.reglist<-function(object,data=object$Lineal$model[2:1], f.lty=c(1:6,1:6)[1:length(object)],f.col=rainbow(length(object)),f.lwd=rep(2,length(object)),...)
{
  plot(data,...)
  
  x.val<-seq(min(data[,1]),max(data[,1]),length.out=100)
  dt<-data[rep(1,100),1,drop=FALSE]
  dt[,1]<-x.val
  y.val<-predict(object,dt)

  for (i in 1:ncol(y.val))
  {
    lines(x.val,y.val[,i],col=f.col[i],lty=f.lty[i],lwd=f.lwd[i])
  }
  
  if (object$Lineal$coefficients[2]>=0) legend("topleft",col=f.col,lty=f.lty,lwd=f.lwd,legend=names(object))
  else legend("topright",col=f.col,lty=f.lty,lwd=f.lwd,legend=names(object))
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Bar plot with associated probability distribution (discrete variables)
# data: sample data
# values: different values of the discrete random variable. Example: values=c(1,2,3)
# prob: associated probabilities for each value. Example: prob=c(0.2,0.5,0.3)
barplot.distribution<-function(data,values,prob,...)
{
  if (is.data.frame(data)) stop("need numeric data")
  if (is.factor(data)) f<-data else f<-factor(data,values)
  f<-f[!is.na(f)]
  t<-table(f)
  h<-matrix(c(prob,t/length(data)),nrow=2,ncol=length(values),byrow=TRUE)
  rownames(h)<-c("Probability","Frequency")
  colnames(h)<-values
  barplot(h, beside=TRUE, legend.text=TRUE,...)
}
# ------------------------------------------------------------------


#plot.continuousdensfreq<-function(data,grid,density,numint="Sturges")
#{
#  hist(data,freq=FALSE,xlim=c(min(grid),max(grid)),breaks=numint)
#  lines(grid,density,col="red")
#}

if (!exists("prop.test.old")) prop.test.old<-prop.test

prop.test<-function (x, n, p = NULL, alternative = c("two.sided", "less", "greater"), conf.level = 0.95, correct = TRUE)
{
    DNAME <- deparse(substitute(x))
    if (is.table(x) && length(dim(x)) == 1L) {
        if (dim(x) != 2L) 
            stop("table 'x' should have 2 entries")
        l <- 1
        n <- sum(x)
        x <- x[1L]
    }
    else if (is.matrix(x)) {
        if (ncol(x) != 2L) 
            stop("'x' must have 2 columns")
        l <- nrow(x)
        n <- rowSums(x)
        x <- x[, 2L]
    }
    else {
        DNAME <- paste(DNAME, "out of", deparse(substitute(n)))
        if ((l <- length(x)) != length(n)) 
            stop("'x' and 'n' must have the same length")
    }
    OK <- complete.cases(x, n)
    x <- x[OK]
    n <- n[OK]
    if ((k <- length(x)) < 1L) 
        stop("not enough data")
    if (any(n <= 0)) 
        stop("elements of 'n' must be positive")
    if (any(x < 0)) 
        stop("elements of 'x' must be nonnegative")
    if (any(x > n)) 
        stop("elements of 'x' must not be greater than those of 'n'")
    if (is.null(p) && (k == 1)) 
        p <- 0.5
    if (!is.null(p)) {
        DNAME <- paste(DNAME, ", null ", if (k == 1) 
            "probability "
        else "probabilities ", deparse(substitute(p)), sep = "")
        if (length(p) != l) 
            stop("'p' must have the same length as 'x' and 'n'")
        p <- p[OK]
        if (any((p <= 0) | (p >= 1))) 
            stop("elements of 'p' must be in (0,1)")
    }
    alternative <- match.arg(alternative)
    if (k > 2 || (k == 2) && !is.null(p)) 
        alternative <- "two.sided"
    if ((length(conf.level) != 1L) || is.na(conf.level) || (conf.level <= 
        0) || (conf.level >= 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    correct <- as.logical(correct)
    ESTIMATE <- x/n
    names(ESTIMATE) <- if (k == 1) 
        "p"
    else paste("prop", 1L:l)[OK]
    NVAL <- p
    CINT <- NULL
    YATES <- if (correct && (k <= 2)) 
        0.5
    else 0
    if (k == 1) {
        z <- qnorm(if (alternative == "two.sided") 
            (1 + conf.level)/2
        else conf.level)
        YATES <- min(YATES, abs(x - n * p))
        z22n <- z^2/(2 * n)
        p.c <- ESTIMATE + YATES/n
        p.u <- if (p.c >= 1) 
            1
        else (p.c + z22n + z * sqrt(p.c * (1 - p.c)/n + z22n/(2 * 
            n)))/(1 + 2 * z22n)
        p.c <- ESTIMATE - YATES/n
        p.l <- if (p.c <= 0) 
            0
        else (p.c + z22n - z * sqrt(p.c * (1 - p.c)/n + z22n/(2 * 
            n)))/(1 + 2 * z22n)
        CINT <- switch(alternative, two.sided = c(max(p.l, 0), 
            min(p.u, 1)), greater = c(max(p.l, 0), 1), less = c(0, 
            min(p.u, 1)))
    }
    else if ((k == 2) & is.null(p)) {
        DELTA <- ESTIMATE[1L] - ESTIMATE[2L]
        YATES <- min(YATES, abs(DELTA)/sum(1/n))
        WIDTH <- (switch(alternative, two.sided = qnorm((1 + 
            conf.level)/2), qnorm(conf.level)) * sqrt(sum(ESTIMATE * 
            (1 - ESTIMATE)/n)) + YATES * sum(1/n))
        CINT <- switch(alternative, two.sided = c(max(DELTA - 
            WIDTH, -1), min(DELTA + WIDTH, 1)), greater = c(max(DELTA - 
            WIDTH, -1), 1), less = c(-1, min(DELTA + WIDTH, 1)))
    }
    if (!is.null(CINT)) 
        attr(CINT, "conf.level") <- conf.level
    METHOD <- paste(if (k == 1) 
        "1-sample proportions test"
    else paste0(k, "-sample test for ", if (is.null(p)) 
        "equality of"
    else "given", " proportions"), if (YATES) 
        "with"
    else "without", "continuity correction")
    if (is.null(p)) {
        p <- sum(x)/sum(n)
        PARAMETER <- k - 1
    }
    else {
        PARAMETER <- k
        names(NVAL) <- names(ESTIMATE)
    }
    names(PARAMETER) <- "df"
    x <- cbind(x, n - x)
    E <- cbind(n * p, n * (1 - p))
    if (any(E < 5)) 
        warning("Chi-squared approximation may be incorrect")
    STATISTIC <- sum((abs(x - E) - YATES)^2/E)
    names(STATISTIC) <- "X-squared"
    if (alternative == "two.sided") 
        PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    else {
        if (k == 1) 
            z <- sign(ESTIMATE - p) * sqrt(STATISTIC)
        else z <- sign(DELTA) * sqrt(STATISTIC)
        PVAL <- pnorm(z, lower.tail = (alternative == "less"))
    }
    RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = as.numeric(PVAL), estimate = ESTIMATE, null.value = NVAL, 
        conf.int = CINT, alternative = alternative, method = METHOD, 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}




# ------------------------------------------------------------------
# Hypothesis tests and confidence intervals - bootstrap 
# ------------------------------------------------------------------

# ------------------------------------------------------------------
#Confidence interval for the proportion based on SCORE method
# x = sample data - numeric or numeric factor
# conf.level = confidence level
# ------------------------------------------------------------------
score.CI<-function(x,conf.level=0.95)
{
  if (is.data.frame(x)) stop("need numeric data")
  if (is.factor(x)) x<-factor.to.numeric(x)	# Two level factor!!!
  x<-x[!is.na(x)]
  p<-mean(x)
  n<-length(x)
  za<-qnorm(conf.level+(1-conf.level)/2)
  za2<-za*za
  center<-p+za2/(2*n)
  spread<-za*sqrt((p*(1-p)+za2/(4*n))/n)
  CI<-c(center-spread, center+spread, center, spread)/(1+za2/n)
  names(CI)<-c("Lim.Inf","Lim.Sup","Mid","Spread")
  return (CI)
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Trick for fixing the seed for lesson evaluation 
# ------------------------------------------------------------------
b.seed<-NA
boot.set.seed<-function(seed=123456789)
{
	b.seed<<-seed
}

# ------------------------------------------------------------------
# Studentized bootstrap test and CI about the mean
# x = sample data - numeric or numeric factor
# alternative = type of alternative hypothesis: "two.sided", "less", "greater"
# mu = expected value under the equality of the null hypothesis
# conf.level = confidence level
# B = number of bootstrap iterations
# boot.sample.size = sample size for bootstrap resamplings 
# ------------------------------------------------------------------
bt.test<-function(x, alternative='two.sided', mu=0.0, conf.level=.95, B=1000, proportion=FALSE, boot.sample.size=n, label="Si", get.sample=FALSE)
{
	if (!is.na(b.seed)) 
	{
		set.seed(b.seed)
		warning("Remember that seed is fixed for bootstrap procedures")
	}
  if (is.data.frame(x)) stop("need numeric data")
  if (is.factor(x))
	{	
		if (proportion) x<-(x==label)*1
		else x<-factor.to.numeric(x)	# Be carefull when using factor. Only numerical ones or two level factors should be used.
	}
  x<-x[!is.na(x)]
  muestra<-x
  n<-length(x)
  den<-sqrt(var(x))
  if (den<1E-10) den<-1
  estad<-(mean(x)-mu)/den 	#Computation of the value of the statistic

  m.muestra<-mean(muestra)  
  
  s<-matrix(sample(muestra, boot.sample.size*B, replace=T),nrow=boot.sample.size,ncol=B)
  dv<-apply(s,2,sd)
  dv[dv<1E-10]<-den
  bootestad<-(apply(s,2,mean)-m.muestra)/dv


#   bootestad<-replicate(B,{
#   s<-sample(muestra, boot.sample.size, replace=T)
#   std<-sqrt(var(s))
#   if (std<1E-10) (mean(s)-m.muestra)/den else (mean(s)-m.muestra)/std
#   })
  
#   xboot<- replicate(B,sample(muestra, boot.sample.size, replace=T))
#   m.muestra<-mean(muestra)
#   ZN<-function(x){(mean(x)-m.muestra)/sqrt(var(x))}	  #Studentized bootstrap both for CI and test
#   bootestad<-apply(xboot,2,ZN)

  alfa<- 1-conf.level
  l.ic<- quantile(bootestad, alfa/2)
  u.ic<- quantile(bootestad, 1-alfa/2)

  IC<-c(m.muestra-u.ic*den, m.muestra-l.ic*den)
  IC<-c(IC,(IC[2]+IC[1])/2,(IC[2]-IC[1])/2)
  names(IC)<-c("Lim.Inf","Lim.Sup","Mid","Spread")

  if (proportion==TRUE) CI.ret<-list(confidence=conf.level, boot.confidence.interval=IC, score.confidence.interval=score.CI(x,conf.level)) 
  else CI.ret<-list(confidence=conf.level, boot.confidence.interval=IC)

  if (get.sample==TRUE) CI.ret$boot.sample=xboot
  
  switch(alternative,
              two.sided=c(estimate=mean(x), statistic=estad, p.value=mean(abs(bootestad)  >= abs(estad)),CI.ret),
              less=c(estimate=mean(x), statistic=estad, p.value=mean(bootestad<=estad),CI.ret),
              greater=c(estimate=mean(x), statistic=estad, p.value=mean(bootestad>=estad),CI.ret))

}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Studentized bootstrap computation of maximum sampling errors in the 
# estimation of the populational mean by means of the sample mean
# x = pilot sample - numerical or numerical factor
# sizes.to.test = sample sizes to be analyzed
# conf.level = confidence level
# B = number of bootstrap iterations
# ------------------------------------------------------------------
bt.errors.for.mean <- function(x, sizes.to.test=length(x), conf.level=0.95, B=1000)
{
	b.copy<-NA
	if (!is.na(b.seed)) 
	{
		set.seed(b.seed)
		warning("Remember that seed is fixed for bootstrap procedures")
		b.copy<-b.seed
		b.seed<<-NA	#Prevent a warning per bt.test execution!
	}
  error<-function(n)
  {
    ret<-bt.test(x, conf.level=conf.level, B=B, proportion=FALSE, boot.sample.size=n)
    return(as.numeric(ret$boot.confidence.interval[4]))
  }
  errors<-data.frame(sample.size=sizes.to.test, maximum.error=sapply(sizes.to.test,error))
	if (!is.na(b.copy)) b.seed<<-b.copy
  return(errors)
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Test and CI for the variance - E((X-mu)^2) - through the bootstrap test and CI for the expectation E((X-mu)^2)
# The sample (X-mu)^2 is approximated by means of (X-mean(X))^2
# x = sample data - numeric or numeric factor
# alternative = type of alternative hypothesis: "two.sided", "less", "greater"
# sigma = standard deviation under the equality of the null hypothesis
# conf.level = confidence level
# B = number of bootstrap iterations
# boot.sample.size = sample size for bootstrap resamplings 
# ------------------------------------------------------------------
bt.vartest<-function(x, alternative='two.sided', sigma=1.0, conf.level=.95, B=1000, boot.sample.size=n)
{
  n<-length(x)
  x<-(n*(x-mean(x))/(n-1))^2
  y<-bt.test(x, alternative=alternative, mu=sigma, conf.level=conf.level, B=B, boot.sample.size=boot.sample.size)
  IC<-sqrt(c(y$boot.confidence.interval[1], y$boot.confidence.interval[2]))
  IC<-c(IC,(IC[2]+IC[1])/2,(IC[2]-IC[1])/2)
  names(IC)<-c("Lim.Inf","Lim.Sup","Mid","Spread")
  y$boot.confidence.interval.sd<-IC
  return(y)
}
# ------------------------------------------------------------------



