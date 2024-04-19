################################################################
###### Simplicial functional principle component analysis ######
################################################################
################### Author: A. Menafoglio ######################
################################################################

rm(list=ls())

###### settings, packages, functions ######
setwd("")

library(robCompositions)
library(fda)
library(fields)
source("NulBaze.R")
source("SmoothingSpline.R")

clr2density <- function(z, z_step, clr) 
{   # Inverse of clr transformation 
    # Input: z = grid of point defining the abscissa 
    #        z_step = step of the grid of the abscissa
    #        clr = grid evaluation of the clr transformed density
    # Output: grid evaluation of the density
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}

trapzc<-function(step,y) 
{   # Numerical integration via trapezoidal formula
    # Input: y = grid evaluation of the function
    #        z_step = step of the grid
  int<-step*(0.5*y[1]+sum(y[2:(length(y)-1)]) + 0.5*y[length(y)])
  return (int)
}

###### preprocessing of the data ######
predictor = read.csv2("predictor_f.csv", header = TRUE, sep = ",", dec=".")#80x60
grain_size = read.csv2("grain_size.csv", header = FALSE, sep = ",", dec=".")#60x1
grain_size=as.matrix(grain_size)


log.grain_size=c(log(2),log(grain_size))
t.fine = seq(min(log.grain_size), max(log.grain_size), length=1000)
t.step = diff(t.fine[1:2])

### widths + volumes of intervals in histogram
width=as.matrix(diff(log.grain_size)) #60x1

centers=NULL
for (i in 1:length(log.grain_size)-1){
  centers[i]=log.grain_size[i]+(log.grain_size[i+1]-log.grain_size[i])/2
}
centers=as.matrix(centers) #60x1

### volume standardization - sum to 1
norm.hist=matrix(nrow=nrow(predictor),ncol=ncol(predictor))
for (j in 1:ncol(predictor)){
  for(i in 1:nrow(predictor)){
    norm.hist[i,j]=predictor[i,j]/sum(predictor[i,]) #proportion
  }
}#80x60
apply(norm.hist,1,sum)

densities=matrix(nrow=nrow(predictor),ncol=ncol(predictor))
for (j in 1:ncol(densities)){
  for(i in 1:nrow(densities)){
    densities[i,j]=norm.hist[i,j]/width[j]
  }
}#80x60

###### clr transformation ######
predictor.clr = cenLR((densities))$x.clr

###### graphs I  ######
### clr transformed densities + (original) densities
options(scipen = 1)
par(mfcol=c(2,1))

matplot(exp(centers),t(predictor.clr),log="x",lty=1:length(centers), type="l",xlab = expression (paste("particle size (",mu,"m)")),ylab="clr density - not smoothed",col="darkblue",pch=16)
abline(h=0,col="darkred")
matplot(exp(centers),t(densities),log="x",lty=1:length(centers), type="l",xlab = expression (paste("particle size (",mu,"m)")),ylab="density- not smoothed",col="darkblue")
dev.off()


###### smoothing of clr densities using smoothing splines ######
#arguments for SmoothingSpline
knots=c(min(t.fine),1.6,2,2.8,4.5,max(t.fine)) #arbitrary - minimization of functional J
w = rep(1,ncol(predictor.clr)) 
k = 4
der = 2
alfa = 0.99
ch = 1     
t =c(t(centers))


par(mfcol=c(2,1))
###### SmoothingSpline ######
#generating z-coeficients for B-spline basis
#one observation
Spline1 = SmoothingSpline0(knots=knots, t=t, f=as.numeric(predictor.clr[1,]), w=w, k=k, der=der, alfa=alfa, ch=ch) #18
abline(v=exp(knots),col="gray",lty=2)

#all observations
J=c()
z_coef=matrix(nrow=length(knots)+1,ncol=nrow(predictor.clr))
for (i in 1:nrow(predictor.clr)){
  z_coef[,i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(predictor.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[2]]
  J[i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(predictor.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[1]]
}
sum(J) 
dev.off()
z_coef=as.matrix(z_coef)

###saving computations
#write.table(predictor.clr, file = "predictor_clr.csv", quote = TRUE, sep = ",", dec = ".")
#write.table(densities, file = "densities.csv", quote = TRUE, sep = ",", dec = ".")
#write.table(centers,file = "centers.csv", quote = TRUE, sep = ",", dec = ".")
#write.table(z_coef,file = "z_coef.csv", quote = TRUE, sep = ",", dec = ".")#7x80

#predictor.clr=read.csv2("predictor_clr.csv", sep = ",", dec = ".") #80x60
#densities=read.csv2("densities.csv", sep = ",", dec = ".") #80x60
#centers=read.csv2("centers.csv", sep = ",", dec = ".") #60x1
#z_coef=read.csv2("z_coef.csv", sep = ",", dec = ".")

###### ZB-spline basis  ######
Z = ZsplineBasis(knots = knots,k)$C0

# clr density
data.l = Z%*%(z_coef) #1000x80

# density in B2
N = nrow(predictor.clr) 
data = NULL
for (i in 1:N)
{
  data = cbind(data, clr2density(t.fine, t.step, data.l[,i])) #1000*80
}

par(mfcol=c(2,1))
matplot(exp(t.fine),data.l,log="x", 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="clr density - smoothed",xlab = expression (paste("log of particle size (",mu,"m)")))
lines(exp(t.fine),t(apply(data.l,1,mean)),col="darkblue",lwd=5)
abline(v=exp(knots),col="lightgray",lty=2)
abline(h=0,col="darkred",lty=1)

matplot(exp(t.fine),data,log="x", 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="density - smoothed",xlab = expression (paste("log of particle size (",mu,"m)")))
lines(exp(t.fine),t(apply(data,1,mean)),col="darkblue",lwd=5)
abline(v=exp(knots),col="lightgray",lty=2)
dev.off()

###### B-spline basis ######
b_coef = t(ZsplineBasis(knots = knots,k)$D)%*%ZsplineBasis(knots = knots,k)$K%*%z_coef
B = create.bspline.basis(range(knots), nbasis = dim(z_coef)[1]+1,norder=k,breaks=knots)
fd.data = fd(b_coef,B) 

par(mfcol=c(2,1))
matplot(exp(t.fine),data.l,log="x", 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="clr density - smoothed",xlab = expression (paste("particle size (",mu,"m)")))
lines(exp(t.fine),t(apply(data.l,1,mean)),col="darkblue",lwd=5)
abline(v=exp(knots),col="lightgray",lty=2)
abline(h=0,col="darkred",lty=1)

matplot(exp(t.fine),data,log="x", 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="density - smoothed",xlab = expression (paste("log of particle size (",mu,"m)")))
lines(exp(t.fine),t(apply(data,1,mean)),col="darkblue",lwd=5)
abline(v=exp(knots),col="lightgray",lty=2)
dev.off()


# Estimation of the mean and of the covariance function
layout(rbind(1,2))
plot.fd(fd.data,xlab="log (particle size)",ylab="clr density")
lines(mean(fd.data),lwd=2)
eval.1 <- eval.fd(t.fine,fd.data)
image.plot(t(cov(t(eval.1))[100:1,]))

##### PCA #####
layout(1)
plot.fd(fd.data,ylab='',xlab="log (particle size)")

# PCA with fda package
pca.l <- pca.fd(fd.data,nharm=5,centerfns=TRUE)

# scree plot
# pca.fd computes all the eigenvalues, only the first 5 are different from zero
layout(rbind(1,2))
plot(pca.l$values[1:6],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca.l$values)[1:6]/sum(pca.l$values),xlab='j',ylab='CPV',ylim=c(0.8,1))
  #cumulative percentage of variability

# first two principal components
layout(rbind(1,2))
plot(pca.l$harmonics[1,],col=1,xlab="log (particle size)",ylab='PC1',lwd=2)
plot(pca.l$harmonics[2,],col=2,xlab="log (particle size)",ylab='PC2',lwd=2)

# plot of the principal components as perturbation of the mean
lmean<-mean(fd.data)

plot(lmean,lwd=2,ylim=c(-2,1.5),ylab='clr(density)',xlab='log (particle size)',main='PC1')
lines(lmean+pca.l$harmonics[1,]*sqrt(pca.l$values[1]), col=2)
lines(lmean-pca.l$harmonics[1,]*sqrt(pca.l$values[1]), col=3)

plot(lmean,lwd=2,ylim=c(-2,1.5),ylab='clr(density)',xlab='log (particle size)',main='PC2')
lines(lmean+pca.l$harmonics[2,]*sqrt(pca.l$values[2]), col=2)
lines(lmean-pca.l$harmonics[2,]*sqrt(pca.l$values[2]), col=3)

dmean=clr2density(t.fine, 1, lmean) 
pca=pca.l
pca$harmonics = harm.std.p = harm.std.m = matrix(0,nrow=length(t.fine), ncol=5)

for(i in 1:5)
{
  pca$harmonics[,i] = clr2density(t.fine,1,pca.l$harmonics[i,])
  # Sign changed according to the choice of representing -harmonics
  harm.std.p[,i] = clr2density(t.fine,1,lmean-pca.l$harmonics[i,]*sqrt(pca.l$values[2]))
  harm.std.m[,i] = clr2density(t.fine,1,lmean+pca.l$harmonics[i,]*sqrt(pca.l$values[2]))
}

dev.off()

#plotting the scores

plot(pca$scores[,1],pca$scores[,2], xlab='PC1', ylab='PC2',main='Scores (SFPCA)')

