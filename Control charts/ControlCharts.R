################################################################
###### Anomaly detection via control charts based on SFPCA #####
################################################################
################### Author: A. Menafoglio ######################
################################################################

rm(list=ls())

############################# LOAD PACKAGES ###############################
library(fda)

############################# LOAD FUNCTIONS ##############################

trapzc = function(t_step,y)
{ # compute the integral of y with step "t_step"
  return(t_step*(0.5*y[1]+sum(y[2:(length(y)-1)]) +0.5*y[length(y)]))
}

clr = function(density, z, z_step)
{ # transform a density to a clr
  return(log(density)-trapzc(z_step,log(density))/(max(z)-min(z)))
}

clr2density <- function(clr, z, z_step)
{ # back-transform a clr to a density
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}

#############################  LOAD DATA  ###############################
######### Example 1
# Gaussian density
# Abscissa
t=seq(-5,5,0.05)

# Data. 
# Sample A: Gaussian densities with mean 0, sd exp(-1+(i-1)/10), i=1,...,21
densities=NULL
for(i in 1:21)
  densities=rbind(densities,
                  dnorm(t,sd=exp(-1+(i-1)/10))/
                    trapzc(t[2]-t[1],dnorm(t,sd=exp(-1+(i-1)/10))))
dA = densities

# Sample B: Gaussian densities with mean 0.5, sd exp(-1+(i-1)/10), i=1,...,21
densities=NULL
for(i in 1:21)
  densities=rbind(densities,
                  dnorm(t,mean = .5, sd=exp(-1+(i-1)/10))/
                    trapzc(t[2]-t[1],dnorm(t,mean = .5, sd=exp(-1+(i-1)/10))))

dB = densities

matplot(t,t(dA),type="l",lty=1,xlab="",lwd=2,
        ylab="Gaussian density",cex.lab=1.2,col=grey(0.1+(0:21)/30), main='Samples A and B')
matlines(t,t(dB),type="l",xlab="",lwd=2,col=grey(0.1+(0:21)/30), lty=2)

# ######### Example 2
# # Metal foam densities (modified from Menafoglio et. al (Techn. 2018))
# # Abscissa
# a=-15
# b=0
# t = seq(a,b,length=302)
# t = t[-c(1,302)]
# 
# dA = t(read.table("dA.txt"))
# dB = t(read.table("dB.txt"))
# 
# matplot(t,t(dA),type="l",lty=1,xlab="",lwd=2,ylim=c(0,0.25),
#         ylab="Density",cex.lab=1.2,col=grey(0.1+(0:21)/30), main='Samples A and B')
# matlines(t,t(dB),type="l",xlab="",lwd=2,col=2, lty=1)

######################## CLR TRANSFORMATION ##############################
t_step = t[2]-t[1]
n=length(t)
clr.dA = lapply(as.list(data.frame(t(dA))), clr, t, t_step)
clr.dB = lapply(as.list(data.frame(t(dB))), clr, t, t_step)

clr.dA=simplify2array(clr.dA)
clr.dB=simplify2array(clr.dB)

matplot(t, clr.dA, type='l', ylim=c(-100,50), col=1, lty=1)
matlines(t, clr.dB, type='l', col=2, lty=1)

######################### SIMPLICIAL FPCA ###############################

# Compute the mean
m1=apply(clr.dA, 1, mean)
m2=apply(clr.dB, 1, mean)

n1=dim(clr.dA)[2]

# Numerical way to get SFPCs
lmod=clr.dA
lmod=lmod*sqrt(t_step)
lmod[1,]=lmod[1,]/sqrt(2)
lmod[n,]=lmod[n,]/sqrt(2)
SS=(n1-1)/n1*cov(t(lmod))
S=SS

pc=eigen(S)
pc$vec[1,]=pc$vec[1,]*sqrt(2)
pc$vec[n,]=pc$vec[n,]*sqrt(2)
pc$vec=pc$vec/sqrt(t_step)

K=Nmax.harm=10

# Compute the centred observations
clr.dA.c=NULL
for(i in 1:n1)
  clr.dA.c=cbind(clr.dA.c, clr.dA[,i]-m1)

# Compute the scores of the sample A
sc.a=matrix(NA, ncol=Nmax.harm, nrow=n1)
for(i in 1:n1)
{
  for(j in 1:Nmax.harm)
    sc.a[i,j]=trapzc(t_step, clr.dA.c[,i]*pc$vec[,j])
}

# Choose the truncation order K: look at the boxplots of the scores and the scree plot 
boxplot(sc.a, las=1, col='gold', main='Principal Components')
plot(1:Nmax.harm, cumsum(pc$val)[1:Nmax.harm]/sum(pc$val), type='b', pch=20, ylim=c(0,1),
     xlab='K', ylab='% Variance')
abline(h=.95, lty=2, col='grey')
abline(h=.98, lty=2, col='grey')
# elbow at K=1 (recall: we are in a 1D exponential family!)
graphics.off()

# or set a threshold to the % of explained variance
thr=.98
K=min(which((cumsum(pc$val)[1:Nmax.harm]/sum(pc$val)>thr)))

# Projected curves  
# Sample A
clr.dA.p=NULL
for(i in 1:n1)
{
  tmp=m1
  for(j in 1:K)
    tmp=tmp+sc.a[i,j]*(pc$vec[,j])
  clr.dA.p=cbind(clr.dA.p, (tmp))
}
dA.p=NULL
for(i in 1:n1)
{
  dA.p=cbind(dA.p,clr2density(clr.dA.p[,i],t,t_step))
}

# Graphics
matplot(t, dA.p, type='l', col=grey(0.1+(0:n1)/30), lty=1, main='Projected sample A')
matplot(t, dA.p, type='l', col=grey(0.1+(0:n1)/30), lty=1, main='Original sample A')

pca=list()
pca$values=pc$values
pca$harmonics=pc$vec
pca$scores=sc.a

######################### PROJECT SAMPLE B ###############################
n2=dim(dB)[1]
clr.dB.c=NULL
for(i in 1:n2)
  clr.dB.c=cbind(clr.dB.c, clr.dB[,i]-m1)

sc.aB=matrix(NA, ncol=Nmax.harm, nrow=n2)
for(i in 1:n2)
{
  for(j in 1:Nmax.harm)
    sc.aB[i,j]=trapzc(t_step, clr.dB.c[,i]*pc$vec[,j])
}

clr.dB.p=NULL
for(i in 1:n2)
{
  tmp=m1
  for(j in 1:K)
    tmp=tmp+sc.aB[i,j]*(pc$vec[,j])
  clr.dB.p=cbind(clr.dB.p, (tmp))
}

dB.p=NULL
for(i in 1:n1)
{
  dB.p=cbind(dB.p,clr2density(clr.dB.p[,i],t,t_step))
}

## Graphics
matplot(t, dB.p, type='l', col=grey(0.1+(0:n1)/30), lty=1, 
        ylab ='Density',main='Projected sample B')
matplot(t, t(dB), type='l', col=grey(0.1+(0:n1)/30), lty=1, 
        ylab='Density',main='Original sample B')

plot(x=sc.a[,1], 
     y=rep(0, n1), 
     pch=19, 
     col=rep(1,n1), ylim=c(0,5), 
     xlab='Scores along SFPC 1', ylab='')
points(x=sc.aB[,1], 
     y=rep(0.5, n2), 
     pch=19, 
     col=rep(2,n2), ylim=c(0,5))

boxplot(sc.a[,1], sc.aB[,1], col=c('grey', 'red'),names = c('Sample A', 'Sample B'))

# Look at other components (if K>1)
# k=2
# plot(x=sc.a[,k], 
#      y=rep(0, n1), 
#      pch=19, 
#      col=rep(1,n1), ylim=c(0,5), 
#      xlab='Scores along SFPC 1', ylab='')
# points(x=sc.aB[,k], 
#        y=rep(0.5, n2), 
#        pch=19, 
#        col=rep(2,n2), ylim=c(0,5))
# boxplot(sc.a[,k], sc.aB[,k], col=c('grey', 'red'),names = c('Sample A', 'Sample B'))

graphics.off()

######################### RESIDUALS ###############################

resA=clr.dA-clr.dA.p
mseA=rep(0,n1)
for(i in 1:n1)
  mseA[i]=trapzc(t_step,resA[,i]^2)

resB=clr.dB-clr.dB.p
mseB=rep(0,n2)
for(i in 1:n2)
{
  mseB[i]=trapzc(t_step,resB[,i]^2)
}

boxplot(mseA, mseB, col=c('grey', 'red'),names = c('Sample A', 'Sample B'))

###################### T2 STATISTIC ##############################

# scores of SFPCA
T2 = matrix(0, 1, (n1+n2))
SCO = rbind(sc.a, sc.aB)
for (i in 1:(n1+n2))
{
  for (k in 1:K)
  {
    T2[i]=T2[i]+(SCO[i,k]^2)/pc$values[k]
  }
}

# Plot of control statistics
# T2
plot(1:(n1+n2),T2, pch=19, cex=.5, type='b', xlab='Sample ID')
abline(v=n2, col='grey', lwd=2)
abline(h=quantile(T2[1:n1], 0.99), col=2, lwd=2)

# MSE
plot(1:(n1+n2),c(mseA, mseB), pch=19, cex=.5, type='b', xlab='Sample ID', ylab='MSE')
abline(v=n2, col='grey', lwd=2)
abline(h=quantile(mseA, 0.99), col=2, lwd=2)

graphics.off()
