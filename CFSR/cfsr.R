#########################################################
###### Compositional function-on-scalar regression ######
#########################################################
###### Authors: R. Talska, A. Menafoglio, I. Pavlu ######
#########################################################

###### settings, packages, functions ######
setwd()

library(fda)
library(robCompositions)
source("SmoothingSpline.R")
source("NulBaze.R")

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
response = read.csv2("response_f.csv", header = FALSE, sep = ";", dec=",")#15x12
x_ax=seq(0.1,0.85,length=12)
#matplot(x_ax,t(response),pch=16)

x=c(0.08,x_ax)
t.fine = seq(min(x), max(x), length=1000)
t.step = diff(t.fine[1:2])

### widths + volumes of intervals in histogram
width=as.matrix(diff(x)) #12

centers=NULL
for (i in 1:length(x)-1){
  centers[i]=x[i]+(x[i+1]-x[i])/2
}
centers=as.matrix(centers) #12x1

### volume standardization - sum to 1
norm.hist=matrix(nrow=nrow(response),ncol=ncol(response))
for (j in 1:ncol(response)){
  for(i in 1:nrow(response)){
    norm.hist[i,j]=response[i,j]/sum(response[i,]) #proportion
  }
}#15x12
apply(norm.hist,1,sum)

densities=matrix(nrow=nrow(response),ncol=ncol(response))
for (j in 1:ncol(densities)){
  for(i in 1:nrow(densities)){
    densities[i,j]=norm.hist[i,j]/width[j]
  }
}#15x12

###### clr transformation ######
response.clr = cenLR((densities))$x.clr

###### graphs I  ######
### clr transformed densities + (original) densities
options(scipen = 1)
par(mfcol=c(1,2))

matplot(centers,t(response.clr),lty=1:length(centers), type="l",xlab = "x",ylab="clr density - not smoothed",col="darkblue",pch=16)
abline(h=0,col="darkred")
matplot(centers,t(densities),lty=1:length(centers), type="l",xlab ="x",ylab="density- not smoothed",col="darkblue")
dev.off()




###### smoothing of clr densities using smoothing splines ######
#arguments for SmoothingSpline
knots=c(min(t.fine),0.15,0.3,0.7,max(t.fine)) #arbitrary - minimization of functional J
w = rep(1,ncol(response.clr)) 
k = 4
der = 2
alfa = 0.99
ch = 1     
t =c(t(centers))


###### SmoothingSpline ######
#generating z-coeficients for B-spline basis
#one observation
Spline1 = SmoothingSpline0(knots=knots, t=t, f=as.numeric(response.clr[1,]), w=w, k=k, der=der, alfa=alfa, ch=ch) #18
abline(v=knots,col="gray",lty=2)

J=c()
z_coef=matrix(nrow=length(knots)+1,ncol=nrow(response.clr))
for (i in 1:nrow(response.clr)){
  z_coef[,i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(response.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[2]]
  J[i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(response.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[1]]
}
sum(J) 
z_coef=as.matrix(z_coef)
dev.off()

#______________________________________________________________________________________
### USE OF PREPROCCESED DATA

clr.raw = response.clr
clr.raw = as.matrix(clr.raw)

dens.raw = densities
dens.raw = as.matrix(dens.raw)

t.raw = read.table("predictor_s.csv", header = FALSE, sep = ";", dec=",")
t.raw = as.matrix(t.raw)


###### construction of B-spline basis ######
# ZB-spline basis 
Z = ZsplineBasis(knots = knots,k)$C0

# clr density
data.l = Z%*%(z_coef) #1000x80

# density in B2
N = nrow(response.clr) 
data = NULL
for (i in 1:N)
{
  data = cbind(data, clr2density(t.fine, t.step, data.l[,i])) #1000*80
}

par(mfcol=c(1,2))
matplot(t.fine,data.l, 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="clr density - smoothed",xlab = "x")
lines(t.fine,t(apply(data.l,1,mean)),col="darkblue",lwd=2)
abline(v=knots,col="lightgray",lty=2)
abline(h=0,col="darkred",lty=1)

matplot(t.fine,data, 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="density - smoothed",xlab = "x")
lines(t.fine,t(apply(data,1,mean)),col="darkblue",lwd=2)
abline(v=knots,col="lightgray",lty=2)
dev.off()

#################################
# compositional splines through B-spline basis
b_coef = t(ZsplineBasis(knots = knots,k)$D)%*%ZsplineBasis(knots = knots,k)$K%*%z_coef
splajn.basis = create.bspline.basis(range(knots), nbasis = dim(z_coef)[1]+1,norder=k,breaks=knots)
fd.data = fd(b_coef,splajn.basis)


#______________________________________________________________________________________
### Regression modeling via B-spline 

# Matrix Y: 
# dim(Y) = N x (g+k+1)
# N observations (B-spline coef) y_i (length g+k+1 = 6); g=3, k=2

Y=as.matrix(t(b_coef))
# design matrix X
# dim(X)=N x 2
x1 = rep(1,N)
x2 = read.csv2("predictor_s.csv", header = FALSE, sep = ";", dec=",")
x2=as.vector(t(x2))
X = as.matrix(cbind(x1,x2))

# Estimate matrix B with LSE 
# dim(B) = 2 x (g+k+1)
B = solve(t(X)%*%X)%*%t(X)%*%Y

# functional estimates of functions beta0, beta1
beta.fd = fd(t(B), splajn.basis)

beta.l = NULL
for (i in 1:2)
{
  beta.l = cbind(beta.l, eval.fd(t.fine, beta.fd)[,i])
}

beta = NULL
for (i in 1:2)
{
  beta = cbind(beta, clr2density(t.fine, t.step,beta.l[,i]))
}

# Compute estimation of B matrix by using command lm
B.reg = lm(Y ~ X[,2])
summary(B.reg)
coef(B.reg)

# plot of the results in clr space and back transformation in Bayes space
#x11(width=21, height=6)
par(mfrow=c(1,3))
plot(beta.fd, xlab="x", ylab="clr(density)", col=c("indianred3","indianred3"), lw=2, cex.lab=1.4, las=1)
legend("topright", lty=c(1,2), col = c("indianred3", "indianred3"),legend = c(expression(paste(clr,"[",hat(beta)[0](t),"]")), expression(paste(clr,"[",hat(beta)[1](t),"]")) ),
       cex = 1.2)

matplot(t.fine,beta[,1],type='l',xlab="x", ylab="density", col="indianred3", lty=1, lw=2,cex.lab=1.4,las=1)
legend("topright", lty=1, col = "indianred3",legend = expression(paste(hat(beta)[0](t))),  cex = 1.2)

matplot(t.fine,beta[,2],type='l',xlab="x", ylab="density", col="indianred3", lty=2, lw=2,cex.lab=1.4,las=1)
legend("topleft", lty=2, col = "indianred3",legend = expression(paste(hat(beta)[1](t))),  cex = 1.2)
dev.off()
# Fit B-spline coef of the response (matrix Y)
Y.est = X %*% B

# functional estimations of functions y(t)
y.fd = fd(t(Y.est), splajn.basis)

y.l = NULL
for (i in 1:N)
{
  y.l = cbind(y.l, eval.fd(t.fine,y.fd)[,i])
}

y = NULL
for (i in 1:N)
{
  y = cbind(y, clr2density(t.fine, t.step, y.l[,i]))
}

#----------------------------------------------------------------------------
### Goodness-of-fit
SSE = rep(0,1000) # var.residuals
SST = rep(0,1000) # var.total
SSF = rep(0,1000) # var.fit

# var.total = var.fit + var.residuals
mean = apply(y.l, 1, mean)

for (i in 1:N)
{
  SSE = SSE + (data.l[,i] - y.l[,i])^2
  SST = SST + (data.l[,i] - mean)^2
  SSF = SSF + (y.l[,i] - mean)^2
}

# Pointwise coefficient of determination
R.t = SSF/SST

# global coefficient of determination
SST.norm = c()
SSF.norm = c()

for (i in 1:N)
{
  SST.norm[i] = trapzc(t.step,(data.l[,i] - mean)^2)
  SSF.norm[i] = trapzc(t.step,(y.l[,i] - mean)^2)
}

R = sum(SSF.norm)/sum(SST.norm) 
R

#prediction & diagnostics
#x11(width=21,height=6)
par(mfrow=c(1,3))
matplot(t.fine,data.l,lty=1,type='l', col="grey", xlab="x", ylab="clr(density)", ylim=c(-3,2),cex.lab=1.4,las=1)
for (i in 1:N){
  matlines(t.fine,y.l[,i], col="indianred3", lty=2)
}
abline(0,0, lty=2)
legend("topright", lty=c(1,2), col = c("grey","indianred3"),legend = c(expression(paste("observed ", y[i](t))),expression(paste("fitted ", hat(y)[i](t)))) ,cex=1)

matplot(t.fine,data,lty=1,type='l', col="grey", xlab="x", ylab="density",cex.lab=1.4,las=1)
for (i in 1:N){
  matlines(t.fine,y[,i], col="indianred3", lty=2)
}
legend("topright", lty=c(1,2), col = c("grey","indianred3"),legend = c(expression(paste("observed ", y[i](t))),expression(paste("fitted ", hat(y)[i](t)))) )

matplot(t.fine,R.t, type='l',xlab="x",
        ylab="Coefficient of determination" , ylim=c(0,1),cex.lab=1.4, las=1)
dev.off()
#----------------------------------------------------------------------------------------------
# Bootstrap
# functional residuals
residua = NULL
for (i in 1:N)
{
  residua = cbind(residua, data.l[,i] - y.l[,i])
}
mean_residua = apply(residua, 1, mean)


matplot(t.fine,residua,col=rainbow(N), type="l" , xlab="x",las=1,ylab="residuals")
matlines(t.fine,mean_residua, col="black",type="l")

# compute bootstrap response Y_boot, R = 200 bootstrap repetitions
weight = x2
weight = rep(weight,200) 
sample_boot = sample(1:15, 3000, replace=T) # sampling with repetition

Y_boot = NULL
for (i in 1:length(sample_boot)){
  Y_boot = cbind(Y_boot,beta.l[,1]+beta.l[,2]*weight[i]+residua[,sample_boot[i]])
}
dim(Y_boot)


matplot(t.fine,Y_boot, type="l",las=1,xlab="x")

# Make a list of Y_boot: each list consist of 15 (=N) Y_boot
d2 = list()
step = 15 
for (i in 1:200){
  if (i==1) {d2[[i]]=Y_boot[,1:step]}
  if (i!=1) {d2[[i]]=Y_boot[,(1+(i-1)*step):((i)*step)]}
}  
dim(d2[[1]])

# Hat matrix
H = solve(t(X)%*%X)%*%t(X)

# Compute bootstrap estimates for betas
Beta.boot.l = NULL
for (i in 1:200)
{
  Beta.boot.l = rbind(Beta.boot.l,H%*%t(d2[[i]]))
}
dim(Beta.boot.l)


matplot(t.fine,t(Beta.boot.l),type="l", xlab="x", ylab="clr(density)", main="beta0, beta1 - bootstrap")
matlines(t.fine,beta.l[,1],type="l", color=1, lwd=4 )
matlines(t.fine,beta.l[,2],type="l", color=1, lwd=4 )

matplot(t.fine,t(Beta.boot.l),type="l", xlab="x", ylab="clr(density)", col="grey",cex.lab=1.4,las=1)
matlines(t.fine,beta.l[,1],type="l", col="black", lwd=1,lty=2 )
matlines(t.fine,beta.l[,2],type="l", col="black", lwd=1)
legend("bottomleft", lty=c(2,1), col = c("black","black"),legend = c(expression(paste(clr,"[",hat(beta)[0](t),"]")),expression(paste(clr,"[",hat(beta)[1](t),"]")) ), cex=1.3)

# Computed beta bootstrap estimates after inverse of clr transformation
Beta.boot = NULL
for (i in 1:400)
{
  Beta.boot = rbind(Beta.boot,clr2density(t.fine,t.step,Beta.boot.l[i,]))
}

matplot(t.fine,t(Beta.boot),type="l", xlab="x", ylab="density")

# Boot-estimates for beta1
Beta.boot2 = NULL
for (i in 1:200)
{
  Beta.boot2 = cbind(Beta.boot2, clr2density(t.fine, t.step, Beta.boot.l[2*i,]))
}

# Boot-estimates for beta0
Beta.boot1 = NULL
for (i in 1:200)
{
  Beta.boot1 = cbind(Beta.boot1, clr2density(t.fine, t.step, Beta.boot.l[2*i-1,]))
}

par(mfrow=c(1,2))
matplot(t.fine,Beta.boot1,type="l", xlab="x", ylab="density" , col="grey",cex.lab=1.4,las=1)
matlines(t.fine,beta[,1], type="l")
legend("topright", lty=1, col = "black",legend = expression(paste(hat(beta)[0](t))),  cex = 1.2)

matplot(t.fine,Beta.boot2,type="l",xlab="x", ylab="density",cex.lab=1.4, col="grey",las=1)
matlines(t.fine, beta[,2], type="l")
legend("topleft", lty=1, col = "black",legend = expression(paste(hat(beta)[1](t))),  cex = 1.2)

