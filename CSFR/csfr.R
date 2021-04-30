##########################################################
###### Compositional scalar-on-function regression  ######
##########################################################
############ Authors: R. Talska, I. Pavlu ################
##########################################################


rm(list=ls())

###### settings, packages, functions ######
setwd()

library(robCompositions)
library(fda)
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
predictor = read.csv2("predictor_raw.csv", header = TRUE, sep = ",", dec=".")#80x60
grain_size = read.csv2("grain_size.csv", header = FALSE, sep = ",", dec=".")#1x60
grain_size=as.matrix(grain_size)


log.grain_size=c(log(2),log(grain_size))
t.fine = seq(min(log.grain_size), max(log.grain_size), length=1000)
t.step = diff(t.fine[1:2])

### widths + volumes of intervals in histogram
width=as.matrix(diff(log.grain_size)) #60

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
}#80*60
apply(norm.hist,1,sum)

densities=matrix(nrow=nrow(predictor),ncol=ncol(predictor))
for (j in 1:ncol(densities)){
  for(i in 1:nrow(densities)){
    densities[i,j]=norm.hist[i,j]/width[j]
  }
} #80x60

###### clr transformation ######
predictor.clr = cenLR((densities))$x.clr

###### graphs I  ######
### clr transformed densities + (original) densities
options(scipen = 1)
par(mfcol=c(1,2))

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


###### SmoothingSpline ######
#generating z-coeficients for B-spline basis
#one observation
Spline1 = SmoothingSpline0(knots=knots, t=t, f=as.numeric(predictor.clr[1,]), w=w, k=k, der=der, alfa=alfa, ch=ch) #18
abline(v=exp(knots),col="gray",lty=2)

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

###### construction of B-spline basis ######
# ZB-spline basis 
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

par(mfcol=c(1,2))
matplot(exp(t.fine),data.l,log="x", 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="clr density - smoothed",xlab = expression (paste("particle size (",mu,"m)")))
lines(exp(t.fine),t(apply(data.l,1,mean)),col="darkblue",lwd=5)
abline(v=exp(knots),col="lightgray",lty=2)
abline(h=0,col="darkred",lty=1)

matplot(exp(t.fine),data,log="x", 
        lty=1:N, type="l",las=1,cex.lab=1,cex.axis=1.2,col=rainbow(N),
        ylab="density - smoothed",xlab = expression (paste("particle size (",mu,"m)")))
lines(exp(t.fine),t(apply(data,1,mean)),col="darkblue",lwd=5)
abline(v=exp(knots),col="lightgray",lty=2)
dev.off()


#################################

# compositional splines through B-spline basis
b_coef = t(ZsplineBasis(knots = knots,k)$D)%*%ZsplineBasis(knots = knots,k)$K%*%z_coef
B = create.bspline.basis(range(knots), nbasis = dim(z_coef)[1]+1,norder=k,breaks=knots)
fd.data = fd(b_coef,B)

# REGRESSION -------------------------------------------------------------
# response
response = read.csv2("response.csv", dec=",", header=T) #60x5
response.clr = cenLR(response)$x.clr

###### cycle for all parts of (response) composition ######
options(scipen = 1)
par(mfcol=c(1,2))
for (z in 1:dim(response.clr)[2]){
   
  
  response.clr_element = response.clr[,z]
  
  #Regression with SFPCA ---------------------------------------------
  ncom = 3 #3
  pca = pca.fd(fd.data, ncom, centerfns = T)
  attach(pca)
  
  # explained variability
  round(100*(values)[1:ncom]/sum(values),4)
  
  # Design matrix
  X = as.matrix(cbind(rep(1,N),scores))
  
  B.pca = solve(t(X)%*%X)%*%t(X)%*%(response.clr_element)
  B.pca
  
  B0 = B.pca[1,]; B1 = B.pca[-1,]
  
  # Beta1(t) estimates
  mat.proj.l = 0
  mat.proj = c()
  for (j in 1:ncom)
  {
    mat.proj.l  =  mat.proj.l + B1[j]*eval.fd(t.fine, harmonics[j,])
  }
  mat.proj = clr2density(t.fine, t.step, mat.proj.l)
 
  # Coeficient of determination
  yhatPCA = X%*%B.pca
  
  # sums of squares
  reziduaPCA  = response.clr_element-yhatPCA
  SS_totalPCA = response.clr_element-mean(response.clr_element)
  SS_regrePCA = yhatPCA-mean(response.clr_element)
  
  SS_totPCA = sum(SS_totalPCA^2)
  SS_rezPCA = sum(reziduaPCA^2)
  SS_regPCA = sum(SS_regrePCA^2)
  
  R2 = 1-(SS_rezPCA/SS_totPCA);R2
  
  # Confidence bands nased on bootstrap procedure
  smpl = sample(1:N, N*100, replace=T)
  score_input = do.call("rbind", rep(list(scores[,1:ncom]), 100))
  
  # compute bootstrap response y_boot:
  y_boot = c()
  for (R in 1:(N*100)){
    y_boot[R] = B.pca[1,1]+sum(B.pca[2:(ncom+1),1]*score_input[R,1:ncom])+reziduaPCA[smpl[R],1]
  }
  
  # Make a matrix of y_boot: each column consist N y_boot
  y_boot_matrix = matrix(y_boot,100,N,byrow=T)
  
  # Hat matrix
  H = solve(t(X)%*%X)%*%t(X)
  
  # Compute bootstrap estimates for betas
  beta.boot.l = NULL
  for (i in 1:100)
  {
    beta.boot.l = cbind(beta.boot.l,H%*%(y_boot_matrix[i,]))
  }
  beta1.boot.l = beta.boot.l[-1,]
  
  # Bootstrap estimates of beta1(t)
  mat.proj.l.boot = matrix(0,length(t.fine),100)
  mat.proj.boot = NULL
  for (R in 1:100)
  {
    for (j in 1:ncom)
    {
      mat.proj.l.boot[,R]  =  mat.proj.l.boot[,R] + beta1.boot.l[j,R]*eval.fd(t.fine, harmonics[j,])
    }
    mat.proj.boot = cbind(mat.proj.boot , clr2density(t.fine, t.step, mat.proj.l.boot[,R]))
  }
  ### Final plots of functional parameter beta_1 and its bootstrap confidence bands
  #clr
  matplot(exp(t.fine),mat.proj.l.boot,col="gray",
          lty=1, type="l",las=1, lwd=2,cex.lab=1.1,cex.axis=1,
          main = paste(names(response)[z],";" ,ncom," FPCs; ",expression(R2),"=",round(R2,3)),
          ylab="clr(density)",xlab = expression (paste("Particle size (",mu,"m)")),log="x")
  matlines(exp(t.fine),mat.proj.l,lty=1,lwd=3,
           col = "darkblue",log="x")
  abline(h=0,col="darkred")
  #density
  matplot(exp(t.fine),mat.proj.boot,col="gray",
          lty=1, type="l",las=1, lwd=2,cex.lab=1.1,cex.axis=1,
          main = paste(names(response)[z],";" ,ncom," FPCs; ",expression(R2),"=",round(R2,3)),
          xlab = expression (paste("Particle size (",mu,"m)")),
          ylab="density",log="x")
  matlines(exp(t.fine),mat.proj,lty=1,lwd=3,
           col = "darkblue")
  abline(h=0,col="darkred")
  
  detach(pca)
}
dev.off()

### Number of components determined by CV:
ncom_cv=5 
pca_cv = pca.fd(fd.data, ncom_cv, centerfns = T)
attach(pca_cv)
criterion = NULL
for (K in 2:ncom_cv)
{
  X1 = as.matrix(cbind(rep(1,N),scores[,1:K]))
  
  e = NULL; e2 = NULL
  for (q in 1:N)
  {
    response.clr_q = response.clr_element[-q]
    X_q = X1[-q,]
    
    B.pca = solve(t(X_q)%*%X_q)%*%t(X_q)%*%response.clr_q
    yhatPCA = X_q%*%B.pca
    
    # predict i-th observation from validing data
    y_predict = sum(X1[q,]*B.pca)
    
    error_i = response.clr_element[q]-y_predict
    
    e = cbind(e, error_i)
    e2 = cbind(e2, error_i^2)
  }
  criterion = cbind(criterion,(1/N)*sum(e2))
}

plot(2:ncom_cv,criterion,las=1,type="b",pch=20,las=1,lwd=2, ylab="CV(K)", xlab="Number of components")
detach(pca_cv)