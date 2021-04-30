#Zero-integral basis
ZsplineBasis = function(knots,k)
{
  library(fda)
  r = length(knots)
  lambda_index = c(0:(r-1)) 
  g = lambda_index[length(lambda_index) - 1]
  
  N = g+(k-1)+1
  
  lambda = c(rep(min(knots),k-1),knots,rep(max(knots),k-1))
  div = seq(min(lambda), max(lambda), length = 1000) 
  
  # standard B-spline basis; collocation matrix := C
  splajn.basis = create.bspline.basis(range(knots),nbasis = N , norder = k, breaks = knots)
  C = eval.basis(div, splajn.basis)
  
  # Matrix D
  differ = lambda[(1+k):(r+2*(k-1))] - lambda[(1:(r+k-2))]
  D = (k)*diag(1/differ)
  
  # Matrix L
  L = array(0, c(N,N-1))
  L[1,1]=1
  L[N,N-1]=-1
  
  for (j in (2:(N-1))){
    L[j,j-1] = (-1)
    L[j,j] = 1
  }
  
  # Spline0 basis: collocation matrix C0
  C0 = C%*%D%*%L
  
  # Matrix M - function for computing integral
  SLP=function(step, c){
    integral = step*(0.5*c[1]+sum(c[2:(length(c)-1)]) +0.5*c[length(c)])
    return (integral)
  }
  
  division = seq(min(lambda), max(lambda),  length = 10000);step=diff(div[1:2])
  CC = eval.basis(division, splajn.basis)
  
  CC0 = CC%*%D%*%L
  
  M=array(0, c(N-1,N-1))
  for (i in 1:(N-1)){
    for (j in 1:(N-1)){
      non_null = c()
      product = CC0[,i]*CC0[,j]
      for (m in 1:length(div)){
        if (product[m] != 0) {non_null[m] = product[m]}
      }
      M[i,j]=SLP(step, product)
    }
  }
  return(list(C0 = C0, M0 = M, K = L, D = D))
}

#ZsplineBasis(knots,k)
