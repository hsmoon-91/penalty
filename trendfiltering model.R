# library(mlbench)
# data = Ozone
# str(data)

rm(list=ls())
# data 
library(GeDS) # install.packages("Rmpfr"); install.packages("GeDS") in R
# data("coalMining")
data = coalMining

# penalized methods
# install.packages("mFilter")
library(mFilter)
library(genlasso)
# install.packages("neariso")
library(neariso)
plot(data)
plot(cbind(data[,1],cumsum(data[,2])))

# p = length(data[,1])
# d.mat = outer(1:p,1:p,function(x,y){1*(x>=y)})
# d.mat = outer(1:p,1:p,function(x,y){1*(x==y)}) # identity matrix
# d.mat = -cbind(diag(p-1),0)+cbind(0,diag(p-1)) # mean
# d.mat = outer(1:(p-1),1:p,function(x,y){1*(x==y)+(-1*(x-y==-1))}) # fused LASSO
# d.mat = outer(1:p,1:p,function(x,y){1*(x==y)+(-1*(x-y==1))}) # mean
# d.mat = outer(1:p,1:p,function(x,y){1*(x==y)+(-2*(x-y==-1))+(1*(x-y==-2))}) # trend

par(mfrow=c(2,2))
# out = trendfilter(data[,2],ord=1) # linear trend
# out = fusedlasso1d(data[,2]) # constant trend
# coef(out,lambda = 1.5)$beta
# softthresh(out,lambda=1.5,gamma=1) # sparse fused LASSO : gamma

## cross validation
y = data[,2]; n = length(y)
cv.mat = cbind(y,rep(0,length.out=n))
for(k in c(2,3,4,5,6)){
  for(i in 2:(n-1)){
    if(i %in% seq(k,n,by=5)){cv.mat[i,2] = k-1}
    # if(i %in% seq(3,n,by=5)){cv.mat[i,2] = 2}
    # if(i %in% seq(4,n,by=5)){cv.mat[i,2] = 3}
    # if(i %in% seq(5,n,by=5)){cv.mat[i,2] = 4}
    # if(i %in% seq(6,n,by=5)){cv.mat[i,2] = 5}
  }
}
y.mat = NULL
# true model : constant trend
tout = fusedlasso1d(data[,2])
lam.vec = tout$lambda
err = 0
# selecting tuning parameter
fold=5
for(fid in 1:fold){
  t.vec = cv.mat[cv.mat[,2]==fid,][,1] 
  v.vec = NULL # (1,3), (2,4), (3,5), (4,6), (5,7)
  for(s in seq(fid,n,by=5)){
    if((s+2)>n) break
    a = cv.mat[s,1]; b = cv.mat[(s+2),1]
    v.vec = c(v.vec,mean(a,b))
  }
  out = fusedlasso1d(t.vec)
  b.mat = NULL
  for(i in 1:length(lam.vec)){
    b.mat = cbind(b.mat,coef(out,lambda = lam.vec[i])$beta)
  }
  err = rbind(err,colMeans((v.vec-b.mat)^2))
}
opt = which.min(colMeans(err))
y.mat = cbind(y.mat,coef(tout,lambda=lam.vec[opt])$beta)

# true model : sparse fused LASSO
tout = fusedlasso1d(data[,2])
lam.vec = tout$lambda
err = 0
gam.vec = c(0.1,0.3,0.5,0.7) #?????????????????
# selecting tuning parameter
fold=5
for(fid in 1:fold){
  t.vec = cv.mat[cv.mat[,2]==fid,][,1] 
  v.vec = NULL # (1,3), (2,4), (3,5), (4,6), (5,7)
  for(s in seq(fid,n,by=5)){
    if((s+2)>n) break
    a = cv.mat[s,1]; b = cv.mat[(s+2),1]
    v.vec = c(v.vec,mean(a,b))
  }
  out = fusedlasso1d(t.vec); b.mat = NULL; gam.opt = NULL
  for(gam in gam.vec){ # sparse fused LASSO
    for(i in 1:length(lam.vec)){
      b.mat = cbind(b.mat,softthresh(out,lambda=lam.vec[i],gamma=gam))
      gam.opt = c(gam.opt,gam)
    }
  }
  # plot(colMeans((v.vec-b.mat)^2))
  err = rbind(err,colMeans((v.vec-b.mat)^2))
  # plot(colMeans(err))
}
opt = which.min(colMeans(err))
y.mat = cbind(y.mat,softthresh(tout,lambda=lam.vec[opt],gamma=gam.opt[opt]))

# true model : linear trend
tout = trendfilter(data[,2],ord=1)
lam.vec = tout$lambda
err = 0
# selecting tuning parameter
fold=5
for(fid in 1:fold){
  t.vec = cv.mat[cv.mat[,2]==fid,][,1] 
  v.vec = NULL # (1,3), (2,4), (3,5), (4,6), (5,7)
  for(s in seq(fid,n,by=5)){
    if((s+2)>n) break
    a = cv.mat[s,1]; b = cv.mat[(s+2),1]
    v.vec = c(v.vec,mean(a,b))
  }
  out = trendfilter(t.vec,ord=1) # linear trend
  b.mat = NULL
  for(i in 1:length(lam.vec)){# linear trend, constant trend
    b.mat = cbind(b.mat,coef(out,lambda = lam.vec[i])$beta)
  }
  plot(colMeans((v.vec-b.mat)^2))
  # plot(colMeans((v.vec-b.mat)^2))
  err = rbind(err,colMeans((v.vec-b.mat)^2))
  # plot(colMeans(err))
}
opt = which.min(colMeans(err))
y.mat = cbind(y.mat,coef(tout,lambda=lam.vec[opt])$beta)

# true model : h-p trend
lam.vec = 0.1*(1:100) # ????????????
tb.mat = NULL
for(i in 1:length(lam.vec)){
  tout = hpfilter(data[,2],freq=lam.vec[i],type="lambda",drift=F)
  tb.mat = cbind(tb.mat,tout$trend)
}
fold=5;err=0
for(fid in 1:fold){
  t.vec = cv.mat[cv.mat[,2]==fid,][,1] 
  v.vec = NULL # (1,3), (2,4), (3,5), (4,6), (5,7)
  for(s in seq(fid,n,by=5)){
    if((s+2)>n) break
    a = cv.mat[s,1]; b = cv.mat[(s+2),1]
    v.vec = c(v.vec,mean(a,b))
  }
  b.mat = NULL
  for(i in 1:length(lam.vec)){
    out = hpfilter(t.vec,freq=lam.vec[i],type="lambda",drift=F)
    b.mat = cbind(b.mat,out$trend)
  }
  # plot(colMeans((v.vec-b.mat)^2))
  err = rbind(err,colMeans((v.vec-b.mat)^2))
  plot(colMeans(err))
}
opt = which.min(colMeans(err))
y.mat = cbind(y.mat,tb.mat[,opt])

par(mfrow=c(1,2))
plot(data[,1],data[,2],xlab="Year",ylab="Number of Coal Mine Diasters")
points(data[,1],y.mat[,1],type="l",lty=1,col=1,lwd=2)
points(data[,1],y.mat[,2],type="l",lty=2,col=2,lwd=2)
points(data[,1],y.mat[,3],type="l",lty=3,col=3,lwd=2)
points(data[,1],y.mat[,4],type="l",lty=4,col=4,lwd=2)
legend("topright",c("meanfiltering","sparse fused LASSO","linearfiltering","H-Pfiltering"),lty = c(1,2,3,4),
       cex=0.6,col=c(1,2,3,4),lwd=0.5)

plot(data[,1],cumsum(data[,2]),cex=0.3,xlab="Year",ylab="Cumulative Number of Diasters")
points(data[,1],cumsum(y.mat[,1]),type="l",lty=1,col=1,lwd=2)
points(data[,1],cumsum(y.mat[,2]),type="l",lty=2,col=2,lwd=2)
points(data[,1],cumsum(y.mat[,3]),type="l",lty=3,col=3,lwd=2)
points(data[,1],cumsum(y.mat[,4]),type="l",lty=4,col=4,lwd=2)
legend("topleft",c("meanfiltering","sparse fused LASSO","linearfiltering","H-Pfiltering"),lty = c(1,2,3,4),
       cex=0.6,col=c(1,2,3,4),lwd=0.5)

# true model : neariso : monotone nondecreasing
library(neariso)
lam.vec = 0.1*(100:1)
lam = lam.vec[1]
tout = neariso(data[,2])

###############################
library(ncvreg)
par(mfrow=c(1,1))
n = dim(data)[1]
x.mat=outer(1:n,1:n,function(x,y){1*(x>=y)})
tout = ncvreg(x.mat,data[,2],family="poisson",nlambda=30)
tout = glmnet(x.mat,data[,2],family="poisson",nlambda=30)

lam.vec = tout$lambda
length(lam.vec)
for(fid in 1:5){
  t.vec = cv.mat[cv.mat[,2]==fid,][,1]
  tn = length(t.vec)
  v.vec = NULL # (1,3), (2,4), (3,5), (4,6), (5,7)
  for(s in seq(fid,n,by=5)){
    if((s+2)>n) break
    a = cv.mat[s,1]; b = cv.mat[(s+2),1]
    v.vec = c(v.vec,mean(a,b))
  }
  tx.mat = outer(1:tn,1:tn,function(x,y){1*(x>=y)})
  out = ncvreg(tx.mat,t.vec,family="poisson",lambda=lam.vec)
  length(out$lambda)
  vn = length(v.vec)
  vx.mat = outer(1:vn,1:vn,function(x,y){1*(x>=y)})
  vy.vec=exp(vx.mat%*%(out$beta)[-1,])
  err = colMeans(v.vec-vy.vec)^2
  plot(err) 
}


plot(data[,1],data[,2])

points(data[,1],ey.vec,type="l")

save.image("C:\\Users\\HSMOON\\Desktop\\do first\\penalty simulation\\trend.RData")

########################################################################################################################################3
# cv.trendfilter(object, k = 5, mode = c("lambda", "df"),
#                approx = FALSE, rtol = 1e-07, btol = 1e-07,
#                verbose = FALSE)
# fit = cv.trendfilter(out,k=5,mode="lambda")
# plot(out,lambda=fit$lambda.min,main="Minimal CV error")
# plot(out,lambda=fit$lambda.1se,main="One standard error rule")
# 
# 
# library(neariso)
# # generate some artificial data
# y <- rnorm(1000) + (1:1000)/3000
# ### run the algorithm as default; will output solution at 100 breakpoints for lambda
# res0 <- neariso(y)
# ### apply function nir and get solution directly
# lambda = 0:10/10
# res <- neariso(y, lambda=lambda)
# ### apply the function and get the solution later
# res2 <- neariso(y, lambda=NULL)
# res2 <- nearisoGetSolution(res2, lambda=lambda)
# ### look at the breakpoints
# lambdaBreaks <- nearisoGetBreakpoints(res2, maxBreaks=1000)
# res3 <- nearisoGetSolution(res2, lambda=lambdaBreaks)