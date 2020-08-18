rm(list=ls())

# install.packages("genlasso")
# install.packages("glmnet")
library(genlasso)
library(glmnet)

data = read.csv("C:/Users/HSMOON/Dropbox/2019 Hyeseong Moon/Thesis/TrimRaw500.csv",header=F)

# install.packages("NormalBetaPrime")
# library(NormalBetaPrime)
# trim32
xx.mat = as.matrix(data)[,-1]; y.vec = drop(data[,1]) 

### sure independent screening 
n = dim(xx.mat)[1]
# d.vec = rep(0,dim(xx.mat)[2])
# for(j in 1:ncol(xx.mat)){ d.vec[j] = abs(coef(glm(y.vec~.,data=data.frame(y.vec,xx.mat[,j]))))[-1] }
# o.vec = order(d.vec,decreasing=T); 

top = round(n/log(n))
# t.vec = nrow(xx.mat)*c(2:4)
# top = t.vec[2]
x.mat = xx.mat[,1:top]; n = dim(x.mat)[1]; b.mat = NULL

# fused lasso
p = dim(x.mat)[2]; val = NULL
for(i in 1:p){for(j in 1:p){ if((j-i)==1) {val = c(val,i,j)} } }
gr = graph(edges=val,directed=F); D = getDgSparse(gr)
t.fit = genlasso::fusedlasso(y.vec,x.mat,D=D,gamma=1)
lam.vec = t.fit$lambda; mval = NULL
for(i in 1:5){
  ## cross validation
  idx = sample(n,size=round(0.3*n),replace=T)
  tx.mat = x.mat[-idx,]; ty.vec = y.vec[-idx]
  vx.mat = x.mat[idx,]; vy.vec = y.vec[idx]
  fit = genlasso::fusedlasso(ty.vec,tx.mat,D=D,gamma=1)
  mse = NULL
  for(pos in 1:length(lam.vec)){
    vpre = predict(fit,lambda=lam.vec[pos],Xnew=vx.mat)
    mse = c(mse,mean((vy.vec-vpre$fit)^2))
    plot(mse)
  }
  mval = rbind(mval,mse)
}
opt = which.min(colMeans(mval)); lam = lam.vec[opt]
b.mat = cbind(b.mat,coef(t.fit,lambda=lam)$beta)

# sparse fused lasso
val = NULL
for(i in 1:p){for(j in 1:p){ if((j-i)==1) {val = c(val,i,j)} } }
gr = graph(edges=val,directed=F); D = getDgSparse(gr)
gam.vec = seq(0.1,0.9,by=0.1); tb.mat=min.val=NULL
# gam = gam.vec[1]
for(gam in gam.vec){
  t.fit = genlasso::fusedlasso(ty.vec,tx.mat,D=D,gamma=gam)
  lam.vec = t.fit$lambda; mval = NULL
  for(i in 1:5){
    ## cross validation
    idx = sample(n,size=round(0.3*n),replace=T)
    tx.mat = x.mat[-idx,]; ty.vec = y.vec[-idx]
    vx.mat = x.mat[idx,]; vy.vec = y.vec[idx]
    fit = genlasso::fusedlasso(ty.vec,tx.mat,D=D,gamma=gam)
    mse = NULL
    for(pos in 1:length(lam.vec)){
      vpre = predict(fit,lambda=lam.vec[pos],Xnew=vx.mat)
      mse = c(mse,mean((vy.vec-vpre$fit)^2))
      plot(mse)
    }
    mval = rbind(mval,mse)
  }
  lam = lam.vec[which.min(colMeans(mval))]; tb.mat=cbind(tb.mat,coef(t.fit,lambda=lam)$beta) 
  min.val = c(min.val,min(colMeans(mval)));
}
b.mat = cbind(b.mat,tb.mat[,which.min(min.val)])

# pairwise fused lasso
# val = NULL
# for(i in 1:p){ for(j in 1:p){ if(i<j){val = c(val,i,j)} } }
# gr = graph(edges=val,directed=F); D = getDgSparse(gr)
# t.fit = genlasso::fusedlasso(ty.vec,tx.mat,D=D,gamma=1)
# lam.vec = t.fit$lambda; mval = NULL
# for(i in 1:5){
#   ## cross validation
#   idx = sample(n,size=round(0.3*n),replace=T)
#   tx.mat = x.mat[-idx,]; ty.vec = y.vec[-idx]
#   vx.mat = x.mat[idx,]; vy.vec = y.vec[idx]
#   fit = genlasso::fusedlasso(ty.vec,tx.mat,D=D,gamma=1)
#   mse = NULL
#   for(pos in 1:length(lam.vec)){
#     vpre = predict(fit,lambda=lam.vec[pos],Xnew=vx.mat)
#     mse = c(mse,mean((vy.vec-vpre$fit)^2))
#     plot(mse)
#   }
#   mval = rbind(mval,mse)
# }
# opt = which.min(colMeans(mval)); lam = lam.vec[opt]
# b.mat = cbind(b.mat,coef(t.fit,lambda=lam)$beta)

# smmooth lasso
m.vec = 10^(-3:1); c.mat = NULL;  minval = rep(0,length(m.vec))
# m = m.vec[1]
for(mid in 1:length(m.vec)){
  cat("mid : ",mid,"\n")
  p = dim(x.mat)[2]; j.mat = -cbind(diag(p-1),0)+cbind(0,diag(p-1))
  dx.mat = rbind(x.mat,m.vec[mid]*j.mat); dy.vec = c(y.vec,c(rep(0,p-1)))
  t.fit = glmnet(dx.mat,dy.vec)
  lam.vec = t.fit$lambda; mval = NULL
  for(i in 1:5){
    cat("cvid : ",i,"\n")
    ## cross validation
    idx = sample(n,size=round(0.3*n),replace=T)
    tx.mat = x.mat[-idx,]; ty.vec = y.vec[-idx]
    vx.mat = x.mat[idx,]; vy.vec = y.vec[idx]
    dtx.mat = rbind(tx.mat,m.vec[mid]*j.mat); dty.vec = c(ty.vec,c(rep(0,p-1)))
    fit = glmnet(dtx.mat,dty.vec,lambda=lam.vec)
    dvx.mat = rbind(vx.mat,m.vec[mid]*j.mat); dvy.vec = c(vy.vec,c(rep(0,p-1)))
    mse = NULL
    for(pos in 1:length(lam.vec)){
      vpre = predict(fit,newx=dvx.mat,s=lam.vec[pos])
      mse = c(mse,mean((dvy.vec-vpre)^2))
      plot(mse)
    }
    mval = rbind(mval,mse)
  }
  opt = which.min(colMeans(mval))
  c.mat = cbind(c.mat,t.fit$beta[,opt])
  minval[mid] = min(colMeans(mval))
}
b.mat=cbind(b.mat,c.mat[,which.min(minval)])
b.mat
round(b.mat,4)

### igraph
# install.packages("igraph")
library(igraph)
par(mfrow=c(1,3))
for(pos in 1:3){
  ad.mat = matrix(0,25,25)
  for(k in 1:25){
    for(i in k:25){
      if((i+1)>25) break
      if(b.mat[k,pos]==b.mat[(i+1),pos]) {ad.mat[k,(i+1)]=1}
    }
  }
  g1 = graph_from_adjacency_matrix(ad.mat,mode="undirected",weighted=NULL,diag=FALSE)
  plot(g1,layout=layout.fruchterman.reingold, vertex.shape= "none",vertex.label.font = 2,vertex.label.cex=1,
       vertex.label.family ="mono")
}

# plot(g1,layout=layout.circle,vertex.shape="none",vertex.label.cex=.8,vertex.label.dist=0.5)
#layout.kamada.kawai
# save.image("C:\\Users\\HSMOON\\Desktop\\do first\\penalty simulation\\fusing.RData")


# adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)
# g1 <- graph_from_adjacency_matrix( adjm )
# 
# adjm <- matrix(sample(0:5, 100, replace=TRUE,
#                       prob=c(0.9,0.02,0.02,0.02,0.02,0.02)), nc=10)
# g2 <- graph_from_adjacency_matrix(adjm, weighted=TRUE)
# graph.edgelist(adjm)
# graph.edgelist(adjm,directed = F)
# plot(g2)



# install.packages("ADMM")
# library(ADMM)
# p = dim(x.mat)[2]
# val = NULL
# for(i in 1:p){for(j in 1:p){ if((j-i)==1) {val = c(val,i,j)} } }
# gr = graph(edges=val,directed=F)
# D = as.matrix(getDgSparse(gr))
# n = dim(x.mat)[1]
# lam.vec = n*seq(0.1,0.9,length.out = 100); tb.mat = NULL
# for(i in lam.vec){
#   t.fit = admm.genlasso(x.mat,y.vec,D,lambda=i)
#   tb.mat = cbind(b.mat,t.fit$x)
# }
# ## cross validation
# m.mat = NULL
# for(fid in 1:5){
#   idx = sample(n,size=round(0.3*n),replace=T)
#   tx.mat = x.mat[-idx,]; ty.vec = y.vec[-idx]
#   vx.mat = x.mat[idx,]; vy.vec = y.vec[idx]
#   b.mat = NULL
#   for(i in lam.vec){
#     fit = admm.genlasso(tx.mat,ty.vec,D,lambda=i)
#     b.mat = cbind(b.mat,fit$x)
#   }
#   pre = vx.mat%*%b.mat
#   mse = colMeans((vy.vec-pre)^2)
#   plot(mse)
#   m.mat = rbind(m.mat,mse)
# }
# opt = which.min(colMeans(m.mat))
# opt
# lam = lam.vec[opt]
# b.vec = tb.mat[,opt]
# round(b.vec,2)
