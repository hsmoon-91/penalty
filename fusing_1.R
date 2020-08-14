rm(list=ls())
library(MASS)
library(glmnet)
library(genlasso)
# lasso
l.fun = function(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec){
  fit = glmnet(tx.mat,ty.vec); lam.vec = fit$lambda
  vfit = glmnet(vx.mat,vy.vec,lambda=lam.vec)
  mse = NULL
  for(pos in 1:length(lam.vec)){
    vpre = predict(vfit,newx=tx.mat,s=lam.vec[pos])
    mse = c(mse,mean((ty.vec-vpre)^2))
  }
  plot(mse)
  l.coef = fit$beta[,which.min(mse)]
  pre = nx.mat %*% l.coef
  return(list(mse=mean((ny.vec-pre)^2),coef=l.coef))
}
# l.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)

# smooth lasso
s.fun = function(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec){
  m.vec = 10^(-3:1); mval = NULL; b.mat = NULL
  # m = m.vec[1]
  for(m in m.vec){
    p = dim(tx.mat)[2]
    j.mat = -cbind(diag(p-1),0)+cbind(0,diag(p-1))
    
    ttx.mat = rbind(tx.mat,m*j.mat)  
    tty.vec = c(ty.vec,c(rep(0,p-1)))
    
    fit = glmnet(ttx.mat,tty.vec); lam.vec = fit$lambda
    
    vvx.mat = rbind(vx.mat,m*j.mat)  
    vvy.vec = c(vy.vec,c(rep(0,p-1)))
    
    vfit = glmnet(vvx.mat,vvy.vec,lambda=lam.vec)
    
    mse = NULL
    for(pos in 1:length(lam.vec)){
      vpre = predict(vfit,newx=ttx.mat,s=lam.vec[pos])
      mse = c(mse,mean((tty.vec-vpre)^2))
    }
    b.mat = cbind(b.mat,fit$beta)
    mval = c(mval,mse)
  }
  plot(mval)
  coef = b.mat[,which.min(mval)]
  pre = nx.mat %*% coef
  return(list(mse=mean((ny.vec-pre)^2),coef=coef))
}
# s.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)

# elastic net
# tvx.mat=rbind(tx.mat,vx.mat); tvy.vec=c(ty.vec,vy.vec)
e.fun = function(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec){
  alp.vec = round(seq(0.1,0.9,length=10),1); aval=NULL; b.mat=NULL
  for(apos in 1:length(alp.vec)){
    fit = glmnet(tx.mat,ty.vec,alpha=alp.vec[apos])
    lam.vec = fit$lambda; vfit = glmnet(vx.mat,vy.vec,alpha=alp.vec[apos],lambda=lam.vec)
    mse = NULL
    for(pos in 1:length(lam.vec)){
      vpre = predict(vfit,newx=tx.mat,s=lam.vec[pos])
      mse = c(mse,mean((ty.vec-vpre)^2))
    }
    b.mat = cbind(b.mat,fit$beta)
    aval = c(aval,mse)
  }
  plot(aval)
  # length(aval);which.min(aval)
  coef = b.mat[,which.min(aval)]
  pre = nx.mat %*% coef
  return(list(mse=mean((ny.vec-pre)^2),coef=coef))
}
# e.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)

# fused lasso
f.fun = function(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec){
  p = dim(tx.mat)[2]
  val = NULL
  for(i in 1:p){for(j in 1:p){ if((j-i)==1) {val = c(val,i,j)} } }
  gr = graph(edges=val,directed=F)
  fit = genlasso::fusedlasso(ty.vec,tx.mat,graph=gr)
  lam.vec = fit$lambda
  mse = NULL
  for(pos in 1:length(lam.vec)){
    vpre = predict(fit,lambda=lam.vec[pos],Xnew=vx.mat)
    mse = c(mse,mean((vy.vec-vpre$fit)^2))
  }
  plot(mse)
  # length(mse); which.min(mse)
  coef = coef(fit)$beta[,which.min(mse)]
  pre = predict(fit,lambda=lam.vec[which.min(mse)],Xnew=nx.mat)
  return(list(mse=mean((ny.vec-pre$fit)^2),coef=coef))
}
# f.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)

# pairwise
p.fun = function(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec){
  p = dim(tx.mat)[2]
  val = NULL
  for(i in 1:p){ for(j in 1:p){ if(i<j){val = c(val,i,j)} } }
  gr = graph(edges=val,directed=F); # D = getDgSparse(gr); plot(gr)
  fit = genlasso::fusedlasso(ty.vec,tx.mat,graph=gr)
  lam.vec = fit$lambda
  mse = NULL
  for(pos in 1:length(lam.vec)){
    vpre = predict(fit,lambda=lam.vec[pos],Xnew=vx.mat)
    mse = c(mse,mean((vy.vec-vpre$fit)^2))
  }
  pre = predict(fit,lambda=lam.vec[which.min(mse)],Xnew=nx.mat)
  coef = coef(fit)$beta[,which.min(mse)]
  return(list(mse=mean((ny.vec-pre$fit)^2), coef=coef))
}
# p.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)

##################################setting1########################################################
p = 8; b.vec = c(3,1.5,0,0,0,2,0,0)
iter.max = 50
mse1.mat = NULL; l.vec=e.vec=s.vec=p.vec=f.vec=NULL
for(iter in 1:iter.max){
  print(iter)
  corr.mat=outer(1:p, 1:p, function(x,y){.9^abs(x-y)}); diag(corr.mat) = 3; mu=rep(0,p)
  tn = 100
  tx.mat = mvrnorm(tn, mu=mu, Sigma=corr.mat)
  ty.vec = tx.mat%*%b.vec + rnorm(tn)
  vn = 100
  vx.mat = mvrnorm(vn, mu=mu, Sigma=corr.mat)
  vy.vec = vx.mat%*%b.vec + rnorm(vn)
  nn = 200
  nx.mat = mvrnorm(nn, mu=mu, Sigma=corr.mat)
  ny.vec = nx.mat%*%b.vec + rnorm(nn)
  
  lasso = l.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  elastic = e.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  smooth = s.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  fused = f.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  pairwise = p.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  
  mse1.mat = rbind(mse1.mat, c(lasso$mse,elastic$mse,smooth$mse,fused$mse,pairwise$mse))
  l.vec = colMeans(rbind(l.vec,lasso$coef))
  e.vec = colMeans(rbind(e.vec,elastic$coef))
  s.vec = colMeans(rbind(s.vec,smooth$coef))
  f.vec = colMeans(rbind(f.vec,fused$coef))
  p.vec = colMeans(rbind(p.vec,pairwise$coef))
}
coef1.mat = rbind(l.vec,e.vec,s.vec,f.vec,p.vec)
coef1.mat; b.vec
boxplot(mse1.mat)
##################################setting2########################################################
p=20; b.vec = c(rep(0,5),rep(2,5),rep(0,5),rep(2,5))
mse2.mat = NULL; l.vec=e.vec=s.vec=p.vec=f.vec=NULL
for(iter in 1:iter.max){
  print(iter)
  corr.mat=outer(1:p, 1:p, function(x,y){.5^abs(x-y)}); diag(corr.mat) = 15; mu=rep(0,p)
  tn = 200
  tx.mat = mvrnorm(tn, mu=mu, Sigma=corr.mat)
  ty.vec = tx.mat%*%b.vec + rnorm(tn)
  vn = 200
  vx.mat = mvrnorm(vn, mu=mu, Sigma=corr.mat)
  vy.vec = vx.mat%*%b.vec + rnorm(vn)
  nn = 400
  nx.mat = mvrnorm(nn, mu=mu, Sigma=corr.mat)
  ny.vec = nx.mat%*%b.vec + rnorm(nn)
 
  lasso = l.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  elastic = e.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  smooth = s.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  fused = f.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  pairwise = p.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  
  mse2.mat = rbind(mse2.mat, c(lasso$mse,elastic$mse,smooth$mse,fused$mse,pairwise$mse))
  l.vec = colMeans(rbind(l.vec,lasso$coef))
  e.vec = colMeans(rbind(e.vec,elastic$coef))
  s.vec = colMeans(rbind(s.vec,smooth$coef))
  f.vec = colMeans(rbind(f.vec,fused$coef))
  p.vec = colMeans(rbind(p.vec,pairwise$coef))
}
coef2.mat = rbind(l.vec,e.vec,s.vec,f.vec,p.vec)
coef2.mat; b.vec
boxplot(mse2.mat)
##################################setting3########################################################
p = 20; b.vec = c(rep(5,3),rep(2,3),rep(10,3),rep(0,11))
mse3.mat = NULL; l.vec=e.vec=s.vec=p.vec=f.vec=NULL
for(iter in 1:iter.max){
  print(iter)
  tn = 200
  tx.mat=matrix(rep(0,tn*p),nrow=tn,ncol=p)
  Z1 = Z2 = Z3 = rnorm(tn,0,1)
  for(pos in 1:3){tx.mat[,pos]=Z1+rnorm(tn,0,0.01)}
  for(pos in 4:6){tx.mat[,pos]=Z2+rnorm(tn,0,0.01)}
  for(pos in 7:8){tx.mat[,pos]=Z3+rnorm(tn,0,0.01)}
  for(pos in 9:20){tx.mat[,pos]=rnorm(tn,0,1)}
  diag(tx.mat) = 15
  ty.vec = tx.mat%*%b.vec + rnorm(tn)
  
  vn = 200
  vx.mat=matrix(rep(0,vn*p),nrow=vn,ncol=p)
  Z1 = Z2 = Z3 = rnorm(vn,0,1)
  for(pos in 1:3){vx.mat[,pos]=Z1+rnorm(vn,0,0.01)}
  for(pos in 4:6){vx.mat[,pos]=Z2+rnorm(vn,0,0.01)}
  for(pos in 7:8){vx.mat[,pos]=Z3+rnorm(vn,0,0.01)}
  for(pos in 9:20){vx.mat[,pos]=rnorm(vn,0,1)}
  diag(vx.mat) = 15
  vy.vec = vx.mat%*%b.vec + rnorm(vn)
  
  nn = 400
  nx.mat=matrix(rep(0,nn*p),nrow=nn,ncol=p)
  Z1 = Z2 = Z3 = rnorm(nn,0,1)
  for(pos in 1:3){nx.mat[,pos]=Z1+rnorm(nn,0,0.01)}
  for(pos in 4:6){nx.mat[,pos]=Z2+rnorm(nn,0,0.01)}
  for(pos in 7:8){nx.mat[,pos]=Z3+rnorm(nn,0,0.01)}
  for(pos in 9:20){nx.mat[,pos]=rnorm(nn,0,1)}
  diag(nx.mat) = 15
  ny.vec = nx.mat%*%b.vec + rnorm(nn)
  
  lasso = l.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  elastic = e.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  smooth = s.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  fused = f.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  pairwise = p.fun(tx.mat,ty.vec,vx.mat,vy.vec,nx.mat,ny.vec)
  
  mse3.mat = rbind(mse3.mat, c(lasso$mse,elastic$mse,smooth$mse,fused$mse,pairwise$mse))
  l.vec = colMeans(rbind(l.vec,lasso$coef))
  e.vec = colMeans(rbind(e.vec,elastic$coef))
  s.vec = colMeans(rbind(s.vec,smooth$coef))
  f.vec = colMeans(rbind(f.vec,fused$coef))
  p.vec = colMeans(rbind(p.vec,pairwise$coef))
}
coef3.mat = rbind(l.vec,e.vec,s.vec,f.vec,p.vec)
rownames(coef3.mat) = c("lasso","elastic","smooth","fused","pairwise")
round(coef3.mat,2)
coef3.mat; b.vec
boxplot(mse3.mat)


#real data
library(MASS)
library(dplyr)
library(spatstat)
data(birthwt)
dim(birthwt); str(birthwt)
birthwt = birthwt %>% mutate(race=as.factor(race),smoke=as.factor(smoke),ptl=as.factor(ptl),
                             ht=as.factor(ht),ui=as.factor(ui),ftv=as.factor(ftv),low=as.factor(low),
                             age=as.numeric(age),lwt=as.numeric(lwt),bwt=as.numeric(bwt))
df = birthwt
df$race = sapply(df$race,function(x) as.factor(ifelse(x==1,"white",ifelse(x==2,"black","other"))))
df$smoke = sapply(df$smoke,function(x) as.factor(ifelse(x==1,"smoke","nonsmoke")))
df$ptl = sapply(df$ptl,function(x) as.factor(ifelse(x==0,"0",ifelse(x==1,"1","more")))) 
df$ht = sapply(df$ht,function(x) as.factor(ifelse(x==0,"no","yes"))) 
df$ui = sapply(df$ui,function(x) as.factor(ifelse(x==0,"no","yes")))
df$ftv = sapply(df$ftv,function(x) as.factor(ifelse(x==0,"0",ifelse(x==1,"1",ifelse(x==2,"2","over"))))) 
x.mat = df[,2:9] 
x.mat.num = x.mat[,1:2]
x.mat.fac = x.mat[,3:8]
x.mat.dum = dummify(x.mat.fac)
colnames(x.mat.dum)
x.mat = cbind(x.mat.num,x.mat.dum[,-c(2,4,6,9,12,13)]) # baseline
y.vec = df[,10] #numeric
mse4.mat = NULL; l.vec=e.vec=s.vec=p.vec=f.vec=NULL
for(iter in 1:iter.max){
  print(iter)
  idx = sample(c(1,2,3),size=length(y.vec), prob=c(.6,.2,.2),replace=T)
  # train
  tr.y.vec = drop(y.vec[idx==1]); tr.x.mat = as.matrix(x.mat[idx==1,])
  # valid
  vd.y.vec = drop(y.vec[idx==2]); vd.x.mat = as.matrix(x.mat[idx==2,])
  # test
  tn.y.vec = drop(y.vec[idx==3]); tn.x.mat = as.matrix(x.mat[idx==3,])
  
  
  lasso = l.fun(tr.x.mat,tr.y.vec,vd.x.mat,vd.y.vec,tn.x.mat,tn.y.vec)
  elastic = e.fun(tr.x.mat,tr.y.vec,vd.x.mat,vd.y.vec,tn.x.mat,tn.y.vec)
  smooth = s.fun(tr.x.mat,tr.y.vec,vd.x.mat,vd.y.vec,tn.x.mat,tn.y.vec)
  fused = f.fun(tr.x.mat,tr.y.vec,vd.x.mat,vd.y.vec,tn.x.mat,tn.y.vec)
  pairwise = p.fun(tr.x.mat,tr.y.vec,vd.x.mat,vd.y.vec,tn.x.mat,tn.y.vec)
  
  mse4.mat = rbind(mse4.mat, c(log(lasso$mse),log(elastic$mse),log(smooth$mse),log(fused$mse),log(pairwise$mse)))
  l.vec = colMeans(rbind(l.vec,lasso$coef))
  e.vec = colMeans(rbind(e.vec,elastic$coef))
  s.vec = colMeans(rbind(s.vec,smooth$coef))
  f.vec = colMeans(rbind(f.vec,fused$coef))
  p.vec = colMeans(rbind(p.vec,pairwise$coef))
}


real.mse = c(lasso$mse,elastic$mse,smooth$mse,fused$mse,pairwise$mse)
real.vec = rbind(lasso$coef,elastic$coef,smooth$coef,fused$coef,pairwise$coef)

save.image("C:/Users/HSMOON/Desktop/seminar/fusing_1.RData")
load("C:/Users/HSMOON/Desktop/seminar/fusing_1.RData")
