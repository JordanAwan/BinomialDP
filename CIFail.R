source("./Tulap.R") 

pvalTulap=function(Tulap,size,theta,lambda,cut){ 
    reps = length(Tulap)
    pval=rep(0,reps)
    values = seq(0,size)

    B = dbinom(values, size=size,prob=theta)

    for(r in 1:reps){
        F = ptulap(values-Tulap[r],lambda,cut)
        pval[r]=t(F)%*%B
    }
    
    return(pval)
}


theta = seq(.01,.99,by=.01)
pvals = rep(0,length(theta))
for(i in 1:length(theta)){
    pvals[i] = pvalTulap(Tulap=5,size=50,theta=theta[i],lambda=lambda,cut=cut)
}
plot(theta,pvals,type="l")





derivPval = function(Tulap,size,theta,lambda,cut){
    if(theta<=0 || theta>=1)
        print("Invalid Theta in derivPval")
    
    reps = length(Tulap)
    deriv=rep(0,reps)
    values = seq(0,size)

    B = dbinom(values, size=size,prob=theta)
    Weight = (values - size*theta)/(theta*(1-theta))

    for(r in 1:reps){
        F = Weight*ptulap(values-Tulap[r],lambda,cut)
        deriv[r]=t(F)%*%B
    }
    return(deriv)
}

theta = seq(.0001,.99,by=.01)
deriv = rep(0,length(theta))
for(i in 1:length(theta)){
    deriv[i] = derivPval(Tulap=1,size=50,theta=theta[i],lambda=lambda,cut=cut)
}
plot(theta,deriv,type="l")




objective = function(V,Tulap,alpha,size,lambda,cut){
    theta1 = exp(V[1])/(1+exp(V[1]))
    theta2 = exp(V[2])/(1+exp(V[2]))
    
    first = (derivPval(Tulap=Tulap,size=size,theta=theta2,lambda=lambda,cut=cut)-
                                            derivPval(Tulap=Tulap,size=size,theta=theta1,lambda=lambda,cut=cut))^2
    #first = (theta1-theta2)^2
    second = 10000*(pvalTulap(Tulap=Tulap,size=size,theta=theta2,lambda=lambda,cut=cut)-
              pvalTulap(Tulap=Tulap,size=size,theta=theta1,lambda=lambda,cut=cut)-(1-alpha))^2
    #print(first)
    #print(second)
    return(first+second)
}

getCI = function(Tulap,alpha,size,lambda,cut){
    initial = c(-1,1)
    answer = optim(par=initial,fn=objective,Tulap=Tulap,alpha=alpha,size=size,lambda=lambda,cut=cut)
    V = answer$par
    print(answer$value)
    CI = exp(V)/(1+exp(V))
    return(CI)
}


n=50
theta_0=.1
ep=1
delta=.01

X = rbinom(n=1,size=n,prob=theta_0)
lambda = exp(-ep)
b=lambda
cut = 2*delta*b/(1-b+2*delta*b)

Z = rtulap(1,median=X,lambda=lambda,cut = cut)
Z
Z/n
alpha=.05

CI=getCI(Tulap=Z,alpha=alpha,size=n,lambda=lambda,cut=cut)
    CI

