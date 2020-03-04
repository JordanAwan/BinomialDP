source("./Tulap.R")
library(tictoc)
#library(rootSolve)

pvalTulap=function(Tulap,size,theta,lambda,cut){ 
    reps = length(Tulap)
    pval=rep(0,reps)
    values = seq(0,size)

    B = dbinom(values, size=size,prob=theta)

    for(r in 1:reps){
        F = ptulap(t=values-Tulap[r],median=0,lambda,cut)
        pval[r]=t(F)%*%B
    }
    
    return(pval)
}


nearlyUnbiasedPval = function(Tulap,size,theta,lambda,cut){
    reps = length(Tulap)
    pval=rep(0,reps)
    values = seq(0,size)

    T = abs(Tulap-n*theta)

    return(pvalTulap(Tulap=T+n*theta,size,theta,lambda,cut) + 1-pvalTulap(Tulap=n*theta-T,size,theta,lambda,cut))
}


oneSide = function(theta,n,alpha,epsilon,delta){
    b=exp(-ep)
    q=2*de*b/(1-b+2*de*b)
    values = seq(0,n)
    B = dbinom(values,size=n,prob=theta)
    #BX = B*(values - n*theta)
    obj = function(m){
        #k = km[1]
        #m = km[2]
        #greaterK = (values>=k)
        phi = ptulap(t=values -m,median=0,lambda=b,cut=q)
        #F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
        #phi = F1*greaterK + F2*(1-greaterK)

        return( (B%*%phi - alpha))    
    }

    lower=-1
    while(obj(lower)<0)
        lower=lower*2
    upper=1
    while(obj(upper)>0)
        upper=upper*2
    #result = optim(par = 0,
    #               fn=obj,method="BFGS")
    result = uniroot(f=obj,interval=c(lower,upper))                                   #print(result)
    #if(result$convergence!=0)
    #    print("UMPU not Converge")
    #k = result$par[1]
    m = result$root
    phi = ptulap(t=values -m,median=0,lambda=b,cut=q)

    return(phi)
}



UMPU = function(theta,n,alpha,epsilon,delta){
    b=exp(-ep)
    q=2*de*b/(1-b+2*de*b)
    values = seq(0,n)
    B = dbinom(values,size=n,prob=theta)
    BX = B*(values - n*theta)
    m=0
    
    ###   Search over k for unbiasedness
    obj = function(k){
        greaterK = (values>=k)

        ###   Search over m for alpha level
        miniObj = function(m){
            F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
            F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
            phi = F1*greaterK + F2*(1-greaterK)

            return(B%*%phi - alpha)
        }###  end miniObj
        
        lower=-1
        while(miniObj(lower)<0)
            lower=lower*2
        upper=1
        while(miniObj(upper)>0)
            upper=upper*2
        
        miniResult = uniroot(f = miniObj,interval = c(lower,upper))
        m = miniResult$root
        #print(m)
        #print(miniObj(m))

        F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
        F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
        phi = F1*greaterK + F2*(1-greaterK)

        return((BX%*%phi))    
    }
    #print(m)
    
    lower = -n
    while(obj(lower)<0)
            lower=lower*2
    upper = 2*n
    while(obj(upper)>0)
        upper = upper*2

    result = uniroot(f = obj,interval = c(lower,upper))
    k = result$root
    greaterK = (values>=k)

    miniObj = function(m){
            F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
            F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
            phi = F1*greaterK + F2*(1-greaterK)

            return(B%*%phi - alpha)
        }###  end miniObj
        
        lower=-1
        while(miniObj(lower)<0)
            lower=lower*2
        upper=1
        while(miniObj(upper)>0)
            upper=upper*2
        
        miniResult = uniroot(f = miniObj,interval = c(lower,upper))
        m = miniResult$root

    F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
    F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
    phi = F1*greaterK + F2*(1-greaterK)
    return(phi)
}


nearlyUMPU =  function(theta,n,alpha,epsilon,delta){
    b=exp(-ep)
    q=2*de*b/(1-b+2*de*b)
    values = seq(0,n)
    B = dbinom(values,size=n,prob=theta)
    BX = B*(values - n*theta)
    k = n*theta
    greaterK = (values>=k)

    obj = function(m){
        #k = n*theta
        #m = km[2]
        F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
        F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
        phi = F1*greaterK + F2*(1-greaterK)

        return((B%*%phi - alpha))    
    }
     lower=-1
    while(obj(lower)<0)
        lower=lower*2
    upper=1
    while(obj(upper)>0)
        upper=upper*2
    #result = optim(par = 0,
    #               fn=obj,method="BFGS")
    result = uniroot(f=obj,interval=c(lower,upper))                                   #print(result)

    m = result$root

    F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
    F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
    phi = F1*greaterK + F2*(1-greaterK)
    return(phi)
}


############################################################

CIobj = function(theta,alpha,Tulap,size,lambda,cut){
    return((pvalTulap(Tulap=Tulap,size=size,theta=theta,lambda=lambda,cut=cut)-alpha)^2)
}


CIoneSide = function(alpha,Tulap,size,lambda,cut){
    L = optim(par=.5,fn=CIobj,
              alpha=alpha,Tulap=Z,size=n,lambda=b,cut=q,
              method="Brent",lower=0,upper=1)

    #U = optim(par=.5,fn=CIobj,alpha=1-alpha/2,Tulap=Z,size=n,lambda=b,cut=q,method="Brent",lower=0,upper=1)
                                        # CI = c(L$par,U$par)
 return(L$par)
}

CItwoSideSimple = function(alpha,Tulap,size,lambda,cut){
    L = optim(par=.5,fn=CIobj,
              alpha=alpha/2,Tulap=Tulap,size=size,lambda=lambda,cut=cut,
              method="Brent",lower=0,upper=1)

    U = optim(par=.5,fn=CIobj,
              alpha=1-alpha/2,Tulap=Tulap,size=size,lambda=lambda,cut=cut,
              method="Brent",lower=0,upper=1)
    CI = c(L$par,U$par)
 return(CI)
}


CIobj2 = function(theta,alpha,Tulap,size,lambda,cut){
    return((nearlyUnbiasedPval(Tulap=Tulap,theta=theta,size=size,lambda=lambda,cut=cut)-alpha)^2)
}

CINearlyUnbiased = function(alpha,Tulap,size,lambda,cut){
    mle = Tulap/size
    mle = max(min(mle,1),0)
    

    if(mle>0){
    L = optim(par=mle/2,fn=CIobj2,
              alpha=alpha,Tulap=Tulap,size=size,lambda=lambda,cut=cut,
              method="Brent",lower=0,upper=mle)
    L = L$par
    }
    else
        L=0
    if(mle<1){
    U = optim(par=(1-mle)/2,fn=CIobj2,
              alpha=alpha,Tulap=Tulap,size=size,lambda=lambda,cut=cut,
              method="Brent",lower=mle,upper=1)
    U = U$par
    }
    else
        U=1
    
    CI = c(L,U)
    return(CI)
        
}




############################################################
###
############################################################
n=30
#truth=.4
null=.5
reps=10000
ep=.1
de=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
numThetas=101
thetas = seq(0,1,length=numThetas)

powerUMPU = rep(0,numThetas)
powerNearly = rep(0,numThetas)
powerLeft = rep(0,numThetas)
powerRight = rep(0,numThetas)
powerSimple = rep(0,numThetas)
powerUpper = rep(0,numThetas)

tic()
for(t in 1:numThetas){
    truth = thetas[t]
    print(t*100/numThetas)

    phiUMPU = UMPU(theta=null, n=n, alpha=alpha, epsilon=epsilon, delta=delta)
    phiNearly = nearlyUMPU(theta=null, n=n, alpha=alpha, epsilon=epsilon, delta=delta)
    phiLeft = oneSide(theta=null, n=n, alpha=alpha, epsilon=epsilon, delta=delta)
    phiRight = 1-oneSide(theta=null, n=n, alpha=1-alpha, epsilon=epsilon, delta=delta)

    phiSimple = #pmax(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideTest_5.pdf",width=7,height=7)
plot(thetas,powerUMPU,type="l",col="red",ylim=c(0,1),ylab="power",lwd=2)
lines(thetas,powerNearly,type = "l", col="purple",lty=2,lwd=2)
lines(thetas,powerSimple,lty=3,lwd=2)
lines(thetas,powerLeft,col="darkgreen",lty=4,lwd=2)
lines(thetas,powerRight,col="blue",lty=6,lwd=2)
legend(.3,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU", "Naive"),col=c("darkgreen","blue","red","purple","black"),lty=c(4,6,1,2,3),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()



############################################################
############################################################
n=30
#truth=.4
#null=.4
reps=1000
ep=1
de=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
numThetas=11
thetas = seq(0,1,length=numThetas)

 
tic()
coverage1 = coverage2 = 0


#width0 = rep(0,numThetas)
width1 = rep(0,numThetas)
width2 = rep(0,numThetas)
for(t in 1:numThetas){
    truth = thetas[t]
    print(t*100/numThetas)
    for(r in 1:reps){
       # if(r%%1000==0)
       #     print((r*100)/reps)
        X = rbinom(n=1,size=n,prob=truth)
        Z = X+rtulap(n=1,median=0,lambda=b,cut=q)

        #CI0=CIoneSide(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        #width0[t] = width0[t] + (1-CI0)/reps
        
        CI1=CItwoSideSimple(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        width1[t] = width1[t] + (CI1[2]-CI1[1])/reps

        #CI0 = CItwoSideSimple(Tulap=Z,size=n,alpha=.05*2,lambda=b,cut=q)
        #width0[t] = width0[t] + (CI0[2]-CI0[1])/reps
        
        CI2=CINearlyUnbiased(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        width2[t] = width2[t] + ( CI2[2]-CI2[1])/reps

    }
}


toc()
pdf(file="CI_30_1.pdf",width=7,height=7)
plot(thetas,(width1),type="l",col="blue",xlab="Theta",ylab="Average Width",ylim=c(.1,.4))
lines(thetas,(width2),col="red",lty=2)
legend(.2,.2,c("Approx UMPU", "Naive"),col=c("red","blue"),lty=c(2,1))
#lines(thetas,width0,col="black",lty=3)
dev.off()
summary(width1)
summary(width2)


plot(thetas,width2/width1,type="l")


############################################################
Check the bias.
############################################################

n=300
#truth=.4
null=.05
reps=100000
ep=1
de=0.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
numThetas=11

tic()
bias=0



for(r in 1:reps){
    if(r%%1000==0)
        print(r*100/reps)
    X = rbinom(n=1,size=n,prob=null)
    Z = X+rtulap(n=1,median=0,lambda=b,cut=q)

    if(nearlyUnbiasedPval(Tulap=Z,size=n,theta=null,lambda=b,cut=q)<=alpha)
        bias = bias + (X-n*null)/reps
}
bias/n
toc()
