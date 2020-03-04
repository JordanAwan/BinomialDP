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

simplePval = function(Tulap,size,theta,lambda,cut){
    return(2*pmin(pvalTulap(Tulap,size,theta,lambda,cut),
                  1-pvalTulap(Tulap,size,theta,lambda,cut)))
}



nearlyUnbiasedPval = function(Tulap,size,theta,lambda,cut){
    #reps = length(Tulap)
    #pval=rep(0,reps)
    #values = seq(0,size)

    T = abs(Tulap-size*theta)

    return(pvalTulap(Tulap=T+size*theta,
                     size=size,theta=theta,lambda=lambda,cut=cut) +
           1-
           pvalTulap(Tulap=size*theta-T,
                     size=size,theta=theta,lambda=lambda,cut=cut))
}


oneSide = function(theta,n,alpha,epsilon,delta){
    b=exp(-epsilon)
    q=2*delta*b/(1-b+2*delta*b)
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
    b=exp(-epsilon)
    q=2*delta*b/(1-b+2*delta*b)
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
    b=exp(-epsilon)
    q=2*delta*b/(1-b+2*delta*b)
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
              alpha=alpha,Tulap=Z,size=size,lambda=b,cut=q,
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



CItwoSideSimple2 = function(alpha,Tulap,size,lambda,cut){
    mle = Tulap/size
    mle = max(min(mle,1),0)
    
    obj = function(theta){
        return(sapply(X=theta,FUN=simplePval,Tulap=Tulap,size=size,lambda=lambda,cut=cut)-alpha)
    }

    if(obj(mle)<0){
        return(c(0,1))

    }
    
    else{
        if(obj(0)<0 & mle>0){
            L = uniroot(f=obj,interval=c(0,mle))
                                        #print(result)
            L=L$root
        }
        else
            L=0
        if(obj(1)<0 & mle<1){

            U = uniroot(f=obj,interval = c(mle,1))
                                        #}
            U = U$root
        }
        else
            U=1
        
        CI = c(L,U)
        return(CI)
    }
    

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


CINearlyUnbiased2= function(alpha,Tulap,size,lambda,cut){
    mle = Tulap/size
    mle = max(min(mle,1),0)
    
    obj = function(theta){
        return(sapply(X=theta,FUN=nearlyUnbiasedPval,Tulap=Tulap,size=size,lambda=lambda,cut=cut)-alpha)
    }

    if(obj(mle)<0){
        return(c(0,1))

    }
    
    else{
        if(obj(0)<0 & mle>0){
            L = uniroot(f=obj,interval=c(0,mle))
                                        #print(result)
            L=L$root
        }
        else
            L=0
        if(obj(1)<0 & mle<1){

            U = uniroot(f=obj,interval = c(mle,1))
                                        #}
            U = U$root
        }
        else
            U=1
        
        CI = c(L,U)
        return(CI)
     }
    
}



############################################################
###   FIGURE 3
############################################################
n=30
#truth=.4
null=.1
#reps=10000
epsilon=.1
delta=0#.01
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

    phiSimple = #pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)#,1)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideTest_1.pdf",width=7,height=7)
plot(thetas,powerUMPU,type="l",col="red",ylim=c(0,1),ylab="power",lwd=2,xlab="theta")
lines(thetas,powerNearly,type = "l", col="purple",lty=2,lwd=2)
lines(thetas,powerSimple,lty=3,lwd=2)
lines(thetas,powerLeft,col="darkgreen",lty=4,lwd=2)
lines(thetas,powerRight,col="blue",lty=6,lwd=2)
legend(.3,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU", "Bonferroni"),col=c("darkgreen","blue","red","purple","black"),lty=c(4,6,1,2,3),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()



############################################################
###   FIGURE 3.5
############################################################
n=30
#truth=.4
null=.1
#reps=10000
epsilon=.1
delta=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
numThetas=101
thetas = seq(0,.2,length=numThetas)

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

    phiSimple = #pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)#,1)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideTest_12.pdf",width=7,height=7)
plot(thetas,powerUMPU,type="l",col="red",ylim=c(0.03,.08),ylab="power",lwd=2,xlab="theta")
lines(thetas,powerNearly,type = "l", col="purple",lty=2,lwd=2)
lines(thetas,powerSimple,lty=3,lwd=2)
lines(thetas,powerLeft,col="darkgreen",lty=4,lwd=2)
lines(thetas,powerRight,col="blue",lty=6,lwd=2)
legend(.05,.08,c("UMP Left","UMP Right","UMPU","Approx UMPU", "Bonferroni"),col=c("darkgreen","blue","red","purple","black"),lty=c(4,6,1,2,3),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()



############################################################
###   FIGURE 4
############################################################
n=100
#truth=.4
null=.5
#reps=10000
epsilon=.1
delta=0#.01
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

    phiSimple = #pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)#,1)
        

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
legend(.3,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU", "Bonferroni"),col=c("darkgreen","blue","red","purple","black"),lty=c(4,6,1,2,3),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()



############################################################
###   FIGURE 5
############################################################
n=100
truth=.7
null=.6
#reps=10000
ep_vect = c(.01,.1,1,10)
numEp = length(ep_vect)
delta=0#.01
#b=exp(-ep)
#q = 2*de*b/(1-b+2*de*b)
alpha=.05
#numThetas=101
#thetas = seq(0,1,length=numThetas)

powerUMPU = rep(0,numEp)
powerNearly = rep(0,numEp)
powerLeft = rep(0,numEp)
powerRight = rep(0,numEp)
powerSimple = rep(0,numEp)
powerUpper = rep(0,numEp)

tic()
for(t in 1:length(ep_vect)){
    #truth = thetas[t]
    epsilon = ep_vect[t]
    #print(t*100/numThetas)

    phiUMPU = UMPU(theta=null, n=n, alpha=alpha, epsilon=epsilon, delta=delta)
    phiNearly = nearlyUMPU(theta=null, n=n, alpha=alpha, epsilon=epsilon, delta=delta)
    phiLeft = oneSide(theta=null, n=n, alpha=alpha, epsilon=epsilon, delta=delta)
    phiRight = 1-oneSide(theta=null, n=n, alpha=1-alpha, epsilon=epsilon, delta=delta)

    phiSimple = pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta),1)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    #print(dbinom(seq(0,n),size=n,prob=null)%*%phiNearly)
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideEpsilon.pdf",width=7,height=7)
plot(powerUMPU,type="l",col="red",ylim=c(0,1),ylab="power",lwd=2,xaxt="n",xlab="epsilon")
axis(1, labels=ep_vect, at = 1:numEp)
lines(powerNearly,type = "l", col="purple",lty=2,lwd=2)
lines(powerSimple,lty=3,lwd=2)
lines(powerLeft,col="darkgreen",lty=4,lwd=2)
lines(powerRight,col="blue",lty=6,lwd=2)
legend(1,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU", "Bonferroni"),col=c("darkgreen","blue","red","purple","black"),lty=c(4,6,1,2,3),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()
############################################################
###   FIGURE 6
############################################################
#n=100
truth=.75
null=.8
nVect = 2^(4:12)
numN = length(nVect)

#reps=10000
epsilon=.1
de=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
#numThetas=101
#thetas = seq(0,1,length=numThetas)

powerUMPU = rep(0,numN)
powerNearly = rep(0,numN)
powerLeft = rep(0,numN)
powerRight = rep(0,numN)
powerSimple = rep(0,numN)
powerUpper = rep(0,numN)

tic()
for(t in 1:numN){
    n = nVect[t]
    #print(t*100/numThetas)

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

pdf(file="TwoSideN.pdf",width=7,height=7)
plot(powerUMPU,type="l",col="red",ylim=c(0,1),ylab="power",lwd=2,xaxt="n",xlab="n")
axis(1, labels=nVect, at = 1:numN)
lines(powerNearly,type = "l", col="purple",lty=2,lwd=2)
lines(powerSimple,lty=3,lwd=2)
lines(powerLeft,col="darkgreen",lty=4,lwd=2)
lines(powerRight,col="blue",lty=6,lwd=2)
legend(1,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU", "Bonferroni"),col=c("darkgreen","blue","red","purple","black"),lty=c(4,6,1,2,3),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()


############################################################
###   CI FIGURE ONE
############################################################
set.seed(100)
n=30
#truth=.4
#null=.4
reps=1000
ep=1
de=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
numThetas=21
thetas = seq(1/(numThetas),1-1/(numThetas),length=numThetas)

 
tic()
coverage1 = coverage2 = rep(0,numThetas)


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

        if(truth>=CI1[1] & truth<=CI1[2])
            coverage1[t] = coverage1[t] + 1/reps
        #CI0 = CItwoSideSimple(Tulap=Z,size=n,alpha=.05*2,lambda=b,cut=q)
        #width0[t] = width0[t] + (CI0[2]-CI0[1])/reps
        
        CI2=CINearlyUnbiased(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        width2[t] = width2[t] + ( CI2[2]-CI2[1])/reps

        if(truth>=CI2[1] & truth<=CI2[2])
            coverage2[t] = coverage2[t] + 1/reps

    }
}


toc()
pdf(file="CI_30_1.pdf",width=7,height=7)
plot(thetas,(width1),type="l",col="blue",xlab="Theta",ylab="Average Width",ylim=c(.1,.4))
lines(thetas,(width2),col="red",lty=2)
legend(.2,.2,c("Approx UMPU", "Bonferroni"),col=c("red","blue"),lty=c(2,1))
#lines(thetas,width0,col="black",lty=3)
dev.off()
summary(width1)
summary(width2)

width2/width1
width1/width2

coverage1
coverage2
plot(coverage1)
points(coverage2,pch="+",col="red")
summary(c(coverage1,coverage2))
summary(coverage2)

sqrt(.95*(1-.95)/(1000))
hist(coverage1)
lines(hist(coverage2))


#plot(thetas,width2/width1,type="l")


############################################################
###   CI FIGURE TWO
############################################################
###   over epsilon??
n=100
truth=.5
#null=.4
reps=1000
#ep=1
ep_vect = c(.01,.1,1,10)

numEp=length(ep_vect)
#de=0#.01
#b=exp(-ep)
#q = 2*de*b/(1-b+2*de*b)
alpha=.05
#numThetas=11
thetas = seq(0,1,length=numThetas)

 
tic()
coverage1 = coverage2 = 0


#width0 = rep(0,numThetas)
width1 = rep(0,numEp)
width2 = rep(0,numEp)
for(t in 1:numEp){
    print(t)
    ep = ep_vect[t]
    de=0#.01
    b=exp(-ep)
    q = 2*de*b/(1-b+2*de*b)
    #print(t*100/numThetas)
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
#pdf(file="CI_30_1.pdf",width=7,height=7)
plot((width1),type="l",col="blue",xlab="Theta",ylab="Average Width",ylim=c(0,1),xaxt="n")
axis(1, labels=ep_vect, at = 1:numEp)
lines((width2),col="red",lty=2)
legend(.2,.2,c("Approx UMPU", "Bonferroni"),col=c("red","blue"),lty=c(2,1))
#lines(thetas,width0,col="black",lty=3)
#dev.off()
#summary(width1)
#summary(width2)


plot(width2/width1,type="l")
summary(width2/width1)


############################################################
###   CI FIGURE THREE
############################################################
###   over n??
#n=100
truth=.5
#null=.4
reps=1000
#ep=1
n_vect = 2^(2:7)
ep=1
numN = length(n_vect)
#numEp=length(ep_vect)
de=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
#numThetas=11
#thetas = seq(0,1,length=numThetas)

 
tic()
coverage1 = coverage2 = rep(0,numN)


#width0 = rep(0,numThetas)
width1 = rep(0,numN)
width2 = rep(0,numN)
for(t in 1:numN){
    print(t)
    n = n_vect[t]
    #ep = ep_vect[t]
    #de=0#.01
    #b=exp(-ep)
    #q = 2*de*b/(1-b+2*de*b)
    #print(t*100/numThetas)
    for(r in 1:reps){
       # if(r%%1000==0)
       #     print((r*100)/reps)
        X = rbinom(n=1,size=n,prob=truth)
        Z = X+rtulap(n=1,median=0,lambda=b,cut=q)

        #CI0=CIoneSide(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        #width0[t] = width0[t] + (1-CI0)/reps
        
        CI1=CItwoSideSimple2(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        width1[t] = width1[t] + (CI1[2]-CI1[1])/reps

        if(truth<=CI1[2] & truth>=CI1[1])
            coverage1[t] = coverage1[t]+1/reps

        #CI0 = CItwoSideSimple(Tulap=Z,size=n,alpha=.05*2,lambda=b,cut=q)
        #width0[t] = width0[t] + (CI0[2]-CI0[1])/reps
        
        CI2=CINearlyUnbiased2(Tulap=Z,size=n,alpha=.05,lambda=b,cut=q)
        width2[t] = width2[t] + (CI2[2]-CI2[1])/reps
         if(truth<=CI2[2] & truth>=CI2[1])
            coverage2[t] = coverage2[t]+1/reps

    }
}


toc()
pdf(file="CI_N.pdf",width=7,height=7)
plot((width1),type="l",col="blue",xlab="n",ylab="Average Width",ylim=c(0,1),xaxt="n")
axis(1, labels=n_vect, at = 1:numN)
lines((width2),col="red",lty=2)
legend(2,1,c("Approx UMPU", "Bonferroni"),col=c("red","blue"),lty=c(2,1))
#lines(thetas,width0,col="black",lty=3)
dev.off()
#summary(width1)
                                        #summary(width2)
summary(width2/width1)
which(width2/width1<.9754)
n_vect[3]
coverage1
coverage2


#plot(thetas,width2/width1,type="l")




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


############################################################

############################################################
theta = .3
n=100
alpha=.05
epsilon=1
delta=0
test = UMPU(theta,n,alpha,epsilon,delta)
values = seq(0,n)
pdf("phi1_0.pdf",width=5,height=5)
plot(values, test,type="l",xlab="x",ylab="phi")
dev.off()
#B = dbinom(values, size=n,prob=theta)
#points(values, B)
#points(values,test)




############################################################
############################################################
############################################################
###   MODIFIED PLOTS
############################################################
############################################################

############################################################
###   FIGURE 3
############################################################
n=30
#truth=.4
null=.1
#reps=10000
epsilon=.1
delta=0#.01
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

    phiSimple = #pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)#,1)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideTest_1Mod.pdf",width=7,height=7)
plot(thetas,powerUMPU,type="l",col="red",ylim=c(0,1),ylab="power",lwd=2,xlab="theta")
lines(thetas,powerNearly,type = "l", col="purple",lty=2,lwd=2)
#lines(thetas,powerSimple,lty=3,lwd=2)
lines(thetas,powerLeft,col="darkgreen",lty=4,lwd=2)
lines(thetas,powerRight,col="blue",lty=6,lwd=2)
legend(.3,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU"),col=c("darkgreen","blue","red","purple"),lty=c(4,6,1,2),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()



############################################################
###   FIGURE 3.5
############################################################
n=30
#truth=.4
null=.1
#reps=10000
epsilon=.1
delta=0#.01
b=exp(-ep)
q = 2*de*b/(1-b+2*de*b)
alpha=.05
numThetas=101
thetas = seq(0,.2,length=numThetas)

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

    phiSimple = #pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)#,1)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideTest_12Mod.pdf",width=7,height=7)
plot(thetas,powerUMPU,type="l",col="red",ylim=c(0.03,.08),ylab="power",lwd=2,xlab="theta")
lines(thetas,powerNearly,type = "l", col="purple",lty=2,lwd=2)
#lines(thetas,powerSimple,lty=3,lwd=2)
lines(thetas,powerLeft,col="darkgreen",lty=4,lwd=2)
lines(thetas,powerRight,col="blue",lty=6,lwd=2)
legend(.05,.08,c("UMP Left","UMP Right","UMPU","Approx UMPU"),col=c("darkgreen","blue","red","purple"),lty=c(4,6,1,2),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()



############################################################
###   FIGURE 4
############################################################
n=100
#truth=.4
null=.5
#reps=10000
epsilon=.1
delta=0#.01
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

    phiSimple = #pmin(
        oneSide(theta=null, n=n, alpha=alpha/2, epsilon=epsilon, delta=delta)+
        1-oneSide(theta=(null), n=n, alpha=1-alpha/2, epsilon=epsilon, delta=delta)#,1)
        

    powerUMPU[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiUMPU
    powerNearly[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiNearly
    powerLeft[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiLeft
    powerRight[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiRight
    powerSimple[t] = dbinom(seq(0,n),size=n,prob=truth)%*%phiSimple
}

toc()

pdf(file="TwoSideTest_5Mod.pdf",width=7,height=7)
plot(thetas,powerUMPU,type="l",col="red",ylim=c(0,1),ylab="power",lwd=2)
lines(thetas,powerNearly,type = "l", col="purple",lty=2,lwd=2)
#lines(thetas,powerSimple,lty=3,lwd=2)
lines(thetas,powerLeft,col="darkgreen",lty=4,lwd=2)
lines(thetas,powerRight,col="blue",lty=6,lwd=2)
legend(.3,.9,c("UMP Left","UMP Right","UMPU","Approx UMPU"),col=c("darkgreen","blue","red","purple"),lty=c(4,6,1,2),lwd=2)
#lines(thetas,powerRight,col="green")
#lines(thetas,powerUpper,lty=3,col="blue")
dev.off()

