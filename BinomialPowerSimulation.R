library(rmutil) 
library(NormalLaplace)
library(DiscreteLaplace)
set.seed(100)
FU = function(u){# domain is -1/2 to 1/2
    return(u+1/2)
}


ptulap = function(t,b){
    cdf = function(t,b){
        ifelse(t<=0,b^(-round(t))/(1+b)*(b+FU(t-round(t))*(1-b)),1-b^(round(t))/(1+b)*(b+FU(round(t)-t)*(1-b)))
    }
    return(sapply(t,cdf,b=b))
}

C = function(M,b){
    return(b^(round(M))*(b+FU(round(M)-M)*(1-b))/(1+b))
}













############################################################
###   SIMULATION
############################################################
epvect = 2^seq(-6,1,by=1)
size=30

null=.8
alt=.9
reps=10000



powerTulap=powerNN=powerNP=rep(0,length(epvect))
avgpTulap=avgpNN=avgpNP=rep(length(epvect))


for(i in 1:length(epvect)){
    print(i)
    ep=epvect[i]

    ###DATA GENERATION
    X=rbinom(n=reps,size=size,prob=alt)

    Z = X+runif(n=reps,min=-1/2,max=1/2)

    pvalNP = rep(0,reps)
    values = seq(0,size)
    B = dbinom(values, size=size,prob=null)
    for(r in 1:reps){
        F = punif(values-Z[r],min=-1/2,max=1/2)
        pvalNP[r] = t(F)%*%B
    }
    avgpNP[i] = mean(pvalNP)
    powerNP[i] = mean(pvalNP<.05)
    


    ###TULAP
    U = runif(n=reps,min=-1/2,max=1/2)
    G1 = rgeom(n=reps,prob=1-exp(-ep))
    G2 = rgeom(n=reps,prob=1-exp(-ep))
    Tulap = X+U+G1-G2

    ###LAPLACE
    E1 = rexp(n=reps,rate=ep)
    E2 = rexp(n=reps,rate=ep)
    Laplace = X+E1-E2

    ###TULAP   Pvalue
    #values=seq(0,size)
    #B = dbinom(values,size=size,prob=null)
    pvalTulap=rep(0,reps)
    for(r in 1:reps){
        F = ptulap(values-Tulap[r],exp(-ep))
        pvalTulap[r]=t(F)%*%B
    }
    avgpTulap[i] = mean(pvalTulap)
    powerTulap[i] = mean(pvalTulap<.05)

###   Playing with pvals
   # plot(ecdf(pvalTulap))
   # values=seq(0,1,by=.1)
   # lines(values,values)


    
    
    

    
    ###NORMAL LAPLCE Pvalue
    #pvalNL = (1-pnl(Laplace,param = c(null*size,sqrt(size*null*(1-null)),ep,ep)))
    #powerNL[i] = mean(pvalNL<.05)

    
    ###NORMAL NORMAL Pvalue
    pvalNN = (1-pnorm(Laplace,m=null*size,s=sqrt(size*null*(1-null)+2/ep^2)))
    avgpNN[i] = mean(pvalNN)
    powerNN[i] = mean(pvalNN<.05)

    #lines(ecdf(pvalNN),col="blue")
    


    
}

#plot(avgpTulap/avgpNN,type="l")
#plot(powerTulap/powerNN,type="l")

#pdf("PValue.pdf",width=10,height=4)
plot((avgpTulap),col="red",pch=3,ylab="average p-value",xaxt="n",xlab="epsilon")
lines((avgpTulap),col="red")
axis(1,at=seq(1,length(epvect)),labels=epvect)
points((avgpNN),col="blue")
lines((avgpNN),col="blue",lty=2)
lines((avgpNP),col="black",lty=1)
legend(5,.4,c("UMP","Normal Approximation"),
       pch=c(3,1),lty=c(1,2),col=c("red","blue"))
#dev.off()

pdf("Power.pdf", width=7, height=7)
plot(powerTulap,col="red",pch=3,xaxt="n",ylab="empirical power",xlab="epsilon",ylim=c(0.06,.44))
lines(powerTulap,col="red",lty=4)
axis(1,at=seq(1,length(epvect)),labels=c("1/64","1/32","1/16","1/8","1/4","1/2","1","2"))
points(powerNN,col="blue",pch=2)
lines(powerNN,col="blue",lty=2)
points(powerNP,col="black",pch=1)
lines(powerNP,col="black",lty=1)
legend(1,.3, c("DP UMP","DP Normal Approximation","Non-private UMP"),
       pch=c(3,2,1),lty=c(4,2,1), col=c("red","blue","black"))
dev.off()


############################################################
###   TYPE I ERROR
############################################################

epvect = 2^seq(-6,1,by=1)
size=30

ep=1

reps=100000

tavect = seq(.05,.95,by=.05)



errorTulap=errorNN=rep(0,length(tavect))
#avgpTulap=avgpNL=avgpNP=rep(length(epvect))


for(i in 1:length(tavect)){
    print(i)
    #ep=epvect[i]
    ta=tavect[i]
    alt=null=ta
    
    ###DATA GENERATION
    X=rbinom(n=reps,size=size,prob=alt)


    ###TULAP
    U = runif(n=reps,min=-1/2,max=1/2)
    G1 = rgeom(n=reps,prob=1-exp(-ep))
    G2 = rgeom(n=reps,prob=1-exp(-ep))
    Tulap = X+U+G1-G2

    ###LAPLACE
    E1 = rexp(n=reps,rate=ep)
    E2 = rexp(n=reps,rate=ep)
    Laplace = X+E1-E2

    ###TULAP   Pvalue
    values=seq(0,size)
    B = dbinom(values,size=size,prob=null)
    pvalTulap=rep(0,reps)
    for(r in 1:reps){
        F = ptulap(values-Tulap[r],exp(-ep))
        pvalTulap[r]=t(F)%*%B
    }
    #avgpTulap[i] = mean(pvalTulap)
    errorTulap[i] = mean(pvalTulap<.05)

###   Playing with pvals
   # plot(ecdf(pvalTulap))
   # values=seq(0,1,by=.1)
   # lines(values,values)


    
    
    

    
    ###NORMAL LAPLCE Pvalue
    #pvalNL = (1-pnl(Laplace,param = c(null*size,sqrt(size*null*(1-null)),ep,ep)))
    #powerNL[i] = mean(pvalNL<.05)

    
    ###NORMAL NORMAL Pvalue
    pvalNN = (1-pnorm(Laplace,m=null*size,s=sqrt(size*null*(1-null)+2/ep^2)))
    #avgpNN[i] = mean(pvalNN)
    errorNN[i] = mean(pvalNN<.05)

    #lines(ecdf(pvalNN),col="blue")
    


    
}


#plot(avgpTulap,col="red",pch=3,ylab="average p-value",xaxt="n",xlab="epsilon")
#lines(avgpTulap,col="red")
#axis(1,at=seq(1,length(epvect)),labels=epvect)
#points(avgpNN,col="blue")
#lines(avgpNN,col="blue",lty=2)
#legend(6,.4,c("UMP","Normal Approximation"),
#pch=c(3,1),lty=c(1,2),col=c("red","blue"))

pdf("./TypeIerror.pdf",width=7, height=7)
plot(errorTulap,col="red",pch=3,xaxt="n",ylab="empirical type I error",xlab="theta",ylim=c(.04,.06))
#lines(errorTulap,col="red")
abline(a=.05,b=0)
axis(1,at=seq(1,length(tavect)),labels=tavect)
points(errorNN,col="blue")
#lines(errorNN,col="blue",lty=2)
legend(12,.06, c("UMP","Normal Approximation"),
       pch=c(3,1),lty=c(4,2), col=c("red","blue"))


ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
#ma(errorTulap)
#lines(ma(errorTulap),col="red")
#lines(ma(errorNN),col="blue",lty=2)
library(zoo)
lines(c(errorTulap[1],rollmean(errorTulap,3),errorTulap[length(errorTulap)]),col="red",lty=4)
lines(c(errorNN[1],rollmean(errorNN,3),errorNN[length(errorTulap)]),col="blue",lty=2)
dev.off()


############################################################
###   SIMULATION VARY N instead of EP
############################################################

set.seed(1000)
ep=1
reps=10000



nvect = 2^seq(4,9)
#size=30

null=.9
alt=.95



powerTulap=powerNN=powerNP=rep(0,length(nvect))
avgpTulap=avgpNN=avgpNP=rep(length(nvect))


for(i in 1:length(nvect)){
    print(i)
    size=nvect[i]

    ###DATA GENERATION
    X=rbinom(n=reps,size=size,prob=alt)

    Z = X+runif(n=reps,min=-1/2,max=1/2)

    pvalNP = rep(0,reps)
    values = seq(0,size)
    B = dbinom(values, size=size,prob=null)
    for(r in 1:reps){
        F = punif(values-Z[r],min=-1/2,max=1/2)
        pvalNP[r] = t(F)%*%B
    }
    avgpNP[i] = mean(pvalNP)
    powerNP[i] = mean(pvalNP<.05)
    


    ###TULAP
    U = runif(n=reps,min=-1/2,max=1/2)
    G1 = rgeom(n=reps,prob=1-exp(-ep))
    G2 = rgeom(n=reps,prob=1-exp(-ep))
    Tulap = X+U+G1-G2

    ###LAPLACE
    E1 = rexp(n=reps,rate=ep)
    E2 = rexp(n=reps,rate=ep)
    Laplace = X+E1-E2

    ###TULAP   Pvalue
    #values=seq(0,size)
    #B = dbinom(values,size=size,prob=null)
    pvalTulap=rep(0,reps)
    for(r in 1:reps){
        F = ptulap(values-Tulap[r],exp(-ep))
        pvalTulap[r]=t(F)%*%B
    }
    avgpTulap[i] = mean(pvalTulap)
    powerTulap[i] = mean(pvalTulap<.05)

    
    ###NORMAL NORMAL Pvalue
    pvalNN = (1-pnorm(Laplace,m=null*size,s=sqrt(size*null*(1-null)+2/ep^2)))
    avgpNN[i] = mean(pvalNN)
    powerNN[i] = mean(pvalNN<.05)

    #lines(ecdf(pvalNN),col="blue")
    


    
}
powerNN/powerTulap
powerTulap-powerNN

#plot(avgpTulap/avgpNN,type="l")
#plot(powerTulap/powerNN,type="l")

#pdf("PValue.pdf",width=10,height=4)
#plot((avgpTulap),col="red",pch=3,ylab="average p-value",xaxt="n",xlab="n")
#lines((avgpTulap),col="red")
#axis(1,at=seq(1,length(nvect)),labels=nvect)
#points((avgpNN),col="blue")
#lines((avgpNN),col="blue",lty=2)
#lines((avgpNP),col="black",lty=1)
#legend(5,.4,c("UMP","Normal Approximation"),
#       pch=c(3,1),lty=c(1,2),col=c("red","blue"))
#dev.off()

pdf("Power.pdf", width=7, height=7)
plot(powerTulap,col="red",pch=3,xaxt="n",ylab="empirical power",xlab="n",ylim=c(0.05,1))
abline(a=1,b=0)
lines(powerTulap,col="red",lty=4)
axis(1,at=seq(1,length(nvect)),labels=nvect)
points(powerNN,col="blue",pch=2)
lines(powerNN,col="blue",lty=2)
points(powerNP,col="black",pch=1)
lines(powerNP,col="black",lty=1)
legend(1,.9, c("DP UMP","DP Normal Approximation","Non-private UMP"),
       pch=c(3,2,1),lty=c(4,2,1), col=c("red","blue","black"))
dev.off()


