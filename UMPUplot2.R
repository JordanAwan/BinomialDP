
source("./Tulap.R")
library(tictoc)



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

    mk = c(m,k)
    return(mk)
}

n=10
theta = .75
alpha = .1
epsilon=1
delta=0
thick=3

mk=UMPU(theta=theta,n=n,alpha=alpha,epsilon=epsilon,delta=delta)
m = mk[1]
k = mk[2]
values = seq(0,n,by=.001)
b=exp(-epsilon)
q=2*delta*b/(1-b+2*delta*b)
F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
greaterK = (values>=k)
phi = F1*greaterK + F2*(1-greaterK)

pdf("UMPUplot.pdf",width=9,height=4)
plot(values,phi,type="l",lty=2,col="red",lwd=thick,xlab="x",ylab="test",cex.axis=1.5,cex.lab=1.5)
abline(v=k,lty=2,col="red",lwd=thick)


n=10
theta = .75
alpha = .01
epsilon=1
delta=0

mk=UMPU(theta=theta,n=n,alpha=alpha,epsilon=epsilon,delta=delta)
m = mk[1]
k = mk[2]
values = seq(0,n,by=.001)
b=exp(-epsilon)
q=2*delta*b/(1-b+2*delta*b)
F1 = ptulap(t=values - k-m,median=0,lambda=b,cut=q)
F2 = ptulap(t=k-values-m,median=0,lambda=b,cut=q)
greaterK = (values>=k)
phi = F1*greaterK + F2*(1-greaterK)

lines(values,phi,lty=4,col="blue",lwd=thick)
abline(v=k,lty=4,col="blue",lwd=thick)
abline(v=n*theta,lwd=thick)
legend("topright",lty=c(2,4),col=c("red","blue"),c(expression(paste(alpha,"=.1")),expression(paste(alpha,"=.01"))),lwd=thick,cex=1.5)
dev.off()
