##########ST502 Project 1: Liam Flaherty, Zixuang Huang, and Austin Larson###########
#######Goal is estimate parameter p from binomial random samples#######
#####1. Generate Random Samples#####
samplesize=1500                                                     #number of observations of n trials#
p=seq(from=0.01, to=0.99, by=(0.99-.01)/29)                         #various probabilities#
n=c(15, 100, 200)                                                   #number of trials#
mylist=list(); set.seed(502)                                        #initialize and replicate pseudo-randomness#
  
for (i in 1:length(n)) {                                            #for each value of n...#
  for (j in 1:length(p)) {                                          #test many values of p#
    mylist[[length(p)*(i-1)+j]]=rbinom(samplesize, n[i], p[j])
  }
}

a=52                                                                #test example#
hist(mylist[[a]], xlab="Number Successes", breaks=40, 
     main=paste0("1500 Random Samples of Binom(n=", 
                 n[1+trunc(a/length(p))], ", p=", 
                 round(p[a-(length(p)*trunc(a/length(p)))],2),
                 ")"))
summary(mylist[[a]])                                                #min=56, max=86#





#####2. Different Methods Of Generating CI's#####
waldCI=function(y,n,alpha) {
  p=y/n                                                             #estimate for p#
  c(p-qnorm(1-(alpha/2))*sqrt((p*(1-p))/n),             
    p+qnorm(1-(alpha/2))*sqrt((p*(1-p))/n))
}

adjustedwaldCI=function(y,n, alpha) {
  p=(y+2)/(n+4)                                                     #add 2 successes and 2 failures#
  c(p-qnorm(1-(alpha/2))*sqrt((p*(1-p))/n),             
    p+qnorm(1-(alpha/2))*sqrt((p*(1-p))/n))
}

scoreCI=function(y,n,alpha){
  p=y/n; z=qnorm(1-(alpha/2))
  c((p+((z^2)/(2*n))-z*sqrt((p*(1-p)+((z^2)/(4*n)))/n))/(1+((z^2)/n)),
  (p+((z^2)/(2*n))+z*sqrt((p*(1-p)+((z^2)/(4*n)))/n))/(1+((z^2)/n)))
}

clopperCI=function(y,n,alpha){                                      
  if(y==0) {
    c(0,0)
  } else if (y==n) {
    c(1,1)
  } else {
    c((1+((n-y+1)/(y*qf((alpha/2), 2*y, 2*(n-y+1)))))^(-1),         #qf_{c,a,b} is (1-c) from F with df a and b#
      (1+((n-y)/((y+1)*qf(1-(alpha/2), 2*(y+1), 2*(n-y)))))^(-1))
  }
}

#rawpercentileCI#                                                   #still need to do#

#bootstraptCI#                                                      #still need to do#





#####3. Generate the CI's#####
waldCIs=list(); adjustedwaldCIs=list(); scoreCIs=list(); clopperCIs=list()

for (i in 1:length(mylist)) {                                      #pick one n,p combo#
  number=n[1+trunc((i-1)/length(p))]
  alpha=0.05
  a=list(); b=list(); c=list(); d=list()
  
  for (j in 1:samplesize) {                                        #then examine each of 1500 observations#
    observation=mylist[[i]][j]
    a[[j]]=waldCI(observation,number,alpha)                        #create CI based on observation#
    b[[j]]=adjustedwaldCI(observation,number,alpha) 
    c[[j]]=scoreCI(observation,number,alpha) 
    d[[j]]=clopperCI(observation,number,alpha) 
  }
  
  waldCIs[[i]]=a                                                   #store the 1500 CI's as one element in list#
  adjustedwaldCIs[[i]]=b
  scoreCIs[[i]]=c
  clopperCIs[[i]]=d
}

waldCIs[[50]][[1000]]                                              #example for how to pull single CI#
waldCIs[[50]][[1000]][2]                                           #or just the upper end#








#####4. Evaluate the different CIs#####
###A) Average length of interval###
nvector=c(rep(15, 30), rep(100, 30), rep(200,30))
waldavglength=vector(); adjwaldavglength=vector(); scoreavglength=vector(); clopperavglength=vector()

for (i in 1:length(waldCIs)) {
  a=vector(); b=vector(); c=vector(); d=vector()
  
  for (j in 1:samplesize) {
    a[j]=waldCIs[[i]][[j]][2]-waldCIs[[i]][[j]][1]
    b[j]=adjustedwaldCIs[[i]][[j]][2]-adjustedwaldCIs[[i]][[j]][1]
    c[j]=scoreCIs[[i]][[j]][2]-scoreCIs[[i]][[j]][1]
    d[j]=clopperCIs[[i]][[j]][2]-clopperCIs[[i]][[j]][1]
  }
  
  waldavglength[i]=mean(a)
  adjwaldavglength[i]=mean(b)
  scoreavglength[i]=mean(c)
  clopperavglength[i]=mean(d)
}

lengthdf=data.frame(nvector,round(p,2), waldavglength, adjwaldavglength, scoreavglength, clopperavglength)
colnames(lengthdf)=c("n", "p", "wald_length", "adj_wald_length", "score_length", "clopper_length")
lengthdf


###B) Proportion in interval###
waldprop=vector(); adjwaldprop=vector(); scoreprop=vector(); clopperprop=vector()

for (i in 1:length(clopperCIs)) {
  
  a=vector(); b=vector(); c=vector(); d=vector()
  truep=p[((i-1)%%length(p))+1]
  
  for (j in 1:samplesize) {
    
    if(truep >= waldCIs[[i]][[j]][1] && truep <= waldCIs[[i]][[j]][2]) {
      a[j]=1
    } else {
      a[j]=0
    }
    
    if(truep >= adjustedwaldCIs[[i]][[j]][1] && truep <= adjustedwaldCIs[[i]][[j]][2]) {
      b[j]=1
    } else {
      b[j]=0
    }
    
    if(truep >= scoreCIs[[i]][[j]][1] && truep <= scoreCIs[[i]][[j]][2]) {
      c[j]=1
    } else {
      c[j]=0
    }
    
    if(truep >= clopperCIs[[i]][[j]][1] && truep <= clopperCIs[[i]][[j]][2]) {
      d[j]=1
    } else {
      d[j]=0
    }
  }
  
  waldprop[i]=sum(a)/samplesize
  adjwaldprop[i]=sum(b)/samplesize
  scoreprop[i]=sum(c)/samplesize
  clopperprop[i]=sum(d)/samplesize
}

propdf=data.frame(nvector,round(p,2), waldprop, adjwaldprop, scoreprop, clopperprop)
colnames(propdf)=c("n", "p", "wald_prop", "adj_wald_prop", "score_prop", "clopper_prop")
propdf





#####5. Visuals#####
length15=lengthdf[lengthdf$n==15,]
length100=lengthdf[lengthdf$n==100,]
length200=lengthdf[lengthdf$n==200,]

prop15=propdf[propdf$n==15,]
prop100=propdf[propdf$n==100,]
prop200=propdf[propdf$n==200,]

mycolors=rainbow(6)

###A) Average length, n=15###
plot(length15$p, length15$wald_length, 
     ylim=c(0,0.01+max(length15[, (ncol(length15)-2):ncol(length15)])),
     main="Average Length Of CI's For Binomial(15,p)",
     ylab="Average Length Of CI", xlab="Probability Of Success",
     type="l", lty = 1, lwd = 2, col=mycolors[1])
lines(length15$p, length15$adj_wald_length, 
      type="l", lty = 1, lwd = 2, col=mycolors[2])
lines(length15$p, length15$score_length,
      type="l", lty = 1, lwd = 2, col=mycolors[3])
lines(length15$p, length15$clopper_length,
      type="l", lty = 1, lwd = 2, col=mycolors[4])



###B) Average length, n=100###
plot(length100$p, length100$wald_length, 
     ylim=c(0,0.01+max(length100[, (ncol(length100)-2):ncol(length100)])),
     main="Average Length Of CI's For Binomial(100,p)",
     ylab="Average Length Of CI", xlab="Probability Of Success",
     type="l", lty = 1, lwd = 2, col=mycolors[1])
lines(length100$p, length100$adj_wald_length, 
      type="l", lty = 1, lwd = 2, col=mycolors[2])
lines(length100$p, length100$score_length,
      type="l", lty = 1, lwd = 2, col=mycolors[3])
lines(length100$p, length100$clopper_length,
      type="l", lty = 1, lwd = 2, col=mycolors[4])



###C) Average length, n=200###
plot(length200$p, length200$wald_length, 
     ylim=c(0,0.01+max(length200[, (ncol(length200)-2):ncol(length200)])),
     main="Average Length Of CI's For Binomial(200,p)",
     ylab="Average Length Of CI", xlab="Probability Of Success",
     type="l", lty = 1, lwd = 2, col=mycolors[1])
lines(length200$p, length200$adj_wald_length, 
      type="l", lty = 1, lwd = 2, col=mycolors[2])
lines(length200$p, length200$score_length,
      type="l", lty = 1, lwd = 2, col=mycolors[3])
lines(length200$p, length200$clopper_length,
      type="l", lty = 1, lwd = 2, col=mycolors[4])


###D) Proportion, n=15###
#can use par(mfrow=c(2,2)) for multiple plots#

plot(prop15$p, prop15$wald_prop, 
     main="Wald Proportion Containing With n=15",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[1], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop15$p, prop15$adj_wald_prop, 
     main="Adj. Wald Proportion Containing With n=15",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[2], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop15$p, prop15$score_prop, 
     main="Score Proportion Containing With n=15",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[3], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop15$p, prop15$clopper_prop, 
     main="Clopper Proportion Containing With n=15",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[4], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)



###E) Proportion, n=100###
#can use par(mfrow=c(2,2)) for multiple plots#

plot(prop100$p, prop100$wald_prop, 
     main="Wald Proportion Containing With n=100",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[1], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop100$p, prop100$adj_wald_prop, 
     main="Adj. Wald Proportion Containing With n=100",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[2], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop100$p, prop100$score_prop, 
     main="Score Proportion Containing With n=100",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[3], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop100$p, prop100$clopper_prop, 
     main="Clopper Proportion Containing With n=100",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[4], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)



###F) Proportion, n=200###
#can use par(mfrow=c(2,2)) for multiple plots#

plot(prop200$p, prop200$wald_prop, 
     main="Wald Proportion Containing With n=200",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[1], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop200$p, prop200$adj_wald_prop, 
     main="Adj. Wald Proportion Containing With n=200",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[2], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop200$p, prop200$score_prop, 
     main="Score Proportion Containing With n=200",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[3], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

plot(prop200$p, prop200$clopper_prop, 
     main="Clopper Proportion Containing With n=200",
     ylab="Proportion", xlab="True Parameter p",
     type="l", lty=1, lwd=2, col=mycolors[4], ylim=c(0.25,1))
abline(h=0.95, col="black", lty=2, lwd=1)

