###################################################################################################
## TITLE:   Simulating the Tail Conditional Covariance
## AUTHOR:  Research project by Dave Jansz 
## PURPOSE: Based on the paper "On Conditional Variance and Tail Covariance" by E.A. Valdez, 2004.
##          and submitted as a research project for the course STA3134 Monte Carlo Methods 
##          http://probability.ca/jeff/teaching/1011/sta3431/
## NOTES:   Install the following Packages before running this code
##          install.packages(c("cubature","moments","mvtnorm","VGAM"))
##          library(cubature); library(moments);library(mvtnorm);library(VGAM);
###################################################################################################
## Program: TCC00_BiNormPic.r   
## Script purpose: Plots the 3D graph of the Distribution being Analyzed
##       and displays the formula for the Standard Bivariate Normal Distribution with Corr(X,Y)=0.5
##       The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
mu1<-0 # setting the expected value of x1
mu2<-0 # setting the expected value of x2
s11<-1 # setting the variance of x1
s12<-0.50 # setting the covariance between x1 and x2
s22<-1 # setting the variance of x2
rho<-0.50 # setting the correlation coefficient between x1 and x2
x1<-seq(-3,3,length=40) # generating the vector series x1
x2<-x1 # copying x1 to x2
#
f<-function(x1,x2){
	term1 <- 1/(2*pi*sqrt(s11*s22*(1-rho^2)))
	term2 <- -1/(2*(1-rho^2))
	term3 <- (x1-mu1)^2/s11
	term4 <- (x2-mu2)^2/s22
	term5 <- -2*rho*((x1-mu1)*(x2-mu2))/(sqrt(s11)*sqrt(s22))
	term1*exp(term2*(term3+term4-term5))
} # setting up the function of the multivariate normal density
#
z<-outer(x1,x2,f) # calculating the density values

persp(x1, x2, z, 
      main="The Bivariate Normal Distribution with Cov = 0.50", 
      # mathematical typesetting for the pdf of the distribution
      sub=expression(italic(f)~(bold(x))==frac(1,2~pi~sqrt(sigma[11]~
                     sigma[22]~(1-rho^2)))~phantom(0)~exp~bgroup("{",
	             list(-frac(1,2(1-rho^2)), 
	             bgroup("[", frac((x[1]~-~mu[1])^2, sigma[11])~-~2~rho~frac(x[1]~-~mu[1],
	             sqrt(sigma[11]))~ frac(x[2]~-~mu[2],sqrt(sigma[22]))~+~ 
	             frac((x[2]~-~mu[2])^2, sigma[22]),"]")),"}")),
      col="grey", 
      theta=30, phi=20, 
      r=50, 
      d=0.1, 
      expand=0.5, 
      ltheta=90, lphi=180, 
      shade=0.75, 
      ticktype="detailed", 
      nticks=5) # produces the 3-D plot
# text below the title
mtext(expression(list(mu[1]==0,mu[2]==0,sigma[11]==1,sigma[22]==1,sigma[12]==0.50,rho==0.50)), side=3) 
###################################################################################################
## Program: TCC01_Closed_Form_Solution.r   
## Script purpose: show the Closed Form Solution of the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# 1. The TCC95 of the Bivariate Normal Distribution with Cov(X,Y) = 0.50 - Closed Form Analytic Solution
TCC95 = 0.50*(1.645*dnorm(1.645)/(1-pnorm(1.645)) + 1)
TCC95
TCC95 = 0.50*(qnorm(0.95)*dnorm(qnorm(0.95))/(1-pnorm(qnorm(0.95))) + 1)   # slightly more accurate results using the qnorm() function
TCC95
###################################################################################################
## Program: TCC02_Direct_Integration.r   
## Script purpose: show the Closed Form Solution of the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# The mean, covariance, standard deviation and correlation 
meanx = 0; meany = 0;
sdxy = 0.5; sdx = 1; sdy = 1;
rho = sdxy/(sdx*sdy);

covBivariateN = function(z){
  x = z[1];
  y = z[2];
  C = 1/(2*pi*sdx*sdy*sqrt(1-rho^2));
  f = (x-0)*(y-0)*C*exp(-1/(1-rho^2)*((x-meanx)^2/(2*sdx^2)-rho*(x-meanx)*(y-meany)/(sdx*sdy)+(y-meany)^2/(2*sdy^2)));
}

TCC95 = 1/(1-.95)*adaptIntegrate(covBivariateN, lowerLimit=c(-10,qnorm(0.95)),upperLimit=c(10,10),maxEval = 10000)[[1]]
TCC95

# Integration using the dmvnorm function in the mvtnorm package
library(mvtnorm);  # Load the package that has functions for the Multivariate Normal distribution

covBivariateN = function(x){x[1]*x[2]*dmvnorm(x,mean=c(0,0), sigma=matrix(c(1,0.50,0.50,1),nrow=2,ncol=2))}

TCC95 = 1/(1-.95)*adaptIntegrate(covBivariateN, lowerLimit=c(-10,qnorm(0.95)),upperLimit=c(10,10),maxEval = 10000)[[1]]
TCC95
###################################################################################################
## Program: TCC02b_Reimann_Sum.r   
## Script purpose: Reimann Sum to estiamte the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
meanx = 0; meany = 0;
sdxy = 0.5; sdx = 1; sdy = 1;
rho = sdxy/(sdx*sdy);

m = 100
g = round(sqrt(m))                 		# no. of grid pts on each axis
x1 = rep(((-10*g):(10*g))/g, times=201)   	# these two lines give
x2 = rep(((-10*g):(10*g))/g, each=201)    #   coordinates of grid points
z = cbind(x1,x2)

hlist2 = c(0)
hlist = rep(0,m)
for (i in 1:length(x1)){
  x = z[i,1];
  y = z[i,2];
  C = 1/(2*pi*sdx*sdy*sqrt(1-rho^2));
  f = .1*.1*(x-0)*(y-0)*C*exp(-1/(1-rho^2)*((x-meanx)^2/(2*sdx^2)-rho*(x-meanx)*(y-meany)/(sdx*sdy)+(y-meany)^2/(2*sdy^2)));
  hlist[i] = f
  if (y > 1.645) {
    hlist2 = c(hlist2,f)
  }
}

TCC95 = 1/(1-.95)*sum(hlist2)
TCC95

###################################################################################################
## Program: TCC03_Classic_Monte_Carlo_Solution.r   
## Script purpose: show the Classic Monte Carlo Solution to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
library(mvtnorm);  # Load the package that has functions for the Multivariate Normal distribution

# The mean, covariance, standard deviation and correlation 
meanx = 0; meany = 0;
sdxy = 0.5; sdx = 1; sdy = 1;
rho = sdxy/(sdx*sdy);

# Exporting a dataframe to a comma-seperated values (csv) file
n = 100
meanList = c(meanx,meany)
sigmaList = matrix(c(sdx,sdxy,sdxy,sdy),nrow=2,ncol=2)
biNorm<-rmvnorm(n,mean=meanList, sigma=sigmaList)
write.table(biNorm, "C:/Users/Dave/Downloads/biNorm.csv", sep=",", 
  col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
###################################################################################################
## Program: TCC04_Monte_Carlo_Integration.r
## Script purpose: show Monte Carlo Integration with a Shifted Laplace Distribution
##                 designed to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# Standard Monte Carlo Integration to Calculate the TCC95 of a Biv Normal with Cov = 0.50
h = function(x,y) {x*y*1/(2*pi*sqrt(1-0.5^2))*exp(-1/(1-0.50^2)*(x^2/2-0.5*x*y+y^2/2))}
M = 10000000
xlist = runif(M,min = -10, max = 10)
ylist = runif(M,min = 1.645, max = 10)
funclist = 1/(1-0.95)*(10-1.645) *(10 - (-10)) * h(xlist,ylist)
TCC95 = mean(funclist)
TCC95
stdErr = sd(funclist) / sqrt(M)
stdErr
###################################################################################################
## Program: TCC05_Metropolis_Normal.r
## Script purpose: use the Metropolis Algorithm with proposals following a Normal Distribution
##                 in order to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################

# Metropolis Algorithm to calculate the Covariance and TCC95 of the bivariate normal

g = function(x) {

	return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )
}

h = function(x) { return( x[1]*x[2] ) }

M = 100000  # run length
B = 000  # amount of burn-in
X = rnorm(2)  # random initial value for X
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0;

for (i in 1:M) {
    Y = X + sigma * rnorm(2)  # MVN proposal value (dim=2) 
    U = runif(1)  # for accept/reject
    alpha = g(Y) / g(X)  # for accept/reject
    if (U < alpha) {
	X = Y  # accept proposal
        numaccept = numaccept + 1;
    }
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}

cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", numaccept/M, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")
###################################################################################################
## Program: TCC05b_Metropolis_Uniform.r
## Script purpose: use the Metropolis Algorithm with proposals following a Uniform Distribution
##                 in order to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################

# Metropolis Algorithm to calculate the Covariance and TCC95 of the bivariate normal

g = function(x) {return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )}

h = function(x) { return( x[1]*x[2] ) }

M = 20000  # run length
B = 5000  # amount of burn-in
X = rnorm(2)  # random initial value for X
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0;

for (i in 1:M) {
    Y = X + runif(n=2,min=-1,max=1)  # MVN proposal value (dim=2) 
    U = runif(1)  	# for accept/reject
    alpha = g(Y) / g(X)  # for accept/reject
    if (U < alpha) {
	X = Y  # accept proposal
        numaccept = numaccept + 1;
    }
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}

cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", numaccept/M, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")




par(mfrow = c(1,2), pty="s")
plot(x1list[5000:15000], x2list[5000:15000], xlim=c(-4,4), ylim=c(-4,4), type="l")
plot(x1list, x2list, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")
  
bvn<-cbind(x1list,x2list)

plot(bvn,col=1:1000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))
###################################################################################################
## Program: TCC06_MH_Normal.r   
## Script purpose: show the Metropolis-Hastings algorithm with Normal Proposals
##                 as an algorithm to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
#  a Metropolis-Hastings algorithm to calculate the Covariance and TCC95 of the bivariate normal
# Proposal distribution: MVN(X_{n-1}, sigma^2 (1+|X_{n-1}|^2) I)

g = function(x) {

	return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )
}

h = function(x) { return( x[1]*x[2] ) }

qq = function(x,y) {
    return( 1/(1+sum(x^2))^2 *
		exp( - sum((y-x)^2) / ( 2 * sigma^2 * (1+sum(x^2))^2 ) ) )
}

M = 50000  # run length
B = 5000  # amount of burn-in
X = rnorm(2) # overdispersed starting distribution (dim=2)
sigma = 0.1  # proposal scaling
x1list = x2list = hlist = rep(0,M)  # for keeping track of values

for (i in 1:M) {
    Y = X + sigma * (1 + sum(X^2)) * rnorm(2)  # proposal value is MVN, variance dependent upon X
    U = runif(1)  # for accept/reject
    A = ( g(Y) * qq(Y,X) ) / ( g(X) * qq(X,Y) )   # for accept/reject
    if (U < A){
	X = Y  # accept proposal
        numaccept = numaccept + 1;
    }
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}

cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", numaccept/M, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " which is close to the true covariance of 0.50")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]))
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")
###################################################################################################
## Program: TCC06b_MH_Uniform.r   
## Script purpose: show the Metropolis-Hastings algorithm with Uniform Proposals
##                 as an algorithm to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# A Metropolis-Hastings algorithm to calculate the Covariance and TCC95 of the bivariate normal
# Proposals from the jump distribution Yn ~ Unifrom(Xn-1 - 2, Xn-1 + 4)

g = function(x) {return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )}

h = function(x) { return( x[1]*x[2] ) }

qq = function(x,y) { 
  x1 = x[1]; x2 = x[2]
  y1 = y[1]; y2 = y[2]
  dunif(y1, x1-2, x1+4)*dunif(y2,x2-2,x2+4)
}

M = 30000; B = 10000; algorithm = "Metropolis-Hastings with a Uniform Jump Distribution" # run length, burn-in and algorithm name
X = runif(2,min=-1.96,max=1.96) # initial state
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0

for (i in 1:M) {
    Y = X + runif(2,min=-2, max = 4)  # assymetrical jump distribution Yn ~ Unifrom(Xn-1 - 2, Xn-1 + 4)
    U = runif(1)  # for accept/reject
    A = min(( g(Y) * qq(Y,X) ) / ( g(X) * qq(X,Y) ),1)   # qq(Y,X)/qq(X,Y) is the Metropolis-Hastings Adj. to allow for an assym. jump dist.
    if (U < A){
	X = Y  # accept proposal
        numaccept = numaccept + 1;
    }
    x1list[i] = X[1];  x2list[i] = X[2];  hlistCov[i] = h(X);  
}

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}

varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }

library(moments)
acceptanceRate = numaccept/M
EX = mean(x1list[(B+1):M]); sdX = sd(x1list[(B+1):M]); rangeX = max(x1list[(B+1):M]) - min(x1list[(B+1):M]) 
skewnessX = skewness(x1list[(B+1):M]); exKurtX = kurtosis(x1list[(B+1):M])-3; 
EY = mean(x2list[(B+1):M]); sdY = sd(x2list[(B+1):M]); rangeY = max(x2list[(B+1):M]) - min(x2list[(B+1):M]) 
skewnessY = skewness(x2list[(B+1):M]); exKurtY = kurtosis(x2list[(B+1):M])-3;
CovXY_sim = mean(hlistCov[(B+1):M]); CovXY_actual = 0.50
iidStdErr_Cov =  sd(hlistCov[(B+1):M]) / sqrt(length(hlistCov[(B+1):M])); trueStdErr_Cov = iidStdErr_Cov  * sqrt( varfact(hlistCov[(B+1):M]) )
TCC95_sim = mean(hlistTCC[-1]); TCC95_actual = 2.196430; 
iidStdErr_TCC =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1])); trueStdErr_TCC = iidStdErr_TCC * sqrt( varfact(hlistTCC[-1]) )

# Create a table that will display some common distributional measures and their standard error
X <- rbind(round(EX,4), round(sdX,4), round(rangeX,4), round(skewnessX,4), round(exKurtX,4))
Y <- rbind(round(EY,4), round(sdY,4), round(rangeY,4), round(skewnessY,4), round(exKurtY,4))
measures <- cbind(X, Y)
table <- data.frame(measures)
dimnames(table) <- list(c("Sample Mean","SD","Range", "Skewness","Ex. Kurtosis"),c(" X ", " Y "))
# Create a table that will display the common risk metrics: VaR, CTE and TCV
h <- rbind(round(CovXY_sim,6), round(TCC95_sim,6))
actuals <- rbind(round(CovXY_actual,6), round(TCC95_actual,6))
diff <- rbind(round(CovXY_sim-CovXY_actual,6), round(TCC95_sim-TCC95_actual,6))
iidSE_h <- rbind(round(iidStdErr_Cov,6), round(iidStdErr_TCC,6))
trueSE_h <- rbind(round(trueStdErr_Cov,6), round(trueStdErr_TCC,6))
measures2 <- cbind(h,actuals, diff, iidSE_h,trueSE_h)
table2 <- data.frame(measures2)
dimnames(table2) <- list(c("Cov(X,Y) ","TCC95(X|Y) "),c("  Est.  ","  Actual   ", "    Diff.  ",  "   iid SE   ","   true SE   "))
# Output the tables to the console window
sampleStats <-structure(table, heading = c("\nSimulated Bivariate Normal Measures","\n"),class = c("anova", "data.frame"))
covTable<-structure(table2, heading = c("\nSimulated Covariance & Tail Conditional Covariance","\n"),class = c("anova", "data.frame"))
cat("\n", algorithm,"\n", " (the Markov Chain has ", M , " successive steps, of which ", B, " were used as burn-in) ", "\n"); 
sampleStats ; covTable

par(mfrow = c(1,2), pty="s")
plot(x1list[(B+1):M], x2list[(B+1):M], xlim=c(-4,4), ylim=c(-4,4), type="l")
plot(x1list, x2list, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")
  
bvn<-cbind(x1list,x2list)

plot(bvn,col=1:1000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))
###################################################################################################
## Program: TCC07_Independence_Sampler.r   
## Script purpose: show the Independce Sampler algorithm to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# The Independence Sampler - a Special Metropolis-Hastings algorithm to calculate the Covariance and TCC95 of the bivariate normal
# Proposal distribution: MVN(0, I)

g = function(x) {return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )}

h = function(x) { return( x[1]*x[2] ) }

qq = function(x,y) {return( exp( - sum((y)^2) /  2  ) )}

M = 100000  # run length
B = 10000  # amount of burn-in
X = rnorm(2) # overdispersed starting distribution (dim=2)
sigma = 0.1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0

for (i in 1:M) {
    Y = rnorm(2)  # proposal value is MVN(0,I) 
    U = runif(1)  # for accept/reject
    A = ( g(Y) * qq(Y,X) ) / ( g(X) * qq(X,Y) )   # for accept/reject
    if (U < A){
	X = Y  # accept proposal
        numaccept = numaccept + 1;
    }
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}

cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", numaccept/M, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")
###################################################################################################
## Program: TCC08_MWG_Systematic.r   
## Script purpose: show how the Metropolis-Within-Gibbs algorithm with Systematic Scanning
##                 can be used to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# Metropolis-within-Gibbs algorithm (systematic scan) to Solve for the TCC95 

g = function(x) {return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )}

h = function(x) { return( x[1]*x[2] ) }


M = 10000  # run length (number of scans)
B = 1000  # amount of burn-in
X = rnorm(2)  # starting distribution
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,2*M);
numaccept = 0;

for (i in 1:M) {
  for (coord in 1:2) {
    Y = X
    Y[coord] = X[coord] + sigma * rnorm(1)  # propose in direction "coord"
    U = runif(1)  # for accept/reject
    alpha = g(Y) / g(X)  # for accept/reject
    if (U < alpha){
	X = Y  # accept proposal
	numaccept = numaccept + 1
    }
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
  }
}

cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", numaccept/(2*M), "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) , using numbers greater than 1.645, or the top ", length(hlistTCC)/(M-B))
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")


###################################################################################################
## Program: TCC09_MWG_Random.r   
## Script purpose: show how the Metropolis-Within-Gibbs algorithm with Random Scanning 
##                 can be used to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# Metropolis-within-Gibbs algorithm (random scan) to Solve for the TCC95 

g = function(x) {return( exp(-1/(1-0.50^2)*(x[1]^2/2-0.5*x[1]*x[2]+x[2]^2/2)) )}

h = function(x) { return( x[1]*x[2] ) }


M = 20000  # run length (number of scans)
B = 2000  # amount of burn-in
X = rnorm(2)  # starting distribution
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,2*M);
numaccept = 0;

for (i in 1:M) {
  coord = floor( runif(1,1,3) )  # uniform on {1,2}
  Y = X
  Y[coord] = X[coord] + sigma * rnorm(1)  # propose in direction "coord"
  U = runif(1)  # for accept/reject
  alpha = g(Y) / g(X)  # for accept/reject
  if (U < alpha){
    X = Y  # accept proposal
    numaccept = numaccept + 1
  }
  x1list[i] = X[1];
  x2list[i] = X[2];
  hlistCov[i] = h(X)  
}


cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", numaccept/M, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")
###################################################################################################
## Program: TCC10_Gibbs_Sampler.r   
## Script purpose: show how the Simple Gibbs Sampler algorithm
##                 can be used to estimate the 95th percentile of the 
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# A Simple Gibbs Sampler to Generate the Bivariate Normal 
# and then Use Monte Carlo to solve for the TCC95

h = function(x) { return( x[1]*x[2] ) }

M = 30000;   # Number of successive steps of the Markov chain 
B = 10000;	# Burn-in
X = rnorm(2)  # random initial value for X
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0;

for (i in 2:M) {
    x <- rnorm(1, 0.5 * y, sqrt(1 - 0.5^2))
    y <- rnorm(1, 0.5 * x, sqrt(1 - 0.5^2))
    X = c(x,y)
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}


cat("ran Metropolis algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", 1, " (a Gibbs sampler is a Metropolis-Hastings algorithm in which all proposals are accepted)", "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")





par(mfrow = c(1,2), pty="s")
plot(x1list[1:1000], x2list[1:1000], xlim=c(-4,4), ylim=c(-4,4), type="l")
plot(x1list, x2list, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")
  
bvn<-cbind(x1list,x2list)

plot(bvn,col=1:1000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))
###################################################################################################
## Program: TCC11_Cholesky_Decomposition.r   
## Script purpose: uses the Cholesky Decomposition Factor to simulate the  
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# Using the Marginal of X and then the Conditional of Y given X to generate the Variates
# and then Use Monte Carlo to solve for the TCC95

h = function(x) { return( x[1]*x[2] ) }

M = 40000;   # Number of successive steps of the Markov chain 
B = 0;	# Burn-in
X = rnorm(2)  # random initial value for X
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0;
meanx = 0; meany = 0;
sdxy = 0.5; sdx = 1; sdy = 1;
rho = sdxy/(sdx*sdy);
covMatrix = matrix(c(sdx,sdxy,sdxy,sdy),nrow=2,ncol=2)

L<-matrix(c(0,.8660254,1,0.5),nrow=2,ncol=2)


for (i in 1:M) {
    Z<-rbind(rnorm(1),rnorm(1))
    X = L%*%Z
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}


cat("ran  algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", 1, " (a Gibbs sampler is a Metropolis-Hastings algorithm in which all proposals are accepted)", "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")


# Using the Marginal of X and then the Conditional of Y given X to generate the Variates
# and then Use Monte Carlo to solve for the TCC95
# Equivalent method to using the for loop to create x1list and x2list
M = 1000000; B = 500000;

x1list <- rnorm(n = M, mean = 0, sd = 1)
x2list <- rnorm(n = M, mean = 0.50* x1list, sd = sqrt(1 - 0.50^2))
hlistCov = x1list*x2list      

cat("ran  algorithm for", M, "iterations, with burn-in", B, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
   
###################################################################################################
## Program: TCC11_Simulation_of_Conditional.r   
## Script purpose: Uses the Marginal of X and then the Conditional of Y given X to 
##                 generate the Variates and then Uses Monte Carlo to solve for
##                 the TCC95 of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################


h = function(x) { return( x[1]*x[2] ) }

M = 40000;   # Number of successive steps of the Markov chain 
B = 10000;	# Burn-in
X = rnorm(2)  # random initial value for X
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0;

for (i in 2:M) {
    x <- rnorm(1,0,1)
    y <- rnorm(1, 0.5 * x, sqrt(1 - 0.5^2))
    X = c(x,y)
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}


cat("ran  algorithm for", M, "iterations, with burn-in", B, "\n");
cat("acceptance rate =", 1, " (a Gibbs sampler is a Metropolis-Hastings algorithm in which all proposals are accepted)", "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")


# Using the Marginal of X and then the Conditional of Y given X to generate the Variates
# and then Use Monte Carlo to solve for the TCC95
# Equivalent method to using the for loop to create x1list and x2list
M = 1000000; B = 500000;

x1list <- rnorm(n = M, mean = 0, sd = 1)
x2list <- rnorm(n = M, mean = 0.50* x1list, sd = sqrt(1 - 0.50^2))
hlistCov = x1list*x2list      

cat("ran  algorithm for", M, "iterations, with burn-in", B, "\n");
cat("mean of h is about", mean(hlistCov[(B+1):M]), " (Covariance = 0.50)")

hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
   
###################################################################################################
## Program: TCC12_Linear_Transformation.r
## Script purpose: Simulate the bivariate normal distribution as a linear transformation of 
##                 independent random variables then do Monte Carlo to solve for the 95th
##                 Percentile of the Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# 
# then do Monte Carlo to solve for the TCC95

M = 100000
u1list = rnorm(M,0,sqrt(0.5))
u2list = rnorm(M,0,sqrt(0.5))
x1list = u1list*1.37 + u2list*.36
x2list = u1list*.36 + u2list*1.37
rbvnorm = cbind(x1list,x2list)
cov(rbvnorm)
hlistCov = x1list*x2list

hlistTCC = c(0)
for (i in 1:M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}
cat("The TCC95 was calculated to be: ", mean(hlistTCC[-1]), " (TCC95 = 2.19643) ")
    
se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")

se =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1]))
cat("iid standard error for the TCC estimator ", se, "\n")


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }
cat("true standard error is about", se * sqrt( varfact(hlistTCC[-1]) ), "\n")

###################################################################################################
## Program: TCC13_Cholesky_MC.r   
## Script purpose: Use the Cholesky factor then do Monte Carlo to find the
##                 Tail Conditional Covariance of a Bivariate Normal with Corr=.5 
##                 The Exact Value of this Risk Measure in this Scenario is TCC95=2.19643
## Author: Dave Jansz
###################################################################################################
# Using the Cholesky factor
# and then Use Monte Carlo to solve for the TCC95

h = function(x) { return( x[1]*x[2] ) }

algorithm = " Simulating the bivariate Normal using the Cholesky factor "
M = 20000; B=0;   
X = rnorm(2)  # random initial value for X
sigma = 1  # proposal scaling
x1list = x2list = hlistCov = rep(0,M)  # for keeping track of values
numaccept = 0;
rho = sdxy/(sdx*sdy);
L<-matrix(c(0,.8660254,1,0.5),nrow=2,ncol=2)
L%*%t(L)   # show that the Cholesky factor is correct

for (i in 1:M) {
    Z<-rbind(rnorm(1),rnorm(1))
    X = L%*%Z
    x1list[i] = X[1];
    x2list[i] = X[2];
    hlistCov[i] = h(X)  
}


hlistTCC = c(0)
for (i in (B+1):M){
  if (x2list[i] > 1.644854) 
    hlistTCC = c(hlistTCC,hlistCov[i]) ;
}


varfact <- function(xxx) { 2 * sum(acf(xxx, plot=FALSE)$acf) - 1 }

library(moments)
acceptanceRate = numaccept/M
EX = mean(x1list[(B+1):M]); sdX = sd(x1list[(B+1):M]); rangeX = max(x1list[(B+1):M]) - min(x1list[(B+1):M]) 
skewnessX = skewness(x1list[(B+1):M]); exKurtX = kurtosis(x1list[(B+1):M])-3; 
EY = mean(x2list[(B+1):M]); sdY = sd(x2list[(B+1):M]); rangeY = max(x2list[(B+1):M]) - min(x2list[(B+1):M]) 
skewnessY = skewness(x2list[(B+1):M]); exKurtY = kurtosis(x2list[(B+1):M])-3;
CovXY_sim = mean(hlistCov[(B+1):M]); CovXY_actual = 0.50
iidStdErr_Cov =  sd(hlistCov[(B+1):M]) / sqrt(length(hlistCov[(B+1):M])); trueStdErr_Cov = iidStdErr_Cov  * sqrt( varfact(hlistCov[(B+1):M]) )
TCC95_sim = mean(hlistTCC[-1]); TCC95_actual = 2.196430; 
iidStdErr_TCC =  sd(hlistTCC[-1]) / sqrt(length(hlistTCC[-1])); trueStdErr_TCC = iidStdErr_TCC * sqrt( varfact(hlistTCC[-1]) )

# Create a table that will display some common distributional measures and their standard error
X <- rbind(round(EX,4), round(sdX,4), round(rangeX,4), round(skewnessX,4), round(exKurtX,4))
Y <- rbind(round(EY,4), round(sdY,4), round(rangeY,4), round(skewnessY,4), round(exKurtY,4))
measures <- cbind(X, Y)
table <- data.frame(measures)
dimnames(table) <- list(c("Sample Mean","SD","Range", "Skewness","Ex. Kurtosis"),c(" X ", " Y "))
# Create a table that will display the common risk metrics: VaR, CTE and TCV
h <- rbind(round(CovXY_sim,6), round(TCC95_sim,6))
actuals <- rbind(round(CovXY_actual,6), round(TCC95_actual,6))
diff <- rbind(round(CovXY_sim-CovXY_actual,6), round(TCC95_sim-TCC95_actual,6))
iidSE_h <- rbind(round(iidStdErr_Cov,6), round(iidStdErr_TCC,6))
trueSE_h <- rbind(round(trueStdErr_Cov,6), round(trueStdErr_TCC,6))
measures2 <- cbind(h,actuals, diff, iidSE_h,trueSE_h)
table2 <- data.frame(measures2)
dimnames(table2) <- list(c("Cov(X,Y) ","TCC95(X|Y) "),c("  Est.  ","  Actual   ", "    Diff.  ",  "   iid SE   ","   true SE   "))
# Output the tables to the console window
sampleStats <-structure(table, heading = c("\nSimulated Bivariate Normal Measures","\n"),class = c("anova", "data.frame"))
covTable<-structure(table2, heading = c("\nSimulated Covariance & Tail Conditional Covariance","\n"),class = c("anova", "data.frame"))
cat("\n", algorithm,"\n"); 
sampleStats ; covTable

par(mfrow = c(1,2), pty="s")
plot(x1list[(B+1):M], x2list[(B+1):M], xlim=c(-4,4), ylim=c(-4,4), type="l")
plot(x1list, x2list, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")
  
bvn<-cbind(x1list,x2list)

plot(bvn,col=1:1000)
plot(bvn,type="l")
plot(ts(bvn[,1]))
plot(ts(bvn[,2]))
hist(bvn[,1],40)
hist(bvn[,2],40)
par(mfrow=c(1,1))