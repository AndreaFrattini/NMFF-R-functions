# Given a Dataset these codes solve the following problems:
# 1) Remove from the dataset of the (PUT-CALL) options those that do not satisfy the Merton's constraints
# 2) With the remaining out-of-the Money (PUT-CALL) options, draw the volatility surface and then
# 3) Calibrate he volatility in the Black and Scholes model

library(readxl)
dataset <- read_excel("Add here the path")
df <- data.frame(dataset)

Maturity <- 51
TimeToMat <- Maturity/365
# POINT 1
mertConst <- function(PriceOpt, S, K, r, TimeToMat, TypeOpt){
  dimens <-length(PriceOpt)
  resCheck <- logical(dimens)
  for(i in c(1:length(PriceOpt))){
    if(TypeOpt[i]=='call'){
      condLow <- PriceOpt[i] >= max(c(S[i] - K[i] * exp(-r[i] * (TimeToMat)), 0))
      condUpp <- PriceOpt[i] <= S[i]
      resCheck[i] <- condLow & condUpp
    }else{
      condLow <- PriceOpt[i] >= max(c(K[i] * exp(-r[i] * TimeToMat)-S[i], 0))
      condUpp <- PriceOpt[i] <= K[i]*exp(-r[i] * TimeToMat)
      resCheck[i] <- condLow & condUpp
    }
  }
  return(resCheck)
}

check <- mertConst(PriceOpt=df$Price.Option, S=df$Underlyng, K=df$Strike, r=df$RiskFreeRate_ccr, 
                   TimeToMat=TimeToMat, TypeOpt = df$Type.Option)
View(check) #All the options satisfy the Merton constraints, so I can move on

# POINT 2
OutOfMoney <- function(S, K, TypeOpt){
  dimens <-length(TypeOpt)
  resCheck <- logical(dimens)
  for(i in c(1:length(TypeOpt))){
    if(TypeOpt[i]=='call'){
      cond <- S[i] < K[i]
      resCheck[i] <- cond
    }else{
      cond <- S[i] > K[i]
      resCheck[i] <- cond
    }
  }
  return(resCheck)
}

check2 <- OutOfMoney( S=df$Underlyng, K=df$Strike, TypeOpt=df$Type.Option)
View(check2) # I have to delete the first 5 and the last 5 elements

df <- df[-c(1:5,39:43),] # eliminating the false values output by "OutOfMoney"
df1 <- df[c(1:16),] # dataframe with call
df2 <- df[c(17:33),] # dataframe with put
S0 <- 52.67
r <- 0.001312074
#call
BScall <- function(S0, K, TimeToMat, r, PriceOpt, sigma){
  dimens <-length(PriceOpt)
  error <- numeric(dimens)
  for(i in c(1:length(PriceOpt))){
      d1 <- (log(S0/K[i])+(r+0.5*sigma^2)*(TimeToMat))/(sigma*sqrt(TimeToMat))
      d2 <- d1-sigma*sqrt(TimeToMat)
      error[i] <- (S0*pnorm(d1)- K[i]*exp(-r*TimeToMat)*pnorm(d2))-PriceOpt[i]
  }
  return(error^2)
}

#put
BSput <- function(S0, K, TimeToMat, r, PriceOpt, sigma){
  dimens <-length(PriceOpt)
  error <- numeric(dimens)
  for(i in c(1:length(PriceOpt))){
    d1 <- (log(S0/K[i])+(r+0.5*sigma^2)*(TimeToMat))/(sigma*sqrt(TimeToMat))
    d2 <- d1-sigma*sqrt(TimeToMat)
    error[i] <- (K[i]*exp(-r*TimeToMat)*pnorm(-d2) - S0*pnorm(-d1))-PriceOpt[i]
  }
  return(error^2)
}


ImpliedVolCall <- numeric(length=length(df1$Strike))
ImpliedVolPut <- numeric(length=length(df2$Strike))

#Implied Volatility for Call
for(i in c(1:length(df1$Strike))){
  dummy<- optim(par=0.1, fn=BScall, lower=0, method ="L-BFGS-B",
                S0=S0 ,K=df1$Strike[i] , TimeToMat=TimeToMat, 
                r=r , PriceOpt=df1$Price.Option[i] )
  if(dummy$convergence==0){
    cat("\n", "successful completion")
  }else{
    cat("\n Check the prob numb. ", c(i, dummy$value))
  }
  ImpliedVolCall[i]<-dummy$par
}

#Implied Volatility for Put
for(i in c(1:length(df2$Strike))){
  dummy<- optim(par=0.1, fn=BSput, lower=0, method ="L-BFGS-B",
                S0=S0 ,K=df2$Strike[i] , TimeToMat=TimeToMat, 
                r=r , PriceOpt=df2$Price.Option[i] )
  if(dummy$convergence==0){
    cat("\n", "successful completion")
  }else{
    cat("\n Check the prob numb. ", c(i, dummy$value))
  }
  ImpliedVolPut[i]<-dummy$par
}

StrikeCall <- df1$Strike
StrikePut <- df2$Strike
plot(x=StrikeCall, y=ImpliedVolCall, type = "l")
plot(x=StrikePut, y=ImpliedVolPut, type = "l")

# POINT 3
#call
BScallmean <- function(S0, K, TimeToMat, r, PriceOpt, sigma){
  dimens <-length(PriceOpt)
  error <- numeric(dimens)
  for(i in c(1:length(PriceOpt))){
    d1 <- (log(S0/K[i])+(r+0.5*sigma^2)*(TimeToMat))/(sigma*sqrt(TimeToMat))
    d2 <- d1-sigma*sqrt(TimeToMat)
    error[i] <- (S0*pnorm(d1)- K[i]*exp(-r*TimeToMat)*pnorm(d2))-PriceOpt[i]
  }
  return(mean(error^2))
}

#put
BSputmean <- function(S0, K, TimeToMat, r, PriceOpt, sigma){
  dimens <-length(PriceOpt)
  error <- numeric(dimens)
  for(i in c(1:length(PriceOpt))){
    d1 <- (log(S0/K[i])+(r+0.5*sigma^2)*(TimeToMat))/(sigma*sqrt(TimeToMat))
    d2 <- d1-sigma*sqrt(TimeToMat)
    error[i] <- (K[i]*exp(-r*TimeToMat)*pnorm(-d2) - S0*pnorm(-d1))-PriceOpt[i]
  }
  return(mean(error^2))
}

calibrationCall=optim(par=0.1, fn=BScallmean, lower=0, method ="L-BFGS-B",
                     S0=S0 ,K=df1$Strike , TimeToMat=TimeToMat, 
                     r=r , PriceOpt=df1$Price.Option )$par
calibrationPut=optim(par=0.1, fn=BSputmean, lower=0, method ="L-BFGS-B",
                     S0=S0 ,K=df2$Strike , TimeToMat=TimeToMat, 
                     r=r , PriceOpt=df2$Price.Option )$par
calibrationCall
calibrationPut

# Starting from a Vasicek model, these codes solves the following problems:
# 1) Create a function that returns the exact conditional probability given the information at time t0 of the event{XT<a ∪ XT>b}with a<b
# 2) Write a function that estimates the event of point 1 using the Monte Carlo Approach based on the exact simulation scheme reporting also the Upper and Lower Bounds for the MC price at the level 95%.
# 3) Write a function that repeats the exercise in point 2 using the exact simulation scheme
# 4) Test the aforementioned points with numerical examples
#Inputs
X0=0
alpha=0.2
mu=0.05
sigma=0.2
N=100 
M=1000
t0=0
FinalT=1
# POINT 1
# Exact Conditional Probability
a=0.01
b=0.07

ExConPr <- function(X0,mu, alpha, sigma, t0, FinalT, a, b){
  standDevOfXT<-sigma/sqrt(2*alpha)*sqrt(1-exp(-2*alpha*FinalT-t0))
  d_bar1 <- (a-X0*exp(-alpha*(FinalT-t0))-mu*(1-exp(-alpha*(FinalT-t0))))/standDevOfXT
  d_bar2 <- (b-X0*exp(-alpha*(FinalT-t0))-mu*(1-exp(-alpha*(FinalT-t0))))/standDevOfXT
  ExactProb <- 1+pnorm(d_bar1)-pnorm(d_bar2)
  return(ExactProb)
}
ExactConditionalProb <- ExConPr(X0=X0,mu=mu, alpha=alpha, sigma=sigma, t0=t0, FinalT=FinalT,
                                a=a, b=b)
ExactConditionalProb

# POINT 2
set.seed(1)
SamplePath <- matrix(0,M,N+1)
Grid <- seq(t0,FinalT, length.out = N+1)
Delta<-(FinalT-t0)/N
SamplePath[,0] <- X0 #Inizialization 
for (i in c(2:(N+1))){
  SamplePath[,i] <- SamplePath[,i-1]*exp(-alpha*Delta)+mu*(1-exp(-alpha*Delta))+sigma/sqrt(2*alpha)*sqrt(1-exp(-2*alpha*Delta))*rnorm(M)
  colnames(SamplePath)<- paste0("t = ",Grid)  
}

X<- SamplePath[,N]
MCExactProb <- mean(X>a & X<b)
MCExactProb
UB <- mean(X) + 1.96*sd(X)
LB <- mean(X) - 1.96*sd(X)
Intervals <- c(LB,UB)
Intervals

# POINT 3
MCEuler <- function(X0,mu, alpha, sigma, t0, FinalT, a, b){
  set.seed(1)
  EulerPath <- matrix(0, nrow = M, ncol = N+1)
  EulerPath[,0]<-X0
  for(i in c(2: (N+1))){
    EulerPath[,i] < EulerPath[,i-1]+alpha*(mu-EulerPath[,i-1])*Delta+sigma*sqrt(Delta)*rnorm(M)
    EulerXT<-EulerPath[,N+1] 
    ProbB <- mean(EulerXT>b)
    ProbA <- mean(EulerXT<a)
    MCProbEuler <- ProbA+ProbB
    return(MCProbEuler)
  }
}
MCEulerProb <- MCEuler(X0,mu, alpha, sigma, t0, FinalT, a, b)
MCEulerProb

# These codes, after having downloaded the closing prices of the VIX index from January 2nd 2019 to December 30th 2020,
# solve the following problems under the assumption that the log-prices Xt follow a Vasicek model:
# 1) Estimate the Vasicek parameters using the MLE approach. Construct the Log-Likelihood function using the transition density of the log- returns.
# 2) Repeat the point 1 where the likelihood function is constructed using the transition density of the log-prices.
# 3) Estimate the Vasicek parameters using the QMLE approach and construct the Quasi Log-Likelihood function using the approximated transition density of the log-returns.
# POINT 1
library(quantmod)
getSymbols("^VIX",  from="2019-01-02", to="2020-12-31")
Price <- na.omit(VIX[,4]) 
daysInOneYear<-365
timeOfObervation <-as.numeric(index(Price))
deltati<-diff(timeOfObervation)/daysInOneYear
length(deltati)
logprice <- log(as.numeric(Price)) 
logret <- diff(logprice)
logprice_old <- log(as.numeric(Price))[-length(as.numeric(Price))]
length(logret)
length(logprice_old)


minuslogLik<- function(par,logprice_old,logret,deltati){
  alpha<-par[1]
  mu<-par[2]
  sig<- par[3]
  vecMean <-  logprice_old*exp(-alpha*deltati)+mu*(1-exp(-alpha*deltati))-logprice_old
  vecsd <- sig/sqrt(2*alpha)*sqrt(1-exp(-2*alpha*deltati))
  -sum(dnorm(logret,mean=vecMean,sd=vecsd,log=TRUE))
}

minuslogLik(par=c(0.1,0.1,0.1),logprice_old=logprice_old,logret=logret,deltati=deltati)

res <- optim(par=c(0.1,0.1,0.1),fn = minuslogLik,method ="L-BFGS-B", 
             lower=c(0.00001, -Inf, 0.0000001), logprice_old=logprice_old,
             logret=logret,deltati=deltati)
res$par
-1*res$value
res$convergence

# POINT 2
deltati<-diff(timeOfObervation)/daysInOneYear
logret <- diff(log(as.numeric(Price)))
logprice_old <- log(as.numeric(Price))[-length(as.numeric(Price))]

minuslogLik1<- function(par,logprice,logprice_old,deltati){
  alpha<-par[1]
  mu<-par[2]
  sig<- par[3]
  vecMean <-  logprice_old*exp(-alpha*deltati)+mu*(1-exp(-alpha*deltati))
  vecsd <- sqrt((sig^2)/(2*alpha)*(1-exp(-2*alpha*deltati)))
  -sum(dnorm(logprice,mean=vecMean,sd=vecsd,log=TRUE))
}

minuslogLik1(par=c(1,0.5,0.2),logprice=logprice,logprice_old=logprice_old,deltati=deltati)

res <- optim(par=c(1,0.5,0.2),fn = minuslogLik1,method ="L-BFGS-B",
             logprice=logprice, logprice_old=logprice_old,deltati=deltati)
res$par
-1*res$value
res$convergence

# POINT 3

minusQloglik_Vasicek <- function(par,logret,logprice_old,deltati){
  alpha <- par[1]
  mu <- par[2]
  sig <- par[3]
  logdens <- dnorm(logret, mean = alpha*(mu-logprice_old)*deltati, sd= sig*sqrt(deltati),log=T )
  -sum(logdens)
}

minusQloglik_Vasicek(par=c(0.7,0.5,0.2),logret=logret,
                     logprice_old=logprice_old,deltati=deltati)

res <- optim(par=c(0.7,0.5,0.2), fn =minusQloglik_Vasicek,
             logret=logret, logprice_old=logprice_old, 
             deltati=deltati, lower=c(-Inf,-Inf,0), method ="L-BFGS-B")

res$par
res$value
res$convergence


# These codes, after having downloaded the closed price of the S&P500 index from January 2nd 2019 to December 30th 2020 and solve the following problems:
# 1) Estimate the volatility (on yearly basis) using the MLE approach
# 2) Using the estimated volatility in point 1 compute the (exact) Black and Scholes price of an At-the-Money European Put option for the following inputs: t0 (today) is December 30th 2020,
# Maturity of the option is February 27th 2020, the interest rate in c.c.r is 1.5% on yearly basis.

# POINT 1
# Estimate the volatility (on yearly basis) using the MLE approach.
library(quantmod)
getSymbols("^GSPC",  from="2019-01-02", to="2020-12-31")
Price <- na.omit(GSPC[,4]) 
daysInOneYear<- 365 
timeOfObervation <-as.numeric(index(Price))
deltati<-diff(timeOfObervation)/daysInOneYear
Xt_i <- diff(log(as.numeric(Price)))

minuslogLik<- function(par,Xt_i,deltati){
  mu<-par[1]
  sig<-par[2]
  vecMean <-  (mu-0.5*sig^2)*deltati
  vecsd <-sig*sqrt(deltati)
  -sum(dnorm(Xt_i,mean=vecMean,sd=vecsd,log=TRUE))
}

minuslogLik(par=c(0.5,0.2),Xt_i=Xt_i,deltati=deltati)

res <- optim(par=c(0.5,0.2),fn = minuslogLik, lower=c(-Inf,0),method ="L-BFGS-B", 
             Xt_i=Xt_i,deltati=deltati)
res$par
-1*res$value
res$convergence
volatility <- res$par[2]
volatility

# POINT 2

#Inputs
#t0 (today) is December 30th 2019
#Maturity of the option is February 27th 2020
sigma=volatility
S0=as.numeric(Price[504,])
k=S0
r=0.015
TimeToMat=60 
daysInOneYear <- 365
deltati = TimeToMat/daysInOneYear 
d1 <- (log(S0/k)+(r+0.5*sigma^2)*(deltati))/(sigma*sqrt(deltati))
d2 <- d1-sigma*sqrt(deltati)

BSput=k*exp(-r*(deltati))*pnorm(-d2)- S0*pnorm(-d1)
BSput


