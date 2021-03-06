# Function in R that returns the final payoff of the following strategy:
# A long position of j units of a call with Strike K1 and Maturity T
# A long position of h units of a call with Strike K2 and Maturity T
# A short position of h + j units of a call with Strike K and Maturity T
# where K is a convex linear combination of K1 and K2:  K=ωK1+(1−ω)K2 and ω=j/h+j

j=1
h=1
k1=95
k2=85
s_t= seq(75, 120, by = 5)

butterfly <-  function(j,h,k1,k2,Tn,s_t){
  payoff1 <- pmax((s_t-k1)*j,0)
  payoff2 <- pmax((s_t-k2)*h,0)
  omega <- (j/h+j)
  kbar <- omega*k1+(1-omega)*k2
  payoff3 <- -pmax((s_t-kbar)*(h+j),0)
  par(mfcol=c(2,2))
  plot1 <- plot(s_t, payoff1, type = "l", col= "green", main = "strategy 1")
  plot2 <- plot(s_t, payoff2, type = "l", col= "blue", main = "strategy 2")
  plot3 <- plot(s_t, payoff3, type = "l", col= "red", main = "strategy 3")
  return(list("payoff1"=payoff1,"payoff2"=payoff2,"payoff3"=payoff3))
  
}


butterfly(j=j, h=h, k1=k1, k2=k2, s_t=s_t)

# function in R that uses the uniperiodal Binomial pricing formula to evaluate 
# an European Put Option.
# The function returns a list containing the price of the Put at time t0, 
# the Martingale probability measure Q = (qu, qd) and the units y⋆ of the risky asset 
# in the replicating strategy
s0=100
k=80
r=0.25
T1=1
t0=0
u=5
d=0.1

UniBin <- function(s0,k,t0,T1,r,u,d){
  qu <- ((1+r)*exp(T1-t0)-d)/(u-d)
  qd <- 1-qu
  I_s0u = pmax(k-s0*u,0)
  I_s0d = pmax(k-s0*d,0)
  p <- (1/(1+r)^(T1-t0))*(qu*I_s0u+qd*I_s0d)
  y <- (1/s0)*((I_s0u-I_s0d)/u - d)
  return(list("price"=p,"y*"=y,"qu"=qu,"qd"=qd))
}
UniBin(s0=s0,k=k,t0=t0,T1=T1,r=r,u=u,d=d)

# function in R that internally calls the previous function and computes the price
# at time t0 of the European Call Option using the put-call parity.

parity <- function(s0,t0,T1,r,k,u,d){
  I_s0u = pmax(k-s0*u,0)
  I_s0d = pmax(k-s0*d,0)
  value <- UniBin(s0,k,t0,T1,r,u,d)
  q_u <- value$qu
  q_d <- value$qd
  C <- s0-k*exp(-r*(T1-t0))+ (1/(1+r)^(T1-t0))*(q_u*I_s0u+q_d*I_s0d)
  return(list(value,"Price of the Call"=C))
}
parity(s0=s0,t0=t0,T1=T1,r=r,k=k,u=u,d=d)

# Function that internally calls the two previous functions that check if the Merton's
# constraints are staisfied.
merton <- function(s0,t0,T1,r,k,u,d){
  PriceOfCall <- parity(s0=s0,t0=t0,T1=T1,r=r,k=k,u=u,d=d)
  CallPrice <- PriceOfCall$`Price of the Call`
  res <- logical(length(CallPrice))
  check <- CallPrice >= max(s0-k*exp(-r*(T1-t0)),0)
  res <- check
  return(res)
}
merton(s0=s0,t0=t0,T1=T1,r=r,k=k,u=u,d=d)


# Function in R that, using the Backward Induction, returns the no-arbitrage price of
# a Put option (European and American) in the Multiperiodal Binomial Model with (N = 3)
# periods and establish which one is more convenient.
s0 = 10
Tn = 1 
u = 1.25
d = 0.8 
r = 0.04
k = 11
N = 3 
European = TRUE

PutOption <- function(s0,Tn,u,d,r,k,N,European = TRUE){
  TreeNodes <- matrix(0,N+1,N+1)
  TreeNodes[1,1] <- s0
  for(t in c(2:(N+1))){
    NumberOfSteps <- t-1
    for(j in c(1:t)){
      NumberOfUp <- t-j #To count the number of up movement 
      #using the index t (column index) and the index j (row index)
      
      TreeNodes[j,t] <- s0*u^NumberOfUp*d^(NumberOfSteps-NumberOfUp)
    }
  }
  #Backward Induction 
  #European Final Payoff 
  if(European==TRUE){
    Delta <- Tn/N
    EuFinalPayoff <- pmax(k-TreeNodes[,N+1],0)
    qu <- ((1+r)^Delta-d) / (u-d)
    qd <- 1-qu
    
    TreeNodePut <- matrix(0, N+1,N+1)
    TreeNodePut[,N+1] <- EuFinalPayoff
    for(t in c((N+1):2)){
      for(j in c(1:(t-1))){
        TreeNodePut[j,t-1] <- 1/(1+r)^(Delta)*(TreeNodePut[j,t]*qu+TreeNodePut[j+1,t]*qd)
        
        
      }
    }
    
  } 
  #Backward Induction 
  #American put option 
  else{
    Delta <- Tn/N
    USFinalpayoff <- pmax(k-TreeNodes[,N+1],0)
    qu <- ((1+r)^Delta-d) / (u-d)
    qd <- 1-qu
    
    TreeNodePut <- matrix(0, N+1,N+1)
    TreeNodePut[,N+1] <- USFinalpayoff
    for(t in c((N+1):2)){
      for(j in c(1:(t-1))){
        TreeNodePut[j,t-1] <- pmax(1/(1+r)^(Delta)*(TreeNodePut[j,t]*qu+TreeNodePut[j+1,t]*qd),pmax(k-TreeNodes[j,t-1],0))
        
      }
    }
  }
  return(TreeNodePut[j,t-1])
  
}

American <- PutOption(10,1,1.25,0.8,0.04,11,3,European = FALSE)
European <- PutOption(10,1,1.25,0.8,0.04,11,3,European = TRUE)

MoreConvenient <- function(American,Europian,k){
  AmericanValue <- max(k-American,0)
  EuropeanValue <- max(k-European,0)
  if (EuropeanValue > AmericanValue){
    print("European option is more convenient than American")
  }else{
    print("American option is more convenient than European")
  }
  return()
}
MoreConvenient(American,Europian,11)

# Function in R that finds in a Biperiodal Binomial model the no arbitrage price of
# an Asian Option with the following final payoff: max((S_t - (1/N+1)*sum(S_t_i));0)
s0 <- 1
u <- 2
d <- 0.5
r <- 0
T1 <- 1
N <- 2

delta <- T1/N

qu <- ((1+r)^delta-d)/(u-d)
qd <- 1-qu

AsianOption <- function(s0,u,d,r,qu,qd,T1,N){
  # initialization of a matrix that will contain the values of up and down movements
  matrixUD<-matrix(NA,2^N,N)
  # we use the for cycle to highlight our values of u and d
  for(i in c(1:N)){
    matrixUD[,i]<-c(rep(u,2^(N-i)),rep(d,2^(N-i)))
  }
  # we add a new column containing the probability that at time t0 the price value
  # is s0 (prob is 1)
  matrixUD <- cbind(1,matrixUD)
  # initialization of the underlying tree, final payoff and number of up movements 
  Tree <- matrix(NA,2^N,N+1)
  Payoff <- matrix(NA,2^N,1)
  NumU <- matrix(NA,2^N,1)
  for(j in c(1:2^N)){
    Tree[j,] <- s0*cumprod(matrixUD[j,])
    Payoff[j] <- max((Tree[j,N+1]) - 1/(N+1)* sum(Tree[j,]),0)
    NumU[j]<-sum(matrixUD[j,-1]==u)
  }
  # we use the analytic formula to calculate the no arbitrage price of an Asian Option
  NoArbPrice <- 1/(1+r)^T1*sum(Payoff*qu^NumU*qd^(N-NumU))
  return( NoArbPrice)
}

AsianOption(s0=s0,u=u,d=d,r=r,qu=qu,qd=qd,T1=T1,N=N)



