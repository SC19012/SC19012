## -----------------------------------------------------------------------------
prof <- c(2119,3196,4423,5472,6068,6588,6914,6710,6819,8344,8888,9319,10979,11999,12597,13097)
num <- c(40,60,80,100,110,120,125,120,120,150,160,170,200,220,230,240)
ad.cost <- c(100,130,200,220,250,260,280,300,330,350,350,340,380,400,420,430)
years <- c(2001:2016)
lm.sol <- lm(prof~num+ad.cost)
summary(lm.sol)$coef
knitr::kable(data.frame(years,prof,num,ad.cost))
plot(lm.sol)
lm.pred <- predict(lm.sol,data.frame(num =240,ad.cost = 440),interval = "prediction",level = 0.95)
lm.pred


## -----------------------------------------------------------------------------
set.seed(3.5)
n <- 1000
sigma <- c(1,2,3)
u <- runif(n)
i <- 1
while (i <= 3) {
  x <- {-2*(sigma[i]^2)*log(1-u)}^{1/2}
  hist(x,prob = TRUE,main = expression(f(x) ==     (x/{sigma^2})*exp(-{x^2}/{2*sigma^2})),col = "YELLOW")
  y <-seq(0,10,.0001,)
  lines(y,{y/{sigma[i]^2}}*exp(-{y^2}/{2*sigma[i]^2}),col = "RED")
  i <- i+1
}


## -----------------------------------------------------------------------------
set.seed(3.11)
n <- 1e3
X1 <- rnorm(n)
X2 <- rnorm(n,3,1)
i <- 1
p <- seq(0,1,.05)
while (i<= length(p)) {
  r <- sample(c(0,1),n,TRUE,c(1-p[i],p[i]))
  i <- i+1
  Z <- r*X1+(1-r)*X2
 par(mfrow = c(1,3))
 q <- p[i]
  hist(Z,col = "dark RED" )

}

## -----------------------------------------------------------------------------
set.seed(3.18)
sigma <- matrix(c(1,.8,.5,.8,2,.6,.5,.6,.3), nrow = 3, ncol = 3,byrow = TRUE)
L <- chol(sigma)
d <- 3
i <- 1
b <- numeric(length(d))
while (i <= 3) {
  b[i] <- (rchisq(1,d-i+1))^(1/2)
  i <- i+1
}
a <- rnorm(3)
A <- matrix(c(b[1],0,0,a[1],b[2],0,a[2],a[3],b[3]), nrow = 3, ncol = 3,byrow = TRUE)
X <- L %*% A %*% t(A) %*% t(L)
cov(X)

## -----------------------------------------------------------------------------
set.seed(5.1)
n <- 1e4
t <- runif(n,0,pi/3)
theta.hat <- (pi/3)*mean(sin(t))
print(c(theta.hat,1-cos(pi/3)))

## -----------------------------------------------------------------------------
set.seed(5.10)
f <- function(u)
  exp(u)/(1+u^2)
m <- 1e3
x <- runif(m)
y <- x[1:(m/2)]
theta.ant <- mean((f(y)+f(1-y))/2)
theta.ant
mean(f(x))
sd.ant <- var((f(y)+f(1-y))/2)
sd.without <- var(f(x))
print(c(sd.ant,sd.without,(sd.without-sd.ant)/sd.without))

## -----------------------------------------------------------------------------
set.seed(5.15)
M <- 1e4
k <- 5
N <- 50
T2 <- numeric(k)
est <- matrix(0,N,2)
x <- matrix(0,k,M/k)
f <- function(u, i)
  exp(-u)/(exp(-(i-1)/k)-exp(-i/k))
g <- function(u)
  exp(-u)/(1+u^2)
for(i in 1:N) {
  est[i,1] <- mean(g(runif(M)))
  for (j in 1:k) {
    u <- runif(M/k)
    x[j,] <- -log(exp(-(j-1)/k)-(exp(-(j-1)/k)-exp(-j/k))*u)
    T2[j] <- mean(g(x[j,])/f(x[j,], j))
  }
  est[i,2] <- sum(T2)
}
est
colMeans(est)
apply(est,2,sd)
#example 5.10
m <- 10000
g <- function(x) {
    exp(-x - log(1+x^2))*(x>0)*(x<1)
    }
u <- runif(m) #f3, inverse transform method 
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat <- mean(fg)
se <- sd(fg)
print(c(theta.hat,se))

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu <- 2
UCL <- numeric(1e4)
for (i in 1:1e4) {
  x <- rchisq(n,2)
  UCL[i] <- as.integer(abs((mean(x) - mu)/(sd(x)/sqrt(n))) < qt(.975, n-1))}
print(mean(UCL))


## -----------------------------------------------------------------------------
n <- 20
qua <- c(0.025, 0.05, 0.95, 0.975)
v.qua <- numeric(length(qua))
for (i in 1:length(qua)) {
  v.qua[i] <- qnorm(qua[i],0,sqrt(6*(n-2) / ((n+1)*(n+3))))
  }
v.qua
x <- sort(rnorm(1e5,0,sqrt(6*(n-2) / ((n+1)*(n+3)))))
x[qua * 1e5]
##
se.hat <- numeric(length(n))
for(i in 1:length(qua)){
  se.hat[i] <- (qua[i]*(1-qua[i])/(n*(dnorm(v.qua[i]))^2))^(1/2)
}
se.hat

## -----------------------------------------------------------------------------
##beta
set.seed(6.7)
alpha <- c(seq(2,50,1))#Beta distribution's parameter
alpha1 <- 0.05 #significance level
power <- numeric(length(alpha))
m <- 1e3
p.val <- numeric(m)
n <- 20
sk <- function(x) { 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) } 
for (i in 1:length(alpha)) {
  for (j in 1:m) {
    x <- rbeta(n,alpha[i],alpha[i])
    p.val[j] <- 2*(1 - pnorm(abs(sk(x)/sqrt(6*(n-2) / ((n+1)*(n+3))))))
  }
  power[i] <- mean(p.val <= alpha1)
  
}
plot(alpha, power, type = 'l')
points(alpha, power)

##t
v <- c(seq(2,50,1))#Beta distribution's parameter
alpha1 <- 0.05 #significance level
power <- numeric(length(v))
p.val <- numeric(m)
m <- 1e3
n <- 20
sk <- function(x) { 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) } 
for (i in 1:length(v)) {
  for (j in 1:m) {
    x <- rt(n,v[i])
    p.val[j] <- 2*(1 - pnorm(abs(sk(x)/sqrt(6*(n-2) / ((n+1)*(n+3))))))
  }
  power[i] <- mean(p.val <= alpha1)
  
}
plot(v, power, type = 'l')
points(v, power)


## -----------------------------------------------------------------------------
set.seed(6.001)
##χ2(1)
n <- 20 
alpha <- .05 
mu0 <- 1 
m <- 10000
p <- numeric(m) 
for (j in 1:m) {
  x <- rchisq(n, 1)
  ttest <- t.test(x, alternative = "greater",mu = mu0) 
  p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
print(p.hat)

##Uniform(0,2)
for (j in 1:m) {
  x <- runif(n, 0,2)
  ttest <- t.test(x, alternative = "greater",mu = mu0) 
  p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
print(p.hat)

##Exponential(1)
for (j in 1:m) {
  x <- rexp(n, 1)
  ttest <- t.test(x, alternative = "greater",mu = mu0) 
  p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
print(p.hat)

## -----------------------------------------------------------------------------
set.seed(7.6)
library(boot)
library(bootstrap)
pairs(scor)#matrix of scatter plots
cor(scor)#sample correlation matrix

b.cor <- function(x,i)
  cor(x[i,1],x[i,2])
obj.12 <- boot(data = scor,statistic = b.cor,R = 2000)
round(c(se.12 = sd(obj.12$t)),3)

b.cor <- function(x,i)
  cor(x[i,3],x[i,4])
obj.34 <- boot(data = scor,statistic = b.cor,R = 2000)
round(c(se.34 = sd(obj.34$t)),3)

b.cor <- function(x,i)
  cor(x[i,3],x[i,5])
obj.35 <- boot(data = scor,statistic = b.cor,R = 2000)
round(c(se.35 = sd(obj.35$t)),3)

b.cor <- function(x,i)
  cor(x[i,4],x[i,5])
obj.45 <- boot(data = scor,statistic = b.cor,R = 2000)
round(c(se.45 = sd(obj.45$t)),3)


## -----------------------------------------------------------------------------
set.seed(7.01)
#normal population
library(boot)
library(bootstrap)
n <- 20
sk <- function(x,i) { 
  xbar <- mean(x[i]) 
  m3 <- mean((x[i] - xbar)^3) 
  m2 <- mean((x[i] - xbar)^2) 
  return( m3 / m2^1.5 ) } 
m <- 1000
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
mu <- 0#expection of N(0,2) dstribution's skewness
for (i in 1:m) {
  x <- rnorm(n,0,2)
  boot.obj <- boot(data = x,statistic = sk,R = 2000)
  ci <- boot.ci(boot.obj,type = c("basic","norm","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>= mu),
    'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
    'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu))
left.norm <- mean(ci.norm[,1] > mu)
right.norm <- mean(ci.norm[,2] < mu)
left.basic <- mean(ci.basic[,1] > mu)
right.basic <- mean(ci.basic[,2] < mu)
left.perc <- mean(ci.perc[,1] > mu)
right.perc <-  mean(ci.perc[,2] < mu)
print(c('left.norm'=left.norm,'right.norm'=right.norm,'left.basic'=left.basic,'right.basic'=right.basic,'left.perc'=left.perc,'right.perc'=right.perc))

#χ2(5)
sk <- function(x,i) { 
  xbar <- mean(x[i]) 
  m3 <- mean((x[i] - xbar)^3) 
  m2 <- mean((x[i] - xbar)^2) 
  return( m3 / m2^1.5 ) } 
m <- 1000
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
mu <- sqrt(8/5)#expection of χ2(5) dstribution's skewness
for (i in 1:m) {
  x <- rchisq(n,5)
  boot.obj <- boot(data = x,statistic = sk,R = 2000)
  ci <- boot.ci(boot.obj,type = c("basic","norm","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>= mu),
    'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
    'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu))
left.norm <- mean(ci.norm[,1] > mu)
right.norm <- mean(ci.norm[,2] < mu)
left.basic <- mean(ci.basic[,1] > mu)
right.basic <- mean(ci.basic[,2] < mu)
left.perc <- mean(ci.perc[,1] > mu)
right.perc <-  mean(ci.perc[,2] < mu)
print(c('left.norm'=left.norm,'right.norm'=right.norm,'left.basic'=left.basic,'right.basic'=right.basic,'left.perc'=left.perc,'right.perc'=right.perc))


## -----------------------------------------------------------------------------
set.seed(7.7)
library(bootstrap)
n <- nrow(scor) 
mat <- as.matrix(scor) ##define matrix
value <- val.hat <- numeric(5)
val.jack <- matrix(0,n,5)
sita.j <- numeric(n)
for (i in 1:n) {
  value <- eigen(cor(mat[-i,]))$value
  val.jack[i,] <- sort(value,decreasing = T)
  sita.j[i] <- val.jack[i,1]/sum(val.jack[i,])
}##obtain sita.hat by jackknife 
sita.jack <- mean(sita.j)
lambda <- sort(eigen(cor(mat))$value,decreasing = T)
sita.hat <- lambda[1]/sum(lambda)
bias.jack <- (n-1)*(mean(sita.jack)-sita.hat) ##bias
se.jack <- sqrt((n-1)*mean((sita.j-sita.hat)^2)) ##se
round(c(bias.jack = bias.jack,se.jack = se.jack),5)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- R1 <- R2 <- R3 <- R4 <- T1 <- T2 <- T3 <- T4 <- yhat1 <- yhat2 <- yhat3 <- yhat4 <- numeric(n)
for (k in 1:n) { 
  y <- magnetic[-k]
  x <- chemical[-k]
  ##linear
  J1 <- lm(y ~ x) 
  yhat1[k] <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1[k]
  R1[k] <- (yhat1[k]-mean(magnetic))^2
  T1[k] <- (magnetic[k]-mean(magnetic))^2
  ##quadratic polynomial
  J2 <- lm(y ~ x + I(x^2)) 
  yhat2[k] <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2[k]
  R2[k] <- (yhat2[k]-mean(magnetic))^2
  T2[k] <- (magnetic[k]-mean(magnetic))^2
 ##log-log
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3[k] <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3[k]
  R3[k] <- (yhat3[k]-mean(magnetic))^2
  T3[k] <- (magnetic[k]-mean(magnetic))^2
  ## cubic polynomial model
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4[k] <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4[k] 
  R4[k] <- (yhat4[k]-mean(magnetic))^2
  T4[k] <- (magnetic[k]-mean(magnetic))^2
} 
error.square <- c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
R.square <- c(sum(R1)/sum(T1)*51/50, sum(R2)/sum(T2)*51/49, sum(R3)/sum(T3)*51/50, sum(R4)/sum(T4)*41/48)## adjusted R^2 = (SSE/n-p-1)/(SST/n-1) 

models <- c("linear","quadratic polynomial","log-log","cubic polynomial")
data.frame(models,error.square,R.square)





## -----------------------------------------------------------------------------
x <- c(158,171,193,199,230,243,248,248,250,267,271,316,327,329)
y <- c(141,148,169,181,203,213,229,244,257,260,271,309)
R <- 999
z <- c(x,y)
K <- 1:length(z)
n <- length(x)
reps <- numeric(R)
t0 <- t.test(x,y)$statistic
for (i in 1:R) {
  k <- sample(K, size = n,replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  reps[i] <- t.test(x1,y1)$statistic
  
}
p <- mean(abs(c(t0,reps)) >= abs(t0))
round(c(p,t.test(x,y)$p.value),3)


## -----------------------------------------------------------------------------
library(Ball)
library(boot)
dCov <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x); m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y))))) stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d); M <- mean(d)
    a <- sweep(d, 1, m); b <- sweep(a, 2, m)
    b + M
  }
  A <- Akl(x); B <- Akl(y)
  sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y
  return(nrow(z) * dCov(x, y)^2)
}
ns <- seq(20, 200, 20); M = 100;
P <- matrix(0, length(ns), 2)
for(m in 1:M) {
  p <- matrix(0, length(ns), 2)
  for(i in 1:length(ns)) {
    # data generate
    n <- ns[i]
    x1 <- rnorm(n); x2 <- rnorm(n); x <- cbind(x1, x2)
    e1 <- rnorm(n); e2 <- rnorm(n); e <- cbind(e1, e2)
    # Model2
    y <- x / 4 + e
    z <- cbind(x, y)
    # Ball test
    p[i, 1] <- bcov.test(z[,1:2],z[,3:4],num.permutations=99,seed=i*m*42)$p.value
    # permutatin: resampling without replacement
    boot.obj <- boot(data = z, statistic = ndCov2, R = 99, 
                     sim = "permutation", dims = c(2, 2))
    tb <- c(boot.obj$t0, boot.obj$t)
    p[i, 2] <- mean(tb>=tb[1])
  }
  P <- P + (p <= 0.05)
}
P <- P / M
plot(ns, P[, 1], type = 'l', col = 'red', ylim = c(0, 1), xlab = "N", ylab = "power")
lines(ns, P[, 2], col = 'blue')
title(main = "Model1: Power")
legend(20, 1, legend = c("Ball", "Dist-Corr"), lty = 1, col = c("red", "blue"))

## -----------------------------------------------------------------------------
P <- matrix(0, length(ns), 2)
for(m in 1:M) {
  p <- matrix(0, length(ns), 2)
  for(i in 1:length(ns)) {
    # data generate
    n <- ns[i]
    x1 <- rnorm(n); x2 <- rnorm(n); x <- cbind(x1, x2)
    e1 <- rnorm(n); e2 <- rnorm(n); e <- cbind(e1, e2)
    # Model2
    y <- x / 4 * e
    z <- cbind(x, y)
    # Ball test
    p[i, 1] <- bcov.test(z[,1:2],z[,3:4],num.permutations=99,seed=i*m*42)$p.value
    # permutatin: resampling without replacement
    boot.obj <- boot(data = z, statistic = ndCov2, R = 99, 
                     sim = "permutation", dims = c(2, 2))
    tb <- c(boot.obj$t0, boot.obj$t)
    p[i, 2] <- mean(tb>=tb[1])
  }
  P <- P + (p <= 0.05)
}
P <- P / M
plot(ns, P[, 1], type = 'l', col = 'red', ylim = c(0, 1), xlab = "N", ylab = "power")
lines(ns, P[, 2], col = 'blue')
title(main = "Model2: Power")
legend(20, 0.2, legend = c("Ball", "Dist-Corr"), lty = 1, col = c("red", "blue"))

## -----------------------------------------------------------------------------
set.seed(7.4)
library(distr)  # distr package can generate Laplace distribution and quantile
D <- DExp(1) # generate the target function
rw.Metropolis <- function(n, sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
            if (u[i] <= (d(D)(y) / d(D)(x[i-1])))
                x[i] <- y  else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
        }
    n <- 4  #degrees of freedom for target Student t dist.
    N <- 2000
    sigma <- c(.05, .2, 1, 4,  16) #different sigma

    x0 <- 25
    rw1 <- rw.Metropolis(n, sigma[1], x0, N)
    rw2 <- rw.Metropolis(n, sigma[2], x0, N)
    rw3 <- rw.Metropolis(n, sigma[3], x0, N)
    rw4 <- rw.Metropolis(n, sigma[4], x0, N)
    rw5 <- rw.Metropolis(n, sigma[5], x0, N)

    
    print(c(rw1$k, rw2$k, rw3$k, rw4$k,rw5$k))#number of candidate points rejected
    print((N-c(rw1$k, rw2$k, rw3$k, rw4$k,rw5$k))/N) #acceptance rates of candidate points

    par(mfrow=c(2,3))  #display 4 graphs together
    refline <- q(D)(c(.025, .975))
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x,rw5$x)
    for (j in 1:5) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
    par(mfrow=c(1,1)) #reset to default



## -----------------------------------------------------------------------------
n <- 1e3
x <- runif(n,0,100)
y <- log(exp(x))
z <- exp(log(x))
isTRUE(y == z)
isTRUE(y == x)
isTRUE(z == x)
isTRUE(all.equal(y,z))
isTRUE(all.equal(y,x))
isTRUE(all.equal(x,z))



## -----------------------------------------------------------------------------
#11.4
p.b <- numeric(24)
for (k in c(5:25,100,500,1000)) {
  #coefficient
  fgma <- function(x){
  return(exp(log(2)+lgamma(x/2)-(1/2)*(log(pi * (x-1)))-lgamma((x-1)/2)))
}
func <- function(y)
  #function
  (fgma(k)*integrate(function(x)(1 + (x^2)/(k-1))^(-k/2),lower = sqrt((y^2)*(k-1)/(k-y^2)),upper = Inf)$value)-(fgma(k+1)*integrate(function(x)(1 + (x^2)/k)^(-(k+1)/2),lower = sqrt((y^2)*k/(k+1-y^2)),upper = Inf)$value)
p.b[k-4] <- uniroot(func,lower = 0.1,upper = 1+sqrt(k)/2)$root
}
plot(p.b[-1],main = "A(k) in 11.4",col = "RED",xlab = "k",ylab = "A(k)")

#11.5
p.a <- numeric(24)
for (k in c(4:25,100)) {
  fgma <- function(x){
  return(exp(log(2)+lgamma(x/2)-(1/2)*(log(pi * (x-1)))-lgamma((x-1)/2)))
}
func <- function(y)
  (fgma(k)*integrate(function(x)(1 + (x^2)/(k-1))^(-k/2),lower = 0,upper = sqrt((y^2)*(k-1)/(k-y^2)))$value)-(fgma(k+1)*integrate(function(x)(1 + (x^2)/k)^(-(k+1)/2),lower = 0,upper = sqrt((y^2)*k/(k+1-y^2)))$value)
p.a[k-3] <- uniroot(func,lower = 0.1,upper = 1+sqrt(k)/2)$root
}
plot(p.a,main = "A(k) in 11.5",col = "RED",xlab = "k",ylab = "A(k)")
print(p.b[which(p.b>0)])
print(p.a[which(p.a>0)])

## -----------------------------------------------------------------------------
n.A <- 28
n.B <- 24
n.OO <- 41
n.AB <- 70
#Log - likelihood function
llk <- function(theta ,n.AO ,n.BO){
  p <- theta[1]
  q <- theta[2]
  return(-(2*n.A*log(p) + 2*n.B*log(q) + 2*n.OO*log(1-p-q) + n.AO*log((1-p-q)/p) + n.BO*log((1-p-q)/q) +n.AB*log(p*q)))
}
st <- 1#step
er <- 1#error
N <- 1000
p1 <- q1 <- numeric(N)
#initial value (random value and less than 1)
p1[1] <- 0.1
q1[1] <- 0.1
#EM algorithm
while (er > 1e-5 && st < N ) {
  p <- p1[st]
  q <- q1[st]
  #n.AO and n.BO subject to B distribution(parameter is below)
  theta <- c(p,q)
  n.AO <- n.A * (1 - (p/(2-p-2*q)))
  n.BO <- n.B * (1 - (q/(2-2*p-q)))
  #p and q mle
  mle <- optim(theta,llk, n.AO = n.AO, n.BO = n.BO)
  theta <- mle$par
  st <- st + 1
  p1[st] <- theta[1]
  q1[st] <- theta[2]
  #compute error
  er <- min(abs(p1[st] - p1[st - 1]),abs(q1[st] - q1[st - 1]))
}
c(p,q,1-p-q)


## -----------------------------------------------------------------------------
rm(list = ls())#clean environment
attach(mtcars)#attach data
#lapply
x <- list(disp,I(1/disp),disp+wt,I(1/disp)+wt)# def list
func <- function(x)#def func
  lm(mpg ~ x)
reg1 <- lapply(x,func)
reg1
#loops
out1 <- vector("list",length(x))
for (i in seq_along(x)) {
  out1[[i]] <- func(x[[i]])
  
}# def loops
out1


## -----------------------------------------------------------------------------
rm(list = ls())
bootstraps <- lapply(1:10, function(i) {
         rows <- sample(1:nrow(mtcars), rep = TRUE)
         mtcars[rows, ]
})#produce bootstraps
boot.list <- as.list(bootstraps)#change bootstraps to list
func <- function(x)
  lm(mpg ~ disp)
#lapply
reg2 <- lapply(boot.list,func)
reg2
#loops
out2 <- vector("list",length(boot.list))
for (i in seq_along(boot.list)) {
  out2[[i]] <- func(boot.list[[i]])
  
}
out2

## -----------------------------------------------------------------------------
rm(list = ls())
attach(mtcars)
#question 3
rsq <- function(mod)
  summary(mod)$r.squared#def rsq
x <- list(disp,I(1/disp),disp+wt,I(1/disp)+wt)# def list
func <- function(x)#def func
  lm(mpg ~ x)
reg1 <- lapply(x,func)
lapply(reg1,rsq)
#question4
bootstraps <- lapply(1:10, function(i) {
         rows <- sample(1:nrow(mtcars), rep = TRUE)
         mtcars[rows, ]
})
boot.list <- as.list(bootstraps)
func <- function(x)
  lm(mpg ~ disp)
reg2 <- lapply(boot.list,func)
lapply(reg2,rsq)


## -----------------------------------------------------------------------------
rm(list = ls())
trials1 <- replicate(
         100,
         t.test(rpois(10, 10), rpois(7, 10)),
         simplify = FALSE
       )#t.test
p.value1 <- sapply(trials1,function(x) x$p.value)#produce p.value from each t.test using anonymous function

trials2 <- replicate(
         100,
         t.test(rpois(10, 10), rpois(7, 10)),
         simplify = FALSE
       )
p.value2 <- sapply(trials2,'[[',3)
data.frame(p.value1 = p.value1,p.value2 = p.value2)#get rid of the anonymous function by using [[ directly.



## -----------------------------------------------------------------------------
rm(list = ls())
library(parallel)
#mcsapply
mcsapply <- function(x,func,n){
  unlist(mclapply(x,func,mc.cores = n))
}#def mcsapply(by unlist(mclapply))
#e.g from advanced R
boot_df <- function(x) 
  x[sample(nrow(x), rep = T), ] 
rsquared <- function(mod) summary(mod)$r.square 
boot_lm <- function(i) {
  rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
}
#by different numbers of core
system.time(sapply(1:1e4, boot_lm))
system.time(mcsapply(1:1e4, boot_lm,2))
system.time(mcsapply(1:1e4, boot_lm,4))


## -----------------------------------------------------------------------------
library(distr)  # distr package can generate Laplace distribution and quantile
D <- DExp(1) 
# generate the target function
rw.MetropolisR <- function(sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
       for (i in 2:N) {
         y <- rnorm(1, x[i-1], sigma)
         if (u[i] <= exp(-(abs(y) - abs(x[i-1]))))
           x[i] <- y else {
             x[i] <- x[i-1]
             k <- k + 1
      }
  }
        return(x)
        }

## -----------------------------------------------------------------------------
library(Rcpp)
#generate the target function using Cpp
cppFunction('NumericVector rw_MetropolisC(double sigma, double x0, int N) 
{
  //Metropolis Randomwalk using C
  NumericVector x(N);
  x[0] = x0;
  double u, y;
  int k = 0;
  for (int i = 1; i < N; i++) 
  {
    y = rnorm(1, x[i-1], sigma)[0];
    u = runif(1)[0];
    if (u <= exp(-(abs(y) - abs(x[i-1])))) 
    {
      x[i] = y; 
    }
    else 
    {
      x[i] = x[i-1];
      k++;
    }
  }
  return (x);
}')

## -----------------------------------------------------------------------------
library(microbenchmark)
N <-2000#sample number
sigma <- 4
x0 <- 25#initial value
ts <- microbenchmark(rw.R = rw.MetropolisR(sigma,x0,N), rw.C = rw_MetropolisC(sigma,x0,N))#time spent by both functions
summary(ts)[,c(1,3,5,6)]
rw.R <- rw.MetropolisR(sigma,x0,N)
rw.C <- rw_MetropolisC(sigma,x0,N)
Q <- ppoints(100)#100 quantiles
Q.T <- q(D)(Q)#TRUE Quantile
Q.RWR <- quantile(rw.R,Q)#Quantile by rw.MetropolisR
qqplot(Q.T, Q.RWR, main="Q-Q PLOT for R", xlab="Laplace Quantiles", ylab="Sample Quantiles")
abline(a = 0,b = 1,col = "RED")#line y=x

Q.RWC <- quantile(rw.C,Q)#Quantile by rw_MetropolisC
qqplot(Q.T, Q.RWC, main="Q-Q PLOT for C", xlab="Laplace Quantiles", ylab="Sample Quantiles")
abline(a = 0,b = 1,col = "RED")


