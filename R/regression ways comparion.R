#' @title cross validation to lasso 
#' @description cross validation to lasso BIC/AIC/CV/GCV
#' @param x observed value
#' @param y response value
#' @param lambdas a pre-specified sequence of lambdas
#' @return Nonzero parameter R/zero parameter zero and MSE
#' @examples
#' \dontrun{
#' library(MASS)
#' library(glmnet)
#' n = 20
#' beta <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
#' Sigma = matrix(0, 8, 8)
#' for(i in 1:8)
#'   for(j in 1:i) {
#'     Sigma[i,j] = Sigma[j,i] = 0.5^(i-j)
#'   }
#' noise.s <- 3
#' sgn.ratio <- 5.7
#' signal.var <- t(beta) %*% Sigma %*% t(t(beta))
#' noise.var <- noise.s^2
#' signal.s <- as.numeric(sgn.ratio*noise.var/signal.var)
#' d <- length(beta)
#' x <- mvrnorm(n = n, mu = rep(0, d), Sigma = sqrt(signal.s)*Sigma)
#' epsilon <- rnorm(n)
#' y <- t(beta %*% t(x) + noise.s * epsilon)
#' # LASSO AIC
#' bestlam <- lambda.GCV(x, y, 10^(seq(-2, 1, length.out = 50)), type = "AIC")
#' out <- glmnet(x, y, alpha = 1, lambda = bestlam)
#' lasso.coef <- predict(out, type = "coefficients", s = bestlam)[2:9, ]
#' R1 = as.numeric(lasso.coef!=0)
#' lasso.pred = predict(out, s = bestlam, newx = x)
#' mse1 <- mean((lasso.pred - y)^2)
#' zeros1 <- sum(lasso.coef==0)
#' # LASSO BIC
#' bestlam <- lambda.GCV(x, y, 10^(seq(-2, 1, length.out = 50)))
#' out <- glmnet(x, y, alpha = 1, lambda = bestlam)
#' lasso.coef <- predict(out, type = "coefficients", s = bestlam)[2:9, ]
#' R2 = as.numeric(lasso.coef!=0)
#' lasso.pred = predict(out, s = bestlam, newx = x)
#' mse2 <- mean((lasso.pred - y)^2)
#' zeros2 <- sum(lasso.coef==0)
#' # LASSO GCV
#' bestlam <- lambda.GCV(x, y, 10^(seq(-2, 1, length.out = 50)), type = "GCV")
#' out <- glmnet(x, y, alpha = 1, lambda = bestlam)
#' lasso.coef <- predict(out, type = "coefficients", s = bestlam)[2:9, ]
#' R3 = as.numeric(lasso.coef!=0)
#' lasso.pred = predict(out, s = bestlam, newx = x)
#' mse3 <- mean((lasso.pred - y)^2)
#' zeros3 <- sum(lasso.coef==0)
#' mse1
#' mse2
#' mse3
#' }
#' @export
lambda.GCV <- function(x, y, lambdas, type = "BIC")
  # lambdas: a pre-specified sequence of lambdas
{
  require(glmnet)
  n = nrow(x)
  d = ncol(x)
  error.gcv <- numeric(length(lambdas))
  for(j in 1 : length(lambdas))
  {
    # Refit the lasso model for each lambda and calculate predict errors
    lambda <- lambdas[j]
    fit.train <- glmnet(x, y, alpha = 1, lambda = lambda)
    pred.test <- predict(fit.train, newx = x, s = lambda)
    coef.test <- predict(fit.train, type = "coefficients", s = lambda)[2:(1+d), ]
    rss <- sum((y - pred.test)^2)
    k <- sum(coef.test != 0)
    beta <- as.vector(fit.train$beta)
    w = 1/beta
    w[w == Inf] = 0
    Winv = diag(w)
    M <- x %*% solve(t(x) %*% x + lambda * Winv, tol = 1e-20) %*% t(x)
    p <- sum(diag(M))
    if(type == "BIC") {
      error.gcv[j] <- k*log(n) + rss/4
    }
    else if(type == "AIC") {
      error.gcv[j] <- k*2 + rss/4
    }
    else {
      error.gcv[j] <- rss/(1-p/n)^2
    }
  }
  # Compute overall predict errors for each lambda
  index <- which.min(error.gcv)
  lambda.optimal <- lambdas[index]
  return(lambda.optimal)
}
