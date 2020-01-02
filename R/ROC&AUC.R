#' @title Binary classification error rate
#' @description Binary classification error rate by AUC
#' @param ypred predicte value
#' @param ytest TRUE value
#' @return ROC and AUC
#' @examples
#' \dontrun{
#' ypred <- runif(1000)
#' ytest1 <- ifelse(runif(1000) >= 0.5, 1, 0)
#' AUC(ypred, ytest2)
#' }
#' @export
AUC <- function(ypred, ytest, acc = 0.01, plot = TRUE) {
  lin <- seq(0, 1, acc)[-1]
  tpr <- c(1)
  fpr <- c(1)
  res <- 0
  for(prob in lin) {
    a1 <- tail(tpr, 1)
    b1 <- tail(fpr, 1)
    tpr = c(tpr, sum(ypred > prob & ytest == 1)/sum(ytest == 1))
    fpr = c(fpr, sum(ypred > prob & ytest == 0)/sum(ytest == 0))
    a2 <- tail(tpr, 1)
    b2 <- tail(fpr, 1)
    res <- res + (a1 + a2)*(b1 - b2)/2
  }
  plot(fpr, tpr, type = 'l')
  return(res)
}

