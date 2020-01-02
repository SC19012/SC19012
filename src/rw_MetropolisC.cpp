#include <Rcpp.h>
using namespace Rcpp;
//' @title Implement a randomwalk Metropolis sampler for generating the standard Laplace distribution
//' @description Implement a randomwalk Metropolis sampler for generating the standard Laplace distribution
//' @param sigma variance
//' @param x0 initial value
//' @param N sample size
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' rw_MetropolisC(1,25,1000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rw_MetropolisC(double sigma, double x0, int N) 
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
}
