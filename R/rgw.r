#' Genaralizing Waring
#'
#' Random generation for the Gamma distribution with parameters a, k, ro
#'
#' @param n number of generated values
#' @param a a parameter
#' @param k k parameter
#' @param ro ro parameter
#'
#' @return vector of generated values
#'
#' @importFrom stats rbeta rgamma rpois
#'
#' @examples
#' rgw(10,2,2,2)
#'
#' @export

rgw <- function(n, a, k, ro){
  if (sum(a <= 0) | sum(k <= 0) | sum(ro <= 0)) stop('Parameters must be positive')
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if ((any(!is.wholenumber(n)==TRUE) | any(n <= 0))) stop('n must be a positive integer')
  u <- rbeta(n, ro, k)
  v <- u / (1 - u)
  lam <- rgamma(n, a, v)
  salida <- rpois(n, lam)
  return(salida)
}
