#' Extract Log-Likelihood
#'
#' This function is generic; method functions can be written to handle specific classes of objects. Classes which have methods for this function include: "gw"
#' 
#' @param object	a fitted object of class inheriting from "gw".
#' @param ... 	further arguments passed to or from other methods.
#'
#' @return Returns an object of class logLik. This is a number with at least one attribute, "df" (degrees of freedom), giving the number of (estimated) parameters in the model.
#' There is a simple print method for "logLik" objects.
#'
#' @examples
#' object<-gw(goals~played,data=goals)
#' logLik(object)
#'
#' @export
#' 
logLik.gw <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$nobs
  attr(val, "df") <- length(object$coefficients)
  class(val) <- "logLik"
  val
}