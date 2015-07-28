#' Partition of the variance in a GWRM model
#'
#' In a GWRM model, the variance may be split into three terms. The first component of this decomposition represents the variability due to randomness and it comes from the underlying Poisson model. The other two components refer to the variability that is not due to randomness but is explained by the presence of liability and proneness, respectively.
#'
#' @param object	an object class "gw" for which the partition is desired.
#' @param newdata	optionally, a data frame in which to look for variables with which to obtain the partition. If omitted, all the cases are used.
#' @param ... 	further arguments passed to or from other methods.
#'
#' @return Two data frames, with ratio of sources of variation and sources of variation in which variance is splitted.
#'
#' @examples
#'
#' data(goals)
#' fit <- gw(goals ~ position, data = goals)
#' pos <- factor(c("Defender", "Midfielder"), levels = c("Defender", "Midfielder", "Forward"))
#' lev <- data.frame(position = pos, played = c(17, 21))
#'
#' partvar(fit, newdata = lev)
#'
#' @export

partvar <- function(object, newdata = NULL, ...){
  UseMethod("partvar",object)
}


#' partvar gw object.
#'
#' partvar gw object
#'

#' @param object	an object class "gw" for which the partition is desired.
#' @param newdata	optionally, a data frame in which to look for variables with which to obtain the partition. If omitted, all the cases are used.
#' @param ... 	further arguments passed to or from other methods.
#' @importFrom stats delete.response .checkMFClasses
#' @export

partvar.gw <- function(object, newdata = NULL, ...){
  tt <- terms(object)
  if (!inherits(object, "gw"))
    warning("calling predict.gw(<fake-gw-object>) ...")

  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    nobs<-nrow(as.matrix(m))
    X <- model.matrix(Terms, m, offset <- rep(0, nobs))
    if (!is.null(off.num <- attr(tt, "offset")))
      for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }

  ncovars <- ncol(X)
  beta <- object$coefficients[1:(ncovars)]
  if(!object$kBool){
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
  }
  else{
    k <- object$k
    ro <- object$betaIIpars[2]
  }


  if (is.null(offset)){
    mus <- exp( X %*% beta)
  }
  else{
    mus <- exp(offset +X %*% beta)
  }

  a <- mus * (ro - 1) / k
  var <- mus * ((a + ro - 1) * (k + ro - 1)) / ((ro - 1) * (ro - 2))
  prand <- mus / var
  pliabi <- ((ro - 1) * (k + 1)) / ((a + ro - 1) * (k + ro - 1))
  pprone <- a / (a + ro - 1)
  rand <- mus
  liabi <- pliabi * var
  prone <- pprone * var
  partvar_rate = data.frame(Levels = X, Randomness = prand, Liability = pliabi, Proneness = pprone)
  partvar = data.frame(Levels = X, Randomness = rand, Liability = liabi, Proneness = prone)

  ans <- list(Prop.Variance.Components = partvar_rate[, c("Randomness", "Liability", "Proneness")], Variance.Components = partvar[, c("Randomness", "Liability", "Proneness")])
  return(ans)
}
