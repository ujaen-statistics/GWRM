#' Extract and Visualize Model Residuals
#'
#' residuals is a method which extracts model residuals from "gw", commonly returned
#' by gw function. Optionally, it produces a normal plot with simulated envelope of the residual.
#'
#' @param object	object of class "gw" holding the fitted model
#' @param type type of residuals to be extracted. Default is pearson. response and deviance are also available.
#' @param rep	number of replications for envelope construction. Default is 19, that is the smallest 95 percent band that can be build
#' @param envelope	a logical value to specify if the envelope is required.
#' @param title	a title for the envelope.
#' @param trace	if TRUE a sort of information is printed during the running time.
#' @param parallel if TRUE use parallel executation
#' @param ncores is the number of cores that we use if parallel is TRUE
#' @param ... 	further arguments passed to or from other methods.
#'
#' @return residuals values and plot
#'
#' @examples
#' data(goals)
#' fit0 <- gw(goals ~ position, data = goals)
#' residuals(fit0, type = "pearson", rep = 19, envelope = TRUE, trace = FALSE, ncores=2)
#'
#' @importFrom stats qnorm
#' @importFrom graphics lines plot polygon
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterExport detectCores
#' @importFrom foreach foreach %dopar% %do%
#' @export

residuals.gw <- function(object, type = "pearson", rep = 19, envelope = TRUE, title = "Simulated Envelope of Residuals", trace = FALSE, parallel=TRUE, ncores=2,  ...){

  resid.gw <- function(object, type){
    mu <- object$fitted.values
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    if (type == 'pearson'){
      if (ro < 2) stop('Variance is infinite')
      else{
        varianza <- ((k + ro - 1) / (ro - 2)) * (mu + (mu ^ 2) / k)
        residuos <- sort(object$residuals / sqrt(varianza))
        return(residuos)
      }
    }
    if (type == 'deviance'){
      a <- mu * (ro - 1) / k
      gama <- a + k + ro
      y <- object$Y
      residuos <- sort(2 * (lgamma(a + y) + lgamma(k + y) - lgamma(gama + y) - (lgamma(a + mu) + lgamma(k + mu) - lgamma(gama + mu))))
      return(residuos)
    }
    if (type == 'response'){
      residuos <- sort(object$residuals)
      return(residuos)
    }
  }

  if (envelope == FALSE) return(resid.gw(object, type))

  else{
    st <- proc.time()

    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    mu <- object$fitted.values
    a <- mu * (ro - 1) / k

    n<-sum(object$W)

    if (parallel){
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      clusterExport(cl, list("object","rgw","type"),
                    envir=environment())

      residuos.sim <- foreach(j=1:(rep-1), .combine = cbind,.multicombine = TRUE,.inorder=FALSE,.packages=c("GWRM"),.verbose=as.logical(trace)) %dopar%{
        converged<-FALSE
        varResponse<-getResponse(object$formula)
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
        while(!converged){
          fit <- try(GWRM::gw(object$formula, data = datos, k = object$k, method = object$methodCode), silent = TRUE)
          if(fit$aic>0)
            converged<-TRUE
          else{ ##Generate new response values
            datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
          }
        }
        as.matrix(resid.gw(fit,type))
      }
      stopCluster(cl)
    }else{
      residuos.sim <- foreach(j=1:(rep-1), .combine = cbind,.multicombine = TRUE,.inorder=FALSE,.packages=c("GWRM"),.verbose=as.logical(trace)) %do%{
        converged<-FALSE
        varResponse<-getResponse(object$formula)
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
        while(!converged){
          fit <- try(GWRM::gw(object$formula, data = datos, k = object$k, method = object$methodCode), silent = TRUE)
          if(fit$aic>0)
            converged<-TRUE
          else{ ##Generate new response values
            datos[varResponse]<-as.matrix(rgw(n, a, k, ro))
          }
        }
        as.matrix(resid.gw(fit,type))
      }
    }
    residuos <- resid.gw(object, type)
    residuos.sim<-cbind(as.matrix(rep(residuos,object$W)),residuos.sim)


    minimos <- apply(residuos.sim[, 2:rep], 1, min)
    maximos <- apply(residuos.sim[, 2:rep], 1, max)

    t <- 1:n
    normal.score <- qnorm(t / (n + 1))
    xx <- c(normal.score, rev(normal.score))
    yy <- c(minimos, rev(maximos))
    plot(normal.score, residuos, type = "l", xlab = "Standard normal quantiles", ylab = paste("Residuals ","(",type[1],")", sep = ""), main = title)
    polygon(xx, yy, col = "gray", border = NA)
    lines(normal.score, residuos)
  }

  et <- proc.time()
  if (trace > 0)
    cat(paste("\nOverall envelope simulation process took", round((et - st)[3], 2), 'seconds'))
  else cat("\n")

  ans <-  list(type = type, residuals = residuos, sim.residuals = residuos.sim)
  class(ans) <- "residuals.gw"
  return(ans)
}
