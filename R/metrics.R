#' Typical metrics for regression
#'
#' @param sim simulated values
#' @param obs observed values
#' @param use everything or pairwise
#' @param method a character string indicating which correlation coefficient
#'
#' @return a vector of metric
#' @export
#'
#' @import Metrics
#'
#' @examples
#' NULL
metric <- function(sim, obs, use="pairwise", method="pearson") {

  n <- length(obs)
  # Check length of sim
  if (length(sim) < n) {
    stop("Sim and obs data match.")
  }

  error <- sim - obs

  me <- mean(error, na.rm = TRUE)
  mae <- mean(abs(error), na.rm = TRUE)

  mse <- mean(error^2, na.rm = TRUE)

  pe <- error / obs * 100
  mape <- mean(abs(pe), na.rm = TRUE)
  mpe <- mean(pe, na.rm = TRUE)

  rho <- cor(obs,sim, use=use, method=method)

  out <- c(me, sqrt(mse), rho, mae, mpe, mape)
  names(out) <- c("ME", "RMSE", "R", "MAE", "MPE", "MAPE")

  return(out)
}

#' Nash-Sutcliffe model efficiency
#'
#' @param sim simulated values
#' @param obs observed values
#' @param j the exponent to be used in the computation of the modified Nash-Sutcliffe effciency. The default value is j=1
#'
#' @return  a vector of metric
#' @export
#'
#' @examples
#' NULL
NSE <- function(sim, obs,j=1) {

  # Nash-Sutcliffe model efficiency (NSE)
  denominator <- sum( abs(obs - mean(obs))^2 )
  nse <- 1 - ( sum( (obs - sim)^2 ) / denominator )

  # modified Nash-Sutcliffe model efficiency (NSE)
  denominator <- sum( abs(obs - mean(obs))^j )
  mnse <- 1 - ( sum( abs(obs - sim)^j ) / denominator )

  out <- c(nse, mnse)
  names(out) <- c("NSE", "mNSE")

  return(out)
}
