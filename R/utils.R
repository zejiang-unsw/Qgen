#' Predictor format
#'
#' @param pred.mat a matrix
#' @param lag max no of lag for predictors
#' @param mv moving window
#'
#' @return a matrix
#' @export
#'
#' @import zoo
#'
pred_mat_lag <- function(pred.mat, lag=3, mv=12){

  # 4th predictor - moving average with window size (n=3 for season/n=12 for month) of each predictor
  n <- nrow(pred.mat)
  n_pred <- ncol(pred.mat)
  Xavg = matrix(NA, nrow=n,ncol=n_pred)
  for(i in 1:ncol(pred.mat)){
    Xavg[, i] <- rollmean(pred.mat[,i],mv,fill=NA,align = 'right')
  }

  ## obtain lagged predictors
  pred.all <- NULL
  for(j in 1:n_pred){
    pred.lag = matrix(NA, nrow=n, ncol=lag)
    for(i in 1:lag){
      pred.lag[(1+i):n,i] <- pred.mat[1:(n-i),j]
    }
    pred.all <- cbind(pred.all, pred.mat[,j], pred.lag, as.numeric(Xavg[,j]))
  }
  head(pred.all)

  vars <- colnames(pred.mat)
  colnames(pred.all) <- c(vars[1], paste0(vars[1],"t-",c(1:lag,"avg")),
                            vars[2], paste0(vars[2],"t-",c(1:lag,"avg")),
                            vars[3], paste0(vars[3],"t-",c(1:lag,"avg")),
                            vars[4], paste0(vars[4],"t-",c(1:lag,"avg")))

  return(pred.all)
}


#' Convert date to season format
#'
#' @param x date
#' @param out.fmt output format
#' @param type season type
#'
#' @return a season
#' @export
#'
time2season <- function (x, out.fmt = "months", type = "default")
{
  if (!(is(x, "Date") | is(x, "POSIXct") | is(x, "POSIXt")))
    stop("Invalid argument: 'x' must be in c('Date', 'POSIXct', 'POSIXt') !")
  if (is.na(match(out.fmt, c("seasons", "months"))))
    stop("Invalid argument: 'out.fmt' must be in c('seasons', 'months')")
  valid.types <- c("default", "FrenchPolynesia","calendar")
  if (length(which(!is.na(match(type, valid.types)))) <= 0)
    stop("Invalid argument: 'type' must be in c('default', 'FrenchPolynesia')")
  months <- format(x, "%m")
  if (type == "default") {
    winter <- which(months %in% c("12", "01", "02"))
    spring <- which(months %in% c("03", "04", "05"))
    summer <- which(months %in% c("06", "07", "08"))
    autumm <- which(months %in% c("09", "10", "11"))
  }
  else if (type == "FrenchPolynesia") {
    winter <- which(months %in% c("12", "01", "02", "03"))
    spring <- which(months %in% c("04", "05"))
    summer <- which(months %in% c("06", "07", "08", "09"))
    autumm <- which(months %in% c("10", "11"))
  } else if (type == "calendar") {
    winter <- which(months %in% c("01", "02", "03"))
    spring <- which(months %in% c("04", "05", "06"))
    summer <- which(months %in% c("07", "08", "09"))
    autumm <- which(months %in% c("10", "11", "12"))
  }
  seasons <- rep(NA, length(x))
  if (out.fmt == "seasons") {
    seasons[winter] <- "winter"
    seasons[spring] <- "spring"
    seasons[summer] <- "summer"
    seasons[autumm] <- "autumm"
  }
  else {
    if (type == "default") {
      seasons[winter] <- "DJF"
      seasons[spring] <- "MAM"
      seasons[summer] <- "JJA"
      seasons[autumm] <- "SON"
    }
    else if (type == "FrenchPolynesia") {
      seasons[winter] <- "DJFM"
      seasons[spring] <- "AM"
      seasons[summer] <- "JJAS"
      seasons[autumm] <- "ON"
    }
    else if (type == "calendar") {
      seasons[winter] <- 1
      seasons[spring] <- 2
      seasons[summer] <- 3
      seasons[autumm] <- 4
    }
  }
  return(seasons)
}

#' Deseason and detrend
#'
#' @param x a vector
#' @param win.len window length
#' @param do.plot a logic value
#'
#' @return a vector
#' @export
#' @import Ecdat stats
#'
deseason_trend <- function(x, win.len=12, do.plot=FALSE){

  trend_x = stats::filter(x, rep(1 / win.len, win.len), sides = 2,
                          circular=T)

  detrend_x = x / trend_x

  m_x = t(matrix(data = detrend_x, nrow = win.len))
  seasonal_x = colMeans(m_x, na.rm = T)

  random_x = x / (trend_x * seasonal_x)

  if(do.plot){
    par(mfrow=c(2,1),mar=c(3,3,1,1), mgp=c(1.5, 0.6, 0))
    plot(as.ts(x), xlab=NA, ylab="x")
    lines(trend_x, col=2)

    plot.ts(cbind(random_x, detrend_x), col=1:2,
            xlab=NA, ylab="x'",
            plot.type = 'single')
  }

  return(random_x)

}
