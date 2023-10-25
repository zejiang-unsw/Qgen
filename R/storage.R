#' Storage
#'
#' @param Qobs a vector of monthly Q
#' @param Qsim a matrix of Q
#'
#' @return
#' @export
storage <- function(Qobs, Qsim) {
  nyr_obs <- length(Qobs)/12
  nyr_sim <- nyr_obs

  nmo_obs <- length(Qobs)
  nmo_sim <- nmo_obs

  nrea <- ncol(Qsim)

  VDmax1 <- numeric(length=nrea)  # max deficit volume (1)
  DDmax1 <- numeric(length=nrea)  # max deficit duration (2)
  VSmax1 <- numeric(length=nrea)  # max surplus volume (3)
  DSmax1 <- numeric(length=nrea)  # max surplus duration (4)
  SM1 <- numeric(length=nrea)     # storage capacity for full compensation (5)
  LDP1 <- numeric(length=nrea)    # longest dry period for full compensation (6)
  Diff1 <- array(NA, dim=c(nyr_sim*12,nrea))   # mass curve cumsum(Q-MQ); 'nrea' realis;

  Qlow <- 0.5*mean(Qobs)   # absolute threshold fixed from OBS !!!
  for (j in 1:nrea) {
    VD <- 0.0; VDmax1[j] <- 0.0
    DD <-0; DDmax1[j] <- 0
    for (i in 1:nmo_sim)
    {
      if (Qsim[i,j] < Qlow)   {VD <- VD + (Qlow-Qsim[i,j]); DD <- DD+1}
      else {VDmax1[j] <- max(VD,VDmax1[j]); VD <- 0; DDmax1[j] <- max(DD,DDmax1[j]); DD <- 0}
    }
    VDmax1[j] <- VDmax1[j]*2.628   # scale factor to obtain hm^3 from sum of m^3/s over several months (2.628*10^6 s/mon)
  }

  # (3,4) Maximum surplus volume (VSmax) and maximum surplus duration (DSmax) for Q > Qhigh
  Qhigh <- 2.0*mean(Qobs)   # absolute threshold fixed from OBS !!!
  for (j in 1:nrea) {
    VS <- 0.0; VSmax1[j] <- 0.0
    DS <-0; DSmax1[j] <- 0
    for (i in 1:nmo_sim)
    {
      if (Qsim[i,j] > Qhigh)   {VS <- VS + (Qsim[i,j]-Qhigh); DS <- DS+1}
      else {VSmax1[j] <- max(VS,VSmax1[j]); VS <- 0; DSmax1[j] <- max(DS,DSmax1[j]); DS <- 0}
    }
    VSmax1[j] <- VSmax1[j]*2.628   # scale factor to obtain hm^3 from sum of m^3/s over several months (2.628*10^6 s/mon)
  }

  # (5,6) Storage capacity for total compensation of flows (SM) and longest dry period (LDP)
  for (j in 1:nrea) {
    Diff1[,j] <- cumsum(Qsim[,j] - mean(Qsim[,j]))
    SM1[j] <- max(Diff1[,j])-min(Diff1[,j])
    SM1[j] <- SM1[j]*2.628                                      # scale factor to obtain hm^3
    LDP1[j] <- abs(which.max(Diff1[,j])-which.min(Diff1[,j]))   # get index for max and min
  }

  out <- list(VD=VDmax1,
              DD=DDmax1,
              VS=VSmax1,
              DS=DSmax1,SM=SM1,LDP=LDP1)
  return(out)

}
