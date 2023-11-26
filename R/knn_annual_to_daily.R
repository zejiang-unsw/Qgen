#' Annual to daily disaggregation
#'
#' Disaggregate Q1 by sampling from Q2 using (modified) k-nearest-neighbour resampling (Nowalk et al, 2010).
#'
#' @param Q1 An annual dataframe (year, Qa).
#' @param Q2 A daily dataframe (year, month, day, Qd).
#' @param K Number of nearest neighbours.
#' @param flag logic value if same year can be sampled for leave out cross validation.
#'
#' @return A data frame (year, nn, month, day, Qd) (nn for the sampled one of the k nearest neighbours).
#' @references Nowak, K., J. Prairie, B. Rajagopalan, and U. Lall (2010), A nonparametric stochastic approach for multisite disaggregation of annual to daily streamflow, Water Resour. Res., 46, W08529, doi:10.102/2009WR008530.
#' @export
#' @import lubridate
knn_annual_to_daily <- function(Q1, Q2, K, flag=TRUE) {
    colnames(Q1) <- c('year', 'Qa')
    # Incomplete years are not sampled
    Q2a <- Q2 %>% group_by(year) %>% summarise(Qa = mean(Qd)) %>% ungroup() %>% dplyr::filter(!is.na(Qa))
    Q2  <- Q2 %>% group_by(year) %>% mutate(index = Qd / mean(Qd)) %>% ungroup()
    if (missing(K)) K <- ceiling(sqrt(nrow(Q2a))) # K = sqrt(number of years in observed record)

    knn1 <- function(q) {
        # Return one of the K nearest neighbours randomly
        # Argument
        #   q      : an annual flow value
        #   output : the year corresponding to the nearest neighbour

        distance <- abs(Q2a$Qa - q) # Eucledian distance

        #  sample() needs a probability vector of length > 1, so if K = 1 we don't sample
        if (K > 1) {
            knn      <- order(distance)[1:K]             # Take the K nearest neighbours
            w        <- (1 / 1:K) / sum(1 / 1:K)         # Calculate weights
            row      <- sample(knn, size = 2, prob = w)  # Sample according to weights
        } else row <- order(distance)[1:2]

        return(Q2a$year[row])
    }

    map_flow <- function(v) {
        # Map an annual flow to a vector of monthly flow
        # Arguments
        #   v: a vector with 2 elements (year, Qa)
        # Returns a data frame (year, nn, month, Qm)

        y.nn     <- knn1(v['Qa'])       # nearest neighbour selected
        #cat(y.nn)
        if(flag){
          y <- y.nn[1]
        } else {
          if(v['year']!=y.nn[1]) y <- y.nn[1] else y <- y.nn[2]
        }
        Q.sub <- dplyr::filter(Q2, year == y)
        #print(Q.sub)


        # check if leap year
        if(leap_year(v['year']) && !leap_year(y)) {
            doy.leap = 31+28
            Q.sub <- Q.sub %>% add_row(day=29,month=2,year=v['year'],Qd=NA,
                                       index=Q.sub$index[sample((doy.leap-7):(doy.leap+7),1)],
                                       .after=doy.leap)
        }

        if(!leap_year(v['year']) && leap_year(y)) {
          Q.sub$index[31+28] <- Q.sub$index[31+28] + Q.sub$index[31+29]
          Q.sub <- Q.sub[-(31+29),]
        }

        n     <- nrow(Q.sub)
        #cat(n)
        q1    <- v['Qa'] * Q.sub$index # calculate monthly flow
        #print(length(q1))

        return(data.frame(year = rep(v['year'], n),
                          nn = rep(y, n),
                          month = Q.sub$month,
                          day = Q.sub$day,
                          Qd = q1
                          ))
    }

    # Apply map_flow to each row of Q1 (i.e. each year in Q1)
    bind_rows(apply(Q1, 1, map_flow))
}
