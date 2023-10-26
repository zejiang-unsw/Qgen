#' Monthly to daily disaggregation
#'
#' Disaggregate Q1 by sampling from Q2 using (modified) k-nearest-neighbour resampling (Nowalk et al, 2010).
#' @param Q1 An annual dataframe (year, Qa).
#' @param Q2 A daily dataframe (year, month, day, Qd).
#' @param K Number of nearest neighbours.
#' @return A data frame (year, nn, month, day, Qd) (nn for the sampled one of the k nearest neighbours).
#' @references Nowak, K., J. Prairie, B. Rajagopalan, and U. Lall (2010), A nonparametric stochastic approach for multisite disaggregation of annual to daily streamflow, Water Resour. Res., 46, W08529, doi:10.102/2009WR008530.
#' @export
#' @import lubridate
knn_monthly_to_daily <- function(Q1, Q2, K) {
    colnames(Q1) <- c('year', 'month', 'Qm')
    # Incomplete years are not sampled
    Q2m <- Q2 %>% group_by(year, month) %>% summarise(Qm = mean(Qd)) %>% ungroup() #%>% dplyr::filter(!is.na(Qa))

    Q2  <- Q2 %>% group_by(year,month) %>% mutate(index = Qd / mean(Qd)) %>% ungroup()
    if (missing(K)) K <- ceiling(sqrt(nrow(Q2m))) # K = sqrt(number of years in observed record)

    knn1 <- function(v) {
        # Return one of the K nearest neighbours randomly
        # Argument
        #   q      : an annual flow value
        #   output : the year corresponding to the nearest neighbour

        ind_m <- which(Q2m$month %in% v['month'])
        distance <- abs(Q2m$Qm[ind_m] - v['Qm']) # Eucledian distance

        #  sample() needs a probability vector of length > 1, so if K = 1 we don't sample
        if (K > 1) {
            knn      <- order(distance)[1:K]             # Take the K nearest neighbours
            w        <- (1 / 1:K) / sum(1 / 1:K)         # Calculate weights
            row      <- sample(knn, size = 2, prob = w)  # Sample according to weights
        } else row <- order(distance)[1:2]

        #print(ind_m[row])
        return(Q2m[ind_m[row],])
    }

    map_flow <- function(v) {
        # Map an annual flow to a vector of monthly flow
        # Arguments
        #   v: a vector with 2 elements (year, Qa)
        # Returns a data frame (year, nn, month, Qm)

        # v <- Q1[1,]
        y.nn     <- knn1(v)       # nearest neighbour selected
        #print(y.nn)
        if(v['year']==y.nn$year[1] & v['month']==y.nn$month[1]) {
          y <- y.nn$year[2]
          m <- y.nn$month[2]
        }else {
          y <- y.nn$year[1]
          m <- y.nn$month[1]
        }

        Q.sub <- dplyr::filter(Q2, year == y & month==m)
        #print(Q.sub)

        # check if leap year
        if(leap_year(v['year']) & v['month']==2 & !leap_year(y)) {
            doy.leap = 28
            Q.sub.y <- dplyr::filter(Q2, year == y)
            Q.sub <- Q.sub %>% add_row(day=29,month=2,year=v['year'],Qm=Q.sub.y$Qm,
                                       index=Q.sub.y$index[sample((doy.leap-7):(doy.leap+7),1)],
                                       .after=doy.leap)
        }

        if(!leap_year(v['year']) & v['month']==2 & leap_year(y)) {
          Q.sub$index[28] <- Q.sub$index[28] + Q.sub$index[29]
          Q.sub <- Q.sub[-29,]
        }

        n     <- nrow(Q.sub)
        #cat(n)
        q1    <- v['Qm'] * Q.sub$index # calculate monthly flow
        #print(length(q1))

        return(data.frame(year = rep(v['year'], n),
                          nn = rep(y, n),
                          nnm = rep(m, n), # can do near by month -2 ~ +2
                          month = Q.sub$month,
                          day = Q.sub$day,
                          Qd = q1
                          ))
    }

    # Apply map_flow to each row of Q1 (i.e. each year in Q1)
    bind_rows(apply(Q1, 1, map_flow))
}
