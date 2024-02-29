## code to prepare `DATASET` dataset goes here
# Harz mountains in Germany comprising five streamflow gauges with long daily observations.
# by Ze Jiang, UNSW
# 01 March 2024

### Update:
# Extract observed Q
# Extract observed climate: P, T, and E
# Extract modelled climate: Historical/RCP85

# Period: A: 1971-2000, B: 2021-2050 and C: 2071-2100
# Resolution: daily & monthly
rm(list=ls()); graphics.off() # clear environment and graphics
library(dplyr)
library(tidyr)

set.seed(2023-10-30) # to form a composite of GCM realizations
flag.save <- T
path.dat <- "data-raw/Harz_Data/"
stn_list <- paste0("Q",c(1,4,5))

## Obs response - q----
Harz_obs_Q <- list()
for(i_stn in stn_list){
  Qday <- read.table(paste0(path.dat,"Flow/",i_stn,"day.dat"),header = T) %>%
    arrange(iyear,imon,iday)

  head(Qday[,1:2]) #%>% print()
  tail(Qday[,1:2]) #%>% print()

  Harz_obs_Q[[length(Harz_obs_Q)+1]] <- Qday
}
names(Harz_obs_Q) <- stn_list

if(flag.save) save(Harz_obs_Q, file=paste0("data/Harz_obs_Q.Rdat"))
#usethis::use_data(Harz_obs_Q, overwrite = TRUE)

## Conditional predictor----
for(i_stn in stn_list){
  # i_stn <- "Q1"
  i_stn %>% print()

  ### Observed - climate----
  etp1 <- get(load(paste0(path.dat,"cli_obs/",i_stn, "sub/Emon.RDat")))
  pcp1 <- get(load(paste0(path.dat,"cli_obs/",i_stn, "sub/Pmon.RDat")))
  tav1 <- get(load(paste0(path.dat,"cli_obs/",i_stn, "sub/Tmon.RDat")))
  glo1 <- get(load(paste0(path.dat,"cli_obs/",i_stn, "sub/Gmon.RDat")))
  head(pcp1[,1:2]) #%>% print()
  tail(pcp1[,1:2]) #%>% print()

  Harz_obs_cli <- list(Pmon=pcp1, Tmon=tav1, Emon=etp1, Gmon=glo1)
  if(flag.save) save(Harz_obs_cli, file=paste0("data/Harz_",i_stn,"_obs_cli.Rdat"))
  #usethis::use_data(Harz_obs_cli, overwrite = TRUE)

  ## GCM hist&fut climate ----
  for(model in c("hist","rcp85")){
    # model <- "hist"
    path.hist <- paste0(path.dat,"cli_mod_",model,"/",i_stn, "sub/") #%>% print()

    Harz_gcm <- NULL
    for(i_mod in 1:6){

      #i_mod <- 1
      etp1 <- get(load(paste0(path.hist,"Emon_mod_",i_mod,".RDat")))
      pcp1 <- get(load(paste0(path.hist,"Pmon_mod_",i_mod,".RDat")))
      tav1 <- get(load(paste0(path.hist,"Tmon_mod_",i_mod,".RDat")))
      glo1 <- get(load(paste0(path.hist,"Gmon_mod_",i_mod,".RDat")))
      head(pcp1[,1:2]) #%>% print()
      tail(pcp1[,1:2]) #%>% print()

      Harz_gcm[[i_mod]] <- list(Pmon=pcp1,Tmon=tav1,Emon=etp1,Gmon=glo1)

    }


    # composite
    samp_ind <- sapply(1:nrow(pcp1), function(i) sample(1:length(Harz_gcm),size=1))

    mod_list <- list()
    var_list <- c("Pmon","Tmon","Emon","Gmon")
    for(i in 1:length(var_list)){
      tmp_mat <- sapply(Harz_gcm, function(ls) ls[[i]][[var_list[i]]])
      tmp <- sapply(1:nrow(pcp1), function(j) tmp_mat[j, samp_ind[j]])

      mod_list[[var_list[i]]] <- data.frame(year=pcp1$year, mon=pcp1$mon, tmp)
    }

    Harz_gcm[[i_mod+1]] <- mod_list


    if(flag.save) save(Harz_gcm, file=paste0("data/Harz_",i_stn,"_",model,"_mon.Rdat"))
    #usethis::use_data(Harz_gcm, overwrite = TRUE)
  }

}
