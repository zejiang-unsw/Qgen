## code to prepare `DATASET` dataset goes here
# Harz mountains in Germany comprising five streamflow gauges with long daily observations.
# by Ze Jiang, UNSW
# 07.07.2023

### Update:
# Extract observed Q
# Extract observed climate: P, T, and E
# Extract modelled climate: Historical/RCP85

# Period: A: 1971-2000, B: 2021-2050 and C: 2071-2100
# Resolution: daily & monthly
rm(list=ls()); graphics.off() # clear environment and graphics
library(dplyr)
library(tidyr)

flag.save <- T

## Observed data----
path.dat <- "data-raw/Harz-Data/"
### predictor - climate----
etp1 <- get(load(paste0(path.dat,"cli_obs/Emon.RDat")))
pcp1 <- get(load(paste0(path.dat,"cli_obs/Pmon.RDat")))
tav1 <- get(load(paste0(path.dat,"cli_obs/Tmon.RDat")))
glo1 <- get(load(paste0(path.dat,"cli_obs/Gmon.RDat")))
head(pcp1[,1:2]) %>% print()
tail(pcp1[,1:2]) %>% print()

Harz_obs_cli <- list(Pmon=pcp1, Tmon=tav1, Emon=etp1, Gmon=glo1)
#if(flag.save) save(Harz_obs_cli, file=paste0("Harz_obs_cli.Rdat"))
usethis::use_data(Harz_obs_cli, overwrite = TRUE)

### response - q----
Harz_obs_Q <- list()
stn_list <- paste0("Q",c(1,4,5))
for(i_stn in stn_list){
  Qday <- read.table(paste0(path.dat,"Flow/",i_stn,"day.dat"),header = T) %>%
    arrange(iyear,imon,iday)

  head(Qday[,1:2]) %>% print()
  tail(Qday[,1:2]) %>% print()

  Harz_obs_Q[[length(Harz_obs_Q)+1]] <- Qday
}
names(Harz_obs_Q) <- stn_list

#if(flag.save) save(Harz_obs_Q, file=paste0("data/Harz_obs_Q.Rdat"))
usethis::use_data(Harz_obs_Q, overwrite = TRUE)

## GCM hist----
# conditional predictor
path.hist <- "data-raw/Harz-Data/cli_mod_hist/"

Harz_gcm_hist <- NULL
for(i_mod in 1:6){

  #i_mod <- 1
  etp1 <- get(load(paste0(path.hist,"Emon_mod_",i_mod,".RDat")))
  pcp1 <- get(load(paste0(path.hist,"Pmon_mod_",i_mod,".RDat")))
  tav1 <- get(load(paste0(path.hist,"Tmon_mod_",i_mod,".RDat")))
  glo1 <- get(load(paste0(path.hist,"Gmon_mod_",i_mod,".RDat")))
  head(pcp1[,1:2]) %>% print()
  tail(pcp1[,1:2]) %>% print()

  Harz_gcm_hist[[i_mod]] <- list(Emon=etp1,Pmon=pcp1,Tmon=tav1,Gmon=glo1)

  if(flag.save) save(Harz_gcm_hist, file=paste0("data/Harz_GCM",i_mod,"_mon_1971-2005.Rdat"))
  #usethis::use_data(Harz_gcm_hist, overwrite = TRUE)
}




## GCM fut----
# conditional predictor
path.rcp85 <- "data-raw/Harz-Data/cli_mod_rcp85/"

Harz_gcm_rcp85 <- NULL
for(i_mod in 1:6){

  #i_mod <- 1

  etp1 <- get(load(paste0(path.rcp85,"Emon_mod_",i_mod,".RDat")))
  pcp1 <- get(load(paste0(path.rcp85,"Pmon_mod_",i_mod,".RDat")))
  tav1 <- get(load(paste0(path.rcp85,"Tmon_mod_",i_mod,".RDat")))
  glo1 <- get(load(paste0(path.rcp85,"Gmon_mod_",i_mod,".RDat")))
  head(pcp1[,1:2]) %>% print()
  tail(pcp1[,1:2]) %>% print()

  Harz_gcm_rcp85[[i_mod]] <- list(Emon=etp1, Pmon=pcp1,Tmon=tav1,Gmon=glo1)

  if(flag.save) save(Harz_gcm_rcp85, file=paste0("data/Harz_GCM",i_mod,"_mon_2006-2100.Rdat"))
  #usethis::use_data(Harz_gcm_rcp85, overwrite = TRUE)
}


