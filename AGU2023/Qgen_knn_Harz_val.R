#==============================================================================#
# By Ze Jiang (ze.jiang@unsw.edu.au)
# Validation using observed climate
# Period A include A period
# Period D exclude D period

rm(list=ls()); graphics.off() # clear environment and graphics
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path));

# load package----
library(scales)
library(zoo)
library(NPRED)
library(WASP)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(data.table)

#devtools::install_github("hzambran/hydroTSM")
#devtools::install_github("itsoukal/anySim")
# library(hydroTSM)
# library(Ecdat)
# library(ExtDist)
# library(anySim)
library(Qgen)

#==============================================================================#
# global variables----
set.seed(20230607)
#set.seed(20231116)
mod_id <- 7 # GCM realization: 1-6 is raw, and 7 is the sampled scenario
station_id <- "Q5" # station
nensemble <- 100 # must larger than 1 for knn bootstrap

flag.save <- T # save the Qsim or not

flag.ver <- switch(2, "","_v1","_v2","_v3")
flag.sel <- switch(2, "ALL","NPRED","WASP") # predictor selection method

if(TRUE){
  flag.val <- switch (2,"_A","_D")
#for(flag.val in c("_A","_D")){

  # output----
  path.out <- paste0(station_id,"/cli_val",flag.val,"/")#,"/mod_",mod_id
  dir.create(path.out, recursive = T)

  if(flag.val=="_A"){
    val.st <- 1971; val.end <- 1990
  } else {
    val.st <- 1998; val.end <- 2017
  }

  time.st <- Sys.time()
  #==============================================================================#
  #load data----
  # obsvered Q
  data("Harz_obs_Q")

  # 1925.11 to 2020.12
  Qday <- Harz_obs_Q[[station_id]] %>% rename("year"="iyear","month"="imon","day"="iday","Qd"="Qval") %>%
    subset(year>=1926&year<=2020) # make sure the complete year for resamping
  Qmon <- aggregate(Qd~year+month, Qday, FUN = mean) %>% rename("Qm"="Qd") %>% arrange(year,month)
  Qyr <- aggregate(Qd~year+month, Qday, FUN = mean) %>% rename("Qa"="Qd") %>% arrange(year)
  head(Qmon[,1:2])
  tail(Qmon[,1:2])
  date_Qobs <- as.Date(paste(Qmon[,1],Qmon[,2],"1", sep="-"))
  range(date_Qobs)

  # observed climate
  data("Harz_obs_cli")
  pred.mat <- sapply(Harz_obs_cli, function(ls) ls[[3]])
  #pred.mat <- pred.mat[,-4] # remove Gmon
  head(pred.mat)

  # add Water Balance here
  pred.mat <- pred.mat %>% data.frame() %>% mutate(Bmon = Pmon-Emon)
  colnames(pred.mat)
  head(pred.mat)

  n_lag <- 3; mv <- 12
  pred.all <- pred_mat_lag(pred.mat, lag=n_lag, mv=mv)
  colnames(pred.all)
  summary(pred.all)

  Pmon <- Harz_obs_cli[[1]]
  date_cli <- as.Date(paste(Pmon[,1],Pmon[,2],"1", sep="-"))
  range(date_cli)

  # ### combine GCM data
  # hist_mon <- get(load(paste0("data/Harz_hist_mon.Rdat")))
  # rcp85_mon <- get(load(paste0("data/Harz_rcp85_mon.Rdat")))
  #
  # Harz_gcm <- lapply(1:4, function(i) rbind(hist_mon[[mod_id]][[i]],
  #                                           rcp85_mon[[mod_id]][[i]]))
  # names(Harz_gcm) <- names(hist_mon[[mod_id]])
  # pred.mat <- sapply(Harz_gcm, function(ls) ls[[3]])
  #
  # Pmon <- Harz_gcm$Pmon
  # date_gcm <- as.Date(paste(Pmon[,1],Pmon[,2],"1", sep="-"))
  # range(date_gcm)
  #
  # pred.mat <- pred.mat %>% data.frame() %>% mutate(Bmon = Pmon-Emon)
  # pred.all.cli <- pred_mat_lag(pred.mat, lag=n_lag, mv=mv)
  # colnames(pred.all.cli)
  # head(pred.all.cli)

  # dates----
  date_year <- val.st:val.end

  date_val <- seq(as.Date(paste(val.st,"1","1",sep="-")),
                  as.Date(paste(val.end,"12","1",sep="-")),
                  by="month")
  range(date_val)
  range(date_Qobs); range(date_cli)
  ind_Qobs <- which(date_Qobs %in% date_val)

  # train
  date_cli1 <- date_cli[! (date_cli %in% date_val)]; range(date_cli1)

  ind_Qobs1 <- which(date_Qobs %in% date_cli1); range(date_Qobs[ind_Qobs1])
  ind_cli1 <- which(! (date_cli %in% date_val)); range(date_cli[ind_cli1])

  ind_cli <-  which(date_cli %in% date_val); range(date_cli[ind_cli])
  #ind_gcm <-  which(date_gcm %in% date_val)

  #==============================================================================#
  # knn conditional bootstrap ----
  out_sel <- readRDS(paste0("Harz_",station_id,"_",flag.sel,".rds"))
  cpy <- out_sel[[1]]
  pw <- out_sel[[2]]

  if(TRUE){
    x <- Qmon$Qm[ind_Qobs1]
    z <- pred.all[ind_cli1, cpy]
    zout <- pred.all[ind_cli, cpy]

    x.boot <- Qgen::knn(x, z, zout, pw = pw,reg = FALSE,nensemble = nensemble)
  }

  ## overview----
  summary(x); summary(z)
  par(pty="s")
  x_obs <- Qmon$Qm[which(date_Qobs %in% date_val)]
  x_hat <- x.boot[,1]
  plot(x_obs, x_hat, xlim=c(0,2),ylim=c(0,2))
  abline(a=0,b=1)

  #==============================================================================#
  # Monthly Disaggregation----
  x.sim <- data.frame(Date=date_val, x.boot) %>%
    gather(group, sim,2:(nensemble+1)) %>% mutate_if(is.character,as.factor) %>% data.table()

  x.sim[,month:=month(Date),]
  x.sim[,season := time2season(Date,type="calendar"),]
  x.sim[,year := year(Date),]

  Qmon_sim <- x.sim
  Qsim_year <- Qmon_sim[,.(value=mean(sim)),by=c("year","group")] %>% spread(group, value)
  summary(Qsim_year)

  Qmon1 <- Qmon %>% subset(!(year %in% date_year)); unique(Qmon1$year)
  k <- floor(0.5 + 3 * sqrt(length(unique(Qmon1$year))));k
  Qsim_mon <- lapply(1:nensemble, function(i) knn_annual_to_monthly(Qsim_year[,c(1,i+1)], Qmon1, K=k))
  #summary(Qsim_mon)

  #==============================================================================#
  # Daily Disaggregation----
  Qsim_year <- Qmon_sim[,.(value=mean(sim)),by=c("year","group")] %>% spread(group, value) %>% data.frame()
  summary(Qsim_year)

  Qday1 <- Qday %>% subset(!(year %in% date_year)); unique(Qday1$year)
  k <- floor(0.5 + 3 * sqrt(length(unique(Qday1$year))));k
  Qsim_day <- lapply(1:nensemble, function(i) knn_annual_to_daily(Qsim_year[,c(1,i+1)], Qday1, K=k))
  #summary(Qsim_day)

  #==============================================================================#
  (Sys.time() - time.st) %>% print()

  # Time difference of 52.70773 secs - NPRED
  # Time difference of 23.30932 secs - ALL

  ## Save output----
  if(flag.save){

    save(Qsim_mon,file=paste0(path.out, "Harz_Qmon_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))
    save(Qsim_day,file=paste0(path.out, "Harz_Qday_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))

  }

}
