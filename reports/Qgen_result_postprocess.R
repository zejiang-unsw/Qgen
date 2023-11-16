rm(list=ls()); graphics.off() # clear environment and graphics

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path)); #setwd("../")

library(dplyr)
library(tidyr)
library(data.table)
library(lubridate) # leap_year

flag.save <- F

nensemble <- 100
flag.sel <- switch(2, "ALL","NPRED","WASP"); flag.sel
flag.ver <- switch(1, "","_v1","_v2","_v3")

# Demo station
station_id <- "Q5"
mod_id <- 7 # GCM realization: 1-6 is raw, and 7 is the sampled scenario

## validation----
if(T){
path.out <- paste0(station_id,"/cli_obs/")
print("validation")
load(paste0(path.out, "Harz_Qmon_",flag.sel,"_r",nensemble,"_val",flag.ver,".Rdat")) # Qsim_mon
load(paste0(path.out, "Harz_Qday_",flag.sel,"_r",nensemble,"_val",flag.ver,".Rdat")) # Qsim_day
summary(Qsim_mon[[1]])
summary(Qsim_day[[1]])

# overview
year.target <- Qsim_mon[[1]]$year
year.sample <- Qsim_mon[[1]]$nn

sum(year.sample==year.target) %>% print()

year.target <- Qsim_day[[1]]$year
year.sample <- Qsim_day[[1]]$nn

sum(year.sample==year.target) %>% print()

Qmon_sim <- rbindlist(Qsim_mon, idcol="group") %>% select(!nn) %>% spread(group, Qm)
Qday_sim <- rbindlist(Qsim_day, idcol="group") %>% select(!nn) %>% spread(group, Qd)

### remove 29 Feb
ind_29 <- which(Qday_sim$month==2&Qday_sim$day==29)
Qday_sim[ind_29, 4];Qday_sim[ind_29-1, 4]
Qday_sim[ind_29-1, -c(1:3)] <- Qday_sim[ind_29-1, -c(1:3)] + Qday_sim[ind_29, -c(1:3)]

Qday_sim <- Qday_sim[-ind_29,]

length(unique(Qmon_sim$year))*365==nrow(Qday_sim)

if(flag.save){

  path.out <- paste0(station_id,"_v2/cli_obs/")

  if (!dir.exists(path.out)){
    dir.create(path.out, recursive = T)
  } else {
    print("Dir already exists!")
  }

  save(Qmon_sim,file=paste0(path.out, "Harz_Qmon_",flag.sel,"_r",nensemble,"_val",flag.ver,".Rdat"))
  save(Qday_sim,file=paste0(path.out, "Harz_Qday_",flag.sel,"_r",nensemble,"_val",flag.ver,".Rdat"))

}
}
## future simulation----
out.Qmon <- list()
gcms <- c("hist","rcp85_near", "rcp85_far","pseudo")[4:4]
for(i in 1:length(gcms)) {
  # i <- 1
  flag.gcm <- gcms[i]
  print(flag.gcm)
  path.out <- paste0(station_id,"/cli_mod_",flag.gcm,"/mod_",mod_id,"/")
  Qsim_mon <- get(load(paste0(path.out, "Harz_Qmon_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))) # Qsim_mon
  Qsim_day <- get(load(paste0(path.out, "Harz_Qday_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))) # Qsim_day

  # overview
  year.target <- Qsim_mon[[1]]$year
  year.sample <- Qsim_mon[[1]]$nn
  sum(year.sample==year.target) %>% print()

  year.target <- Qsim_day[[1]]$year
  year.sample <- Qsim_day[[1]]$nn
  sum(year.sample==year.target) %>% print()

  Qmon_sim <- rbindlist(Qsim_mon, idcol="group") %>% select(!nn) %>% spread(group, Qm)
  Qday_sim <- rbindlist(Qsim_day, idcol="group") %>% select(!nn) %>% spread(group, Qd)

  ### remove 29 Feb
  ind_29 <- which(Qday_sim$month==2&Qday_sim$day==29)
  Qday_sim[ind_29, 4];Qday_sim[ind_29-1, 4]
  Qday_sim[ind_29-1, -c(1:3)] <- Qday_sim[ind_29-1, -c(1:3)] + Qday_sim[ind_29, -c(1:3)]

  Qday_sim <- Qday_sim[-ind_29,]

  length(unique(Qmon_sim$year))*365==nrow(Qday_sim)

  if(flag.save){

    #path.out <- paste0(station_id,"_v2/cli_mod_",flag.gcm,"/mod_",mod_id,"/")
    path.out <- paste0(station_id,"_v2/cli_mod_",flag.gcm,"/")
    if (!dir.exists(path.out)){
      dir.create(path.out, recursive = T)
    } else {
      print("Dir already exists!")
    }

    save(Qmon_sim,file=paste0(path.out, "Harz_Qmon_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))
    save(Qday_sim,file=paste0(path.out, "Harz_Qday_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))

  }
}

