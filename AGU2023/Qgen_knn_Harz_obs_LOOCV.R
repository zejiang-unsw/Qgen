#==============================================================================#
# By Ze Jiang (ze.jiang@unsw.edu.au)
# Find predictor
# LOOCV

rm(list=ls()); graphics.off() # clear environment and graphics

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
library(hydroTSM)
library(Ecdat)
library(ExtDist)
library(anySim)
library(Qgen)

#==============================================================================#
# global variables----
flag.save <- T # save the Qsim or not
op <- par()
set.seed(20230607)
nensemble <- 100 # must larger than 1 for knn bootstrap

#v0 only val for  # v1 out of val # v2 number of k in disaggregation
flag.ver <- switch(1,"","_v1","_v2")

flag.sel <- switch(2, "ALL","NPRED","WASP") # predictor selection method
dpi_fig <- switch(1, 300, 500) # result quick check
#val.st <- 1971; val.end <- 2000 # validation period
val.st <- 1951; val.end <- 2017 # validation period


time.st <- Sys.time()
#==============================================================================#
#load data----
# obsvered Q
station_id <- "Q5"
data("Harz_obs_Q")

Qday <- Harz_obs_Q[[station_id]] %>% rename("year"="iyear","month"="imon","day"="iday","Qd"="Qval") %>%
  subset(year>=1926&year<=2020) # make sure the complete year for resamping
head(Qday)
tail(Qday)
Qmon <- aggregate(Qd~year+month, Qday, FUN = mean) %>% rename("Qm"="Qd") %>% arrange(year,month)
Qyr <- aggregate(Qd~year+month, Qday, FUN = mean) %>% rename("Qa"="Qd") %>% arrange(year)
head(Qmon[,1:2])
tail(Qmon[,1:2])
date_Qobs <- as.Date(paste(Qmon[,1],Qmon[,2],"1", sep="-"))
range(date_Qobs)

# obsvered climate
data("Harz_obs_cli")
pred.mat <- sapply(Harz_obs_cli, function(ls) ls[[3]])
#pred.mat <- pred.mat[,-4] # remove Gmon
head(pred.mat)

# add Water Balance here
pred.mat <- pred.mat %>% data.frame() %>% mutate(Bmon = Pmon-Emon)
colnames(pred.mat)
head(pred.mat)

Pmon <- Harz_obs_cli[[1]]
date_cli <- as.Date(paste(Pmon[,1],Pmon[,2],"1", sep="-"))
range(date_cli)

# dates----
date_year <- val.st:val.end

date_val <- seq(as.Date(paste(val.st,"1","1",sep="-")),
                  as.Date(paste(val.end,"12","1",sep="-")),
                  by="month")

ind_Qobs <- which(date_Qobs %in% date_val)
ind_Qobs1 <- which(date_Qobs %in% date_cli)

ind_cli <-  which(date_cli %in% date_val)

# output----
path.out <- paste0(station_id,"/cli_obs/")

if (!dir.exists(path.out)){
  dir.create(path.out, recursive = T)
} else {
  print("Dir already exists!")
}
#==============================================================================#
# predictor selection----
n_lag <- 3; mv <- 12
pred.mat.n <- pred_mat_lag(pred.mat, lag=n_lag, mv=mv)
summary(pred.mat.n)

mat <- cbind(Q=Qmon$Qm[ind_Qobs1], pred.mat.n) %>% na.omit()
#mat <- cbind(Q=Qmon$Qm[ind_Qobs], pred.mat.n[ind_cli, ]) %>% na.omit()
summary(mat)
colnames(mat)

## check data----
plot.ts(mat[,1:10], main=NA, xlab=NA)
pairs(mat[,c(1, seq(2,by=5,length.out=ncol(pred.mat)))])

x <- mat[,1]
py<- mat[,-1]

if(flag.sel=="NPRED"){
  out <- NPRED::stepwise.PIC(x=x, py=py)
} else if(flag.sel=="WASP") {
  out <- WASP::stepwise.VT(list(x=x, dp=py), wf="haar", mode="MODWT",J=5)
  #out <- WASP::stepwise.VT(list(x=x, dp=py), wf="d4", mode="MRA", J=5)
} else if(flag.sel=="ALL"){
  #out <- list(cpy=1:ncol(py), wt=rep(1, ncol(py)))

  cpy <- seq(1,by=5,length.out=ncol(pred.mat))
  wt <- rep(1, ncol(pred.mat))
  out <- list(cpy=cpy, wt=wt)
}

if(is.null(out$cpy)){
  cpy <- seq(1,by=5,length.out=ncol(pred.mat))
  wt <- rep(1, ncol(pred.mat))
} else {
  cpy <- out$cpy
  wt <- out$wt
}
print(cpy); print(wt)
colnames(pred.mat.n)[cpy] %>% print()

## Save NPERD----
out <- list(cpy, wt)
saveRDS(out,file=paste0("Harz_",station_id,"_",flag.sel,".rds"))

#==============================================================================#
# knn conditional bootstrap ----
if(TRUE){ # cross vadidation for target period
x <- Qmon$Qm[ind_Qobs1]
z <- pred.mat.n[, cpy]
zout <- pred.mat.n[ind_cli, cpy]

x.boot <- t(sapply(1:length(date_val), function(i) knn(x[-ind_cli[i]],
                                                         z[-ind_cli[i],], zout[i,],
                                                         pw = wt,reg = FALSE,nensemble = nensemble)))
} else { # cross validation for all period
x <- Qmon$Qm[ind_Qobs]
z <- pred.mat.n[ind_cli, cpy]

x.boot <- t(sapply(1:length(date_val), function(i) knn(x[-i],z[-i,], z[i,],
                                                       pw = wt,reg = FALSE,nensemble = nensemble)))
}

## Visual check----
plot(x[ind_cli],x.boot[,1], xlab="obs",ylab="sim")
for(i in 2:nensemble){
  points(x[ind_cli], x.boot[,i],col=i)
}

plot.ts(cbind(x.boot[,1:min(9,nensemble)],x), plot.type = "single",
        col=c(rep("grey",min(9,nensemble)),"red"),
        lwd=c(2,rep(1,min(9,nensemble))),xlab=NA, ylab="Qmon")
legend("topright",legend=c("Sim", "Obs"),
       lwd=c(1,1),bg="transparent",bty = "n",
       col=c("grey","red"),cex=0.8,horiz=TRUE)

#==============================================================================#
# Monthly----
x.sim <- data.frame(Date=date_val, obs=Qmon$Qm[ind_Qobs], x.boot) %>%
  gather(group, sim,3:(nensemble+2)) %>% mutate_if(is.character,as.factor) %>% data.table()

x.sim[,month:=month(Date),]
x.sim[,season := Qgen::time2season(Date,type="calendar"),]
x.sim[,year := year(Date),]

summary(x.sim)

## monthly climatology----
x.mon.mean <- x.sim[,.(sim=mean(sim),
                       obs=mean(obs)),by=c("month","group")]
summary(x.mon.mean)

p.mon.clim1 <- ggplot(x.mon.mean,aes(x=factor(month), y=sim)) +

  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black") +

  geom_point(aes(x=factor(month), y=obs, color="Qobs"),shape=17) +

  #scale_x_discrete(labels=1:lag.max) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values="red") +

  labs(x="Month",y="Streamflow",title="Monthly Climatology",
       color=NULL, fill=NULL) +
  guides(fill="none") +
  theme_bw() +
  theme(text = element_text(size = 16,face='bold',family="serif"),
        plot.margin = unit(c(0.1,0.1,0.1, 0.1), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        # strip.background = element_blank(),
        # strip.placement = "outside",

        legend.position = "right",
        #legend.position = c(.1,.85),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.title=element_blank()
  )

p.mon.clim1 %>% print()

### save plot----
ggsave(paste0(path.out, "Figure_",flag.sel,"_mon_clim",flag.ver,".jpg"),p.mon.clim1,
       width = 9, height = 7, dpi=dpi_fig)

## seasonal climatology ----
Qsea <- Qmon %>% mutate(season=cut(Qmon$month, breaks=c(0, 3, 6, 9, 12), right = T,
                                  labels = c(1,2,3,4))) %>%
  aggregate(Qm~year+season, FUN=mean) %>% data.table()
summary(Qsea)

Qsea.mean <- Qsea[,.(obs=mean(Qm)),by=c("season")]
summary(Qsea.mean)

summary(x.sim)
x.sea.mean_r <- x.sim[,.(value=mean(sim)),by=c("season","year","group")]
x.sea.mean_r1 <- x.sea.mean_r[,.(mean=mean(value)),by=c("season","group")]
summary(x.sea.mean_r1)

p.sea.clim1 <- ggplot(x.sea.mean_r1,aes(x=factor(season), y=mean)) +

  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black") +

  geom_point(data=Qsea.mean,aes(x=season, y=obs, color="Qobs"),shape=17) +

  #scale_x_discrete(labels=label_sea) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_color_manual(values="red") +

  labs(x="Season",y="Streamflow",title="Seasonal Climatology",
       color=NULL, fill=NULL) +
  guides(fill="none") +
  theme_bw() +
  theme(text = element_text(size = 16,face='bold',family="serif"),
        plot.margin = unit(c(0.1,0.1,0.1, 0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        # strip.background = element_blank(),
        # strip.placement = "outside",

        legend.position = "right",
        #legend.position = c(.1,.85),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.title=element_blank()
  )

p.sea.clim1 %>% print()

### save plot----
ggsave(paste0(path.out,"Figure_",flag.sel,"_sea_clim",flag.ver,".jpg"),
       width = 9, height = 7, dpi=dpi_fig)

## detrend and deseason----
x_n <- deseason_trend(x, 12, do.plot = F)
x.boot_n <- apply(x.boot, 2, function(vec) deseason_trend(vec, 12, do.plot = F))

## month-to-month correlation----
x.mon.mat <- lapply(1:nensemble, function(r) matrix(x.boot[,r], ncol=12,byrow=T))
seasonCor_mat <- sapply(1:nensemble, function(r) s2scor(x.mon.mat[[r]]))
seasonCor_obs <- s2scor(matrix(x, ncol=12, byrow=T));seasonCor_obs

mon2mon_cor <- data.frame(month=factor(1:12), obs=seasonCor_obs, seasonCor_mat) %>%
  gather(group,value, 3:(nensemble+2))
summary(mon2mon_cor)

p.mon.cor <- ggplot(mon2mon_cor,aes(x=month, y=value)) +

  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black") +

  geom_point(data=mon2mon_cor, aes(x=month, y=obs, color="obs"),shape=17) +

  scale_y_continuous(limits=c(NA,1), breaks = pretty_breaks()) +
  scale_color_manual(values="red") +

  labs(x="Month",y="Correlation",title="Month to Month", color=NULL, fill=NULL) +
  guides(fill="none") +
  theme_bw() +
  theme(text = element_text(size = 16,face='bold',family="serif"),
        plot.margin = unit(c(0.1,0.1,0.1, 0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        # strip.background = element_blank(),
        # strip.placement = "outside",

        legend.position = "right",
        #legend.position = c(.1,.85),
        legend.direction = c("horizontal","vertical")[2],
        legend.background = element_rect(fill="transparent"),
        legend.title=element_blank()
  )

p.mon.cor %>% print()

### save plot----
ggsave(paste0(path.out,"Figure_",flag.sel,"_mon_cor",flag.ver,".jpg"), p.mon.cor,
       width = 9, height = 7, dpi=dpi_fig)

## Lag l autocorrelation----
flag.acf <- "ACF"; lag.max <- 12
if(flag.acf=="ACF"){
  acf.x <-  acf(x_n, lag.max = lag.max, plot = F)$acf[-1,1,1]
  acf.xboot <- apply(x.boot_n, 2, function(vec) acf(vec, lag.max = lag.max, plot = F)$acf[-1,1,1])

} else {
  acf.x <-  pacf(x_n, lag.max = lag.max, plot = F)$acf[,1,1]
  acf.xboot <- apply(x.boot_n, 2, function(vec) pacf(vec, lag.max = lag.max, plot = F)$acf[,1,1])

}

df.acf.mon <- data.frame(Lag=paste0("lag",1:lag.max), obs=acf.x, acf.xboot) %>%
  gather(nensemble, sim, 3:(ncol(acf.xboot)+2)) %>%
  mutate_if(is.character, as.factor)
summary(df.acf.mon)

df.acf.mon$Lag <- factor(df.acf.mon$Lag, levels=paste0("lag",1:lag.max))
lim_y_acf <- c(-0.5,0.5)

p.acf.mon <- ggplot(df.acf.mon, aes(x=Lag, y=sim)) +

  geom_boxplot(alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black", fill="black") +

  geom_point(aes(x=Lag, y=obs, color="Qobs"),shape=17) +

  scale_x_discrete(labels=1:lag.max) +
  scale_y_continuous(limits=lim_y_acf, breaks=pretty_breaks()) +
  scale_color_manual(values="red") +

  labs(x="Lag (month)",y=flag.acf,title="Monthly autocorrelation",
       color=NULL, fill=NULL) +
  guides(fill="none") +

  theme_bw() +
  theme(text = element_text(size = 16,face='bold',family="serif"),
        plot.margin = unit(c(0.1,0.1,0.1, 0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        # strip.background = element_blank(),
        # strip.placement = "outside",

        #legend.position = "right",
        legend.position = c(.8,.9),
        #legend.direction = c("horizontal","vertical")[2],
        #legend.background = element_rect(fill="transparent"),
        legend.title=element_blank()
  )
p.acf.mon %>% print()

### save plot----
ggsave(paste0(path.out,"Figure_",flag.sel,"_mon_",flag.acf,flag.ver,".jpg"),p.acf.mon,
       width = 9, height = 7, dpi=dpi_fig)

## Distribution----
summary(x.sim)
p.dis.mon <- ggplot(x.sim) +
  # stat_qq(aes(sample=x.sim, fill=group, color="Simulation")) +
  # #stat_qq_line(aes(sample=x.sim, fill=group, color="Simulation")) +
  #
  # stat_qq(aes(sample=x, color="Observed")) +
  # stat_qq_line(aes(sample=x, color="Observed")) +

  stat_ecdf(aes(x=sim, fill=group,color="Simulation"), pad = TRUE) +
  stat_ecdf(aes(x=obs, color="Observed"),pad = TRUE) +

  #geom_point(aes(color = id), size = 1.5, alpha = 0.3, show.legend = FALSE) +
  #geom_smooth(method = "lm", formula = y~x) +

  scale_color_manual(values=c("black","red")) +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks=pretty_breaks(), expand = c(0,0)) +

  labs(x="Streamflow", y="CDF", title="Monthly",color=NULL) +
  guides(fill='none') +
  theme_bw(base_size = 12) +
  theme(text = element_text(size = 16),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(1.2, "lines"),

        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position="right",
        legend.position = c(.8,.25),
        legend.background = element_rect(fill="transparent")
  )

p.dis.mon %>% print()

### save plot----
ggsave(paste0(path.out,"Figure_",flag.sel,"_mon_cdf",flag.ver,".jpg"),
       width = 9, height = 7, dpi=dpi_fig)

#==============================================================================#
# Monthly Disaggregation----
Qmon_sim <- x.sim
Qsim_year <- Qmon_sim[,.(value=mean(sim)),by=c("year","group")] %>% spread(group, value)
summary(Qsim_year)

#Qmon1 <- Qmon
Qmon1 <- Qmon %>%  subset(year>=val.st&year<=val.end) # use only the validation period can improve general results but the min/max might not
summary(Qmon1)

k <- floor(0.5 + 3 * sqrt(length(unique(Qmon1$year))));k
Qsim_mon <- lapply(1:nensemble, function(i) knn_annual_to_monthly(Qsim_year[,c(1,i+1)], Qmon1, K=k, flag=F))
#summary(Qsim_mon)

# overview
year.target <- Qsim_mon[[1]]$year
year.sample <- Qsim_mon[[1]]$nn

#plot(year.target, year.sample)
sum(year.sample==year.target) %>% print()

#==============================================================================#
# Daily Disaggregation----
Qsim_year <- Qmon_sim[,.(value=mean(sim)),by=c("group","year")] %>% spread(group, value) %>% data.frame()
summary(Qsim_year)


Qday1 <- Qday %>%  subset(year>=val.st&year<=val.end) # use only the validation period can improve results but the min/max might not

summary(Qday1)

if(TRUE){ # disaggredate from yearly
  k <- floor(0.5 + 3 * sqrt(length(unique(Qday1$year))));k
  Qsim_day <- lapply(1:nensemble, function(i) knn_annual_to_daily(Qsim_year[,c(1,i+1)], Qday1, K=k, flag=F))
  summary(Qsim_day)

} else { # disaggredate from monthly
# Two option here: option 2 is better than the option 1 but both worse than the above one.
# option 1: use predicted Qmon_sim; option 2: used predicted and rescaled Qsim_mon
# Qsim_mon1 <- Qmon_sim[,.(value=mean(sim)),by=c("group","year","month")] %>% spread(group, value) %>% data.frame()
# summary(Qsim_mon1)

  Qsim_mon_df <- rbindlist(Qsim_mon, idcol="group") %>% rename("sim"="Qm")
  Qsim_mon1 <- Qsim_mon_df[,.(value=mean(sim)),by=c("group","year","month")] %>% spread(group, value) %>% data.frame()
  summary(Qsim_mon1)

  k <- floor(0.5 + 3 * sqrt(length(unique(Qday1$year))));k
  Qsim_day <- lapply(1:nensemble, function(i) knn_monthly_to_daily(Qsim_mon1[,c(1:2,i+2)], Qday1, K=k))

}

# overview
year.target <- Qsim_day[[1]]$year
year.sample <- Qsim_day[[1]]$nn

sum(year.sample==year.target) %>% print()

#==============================================================================#
(Sys.time() - time.st) %>% print()

# Time difference of 52.70773 secs - NPRED
# Time difference of 23.30932 secs - ALL

## Save output----
if(flag.save){

  save(Qsim_mon,file=paste0(path.out, "Harz_Qmon_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))
  save(Qsim_day,file=paste0(path.out, "Harz_Qday_",flag.sel,"_r",nensemble,flag.ver,".Rdat"))

}

