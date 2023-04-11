# driver script to set up the system to generate (retrospective) forecast for select US states
# 7/28/22: using variant specific case/death counts (apportioned per GISAID data)

dir_data = './data/'
dir_code = './code/'
dir_res = paste0('./results/test/')
if(! file.exists(dir_res))  dir.create(dir_res,recursive = T)

source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code,'SEIRS.R'))
source(paste0(dir_code,'EAKF_rproj_scenarios.R'))
source(paste0(dir_code, 'get_fcastProbDist.R'))
source(paste0(dir_code,'get_relR0.R'))

# DUMMY VARIABLES TO SET FOR DIFFERENT MODEL/LOCATION/FORECAST WEEK, BEFORE RUNNING 
if(F){
  # for test run
  dummy_isn = 1
  dummy_iloc = 1
  dummy_tno = 1
  dummy_idf = 1 # for deflation
  
  iwk_fcast = 1
  iwk_fcast = 90
}

dtag = '_upto20221002' # time stamp for data used for the study


fcast.start_vec_non.Omicron = c(seq(as.Date('2020/7/1'), as.Date('2021/8/15') , by = '1 week')
)
fcast.start_vec_Omicron =  seq(as.Date('2021/12/1'), as.Date('2022/9/30') , by = '1 week') # t32
fcast.start_vec = c(fcast.start_vec_non.Omicron, fcast.start_vec_Omicron)
length(fcast.start_vec)


vtags = c('non.Omicron', 'Omicron')
fcast.deflat_vec = seq(.9, 1, by = .025)
v.deflat = c('S1', 'E1', 'I1') # only adjust the state variables

{
  # do not use
  # tested and found these were not good 
  exclExtreme = F # exclude the extreme values before projection, for long term outcome
  perc2excl = .1 # fraction of extreme values to exclude
  useMean = F # use the mean for projection
  useMean.wd = .2 # use the mean for projection but widen it by +/- SD * useMean.wd
  
  tag.proj = ifelse(useMean, paste0('useMean',useMean.wd,'SD'),
                    ifelse(exclExtreme, paste0('exclExtreme',perc2excl * 100, 'perc'), 
                           'useEns'))
  
}

tag.model = ''
num_runs = 100  # run more


# seasonality options
tab_sn = data.table(seasonality = c(T, T, F), 
                    useTranRelR0 = c(T, F, F),
                    sn.tag = c('w.sn.tran','w.sn', 'no.sn'))

# specify dummy_isn for type of seasonality to use
seasonality = tab_sn[dummy_isn]$seasonality; 
sn.tag = tab_sn[dummy_isn]$sn.tag
useTranRelR0 = tab_sn[dummy_isn]$useTranRelR0 # whether to use transformed seasonal trend (relR0)

tda = NULL
for(v.tag in c('non.Omicron','Omicron')){
  tmp = read.csv(paste0(dir_data, 'da_case_death_us_nyt_',v.tag, dtag,'.csv')) %>% data.table()
  tda = rbind(tda, tmp, fill=T)
}
week.starts = tda$date %>% as.Date %>% unique
fcast.start_vec_non.Omicron = week.starts[week.starts >= as.Date('2020/7/1') & week.starts <= as.Date('2021/8/15')]
fcast.start_vec_Omicron = week.starts[week.starts >= as.Date('2021/12/1')] 
fcast.start_vec_Omicron = append(fcast.start_vec_Omicron, max(fcast.start_vec_Omicron)+7)
fcast.start_vec = c(fcast.start_vec_non.Omicron, fcast.start_vec_Omicron)
length(fcast.start_vec)

fcast.start.month_vec = format(fcast.start_vec, '%Y%m%d')
fcast.start.month = fcast.start.month_vec[iwk_fcast]  # this may cause problem at the end when the last week of data match with multiple fcast.start from the which.min function
fcast.start = as.Date(fcast.start_vec[iwk_fcast]) #  %>% as.Date #  as.Date('2020/09/01')

# fcast.start_vec
fcast.start

if(fcast.start < as.Date('2021/11/21')){
  variant.tag = 'non.Omicron'
} else {
  variant.tag = 'Omicron'
}

print(variant.tag)


vax.start = as.Date('2020/12/14')
date.start = as.Date('2020/03/01') # start from 3/1/20
mob.type = 'business'

# for setting seeding
case.us = read.csv(paste0(dir_data, 'da_case_death_us_national.csv')) %>% data.table() # us national level case data, 7-day moving average
case.us$date = case.us$date %>% as.Date

DAT.EPI = read.csv(paste0(dir_data, 'da_case_death_us_nyt_',variant.tag, dtag,'.csv')) %>% data.table()
DAT.EPI$date = DAT.EPI$date %>% as.Date


# read mobility data
DAT.MOB = read.csv(paste0(dir_data, 'da_mobility_us.csv')) %>% data.table()
DAT.MOB = DAT.MOB[data.type == mob.type]

# read vaccination data
DAT.VAC = read.csv(paste0(dir_data,'da_vx_perM_us_lagged_booster.csv')) %>% data.table()

# read seasonal trend
DAT.SN = read.csv(paste0(dir_data, 'est.sn.2000t2020_us.csv')) %>% data.table()
DAT.SNraw = read.csv(paste0(dir_data, 'est.raw.sn.2000t2020_us.csv')) %>% data.table()

# read variant data, for detecting rising new variants
DAT.VARIANT = read.csv(paste0(dir_data, 'perc.variant.wkly.smoothed_by_state.csv')) %>% data.table()
DAT.VARIANT$week = DAT.VARIANT$week %>% as.Date
# read the dates key variants became predominant in different states - for setting VE per variant
vdates = read.csv(paste0(dir_data, 'dates_variant_by_state.csv')) %>% data.table()

# the list of states 
locs = DAT.EPI$state %>% unique()

N = 1e6; # per 1 M

num_gr = num_obs = length(N); # no age structure
num_ens = 500

# epi.model = 'SEIRSV' # susceptible-exposed-infectious-recovered-susc
epi.model = 'SEIRSVimmLoss'  # track immune loss of those vaccinated seperately

stoch = T

source(paste0(dir_code,'set_tm2event.R'))


seed = .1

# set parms for space reprobing (SR)
doSR = T
percSRmajor = .1
SR.perc = .05
SR.perc.local = .03;
# SR.perc.full = .07 -> large bounds t13
SR.perc.full = .03; # 0.03;
SR.perc.imm = .06; # for immune evasion, probing on S
SR.perc.extra = .15; 
SR.perc.Extra = .33; # for decrease in detection 
SR.perc.Extra1 = .15; # for increase in detection rate
SR.perc.Extra2 = .9; # for decrease in detection 
SR.perc.Extra3 = .6;
SR.perc.EXTRA = .95; # RE-NEW THE SYSTEM
SR.perc.extraS = .25 # for super-spreading? 

SR.var.local= c('beta','Tei','Tir','Trs','Td.mean','Td.sd','p.mob','alpha','ifr')
# SR.var.full= c('beta','Tir','ifr','alpha') 
# 3/4/22 - add p.mob in SR.var.full
SR.var.full= c('beta','Tir','ifr','alpha','p.mob') # ,'alpha','ifr' 
SR.var.tx = c('beta','Tir')

donotUpdateS = T # do not allow the filter to update S during first wave

rednUpdateEI = T # do not allow or reduce the level allowed, the filter to update E or I - it takes OEVr all efforts

redn.priority = 1  # no prioritization of vac

# for omicron
tm2death_adj_omicron = 10 # t32 used 10 days 5 # assume the average time to death is 5 days longer for omicron
tm.to.death.max_omicron = 90 
# NYC data: mean = 14.41, sd = 17.15, excluded posthumous (25 - 40%), so likely shorter
# check time lags

# parameter prior ranges
PRAMBOUNDS = read_excel(paste0(dir_code, 'parm.bounds_pop.den_lowerbeta.xlsx'), sheet = 1) %>% data.table() # adjust beta per pop density slightly 
PRAMBOUNDS$date.start = PRAMBOUNDS$date.start %>% as.Date()
PRAMBOUNDS$date.end = PRAMBOUNDS$date.end %>% as.Date()
PRAMBOUNDS$lwr = PRAMBOUNDS$lwr %>% as.numeric()
PRAMBOUNDS$upr = PRAMBOUNDS$upr %>% as.numeric()

# intended for adjusting initial OEV but not used
# date_ini = as.Date('2020/4/1') 
# date_omicron = as.Date('2021/11/15')

# estimated from VE data
est_parmVimmLoss = read.csv(paste0(dir_code, 'est_parmVimmLoss.csv')) %>% data.table()

# best parameter settings for the transformed seasonality
sn.tran.parms = read.csv(paste0(dir_code,'best_sn.tran.parms.csv')) %>% data.table()


loc.t = 'New York'
loc.t = 'Florida'
loc.t = 'Michigan'
loc.t = "California" 
loc.t = 'Massachusetts'
loc.t = 'Texas'
loc.t = 'Wyoming'
loc.t = 'Iowa'
loc.t = 'Alaska'
locs.ab.t = c('MA', 'NY', 'PA', 'MI', 'TX', 'IA', 'WY', 'CA', 'WA', 'FL')
locs = c(state.name, "District of Columbia")
names(locs) = c(state.abb, 'DC')
locs.t = locs[names(locs) %in% locs.ab.t]

week.starts = DAT.EPI[date >= date.start]$date %>% unique %>% as.Date %>% sort


# fcast.start = week.starts[which.min(abs(week.starts - fcast.start))]
nfcast = 26 # number weeks forecast done
date.end = fcast.start + (nfcast - 1) * 7
iwk_fcast
fcast.start.month_vec[iwk_fcast]
fcast.start


for(loc.t in locs.t[dummy_iloc]){  # dummy_iloc: for running on a cluster
  
  if(length(fcast.start) < 1) 
    break
  
  print(loc.t)
  print(paste0('seasonality: ', seasonality), quote=F)
  print(paste0('fcast.start: ', fcast.start), quote=F)
  
  
  # get the best parm settings for the transformed seasonality
  sn.tran.parm = sn.tran.parms %>% dplyr::filter(loc == loc.t)
  sn.relR0min.adj = sn.tran.parm$sn.relR0min.adj
  sn.dur0 = sn.tran.parm$sn.dur0;
  sn.pshift = sn.tran.parm$sn.pshift;
  
  
  da.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
  
  
  if(fcast.start - max(da.t$date) > 7)
    break
  
  da.mob.t = DAT.MOB[state == loc.t & data.type == mob.type] %>% setnames('value', 'mob')
  da.t$date = da.t$date %>% as.Date()
  da.mob.t$date = da.mob.t$date  %>% as.Date()
  da.t = merge(da.t, da.mob.t, x.all = T, by = c('date', 'year', 'week'))
  da.t = da.t[date >= date.start]
  da.t = da.t[date <= date.end]
  da.t = da.t[order(date)]
  
  rel.mob = da.t$mob %>% as.matrix()
  
  da.vacc = DAT.VAC[state == loc.t]
  
  # seasonality
  relR0 = DAT.SN[state == loc.t] %>% .[order(week)]
  relR0 = matrix(relR0$value, nrow = nrow(relR0), ncol = num_ens)
  
  # 4/8/22 get the transformed relR0
  relR0tran = fn_getR0trans(relR0=DAT.SN[state == loc.t] %>% .[order(week)] %>% .$value, # 8/10/22 used the capped relR0 so that it won't be too extreme
                            num_ens = num_ens,
                            relR0uncapped = DAT.SNraw[state == loc.t] %>% .[order(week)] %>% .$value,
                            bound.relR0min = sn.relR0min.adj+c(.8, 1.2),
                            bound.dur0 = sn.dur0 + c(-4, 4),
                            bound.pshift = sn.pshift + c(-2, 2))
  
  if(useTranRelR0){
    relR0 = relR0tran
  }
  
  VE1wt = .85
  VE2wt = .95  # higher b/c we are using mortality data too
  
  VE1delta = .5
  VE2delta = .8  # higher b/c we are using mortality data too
  
  VE1omicron = .1 # .5  # 1st dose
  VE2omicron = .7  # 2nd and 3rd dose
  
  # date a given variant became predominant - for setting VE
  # date.delta = vdates[state == loc.t & variant == 'Delta']$date  %>% as.Date()
  # date.omicron = vdates[state == loc.t & variant == 'Omicron']$date  %>% as.Date()
  
  date.delta = vdates[state == loc.t & variant == 'Delta']$date.dominant  %>% as.Date()
  date.omicron = vdates[state == loc.t & variant == 'Omicron']$date.dominant  %>% as.Date()
  
  # for tuning the vaccine-induced immunity against infection
  # for wildtype, alpha
  tm.imm.wt = est_parmVimmLoss[variant == 'wt']$tm.imm; # during of vaccine-incuded protection against infection
  # tm.ini.imm.wt = 4 * 30; # initial period with near 0 imm loss
  p.imm.wane.max.wt =  est_parmVimmLoss[variant == 'wt']$p.imm.wane.max; # maximal level of immunity loss (=1 or lower)
  k.wt =  est_parmVimmLoss[variant == 'wt']$k; 
  
  tm.imm.delta =  est_parmVimmLoss[variant == 'delta']$tm.imm; # during of vaccine-incuded protection against infection
  # tm.ini.imm.delta = 3 * 30; # initial period with near 0 imm loss
  p.imm.wane.max.delta =  est_parmVimmLoss[variant == 'delta']$p.imm.wane.max; # maximal level of immunity loss (=1 or lower)
  k.delta =  est_parmVimmLoss[variant == 'delta']$k;
  
  tm.imm.omicron =  est_parmVimmLoss[variant == 'omicron']$tm.imm; # during of vaccine-incuded protection against infection
  # tm.ini.imm.omicron = 30; # initial period with near 0 imm loss
  p.imm.wane.max.omicron = est_parmVimmLoss[variant == 'omicron']$p.imm.wane.max; # maximal level of immunity loss (=1 or lower)
  k.omicron = est_parmVimmLoss[variant == 'omicron']$k;
  
  
  seed_max = (da.t$case[1] + .1) * 100 # 10000
  tm_largerVar = 5 # number of initial weeks to have larger OEV
  pOEV = 1
  

  p.mob_bounds = c(.5, 1.5); # scaling for mobility
  Td.mean_bounds =  c(5,8) # mean Td: reporting delay
  Td.mean_SRbounds = c(5, 7)
  Td.sd_bounds = c(1,3) # Td, sd: reporting delay sd 
  
  # imm_bounds = c(2, 3) * 365
  if(variant.tag == 'non.Omicron'){
    imm_bounds = c(2, 3) * 365
  } else if (variant.tag == 'Omicron'){
    imm_bounds = c(1, 3) * 365
  }
  
  
  # read the parm bounds for different stage
  parm.bound_vec = PRAMBOUNDS[location == loc.t] %>% data.table()
  if(variant.tag == 'non.Omicron'){
    beta_bounds = parm.bound_vec[parm == 'beta' & type == 'initialization'] %>% .[,c('lwr','upr')] %>% unlist 
    alpha_bounds = parm.bound_vec[parm == 'alpha' & type == 'initialization'] %>% .[,c('lwr','upr')] %>% unlist 
    ifr_bounds = parm.bound_vec[parm == 'ifr' & type == 'initialization'] %>% .[,c('lwr','upr')] %>% unlist
  } else if (variant.tag == 'Omicron'){
    # print('line 412')
    beta_bounds = parm.bound_vec[parm == 'beta' & type == 'SR' & period == 'Omicron_BA.1'] %>% .[,c('lwr','upr')] %>% unlist 
    alpha_bounds = parm.bound_vec[parm == 'alpha' & type == 'SR' & period == 'Omicron0'] %>% .[,c('lwr','upr')] %>% unlist 
    # print('line 415')
    ifr_bounds = parm.bound_vec[parm == 'ifr' & type == 'SR' & period == 'Omicron_BA.1'] %>% .[,c('lwr','upr')] %>% unlist
    
  }
  
  # do not use the bounds for Omicron if the data are for non-Omicron variants
  if(variant.tag == 'non.Omicron'){
    parm.bound_vec = parm.bound_vec %>% dplyr::filter(!grepl('omicron', tolower(type)) & 
                                                  !grepl('omicron', tolower(parm)) &
                                                !grepl('omicron', tolower(variant)) &
                                                  !grepl('omicron', tolower(period))) %>% data.table()
    # also update the end time stamp for the Delta wave
    parm.bound_vec %>% dplyr::filter(grepl('delta', tolower(type)) | 
                                grepl('delta', tolower(parm)) |
                                grepl('delta', tolower(variant)) |
                                grepl('delta', tolower(period)))
    parm.bound_vec %>% dplyr::filter(parm=='alpha' & type == 'SR' & period == 'Delta') %>% 
      filter(date.end == max(date.end)) -> toupdate; # = as.Date('2022/12/31')
    # toupdate$date.end = as.Date('2022/12/31')
    parm.bound_vec = parm.bound_vec %>% mutate(date.end = case_when(parm=='alpha' & type == 'SR' & date.end == toupdate$date.end ~ as.Date('2022/12/31'),
                                                   T ~ date.end))
    
    parm.bound_vec %>% filter(parm=='beta' & type == 'SR' & period == 'Delta') %>% 
      filter(date.end == max(date.end)) -> toupdate; # = as.Date('2022/12/31')
    # toupdate$date.end = as.Date('2022/12/31')
    parm.bound_vec = parm.bound_vec %>% mutate(date.end = case_when(parm=='beta' & type == 'SR' & date.end == toupdate$date.end ~ as.Date('2022/12/31'),
                                                                    T ~ date.end))
    parm.bound_vec %>% filter(parm=='beta' & type == 'DA' & period == 'Delta') %>% 
      filter(date.end == max(date.end)) -> toupdate; # = as.Date('2022/12/31')
    # toupdate$date.end = as.Date('2022/12/31')
    parm.bound_vec = parm.bound_vec %>% mutate(date.end = case_when(parm=='beta' & type == 'DA' & date.end == toupdate$date.end ~ as.Date('2022/12/31'),
                                                                    T ~ date.end)) 
    # ifr
    parm.bound_vec %>% filter(parm=='ifr' & type == 'SR' & period == 'Delta') %>% 
      filter(date.end == max(date.end)) -> toupdate; # = as.Date('2022/12/31')
    # toupdate$date.end = as.Date('2022/12/31')
    parm.bound_vec = parm.bound_vec %>% mutate(date.end = case_when(parm== 'ifr' & type == 'SR' & date.end == toupdate$date.end ~ as.Date('2022/12/31'),
                                                                    T ~ date.end)) %>% data.table()
    
    
  } else if (variant.tag == 'Omicron'){
    
    parm.bound_vec = parm.bound_vec %>% dplyr::filter(grepl('omicron', tolower(type)) | 
                                grepl('omicron', tolower(parm)) |
                                grepl('omicron', tolower(variant)) |
                                grepl('omicron', tolower(period)))  %>% data.table()
    
    # adjust the start
    # print('line 463')
    parm.bound_vec[type == 'SR' & period == 'Omicron_BA.1']$date.start = as.Date('2021/10/15')
    # print('line 465')
    parm.bound_vec[type == 'DA' & period == 'Omicron_BA.1']$date.start = as.Date('2021/10/15') %>% data.table()
    
    # if the first week, no cases were detected, set the initial to lower level
    if(da.t[1]$case == 0){
      # none get detected
      alpha_bounds = c(0, .01)
    }
    # if no cases during first week but much more cases 2nd week, increase level of SR for detection rate (alpha)
    if(da.t[1]$case == 0 & da.t[2]$case >10){
      tmp = parm.bound_vec[type == 'SR' & parm == 'alpha' & date.end >= da.t$date[2] & date.start < da.t$date[2]]
      tmp0 = tmp1 = tmp
      tmp1$SR.level = 'Extra3.first1wk'
      tmp1$date.start = da.t[2]$date
      # tmp1$lwr = pmax(tmp1$lwr, .01)
      tmp1$upr = pmax(tmp1$upr, .1)
      tmp0$date.end = da.t[2]$date -1
      parm.bound_vec = rbind(tmp0, tmp1, parm.bound_vec[!(type == 'SR' & parm == 'alpha' & date.end >= da.t$date[2] & date.start < da.t$date[2])])
    }
  }
  
  # time period need filtering restriction on S, E, or I
  tm_rednUpdateS = NULL
  tm_redn = parm.bound_vec[parm == 'S' & type == 'filtering' & period == 'no.imm.escape'] %>% .[,c('date.start','date.end')] 
  tm_redn = tm_redn[complete.cases(tm_redn)]
  if(any(!is.na(tm_redn$date.start))){
    for(i in 1:nrow(tm_redn)){
      tm_rednUpdateS = append(tm_rednUpdateS, seq(as.Date(tm_redn[i]$date.start), as.Date(tm_redn[i]$date.end), by = 'day'))
    }
  }
    
  tm_rednUpdateEI = NULL
  tm_redn = parm.bound_vec[parm == 'EI' & type == 'filtering' & period == 'voc.ini'] %>% .[,c('date.start','date.end')] 
  tm_redn = tm_redn[complete.cases(tm_redn)]
  if(any(!is.na(tm_redn$date.start))){
    for(i in 1:nrow(tm_redn)){
      tm_rednUpdateEI = append(tm_rednUpdateEI, seq(as.Date(tm_redn[i]$date.start), as.Date(tm_redn[i]$date.end), by = 'day'))
    }
  }
  
  tm_reset_cntSR.S = NULL
  tm_redn = parm.bound_vec[parm == 'S' & type == 'filtering' & period == 'reset cntSR.S'] %>% .[,c('date.start','date.end')] 
  tm_redn = tm_redn[complete.cases(tm_redn)]
  if(any(!is.na(tm_redn$date.start))){
    for(i in 1:nrow(tm_redn)){
      tm_reset_cntSR.S = append(tm_reset_cntSR.S, seq(as.Date(tm_redn[i]$date.start), as.Date(tm_redn[i]$date.end), by = 'day'))
    }
  }
  
  obs_i = (da.t$case) %>% as.matrix() 
  obs_vars_i = obs_i
  for(j in 1:num_obs){
    tmp=rep(0,nrow(da.t))
    for (i in 3:nrow(da.t)){
      tmp[i]=mean(obs_i[(i-2):(i-0),j]);
    }
    
    # obs_vars_i[,j] = N/1000 + pmin(tmp * 9, pmax(tmp^2 /50, tmp * 6))
    if(variant.tag == 'non.Omicron'){
      obs_vars_i[,j] = N/1000 + pmin(tmp * 6, pmax(tmp^2 /50, tmp * 4))
    } else {
      # higher infection rate / uncertainty
      obs_vars_i[,j] = N/500 + pmin(tmp * 6, pmax(tmp^2 /50, tmp * 4))
    }
  }
  
  # 8/30/22 - some weeks have very large jumps, e.g. CA 2021/6/27 
  # increase OVE for those?
  p.adj = matrix(1, nrow = nrow(da.t), ncol = num_obs)
  # seach to see if there are such jumps
  for(j in 1:num_obs){
    for (i in 2:(nrow(da.t)-1)){
      if(obs_i[i,j] > obs_i[i-1,j] * 5 & obs_i[i+1, j] < obs_i[i, j] * .9)
        p.adj[i, j] = pmax(2, pmin(obs_i[i,j] / (obs_i[i-1,j] * 2), 3))
    }
  }
  obs_vars_i = obs_vars_i * p.adj
  
  
  # obs_vars_i = ((obs_i + 10) * 4)  %>% as.matrix() 
  # for weeks with very low case rate, increase OVE 
  if(F){
    low.t = quantile(da.t[date > date_ini & date < date_omicron]$case, .1)
    idx.low = which(da.t$case < low.t)
    obs_vars_i[idx.low,] = pmax(2000, obs_vars_i[idx.low,])
  }

  
  obs_d = (da.t$death) %>% as.matrix() 
  obs_vars_d = obs_d
  for(j in 1:num_obs){
    tmp=rep(0,nrow(da.t))
    for (i in 3:nrow(da.t)){
      tmp[i]=mean(obs_d[(i-2):(i-0),j]);
    }
    
    obs_vars_d[,j]= (c(rep(N/1e4,tm_largerVar),rep(N/1e4,nrow(da.t)-tm_largerVar)) + 
                       pmin((tmp^2)/10, tmp*20)
    ) * pOEV;
  }
  

  weeks = da.t$week # for seasonality if applicable
  Week.starts = da.t$date
  
  fcast.wk.starts = as.Date(fcast.start) + seq(0, length.out = nfcast, by = 7)
  weeks.fcast = MMWRweek(fcast.wk.starts)['MMWRweek'] %>% unlist
  
  # if there are fewer than 4 weeks's data, skip
  if(which(Week.starts == fcast.wk.starts[1] - 7) < 4){
    next
  }
  
  # do n runs
  for(ir in dummy_tno){ # tno
    
    print(paste('run', ir))
    
    source(paste0(dir_code,'set_tm2event.R'))
    
    if(variant.tag == 'non.Omicron'){
      So=t(lhs(num_ens,rect = rbind(cbind(.99, 1) * N, # S0
                                    cbind(seed_max/20,seed_max/2), # E0
                                    cbind(seed_max/20,seed_max/2), # I0
                                    cbind(0,seed_max/100) # deaths0
      )))
    } else if(variant.tag == 'Omicron'){
      iniSomicron = c(.5, .9)
      if(Week.starts[1] < as.Date('2021/11/21')){ # f2021/11/21 = irst week of detection in NYC
        # start early, likely low
        # tmp = da.t[date <= as.Date('2021/11/21')]$case
        # seed_max = (weighted.mean(tmp, w = c(.8, .2/length(tmp))) + .1) * 100 # 10000
        if(da.t[1]$case == 0 & da.t[2]$date == as.Date('2021/11/21') & da.t[2]$case > 0){
          # assume ~2 day dubbling time: c1 = c0 * 2^(7/2) -> c0 = c1 / 2^(7/2)
          seed_max = (da.t[2]$case / 2^(7/2) + .1) * 100 # and assume only 1% get detected
        }
        So=t(lhs(num_ens,rect = rbind(iniSomicron * N, # S0
                                      cbind(seed_max/20,seed_max/2), # E0, more uncertainty
                                      cbind(seed_max/20,seed_max/2), # I0
                                      cbind(0,seed_max/100/100) # deaths0
        )))
      } else {
        So=t(lhs(num_ens,rect = rbind(iniSomicron * N, # S0
                                      cbind(seed_max/20,seed_max/2), # E0
                                      cbind(seed_max/20,seed_max/2), # I0
                                      cbind(0,seed_max/100/100) # deaths0
        )))
      }
      
    }
    
    
    
    S0 = So[1:num_gr,,drop=F]
    E0 = So[1:num_gr+num_gr,,drop=F]
    I0 = So[1:num_gr+num_gr*2,,drop=F]
    D0 = So[1:num_gr+num_gr*3,,drop=F]
    
    newItot = I0; newIobs = I0; 
    rownames(S0)=paste0('S',1:num_gr); 
    rownames(E0)=paste0('E',1:num_gr); 
    rownames(I0)=paste0('I',1:num_gr); 
    rownames(D0)=paste0('death',1:num_gr);
    rownames(newItot)=paste0('newItot',1:num_gr); 
    rownames(newIobs)=paste0('newIobs',1:num_gr); 
    
    
    parm.bounds = rbind(beta_bounds, # beta for all loc's
                        c(2,5), # Tei: time from exposed to infectious: incubation time mean = 4
                        c(2,5), # Tir: time from infectous to not (remOEVd)
                        imm_bounds, # immunity period, Trs
                        # c(3,8), # mean Td: reporting delay 
                        # make the range smaller to constrain the model better
                        Td.mean_bounds, # c(5,7), # mean Td: reporting delay
                        Td.sd_bounds, # c(1,3), # Td, sd: reporting delay sd
                        p.mob_bounds, # scaling for mobility
                        alpha_bounds, # reporting rate
                        ifr_bounds # infection fatality risk
    )
    parm.names = c('beta','Tei','Tir','Trs', 'Td.mean', 'Td.sd', 'p.mob','alpha', 'ifr')
    
    rownames(parm.bounds) = parm.names
    parm.bounds
    
    parm0=t(lhs(num_ens,parm.bounds)); rownames(parm0)=rownames(parm.bounds)
    
    
    STATE0=rbind(S0, E0, I0, D0, newIobs, newItot, parm0)
    state.names=rownames(STATE0)
    idx.obs_i= which(state.names == 'newIobs1')  # the random tests are testing the prevalence of infectious - I
    idx.obs_d= which(state.names == 'death1')  # the random tests are testing the prevalence of infectious - I
    idx.newItot = which(state.names == 'newItot1') 
    idx.e = which(state.names == 'E1') 
    idx.i = which(state.names == 'I1') 
    
    num_state = 4 + 2
    
    DAbounds = rbind(matrix(c(rep(0,num_state * num_gr), rep(N,num_state)),num_state * num_gr,2),
                      cbind(parm.bounds[,1]*.5, parm.bounds[,2]*1.5)) # cbind(parm.bounds[,1]*.5, parm.bounds[,2]*1.5)
    rownames(DAbounds)=state.names
    DAbounds[c('E1','I1'),2] = N * .15 # / 20
    DAbounds['death1',2] = N / 100 # 200
    DAbounds['S1',1] = N / 10
    DAbounds[c('newItot1','newIobs1'),2] = N * .2
    
    DAbounds['Td.mean',1] = parm.bounds['Td.mean',1] * .6
    DAbounds['Td.mean',2] = parm.bounds['Td.mean',2] * 1.1
    DAbounds['Trs',1] = ifelse(variant.tag == 'Omicron', 100, 200)
    
    DAbounds['p.mob',] =  parm.bounds['p.mob',]
    
    # SRbounds = cbind(parm.bounds[,1]*.75, parm.bounds[,2]*1.25)
    SRbounds = parm.bounds # cbind(parm.bounds[,1]*.9, parm.bounds[,2]*1.1) # 
    rownames(SRbounds)=parm.names
    SRbounds['Td.mean',] = Td.mean_SRbounds
    
    # fcast.deflat = .9
    tm.ini=1; tmstep=7; newI.previous = NULL; inflat=1.03; state0=STATE0
    fcast.deflat = fcast.deflat_vec[dummy_idf]
    
    severity['death',] = STATE0['ifr',]
    
    # model training and forecast
    {  
      print('running...')
      train.proj = EAKF_rproj(epi.model=epi.model, num_ens=num_ens,inflat=1.03, 
                              fcast.deflat = fcast.deflat, # shrink the distribution during fcast period
                              obs_i=obs_i, obs_vars_i=obs_vars_i, # case
                              obs_d=obs_d, obs_vars_d=obs_vars_d,
                              weeks=weeks,Week.starts=Week.starts,
                              parm.bounds=parm.bounds, DAbounds=DAbounds, SRbounds=SRbounds, 
                              parm.names = rownames(parm.bounds), rel.mob = rel.mob,
                              state0=STATE0, state.names=rownames(STATE0),
                              severity = severity,
                              tm.ini=1, tmstep=7,
                              newI.previous = NULL,
                              parm.bound_vec = parm.bound_vec,
                              # for the projection
                              weeks.fcast, # week of the year to get seasonality 
                              fcast.wk.starts,
                              incl.sce.newV = T,
                              save.cumIetc = F
      )
      
      train.proj$fcast.start.week = fcast.start
      train.proj$fcast.start.month = fcast.start.month
      train.proj$loc = loc.t
      train.proj$fcast.deflat = fcast.deflat
      train.proj$fcast.type = tag.proj
      train.proj$variant = variant.tag
      train.proj$seasonality = sn.tag
      save(train.proj, file = paste0(dir_res, gsub(' ','',loc.t),'_',variant.tag,'_train.proj_',tag.proj,'_df',fcast.deflat,'_',sn.tag,'_',
                                     format(fcast.start,'%Y%m%d'),
                                     '_r',ir,'.RData'))
      
    } 
    
  } # end this run
} # end this location
  
if(F){
  # check outputs
  source(paste0(dir_code,'getPlot.R'))
  # asIs
  da.full.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= as.Date(date.start) & date <= as.Date(tail(fcast.wk.starts,1))]
  train.t = train.proj$states_stats[state == 'newIobs1' & Week.start < fcast.start]
  train.t = merge(train.t, da.t[,c('date', 'case'), with=F] %>% setnames(c('date', 'case'), c('Week.start', 'obs')), 
                  by = 'Week.start')
  proj.t = train.proj$fcast_stats[measure == 'Cases' & scenario == 'asIs']
  proj.t = merge(proj.t, da.full.t %>% setnames(c('date', 'case'), c('Week.start', 'obs')), by = 'Week.start', all.x=T)
  train.t$loc = loc.t
  proj.t$loc = loc.t
  getPlotProj80(tail(train.t, 20), proj.t)
  
  # newV
  da.full.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= as.Date(date.start) & date <= as.Date(tail(fcast.wk.starts,1))]
  train.t = train.proj$states_stats[state == 'newIobs1' & Week.start < fcast.start]
  train.t = merge(train.t, da.t[,c('date', 'case'), with=F] %>% setnames(c('date', 'case'), c('Week.start', 'obs')), 
                  by = 'Week.start')
  proj.t = train.proj$fcast_stats[measure == 'Cases' & scenario == 'newV']
  proj.t = merge(proj.t, da.full.t %>% setnames(c('date', 'case'), c('Week.start', 'obs')), by = 'Week.start', all.x=T)
  train.t$loc = loc.t
  proj.t$loc = loc.t
  getPlotProj80(tail(train.t, 20), proj.t)
  
  da.full.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= as.Date(date.start) & date <= as.Date(tail(fcast.wk.starts,1))]
  train2.t = train.proj$states_stats[state == 'death1' & Week.start < fcast.start]
  train2.t = merge(train2.t, da.t[,c('date', 'death'), with=F] %>% setnames(c('date', 'death'), c('Week.start', 'obs')), 
                  by = 'Week.start')
  proj.t = train.proj$fcast_stats[measure == 'Deaths'  & scenario == 'asIs']
  proj.t = merge(proj.t, da.full.t %>% setnames(c('date', 'death'), c('Week.start', 'obs')), by = 'Week.start', all.x=T)
  train2.t$loc = loc.t
  proj.t$loc = loc.t
  getPlotProj80(tail(train2.t, 20), proj.t)
  
  da.full.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= as.Date(date.start) & date <= as.Date(tail(fcast.wk.starts,1))]
  train2.t = train.proj$states_stats[state == 'death1' & Week.start < fcast.start]
  train2.t = merge(train2.t, da.t[,c('date', 'death'), with=F] %>% setnames(c('date', 'death'), c('Week.start', 'obs')), 
                   by = 'Week.start')
  proj.t = train.proj$fcast_stats[measure == 'Deaths'  & scenario == 'newV']
  proj.t = merge(proj.t, da.full.t %>% setnames(c('date', 'death'), c('Week.start', 'obs')), by = 'Week.start', all.x=T)
  train2.t$loc = loc.t
  proj.t$loc = loc.t
  getPlotProj80(tail(train2.t, 20), proj.t)
  
  
}


