# generate the initial conditions, SR/DA bounds based on timing of variant circulation
# 1/28/22
# 7/28/22 double check variant timing for reprobing S and 
# add SR bound for detection rate during the Omicron wave (much lower later on, per estimates for NYC after Feb 2022, only ~15% detected)
# 8/26/22 - further consider for some states (CA), co-circulation of multiple variants (alpha, epsilon, and gamma starting april 2021)
# some some immune erosion may occur earlier
# 10/31/22 - clearer labeling

tag = '_lowerbeta'

library(tidyverse)
library(data.table)
library(writexl)
library(readxl)

dir_data = './data/'
dir_code = './code/'

# population density data from https://www.statista.com/statistics/183588/population-density-in-the-federal-states-of-the-us/
pop.den = read_xlsx(paste0(dir_data,'statistic_id183588_population-density-in-the-us-by-state-2020.xlsx'), sheet='Data', skip = 5, col_names = F) %>% data.table()
colnames(pop.den) = c('loc','density')
pop.den$density = pop.den$density %>% as.numeric()
pop.den.us = pop.den[loc == 'United States']$density

tm_variants = read.csv(paste0(dir_data, 'tm.variant_by_state_smoothed.wkly.csv')) %>% data.table()
tm_variants$date.start = tm_variants$date.start %>% as.Date
tm_variants$date.end = tm_variants$date.end %>% as.Date

DAT.EPI = read.csv(paste0(dir_data, 'da_case_death_us_nyt.csv')) %>% data.table()
DAT.EPI$date = DAT.EPI$date %>% as.Date
# the list of states 
locs = DAT.EPI$state %>% unique()

da.vac = read.csv(paste0(dir_data, 'da_vx_perM_us_lagged.csv')) %>% data.table()
 
da.variant = read.csv(paste0(dir_data, 'perc.variant.wkly.smoothed_by_state.csv'))  %>% data.table()
tmp = da.variant %>% filter(location == 'Texas')  

DAT.SNraw = read.csv(paste0(dir_data, 'est.raw.sn.2000t2020_us.csv')) %>% data.table()

# transmission rate
beta_bounds_ini = c(.5, .8) # too high?
beta_bounds_ini = c(.5, .8) - .05 # try slightly lower
beta_bounds_Iota = beta_bounds_ini * c(1.2, 1.1)
beta_bounds_Epsilon = beta_bounds_ini * c(1.2, 1.1)
beta_bounds_mixed = beta_bounds_ini * c(1, 1.1)
beta_bounds_Alpha = beta_bounds_ini * c(1.3, 1.5) # c(1.4, 1.5)
beta_bounds_Beta = beta_bounds_ini * c(1.3, 1.4)
beta_bounds_Gamma = beta_bounds_ini * c(1.3, 1.6)
beta_bounds_Delta = beta_bounds_ini * c(1.2, 1.7)  # allow lower bounds b/c mass vaccination
beta_bounds_Omicron_BA.1 = beta_bounds_ini * c(1.5, 2) # allow lower bounds b/c mass vaccination
beta_bounds_Omicron_nonBA.1 = beta_bounds_Omicron_BA.1 * c(1.2, 1.4)
beta_bounds_Omicron_BA.2 = beta_bounds_Omicron_BA.1 * c(1.2, 1.4)
beta_bounds_Omicron_BA.2.12.1 = beta_bounds_Omicron_BA.1  * c(1.2, 1.4)
beta_bounds_Omicron_nonBA.1o2 = beta_bounds_Omicron_BA.1  * c(1.2, 1.4)
beta_bounds = rbind(
  data.table(SR.level = '', period = 'start', beta_bounds_ini %>% t, date.start = '2020/01/01', date.end = '2020/9/30'), # initial detection rate
  data.table(SR.level = 'normal', period = 'wave1', beta_bounds_ini %>% t, date.start = '2020/01/01', date.end = '2020/9/30'), # initial detection rate
  data.table(SR.level = 'normal', period = 'wave2', beta_bounds_mixed %>% t, date.start = '2020/10/01', date.end = '2021/2/28'), # fall/winter wave, depending on location
  data.table(SR.level = 'extra.first2wk', period = 'Alpha', beta_bounds_Alpha %>% t, date.start = '2021/03/01', date.end = '2021/6/30'),# beta wave more severe, less vaccination at the time
  data.table(SR.level = 'normal', period = 'Delta', beta_bounds_Delta %>% t, date.start = '2021/07/01', date.end = '2021/11/30'), # delta wave, more severe
  data.table(SR.level = 'extra.first2wk', period = 'Omicron_BA.1', beta_bounds_Omicron_BA.1 %>% t, date.start = '2021/12/1', date.end = '2022/2/28'), # omicron wave, awareness/holiday pre-gathering testing could increase detection
  # data.table(SR.level = 'normal', period = 'Omicron_nonBA.1', beta_bounds_Omicron_nonBA.1 %>% t, date.start = '2022/3/1', date.end = '2022/12/31')
  data.table(SR.level = 'extra.first2wk', period = 'Omicron_BA.2', beta_bounds_Omicron_BA.2 %>% t, date.start = '2022/2/15', date.end = '2022/3/31'),
  data.table(SR.level = 'extra.first2wk', period = 'Omicron_BA.2.12.1', beta_bounds_Omicron_BA.2 %>% t, date.start = '2022/4/1', date.end = '2022/6/15'),
  data.table(SR.level = 'extra.first2wk', period = 'Omicron_nonBA.1o2', beta_bounds_Omicron_nonBA.1o2 %>% t, date.start = '2022/6/16', date.end = '2022/12/31')
) 
colnames(beta_bounds)[3:4] = c('lwr', 'upr')
beta_bounds$parm = 'beta'
beta_bounds$type = 'SR'
beta_bounds[SR.level == '' & period=='start']$type = 'initialization'
beta_bounds$date.start = beta_bounds$date.start %>% as.Date
beta_bounds$date.end = beta_bounds$date.end %>% as.Date

beta_DAbounds = rbind(
  data.table(SR.level = '', period = 'waves', beta_bounds_ini * c(.2, 1.7) %>% t, date.start = '2020/01/01', date.end = '2021/1/31'), # initial detection rate
  data.table(SR.level = '', period = 'Alpha', beta_bounds_Alpha * c(.2, 1.5) %>% t, date.start = '2021/02/01', date.end = '2021/6/30'),# beta wave more severe, less vaccination at the time
  data.table(SR.level = '', period = 'Delta', beta_bounds_Delta * c(.2, 1.5) %>% t, date.start = '2021/07/01', date.end = '2021/11/30'), # delta wave, more severe
  data.table(SR.level = '', period = 'Omicron_BA.1', beta_bounds_Omicron_BA.1 * c(.2, 1.5) %>% t, date.start = '2021/12/1', date.end = '2022/2/28'), # omicron wave, awareness/holiday pre-gathering testing could increase detection
  data.table(SR.level = '', period = 'Omicron_nonBA.1o2', beta_bounds_Omicron_nonBA.1o2 * c(.2, 1.5) %>% t, date.start = '2022/3/1', date.end = '2022/12/31') # subvariants
) 
colnames(beta_DAbounds)[3:4] = c('lwr', 'upr')
beta_DAbounds$parm = 'beta'
beta_DAbounds$type = 'DA'
beta_DAbounds$date.start = beta_DAbounds$date.start %>% as.Date
beta_DAbounds$date.end = beta_DAbounds$date.end %>% as.Date

# detection rate
alpha_bounds_ini = c(.01, .05) # initial detection rate
alpha_bounds_w1 = c(.08, .35) # gradually increases after a few weeks as testing ramped up
# w1 and w2 are given similar detection bounds, b/c even though in general detection rate in w1 is likely lower, 
# the infection demographics leaned toward older adults, which could lead to overall higher detection rate due to more severe disease
alpha_bounds_2020summer = c(.05, .25) # summer lower detection in general
alpha_bounds_w2 = c(.12, .4) # fall/winter wave, depending on location
alpha_bounds_Alpha = c(.2, .4) # alpha wave more severe, less vaccination at the time
alpha_bounds_2021vx = c(.05, .3) # post wave 3, more ppl getting vaccinated
alpha_bounds_2021summer = c(.03, .25) # summer, prior to delta, more ppl getting vaccinated
alpha_bounds_Delta = c(.1, .5) # delta wave, more severe
alpha_bounds_Omicron0 = c(0.001, .05) # initial very low detection
alpha_bounds_Omicron1a = c(.1, .5) # omicron wave, make it higher b/c it is massively replaced, awareness/holiday pre-gathering testing could increase detection
alpha_bounds_Omicron1b = c(.1, .6) # omicron wave, make it higher b/c it is massively replaced, awareness/holiday pre-gathering testing could increase detection
alpha_bounds_Omicron2 = c(.05, .2) # omicron wave, after the first couple weeks, test capacity max-out, milder disease -> lower detection
alpha_bounds_Omicron3 = c(.02, .18) # omicron wave, much lower detection rate later on
alpha_bounds = rbind(
  data.table(SR.level = '', period = 'start', alpha_bounds_ini %>% t, date.start = '2020/01/01', date.end = '2020/3/7'), # initial detection rate
  data.table(SR.level = 'normal', period = 'start', alpha_bounds_ini %>% t, date.start = '2020/01/01', date.end = '2020/3/7'), # initial detection rate
  data.table(SR.level = 'EXTRA.first1wk', period = 'wave1', alpha_bounds_w1 %>% t, date.start = '2020/03/08', date.end = '2020/5/31'), # gradually increases after a few weeks as testing rampped up
  data.table(SR.level = 'extra.first2wk', period = '2020summer', alpha_bounds_2020summer %>% t, date.start = '2020/06/01', date.end = '2020/08/31'),# summer lower detection in general
  data.table(SR.level = 'extra.first6wk', period = 'wave2', alpha_bounds_w2 %>% t, date.start = '2020/09/01', date.end = '2021/2/28'), # fall/winter wave, depending on location
  data.table(SR.level = 'extra.first2wk', period = 'Alpha', alpha_bounds_Alpha %>% t, date.start = '2021/03/01', date.end = '2021/4/30'),# alpha wave more severe, less vaccination at the time
  data.table(SR.level = 'extra.first2wk', period = 'massvax', alpha_bounds_2021vx %>% t, date.start = '2021/05/01', date.end = '2021/5/31'), # post wave 3, more ppl getting vaccinated
  data.table(SR.level = 'Extra.first1wk', period = '2021summer', alpha_bounds_2021summer %>% t, date.start = '2021/06/01', date.end = '2021/6/30'),# summer, prior to delta, more ppl getting vaccinated 
  data.table(SR.level = 'normal', period = 'Delta', alpha_bounds_Delta %>% t, date.start = '2021/07/01', date.end = '2021/11/30'), # delta wave, more severe
  data.table(SR.level = 'extra.first2wk', period = 'Omicron0', alpha_bounds_Omicron0 %>% t, date.start = '2021/11/1', date.end = '2021/11/21'),# omicron wave, awareness/holiday pre-gathering testing could increase detection 
  data.table(SR.level = 'Extra3.first1wk', period = 'Omicron1a', alpha_bounds_Omicron1a %>% t, date.start = '2021/11/22', date.end = '2021/11/28'),# omicron wave, awareness/holiday pre-gathering testing could increase detection
  data.table(SR.level = 'EXTRA.first1wk', period = 'Omicron1b', alpha_bounds_Omicron1b %>% t, date.start = '2021/11/29', date.end = '2021/12/25'),# omicron wave, awareness/holiday pre-gathering testing could increase detection
  data.table(SR.level = 'Extra3.first1wk', period = 'Omicron2', alpha_bounds_Omicron2 %>% t, date.start = '2021/12/26', date.end = '2022/2/28'), # omicron wave, after the first couple weeks, test capacity max-out, milder disease -> lower detection
  data.table(SR.level = 'extra.first4wk', period = 'Omicron3', alpha_bounds_Omicron3 %>% t, date.start = '2022/3/1', date.end = '2022/12/31') # subvariant waves
) 
colnames(alpha_bounds)[3:4] = c('lwr', 'upr')
alpha_bounds$parm = 'alpha'
alpha_bounds$type = 'SR'
alpha_bounds[SR.level == '' & period=='start']$type = 'initialization'
alpha_bounds$date.start = alpha_bounds$date.start %>% as.Date
alpha_bounds$date.end = alpha_bounds$date.end %>% as.Date


alpha_DAbounds = rbind(
  data.table(SR.level = '', period = 'start', c(.01,.3) %>% t, date.start = '2020/01/01', date.end = '2020/3/7'), # initial detection rate
  data.table(SR.level = '', period = 'waves', c(.05,.6) %>% t, date.start = '2020/03/08', date.end = '2020/5/15'), # gradually increases after a few weeks as testing rampped up
  data.table(SR.level = '', period = 'massvax', c(.03,.6) %>% t, date.start = '2020/05/16', date.end = '2022/12/31'), # mass vacc / variants
  data.table(SR.level = '', period = 'omicron0', c(0,.2) %>% t, date.start = '2021/11/1', date.end = '2021/11/21'), # initial phase of Omicron wave
  data.table(SR.level = '', period = 'omicron1', c(.01,.6) %>% t, date.start = '2021/11/22', date.end = '2022/2/28'),
  data.table(SR.level = '', period = 'omicron2', c(.01,.4) %>% t, date.start = '2021/3/1', date.end = '2022/12/31')
  ) 
colnames(alpha_DAbounds)[3:4] = c('lwr', 'upr')
alpha_DAbounds$parm = 'alpha'
alpha_DAbounds$type = 'DA'
alpha_DAbounds$date.start = alpha_DAbounds$date.start %>% as.Date
alpha_DAbounds$date.end = alpha_DAbounds$date.end %>% as.Date

# infection-fatality risk
ifr_bounds_ini = c(.005, .015)
ifr_bounds_w1early = c(.005, .025)  # higher early on for places with high infection rate in nursing homes
ifr_bounds_w1late = c(.001, .015) # lower later
ifr_bounds_2020summer = c(.0001, .005) # summer, more young ages infected, lower overall ifr
ifr_bounds_w2 = c(.0002, .01) # winter wave, increased infections in older ages
ifr_bounds_Alpha = c(.0001, .015) # alpha wave, more severe, but also more vaccinated
ifr_bounds_Delta = c(.0001, .015)  # delta wave more severe, but also more vaccinated?
ifr_bounds_vx = c(.0001, .01) # delta wave, later phase
ifr_bounds_Omicron1 = c(.00004, .004) # ~40% lower
ifr_bounds_Omicron2 = c(.00004, .004) * c(.2, .8) # later phase (~April/May per NYC est), with reinfection/boosters/etc, further decline
ifr_bounds = rbind(
  data.table(SR.level = '', period = 'start', ifr_bounds_ini %>% t, date.start = '2020/01/01', date.end = '2020/3/15'), # initial detection rate
  data.table(SR.level = 'normal', period = 'start', ifr_bounds_ini %>% t, date.start = '2020/01/01', date.end = '2020/3/15'), # initial detection rate
  data.table(SR.level = 'extra.first2wk', period = 'wave1early', ifr_bounds_w1early %>% t, date.start = '2020/03/16', date.end = '2020/4/15'), # gradually increases after a few weeks as testing rampped up
  data.table(SR.level = 'extra.first2wk', period = 'wave1late', ifr_bounds_w1late %>% t, date.start = '2020/04/16', date.end = '2020/5/31'), # gradually increases after a few weeks as testing rampped up
  data.table(SR.level = 'extra.first2wk', period = '2020summer', ifr_bounds_2020summer %>% t, date.start = '2020/06/01', date.end = '2020/9/30'),# summer lower detection in general
  data.table(SR.level = 'extra.first2wk', period = 'wave2', ifr_bounds_w2 %>% t, date.start = '2020/10/01', date.end = '2021/2/28'), # fall/winter wave, depending on location
  data.table(SR.level = 'normal', period = 'Alpha', ifr_bounds_Alpha %>% t, date.start = '2021/03/01', date.end = '2021/4/30'),# ifr wave more severe, less vaccination at the time
  data.table(SR.level = 'extra.first2wk', period = 'massvax', ifr_bounds_vx %>% t, date.start = '2021/05/01', date.end = '2021/7/31'), # post wave 3, more ppl getting vaccinated
  data.table(SR.level = 'normal', period = 'Delta', ifr_bounds_Delta %>% t, date.start = '2021/08/01', date.end = '2021/11/30'), # delta wave, more severe
  data.table(SR.level = 'extra.first2wk', period = 'Omicron_BA.1', ifr_bounds_Omicron1 %>% t, date.start = '2021/12/1', date.end = '2022/3/31'), # omicron wave, after the first couple weeks, test capacity max-out, milder disease -> lower detection
  data.table(SR.level = 'extra.first2wk', period = 'Omicron_nonBA.1', ifr_bounds_Omicron2 %>% t, date.start = '2022/4/1', date.end = '2022/12/31') # after the initial BA.1 wave
) 
colnames(ifr_bounds)[3:4] = c('lwr', 'upr')
ifr_bounds$parm = 'ifr'
ifr_bounds$type = 'SR'
ifr_bounds[SR.level == '' & period=='start']$type = 'initialization'
ifr_bounds$date.start = ifr_bounds$date.start %>% as.Date
ifr_bounds$date.end = ifr_bounds$date.end %>% as.Date

ifr_DAbounds = rbind(
  data.table(SR.level = '', period = '', c(1e-5, .05) %>% t, date.start = '2020/01/01', date.end = '2022/12/31'), # initial detection rate
  data.table(SR.level = '', period = 'Omicron', c(1e-6, .05) %>% t, date.start = '2020/01/01', date.end = '2022/12/31') # initial detection rate
) 
colnames(ifr_DAbounds)[3:4] = c('lwr', 'upr')
ifr_DAbounds$parm = 'ifr'
ifr_DAbounds$type = 'DA'
ifr_DAbounds$date.start = ifr_DAbounds$date.start %>% as.Date
ifr_DAbounds$date.end = ifr_DAbounds$date.end %>% as.Date


Tei_bounds = rbind(
  ini = c(2,5), 
  Omicron = c(1.5, 4)
)


Tei_bounds_ini = c(2,5)
Tir_bounds_ini = c(2,5)
p.mob_bounds_ini = c(.5, 1.5); # scaling for mobility
Td.mean_bounds_ini =  c(5,8) # mean Td: reporting delay
Td.mean_SRbounds_ini = c(5, 7)
Td.sd_bounds_ini = c(1,3) # Td, sd: reporting delay sd 
imm_bounds_ini = c(2, 3) * 365

parm.bounds = rbind(beta_bounds_ini, # beta for all loc's
                    Tei_bounds_ini, # Tei: time from exposed to infectious: incubation time mean = 4
                    Tir_bounds_ini, # Tir: time from infectous to not (remOEVd)
                    imm_bounds_ini, # immunity period, Trs
                    # c(3,8), # mean Td: reporting delay 
                    # make the range smaller to constrain the model better
                    Td.mean_bounds_ini, # c(5,7), # mean Td: reporting delay
                    Td.sd_bounds_ini, # c(1,3), # Td, sd: reporting delay sd
                    p.mob_bounds_ini, # scaling for mobility
                    alpha_bounds_ini, # reporting rate
                    ifr_bounds_ini # infection fatality risk
)
date_ini = as.Date('2020/4/1')
date_omicron = as.Date('2021/11/15')
loc.t = 'New York' # 'Florida' # 'California'
da.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
low.t = quantile(da.t[date > date_ini & date < date_omicron]$case, .05)
idx.low = which(da.t$case < low.t)
da.t$date[idx.low]

loc.t = 'New York'
loc.t = 'California'
loc.t = 'Texas'
# recode weeks with low cases and need larger OVE
BOUNDS = NULL
for(loc.t in locs){
  
  if(T){
    da.t = DAT.EPI[state == loc.t] %>% dcast(., date + year + week ~ data.type, value.var = 'value')
    low.t = quantile(da.t[date > date_ini & date < date_omicron]$case, .05)
    vlow.t = quantile(da.t[date > date_ini & date < date_omicron]$case, .025)
    high.t = quantile(da.t[date > date_ini & date < date_omicron]$case, .8)
    idx.low = which(da.t$case < low.t)
    idx.vlow = which(da.t$case < vlow.t)
    idx.high = which(da.t$case > high.t)
    dates.low.case = da.t$date[idx.low] %>% as.Date # %>% sQuote() %>% paste(.,collapse = ',')
    dates.vlow.case = da.t$date[idx.vlow] %>% as.Date
    dates.high.case = da.t$date[idx.high]  %>% as.Date
    
    
    dates.incr = da.t$date[which(da.t$case[-1] / da.t$case[-nrow(da.t)] > 1.5)+1] %>% as.Date # for detection, so no shift
    # for beta, so shift 1 by 1 week
    dates.vincr = da.t$date[which(da.t$case[-1] / da.t$case[-nrow(da.t)] > 2)] %>% as.Date  # very large increase
    
    
    # pool them together to get the number of waves?
    dates.t = append(dates.low.case, dates.high.case) %>% sort
    names(dates.t) = rep('low',length(dates.t))
    names(dates.t)[dates.t %in% dates.high.case] = 'high'
    
    # def a wave as at least 3 consecutive weeks with case > the 25%th percentile?
    wcut.t = quantile(da.t[date > date_ini & date < date_omicron]$case, .25)
    idx.t = which(da.t$case > wcut.t)
    wdates.t = da.t$date[idx.t]  %>% as.Date
    wdates.t
  }
  
  
  bounds.t = NULL
  
  beta_bounds.t = beta_bounds
  
  # adjust beta based on population density?
  pden.t = pop.den[loc == loc.t]$density
  if(pden.t < pop.den.us / 2.5){
    
    
    # for Wyoming - ~ .7 
    if(tag == '_lowerbeta'){
      L.upr = .25; L.lwr = .1; k = .5
    } else {
      L.upr = .35; L.lwr = .2; k = .5
    }

    # (1 - L.lwr/(1+exp(-k * (log(pop.den.us)-log(pop.den[density < pop.den.us /1.5]$density)))))
    # (1 - L.upr/(1+exp(-k * (log(pop.den.us)-log(pop.den[density < pop.den.us /1.5]$density)))))
    # use a logistic function with a mean of .6?
    beta_bounds.t$lwr = (beta_bounds.t$lwr * (1 - L.lwr/(1+exp(-k * (log(pop.den.us) - log(pden.t)))))) %>% round(2)
    beta_bounds.t$upr = (beta_bounds.t$upr * (1 - L.upr/(1+exp(-k * (log(pop.den.us) - log(pden.t)))))) %>% round(2)
    
    # only adjust the upper bound a bit, otherwise too much uncertainty
    if(F){ # not so good
      if(pden.t < pop.den.us / 10){
        
        beta_bounds.t$lwr = beta_bounds.t$lwr * .9
        beta_bounds.t$upr = beta_bounds.t$upr * .9
      } else {
        # beta_bounds.t$lwr = beta_bounds.t$lwr * .9
        beta_bounds.t$upr = beta_bounds.t$upr * .9
      }
    }
  } else if (pden.t < pop.den.us / 1.2) {
    # for Wyoming - ~ .7 
    if(tag == '_lowerbeta'){
      L.upr = 0;  L.lwr = .05; k = .5
    } else {
      L.upr = .1;  L.lwr = .2; k = .5
    }
    
    # (1 - L.lwr/(1+exp(-k * (log(pop.den.us)-log(pop.den[density < pop.den.us /1.2]$density)))))
    # (1 - L.upr/(1+exp(-k * (log(pop.den.us)-log(pop.den[density < pop.den.us /1.2]$density)))))
    # use a logistic function with a mean of .6?
    beta_bounds.t$lwr = (beta_bounds.t$lwr * (1 - L.lwr/(1+exp(-k * (log(pop.den.us) - log(pden.t)))))) %>% round(2)
    beta_bounds.t$upr = (beta_bounds.t$upr * (1 - L.upr/(1+exp(-k * (log(pop.den.us) - log(pden.t)))))) %>% round(2)
    
  }
  
  alpha_bounds.t = alpha_bounds
  ifr_bounds.t = ifr_bounds
  alpha_DAbounds.t = alpha_DAbounds
  beta_DAbounds.t = beta_DAbounds
  ifr_DAbounds.t = ifr_DAbounds
  # based on the variants and cross-check with case data?
  vda.t = tm_variants[location == loc.t] 
  # adjust timing based on local variant circulation 
  # break it into type of parameters, as some are based on timing, some are based on variant
  for(bbt in c('alpha_bounds.t', 'alpha_DAbounds.t')){  # 'beta_bounds.t','ifr_bounds.t', 'beta_DAbounds.t','ifr_DAbounds.t'
    bb = get(bbt)
    
    # update start of mass vx based on vax data
    da.vac.t = da.vac %>% filter(state == loc.t)
    da.vac.t$cum.v2 = cumsum(da.vac.t$n.v2) / 1e6 * 100
    massvax.date.t = da.vac.t %>% filter(cum.v2 > 40) %>% head(1) %>% .$date %>% as.Date
    if(bbt == 'alpha_bounds.t'){
      bb[period == 'massvax']$date.start = massvax.date.t
      if(bb[period == 'massvax']$date.end < massvax.date.t){
        # if it coincide with the summer, don't need it
        bb = bb %>% filter(period != 'massvax')
      }
    } # only do this for the SR bounds
    
    
    sn.t = DAT.SNraw %>% filter(state == loc.t & week %in% 20:52) 
    if(nrow(sn.t)>0){
      # adjust for wave2, when det rate started to increase
      sn.t$date = MMWRweek::MMWRweek2Date(MMWRyear = rep(2020, nrow(sn.t)), MMWRweek = sn.t$week, MMWRday = 1)
      sn.t = sn.t %>% filter(date >= alpha_bounds.t[period == 'wave2']$date.start)
      dates.w2 = wdates.t[wdates.t >= as.Date('2020/9/1') & wdates.t <= as.Date('2021/02/28')]
      alpha_bounds.t[period == 'wave2']$date.start = min(sn.t[value > 1] %>% head(1)  %>% .$date,
                                                         dates.w2[1])
      alpha_bounds.t[period == '2020summer']$date.end = alpha_bounds.t[period == 'wave2']$date.start -1 
      
    }
    # check if there is a summer wave, if so, allow higher det rate
    nwk.high.2020summer = length(wdates.t[wdates.t >= alpha_bounds.t[period == '2020summer']$date.start & wdates.t <= alpha_bounds.t[period == '2020summer']$date.end])
    flag.2020summer.wave = nwk.high.2020summer >= 4
    if(flag.2020summer.wave){
      alpha_bounds.t[period == '2020summer']$upr = alpha_bounds.t[period == '2020summer']$upr * ifelse(nwk.high.2020summer>=8, 1.5, 1.25)
    }
    
    for(ii in 1:nrow(bb)){
      for(v.t in c('Alpha','Delta', 'Omicron0','Omicron1a', 'Omicron1b', 'Omicron1', 'Omicron2', 'Omicron3')){  # , 'Omicron_BA.1','Omicron_nonBA.1'
        if(bb[ii]$period== v.t){ # gsub('1','',v.t)
          
          # vt = vda.t[which(v.main==gsub('1','',v.t))]
          # use the 10% mark for detection rate
          if(v.t %in% c('Omicron0', 'Omicron1', 'Omicron1a')){ # ,'Omicron2'
            vt = vda.t %>% filter(type == '10perc' & 
                                    variant == 'Omicron_BA.1' # use the timing of BA.1
                                  # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
            )
          } else if(v.t %in% c('Omicron1b', 'Omicron2')){ # ,'Omicron2'
            vt = vda.t %>% filter(type == '25perc' & 
                                    variant == 'Omicron_BA.1' # use the timing of BA.1
                                  # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
            )
          } else if(v.t %in% c('Omicron3')){
            vt = vda.t %>% filter(type == '25perc' & 
                                    variant == 'Omicron_nonBA.1o2' # use the timing of BA.1
                                  # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
            )
          } else {
            vt = vda.t %>% filter(type == '25perc' & 
                                    variant == v.t
                                  # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
            )
          }
          
          if(nrow(vt) < 1) 
            next
          
          if(!grepl('Omicron', v.t)){
            bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
            bb[ii]$date.end = (vt$date.end %>% as.Date) # %>% format('%Y-%m-%d')
          }
          
          # also modified adjacent rows
          if(v.t == 'Alpha'){ # it can overlap with summer/mass vacc
            if(bbt %in% c('alpha_bounds.t', 'ifr_bounds.t')){
              if(nrow(bb[period == 'massvax']) > 0){ # maxx vax period exists
                if(bb[ii]$date.end >= bb[period == 'massvax']$date.start)
                  bb[ii]$date.end = (bb[period == 'massvax']$date.start %>% as.Date) - 1 # use the massvx instead
              } else if(nrow(bb[period == '2021summer']) > 0){ # if not, use summer instead
                if(bb[ii]$date.end >= bb[period == '2021summer']$date.start)
                  bb[ii]$date.end = (bb[period == '2021summer']$date.start %>% as.Date) - 1 # use the 2021summer instead
              }
              
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((bb[ii]$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
              
            } 
          } else if(v.t == 'Delta' & bbt %in% c('alpha_bounds.t') ){ # it can overlap with 2021summer
            if(bb[ii]$date.start <= bb[period == '2021summer']$date.end # & bb[ii]$date.start %in% dates.low.case
               ){
              
              tmp = seq(bb[ii]$date.start, bb[period == '2021summer']$date.end, by = 'day')
              tmp = tmp[weekdays(tmp)=='Sunday']
              date.sum.but.high = tmp[which(!tmp %in% dates.low.case)] # it is summer but case is not low
              if(length(date.sum.but.high) > 0){
                bb[ii]$date.start = min((bb[period == '2021summer']$date.end %>% as.Date) + 1,
                                         (date.sum.but.high %>% as.Date) - 6) # use the 2021summer dates instead
                # also need to update end of summer 
                bb[ii - 1]$date.end = bb[ii]$date.start - 1
              } else {
                bb[ii]$date.start = (bb[period == '2021summer']$date.end %>% as.Date) + 1 # use the 2021summer dates instead
                
              }
              
              # bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((bb[ii]$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
            } else {
              
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((bb[ii]$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
            }
            
          } else if(grepl('Omicron', v.t) & bbt %in% c('alpha_bounds.t')) {
            # omicron wave, change by stage (initial detection, lower afterwards)
            # adjust the earlier bounds
            # v.t == 'Omicron0' # only adjust the end as the start of the initial detection, and set the start to 11/1/21
            
            if(v.t == 'Omicron0'){ # initial phase
              # only adjust the end as the start of the initial detection, and set the start to 11/1/21
              bb[ii]$date.start = as.Date('2021/10/15') # make it earlier as it doesn't matter, as long as it is earlier than the initial emergence
              bb[ii]$date.end = (vt$date.start %>% as.Date) -1 # %>% format('%Y-%m-%d') # need this b/c the start may be adjust when checking for delta
            } else if(v.t == 'Omicron1'){ # initial phase
              # only adjust the start
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
            } else if(v.t == 'Omicron1a'){ # initial phase
              # only adjust the start
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
            } else if(v.t == 'Omicron1b'){ # initial phase
              # only adjust the start
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
            } else if (v.t =='Omicron2'){ # BA.1 wave later stage
              # only adjust the end
              bb[ii]$date.end = (vt$date.end %>% as.Date) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((vt$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
            } else if (v.t =='Omicron3'){ # higher transmission due to immune erosive BA.4/BA.5
              # only adjust the start
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
            }
            
          } else {
            bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
            bb[ii+1]$date.start = ((vt$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
          }
          
          # check the start and end
          # in case a wave started late
          if(bb[ii]$date.start < bb[ii-1]$date.end & 
             bb[ii]$period != 'Omicron0' # seperate non.omicron and omicron
             )
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1) # 
          if(bb[ii]$date.start != bb[ii-1]$date.end +1 & 
             bb[ii]$period != 'Omicron0' # seperate non.omicron and omicron
             )
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1)
        }
      }
    }
    
    # check if 2021summer is needed at all
    if(bbt == 'alpha_bounds.t') {
      if(bb[period=='2021summer']$date.end < bb[period=='2021summer']$date.start){
        bb = bb %>% filter(period!='2021summer')
        # then need to re-match the dates
        for(ii in 2:nrow(bb)){
          if(bb[ii]$date.start < bb[ii-1]$date.end & bb[ii]$period != 'Omicron0') # seperate non.omicron and omicron)
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1) # 
          if(bb[ii]$date.start != bb[ii-1]$date.end +1 & bb[ii]$period != 'Omicron0')
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1)
        }
      }
    } 
    
    # wave 2 - was it due to a more severe variant (e.g. epsilon and iota and some alpha)? if so, may have higher detection rate
    da.variant.t = da.variant %>% filter(location == loc.t & as.Date(week) %in% seq(bb[period %in% c('wave2','waves')]$date.start, bb[period %in% c('wave2','waves')]$date.end, by = 'day'))
    if(any(!is.na(da.variant.t[,c('Epsilon','Iota','Alpha')])))
      if(max(da.variant.t$Epsilon + da.variant.t$Iota + da.variant.t$Alpha) > .4){
        tmp.p = ifelse(max(c(da.variant.t$Epsilo, da.variant.t$Iota, da.variant.t$Alpha)) > .6, 1.35, 1.2)
        bb[period %in% c('wave2','waves')]$lwr = bb[period %in% c('wave2','waves')]$lwr * 1.2
        bb[period %in% c('wave2','waves')]$upr = bb[period %in% c('wave2','waves')]$upr * tmp.p # 1.2
        # also increase the det rate for alpha wave a bit, due to likely continued higher awareness
        bb[period == 'Alpha']$upr = bb[period == 'Alpha']$upr * 1.2
      }
      
      
    
    # update it
    bb$location = loc.t
    eval(parse(text = paste(bbt, '= bb')))
  } # this type of bound
  
  for(bbt in c('beta_bounds.t','beta_DAbounds.t')){  # 'ifr_bounds.t', ,'ifr_DAbounds.t'
    bb = get(bbt)
    for(ii in 1:nrow(bb)){
      for(v.t in c('Alpha','Delta','Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.2.12.1','Omicron_nonBA.1o2')){  # , 
        if(bb[ii]$period== v.t){ # gsub('1','',v.t)
          
          # vt = vda.t[which(v.main==gsub('1','',v.t))]
          # use the 25% mark for detection rate
          vt = vda.t %>% filter(type == '25perc' & 
                                  variant == v.t
                                # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
          )
          
          if(nrow(vt) < 1) 
            next
          
          bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
          bb[ii]$date.end = (vt$date.end %>% as.Date) # %>% format('%Y-%m-%d')
          
          {
            bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
            bb[ii+1]$date.start = ((vt$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
          }
          
          # check the start and end
          # in case a wave started late
          if(bb[ii]$date.start < bb[ii-1]$date.end)
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1) # 
          if(bb[ii]$date.start != bb[ii-1]$date.end +1)
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1)
          
        }
      }
    }
    
    
    # update it
    bb$location = loc.t
    eval(parse(text = paste(bbt, '= bb')))
  } # this type of bound
  
  for(bbt in c('ifr_bounds.t', 'ifr_DAbounds.t')){  # 
    bb = get(bbt)
    
    # update start of mass vx based on vax data
    if(bbt == 'ifr_bounds.t'){
      da.vac.t = da.vac %>% filter(state == loc.t)
      da.vac.t$cum.v2 = cumsum(da.vac.t$n.v2) / 1e6 * 100
      massvax.date.t = da.vac.t %>% filter(cum.v2 > 40) %>% head(1) %>% .$date %>% as.Date
      bb[period == 'massvax']$date.start = massvax.date.t
      if(bb[period == 'massvax']$date.end < massvax.date.t){
        # if it coincide with the summer, don't need it
        bb = bb %>% filter(period != 'massvax')
      }
    }
    
    
    for(ii in 1:nrow(bb)){
      for(v.t in c('Alpha','Delta','Omicron_BA.1', 'Omicron_nonBA.1')){  # , 
        if(bb[ii]$period== v.t){ # gsub('1','',v.t)
          
          # vt = vda.t[which(v.main==gsub('1','',v.t))]
          # use the 50% mark for IFR
          vt = vda.t %>% filter(type == '50perc' & 
                                  variant == v.t
                                # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
          )
          if(v.t == 'Omicron_nonBA.1'){ # use Omicron_nonBA.1o2 instead
            vt = vda.t %>% filter(type == '25perc' & 
                                    variant == 'Omicron_nonBA.1o2'
                                  # variant == gsub(pattern='1|2|3|_BA.1|_nonBA.1',replacement='', x = v.t, perl=T)
            )
          }
          
          if(nrow(vt) < 1) 
            next
          
          bb[ii]$date.start = (vt$date.start %>% as.Date) # %>% format('%Y-%m-%d')
          bb[ii]$date.end = (vt$date.end %>% as.Date) # %>% format('%Y-%m-%d')
          
          # also modified adjacent rows
          if(v.t == 'Alpha'){ # it can overlap with summer/mass vacc
            if(bbt %in% c('alpha_bounds.t', 'ifr_bounds.t')){
              if(nrow(bb[period == 'massvax']) > 0){ # maxx vax period exists
                if(bb[ii]$date.end >= bb[period == 'massvax']$date.start)
                  bb[ii]$date.end = (bb[period == 'massvax']$date.start %>% as.Date) - 1 # use the massvx instead
              } else if(nrow(bb[period == '2021summer']) > 0){ # if not, use summer instead
                if(bb[ii]$date.end >= bb[period == '2021summer']$date.start)
                  bb[ii]$date.end = (bb[period == '2021summer']$date.start %>% as.Date) - 1 # use the 2021summer instead
              }
              
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((bb[ii]$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
              
            } 
          } else if(v.t == 'Delta' & bbt %in% c('alpha_bounds.t')){ # it can overlap with 2021summer
            if(bb[ii]$date.start <= bb[period == '2021summer']$date.end # & bb[ii]$date.start %in% dates.low.case
            ){
              bb[ii]$date.start = (bb[period == '2021summer']$date.end %>% as.Date) + 1 # use the 2021summer dates instead
              # bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((bb[ii]$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
            } else {
              bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
              bb[ii+1]$date.start = ((bb[ii]$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
            }
            
          } else {
            bb[ii-1]$date.end = ((vt$date.start %>% as.Date) - 1) # %>% format('%Y-%m-%d')
            bb[ii+1]$date.start = ((vt$date.end %>% as.Date) + 1) # %>% format('%Y-%m-%d')
          }
          
          # check the start and end
          # in case a wave started late
          if(bb[ii]$date.start < bb[ii-1]$date.end)
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1) # 
          if(bb[ii]$date.start != bb[ii-1]$date.end +1)
            bb[ii-1]$date.end = ((bb[ii]$date.start %>% as.Date) - 1)
        }
      }
    }
    
    # wave 2 - was it due to a more severe variant (e.g. epsilon and iota and some alpha)? if so, may have higher detection rate
    if(bbt == 'ifr_bounds.t'){
      da.variant.t = da.variant %>% filter(location == loc.t & as.Date(week) %in% seq(bb[period == 'wave2']$date.start, bb[period == 'wave2']$date.end, by = 'day'))
      if(any(!is.na(da.variant.t[,c('Epsilon','Iota','Alpha')])))
        if(max(da.variant.t$Epsilon + da.variant.t$Iota + da.variant.t$Alpha) > .4){
          # bb[period == 'wave2']$lwr = bb[period == 'wave2']$lwr * 1.2
          bb[period == 'wave2']$upr = bb[period == 'wave2']$upr * 1.25
        }
    }
    
    
    # update it
    bb$location = loc.t
    eval(parse(text = paste(bbt, '= bb')))
  } # this type of bound
  
  # for detection rate, check if cases were very low during 2021 summer prior to the delta surge
  if(nrow(alpha_bounds.t[period == '2021summer']) > 0)
    if(sum(seq(as.Date(alpha_bounds.t[period == '2021summer']$date.start), as.Date(alpha_bounds.t[period == '2021summer']$date.end), by = 'day') %in% dates.low.case) >=3 & 
       sum(seq(as.Date(alpha_bounds.t[period == '2021summer']$date.start), as.Date(alpha_bounds.t[period == '2021summer']$date.end), by = 'day') %in% dates.vlow.case) >=1){
      # the time to the first week in dates.vlow.case
      d.t = dates.vlow.case[which((dates.vlow.case %in% seq(as.Date(alpha_bounds.t[period == '2021summer']$date.start), as.Date(alpha_bounds.t[period == '2021summer']$date.end), by = 'day'))==T)][1]
      tmp0 = tmp1 = alpha_bounds.t[period == '2021summer']
      tmp0$date.end = as.Date(d.t) - 1
      tmp0$SR.level = 'Extra.first1wk'
      tmp1$date.start = as.Date(d.t)
      tmp1$lwr = .02; tmp1$upr = .1
      tmp1$SR.level = 'Extra2.first1wk'
      alpha_bounds.t = rbind(alpha_bounds.t[period != '2021summer'], tmp0, tmp1)
      
      # if this is used, also need special adj to increase it
      d.t = dates.incr[which((dates.incr %in% seq(as.Date(alpha_bounds.t[period == 'Delta']$date.start), as.Date(alpha_bounds.t[period == 'Delta']$date.end), by = 'day'))==T)][1]
      ds.t = dates.incr[which((dates.incr %in% seq(as.Date(alpha_bounds.t[period == 'Delta']$date.start), as.Date(alpha_bounds.t[period == 'Delta']$date.end), by = 'day'))==T)]
      if(as.Date(tail(ds.t,1)) <= as.Date(d.t) + 3*7){
        d2.t = as.Date(tail(ds.t,1))
      } else {
        d2.t = as.Date(d.t) + 3*7
      }
      tmp0 = tmp1 = tmp2 = alpha_bounds.t[period == 'Delta']
      tmp0$date.end = as.Date(d.t) - 1
      tmp1$date.start = as.Date(d.t)
      tmp1$date.end = d2.t
      tmp1$lwr = .3; tmp1$upr = .5; 
      tmp1$SR.level = 'Extra1.first1wk'
      tmp2$date.start = as.Date(d2.t) + 1
      tmp2$SR.level = 'normal'
      alpha_bounds.t = rbind(alpha_bounds.t[period != 'Delta'], tmp0, tmp1, tmp2)
      
      alpha_bounds.t = alpha_bounds.t[order(date.start)]
      
    }
    
  # for beta, check if there is potential super-spreading events?
  wdelta = seq(as.Date(beta_bounds.t[period=='Delta']$date.start), as.Date(beta_bounds.t[period=='Delta']$date.end), by = 'day')
  dholi = seq(as.Date('2021/07/01'),as.Date('2021/07/20'), by = 'day') # ~4th of July
  if(any(wdelta %in% dholi) & any(dates.vincr %in% dholi)){
    ds.t = dates.vincr[(dates.vincr %in% wdelta) & (dates.vincr %in% dholi)]
    d.t = ds.t[1]
    if(as.Date(tail(ds.t,1)) == as.Date(d.t)){ # 1 week
      nwk.t = 1
      d2.t = as.Date(d.t) + 6 + 7 * nwk.t # extend 1 week
      
    } else {
      nwk.t = 2
      d2.t = as.Date(d.t) + 13 + 7 * nwk.t  # 2 weeks at most
      
    }
    tmp0 = tmp1 = tmp2 = beta_bounds.t[period == 'Delta']
    tmp0$date.end = as.Date(d.t) - 1
    tmp1$date.start = as.Date(d.t)
    tmp1$date.end = d2.t
    tmp1$lwr = beta_bounds_ini[1] * 1.9 # increase the lower bound
    tmp1$upr = beta_bounds_ini[2] * 1.7 
    # tmp1$SR.level = paste0('Extra.first1wk') # too much!
    # tmp1$SR.level = paste0('extraS.first',nwk.t,'wk')
    tmp1$SR.level = paste0('extra.first',nwk.t,'wk')
    if(nwk.t==2){
      tmp3 = tmp2 # add another one
      tmp2$date.start = as.Date(d2.t) + 1
      tmp2$date.end = as.Date(d2.t) + (nwk.t+1) * 7
      tmp2$SR.level = paste0('extra.first',nwk.t,'wk')
      tmp2$upr = tmp2$upr * .8
      
      tmp3$date.start = tmp2$date.end + 1
      tmp3$SR.level = 'normal'
      
      tmp2 = rbind(tmp2, tmp3)
    } else {
      tmp2$date.start = as.Date(d2.t) + 1
      tmp2$SR.level = paste0('Extra.first1wk') # allow it to go down
      # tmp2$SR.level = paste0('extraS.first',nwk.t,'wk') # allow it to go down
      # tmp2$SR.level = paste0('extra.first',nwk.t,'wk') # allow it to go down
    }
    
    

    beta_bounds.t = rbind(beta_bounds.t[period != 'Delta'], tmp0, tmp1, tmp2)
    beta_bounds.t = beta_bounds.t[order(date.start)]
    
    # if a super-spreading event is possible, cross-check with detection rate
    if(nrow(alpha_bounds.t[period=='Delta' & SR.level=='Extra1.first1wk'])>0){
      ds1.t = seq(as.Date(ds.t[1]),as.Date(tail(ds.t,1)), by = 'day')
      ds2.t = seq(as.Date(alpha_bounds.t[period=='Delta' & SR.level=='Extra1.first1wk']$date.start),
                 as.Date(alpha_bounds.t[period=='Delta' & SR.level=='Extra1.first1wk']$date.end), by='day')
      if(sum(ds2.t %in% ds1.t)>=7){
        ii = which(alpha_bounds.t$period=='Delta' & alpha_bounds.t$SR.level=='Extra1.first1wk')
        # delay it
        alpha_bounds.t[ii]$date.start = as.Date(d.t[1]+6)
        alpha_bounds.t[ii-1]$date.end = alpha_bounds.t[ii]$date.start - 1
        # not good
        # alpha_bounds.t[ii]$SR.level = 'normal'
        # alpha_bounds.t[ii]$lwr = alpha_bounds[period=='Delta']$lwr
        # alpha_bounds.t[ii]$upr = alpha_bounds[period=='Delta']$upr
      }
    }
    
    
  }
         
  
  # vtimes = vda.t %>% setnames('v.main', 'type')
  # use the 10% mark for probing S
  vtimes = vda.t %>% filter(type == '10perc') # setnames('v.main', 'type')
  vtimes[variant == 'Others']$variant = 'Wildtype'
  vtimes$SR.level = ''; vtimes$parm = 'wave.start'; vtimes$lwr = ''; vtimes$upr = ''; vtimes$period = ''
  # vchar = vtimes
  vchar = vda.t %>% filter(type == 'predominant') # setnames('v.main', 'type')
  vchar[variant == 'Others']$variant = 'Wildtype'
  vchar$SR.level = ''; 
  vchar$lwr = ''; vchar$upr = ''; vchar$period = ''
  vchar$parm = paste0('variant-',tolower(vchar$variant))
  vchar$type = 'variant characterization'
  
  # also set filtering settings for S
  bbS = data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'no.imm.escape', variant = NA,
                   date.start = as.Date('2020/01/01'), 
                   date.end = as.Date('2021/05/31') # as.Date(vtimes[type =='Delta']$date.start) - 30
                   )
  if(nrow(vtimes[variant=='Delta'])>0){
    ii = which(vtimes$variant=='Delta')
    
    # in case there is cocirculation with Gamma
    ii.gamma = which(vtimes$variant=='Gamma')
    
    # bbS$date.end = as.Date(vtimes[ii]$date.start) - 45
    
    # start earlier to allow S update during alpha period too? 
    bbS$date.end = as.Date(min(vtimes[variant %in% c('Iota','Epsilon', 'Alpha','Gamma','Beta')]$date.start))-1 # - 45
    
    # d.reset = as.Date(vtimes[ii]$date.start) # - 14
    # 8/26/22 in case there is cocirculation with Gamma
    d.reset = as.Date(min(vtimes[c(ii, ii.gamma)]$date.start))
    
    if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
      # find the first date cases > the low point or if it's increasing substantially
      d.reset = da.t$date[(! da.t$date %in% dates.low.case & da.t$date >= d.reset) |
                            (da.t$date %in% dates.incr & da.t$date >= d.reset)
                          ] %>% unique %>% sort %>% head(1)
    } else {
      # try 1-wk prior
      d.reset = as.Date(vtimes[ii]$date.start) - 7
      if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
        # find the first date cases > the low point
        # d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
        # find the first date cases > the low point or if it's increasing substantially
        d.reset = da.t$date[(! da.t$date %in% dates.low.case & da.t$date >= d.reset) |
                              (da.t$date %in% dates.incr & da.t$date >= d.reset)
        ] %>% unique %>% sort %>% head(1)
      } else {
        # try 1-wk prior
        d.reset = as.Date(vtimes[ii]$date.start) - 7 # - 14
        if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
          # find the first date cases > the low point
          # d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
          # find the first date cases > the low point or if it's increasing substantially
          d.reset = da.t$date[(! da.t$date %in% dates.low.case & da.t$date >= d.reset) |
                                (da.t$date %in% dates.incr & da.t$date >= d.reset)
          ] %>% unique %>% sort %>% head(1)
        } 
      }
    }
    d.reset = pmax(d.reset - 7, vtimes[ii]$date.start) # 1 week prior
    
    ii.omiBA1 = which(vtimes$variant=='Omicron_BA.1')
    
    bbS = rbind(bbS, 
                data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'imm.escape', variant = NA,
                           date.start = as.Date(vtimes[ii]$date.start) - 45 + 1, 
                           date.end = as.Date(vtimes[ii]$date.start) + 30 * 3), # starting a month before delta, allow the system to adjust first
                data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'no.imm.escape', variant = NA,
                           date.start = as.Date(vtimes[ii]$date.start) + 30 * 3 + 1, # allow to adjust for 3 months
                           date.end = as.Date(vtimes[ii.omiBA1]$date.start) - 30), 
                data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'reset cntSR.S', variant = 'Delta',
                           date.start = d.reset, # 4 weeks before it becomes dominant
                           date.end = d.reset+6))
  }
  if(nrow(vtimes[variant == 'Omicron_BA.2.12.1'])>0){ # omicron BA.1 is used as reference so no need to probe, start with BA.2.12.1
    ii = which(vtimes$variant=='Omicron_BA.2.12.1')
    
    d.reset = as.Date(vtimes[ii]$date.start) # - 14
    if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
      # find the first date cases > the low point
      d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
    } else {
      # try 1-wk prior
      d.reset = as.Date(vtimes[ii]$date.start) - 7
      if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
        # find the first date cases > the low point
        d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
      } else {
        # try 2-wks prior
        d.reset = as.Date(vtimes[ii]$date.start) - 14
        if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
          # find the first date cases > the low point
          d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
        } 
      }
    }
    
    bbS = rbind(bbS, 
                data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'imm.escape', variant = NA,
                           date.start = as.Date(vtimes[ii]$date.start) - 30 + 1,  
                           date.end = max(vtimes$date.end)), # allow the system to adjust throughout afterwards
                data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'reset cntSR.S', variant = 'Omicron_BA.2.12.1',
                           date.start =d.reset, # 4 weeks before it becomes dominant
                           date.end = d.reset+6))
  }
  
  if(nrow(vtimes[variant == 'Omicron_nonBA.1o2'])>0){
    ii = which(vtimes$variant=='Omicron_nonBA.1o2')
    
    d.reset = as.Date(vtimes[ii]$date.start) # - 14
    if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
      # find the first date cases > the low point
      d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
    } else {
      # try 1-wk prior
      d.reset = as.Date(vtimes[ii]$date.start) - 7
      if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
        # find the first date cases > the low point
        d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
      } else {
        # try 2-wks prior
        d.reset = as.Date(vtimes[ii]$date.start) - 14
        if(d.reset %in% dates.low.case){ # cases are very low at the moment, so not accurate
          # find the first date cases > the low point
          d.reset = da.t$date[! da.t$date %in% dates.low.case & da.t$date >= d.reset] %>% head(1)
        } 
      }
    }
    
    bbS = rbind(bbS, 
                data.table(parm = 'S', type = 'filtering', SR.level = '', period = 'reset cntSR.S', variant = 'Omicron_subvariants',
                           date.start =d.reset, # 4 weeks before it becomes dominant
                           date.end = d.reset+6))
  }
  
  
  bbS$location = loc.t; bbS$lwr = ''; bbS$upr = ''
  bbS$date.start = bbS$date.start %>% as.Date # %>% format('%Y-%m-%d')
  bbS$date.end = bbS$date.end %>% as.Date #  %>% format('%Y-%m-%d')
  # bbS$type = NULL
  # assemble it
  vtimes$type = NULL
  setnames(vtimes, 'variant', 'type')
  bounds.t = rbind(vtimes, alpha_bounds.t, beta_bounds.t, ifr_bounds.t, alpha_DAbounds.t, beta_DAbounds.t, ifr_DAbounds.t, 
                   bbS, vchar, fill = T)
  
  BOUNDS = rbind(BOUNDS, bounds.t)
}

write_xlsx(BOUNDS, paste0(dir_code, 'parm.bounds_pop.den',tag,'.xlsx'))

