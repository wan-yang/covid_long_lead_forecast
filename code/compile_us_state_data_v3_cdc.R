# compile COVID data
# 3/23/23 - use CDC data for case/death instead
# b/c NYT data were problematic for some states and weeks later on

data.tag = 'cdc'

dir_code = './code/'
dir_data = './data/'


source(paste0(dir_code, 'Fn_util.R'))

library("RSocrata")
library(data.table); library(magrittr)
library(rworldmap)
library(classInt)
library(RColorBrewer)
library(mapdata)
library(maptools) # for shapefiles
library(scales) # for transparency
library(plotrix); # for color.legend
# library(graphicsQC)
library(rgdal); # for projection system conversion
library(sp);
library(TeachingDemos);
library('readr')
library('readxl')
library('writexl')
library(MMWRweek)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RJSONIO)


# cdc case and death data
# CAUTION!!! THE CDC DATA, EACH WEEK STARTS FROM THURSDAY AND ENDS ON WEDNESDAY
# data source: https://data.cdc.gov/Case-Surveillance/United-States-COVID-19-Cases-and-Deaths-by-State-o/9mfq-cb36
# https://data.cdc.gov/Case-Surveillance/Weekly-United-States-COVID-19-Cases-and-Deaths-by-/pwn4-m3yp
# There are currently 60 public health jurisdictions reporting cases of COVID-19. 
# This includes the 50 states, the District of Columbia, New York City, the U.S. territories of American Samoa, 
# Guam, the Commonwealth of the Northern Mariana Islands, Puerto Rico, and the U.S Virgin Islands
# as well as three independent countries in compacts of free association with the United States, 
# Federated States of Micronesia, Republic of the Marshall Islands, and Republic of Palau.
# New York State’s reported case and death counts do not include New York City’s counts as 
# they separately report nationally notifiable conditions to CDC.
# new_case = new confirmed case + new probable case 
# new_death = new confirmed death + new probable death 
# March 16, 2023: Due to technical difficulties, Arkansas’s aggregate case and death data will be reported as 0. As a result, case and death metrics will appear lower than expected in the March 16, 2023, weekly release. 
# March 16, 2023: Due to potential technical difficulties that are being reviewed with the state, 
#  Texas's aggregate case and death data will be reported as 0. As a result, case and death metrics will appear lower than expected in the March 16, 2023, weekly release.
# March 30, 2023: Due to potential technical issues, Florida was unable to report aggregate case and death data to CDC. As a result, case and death metrics will be reported as 0 in the March 30, 2023, weekly release.

df <- read.socrata(
  # "https://data.cdc.gov/resource/9mfq-cb36.json",
  "https://data.cdc.gov/resource/pwn4-m3yp.json",
  app_token = "YOUR TOKEN",
  email     = "YOUR EMAIL ACCOUNT",
  password  = "YOUR PW"
)

setDT(df)
df = df %>% mutate(tot_cases = tot_cases %>% as.numeric,
                   new_cases = new_cases %>% as.numeric,
                   tot_deaths = tot_deaths %>% as.numeric,
                   new_deaths = new_deaths  %>% as.numeric,
                   new_historic_cases = new_historic_cases  %>% as.numeric,
                   new_historic_deaths = new_historic_deaths %>% as.numeric
)
df.wkly = df %>% dplyr::select(state, start_date, new_cases, new_deaths) %>%
  mutate(date = as.Date(start_date) + 3)
# EACH WEEK STARTS FROM THURSDAY AND ENDS ON WEDNESDAY
weeks = df.wkly %>% .$date %>% MMWRweek() %>% setnames(paste0('MMWR', c('year', 'week', 'day')), c('year', 'week', 'day')) 
df.wkly = df.wkly %>% cbind(weeks, .)  %>% setDT
# combine NY exclude NYC and NYC

da.nys = df.wkly %>% dplyr::select(state, date, year, week, new_cases, new_deaths) %>% 
  dplyr::filter(state %in% c('NY','NYC')) %>%
  reshape2::melt(., id.vars = c('state', 'date', 'year','week')) %>%
  setnames(., c('variable'), 'data.type') %>% 
  dplyr::filter(data.type %in% c('new_cases', 'new_deaths')) %>%
  group_by(date, year, week, data.type) %>% 
  summarise(value = sum(value)) %>%
  mutate(state = 'New York') %>% ungroup()

da_cdc = df.wkly %>% dplyr::select(state, date, year, week, new_cases, new_deaths) %>%
  dplyr::filter(!state %in% c('NY','NYC')) %>%
  mutate(state = factor(state, levels = state.abb, labels = state.name)) %>%
  reshape2::melt(., id.vars = c('state', 'date', 'year','week')) %>%
  setnames(., c('variable'), 'data.type') %>%
  rbind(., da.nys) %>%
  dplyr::filter(data.type %in% c('new_cases', 'new_deaths')) %>%
  mutate(data.type = factor(data.type, levels = c('new_cases', 'new_deaths'), labels = c('case','death'))) %>%
  filter(!is.na(state)) %>%
  setcolorder(., c("state", "date", "year","week","value", "data.type"))

write.csv(da_cdc, paste0(dir_data, 'da_raw_case_death_us_cdc.csv'), row.names = F)

tmp = da_cdc %>% filter(date >= as.Date('2023/3/1'))
  

url.t = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv'
case_hopkins = read_csv(url(url.t)) %>% data.table() 

url.t = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv'
death_hopkins = read_csv(url(url.t)) %>% data.table() 

url.t = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv'
pop = read_csv(url(url.t)) %>% data.table() 

url.t = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv'
vx = read_csv(url(url.t)) %>% data.table()

url.t = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv'
da_nyt = read_csv(url(url.t)) %>% data.table() # cross-checking with new york times data

# variant data
url.t = 'https://raw.githubusercontent.com/hodcroftlab/covariants/master/cluster_tables/USAClusters_data.json'
d.variant = fromJSON(url.t) 
d.variant = d.variant$countries # %>% data.table()
nloc = length(d.variant)
locs = names(d.variant)
res = NULL
for(i in 1:nloc){
  
  da.t = d.variant[[i]]
  cnames = names(da.t)
  da.t = da.t %>% unlist %>% matrix(nrow=length(da.t), byrow=T) %>% t %>% data.table()
  colnames(da.t) = cnames
  # convert to percentage?
  variants = cnames[!cnames %in% c('week','total_sequences')]
  percs = NULL
  for(v in variants){
    percs = cbind(percs, as.numeric(unlist(da.t[,v,with=F])) / as.numeric(da.t$total_sequences))
  }
  colnames(percs) = variants
  percs = percs %>% data.table()
  percs$week = da.t$week
  percs$location = locs[i]
  setcolorder(percs, c('location', 'week', variants))
  
  res = rbind(res, percs, fill = T)
}
# combine different delta and omicron clusters?
col.omi = colnames(res)[grep('Omicron', colnames(res))]
col.omi.nonBA1 = col.omi[!grepl('21K', col.omi)]
col.omi.nonBA1o2 = col.omi[!grepl('21K|21L|22C', col.omi, perl =T)]
# https://covariants.org/variants/21K.Omicron
# Omicron: 21K(Omicron) or BA.1, 21L(Omicron) or BA.2, 22A(Omicron) or BA.4, 22B(Omicron) or BA.5, 22C(Omicron) or BA.2.12.1, 22D(Omicron) or 
res$Delta = rowSums(res[,grep('Delta', colnames(res)),with=F], na.rm = T)
res$Omicron = rowSums(res[,grep('Omicron', colnames(res)),with=F], na.rm = T)
res$Omicron_BA.1 = rowSums(res[,'21K (Omicron)',with=F], na.rm = T)
res$Omicron_BA.2 = rowSums(res[,'21L (Omicron)',with=F], na.rm = T)
res$Omicron_BA.2.12.1 = rowSums(res[,'22C (Omicron)',with=F], na.rm = T)
res$Omicron_nonBA.1o2 = rowSums(res[,col.omi.nonBA1o2,with=F], na.rm = T)
res$Alpha = rowSums(res[,grep('Alpha', colnames(res)),with=F], na.rm = T)
res$Beta = rowSums(res[,grep('Beta', colnames(res)),with=F], na.rm = T)
res$Gamma = rowSums(res[,grep('Gamma', colnames(res)),with=F], na.rm = T)
res$Iota = rowSums(res[,grep('Iota', colnames(res)),with=F], na.rm = T)
res$Epsilon = rowSums(res[,grep('Epsilon', colnames(res)),with=F], na.rm = T)
write.csv(res, paste0(dir_data, 'perc.variant.biwkly_by_state.csv'), row.names = F)

tmp = res %>% filter(location=='California') %>% dplyr:: select(week,Omicron)
tmp$Omicron = tmp$Omicron * 100
tmp$ma2left = stats::filter(tmp$Omicron[-1], filter = rep(1/2,2)) %>% append(NA, values = .)
tmp$ma2right = stats::filter(tmp$Omicron[-nrow(tmp)], filter = rep(1/2,2)) %>% append(x =., values = NA)

# to estimate the weekly proportion, by smoothing
variants_vec = c('Delta', 'Omicron','Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.2.12.1', 'Omicron_nonBA.1o2', 'Alpha', 'Beta', 'Gamma', 'Iota', 'Epsilon')
variants_vec.nooverlap = c('Delta', 'Omicron','Alpha', 'Beta', 'Gamma', 'Iota', 'Epsilon')
states = res$location %>% unique()
res.wkly = NULL
for(loc in states){
  for(v in variants_vec){
    tmp = eval(parse(text = paste0('res %>% filter(location==loc) %>% dplyr:: select(week,', v,')'))) 
    tmp = tmp %>% setnames(., c('week', v), c('week', 'perc'))
    # only use data that the variant is non-zero
    tmp$week = tmp$week %>% as.Date()
    
    # skip if there is no detection at all
    if(nrow(tmp %>% filter(perc > 0)) < 1)
      next
    
    if(F){
      d.start = tmp %>% filter(perc > 0) %>% .$week %>% min
      d.end = tmp %>% filter(perc > 0) %>% .$week %>% max
      tmp = tmp %>% filter(week >= d.start-14 & week <=d.end)
    }
    
    d.start = tmp %>% filter(perc > 0) %>% .$week %>% min
    min.t = tmp %>% filter(perc > 0) %>% .$perc %>% min
    
    
    
    
    # tmp$id = (1: nrow(tmp) -1) *2
    tmp$id = as.numeric(tmp$week - tmp$week[1]) / 7  # for potential missing data
    
    
    m1 = stats::smooth.spline(x = tmp$id, y = tmp$perc)  # , all.knots = F
    tmp.pred = data.table(id = seq(0, max(tmp$id), by = 1), perc = NA)
    tmp.pred$perc = predict(m1, x = tmp.pred$id)$y
    
    if(F){
      plot(tmp$id, tmp$perc)
      lines(tmp.pred$id, tmp.pred$perc)
    }
    

    # it looks ok
    # check DA
    DAneg = tmp.pred[perc < 0]
    DApos = tmp.pred[perc > 1]
    flag = F
    if(any(DAneg$perc < -.2)){
      flag = T
    } else if(any(DApos$perc > 1.2)){
      flag = T
    } 
    if(flag)
      print('need QC')
    
    tmp.pred[perc < 0]$perc = 0
    tmp.pred[perc > 1]$perc = 1
    
    # tmp.pred$week = seq(tmp$week[1] + 7, tmp[nrow(tmp)]$week + 7, by = 'week') # shift it by 1 week to center the data
    tmp.pred$week = tmp$week[1] + 7 + tmp.pred$id * 7 - 1 # make the week starts Sundays
    tmp.pred[perc < min.t / 2 & week < d.start - 14]$perc = 0
    
    # for the initial detection, perc tends to be very low and may be rounded to 0 
    # this affect the initial case/death apportion
    # fill in this gap a bit?
    d.start.est = tmp.pred %>% filter(perc > 0) %>% .$week %>% head(., 1)
    if((d.start - 1) < d.start.est){
      # the observed start (2-wk span) is earlier than the estimated -> all the detection is assigned to the 2nd week
      
      idx0 = which(tmp.pred$week == d.start.est)
      
      # tmp0 = tmp.pred[idx0 + c(-1:1)]
      # assume exponential increase?
      # log(perc1/perc0) = log(perc2/perc1)
      # p2 = p1 * 2^(7/tau)
      # p1 = p0 * 2^(7/tau)
      tau = 7 * log(2) / log(tmp.pred$perc[idx0+1]/ tmp.pred$perc[idx0])
      tmp.pred$perc[idx0-1] = tmp.pred$perc[idx0]/2^(7/tau)
      
      print(paste(loc, 'adjust for', v))
    }
    
    if(F){
      plot(tmp.pred$week, tmp.pred$perc, type = 'l')
      points(tmp$week, tmp$perc)
    }
    
    
    # save it
    res.wkly = rbind(res.wkly, 
                     data.table(location = loc, variant = v, week = tmp.pred$week, perc = tmp.pred$perc, needQC = flag))
    
    rm(tmp, tmp.pred, flag)
    
  }
}

tmpqc = res.wkly %>% filter(needQC == T)
# only 2 instances: Vigin Island & Alpha (neg set 0, not a prob)
# Northern Mariana Islands & Alpha (neg set 0 and >1 set to 1, not a prob)

res.wkly2 = res.wkly %>% dcast(., location + week ~ variant, value.var = 'perc')
omis = res.wkly2 %>% dplyr::select(location, week, Omicron, Omicron_BA.1, Omicron_BA.2, Omicron_BA.2.12.1, Omicron_nonBA.1o2) # check the numbers

# scale to max of 1
ff = data.table(fs = rowSums(res.wkly2[,variants_vec.nooverlap, with=F]))
ff = ff %>% mutate(f = case_when(fs >1 ~ 1/fs,
                                 T ~ 1))

res.wkly2[, variants_vec.nooverlap] = res.wkly2[, variants_vec.nooverlap, with=F] * ff$f

# further scale the omicron subvariants to the total
ff.omi = data.table(fs = rowSums(res.wkly2 %>% dplyr::select(Omicron_BA.1, Omicron_BA.2, Omicron_BA.2.12.1, Omicron_nonBA.1o2)), 
                    omi.total = res.wkly2$Omicron)
ff.omi = ff.omi %>% mutate(f = case_when(fs != omi.total ~ omi.total/fs,
                                 T ~ 1)) %>%
  mutate(f = case_when(is.infinite(f) ~ 1, 
                       T ~ f))

res.wkly2[, c('Omicron_BA.1', 'Omicron_BA.2','Omicron_BA.2.12.1', 'Omicron_nonBA.1o2')] = res.wkly2[, c('Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.2.12.1', 'Omicron_nonBA.1o2'), with=F] * ff.omi$f
# check again
omis = res.wkly2 %>% dplyr::select(location, week, Omicron, Omicron_BA.1, Omicron_BA.2, Omicron_BA.2.12.1, Omicron_nonBA.1o2) # check the numbers


res.wkly2$Others = 1 - rowSums(res.wkly2[,variants_vec.nooverlap, with=F])
any(res.wkly2$Others < 0)
tmp = res.wkly2 %>% filter(Others <0) # rounding err
res.wkly2[Others < 0]$Others = 0 # set to 0

tmp = res.wkly2 %>% filter(location == 'New York' & week > as.Date('2021/11/1'))
  
summary(rowSums(res.wkly2[,c('Delta', 'Omicron_BA.1', 'Omicron_BA.2','Omicron_BA.2.12.1', 'Omicron_nonBA.1o2', 'Alpha', 'Beta', 'Gamma', 'Iota', 'Epsilon','Others'), with=F]))
summary(rowSums(res.wkly2[,c('Delta', 'Omicron', 'Alpha', 'Beta', 'Gamma', 'Iota', 'Epsilon','Others'), with=F]))
# save
write.csv(res.wkly2, paste0(dir_data, 'perc.variant.wkly.smoothed_by_state.csv'), row.names = F)


# identify the dominant variant for each week 
fn_v.main = function(x, variants = c('Alpha','Beta','Gamma','Iota','Epsilon','Delta','Omicron','Others')){
  variants[which.max(x)]
}

res.wkly2$v.main = res.wkly2[,c('Alpha','Beta','Gamma','Iota','Epsilon','Delta','Omicron','Others'),with = F] %>% apply(1, fn_v.main)
res.wkly2$v.main2 = res.wkly2[,c('Alpha','Beta','Gamma','Iota','Epsilon','Delta','Omicron_BA.1','Omicron_BA.2', 'Omicron_BA.2.12.1','Omicron_nonBA.1o2','Others'),with = F] %>% 
  apply(1, fn_v.main, variants = c('Alpha','Beta','Gamma','Iota','Epsilon','Delta','Omicron_BA.1','Omicron_BA.2','Omicron_BA.2.12.1','Omicron_nonBA.1o2','Others'))


# find the critical times
date.start.unk = as.Date('2020/1/1')
date.end.unk = as.Date('2022/12/31')
date.end.unk = Sys.Date() # update 
res4 = NULL
for(loc.t in locs){
  da.t = res.wkly2[location == loc.t]
  for(v in c('Alpha','Beta','Gamma','Iota','Epsilon','Delta','Omicron','Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.2.12.1', 'Omicron_nonBA.1o2','Others')){
    if(v %in% c('Alpha','Beta','Gamma','Iota','Epsilon','Delta','Omicron', 'Others')){
      tmp = da.t[v.main == v]
    } else {
      tmp = da.t[v.main2 == v]
    }
   
    # tmp2 = da.t %>% filter(week > as.Date('2021/11/15'))
    tmp10perc = eval(parse(text = paste('da.t %>% filter(',v, '> .1)')))
    tmp25perc = eval(parse(text = paste('da.t %>% filter(',v, '> .25)')))
    tmp50perc = eval(parse(text = paste('da.t %>% filter(',v, '> .5)')))
    # at least 10%
    if(nrow(tmp10perc) > 0){
      
      if(tail(tmp10perc$week,1) == max(da.t$week)){
        date.end.t = date.end.unk
      } else {
        date.end.t = as.Date(tail(tmp10perc$week,1))+6
      }
      
      if(tmp10perc$week[1] == min(da.t$week)){
        date.start.t = date.start.unk
      } else {
        date.start.t = as.Date(tmp10perc$week[1])
      }
      
      res4 = rbind(res4, data.table(location = loc.t, variant = v, type = '10perc',
                                    date.start = date.start.t, date.end = date.end.t))
    }
    # at least 25%
    if(nrow(tmp25perc) > 0){
      
      if(tail(tmp25perc$week,1) == max(da.t$week)){
        date.end.t = date.end.unk
      } else {
        date.end.t = as.Date(tail(tmp25perc$week,1))+6
      }
      
      if(tmp25perc$week[1] == min(da.t$week)){
        date.start.t = date.start.unk
      } else {
        date.start.t = as.Date(tmp25perc$week[1])
      }
      
      res4 = rbind(res4, data.table(location = loc.t, variant = v, type = '25perc',
                                    date.start = date.start.t, date.end = date.end.t))
    }
    # at least 50%
    if(nrow(tmp50perc) > 0){
      
      if(tail(tmp50perc$week,1) == max(da.t$week)){
        date.end.t = date.end.unk
      } else {
        date.end.t = as.Date(tail(tmp50perc$week,1))+6
      }
      
      if(tmp50perc$week[1] == min(da.t$week)){
        date.start.t = date.start.unk
      } else {
        date.start.t = as.Date(tmp50perc$week[1])
      }
      
      res4 = rbind(res4, data.table(location = loc.t, variant = v, type = '50perc',
                                    date.start = date.start.t, date.end = date.end.t))
    }
    # predominant
    if(nrow(tmp) > 0){
      
      if(tail(tmp$week,1) == max(da.t$week)){
        date.end.t = date.end.unk
      } else {
        date.end.t = as.Date(tail(tmp$week,1))+6
      }
      
      if(tmp$week[1] == min(da.t$week)){
        date.start.t = date.start.unk
      } else {
        date.start.t = as.Date(tmp$week[1])
      }
      
      res4 = rbind(res4, data.table(location = loc.t, variant = v, type = 'predominant',
                                    date.start = date.start.t, date.end = date.end.t))
    }
  }
}

res4 = res4[order(location, variant, date.start)]

write.csv(res4, paste0(dir_data, 'tm.variant_by_state_smoothed.wkly.csv'), row.names = F)

# mobility - no updates from 10/15/22 onwards
Fn_checkNA = function(x){any(is.na(x))}

N = 1e6

# aggregate to weekly intervals
# states = case$Province_State %>% unique %>% .[!(. %in% c("American Samoa", "Diamond Princess","Grand Princess"))]
# states = pop$Province_State  %>% unique
states = da_cdc$state %>% unique
pop.t = pop %>% filter(iso2 == 'US' & Province_State %in% states & is.na(Admin2)) %>% dplyr::select(Province_State, Population) %>% data.table() #  %>% sum(na.rm = T)
# normalized to per 1 M

res = da_cdc %>% 
  left_join(x = ., 
            y = pop.t %>% data.table::copy(.) %>% setnames('Province_State', 'state'), x.all = T, by = 'state') %>%
  mutate(value = value / Population * N) %>%
  mutate(Population = NULL)

write.csv(res, paste0(dir_data, 'da_case_death_us_cdc.csv'), row.names = F)


# compute variant specific case and death
da = read.csv(paste0(dir_data, 'da_case_death_us_cdc.csv'))
da$date = da$date %>% as.Date

# check the time lag from case to death
tda = da %>% reshape2::dcast(., state + date ~ data.type, value.var = 'value')
lags = NULL
for(state.t in (tda$state %>% unique)){
  ttda0 = tda %>% filter(state == state.t & date >= as.Date('2020/6/1') & date < as.Date('2021/12/1'))
  tmp0 = ccf(ttda0$death, ttda0$case, lag.max = 8, plot = F) 
  ttda1 = tda %>% filter(state == state.t & date >= as.Date('2020/3/1') & date < as.Date('2021/12/1'))
  tmp1 = ccf(ttda1$death, ttda1$case, lag.max = 8, plot = F) 
  ttda2 = tda %>% filter(state == state.t & date >= as.Date('2022/1/1'))
  tmp2 = ccf(ttda2$death, ttda2$case, lag.max = 8, plot = F) 
  lags = rbind(lags, 
               data.table(state = state.t, 
                          lag.max.beforeOmicron0 = tmp0$lag[tmp0$acf %>% which.max()],
                          lag.max.beforeOmicron = tmp1$lag[tmp1$acf %>% which.max()],
                          lag.max.Omicron = tmp2$lag[tmp2$acf %>% which.max()]
                          )
               )
}

da.variant = read.csv(paste0(dir_data, 'perc.variant.wkly.smoothed_by_state.csv')) %>% 
  data.table %>% setnames(., c('location','week'), c('state', 'date'))
variant.t = 'Omicron'
res = NULL
res2 = NULL # for the rest of it, i.e. non-omicron
for(d.type in c('case','death')){
  
  
  tda = da %>% filter(data.type == d.type) 
  da.variant.t = eval(parse(text = paste0('da.variant %>% dplyr::select(date, state,', variant.t,')'))) %>% 
    data.table::copy() %>% setnames(., variant.t, 'perc')
  
  da.variant.t$date = da.variant.t$date %>% as.Date
  tmp = da.variant.t %>% reshape2::dcast(., date ~ state, value.var = 'perc')
  tmp2 = data.table(date = seq(max(tmp$date)+7, Sys.Date(), by = 'week'), 
                    tmp[1, -1])
  tmp2 = tmp2 %>% reshape2::melt(., id.vars = c('date')) %>% setnames(., c('variable','value'), c('state', 'perc')) %>%
    mutate(perc = 1) %>% # set all to 1
    reshape2::dcast(., date ~ state, value.var = 'perc')
  
  tmp = rbind(tmp, tmp2)
  tmp = tmp %>% reshape2::melt(., id.vars = c('date')) %>% setnames(., c('variable','value'), c('state', 'perc')) %>%
    mutate(perc = case_when(date > as.Date('2022/12/31') & is.na(perc) ~ 1, 
                            T ~ perc))
  
  da.variant.t = tmp
  
  # 3/23/23 
  # since we are only separating non-Omicron and Omicron, 
  # set those after 2022 to Omicron [UNLESS A NEW NON-OMICRON VARIANT EMERGES IN LATER YEARS]
  
  # exclude last week without variant data
  date.max = da.variant.t %>% group_by(state) %>% summarise(date.max = max(date))
  tda = tda %>% left_join(x = ., y = date.max, by = 'state') %>% 
    filter(date <= date.max) %>% mutate(date.max = NULL)
  
  # for death, shift the time
  if(d.type == 'death'){
    da.variant.t$date = da.variant.t$date + 21
  } 
  
  # merge the case/death data with variant
  
  # for the complement - do it before the time shift? no, should be also afterwards
  tda2 = tda %>% left_join(x = .,
                           y = da.variant.t,
                           by = c('state', 'date')
  ) %>% mutate(perc = case_when(is.na(perc) ~ 0,
                                T ~ perc)) %>% 
    dplyr::mutate(value = value * (1-perc))
  
  # for the variant of interest
  tda = tda %>% left_join(x = .,
                          y = da.variant.t,
                          by = c('state', 'date')
                          ) %>% 
    # dplyr::mutate(value = case_when(value == 0 & perc >0 ~ NA_real_,
    #                                T ~ value)) %>% 
    dplyr::mutate(value = value * perc)
  
  # perc > 0 and case == 0, e.g. Florida happened to have 0 cases but perc > 0
  # for flagging, these are set to NA
  # tda %>% filter(is.na(value) & !is.na(perc))

  
  # only include weeks the variant has been detected
  tda.v = vdates[variant == variant.t] %>% dplyr::select(state, date.1st) %>% data.table::copy(.) # %>% setnames(., 'date.1st', 'date')
  
  
  
  tda = tda %>% left_join(x = ., y = tda.v, by = 'state') %>% 
    filter(date >= date.1st)
  
  res = rbind(res, tda)
  res2 = rbind(res2, tda2)

}

write.csv(res, paste0(dir_data,'da_case_death_us_cdc_', variant.t,'.csv'), row.names = F)
write.csv(res2, paste0(dir_data,'da_case_death_us_cdc_non.', variant.t,'.csv'), row.names = F)


# also update NYT and Hopkins data for comparison
# NYT
res = NULL
for(da.type in c('case', 'death')){
  
  if(da.type == 'case'){
    tda = da_nyt %>% dplyr::select(date, state, cases) %>% 
      reshape2::dcast(., date ~ state, value.var = 'cases') %>% data.table()
  } else if(da.type == 'death'){
    tda = da_nyt %>% dplyr::select(date, state, deaths) %>% 
      reshape2::dcast(., date ~ state, value.var = 'deaths') %>% data.table()
  }
  
  
  tda = data.table(date = tda$date[-1] %>% as.Date, 
                   tda[-1,2:ncol(tda), with=F] - tda[-nrow(tda),2:ncol(tda), with=F]) %>% 
    filter(date < Sys.Date() - (format(Sys.Date(), '%w') %>% as.numeric())) # exclude week with incomplete data
  
  tmp = tda %>% .$date %>% MMWRweek %>% setnames(c('MMWRyear','MMWRweek','MMWRday'),c('year','week','day'))
  tda = cbind(tmp, tda) %>% data.table()
  tda = tda[, lapply(.SD, sum, na.rm=T), .SDcol = 5:ncol(tda), by = c('year', 'week')]
  tda = data.table(date = MMWRweek2Date(MMWRyear = tda$year, MMWRweek = tda$week), tda)
  
  matplot(tda[,4:ncol(tda), with=F], type = 'l')
  tda = tda %>% melt(., id.vars = c('date', 'year', 'week')) %>% setnames('variable', 'state')
  tda[value < 0]  # some <0
  tda[value < 0]$value = 0 # set it to 0
  
  # normalized to per 1 M
  tda = merge(x = tda, y = pop.t %>% data.table::copy(.) %>% setnames('Province_State', 'state'), x.all = T, by = 'state') %>% 
    mutate(value = value / Population * N) %>%
    mutate(Population = NULL)
  
  tda$data.type = da.type
  
  res = rbind(res, tda)
  
  # save(tda2, file = paste0(dir_data, 'da_wkly.',da.type,'_us.RData'))
}
write.csv(res, paste0(dir_data, 'da_case_death_us_nyt.csv'), row.names = F)
# Hopkins
res = NULL
for(da.type in c('case', 'death')){
  da.t = get(paste0(da.type,'_hopkins'))
  
  tda = da.t %>% filter(iso2 == 'US' & Province_State %in% states)  # county-level data
  tmp = tda %>% filter(iso2 == 'US' & Province_State == 'New York') # out of NY, unassigned 
  ids = colnames(tda) %>% as.Date(format = '%m/%d/%y') %>% is.na(.) %>% colnames(tda)[.] # those are not dates
  tda1 = tda %>% melt(data = ., id.vars = ids) %>% setnames('variable', 'date') %>% 
    mutate(value = value %>% as.numeric(), 
           date = date %>% as.Date(format = '%m/%d/%y')) %>% 
    group_by(Province_State, date) %>% summarise(value = sum(value, na.rm = T)) %>% data.table()
  # normalize to per 1 M pop
  tda1 = merge(tda1, pop.t, x.all = T, by = 'Province_State') %>% 
    mutate(value = value / Population * N)
  
  tda2 = tda1 %>% dcast(date ~ Province_State, value.var = 'value')
  # compute the increment for each date
  tda2 = data.table(date = tda2$date[-1], 
                    tda2[-1,2:ncol(tda2), with=F] - tda2[-nrow(tda2),2:ncol(tda2), with=F])
  
  
  # exclude incomplete week
  tda2 = tda2 %>% filter(date < Sys.Date() - (format(Sys.Date(), '%w') %>% as.numeric()))
  
  
  tmp = tda2 %>% .$date %>% MMWRweek %>% setnames(c('MMWRyear','MMWRweek','MMWRday'),c('year','week','day'))
  tda2 = cbind(tmp, tda2) %>% data.table()
  tda2 = tda2[, lapply(.SD, sum, na.rm=T), .SDcol = 5:ncol(tda2), by = c('year', 'week')]
  tda2 = data.table(date = MMWRweek2Date(MMWRyear = tda2$year, MMWRweek = tda2$week), tda2)
  
  matplot(tda2[,4:ncol(tda2), with=F], type = 'l')
  tda2 = tda2 %>% melt(., id.vars = c('date', 'year', 'week')) %>% setnames('variable', 'state')
  tda2[value < 0]  # some <0
  tda2[value < 0]$value = 0 # set it to 0
  
  tda2$data.type = da.type
  
  res = rbind(res, tda2)
  # save(tda2, file = paste0(dir_data, 'da_wkly.',da.type,'_us.RData'))
}
write.csv(res, paste0(dir_data, 'da_case_death_us_hopkins.csv'), row.names = F)

# cp Hopkins and NYT data
d_hop = read.csv(paste0(dir_data,'da_case_death_us_hopkins.csv')) %>% data.table() %>% 
  mutate(date = date %>% as.Date, source = 'Hopkins') 
d_nyt = read.csv(paste0(dir_data,'da_case_death_us_nyt.csv'))  %>% data.table() %>% 
  mutate(date = date %>% as.Date, source = 'NYT')
d_cdc = read.csv(paste0(dir_data,'da_case_death_us_cdc.csv'))  %>% data.table() %>% 
  mutate(date = date %>% as.Date, source = 'CDC')
d = rbind(d_hop, d_nyt, d_cdc) %>%
  mutate(source = factor(source, levels = c('Hopkins','NYT','CDC')))

dates.t = d$date %>% unique %>% as.Date

theme.t = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=8,angle = 45, hjust = 1),
                plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=8), legend.text=element_text(size=8),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                legend.key.size = unit(.2, 'cm'), #change legend key size
                legend.key.height = unit(.5, 'cm'), #change legend key height
                legend.key.width = unit(.2, 'cm')) #change legend key width)

pdf(paste0(dir_data, 'fig_cp_case_death_hopkins_nyt_cdc_',format(Sys.Date(),'%Y%m%d'),'.pdf'), width = 7, height = 4)
for(st in states){
  pp = ggplot(d %>% filter(state == st), aes(x = date, y = value, color = source)) + 
    geom_line() + labs(title = st) +
    facet_wrap(~data.type, scale = 'free_y')
  try(print(pp))
}
dev.off()

states.t = c('California', 'Florida', 'Iowa', 'Massachusetts', 'Michigan', 
             'New York', 'Pennsylvania', 'Texas', 'Washington', 'Wyoming')
pdf(paste0(dir_data, 'fig_cp_case_death_hopkins_nyt_cdc_10states_',format(Sys.Date(),'%Y%m%d'),'.pdf'), width = 10, height = 4)
for(st in states.t){
  pp = ggplot(d %>% filter(state == st), aes(x = date, y = value, color = source)) + 
    geom_line() + labs(title = st, y = 'Number per 1 million people', x = 'Week start') +
    facet_wrap(~data.type, scale = 'free_y') +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 8)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 8)],'%m/%d/%y')) +
    theme_minimal() +  theme.t
  try(print(pp))
}
dev.off()
# SIMILAR FOR MOST STATES BUT NYT DATA SEEMED SMOOTHED FOR THOSE WITH PROBLEMS

# Washington mortality data look off
st = 'Washington'
pdf(paste0(dir_data, 'fig_cp_case_death_hopkins_nyt_cdc_',st,'_',format(Sys.Date(),'%Y%m%d'),'.pdf'), width = 8, height = 5)
pp = ggplot(d %>% filter(state == st & data.type == 'death' & date >= dates.t[2]), aes(x = date, y = value, color = source)) + 
  geom_line() + labs(title = st, y = 'Deaths per 1 million people', x = 'Week start') +
  scale_x_date(breaks = dates.t[seq(2, length(dates.t), by = 4)],
               labels = format(dates.t[seq(2, length(dates.t), by = 4)],'%m/%d/%y')) +
  theme_minimal() +  theme.t
print(pp)
dev.off()

tmp = d %>% filter(data.type == 'death' & state == 'Washington') %>%
  dcast(., date ~ source, value.var = 'value')
which.max(tmp$CDC) - which.max(tmp$NYT)
tmp$date[which.max(tmp$CDC)]
tmp$date[which.max(tmp$NYT)]

# Generate the truth for forecast evaluation
for(variant.tag in c('non.Omicron', 'Omicron')){
  DAT.EPI = read.csv(paste0(dir_data, 'da_case_death_us_', data.tag, '_',variant.tag,'.csv')) %>% data.table()
  DAT.EPI$date = DAT.EPI$date %>% as.Date
  
  nwk.fcast = 26
  wk.fcast.start_vec = DAT.EPI %>% # filter(date > as.Date('2020/5/1') & date <= as.Date('2021/09/30')) %>% 
    .$date %>% unique %>% sort %>% as.character()
  
  states = DAT.EPI$state %>% unique
  truth = NULL; 
  for(st in states){
    for(d.type in c('case', 'death')){
      for(wk.fcast in wk.fcast.start_vec){
        wks.fcast.t = seq(as.Date(wk.fcast), length.out = nwk.fcast, by = 'week')
        da.t = DAT.EPI %>% filter(state == st & data.type == d.type & date %in% wks.fcast.t)
        if(nrow(da.t) < 1) next
        pw = ifelse(nrow(da.t) == nwk.fcast, da.t$date[which.max(da.t$value)], NA)
        pi = ifelse(nrow(da.t) == nwk.fcast, max(da.t$value), NA)
        tot = ifelse(nrow(da.t) == nwk.fcast, sum(da.t$value), NA)
        d.wkly = data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), 
                            target = paste0(1:nwk.fcast,'wk ahead'), value = c(da.t$value, rep(NA, nwk.fcast - nrow(da.t))))
        
        truth = rbind(truth, data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), target = 'peak week', value = pw),
                      data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), target = 'peak intensity', value = pi),
                      d.wkly, 
                      data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), target = 'total', value = tot)
        )
      }
    }
  }
  save(truth, file = paste0(dir_data, 'truth_',data.tag, '_',variant.tag,'.RData'))
  
}


# note re data issue:
# NYT data used for the forecast appeared off for more recent months for some states, likely due to infrequent/irregular reporting period
# as such, use CDC data for evaluation instead, which look more reasonable
# however, for Washington, mortality data from CDC prior to Dec 2022 look off (timing shifted by ~3 months) and NYT data look ok
# so for WA mortality, use NYT data for evaluation instead
# for the 3/30/23 CDC data release, FL case and deaths for the week of 3/26/23 were recorded as 0, due to data issues
# thus, need to exclude the last week for FL

# generate the truth based on a mixture of the most plausible data
# Generate the truth for forecast evaluation
for(variant.tag in c('non.Omicron', 'Omicron')){
  
  
  DAT.EPI_cdc = read.csv(paste0(dir_data, 'da_case_death_us_', 'cdc', '_',variant.tag,'.csv')) %>% data.table()
  DAT.EPI_cdc$date = DAT.EPI_cdc$date %>% as.Date
  DAT.EPI_nyt = read.csv(paste0(dir_data, 'da_case_death_us_', 'nyt' , '_',variant.tag,'.csv')) %>% data.table()
  DAT.EPI_nyt$date = DAT.EPI_nyt$date %>% as.Date
  
  date.t = as.Date('2022/12/1') # data after this week are about the same from both cdc and nyt
  DAT.EPI = rbind(DAT.EPI_cdc %>% filter(!(state == 'Washington' & data.type == 'death' & date <= date.t)),
                  DAT.EPI_nyt %>% filter((state == 'Washington' & data.type == 'death' & date <= date.t)))
  
  nwk.fcast = 26
  wk.fcast.start_vec = DAT.EPI %>% # filter(date > as.Date('2020/5/1') & date <= as.Date('2021/09/30')) %>% 
    .$date %>% unique %>% sort %>% as.character()
  
  states = DAT.EPI$state %>% unique
  truth = NULL; 
  for(st in states){
    for(d.type in c('case', 'death')){
      for(wk.fcast in wk.fcast.start_vec){
        wks.fcast.t = seq(as.Date(wk.fcast), length.out = nwk.fcast, by = 'week')
        da.t = DAT.EPI %>% filter(state == st & data.type == d.type & date %in% wks.fcast.t)
        if(nrow(da.t) < 1) next
        pw = ifelse(nrow(da.t) == nwk.fcast, da.t$date[which.max(da.t$value)], NA)
        pi = ifelse(nrow(da.t) == nwk.fcast, max(da.t$value), NA)
        tot = ifelse(nrow(da.t) == nwk.fcast, sum(da.t$value), NA)
        d.wkly = data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), 
                            target = paste0(1:nwk.fcast,'wk ahead'), value = c(da.t$value, rep(NA, nwk.fcast - nrow(da.t))))
        
        truth = rbind(truth, data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), target = 'peak week', value = pw),
                      data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), target = 'peak intensity', value = pi),
                      d.wkly, 
                      data.table(state = st, data.type = d.type, wk.fcast = as.Date(wk.fcast), target = 'total', value = tot)
        )
      }
    }
  }
  save(truth, file = paste0(dir_data, 'truth_','mixed', '_',variant.tag,'.RData'))
  
}

