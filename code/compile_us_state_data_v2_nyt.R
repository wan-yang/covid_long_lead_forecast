# compile covid data
# 8/3/22 - use NYT data for case/death instead of hopkins

data.tag = 'nyt'

dir_code = './code/'
dir_data = './data/'

source(paste0(dir_code, 'Fn_util.R'))

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

url.t = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv'
# case = read_csv(url(url.t)) %>% data.table() 

url.t = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv'
# death= read_csv(url(url.t)) %>% data.table() 

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

# mobility 
mob.type.bus = c('retail_and_recreation_percent_change_from_baseline',
                 'grocery_and_pharmacy_percent_change_from_baseline',
                 'transit_stations_percent_change_from_baseline',
                 'workplaces_percent_change_from_baseline')  
mob.type.bus0 = c('retail_and_recreation_percent_change_from_baseline',
                 'transit_stations_percent_change_from_baseline',
                 'workplaces_percent_change_from_baseline')
# ,'parks_percent_change_from_baseline',
# 'residential_percent_change_from_baseline'
mob.type.full = c('retail_and_recreation_percent_change_from_baseline',
                  'grocery_and_pharmacy_percent_change_from_baseline',
                  'parks_percent_change_from_baseline',
                  'transit_stations_percent_change_from_baseline',
                  'workplaces_percent_change_from_baseline',
                  'residential_percent_change_from_baseline')

url.mob = 'https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv'
d.mob = read_csv(url(url.mob)) %>% data.table() 

mob = d.mob[country_region %in% c("United States")]
rm(d.mob)
mob = mob[, lapply(.SD, median, na.rm=T), by = c('date', "sub_region_1"), 
            .SD = mob.type.full]
mob$mob.bus = 1 + rowMeans(mob[,mob.type.bus,with=F])/ 100  # take the average for the 
mob$mob.full = 1 + rowMeans(mob[,mob.type.full,with=F])/ 100 
mob$mob.bus0 = 1 + rowMeans(mob[,mob.type.bus0,with=F])/ 100  # take the average for the 

matplot(mob[sub_region_1=="Florida",c('mob.bus0','mob.bus','mob.full')], type='l', lty =1, col = c('orange','blue','black'))
abline(h=0); legend('bottomleft',c('business','all'), lty=1, col = c('blue','black'), bty='n')

# make it weekly
d.mob.bus = dcast(mob, date ~ sub_region_1, value.var = "mob.bus") %>% setnames('NA', 'US')
d.mob.full = dcast(mob, date ~ sub_region_1, value.var = "mob.full") %>% setnames('NA', 'US')

tmp = MMWRweek(d.mob.bus$date) %>% setnames(paste0('MMWR',c('year','week','day')),c('year','week','day'))
d.mob.bus = cbind(tmp, d.mob.bus) %>% data.table()
d.mob.bus = d.mob.bus[,lapply(.SD, mean, na.rm=T), by = c('year', 'week'), 
              .SDcols = 5:ncol(d.mob.bus)]
d.mob.bus$date = MMWRweek2Date(MMWRyear = d.mob.bus$year, MMWRweek = d.mob.bus$week, MMWRday = 1)

tmp = MMWRweek(d.mob.full$date) %>% setnames(paste0('MMWR',c('year','week','day')),c('year','week','day'))
d.mob.full = cbind(tmp, d.mob.full) %>% data.table()
d.mob.full = d.mob.full[,lapply(.SD, mean, na.rm=T), by = c('year', 'week'), 
                      .SDcols = 5:ncol(d.mob.full)]
d.mob.full$date = MMWRweek2Date(MMWRyear = d.mob.full$year, MMWRweek = d.mob.full$week, MMWRday = 1)

# check NA
Fn_checkNA = function(x){any(is.na(x))}
loc.na = d.mob.bus %>% apply(2, Fn_checkNA)
loc.na = loc.na[loc.na] %>% names
d.mob.bus[,loc.na] = fn_segIntrpl(d.mob.bus, loc.na)

loc.na = d.mob.full %>% apply(2, Fn_checkNA)
loc.na = loc.na[loc.na] %>% names
d.mob.full[,loc.na] = fn_segIntrpl(d.mob.full, loc.na)

d.mob.bus = melt(d.mob.bus, id.vars = c("date", "year","week")) %>% setnames('variable', 'state')
d.mob.full = melt(d.mob.full, id.vars = c("date", "year","week")) %>% setnames('variable', 'state')
d.mob.bus$data.type = 'business'
d.mob.full$data.type = 'all'
d.mob = rbind(d.mob.bus, d.mob.full)
write.csv(d.mob, paste0(dir_data, 'da_mobility_us.csv'), row.names = F)


N = 1e6

# aggregate to weekly intervals
# states = case$Province_State %>% unique %>% .[!(. %in% c("American Samoa", "Diamond Princess","Grand Princess"))]
# states = pop$Province_State  %>% unique
states = da_nyt$state %>% unique
pop.t = pop %>% filter(iso2 == 'US' & Province_State %in% states & is.na(Admin2)) %>% dplyr::select(Province_State, Population) %>% data.table() #  %>% sum(na.rm = T)

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






# Florida data are odd with 0 cases for some weeks
zeros = res %>% filter(date > as.Date('2020/6/1') & value == 0 & data.type == 'case')
# 2022-03-13
tmp2 = res %>% filter(date > as.Date('2021/11/1') & state == 'Florida' & data.type == 'case')
plot(tmp2$date, tmp2$value, type = 'l')
tmp3 = res %>% filter(date > as.Date('2021/11/1') & state == 'Florida' & data.type == 'death')
plot(tmp3$date, tmp3$value, type = 'l')

# NYT data look more smoothed




# look at vaccination data
# https://github.com/owid/covid-19-data/tree/master/public/data/vaccinations
# people_vaccinated: total number of people who received at least one vaccine dose. If a person receives the first dose of a 2-dose vaccine, this metric goes up by 1. If they receive the second dose, the metric stays the same.
# people_vaccinated_per_hundred: people_vaccinated per 100 people in the total population of the state.
# people_fully_vaccinated: total number of people who received all doses prescribed by the vaccination protocol. If a person receives the first dose of a 2-dose vaccine, this metric stays the same. If they receive the second dose, the metric goes up by 1.
# people_fully_vaccinated_per_hundred: people_fully_vaccinated per 100 people in the total population of the state.
# booster:
# Aug 13, 2021 (Reuters) - U.S regulators authorized a third dose of COVID-19 vaccines by Pfizer Inc (PFE.N)-BioNTech and Moderna Inc (MRNA.O) on Friday for people with compromised immune systems who are likely to have weaker protection from the two-dose regimens.
# https://www.reuters.com/world/middle-east/us-fda-authorizes-covid-19-vaccine-boosters-immunocompromised-2021-08-13/
# 9/17/21 F.D.A. Advisory Panel Recommends Pfizer Boosters for Older People and Others at High Risk  https://www.nytimes.com/2021/09/17/us/politics/fda-pfizer-booster-covid.html
# 9/24/21 CDC Chief Clears Boosters for Millions in U.S., Overrules Panel
# https://news.bloomberglaw.com/coronavirus/pfizer-boosters-weighed-by-cdc-after-fda-approval-for-seniors
# it looks like in NYC, boosters were administered to some people starting Aug, more so in Sep 2021

colnames(vx)
vx[1000] 
# vx.t = vx %>% filter(location == 'New York State') %>% dplyr::select(date, people_vaccinated_per_hundred, people_fully_vaccinated_per_hundred)
fn_getDailyVx = function(vx.t, vx.start = as.Date('2020/12/14'), boost.start = as.Date('2021/8/13')){
  # mising data for earliest weeks
  vx.t$date = vx.t$date %>% as.Date
  vx.pad = data.table(date = seq(as.Date(vx.start)-1, min(vx.t$date)-1, by = 'day'), people_vaccinated_per_hundred = NA_real_, 
                      people_fully_vaccinated_per_hundred = NA_real_,
                      total_boosters_per_hundred = NA_real_)
  vx.pad[1]$people_vaccinated_per_hundred = 0
  vx.pad[1:(3*7)]$people_fully_vaccinated_per_hundred = 0
  vx.pad$total_boosters_per_hundred = 0
  
  vx.t = rbind(vx.pad, vx.t)
  # for booster - first ~7 months should be all 0
  vx.t[date < as.Date(boost.start)]$total_boosters_per_hundred = 0
  
  vx.t$idx = 1:nrow(vx.t)
  vx.t$daily.vx.t1 = 0
  vx.t$daily.vx.t2 = 0
  vx.t$daily.vx.t3 = 0
  
  # do the two dose separately, plus boosters
  for(dv in c('people_vaccinated_per_hundred', 'people_fully_vaccinated_per_hundred', 'total_boosters_per_hundred')){
    
    idx.all = which(is.na(vx.t[,dv,with=F]))
    
    # idx.all = which(is.na(vx.t$people_vaccinated_per_hundred))
    
    consecutive =  idx.all[-1] - idx.all[-length(idx.all)]
    i.div = which(consecutive >=2) # 1/27/22 check here if err, change '>' to '>='
    w.div = idx.all[which(consecutive >=2)] # 1/27/22 check here if err, change '>' to '>='
    grps = list()
    if(length(i.div) == 1){
      grps[[1]] = idx.all[1]: w.div[1]
      grps[[2]] =idx.all[i.div[1]+1]: tail(idx.all,1)
    } else {
      for(id in 1: length(i.div)){
        if(id == 1){
          grps[[id]] = idx.all[1]: w.div[id]
        } else if (id == length(i.div)){
          # both before and after
          
          grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]]
          grps[[id+1]] = idx.all[i.div[id]+1]: tail(idx.all,1)
        } else {
          grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]] # (idx.all[i.div[id-1]+1]) : w.div[id]
        }
        
      }
    }
    for(ig in 1: length(grps)){
      idx0 = grps[[ig]]
      
      n.add = ifelse(length(idx0) < 2, 1, 1)
      idx = seq((idx0[1] - n.add) %>% pmax(1), 
                (tail(idx0,1)+ n.add) %>% pmin(nrow(vx.t)),
                by = 1
      ) # add an additional week in case there is back adjustment
      
      vx.t0 = vx.t[idx]; 
      vx.t0 = vx.t0[complete.cases(vx.t0[,c('date', dv), with =F])]
      
      if(dv == 'people_vaccinated_per_hundred'){
        fit1 = lm(people_vaccinated_per_hundred ~ idx, data = vx.t0)
        vx.t$daily.vx.t1[idx0] = predict(fit1, newdata = vx.t[idx0])
      } else if (dv == 'people_fully_vaccinated_per_hundred'){
        fit2 = lm(people_fully_vaccinated_per_hundred ~ idx, data = vx.t0)
        vx.t$daily.vx.t2[idx0] = predict(fit2, newdata = vx.t[idx0])
      } else if (dv == 'total_boosters_per_hundred'){
        fit3 = lm(total_boosters_per_hundred ~ idx, data = vx.t0)
        vx.t$daily.vx.t3[idx0] = predict(fit3, newdata = vx.t[idx0])
      }
      
    }
  } # do the two dose separately
  
 
  vx.t[is.na(people_vaccinated_per_hundred)]$people_vaccinated_per_hundred = vx.t[is.na(people_vaccinated_per_hundred)]$daily.vx.t1
  vx.t[is.na(people_fully_vaccinated_per_hundred)]$people_fully_vaccinated_per_hundred = vx.t[is.na(people_fully_vaccinated_per_hundred)]$daily.vx.t2
  vx.t[is.na(total_boosters_per_hundred)]$total_boosters_per_hundred = vx.t[is.na(total_boosters_per_hundred)]$daily.vx.t3
  
  vx.t$n.v1 = c(vx.t$people_vaccinated_per_hundred[1], (vx.t$people_vaccinated_per_hundred[-1] - vx.t$people_vaccinated_per_hundred[-nrow(vx.t)])) / 100 * N  # convert to per Million
  vx.t$n.v2 =c(vx.t$people_fully_vaccinated_per_hundred[1], (vx.t$people_fully_vaccinated_per_hundred[-1] - vx.t$people_fully_vaccinated_per_hundred[-nrow(vx.t)])) / 100 * N
  vx.t$n.v3 =c(vx.t$total_boosters_per_hundred[1], (vx.t$total_boosters_per_hundred[-1] - vx.t$total_boosters_per_hundred[-nrow(vx.t)])) / 100 * N
  
  vx.t[n.v1 < 0]$n.v1 = 0
  vx.t[n.v2 < 0]$n.v2 = 0
  vx.t[n.v3 < 0]$n.v3 = 0
  
  vx.t[,c('date', 'n.v1', 'n.v2', 'n.v3'), with = F]
}

vx[location == "New York State"]$location = 'New York'
res = NULL
for(st in states){
  vx.t = vx %>% filter(location == st) %>% dplyr:: select(date, people_vaccinated_per_hundred, people_fully_vaccinated_per_hundred, total_boosters_per_hundred)
  vx.t = vx.t %>% fn_getDailyVx
  vx.t$state = st
  res = rbind(res, vx.t)
}
write.csv(res, paste0(dir_data, 'da_vx_perM_us_nolag.csv'), row.names = F)

# lag it
lagV1 = 14 # 
lagV2 = 7
lagV3 = 7

res1daily = res[,c('state', 'date','n.v1'), with=F]
res1daily = res1daily[order(date)]
res1daily$date = res1daily$date + lagV1
res2daily = res[,c('state', 'date','n.v2'), with=F]
res2daily = res2daily[order(date)]
res2daily$date = res2daily$date + lagV2
res3daily = res[,c('state', 'date','n.v3'), with=F]
res3daily = res3daily[order(date)]
res3daily$date = res3daily$date + lagV3
res = merge(res1daily, res2daily, all=T,by= c('state', 'date'))
res = merge(res, res3daily, all=T,by= c('state', 'date'))
res[is.na(res)] = 0
write.csv(res, paste0(dir_data,'da_vx_perM_us_lagged.csv'), row.names = F)


# for the omicron variant, VE is much lower so need booster
# to simplify, for period before omicron, treat 1st/2nd dose as is, and 
# after omicron became dominant, treat 2nd dose as 1st (i.e. combine both) and 3rd dose as 2nd
# use variant data from CoVariants (GISAID) to find the date omicron became predominant in each state
d.variant = read.csv(paste0(dir_data, 'perc.variant_by_state.csv')) %>% data.table()
d.variant[location == 'Washington DC']$location = 'District of Columbia'
res2 = NULL
for(st in states){
  da.t = res[state == st]
  vda.t = d.variant[location == st]
  date.omicron = (vda.t$week[which(vda.t$Omicron > .6)] %>% head(1) %>% as.Date)  # use the end, b/c GISAID likely bias toward more new variant sequences
  
  
  if(length(date.omicron)>0){
    
    # now for those before date.omicron: b/c most got their 3nd does by that time protection had waned
    # add n.v3 to n.v2
    d1 = da.t[date < date.omicron,]
    d1$n.v2 = d1$n.v2 + d1$n.v3
    d1$n.v3 = NULL
    
    # for those after date.omicron, treat 2nd dose as 1st (i.e. combine both) and 3rd dose as 2nd
    # if it is more about recency of vaccination, maybe 2nd and 3rd dose when first vax afford similar effect and should be added together?
    
    d2 = da.t[date >= date.omicron]
    setnames(d2, c('n.v1','n.v2','n.v3'), c('n.v0','n.v1','n.v2'))
    d2$n.v0 = NULL
    da.t = rbind(d1, d2)
  } else {
    # break;
    da.t$n.v3 = NULL
  }
  
  res2 = rbind(res2, da.t)
  
}
write.csv(res2, paste0(dir_data,'da_vx_perM_us_lagged_booster.csv'), row.names = F)


# get the dates different variants became predominant by state
d.variant = read.csv(paste0(dir_data, 'perc.variant_by_state.csv')) %>% data.table()
d.variant[location == 'Washington DC']$location = 'District of Columbia'
key.variants = c('Delta','Omicron')
vdates = NULL
for(st in states){
  vda.t = d.variant[location == st]
  
  for(v in key.variants){
    date.t1 = (vda.t$week[which(vda.t[,v,with=F] > .6)] %>% head(1) %>% as.Date)  # use the end, b/c GISAID likely bias toward more new variant sequences
    
    # for some reason, New Jessey detected Omicron during the weeks starting 5/17/21 - could it be an error?
    # to ensure it's real, start with the one with consecutive >0
    # date.t0 = (vda.t$week[which(vda.t[,v,with=F] > 0)] %>% head(1) %>% as.Date) # when first reported
    wk.t = vda.t$week[which(vda.t[,v,with=F] > 0)]  %>% as.Date
    date.t0 = wk.t[which(wk.t[-1] - wk.t[-length(wk.t)] == 14)] %>% head(1)
    
    vdates = rbind(vdates, data.table(state = st, variant = v, 
                                      date.1st = date.t0 - 1,
                                      date.dominant = date.t1 -1)) # shift by 1 day to make it start from Sunday
  }
}
write.csv(vdates, paste0(dir_data, 'dates_variant_by_state.csv'), row.names = F)  

# compute variant specific case and death
da = read.csv(paste0(dir_data, 'da_case_death_us_nyt.csv'))
da$date = da$date %>% as.Date

# check the time lag from case to death
tda = da %>% reshape2::dcast(., state + date ~ data.type, value.var = 'value')
lags = NULL
for(state.t in (tda$state %>% unique)){
  ttda1 = tda %>% filter(state == state.t & date >= as.Date('2020/3/1') & date < as.Date('2021/12/1'))
  tmp1 = ccf(ttda1$death, ttda1$case, lag.max = 8, plot = F) 
  ttda2 = tda %>% filter(state == state.t & date >= as.Date('2021/12/1'))
  tmp2 = ccf(ttda2$death, ttda2$case, lag.max = 8, plot = F) 
  lags = rbind(lags, 
               data.table(state = state.t, 
                          lag.max.beforeOmicron = tmp1$lag[tmp1$acf %>% which.max()],
                          lag.max.Omicron = tmp2$lag[tmp2$acf %>% which.max()]
                          )
               )
}
lags %>% .$lag.max.beforeOmicron %>% summary  # mean = 2.51, median = 2 weeks
lags %>% .$lag.max.Omicron %>% summary # mean = 2.53, median = 2 weeks


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

write.csv(res, paste0(dir_data,'da_case_death_us_nyt_', variant.t,'.csv'), row.names = F)
write.csv(res2, paste0(dir_data,'da_case_death_us_nyt_non.', variant.t,'.csv'), row.names = F)


# cp Hopkins and NYT data
d_hop = read.csv(paste0(dir_data,'da_case_death_us.csv')) %>% data.table() %>% 
  mutate(date = date %>% as.Date, source = 'Hopkins') 
d_nyt = read.csv(paste0(dir_data,'da_case_death_us_nyt.csv'))  %>% data.table() %>% 
  mutate(date = date %>% as.Date, source = 'NYT')
d = rbind(d_hop, d_nyt)

pdf(paste0(dir_data, 'fig_cp_case_death_hopkins_v_nyt.pdf'), width = 7, height = 4)
for(st in states){
  pp = ggplot(d %>% filter(state == st), aes(x = date, y = value, color = source)) + 
    geom_line() + labs(title = st) +
    facet_wrap(~data.type, scale = 'free_y')
  try(print(pp))
}
dev.off()

# SIMILAR FOR MOST STATES BUT NYT DATA SEEMED SMOOTHED FOR THOSE WITH PROBLEMS

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

