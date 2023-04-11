# ARIMA forecast, for comparison with the forecast methods in this study
# 3/28/23

# for the first couple years, not enough data for seasonal arima, but allow the seasonal component for later years?
# note the last forecast was initiated at the end of the week of 2022-10-02
# last forecast period: 2022-10-02 to 2023-03-26
# note that the first forecast week is also in the training set, 
# b/c the forecast was set up to allow the use of retrospectively estimated paramesters (ifr and idr/alpha)

# to mimic the same retrospective forecast using arima, 
# allow the use of all available non-EPI data (x = mobility, + vx? probly not death/case) up to the week of 2022-10-02 
# can't use the EPI during the forecast period, otherwise, it's just fitting
# non-epi data for the ARIMAx model: 
# mobility for each week, 
# cumulative vax in the last 3 months for case / 9 months for deaths
# seasonality

dir_data = './data/'
dir_code = './scripts/'
dir_res = paste0('./results/rproj_arima/')

if(!file.exists(dir_res))
  dir.create(dir_res, recursive = T)

library(forecast); library(ggplot2)
source(paste0(dir_code,'getPlot.R'))
source(paste0(dir_code, 'get_fcastProbDist.R'))

ir = 1 # repeated runs
num_ens = 500
ci_levels = c(50, 80, 90, 95)
tag.proj = 'useEns'
sn.tag = 'mix.sn'
fcast.deflat = 1

models = c('arima','arimax.mob','arimax.sn','arimax.ms','arimax.full')
vars.arima = NULL
vars.arimax.mob = 'mob'
vars.arimax.sn = 'sn'
vars.arimax.ms = c('mob','sn')
vars.arimax.full = c('mob','sn', 'cum.v1', 'cum.v2')

fcast.start_vec_non.Omicron = c(seq(as.Date('2020/7/1'), as.Date('2021/8/15') , by = '1 week')
)
fcast.start_vec_Omicron =  seq(as.Date('2021/12/1'), as.Date('2022/9/30') , by = '1 week') 
fcast.start_vec = c(fcast.start_vec_non.Omicron, fcast.start_vec_Omicron)
length(fcast.start_vec)

da.tag = '_upto20221002' # time stamp for data used for the study
vtags = c('non.Omicron', 'Omicron')


# read mobility data
DAT.MOB = read.csv(paste0(dir_data, 'da_mobility_us.csv')) %>% data.table()
DAT.MOB = DAT.MOB[data.type == mob.type]
DAT.MOB$date = DAT.MOB$date %>% as.Date


# read vaccination data
if(F){
  DAT.VAC = read.csv(paste0(dir_data,'da_vx_perM_us_lagged_booster.csv')) %>% data.table()
  DAT.VAC$date = DAT.VAC$date %>% as.Date
  # compute the cumsum over the recent weeks
  fn_cumsum = function(da, var.t, date.t, dur = 90){
    da %>% filter(date <= date.t & date >= date.t - 90) %>% .[,var.t, with=F] %>% sum
  }
  
  
  DAT.VAC = DAT.VAC[,list(cum.v1.case = fn_cumsum(da = DAT.VAC,  var.t = 'n.v1', date.t = date, dur=90),
                          cum.v2.case = fn_cumsum(da = DAT.VAC,  var.t = 'n.v2', date.t = date, dur=90),
                          cum.v1.death = fn_cumsum(da = DAT.VAC,  var.t = 'n.v1', date.t = date, dur=270),
                          cum.v2.death = fn_cumsum(da = DAT.VAC,  var.t = 'n.v2', date.t = date, dur=270)
  ), by = c('state','date')]
  
  write.csv(DAT.VAC, paste0(dir_data,'DAT.VAC_uptoMarch2023.csv'), row.names = F)
}


DAT.VAC = read.csv(paste0(dir_data,'DAT.VAC_uptoMarch2023.csv')) %>% data.table()
DAT.VAC$date = DAT.VAC$date %>% as.Date

# read seasonal trend
DAT.SN = read.csv(paste0(dir_data, 'est.sn.2000t2020_us.csv')) %>% data.table()


tda = NULL
for(v.tag in c('non.Omicron','Omicron')){
  tmp = read.csv(paste0(dir_data, 'da_case_death_us_nyt_',v.tag, da.tag,'.csv')) %>% data.table()
  tda = rbind(tda, tmp, fill=T)
}
week.starts = tda$date %>% as.Date %>% unique
fcast.start_vec_non.Omicron = week.starts[week.starts >= as.Date('2020/7/1') & week.starts <= as.Date('2021/8/15')]
fcast.start_vec_Omicron = week.starts[week.starts >= as.Date('2021/12/1')] 
fcast.start_vec_Omicron = append(fcast.start_vec_Omicron, max(fcast.start_vec_Omicron)+7)
fcast.start_vec = c(fcast.start_vec_non.Omicron, fcast.start_vec_Omicron)
length(fcast.start_vec)

fcast.start.month_vec = format(fcast.start_vec, '%Y%m%d')

locs.ab.t = c('MA', 'NY', 'PA', 'MI', 'TX', 'IA', 'WY', 'CA', 'WA', 'FL')
locs = c(state.name, "District of Columbia")
names(locs) = c(state.abb, 'DC')
locs.t = locs[names(locs) %in% locs.ab.t]

mob.type = 'business'  # type of mobility data to use

nfcast = 26 # number weeks forecast done

date.start = as.Date('2020/3/1')
# date.end = tail(fcast.start_vec,1) + nfcast * 7
date.end = as.Date('2022/10/02') + nfcast * 7  # last forecast in the study
tab_dates = data.table(date = seq(date.start, date.end, by = 'week'))
tab_dates$year = tab_dates$date  %>% MMWRweek() %>% .$MMWRyear # 
tab_dates$week = tab_dates$date %>% MMWRweek() %>% .$MMWRweek

for(loc.t in locs.t){
  # extend the data to the end of the forecast period
  da.mob.t = DAT.MOB[state == loc.t & data.type == mob.type] %>% setnames('value', 'mob') %>%
    dplyr::select(date, mob) %>%
    right_join(tab_dates, by = 'date') %>%
    arrange(., date)
  
  mob_by_wk = da.mob.t %>% # dplyr::filter(date > as.Date('2020/4/1')) %>% # exclude pre-pandemic
    reshape2::dcast(., week ~ year, value.var = 'mob') # %>% dplyr::select(-`2020`) # exclude year 2020
  
  mob_by_wk$mx = da.mob.t %>% dplyr::filter(date > as.Date('2020/4/1')) %>% 
    reshape2::dcast(., week ~ year, value.var = 'mob') %>% 
    dplyr::select(-week) %>% apply(., 1, max, na.rm=T)
  
  
  # set weeks without data to the historical maximum 
  j0 = which(names(mob_by_wk)=='2022')
  for(j in j0:(ncol(mob_by_wk)-1)){
    for(i in 1:nrow(mob_by_wk)){
      if(is.na(mob_by_wk[i,j]))
        mob_by_wk[i,j] = mob_by_wk$mx[i]
    }
  }
  
  mob_by_wk$mx = NULL
  
  yr_wk_vec = tab_dates[,c('year','week')] %>% apply(1, paste, collapse='-')
  da.mob.t = mob_by_wk %>% reshape2::melt(., id.vars = 'week') %>%
      suppressMessages() %>%
      dplyr::mutate(date = MMWRweek2Date(MMWRyear = variable %>% as.character() %>% as.numeric, MMWRweek = week, MMWRday = 1))
    # exclude dates (i.e., week 53) that do not exist
  da.mob.t$yr_wk = da.mob.t[,c('variable','week')] %>% apply(1, paste, collapse='-') %>% gsub(' ','',.)
  da.mob.t = da.mob.t %>% 
    filter(yr_wk %in% yr_wk_vec) %>%
      dplyr::select(date, value) %>%
      filter(!is.na(value)) %>% 
      setnames(., 'value', 'mob') %>%
      arrange(., date)
  
  # vaccination
  da.vx.t = DAT.VAC[state == loc.t]
  # done in the end of March 2023, so include the full data
  da.vx.t = da.vx.t %>% right_join(., y = tab_dates, by = 'date') %>% 
    mutate(state = NULL, year = NULL, week = NULL) %>% arrange(., date)
  if(!is.na(da.vx.t$cum.v1.case[nrow(da.vx.t)])){
    # include all the needed data
    # simply replace the period w/o vx with 0
    da.vx.t = da.vx.t %>% replace_na(., list(cum.v1.case = 0, cum.v2.case = 0, cum.v1.death = 0, cum.v2.death =0)) %>%
      arrange(., date)
  }
    
 
  for(ir in 1:10){  # 
    
    for(iwk_fcast in 1:103){
      
      print(paste(loc.t, 'run', ir, 'wk', iwk_fcast))
      
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
      ff =  paste0(dir_res, gsub(' ','',loc.t),'_',variant.tag,'_train.proj_',tag.proj,'_df',fcast.deflat,'_',sn.tag,'_',
                   # fcast.start.month, # # this may cause redundancy
                   format(fcast.start,'%Y%m%d'),
                   '_r',ir,'.RData')
      tmp = try(load(ff), silent = T) %>% suppressMessages()
      if(!any(class(tmp)=="try-error"))
        next 
      
      # da = da[date < as.Date(date.stop)]
      # set population size to 20k? make it smaller for stochasticity
      vax.start = as.Date('2020/12/14')
      
      if(F){
        if(variant.tag == 'non.Omicron'){
          date.start = as.Date('2020/03/01') # start from 3/1/20
        } else {
          date.start = as.Date('2021/11/21') # for the Omicron period
        }
      }
      
      DAT.EPI = read.csv(paste0(dir_data, 'da_case_death_us_nyt_',variant.tag,da.tag,'.csv')) %>% data.table()
      DAT.EPI$date = DAT.EPI$date %>% as.Date
      # da = da[date >= date.start]
      
      week.starts = DAT.EPI[date >= date.start]$date %>% unique %>% as.Date %>% sort
      
      # fcast.start = week.starts[which.min(abs(week.starts - fcast.start))]
   
      date.end.epi = fcast.start - 7
      date.end.nonepi = fcast.start + (nfcast - 1) * 7
      iwk_fcast
      fcast.start.month_vec[iwk_fcast]
      fcast.start
      fcast.wk.starts = as.Date(fcast.start) + seq(0, length.out = nfcast, by = 7)
      wk.starts = seq(date.start, tail(fcast.wk.starts,1), by = 'week')
      
      fcastDist = NULL
      fcast_stats = NULL
      states_stats = NULL
      for(m.t in models){
        
        model.fit.case = model.fit.death = fcast.case = fcast.death = NULL # set it to null to begin with
        
        for(mea.t in c('case', 'death')){
          
          da.full.t = DAT.EPI[state == loc.t] %>% 
            dcast(., date + year + week ~ data.type, value.var = 'value')
          da.t = da.full.t %>%
            filter(date >= date.start & date <= date.end.epi) %>%
            arrange(., date)
          
          # mobility - on the top
          # seasonality
          da.sn.t = DAT.SN[state == loc.t] %>% .[order(week)]
          # match to dates
          dates = seq(as.Date(date.start), tail(fcast.wk.starts,1), by = 'week')
          tmp = MMWRweek(dates) %>% setnames(., new = c('year','week','day')) %>%
            mutate(date = dates)
          da.sn.t = da.sn.t %>% right_join(., y = tmp, by = 'week') %>% dplyr::select(date, value) %>%
            setnames(., 'value','sn')
          
          
          # data for the training period
          da.t = da.t %>% left_join(., da.mob.t, by = 'date') %>%
            left_join(., da.sn.t, by = 'date') %>%
            left_join(., da.vx.t, by = 'date') 
          
          # like the filter based forecasts, no forecasts if fewer than 5 weeks of training data
          if(nrow(da.t) < 5){
            print('no forecasts, b/c fewer than 5 weeks of training data')
            break
          }
          # data for the fcast period
          fda.t = da.mob.t %>% 
            left_join(., da.sn.t, by = 'date') %>%
            left_join(., da.vx.t, by = 'date') %>% 
            filter(date %in% fcast.wk.starts) %>%
            data.table()
          
          da.trn.t = da.t[,mea.t, with=F]
          var.vx.t = colnames(da.t)[grepl('cum.v', colnames(da.t))]
          var.vx.t = var.vx.t[grepl(mea.t, var.vx.t)]
          da.t %>% setnames(., var.vx.t, c('cum.v1','cum.v2'))
          fda.t %>% setnames(., var.vx.t, c('cum.v1','cum.v2'))
          vars.t = get(paste0('vars.', m.t))
          if(is.null(vars.t)){
            xreg.trn.t = NULL
            xreg.fcast.t = NULL
          } else {
            xreg.trn.t = da.t[,c(vars.t),with=F] %>% as.matrix()
            xreg.fcast.t = fda.t[,c(vars.t),with=F] %>% as.matrix() # %>% dplyr::select(mob, sn) 
            
          }
          
          tmp = try({
            model.t = da.trn.t %>%
              auto.arima(., xreg = xreg.trn.t)
          })
          # also check likelihood - when death # is 0, can't find a model
          # is.infinite(model.t$loglik)
          
          if(!any(class(tmp) == "try-error")){
            loglik.t = model.t$loglik
          } else {
            loglik.t = NA
          }
          
          if(any(class(tmp) == "try-error") | is.infinite(loglik.t)){
            # also check likelihood - when death # is 0, can't find a model
            # is.infinite(model.t$loglik)
            next
          } else {
            ci.t = sqrt(model.t$sigma2) %>% round(., 4)
            fitted.t = fitted(model.t) %>% as.numeric %>% unname() %>% round(., 4)
            model.fit.t = data.table(Week.start = da.t$date, state = mea.t, 
                                     mean =  fitted.t, median =  fitted.t,
                                     iqr.lwr =  (fitted.t - 2/3*ci.t) %>% round(., 4),
                                     iqr.upr =  (fitted.t + 2/3*ci.t) %>% round(., 4),
                                     ci95.lwr =  (fitted.t - 1.96*ci.t) %>% round(., 4),
                                     ci95.upr =  (fitted.t + 1.96*ci.t) %>% round(., 4)
            )
            model.fit.t[model.fit.t<0] = 0
            
            fcast.t = model.t %>%
              forecast(h=nfcast, xreg = xreg.fcast.t, level = ci_levels)
            ens_fcast.t = matrix(0, nfcast, num_ens)
            for (i in 1:num_ens)
              ens_fcast.t[, i] <- simulate(model.t, nsim = nfcast, bootstrap = TRUE, xreg = xreg.fcast.t)
            ens_fcast.t[ens_fcast.t < 0] = 0 # at least >=0
            
            
            fcast_stats.t = cbind(ens_fcast.t %>% apply(., 1, median), fcast.t$lower, fcast.t$upper) %>% data.table()
            colnames(fcast_stats.t) = c('median', paste0('ci',ci_levels,'.lwr'), paste0('ci',ci_levels,'.upr')) %>% gsub('ci50','iqr',.)
            fcast_stats.t[fcast_stats.t<0] = 0
            fcast_stats.t$measure = mea.t
            fcast_stats.t$Week.start = fcast.wk.starts
            
            
            assign(paste0('fcast.',mea.t), ens_fcast.t)
            assign(paste0('fcast.stat.',mea.t), fcast_stats.t)
            assign(paste0('model.fit.',mea.t), model.fit.t)
          }
          
          
        } # measure
        
        if(!is.null(fcast.case) & !is.null(fcast.death)){
          fcastDist.t = fn_getProbDist(fcast.case, fcast.death, bins.case, bins.death)
          # organize the output similar to the other ensemble forecast
          fcast_stats.t = rbind(fcast.stat.case, fcast.stat.death)
          states_stats.t = rbind(model.fit.case, model.fit.death)
          
          fcastDist.t$fcast.type = m.t # model
          fcast_stats.t$fcast.type = m.t # model
          states_stats.t$fcast.type = m.t
          # sn.tag = ifelse(grepl('sn', vars.t), 'w.sn', 'no.sn')
          # if(length(sn.tag)==0)
          #   sn.tag = 'no.sn'
          
          fcastDist = rbind(fcastDist, fcastDist.t)
          fcast_stats = rbind(fcast_stats, fcast_stats.t)
          states_stats = rbind(states_stats, states_stats.t)
          
        } else {
          print(paste('no model found for', m.t))
        }
        
      }  # diff models
      
      
      train.proj = list(states_stats = states_stats,
                        fcast_stats = fcast_stats, 
                        fcastDist = fcastDist)
      train.proj$fcast.start.week = fcast.start
      train.proj$fcast.start.month = fcast.start.month
      train.proj$loc = loc.t
      train.proj$fcast.deflat = fcast.deflat
      train.proj$scenario = 'asIs'
      train.proj$variant = variant.tag
      train.proj$seasonality = sn.tag
      
      save(train.proj, file = paste0(dir_res, gsub(' ','',loc.t),'_',variant.tag,'_train.proj_',tag.proj,'_df',fcast.deflat,'_',sn.tag,'_',
                                     # fcast.start.month, # # this may cause redundancy
                                     format(fcast.start,'%Y%m%d'),
                                     '_r',ir,'.RData'))
    } # wk
  } # repeated runs
  
}  # loc 


