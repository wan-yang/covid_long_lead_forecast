# combine results from multiple runs - diff eval method


state.names = c('case','death')
sn_levels = c('w.sn', 'w.sn.tran', 'no.sn', 'mix.sn',
              paste0('w.sn.tran',1:8))
sn_labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality', 'Mixed seasonality',
              paste0('w.sn.tran',1:8))

key_measures = c('Cases','Hospitalizations','Deaths') 

# get the training estimates
for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_out, 'res.train_', loc.t,date.tag,'.RData')))
  if(class(tmp) != "try-error"){
    rm(res.train)
    next
  }
  
  res.train <-  lapply(files.t, function(x) {
    
    # print(x)
    try(load(x))
    
    id <- gsub(dir_res,'',x) %>% 
      strsplit('_') %>% unlist
    
    d <- train.proj$states_stats
    if(is.null(d)){
      d = NULL
    } else {
      d <- train.proj$states_stats %>%
        melt(id.vars = c('fcast.type','state','Week.start'), variable.factor = F)
      d$fcast.start.month = train.proj$fcast.start.month
      d$fcast.start.week = train.proj$fcast.start.week #  train.proj$fcast_stats$Week.start %>% min %>% as.Date
      d$loc = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
      
      d$seasonality = id[grepl('sn',id)] # id[5]
      d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
      # d$fcast.type = train.proj$fcast.type  # id[5]
      d$fcast.deflat = train.proj$fcast.deflat
      d$variant = train.proj$variant
    }
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, j = list(value = round(mean(value), 6)), by = list(variant, fcast.type, fcast.deflat, loc, seasonality, state, variable, fcast.start.month, fcast.start.week, Week.start)]) # %>%
  
  
  res.train$state = factor(res.train$state, levels=state.names)
  res.train$seasonality = factor(res.train$seasonality, levels = sn_levels, 
                                 labels = sn_labels)
  
  save(res.train, file = paste0(dir_out, 'res.train_', loc.t, date.tag,'.RData'))
  rm(res.train)
  print(paste(loc.t, 'res.train done'))
}

# get the forecasts
for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_out, 'res.proj_', loc.t,date.tag,'.RData')))
  if(class(tmp) != "try-error"){
    rm(res.proj)
    next
  }
  
  res.proj <-  lapply(files.t, function(x) {
    
    try(load(x))
    
    id <- gsub(dir_res,'',x) %>% 
      strsplit('_') %>% unlist
    
    d <- train.proj$fcast_stats
    if(is.null(d)){
      d = NULL
    } else {
      if(! is.null(train.proj$fcast_stats$fcast.type)){  # this is diff arima model used
        d <- train.proj$fcast_stats %>%
          melt(id.vars = c('fcast.type', 'measure','Week.start'), variable.factor = F)
      } else {
        d <- train.proj$fcast_stats %>%
          melt(id.vars = c('measure','Week.start'), variable.factor = F)
      }
      
      
      d$fcast.start.month = train.proj$fcast.start.month
      d$fcast.start.week = train.proj$fcast.start.week #  d$Week.start %>% min %>% as.Date
      d$loc = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
      d$seasonality = id[grepl('sn',id)] # id[5]
      d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
      d$scenario = train.proj$scenario  # id[5]
      d$fcast.deflat = train.proj$fcast.deflat
      d$variant = train.proj$variant
    }
      
    
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, j = list(value = round(mean(value), 0)), by = list(variant, fcast.type, fcast.deflat, loc, seasonality, scenario, measure, variable, fcast.start.month, fcast.start.week, Week.start)]) # %>%
  
  # res.proj$seasonality = factor(res.proj$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn'), labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality'))
  res.proj$seasonality = factor(res.proj$seasonality, levels = sn_levels, 
                                 labels = sn_labels)
  
  save(res.proj, file = paste0(dir_out, 'res.proj_', loc.t,date.tag,'.RData'))
  rm(res.proj)
  print(paste(loc.t,'res.proj done'))
  
}

# evaluation of the projections
fn_abs.err = function(da, wk.horizon.t, mea.t, sce.t, obs.name, est.name){  # relative absolute error
  # not good, could div by 0
  # mean(unlist(abs(da[,est.name,with=F] - da[,obs.name,with=F]) / da[,obs.name,with=F])) %>% round(3)
  (mean(unlist(abs(da[wk.horizon == wk.horizon.t & measure == mea.t & scenario == sce.t, est.name,with=F] - da[wk.horizon == wk.horizon.t & measure == mea.t & scenario == sce.t,obs.name,with=F]))) / mean(unlist(da[wk.horizon == wk.horizon.t & measure == mea.t & scenario == sce.t,obs.name,with=F]))) %>% round(3)
}
fn_rrmse = function(da, wk.horizon.t, mea.t, sce.t, obs.name, est.name){
  (sqrt(mean(unlist((da[wk.horizon == wk.horizon.t & measure == mea.t & scenario == sce.t,est.name,with=F] - da[wk.horizon == wk.horizon.t & measure == mea.t & scenario == sce.t,obs.name,with=F])^2))) / mean(unlist(da[wk.horizon == wk.horizon.t & measure == mea.t & scenario == sce.t,obs.name,with=F])))  %>% round(3)
}
fn_corr = function(da, wk.horizon.t, mea.t, sce.t, obs.name, est.name){  # relative absolute error
  cor(da[wk.horizon %in% 1:wk.horizon.t & measure == mea.t & scenario == sce.t,est.name,with=F], da[wk.horizon %in% 1:wk.horizon.t & measure == mea.t & scenario == sce.t,obs.name,with=F]) %>% as.numeric %>% round(3)
}
fn_acc.cover = function(da){
  # accuracy base on whether the lower and upper bound cover the obs
  fn_covered = function(x){
    ifelse(x[1] >= x[2] & x[1] <= x[3], 1, 0)
  }
  da %>% apply(1, fn_covered) # %>% mean  %>% round(3)
}

cut0 = .1 # count those with abs rel error within 10% of the observed as accurate
cut1 = .25 # count those with abs rel error within 25% 10% of the observed as accurate
cut2 = .5 # count those with abs rel error within 50% 10% of the observed as accurate

fn_eval_arima = function(d.proj, d.obs){
  # evaluate the accuracy of the projection, cp observation
  # criteria:
  # short-term to longer-term, the first 4 weeks, 8 weeks, 12 weeks, ...
  # mean absolute error (median v obs)
  # correlation
  # coverage: for 50%, 80%, 90%, 95%  CI
  # meas = d.obs$measure %>% unique
  date.f.start = d.proj$Week.start %>% as.Date %>% min
  # horizons = c(seq(4, 26, by = 4), 26)
  # horizons = 1:26
  d.proj.t = d.proj %>% # filter(measure == mea.t) %>% 
    dcast(., measure + fcast.type + Week.start ~ variable, value.var = 'value')
  d.obs.t = d.obs # %>% filter(measure == mea.t)
  da.t = left_join(x = d.proj.t  %>% filter(measure %in% unique(d.obs.t$measure)), y = d.obs.t, by = c('measure', 'Week.start')) %>% 
    mutate(wk.horizon = as.numeric((Week.start - date.f.start)/7 + 1)) %>% data.table()
  # should based on number of obs
  # remove incomplete records
  da.t = da.t[complete.cases(da.t)]
  
  
  tmp = da.t %>% dplyr::select(measure, fcast.type, wk.horizon)
  # 10/5/22
  tmp$abs.relerr = da.t[, c('value', 'median')] %>% apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[1]}) # apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[2]})
  tmp$abs.relerr.nodiv0 = da.t[, c('value', 'median')] %>% apply(., 1, FUN = function(x){abs(x[2]-x[1])/(x[1]+.1)}) # apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[2]})
  tmp$acc95ci = da.t[,c('value','ci95.lwr','ci95.upr')] %>% fn_acc.cover()
  tmp$acc90ci = da.t[,c('value','ci90.lwr','ci90.upr')] %>% fn_acc.cover()
  tmp$acc80ci = da.t[,c('value','ci80.lwr','ci80.upr')] %>% fn_acc.cover()
  tmp$acc50ci = da.t[,c('value','iqr.lwr','iqr.upr')] %>% fn_acc.cover()
  # also recode the accuracy for the n-wk ahead forecast
  tmp = tmp %>% 
    # recode
    mutate(acc.strict = case_when(abs.relerr <= cut0 | abs.relerr.nodiv0 <= cut0 ~ 1, T ~ 0),
           acc.loose1 = case_when(abs.relerr <= cut1 | abs.relerr.nodiv0 <= cut1 ~ 1, T ~ 0),
           acc.loose2 = case_when(abs.relerr <= cut2 | abs.relerr.nodiv0 <= cut2 ~ 1, T ~ 0)) 
  tmp
}
fn_eval.pt_arima = function(d.proj, d.obs){ # for point estimates: peak week, peek intensity
  # evaluate the accuracy of the projection, cp observation
  # criteria:
  # short-term to longer-term, the first 4 weeks, 8 weeks, 12 weeks, ...
  # mean absolute error (median v obs)
  # correlation
  # coverage: for 50%, 80%, 90%, 95%  CI
  # meas = d.obs$measure %>% unique
  date.f.start = d.proj$Week.start %>% as.Date %>% min
  # horizons = c(seq(4, 26, by = 4), 26)
  # horizons = 1:26
  d.proj.t = d.proj %>% # filter(measure == mea.t) %>% 
    dcast(., measure + fcast.type + Week.start ~ variable, value.var = 'value')
  d.obs.t = d.obs # %>% filter(measure == mea.t)
  da.t = left_join(x = d.proj.t  %>% filter(measure %in% unique(d.obs.t$measure)), y = d.obs.t, by = c('measure', 'Week.start')) %>% 
    # mutate(wk.horizon = as.numeric((Week.start - date.f.start)/7 + 1)) %>% 
    data.table()
  
  if(any(is.na(da.t$value))){
    pt = NULL # incomplete data
  } else {
    pt1 = da.t %>% filter(measure %in% key_measures) %>% group_by(fcast.type, measure) %>%  
      summarise(proj.pi = max(median), proj.tot = sum(median),
                obs.pi = max(value), obs.tot = sum(value)) %>%
      melt(., id.var = c('measure','fcast.type')) %>%
      mutate(target = factor(variable, levels = c('proj.pi','proj.tot','obs.pi','obs.tot'),
                             labels = rep(c('peak intensity', 'total'), 2)),
             type = factor(variable, levels = c('proj.pi','proj.tot','obs.pi','obs.tot'),
                           labels = rep(c('proj','obs'), e=2))) %>%
      dcast(., measure + fcast.type + target ~ type, value.var = 'value') %>% data.table()
    tmp1 = pt1 %>% dplyr::select(measure, fcast.type, target)
    tmp1$abs.relerr = pt1[, c('obs', 'proj')] %>% apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[1]}) # apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[2]})
    tmp1$abs.relerr.nodiv0 = pt1[, c('obs', 'proj')] %>% apply(., 1, FUN = function(x){abs(x[2]-x[1])/(x[1]+.1)}) # apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[2]})
    # also recode the accuracy for the n-wk ahead forecast
    tmp1 = tmp1 %>% 
      # recode
      mutate(acc.strict = case_when(abs.relerr <= cut0 | abs.relerr.nodiv0 <= cut0 ~ 1, T ~ 0),
             acc.loose1 = case_when(abs.relerr <= cut1 | abs.relerr.nodiv0 <= cut1 ~ 1, T ~ 0),
             acc.loose2 = case_when(abs.relerr <= cut2 | abs.relerr.nodiv0 <= cut2 ~ 1, T ~ 0)) 
    pt1 = tmp1 %>% melt(., id.var = c('measure','fcast.type','target'))
    
    # for peak week
    pt2 = da.t %>% filter(measure %in% key_measures) %>% group_by(fcast.type, measure) %>%  
      summarise(proj.pw = which.max(median), 
                obs.pw = which.max(value)) %>%
      mutate(acc.strict = case_when(proj.pw == obs.pw ~ 1,
                                    T ~ 0),
             acc.loose1 = case_when(abs(proj.pw - obs.pw) <= 1 ~ 1,
                                    T ~ 0),
             acc.loose2 = case_when(abs(proj.pw - obs.pw) <= 2 ~ 1,
                                    T ~ 0)) %>%
      mutate(target = 'peak week', proj.pw = NULL, obs.pw = NULL) %>% 
      melt(., id.var = c('measure','fcast.type','target'))
    
    pt = rbind(pt1, pt2) %>% setnames('variable','metric')
    
  }
  
  pt
}

median.local = function(x){
  if(all(is.na(x))){
    NA
  } else {
    median(x, na.rm = T)
  }
}


# get the point prediction accuracy evaluation
for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_out, 'res.eval_', loc.t,date.tag,'.RData')))
  if(class(tmp) != "try-error"){
    rm(res.eval)
    next
  }
  
  res.eval <-  lapply(files.t, function(x) {
    
    try(load(x))
    
    # print(x)
    
    id <- gsub(dir_res,'',x) %>% 
      strsplit('_') %>% unlist
    
    
    loc.t = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
    
    d = train.proj$fcast_stats
    if(is.null(d)){
      d = NULL
    } else {
      if(! is.null(train.proj$fcast_stats$fcast.type)){
        d.proj <- train.proj$fcast_stats %>%
          melt(id.vars = c('fcast.type', 'measure','Week.start'), variable.factor = F)
      } else {
        d.proj <- train.proj$fcast_stats %>%
          melt(id.vars = c('measure','Week.start'), variable.factor = F)
      }
      
      d.proj$measure = d.proj$measure %>% factor(., levels = c('case','death','cum.case','cum.death'),
                                                 labels = c('Cases','Deaths',"Cumulative Cases","Cumulative Deaths"))  
      
      
      variant.t = train.proj$variant
      
      # d.obs <- DAT.EPI[state == list.loc.names[list.locs==loc.t]]
      d.obs <- DAT.EPI %>% filter(state == loc.t, variant == variant.t)
      
      # setnames(d.obs, c('data.type','date'),c('measure','Week.start'))
      d.obs = d.obs %>% filter (as.Date(date) >= as.Date(start.date)) %>% dcast(., date ~ data.type, value.var = 'value')
      d.obs$cum.case = d.obs$case %>% cumsum()
      if(! is.null(d.obs$hosp))
        d.obs$cum.hosp = d.obs$hosp %>% cumsum
      d.obs$cum.death = d.obs$death %>% cumsum
      d.obs = d.obs %>% melt(., id.vars = 'date')
      d.obs$measure = d.obs$variable %>% factor(., levels = c('case', 'hosp','death','cum.case','cum.hosp','cum.death'),
                                                labels = c('Cases', 'Hospitalizations','Deaths',"Cumulative Cases","Cumulative Hospitalizations","Cumulative Deaths"))  
      d.obs$Week.start = d.obs$date %>% as.Date
      d.obs$variable = NULL; d.obs$date = NULL
      

      # only run it if there are complete data for evaluation 
      d1 = fn_eval.pt_arima(d.proj, d.obs) # for pt estimates
      
      d2 = fn_eval_arima(d.proj, 
                   d.obs)  %>% melt(., id.vars = c('measure','fcast.type','wk.horizon')) %>% 
        setnames('variable','metric') %>%
        mutate(target = factor(wk.horizon, levels = 1:26, labels = paste0(1:26,'wk ahead')),
               wk.horizon = NULL) %>% data.table()
      setcolorder(d2, c('measure','fcast.type','target','metric','value'))
      d = rbind(d1, d2[complete.cases(d2),]) # exclude NAs, , use.names = T
      
      d$fcast.start.month = train.proj$fcast.start.month
      d$fcast.start.week = train.proj$fcast.start.week #  d$Week.start %>% min %>% as.Date
      d$loc = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
      d$seasonality = id[grepl('sn',id)] # id[5]
      d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
      d$scenario = train.proj$scenario  # id[5]
      d$fcast.deflat = train.proj$fcast.deflat
      d$variant = train.proj$variant
    }
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, list(mean = round(mean(value, na.rm = T), 3), 
                          median  = round(median(value, na.rm = T), 3), 
                          ci50lwr = round(quantile(value, prob = .25, na.rm = T), 3), 
                          ci50upr = round(quantile(value, prob = .75, na.rm = T), 3),
                          ci95lwr = round(quantile(value, prob = .025, na.rm = T), 3), 
                          ci95upr = round(quantile(value, prob = .975, na.rm = T), 3)), 
                   by = list(variant, fcast.type, fcast.deflat, loc, seasonality, fcast.start.month, fcast.start.week, measure, metric, target, scenario)]) # %>%
  
  # res.eval$seasonality = factor(res.eval$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn'), labels = c('Seasonality assumed', 'Transformed seasonality', 'No seasonality'))
  res.eval$seasonality = factor(res.eval$seasonality, levels = sn_levels, 
                                 labels = sn_labels)
  
  save(res.eval, file = paste0(dir_out, 'res.eval_', loc.t,date.tag,'.RData'))
  rm(res.eval)
  print(paste(loc.t,'res.eval done'))
  
}

# prob score
for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_out, 'res.score_', loc.t,date.tag,'.RData')))
  if(class(tmp) != "try-error"){
    rm(res.score)
    next
  }
  
  res.score <-  lapply(files.t, function(x) {
    
    try(load(x))
    # print(x)
    id <- gsub(dir_res,'',x) %>% 
      strsplit('_') %>% unlist
    
    loc.t = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
    # nfcast = 26 # train.proj$fcast_stats$Week.start %>% unique %>% length
    # bins_pwk = 1: nfcast
    
    fcast = train.proj$fcastDist
    
    if(is.null(fcast)){
      d = NULL
    } else {
      fcast$wk.fcast = train.proj$fcast.start.week
      fcast$state = loc.t
      fcast$variant = train.proj$variant
      
      # in case there are multiple scenarios for the projections
      # to use the same scorer, rename 'fcast.type' to scenario
      fcast$scenario = fcast$fcast.type
      fcast$fcast.type = NULL
      
      d = getScore(fcast = fcast, 
                   truth = truth) %>% 
        setnames(c('state', 'wk.fcast', 'data.type'), c('loc', 'fcast.start.week', 'measure'))
      
      d = melt(d, id.vars = c('loc', 'variant', 'fcast.start.week', 'measure', 'target', 'scenario')) %>% 
        setnames(c('variable', 'scenario'), c('metric', 'fcast.type'))
      
      
      d$fcast.start.month = train.proj$fcast.start.month
      d$seasonality = id[grepl('sn',id)] # id[5]
      d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
      d$scenario = train.proj$scenario  # id[5]
      d$fcast.deflat = train.proj$fcast.deflat
      # d$variant = train.proj$variant
    }
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, list(mean = round(mean(value, na.rm = T), 3), 
                          median  = round(median(value, na.rm = T), 3), 
                          ci50lwr = round(quantile(value, prob = .25, na.rm = T), 3), 
                          ci50upr = round(quantile(value, prob = .75, na.rm = T), 3),
                          ci95lwr = round(quantile(value, prob = .025, na.rm = T), 3), 
                          ci95upr = round(quantile(value, prob = .975, na.rm = T), 3)), 
                   by = list(variant, fcast.type, fcast.deflat, loc, seasonality, fcast.start.month, fcast.start.week, measure, target, metric, scenario)]) # %>%
  
  res.score$seasonality = factor(res.score$seasonality, levels = sn_levels, 
                                 labels = sn_labels)
  
  
  save(res.score, file = paste0(dir_out, 'res.score_', loc.t,date.tag,'.RData'))
  rm(res.score)
  print(paste(loc.t, 'res.score done'))
  
}
