# combine results from multiple runs - diff eval method for the EAKF approaches

state.names = c('S1', 'E1','I1',
                'death1', 'newIobs1','newItot1','beta',
                'Tei','Tir','Trs','Td.mean','Td.sd',
                'p.mob','alpha','ifr')
  
files = list.files(dir_res, full.names = T)
files = files[grepl('.RData',files)]
files = files[grepl('train.proj', files)]
# files = files[grepl('Mean', files)] # transformed seasonal trend
# compare performance across runs?
length(files) # 15000 per state
length(files)

# too many!
# locs = DAT.EPI$state %>% unique

for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_res, 'res.train_', loc.t,date.tag,'.RData')))
  if(class(tmp) != "try-error"){
    rm(res.train)
    next
  }
  
  res.train <-  lapply(files.t, function(x) {
    
    try(load(x))
    
    id <- gsub(dir_res,'',x) %>% 
      strsplit('_') %>% unlist
    
    # get %S
    tmp = train.proj$states_stats[state == 'S1', ]
    tmp[, c("mean","median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr")] = tmp[, c("mean","median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr")] / N * 100
    tmp$state = 'Susceptibility'
    
    # cumulative infection rate
    tmp2 = train.proj$cumIperc_stats
    if(is.null(tmp2)){
      
      # cumulative immue loss
      tmp3 = train.proj$immLoss_stats
      
      d <- rbind(
        data.table(state = 'Rt',train.proj$Rt_stats),
        data.table(state = 'Rtx',train.proj$Rtx_stats),
        data.table(state = 'R0',train.proj$R0_stats),
        data.table(tmp),
        # data.table(tmp2),
        data.table(tmp3),
        data.table(train.proj$states_stats),
        fill = T
      ) %>%
        melt(id.vars = c('state','Week.start'), variable.factor = F)
    } else {
      tmp2$state = 'cumItot'; tmp2$sd = NULL;
      
      # cumulative immue loss
      tmp3 = train.proj$immLoss_stats
      
      d <- rbind(
        data.table(state = 'Rt',train.proj$Rt_stats),
        data.table(state = 'Rtx',train.proj$Rtx_stats),
        data.table(state = 'R0',train.proj$R0_stats),
        data.table(tmp),
        data.table(tmp2),
        data.table(tmp3),
        data.table(train.proj$states_stats),
        fill = T
      ) %>%
        melt(id.vars = c('state','Week.start'), variable.factor = F)
    }
    
    
    d$fcast.start.month = train.proj$fcast.start.month
    d$fcast.start.week = train.proj$fcast.start.week #  train.proj$fcast_stats$Week.start %>% min %>% as.Date
    d$loc = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
    
    d$seasonality = id[grepl('sn',id)] # id[5]
    d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
    d$fcast.type = train.proj$fcast.type  # id[5]
    d$fcast.deflat = train.proj$fcast.deflat
    d$variant = train.proj$variant
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, j = list(value = round(mean(value), 6)), by = list(variant, fcast.type, fcast.deflat, loc, seasonality, state, variable, fcast.start.month, fcast.start.week, Week.start)]) # %>%
  
  
  res.train$state = factor(res.train$state, levels=c("IimmLoss", "VimmLoss", 'Rt','R0','Rtx', 'Susceptibility', 'cumItot', state.names), 
                           labels = c('recoveree immune loss', 'vaccinee immune loss', 'Rt','R0','Rtx','Susceptibility', 'Cumulative infection rate', 'Susceptible', 'Exposed','Infectious',
                                      'death', 'case','infection','transmission rate',
                                      'latent period','infectious period','immunity period','Td.mean','Td.sd',
                                      'p.mob','infection detection rate','IFR') 
  )
  res.train$seasonality = factor(res.train$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn',
                                                                   paste0('w.sn.tran',1:8)), 
                                 labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality',paste0('w.sn.tran',1:8)))
  
  save(res.train, file = paste0(dir_res, 'res.train_', loc.t,date.tag,'.RData'))
  rm(res.train)
  print(paste(loc.t, 'res.train done'))
}

for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_res, 'res.proj_', loc.t,date.tag,'.RData')))
  if(class(tmp) != "try-error"){
    rm(res.proj)
    next
  }
  
  res.proj <-  lapply(files.t, function(x) {
    
    try(load(x))
    
    id <- gsub(dir_res,'',x) %>% 
      strsplit('_') %>% unlist
    
    if(! is.null(train.proj$fcast_stats$scenario)){
      d <- train.proj$fcast_stats %>%
        melt(id.vars = c('scenario', 'measure','Week.start'), variable.factor = F)
    } else {
      d <- train.proj$fcast_stats %>%
        melt(id.vars = c('measure','Week.start'), variable.factor = F)
    }
    
    
    d$fcast.start.month = train.proj$fcast.start.month
    d$fcast.start.week = train.proj$fcast.start.week #  d$Week.start %>% min %>% as.Date
    d$loc = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
    d$seasonality = id[grepl('sn',id)] # id[5]
    d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
    d$fcast.type = train.proj$fcast.type  # id[5]
    d$fcast.deflat = train.proj$fcast.deflat
    d$variant = train.proj$variant
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, j = list(value = round(mean(value), 0)), by = list(variant, fcast.type, fcast.deflat, loc, seasonality, scenario, measure, variable, fcast.start.month, fcast.start.week, Week.start)]) # %>%
  
  # res.proj$seasonality = factor(res.proj$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn'), labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality'))
  res.proj$seasonality = factor(res.proj$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn',
                                                                   paste0('w.sn.tran',1:8)), 
                                 labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality',paste0('w.sn.tran',1:8)))
  
  save(res.proj, file = paste0(dir_res, 'res.proj_', loc.t,date.tag,'.RData'))
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

fn_eval = function(d.proj, d.obs){
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
    dcast(., measure + scenario + Week.start ~ variable, value.var = 'value')
  d.obs.t = d.obs # %>% filter(measure == mea.t)
  da.t = left_join(x = d.proj.t  %>% filter(measure %in% unique(d.obs.t$measure)), y = d.obs.t, by = c('measure', 'Week.start')) %>% 
    mutate(wk.horizon = as.numeric((Week.start - date.f.start)/7 + 1)) %>% data.table()
  
  
  tmp = da.t %>% dplyr::select(measure, scenario, wk.horizon)
  # 10/5/22
  tmp$abs.relerr = da.t[, c('value', 'median')] %>% apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[1]}) # apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[2]})
  tmp$abs.relerr.nodiv0 = da.t[, c('value', 'median')] %>% apply(., 1, FUN = function(x){abs(x[2]-x[1])/(x[1]+.1)}) # apply(., 1, FUN = function(x){abs(x[2]-x[1])/x[2]})
  tmp$acc95ci = da.t[,c('value','ci95.lwr','ci95.upr')] %>% fn_acc.cover() # it is by row for each week
  tmp$acc90ci = da.t[,c('value','ci90.lwr','ci90.upr')] %>% fn_acc.cover()
  tmp$acc80ci = da.t[,c('value','ci80.lwr','ci80.upr')] %>% fn_acc.cover()
  tmp$acc50ci = da.t[,c('value','iqr.lwr','iqr.upr')] %>% fn_acc.cover()
  
  tmp
}

for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_res, 'res.eval_', loc.t,date.tag,'.RData')))
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
    
    if(! is.null(train.proj$fcast_stats$scenario)){
      d.proj <- train.proj$fcast_stats %>%
        melt(id.vars = c('scenario', 'measure','Week.start'), variable.factor = F)
    } else {
      d.proj <- train.proj$fcast_stats %>%
        melt(id.vars = c('measure','Week.start'), variable.factor = F)
    }
    
    variant.t = train.proj$variant
    
    # d.obs <- DAT.EPI[state == list.loc.names[list.locs==loc.t]]
    d.obs <- DAT.EPI %>% filter(state == loc.t, variant == variant.t)
    
    # setnames(d.obs, c('data.type','date'),c('measure','Week.start'))
    d.obs = d.obs %>% filter(as.Date(date) >= as.Date(start.date)) %>% dcast(., date ~ data.type, value.var = 'value')
    d.obs$cum.case = d.obs$case %>% cumsum()
    d.obs$cum.death = d.obs$death %>% cumsum
    d.obs = d.obs %>% melt(., id.vars = 'date')
    d.obs$measure = d.obs$variable %>% factor(., levels = c('case','death','cum.case','cum.death'),
                                              labels = c('Cases','Deaths',"Cumulative Cases","Cumulative Deaths"))  
    d.obs$Week.start = d.obs$date %>% as.Date
    d.obs$variable = NULL; d.obs$date = NULL
    
    # in case there are multiple scenarios for the projections
    sce_vec = d.proj$scenario %>% unique()
    if(is.null(sce_vec))
      d.proj$scenario = 'asIs' # just the baseline assuming business as usual
    
    d = fn_eval(d.proj, 
                  d.obs)  %>% melt(., id.vars = c('measure','scenario','wk.horizon')) %>% setnames('variable','metric')
    
    d$fcast.start.month = train.proj$fcast.start.month
    d$fcast.start.week = train.proj$fcast.start.week #  d$Week.start %>% min %>% as.Date
    d$loc = train.proj$loc # id[3] %>% strsplit('\\//') %>% unlist %>% tail(1)
    d$seasonality = id[grepl('sn',id)] # id[5]
    d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
    d$fcast.type = train.proj$fcast.type  # id[5]
    d$fcast.deflat = train.proj$fcast.deflat
    d$variant = train.proj$variant
    
    d
  }) %>%
    rbindlist() %>%
    (function(d) d[, list(mean = round(mean(value, na.rm = T), 3), 
                          median  = round(median(value, na.rm = T), 3), 
                          ci50lwr = round(quantile(value, prob = .25, na.rm = T), 3), 
                          ci50upr = round(quantile(value, prob = .75, na.rm = T), 3),
                          ci95lwr = round(quantile(value, prob = .025, na.rm = T), 3), 
                          ci95upr = round(quantile(value, prob = .975, na.rm = T), 3)), 
                   by = list(variant, fcast.type, fcast.deflat, loc, seasonality, fcast.start.month, fcast.start.week, measure, metric, wk.horizon, scenario)]) # %>%
  
  # res.eval$seasonality = factor(res.eval$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn'), labels = c('Seasonality assumed', 'Transformed seasonality', 'No seasonality'))
  res.eval$seasonality = factor(res.eval$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn',
                                                                   paste0('w.sn.tran',1:8)), 
                                 labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality',paste0('w.sn.tran',1:8)))
  
  save(res.eval, file = paste0(dir_res, 'res.eval_', loc.t,date.tag,'.RData'))
  rm(res.eval)
  print(paste(loc.t,'res.eval done'))
  
}

# prob score
for(loc.t in locs){
  files.t = files[grepl(gsub(' ', '', loc.t), files)]
  if(length(files.t) < 1)
    next
  
  print(loc.t)
  
  tmp = try(load(paste0(dir_res, 'res.score_', loc.t,date.tag,'.RData')))
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
    
    fcast$wk.fcast = train.proj$fcast.start.week
    fcast$state = loc.t
    fcast$variant = train.proj$variant
    
    # in case there are multiple scenarios for the projections
    sce_vec = fcast$scenario %>% unique()
    if(is.null(sce_vec))
      fcast$scenario = 'asIs' # just the baseline assuming business as usual
    
    d = getScore(fcast = fcast, 
                   truth = truth) %>% 
      setnames(c('state', 'wk.fcast', 'data.type'), c('loc', 'fcast.start.week', 'measure'))
    
    d = melt(d, id.vars = c('loc', 'variant', 'fcast.start.week', 'measure', 'target', 'scenario')) %>% setnames('variable', 'metric')
    
    d$fcast.start.month = train.proj$fcast.start.month
    d$seasonality = id[grepl('sn',id)] # id[5]
    d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
    d$fcast.type = train.proj$fcast.type  # id[5]
    d$fcast.deflat = train.proj$fcast.deflat
    # d$variant = train.proj$variant
    
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
  
  # res.score$seasonality = factor(res.score$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn'), labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality'))
  res.score$seasonality = factor(res.score$seasonality, levels = c('w.sn', 'w.sn.tran', 'no.sn',
                                                                   paste0('w.sn.tran',1:8)), 
                                 labels = c('Seasonality assumed','Transformed seasonality', 'No seasonality',paste0('w.sn.tran',1:8)))
  
  
  save(res.score, file = paste0(dir_res, 'res.score_', loc.t,date.tag,'.RData'))
  rm(res.score)
  print(paste(loc.t, 'res.score done'))
  
}

