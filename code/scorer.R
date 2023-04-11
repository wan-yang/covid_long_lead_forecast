# compute the log score of the forecast

# load(paste0(dir_data, 'truth.RData'))

# fcast = train.proj$fcastDist
# fcast$wk.fcast = fcast.start
# fcast$state = loc.t

logScore <- function(x){
  if(length(x))
    log(sum(x, na.rm = T)) %>% round(4)
}


fn_get.bin_incl_loose = function(bin_incl_strict,bin_incl_loose){
  idx = which(bin_incl_strict == 1)
  if(length(idx) > 0 )
    if (idx != length(bin_incl_loose) & idx!=1)
      bin_incl_loose[idx + c(-1:1)] = c(1,1,1)
  else if(idx == 1 & idx != length(bin_incl_loose) )
    bin_incl_loose[idx + c(0:1)] = c(1,1)
  else if(idx != 1 & idx == length(bin_incl_loose) )
    bin_incl_loose[idx + c(-1:0)] = c(1,1)
  bin_incl_loose
}

fn_get.logScore.strict = function(bin_incl_strict,prob,logS.lowest = -10){
  idx = which(bin_incl_strict==1) 
  res = logScore(prob[idx])
  if(is.null(res)){
    res = logS.lowest
  } else if (is.infinite(res)){
    res = logS.lowest
  }
  res
}


fn_get.logScore.loose = function(bin_incl_loose,prob,logS.lowest = -10){
  idx = which(bin_incl_loose==1) 
  res = logScore(prob[idx])
  if(is.null(res)){
    res = logS.lowest
  } else if (is.infinite(res)){
    res = logS.lowest
  }
  
  res
}

getScore = function(fcast, truth){
  truth$wk.fcast = truth$wk.fcast %>% as.Date
  loc_vec = fcast$state %>% unique
  wk.fcast_vec = fcast$wk.fcast %>% unique %>% as.character()
  dtype_vec = fcast$data.type %>% unique
  target_vec = fcast$target %>% unique
  
  da.t = fcast %>% left_join(x = ., y = truth, by = c('data.type', 'state', 'wk.fcast', 'target', 'variant'))
  da.t$wk.fcast = da.t$wk.fcast %>% as.Date
  
  # remove NA: weeks when there are no observations for evaluation (otherwise, would be assigned lowest score)
  da.t = da.t %>% filter(!is.na(value))
  
  da.t = da.t %>% mutate(value = case_when(target == 'peak week' ~ as.numeric((as.Date(value, origin = '1970/1/1')  - wk.fcast)/7+1),
                                           T ~ value))
  
  da.t = da.t %>% mutate(bin_incl_strict = case_when(value >= lwr & value <= upr ~ 1,
                                                     T ~ 0)) %>%
    mutate(bin_incl_loose = bin_incl_strict) %>% data.table()
  
  da.t = da.t[, list(lwr, upr,prob, value, bin_incl_strict,
                     bin_incl_loose = fn_get.bin_incl_loose(bin_incl_strict,bin_incl_loose)),
              by = c('target', 'state', 'data.type', 'variant', 'wk.fcast', 'scenario')] %>% data.table()
  
  score.t = da.t[, list(logScore.strict = fn_get.logScore.strict(bin_incl_strict,prob),
                        logScore.loose = fn_get.logScore.loose(bin_incl_loose,prob)),
                 by = c('target', 'state', 'data.type', 'variant', 'wk.fcast', 'scenario')]
  
  score.t
}



