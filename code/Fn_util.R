# some utility functions for analysis
# 2/5/20

# segment linear interpolation - fill in missing data based on the neighboring data
fn_segIntrpl = function(da.t, var.t){
  
  da.t$idx = 1:nrow(da.t)
  
  
  # do the two dose separately, plus boosters
  for(dv in var.t){
    
    da.t$new.v = 0
    
    idx.all = which(is.na(da.t[,dv,with=F]))
    
    
    consecutive =  idx.all[-1] - idx.all[-length(idx.all)]
    i.div = which(consecutive >=2)
    w.div = idx.all[which(consecutive >=2)]
    grps = list()
    if(length(i.div)==0 & length(idx.all) <=2){ # only 2 NAs and they are adjacent
      grps[[1]] = idx.all
    } else if(length(i.div) == 1){
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
                (tail(idx0,1)+ n.add) %>% pmin(nrow(da.t)),
                by = 1
      ) # add an additional week in case there is back adjustment
      
      da.t0 = da.t[idx]; 
      da.t0 = da.t0[complete.cases(da.t0[,c(dv), with =F])]
      
      eval(parse(text = paste('fit = lm(`',dv,'` ~ idx, data = da.t0)', sep='')))
      da.t$new.v[idx0] = predict(fit, newdata = da.t[idx0])
      
    }
    
    eval(parse(text = paste('da.t[is.na(`',dv,'`)]$`',dv, '`= da.t[is.na(`',dv,'`)]$new.v', sep='')))

  } # do each variable separately
  
  da.t[da.t < 0] = 0
  da.t[, var.t, with = F]
}



fn_find_1st.day.aboveCut=function(x,cut){
  tmp=which(x>=cut)
  if(length(tmp)>0){
    return(tmp[1])
  } else {
    return(NA)
  }
}

fn_get_stats = function(tda){
  # data format: age group, num_ens, days
  mean = apply(tda, c(1,3), mean)
  median = apply(tda, c(1,3), median)
  sd = apply(tda, c(1,3), sd)
  ci.lwr = apply(tda, c(1,3), quantile, .025)
  ci.upr = apply(tda, c(1,3), quantile, .975)
  
  return(list(mean=mean, median=median, ci.lwr=ci.lwr, ci.upr=ci.upr, sd=sd))
}

# combine forecasts
fn_ens.stats=function(means,sds){
  vars=sds^2;
  ens.mean=mean(means);
  ens.var = mean(vars) + (1+1/length(means))*1/(length(means)-1)*sum((means-ens.mean)^2)
  ci95.low=ens.mean-1.96*sqrt(ens.var)
  ci95.up=ens.mean+1.96*sqrt(ens.var)
  ci50.low=ens.mean-.67*sqrt(ens.var)
  ci50.up=ens.mean+.67*sqrt(ens.var)
  ci90.low=ens.mean-1.28*sqrt(ens.var)
  ci90.up=ens.mean+1.28*sqrt(ens.var)
  stats = c(ens.mean,ens.var,ci95.low,ci95.up,ci90.low,ci90.up,ci50.low,ci50.up)
  names(stats) = c('mean','var','CI95lwr','CI95upr','CI90lwr','CI90upr','CI50lwr','CI50upr')
  return(stats)
}


fn_format= function(x,roundit=T,roundigt=0, scientific = T){
  x=unlist(x) %>% as.numeric()
  if(roundit==T){
    if(scientific){
      # only do it for small numbers
      x[x <= 1e-3 & !is.na(x)] =signif(x[x <= 1e-3 & !is.na(x)],roundigt)
      x[x > 1e-3 & !is.na(x)] = round(x[x > 1e-3 & !is.na(x)],roundigt)
    } else {
      x=round(x,roundigt)
    }
  }
  paste0(x[1],' (',x[2],', ',x[3],')')
}
fn_formatCI= function(x,roundit=T,roundigt=0, scientific = T){
  x=unlist(x) %>% as.numeric()
  if(roundit==T){
    if(scientific){
      # only do it for small numbers
      x[x <= 1e-3 & !is.na(x)] =signif(x[x <= 1e-3 & !is.na(x)],roundigt)
      x[x > 1e-3 & !is.na(x)] = round(x[x > 1e-3 & !is.na(x)],roundigt)
    } else {
      x=round(x,roundigt)
    }
  }
  paste0('(',x[1],', ',x[2],')')
}


# key functions for processing the forecast
fn_getVariant = function(week.t, loc.t = loc.t, type.t = 'predominant', vdates.dom = vdates.dom){
  # to find the predomiant variant for week.t
  vdates.t = vdates.dom %>% filter(location == loc.t & type == type.t & 
                                     as.Date(week.t) >= date.start & as.Date(week.t) <= date.end &
                                     ! variant %in% c('Omicron'))
  if(nrow(vdates.t) == 1){
    main.variant = ifelse(vdates.t$variant!='Others', vdates.t$variant %>% as.character(), NA_character_)
  } else {
    main.variant = NA_character_
  }
  main.variant
}
fn_getWave = function(week.t, loc.t = loc.t, type.t = 'predominant', vdates.dom = vdates.dom, variant = variant){
  # to find the predomiant variant for week.t
  # also include variant, b/c there are some overlaps b/w delta and omicron
  # if so base on the variant label
  
  if(is.na(week.t)){
    wave.t = NA_character_
  } else {
    vdates.t = vdates.dom %>% filter(location == loc.t & type == type.t & 
                                       as.Date(week.t) >= date.start & as.Date(week.t) <= date.end &
                                       ! variant %in% c('Omicron'))
    b4.alpha.t = vdates.dom %>% filter(location == loc.t & type == type.t & 
                                         variant == 'Alpha')
    b4.delta.t = vdates.dom %>% filter(location == loc.t & type == type.t & 
                                         variant == 'Delta')
    
    if(nrow(b4.alpha.t) <1){
      # no alpha wave, try delta?
      if(nrow(b4.delta.t) == 1){
        wave2end = min(as.Date('2021/5/31'), as.Date(b4.delta.t$date.start) - 7)
      } else {
        wave2end = as.Date('2021/5/31') # arbitrarily set to end of May 2021
      }
    } else {
      wave2end = as.Date(b4.alpha.t$date.start) - 7
    }
    
    if(as.Date(week.t) <= wave2end){
      wave.t = 'wave2'
    } else {
      if(nrow(vdates.t) == 1){
        wave.t = ifelse(vdates.t$variant!='Others', vdates.t$variant %>% as.character(), NA_character_)
      } else {
        wave.t = NA_character_
      }
    }
  }
  
  if(!is.na(wave.t)){
    # check labeling for overlapping weeks of delta and omicron
    if(wave.t == 'Delta' & variant == 'Omicron')
      wave.t = 'Omicron_BA.1'  # relabel it as BA.1 because only Omicron/Delta specific data are used here
    # if the forecast is made using non-Omicron data (e.g., long horizon, 17-26wk ahead)
    # they could be counted as BA.1 or others, need to refine it as Delta, as it's delta data were used
    if(variant == 'non.Omicron' & grepl('omicron', tolower(wave.t)))
      wave.t = 'Delta'  # relabel it as BA.1 because only Omicron/Delta specific data are used here
  }
  
  wave.t
}


fn_getStats.byloc = function(da.t, loc.t, variant.t, measure.t, metric.t, target2.t, v.cp, v.ref, variable.t = 'mean'){
  setDT(da.t)
  tda = da.t %>% dplyr::filter(loc == loc.t & variant == variant.t & measure == measure.t & metric == metric.t & target2 == target2.t & variable == variable.t) %>% data.table()
  # first test if it is normally distributed
  mean.cp = mean({eval(parse(text = paste(tda[,v.cp,with=F])))}) # %>% signif(digits = 3)
  mean.ref = mean({eval(parse(text = paste(tda[,v.ref,with=F])))}) # %>% signif(digits = 3)
  diff.mean = tda$diff %>% mean() # %>% signif(digits = 3)
  # cumulative relative diff, % relative diff
  # perc.diff = ((sum({eval(parse(text = paste(tda[,v.cp,with=F])))}) - sum({eval(parse(text = paste(tda[,v.ref,with=F])))})) / 
  #   abs(sum({eval(parse(text = paste(tda[,v.ref,with=F])))})) * 100) %>% round(1)
  # if it is log score - exponentiate it first
  if(mean.cp <= 0 & mean.ref <= 0){
    perc.diff = ((exp(mean.cp) - exp(mean.ref)) / 
                   exp(mean.ref) * 100) %>% signif(digits = 3)
  } else {
    perc.diff = ((mean.cp - mean.ref) / 
                   (mean.ref) * 100) %>% signif(digits = 3) # round(1)
  }
  
  n.smp = tda %>% nrow()
  # t.norm = shapiro.test(tda$diff)
  # pvalue.shapiro = t.norm$p.value # %>% format(digits = 2)
  t.wrs = wilcox.test(tda$diff, alternative = 'two.sided', conf.int = T)
  diff.median.wrs = t.wrs$estimate %>% unname() # %>% signif(digits = 3)
  diff.median.wrs.ci95lwr = t.wrs$conf.int[1] # %>% signif(digits = 3)
  diff.median.wrs.ci95upr = t.wrs$conf.int[2] # %>% signif(digits = 3)
  pvalue.wrs.2side = t.wrs$p.value  
  pvalue.wrs.1side = NA_real_
  if(!is.na(pvalue.wrs.2side)){
    if(pvalue.wrs.2side < .05){
      # yes, the median is not 0, test it again with a 1-sided test
      pvalue.wrs.1side = wilcox.test(tda$diff, alternative = ifelse(diff.median.wrs < 0, 'less', 'greater'), conf.int = F)$p.value  # %>% format(digits = 2)
    }
  }
  
  pvalue.wrs = ifelse(is.na(pvalue.wrs.1side), pvalue.wrs.2side, pvalue.wrs.1side)
  # return(list(diff.mean = diff.mean, perc.diff = perc.diff, n.smp = n.smp, 
  #             diff.median.wrs = diff.median.wrs, diff.median.wrs.ci95lwr = diff.median.wrs.ci95lwr, diff.median.wrs.ci95upr = diff.median.wrs.ci95upr,
  #             pvalue.wrs.2side = pvalue.wrs.2side, pvalue.wrs.1side = pvalue.wrs.1side, pvalue.wrs = pvalue.wrs, pvalue.shapiro = t.norm$p.value))
  
  return(c(n.smp, mean.ref %>% signif(digits = 3), mean.cp %>% signif(digits = 3), diff.mean %>% signif(digits = 3), perc.diff,  
           diff.median.wrs %>% signif(digits = 3), diff.median.wrs.ci95lwr %>% signif(digits = 3), diff.median.wrs.ci95upr %>% signif(digits = 3),
           pvalue.wrs.2side %>% signif(digits = 2), 
           pvalue.wrs.1side %>% signif(digits = 2), 
           pvalue.wrs %>% signif(digits = 2))) # , pvalue.shapiro %>% signif(digits = 2)
  
}

# also by respiratory virus season
fn_getStats.byloc.bysn = function(da.t, loc.t, variant.t, resp.sn.t, measure.t, metric.t, target2.t, v.cp, v.ref, variable.t = 'mean'){
  setDT(da.t)
  
  tda = da.t %>% dplyr::filter(loc == loc.t & variant == variant.t & resp.sn == resp.sn.t & measure == measure.t & metric == metric.t & target2 == target2.t & variable == variable.t) %>% data.table()
  n.smp = tda %>% nrow()
  
  # print(paste(loc.t, variant.t, resp.sn.t, measure.t, metric.t, target2.t, n.smp))
  
  if(n.smp < 2){ # need at least 2
    # return(c(n.smp, rep(NA_real_, 7), rep(NA_character_,3))) # , pvalue.shapiro %>% format(digits = 2)
    return(c(n.smp, rep(NA_real_, 10))) 
  } else {
    # first test if it is normally distributed
    mean.cp = mean({eval(parse(text = paste(tda[,v.cp,with=F])))}) # %>% signif(digits = 3)
    mean.ref = mean({eval(parse(text = paste(tda[,v.ref,with=F])))}) # %>% signif(digits = 3)
    diff.mean = tda$diff %>% mean() # %>% signif(digits = 3)
    # cumulative relative diff, % relative diff
    # perc.diff = ((sum({eval(parse(text = paste(tda[,v.cp,with=F])))}) - sum({eval(parse(text = paste(tda[,v.ref,with=F])))})) / 
    #   abs(sum({eval(parse(text = paste(tda[,v.ref,with=F])))})) * 100) %>% round(1)
    # if it is log score - exponentiate it first
    if(mean.cp <= 0 & mean.ref <= 0){
      perc.diff = ((exp(mean.cp) - exp(mean.ref)) / 
                     exp(mean.ref) * 100) %>% signif(digits = 3)
    } else {
      perc.diff = ((mean.cp - mean.ref) / 
                     (mean.ref) * 100) %>% signif(digits = 3) # round(1)
    }
    
    # t.norm = shapiro.test(tda$diff)
    # pvalue.shapiro = t.norm$p.value # %>% format(digits = 2)
    t.wrs = wilcox.test(tda$diff, alternative = 'two.sided', conf.int = T)
    diff.median.wrs = t.wrs$estimate %>% unname() # %>% signif(digits = 3)
    diff.median.wrs.ci95lwr = t.wrs$conf.int[1] # %>% signif(digits = 3)
    diff.median.wrs.ci95upr = t.wrs$conf.int[2] # %>% signif(digits = 3)
    pvalue.wrs.2side = t.wrs$p.value  
    pvalue.wrs.1side = NA_real_
    if(!is.na(pvalue.wrs.2side)){
      if(pvalue.wrs.2side < .05){
        # yes, the median is not 0, test it again with a 1-sided test
        pvalue.wrs.1side = wilcox.test(tda$diff, alternative = ifelse(diff.median.wrs < 0, 'less', 'greater'), conf.int = F)$p.value  # %>% format(digits = 2)
      }
    }
    
    pvalue.wrs = ifelse(is.na(pvalue.wrs.1side), pvalue.wrs.2side, pvalue.wrs.1side)
    # return(list(diff.mean = diff.mean, perc.diff = perc.diff, n.smp = n.smp, 
    #             diff.median.wrs = diff.median.wrs, diff.median.wrs.ci95lwr = diff.median.wrs.ci95lwr, diff.median.wrs.ci95upr = diff.median.wrs.ci95upr,
    #             pvalue.wrs.2side = pvalue.wrs.2side, pvalue.wrs.1side = pvalue.wrs.1side, pvalue.wrs = pvalue.wrs, pvalue.shapiro = t.norm$p.value))
    
    return(c(n.smp, mean.ref %>% signif(digits = 3), mean.cp  %>% signif(digits = 3), diff.mean  %>% signif(digits = 3), perc.diff,  
             diff.median.wrs %>% signif(digits = 3), diff.median.wrs.ci95lwr  %>% signif(digits = 3), diff.median.wrs.ci95upr %>% signif(digits = 3),
             pvalue.wrs.2side %>% signif(digits = 2), 
             pvalue.wrs.1side %>% signif(digits = 2), 
             pvalue.wrs %>% signif(digits = 2))) # , pvalue.shapiro %>% signif(digits = 2)
  }
}
fn_getStats.alllocs = function(da.t, variant.t, measure.t, metric.t, target2.t, v.cp, v.ref, variable.t = 'mean'){
  setDT(da.t)
  tda = da.t %>% dplyr::filter(variant == variant.t & measure == measure.t & metric == metric.t & target2 == target2.t & variable == variable.t) %>% data.table()
  # first test if it is normally distributed
  mean.cp = mean({eval(parse(text = paste(tda[,v.cp,with=F])))}) %>% signif(digits = 3)
  mean.ref = mean({eval(parse(text = paste(tda[,v.ref,with=F])))}) %>% signif(digits = 3)
  diff.mean = tda$diff %>% mean() %>% signif(digits = 3)
  # cumulative relative diff, % relative diff
  # perc.diff = ((sum({eval(parse(text = paste(tda[,v.cp,with=F])))}) - sum({eval(parse(text = paste(tda[,v.ref,with=F])))})) / 
  #   abs(sum({eval(parse(text = paste(tda[,v.ref,with=F])))})) * 100) %>% round(1)
  # if it is log score - exponentiate it first
  if(mean.cp <= 0 & mean.ref <= 0){
    perc.diff = ((exp(mean.cp) - exp(mean.ref)) / 
                   exp(mean.ref) * 100) %>% signif(digits = 3)
  } else {
    perc.diff = ((mean.cp - mean.ref) / 
                   (mean.ref) * 100) %>% signif(digits = 3) # round(1)
  }
  n.smp = tda %>% nrow()
  # t.norm = shapiro.test(tda$diff)
  # pvalue.shapiro = try(t.norm$p.value) # %>% format(digits = 2)
  t.wrs = wilcox.test(tda$diff, alternative = 'two.sided', conf.int = T)
  diff.median.wrs = t.wrs$estimate %>% unname()  %>% signif(digits = 3)
  diff.median.wrs.ci95lwr = t.wrs$conf.int[1] %>% signif(digits = 3)
  diff.median.wrs.ci95upr = t.wrs$conf.int[2] %>% signif(digits = 3)
  pvalue.wrs.2side = t.wrs$p.value  
  pvalue.wrs.1side = NA_real_
  if(pvalue.wrs.2side < .05){
    # yes, the median is not 0, test it again with a 1-sided test
    pvalue.wrs.1side = wilcox.test(tda$diff, alternative = ifelse(diff.median.wrs < 0, 'less', 'greater'), conf.int = F)$p.value  # %>% format(digits = 2)
  }
  pvalue.wrs = ifelse(is.na(pvalue.wrs.1side), pvalue.wrs.2side, pvalue.wrs.1side)
  # return(list(diff.mean = diff.mean, perc.diff = perc.diff, n.smp = n.smp, 
  #             diff.median.wrs = diff.median.wrs, diff.median.wrs.ci95lwr = diff.median.wrs.ci95lwr, diff.median.wrs.ci95upr = diff.median.wrs.ci95upr,
  #             pvalue.wrs.2side = pvalue.wrs.2side, pvalue.wrs.1side = pvalue.wrs.1side, pvalue.wrs = pvalue.wrs, pvalue.shapiro = t.norm$p.value))
  
  return(c(n.smp, mean.ref, mean.cp, diff.mean, perc.diff,  
           diff.median.wrs, diff.median.wrs.ci95lwr, diff.median.wrs.ci95upr,
           pvalue.wrs.2side %>% format(digits = 2), 
           pvalue.wrs.1side %>% format(digits = 2), 
           pvalue.wrs %>% format(digits = 2)))
  
}

fn_getAllStats = function(tda.t, v.cp.t, v.ref.t){
  setDT(tda.t)
  # 1 aggregate by loc, by target; also aggregate by wave, but should do it seperately so won't double count
  stats.t.byloc = tda.t[, as.list(fn_getStats.byloc(da.t = tda.t, loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                    v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                        by = list(loc,variant, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>%
            mutate(variant = wave) %>% data.table() %>% .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                     v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                         by = list(loc,variant, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc))
  colnames(stats.t.byloc)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                   'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 2 aggregate, by target, all locs; also aggregate by wave, but should do it separately so won't double count
  stats.t.alllocs = tda.t %>% mutate(loc = 'all') %>% data.table() %>% # aggregate over all location
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>% # mutate(loc = 'all') %>% 
            mutate(variant = wave, loc = 'all') %>% data.table() %>% .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                                      by = list(loc,variant, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs))
  colnames(stats.t.alllocs)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                     'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 3 aggregate by loc, all targets; also aggregate by wave, but should do it separately so won't double count
  stats.t.byloc.alltargets = tda.t %>% mutate(target2 = 'all') %>% data.table() %>% # aggregate over all location
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>% # mutate(target2 = 'all') %>% 
            mutate(variant = wave, target2 = 'all') %>% data.table() %>% .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                                      v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                                          by = list(loc,variant, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.alltargets))
  colnames(stats.t.byloc.alltargets)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                              'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 4 aggregate by loc, by target; all times
  stats.t.byloc.bylarget.alltimes = tda.t %>% mutate(variant = 'all') %>% data.table() %>% # aggregate over time
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.bylarget.alltimes))
  colnames(stats.t.byloc.bylarget.alltimes)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                     'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 5 aggregate all locs, by target; all times
  stats.t.alllocs.bylarget.alltimes = tda.t %>% mutate(variant = 'all', loc = 'all') %>% data.table() %>% # aggregate over time
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.bylarget.alltimes))
  colnames(stats.t.alllocs.bylarget.alltimes)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                       'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 6 aggregate by loc, all targets; all times
  stats.t.byloc.allltargets.alltimes = tda.t %>% mutate(variant = 'all', target2 = 'all') %>% data.table() %>% # aggregate over time
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.allltargets.alltimes))
  colnames(stats.t.byloc.allltargets.alltimes)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                        'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  
  # 7 aggregate all locs, all targets; also aggregate by wave, but should do it separately so won't double count
  stats.t.alllocs.alltargets = tda.t %>% mutate(loc = 'all', target2 = 'all') %>% data.table() %>% # aggregate over all location
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>% mutate(loc = 'all', target2 = 'all') %>% 
            mutate(variant = wave, target2 = 'all') %>% data.table() %>% .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                                      v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                                          by = list(loc,variant, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.alltargets))
  colnames(stats.t.alllocs.alltargets)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 8 aggregate all locs, all targets; all times
  stats.t.alllocs.alltargets.alltimes = tda.t %>% mutate(loc = 'all', target2 = 'all', variant = 'all')  %>% data.table() %>% # aggregate over all locations, targets, times/variants
    .[, as.list(fn_getStats.byloc(da.t = ., loc.t = loc, variant.t = variant, measure.t = measure, metric.t = metric, target2.t = target2, 
                                  v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.alltargets.alltimes))
  colnames(stats.t.alllocs.alltargets.alltimes)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                         'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # BY RESPIRATORY VIRUS SEASON
  # 9 aggregate by locs, by targets, by respiratory season v off-season, by variant
  stats.t.byloc.bysn = tda.t[, as.list(fn_getStats.byloc.bysn(da.t = tda.t, loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                              v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                             by = list(loc,variant, resp.sn, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>%
            mutate(variant = wave) %>% data.table() %>% .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                          v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                         by = list(loc,variant, resp.sn, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.bysn))
  colnames(stats.t.byloc.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                        'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 10 aggregate, by target, all locs; by respiratory season v off-season, by variant
  stats.t.alllocs.bysn = tda.t %>% mutate(loc = 'all') %>% data.table() %>% # aggregate over all location
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, resp.sn, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>% mutate(loc = 'all') %>%
            mutate(variant = wave) %>% data.table() %>% .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                          v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                         by = list(loc,variant, resp.sn, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.bysn))
  colnames(stats.t.alllocs.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                          'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 11 aggregate by loc, all targets; by respiratory season v off-season, by variant
  stats.t.byloc.alltargets.bysn = tda.t %>% mutate(target2 = 'all') %>% data.table() %>% # aggregate over all location
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant,resp.sn, measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>% mutate(target2 = 'all') %>% 
            mutate(variant = wave) %>% data.table() %>% .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                          v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                         by = list(loc,variant, resp.sn, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.alltargets.bysn))
  colnames(stats.t.byloc.alltargets.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                   'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  
  # 12 aggregate by loc, by target; by.sn, all times
  stats.t.byloc.bylarget.alltimes.bysn = tda.t %>% mutate(variant = 'all') %>% data.table() %>% # aggregate over time
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn,measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, resp.sn, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.bylarget.alltimes.bysn))
  colnames(stats.t.byloc.bylarget.alltimes.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                          'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 13 aggregate all locs, by target; by.sn, all times
  stats.t.alllocs.bylarget.alltimes.bysn = tda.t %>% mutate(variant = 'all', loc = 'all') %>% data.table() %>% # aggregate over time
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, resp.sn, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.bylarget.alltimes.bysn))
  colnames(stats.t.alllocs.bylarget.alltimes.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                            'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 14 aggregate by loc, alls target; by.sn, all times
  stats.t.byloc.allltargets.alltimes.bysn = tda.t %>% mutate(variant = 'all', target2 = 'all') %>% data.table() %>% # aggregate over time
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant,resp.sn, measure, metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.byloc.allltargets.alltimes.bysn))
  colnames(stats.t.byloc.allltargets.alltimes.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                             'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 15 aggregate all locs, all targets; by respiratory season v off-season, by variant
  stats.t.alllocs.alltargets.bysn = tda.t %>% mutate(loc = 'all', target2 = 'all') %>% data.table() %>% # aggregate over all location
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant,resp.sn,  measure,metric,target2)] %>% 
    mutate(wave = NA) %>% 
    rbind(., 
          # aggregate by wave for the pre-Omicron period
          tda.t %>% dplyr::filter(wave %in% c('wave2','Alpha','Delta')) %>% mutate(loc = 'all', target2 = 'all') %>% 
            mutate(variant = wave) %>% data.table() %>% .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                                                          v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
                                         by = list(loc,variant, resp.sn, measure,metric,target2)] %>%
            mutate(wave = variant, variant = 'non.Omicron')
    ) %>% rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.alltargets.bysn))
  colnames(stats.t.alllocs.alltargets.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                     'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # 16 aggregate all locs, all targets; all times
  stats.t.alllocs.alltargets.alltimes.bysn = tda.t %>% mutate(loc = 'all', target2 = 'all', variant = 'all') %>% data.table() %>% # aggregate over all locations, targets, times/variants
    .[, as.list(fn_getStats.byloc.bysn(da.t = ., loc.t = loc, variant.t = variant, resp.sn.t = resp.sn, measure.t = measure, metric.t = metric, target2.t = target2, 
                                       v.cp = v.cp.t, v.ref = v.ref.t, variable.t = 'mean')),
      by = list(loc,variant, resp.sn, measure,metric,target2)] %>% 
    rename(target = target2)
  idx = grep('V',colnames(stats.t.alllocs.alltargets.alltimes.bysn))
  colnames(stats.t.alllocs.alltargets.alltimes.bysn)[idx] = c('n.smp', 'mean.base', 'mean.opt', 'diff.mean', 'perc.diff',  'diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr', 
                                                              'pvalue.wrs.2side', 'pvalue.wrs.1side', 'pvalue.wrs')
  
  # combine all
  stats.t.all = rbind(stats.t.byloc, stats.t.alllocs, stats.t.byloc.alltargets, stats.t.byloc.bylarget.alltimes,
                      stats.t.alllocs.bylarget.alltimes, stats.t.byloc.allltargets.alltimes, stats.t.alllocs.alltargets, stats.t.alllocs.alltargets.alltimes,
                      # by resp.sn
                      stats.t.byloc.bysn, stats.t.alllocs.bysn, stats.t.byloc.alltargets.bysn, stats.t.byloc.bylarget.alltimes.bysn,
                      stats.t.alllocs.bylarget.alltimes.bysn, stats.t.byloc.allltargets.alltimes.bysn, stats.t.alllocs.alltargets.bysn, stats.t.alllocs.alltargets.alltimes.bysn, 
                      fill = T)
  # format 
  stats.t.all$median.diff = stats.t.all[,c('diff.median.wrs', 'diff.median.wrs.ci95lwr', 'diff.median.wrs.ci95upr'), with=F] %>% apply(., 1, fn_format, roundigt = 2)
  setcolorder(stats.t.all, c('loc',"variant", 'wave', 'resp.sn', "measure","metric",
                             "target","n.smp","mean.base","mean.opt","diff.mean","perc.diff","median.diff","pvalue.wrs", 
                             "diff.median.wrs","diff.median.wrs.ci95lwr","diff.median.wrs.ci95upr", "pvalue.wrs.2side","pvalue.wrs.1side"))
  return(stats.t.all)
  
  # tmp = stats.t.all[duplicated(stats.t.all)]
}
fn_format.stat.output = function(stats.t, met.t = 'logScore.loose', loc.t, target.t, mea.t){
  out.stats.t = stats.t %>% dplyr::filter(metric == met.t & loc %in% loc.t & target %in% target.t & measure %in% mea.t) %>% 
    dplyr::select(loc, variant, wave, measure, target, n.smp, mean.base, mean.opt, diff.mean, perc.diff, median.diff, pvalue.wrs) %>%
    rename(State = loc, Period = variant, Wave = wave, Measure = measure, Target = target, `No. prediction` = n.smp, 
           `relative differnce (%)` = perc.diff,  pvalue = pvalue.wrs) %>%
    mutate(Period = factor(Period, levels = c('non.Omicron','Omicron'), labels = c('pre-Omicron','Omicron'))) %>%
    mutate(Wave = factor(Wave, levels = c('wave2','Alpha','Delta'), labels = c('2nd wave','Alpha','Delta'))) %>%
    mutate(Target = factor(Target, levels = c('all', 'peak week', 'peak intensity','total', '1-8wk ahead','9-16wk ahead','17-26wk ahead'))) %>%
    mutate(Measure = factor(Measure, levels = c('case', 'death'), labels = c('Cases','Deaths')))  %>%
    arrange(State, Period, Wave, Measure, Target)
  return(out.stats.t)
}

# plot stats
fn_plotStat = function(da.t){
  
  p.byv.xv = ggplot(da.t %>% filter(wave %in% c('Pre-Omicron','Omicron')), 
                    aes(y = diff.median.wrs %>% as.numeric, x = wave, color = loc)) + 
    geom_pointrange(aes(ymin = diff.median.wrs.ci95lwr %>% as.numeric, ymax = diff.median.wrs.ci95upr %>% as.numeric),position = position_dodge(width = .8), fatten = .5, size = .5) +
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    facet_rep_wrap(~measure + target, scales = 'free_y', ncol = 7,
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') +
    guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t5 # + theme(legend.position = 'none')
  # p3
  
  p.byv.xv.c = ggplot(da.t %>% filter(measure == 'Cases' & wave %in% c('Pre-Omicron','Omicron')), 
                      aes(y = diff.median.wrs %>% as.numeric, x = wave, color = loc)) + 
    geom_pointrange(aes(ymin = diff.median.wrs.ci95lwr %>% as.numeric, ymax = diff.median.wrs.ci95upr %>% as.numeric),position = position_dodge(width = .8), fatten = .5, size = .5) +
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    facet_rep_wrap(~target, ncol = 7, # , scales = 'free_y'
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') + ggtitle(paste0('(A) Cases')) +
    guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t5 # + theme(legend.position = 'none')
  p.byv.xv.d = ggplot(da.t %>% filter(measure == 'Deaths' & wave %in% c('Pre-Omicron','Omicron')), 
                      aes(y = diff.median.wrs %>% as.numeric, x = wave, color = loc)) + 
    geom_pointrange(aes(ymin = diff.median.wrs.ci95lwr %>% as.numeric, ymax = diff.median.wrs.ci95upr %>% as.numeric),position = position_dodge(width = .8), fatten = .5, size = .5) +
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    facet_rep_wrap(~target, ncol = 7, # , scales = 'free_y'
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') + ggtitle(paste0('(B) Deaths')) +
    guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t5 # + theme(legend.position = 'none')
  
  p.byv.xtarget = ggplot(da.t %>% filter(wave %in% c('Pre-Omicron','Omicron')), 
                         aes(y = diff.median.wrs %>% as.numeric, x = target, color = loc)) + 
    geom_pointrange(aes(ymin = diff.median.wrs.ci95lwr %>% as.numeric, ymax = diff.median.wrs.ci95upr %>% as.numeric),position = position_dodge(width = .8), fatten = .5, size = .5) +
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    facet_rep_wrap(~measure + wave, scales = 'free_y', ncol = 2,
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') +
    guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  p.byw.xwave = ggplot(da.t %>% filter(wave %in% c('2nd wave','Alpha','Delta')), 
                       aes(y = diff.median.wrs %>% as.numeric, x = wave, color = loc)) + 
    geom_pointrange(aes(ymin = diff.median.wrs.ci95lwr %>% as.numeric, ymax = diff.median.wrs.ci95upr %>% as.numeric),position = position_dodge(width = .8), fatten = .5, size = .5) +
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    facet_rep_wrap(~measure + target, scales = 'free_y', ncol = 4,
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') +
    guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t5 # + theme(legend.position = 'none')
  # p3
  
  p.byw.xtarget = ggplot(da.t %>% filter(wave %in% c('2nd wave','Alpha','Delta')), 
                         aes(y = diff.median.wrs %>% as.numeric, x = target, color = loc)) + 
    geom_pointrange(aes(ymin = diff.median.wrs.ci95lwr %>% as.numeric, ymax = diff.median.wrs.ci95upr %>% as.numeric),position = position_dodge(width = .8), fatten = .5, size = .5) +
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    facet_rep_wrap(~measure + wave, scales = 'free_y', ncol = 3,
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') +
    guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  
  return(list(p.byv.xv = p.byv.xv, p.byv.xv.c = p.byv.xv.c, p.byv.xv.d = p.byv.xv.d, 
              p.byv.xtarget = p.byv.xtarget, p.byw.xwave = p.byw.xwave, p.byw.xtarget = p.byw.xtarget)) 
}


fn_plot.cp.sce.byloc.multisns = function(da.t, stat.t, sn.opt){
  rm(ymax.t, da.tt, stat.tt, ymax.t2, da.tt2, stat.tt2)
  col.t = 'grey50'; sz = 2.5;  col.tx = 'blue'
  da.tt = da.t %>% filter(variant %in% c('non.Omicron','Omicron')) %>%
    mutate(variant = factor(variant, levels = c('non.Omicron','Omicron'), labels = c('Pre-Omicron','Omicron')))
  stat.tt = stat.t %>% filter(seasonality == sn.opt & wave %in% c('Pre-Omicron','Omicron')) %>%
    mutate(variant = factor(variant, levels = c('non.Omicron','Omicron'), labels = c('Pre-Omicron','Omicron'))) %>%
    mutate(lab = paste0(perc.diff,'%'))
  ymax.t = max(da.tt$diff, na.rm = T)
  pp1 = ggplot(da.tt,
               aes(x = target, y = diff, color = seasonality)) +
    geom_violin(size = .2) + 
    geom_sina(aes(color = seasonality),  size = .05, shape = 20) + 
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    geom_text(aes(x = target, y = ymax.t * 1.1, group = seasonality, label = mean.base %>% as.numeric %>% round(2)), color =  col.tx, size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .95, group = seasonality, label = mean.opt %>% as.numeric %>% round(2)),  color = col.tx, size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .8, group = seasonality, label = lab), color =  col.tx, size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .65, group = seasonality, label = pvalue.wrs),  color =  col.tx, size = sz, data = stat.tt) +
    facet_rep_wrap(~measure + variant, ncol = 2, # scales = 'free_y', 
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'sn') +
    # guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  da.tt2 = da.t %>% filter(wave %in% c('2nd wave','Alpha','Delta')) 
  stat.tt2 = stat.t %>% filter(seasonality == sn.opt &wave %in% c('2nd wave','Alpha','Delta')) %>%
    mutate(lab = paste0(perc.diff,'%'))
  ymax.t2 = max(da.tt$diff, na.rm = T)
  pp2 = ggplot(da.tt2,
               aes(x = target, y = diff, color = seasonality)) +
    geom_violin(size = .2) + 
    geom_sina(aes(color = seasonality),  size = .05, shape = 20) + 
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    geom_text(aes(x = target, y = ymax.t2 * 1.1, label = mean.base %>% as.numeric %>% round(2)), color = col.tx, size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .95, label = mean.opt %>% as.numeric %>% round(2)), color =  col.tx, size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .8, label = lab), color = col.tx, size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .65, label = pvalue.wrs), color =  col.tx, size = sz, data = stat.tt2) +
    facet_rep_wrap(~measure + wave, ncol = 3, # scales = 'free_y', 
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'sn') +
    # guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  return(list(pp.byv = pp1, pp.byw = pp2))
}

fn_plot.cp.sn.byloc = function(da.t, stat.t, cp.t = 'no.v.ytf'){
  
  rm(ymax.t, da.tt, stat.tt, ymax.t2, da.tt2, stat.tt2)
  col.t = 'grey50'; sz = 2.5;  col.tx = 'darkgreen'
  da.tt = da.t %>% filter(variant %in% c('non.Omicron','Omicron')) %>%
    mutate(variant = factor(variant, levels = c('non.Omicron','Omicron'), labels = c('Pre-Omicron','Omicron')))
  stat.tt = stat.t %>% filter(cp == cp.t & wave %in% c('Pre-Omicron','Omicron')) %>%
    mutate(variant = factor(variant, levels = c('non.Omicron','Omicron'), labels = c('Pre-Omicron','Omicron'))) %>%
    mutate(lab = paste0(perc.diff,'%'))
  ymax.t = max(da.tt$diff, na.rm = T)
  # print(ymax.t)
  pp1 = ggplot(da.tt,
               aes(x = target, y = diff, color = cp)) +
    geom_violin(size = .2) + 
    geom_sina(aes(color = cp),  size = .05, shape = 20) + 
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    geom_text(aes(x = target, y = ymax.t * 1.1,  label = mean.base %>% as.numeric %>% round(2)), color =  col.tx, size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .95, label = mean.opt %>% as.numeric %>% round(2)),  color = col.tx, size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .8,  label = lab), color =  col.tx, size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .65,  label = pvalue.wrs),  color =  col.tx, size = sz, data = stat.tt) +
    facet_rep_wrap(~measure + variant, ncol = 2, # scales = 'free_y', 
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'sn') +
    # guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  da.tt2 = da.t %>% filter(wave %in% c('2nd wave','Alpha','Delta')) 
  stat.tt2 = stat.t %>% filter(cp == cp.t & wave %in% c('2nd wave','Alpha','Delta')) %>%
    mutate(lab = paste0(perc.diff,'%'))
  ymax.t2 = max(da.tt2$diff, na.rm = T)
  pp2 = ggplot(da.tt2,
               aes(x = target, y = diff, color = cp)) +
    geom_violin(size = .2) + 
    geom_sina(aes(color = cp),  size = .05, shape = 20) + 
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    geom_text(aes(x = target, y = ymax.t2 * 1.1, label = mean.base %>% as.numeric %>% round(2)), color = col.tx, size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .95, label = mean.opt %>% as.numeric %>% round(2)), color =  col.tx, size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .8, label = lab), color = col.tx, size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .65, label = pvalue.wrs), color =  col.tx, size = sz, data = stat.tt2) +
    facet_rep_wrap(~measure + wave, ncol = 3, # scales = 'free_y', 
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'sn') +
    # guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  return(list(pp.byv = pp1, pp.byw = pp2))
}

fn_plot.dist.byloc = function(da.t, stat.t){
  
  col.t = 'grey50'; sz = 2.5
  da.tt = da.t %>% filter(variant %in% c('non.Omicron','Omicron')) %>%
    mutate(variant = factor(variant, levels = c('non.Omicron','Omicron'), labels = c('Pre-Omicron','Omicron')))
  stat.tt = stat.t %>% filter(wave %in% c('Pre-Omicron','Omicron')) %>%
    mutate(variant = factor(variant, levels = c('non.Omicron','Omicron'), labels = c('Pre-Omicron','Omicron'))) %>%
    mutate(lab = paste0(perc.diff,'%'))
  ymax.t = max(da.tt$diff, na.rm = T)
  pp1 = ggplot(da.tt,
               aes(x = target, y = diff)) +
    geom_violin(color = col.t, size = .2) + 
    geom_sina(color = col.t,  size = .05, shape = 20) + 
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    geom_text(aes(x = target, y = ymax.t * 1.1, label = mean.base %>% as.numeric %>% round(2)), color = 'red', size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .95, label = mean.opt %>% as.numeric %>% round(2)), color = 'red', size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .8, label = lab), color = 'red', size = sz, data = stat.tt) +
    geom_text(aes(x = target, y = ymax.t * .65, label = pvalue.wrs), color = 'red', size = sz, data = stat.tt) +
    facet_rep_wrap(~measure + variant, ncol = 2, # scales = 'free_y', 
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') +
    # guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  da.tt2 = da.t %>% filter(wave %in% c('2nd wave','Alpha','Delta')) 
  stat.tt2 = stat.t %>% filter(wave %in% c('2nd wave','Alpha','Delta')) %>%
    mutate(lab = paste0(perc.diff,'%'))
  ymax.t2 = max(da.tt$diff, na.rm = T)
  pp2 = ggplot(da.tt2,
               aes(x = target, y = diff)) +
    geom_violin(color = col.t, size = .2) + 
    geom_sina(color = col.t,  size = .05, shape = 20) + 
    geom_hline(yintercept = 0, color = 'black', size = .2) +
    geom_text(aes(x = target, y = ymax.t2 * 1.1, label = mean.base %>% as.numeric %>% round(2)), color = 'red', size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .95, label = mean.opt %>% as.numeric %>% round(2)), color = 'red', size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .8, label = lab), color = 'red', size = sz, data = stat.tt2) +
    geom_text(aes(x = target, y = ymax.t2 * .65, label = pvalue.wrs), color = 'red', size = sz, data = stat.tt2) +
    facet_rep_wrap(~measure + wave, ncol = 3, # scales = 'free_y', 
                   repeat.tick.labels = T, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = '', y = 'Diffence in log-score', color = 'State') +
    # guides(color = guide_legend(override.aes = list(size = 0.1))) +
    theme_minimal() +theme.t6 # + theme(legend.position = 'none')
  
  return(list(pp.byv = pp1, pp.byw = pp2))
}
