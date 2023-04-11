# get the prob distribution for the forecast
# 2/23/22

# bins for case: before the omicron wave, most peak < 1% of the population, so bin up to 1%, .05% per bin
N = 1e6
bins.case = c(seq(0, 1, by = .05) / 100 * N, N)
bins.death = c(seq(0, 2, by = .1) / 100 / 100 * N, N)

bins.tot.case = c(seq(0, 10, by = 2)/100, seq(15, 50, by = 5)/100, 1) * N # for the cumulative
bins.tot.death = c(seq(0, 10, by = 2)/100 / 100, seq(15, 50, by = 5)/100 / 100, 1) * N # for the cumulative

fn_getProbDist = function(fcast.case, fcast.death, bins.case, bins.death, fcast.hosp = NULL, bins.hosp = bins.hosp){
  num_wk.fast = fcast.case %>% nrow
  num_ens = fcast.case %>% ncol
  
  meas = c('case', 'death')
  if(! is.null(fcast.hosp)){
    meas = c('case', 'hosp', 'death')
  }
  
  # find prob. distributions for peak weeks, peak intensities and forecasts for the next n weeks
  ProbDist = NULL
  for(ftype in meas){
    
    fcast.t = get(paste0('fcast.', ftype))
    bins.t = get(paste0('bins.', ftype))
    bins.tot.t = get(paste0('bins.tot.', ftype))
    
    peakWeeks.t= fcast.t %>% apply(2, which.max)
    peakIntensities.t = fcast.t %>% apply(2, max)
    
    totals.t = fcast.t %>% apply(2, sum) # sum over the fcast period
    
    peakWeeksDist.t = matrix(NA, nrow=length(unique(peakWeeks.t)), ncol=3)
    row=1
    for(i in sort(unique(peakWeeks.t))){
      peakWeeksDist.t[row, 1] = i;
      peakWeeksDist.t[row, 2] = i;
      peakWeeksDist.t[row, 3] = round(length(peakWeeks.t[peakWeeks.t==i])/length(peakWeeks.t), 4)
      row=row+1
    } 
    
    peakIntensitiesDist.t  = matrix(NA, nrow=length(bins.t)-1, ncol=3)
    row=1
    for(i in 2:length(bins.t)){
      peakIntensitiesDist.t [row, 1] = bins.t[i-1]
      peakIntensitiesDist.t [row, 2] = bins.t[i]
      peakIntensitiesDist.t [row, 3] = round(length(peakIntensities.t[peakIntensities.t >= bins.t[i-1] & peakIntensities.t < bins.t[i]])/length(peakIntensities.t), 4)
      row=row+1    
    }
    
    totalsDist.t  = matrix(NA, nrow=length(bins.tot.t)-1, ncol=3)
    row=1
    for(i in 2:length(bins.tot.t)){
      totalsDist.t [row, 1] = bins.tot.t[i-1]
      totalsDist.t [row, 2] = bins.tot.t[i]
      totalsDist.t [row, 3] = round(length(totals.t[totals.t >= bins.tot.t[i-1] & totals.t < bins.tot.t[i]])/length(totals.t), 4)
      row=row+1    
    }
    
    # calculate prob. distribution for next n week forecasts
    nextDist.t = matrix(NA, nrow=length(bins.t)-1, ncol=nrow(fcast.t)+2)
    row=1
    for(i in 2:length(bins.t)){
      nextDist.t[row, 1] = bins.t[i-1] # lower bound
      nextDist.t[row, 2] = bins.t[i] # upper bound
      for(j in 1:nrow(fcast.t)){
        values = fcast.t[j,]
        values = values[!is.na(values)]
        nextDist.t[row, j+2] = round(length(values[values >= bins.t[i-1] & values < bins.t[i]]) / length(values), 4)
      }
      row=row+1
    }
    
    colnames(peakWeeksDist.t) = c('lwr','upr', 'prob')
    colnames(peakIntensitiesDist.t) = c('lwr','upr', 'prob')
    colnames(totalsDist.t) = c('lwr','upr', 'prob')
    colnames(nextDist.t) = c('lwr','upr', paste0('w', 1:num_wk.fast))
    nextDist.t = melt(nextDist.t %>% data.table(), id.vars = c('lwr','upr')) %>% setnames(c('variable','value'),c('target','prob'))
    nextDist.t$target = factor(nextDist.t$target, levels = paste0('w', 1:num_wk.fast), labels = paste0(1:num_wk.fast,'wk ahead'))
    peakWeeksDist.t = peakWeeksDist.t %>% data.table()
    peakIntensitiesDist.t = peakIntensitiesDist.t %>% data.table()
    totalsDist.t = totalsDist.t %>% data.table()
    peakWeeksDist.t$target = 'peak week'
    peakIntensitiesDist.t$target = 'peak intensity'
    totalsDist.t$target = 'total'
    
    ProbDist = rbind(ProbDist, data.table(data.type = ftype, rbind(peakWeeksDist.t, peakIntensitiesDist.t, totalsDist.t, nextDist.t)))
    
  }
  
  ProbDist
}

