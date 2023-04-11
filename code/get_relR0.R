# use weather data to represent seasonality
# 2/28/21

Rwea_parm.bounds = rbind(
  c(2.34, 2.93), # R0max
  c(.86,1.18), # R0diff - this determines the magnitude of seasonality
  # reduce the leve of seasonality, b/c cp flu corona virus seasonality is less profound
  c(2.2,4.0)/1000, # qmin
  c(17,20)/1000, # qmax
  c(10.2,11.3)/1000, # qmid
  c(20.2,24), # Tc
  c(.4, 5.1), # Tdiff
  c(.95,1.54) # Texp
)
# use the mean instead to reduce uncertainty?
Rwea_parm.bounds = cbind(rowMeans(Rwea_parm.bounds),rowMeans(Rwea_parm.bounds))
rownames(Rwea_parm.bounds) = c('R0max','R0diff','qmin','qmax','qmid','Tc','Tdiff','Texp')


fn_getRelR0 = function(loc.t, ref.wk, Rwea_parm.bounds, smooth = T){
  # loc.t = 'UK'
  da = read.csv(paste0(dir_data,'wea.by.week_',toupper(loc.t),'.csv')) %>% data.table()
  
  # num_ens = 500
  parm_Rwea = lhs(num_ens,rect = Rwea_parm.bounds) 
  colnames(parm_Rwea) = rownames(Rwea_parm.bounds)
  
  # calculate R0
  calc_R0wea <- function(da_wea, parm_Rwea) {
    
    spec.hum = da_wea$spec.hum;
    temp = da_wea$temp
    num_ens = nrow(parm_Rwea)
    
    # Create matrices for storing/returning results:
    res.temp = res.temp.red = matrix(0, nrow(da_wea), num_ens)
    
    # Loop through all ensemble members:
    for (ix in 1:num_ens) {
      
      # Assign parameters:
      q.mn <- parm_Rwea[ix, 'qmin']; q.mx <- parm_Rwea[ix, 'qmax']; q.md <- parm_Rwea[ix, 'qmid']
      R0.max <- parm_Rwea[ix, 'R0max']; R0.diff <- parm_Rwea[ix, 'R0diff']
      Tc <- parm_Rwea[ix, 'Tc']; Tdiff <- parm_Rwea[ix, "Tdiff"]; t.exp <- parm_Rwea[ix, 'Texp']
      
      
      q.mn.cut <- q.mn
      
      # Calculate and correct R0.min
      R0.min <- R0.max - R0.diff
      if (R0.min < 0) {
        R0.min <- 0.1
      }
      Tmin <- Tc - Tdiff
      
      # Calculate parabola params:
      
      # given the symmetry:
      q.mx.left = 2 * q.md - q.mn; 
      b.left <- ((R0.max - R0.min) * (q.mx.left + q.mn)) / ((q.mx.left - q.md) * (q.mn - q.md))
      a.left <- (-1 * b.left) / (q.mx.left + q.mn)
      c.left <- R0.min - a.left * q.md ** 2 - b.left * q.md
      
      q.mn.right = 2 * q.md - q.mx
      b.right <- ((R0.max - R0.min) * (q.mx + q.mn.right)) / ((q.mx - q.md) * (q.mn.right - q.md))
      a.right <- (-1 * b.right) / (q.mx + q.mn.right)
      c.right <- R0.min - a.right * q.md ** 2 - b.right * q.md
      
      fit1 = fit2 =numeric(nrow(da_wea))
      
      # split the data into two sets (those >=q.md, and those <q.md)
      idx.left = which(spec.hum < q.md); idx.right = which(spec.hum >= q.md)
      
      # Full model:
      q1.left <- spec.hum[idx.left]; q1.right = spec.hum[idx.right]
      t1.left <- temp[idx.left]; t1.right = temp[idx.right]
      fit1[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
      fit1[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
      
      # Reduced model:
      q1 <- spec.hum
      q1[q1 < q.mn.cut] <- q.mn.cut; q1[q1 > q.mx] <- q.mx
      t1 <- temp; t1[t1 < Tmin] <- Tmin
      q1.left <- q1[idx.left]; q1.right = q1[idx.right]
      t1.left <- t1[idx.left]; t1.right = t1[idx.right]
      fit2[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
      fit2[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
      
      
      # Store results in matrices:
      res.temp[, ix] <- fit1; res.temp.red[, ix] <- fit2
    }
    
    # Return results:
    return(list(res.temp, res.temp.red))
  }
  
  tmp = calc_R0wea(da_wea = da, parm_Rwea)
  relR0 = tmp[[2]] 
  
  # relR0 = relR0 / matrix(colMeans(relR0[ref.wk + (-1:1),]), nrow(relR0), num_ens, byrow=T) # relative to to spring
  # relative to the mean?
  relR0 = relR0 / matrix(colMeans(relR0[-53,]), nrow(relR0), num_ens, byrow=T) # relative to to spring
  # smooth the curve
  relR0 = relR0 %>% apply(2, stats::filter, filter = rep(1/3,3), circular = T)
  
  relR0
}


# matplot(relR0, type = 'l', lty=1)  
# lines(relR0 %>% rowMeans(), lwd=2)
  
  
# calculate R0
calc_R0wea <- function(da_wea, parm_Rwea, smooth = T) {
  
  spec.hum = da_wea$spec.hum;
  temp = da_wea$temp
  num_ens = nrow(parm_Rwea)
  
  # Create matrices for storing/returning results:
  res.temp = res.temp.red1 = res.temp.red2 = matrix(0, nrow(da_wea), num_ens)
  
  # Loop through all ensemble members:
  for (ix in 1:num_ens) {
    
    # Assign parameters:
    q.mn <- parm_Rwea[ix, 'qmin']; q.mx <- parm_Rwea[ix, 'qmax']; q.md <- parm_Rwea[ix, 'qmid']
    R0.max <- parm_Rwea[ix, 'R0max']; R0.diff <- parm_Rwea[ix, 'R0diff']
    Tc <- parm_Rwea[ix, 'Tc']; Tdiff <- parm_Rwea[ix, "Tdiff"]; t.exp <- parm_Rwea[ix, 'Texp']
    
    
    q.mn.cut <- q.mn
    
    # Calculate and correct R0.min
    R0.min <- R0.max - R0.diff
    if (R0.min < 0) {
      R0.min <- 0.1
    }
    Tmin <- Tc - Tdiff
    
    # Calculate parabola params:
    
    # given the symmetry:
    q.mx.left = 2 * q.md - q.mn; 
    b.left <- ((R0.max - R0.min) * (q.mx.left + q.mn)) / ((q.mx.left - q.md) * (q.mn - q.md))
    a.left <- (-1 * b.left) / (q.mx.left + q.mn)
    c.left <- R0.min - a.left * q.md ** 2 - b.left * q.md
    
    q.mn.right = 2 * q.md - q.mx
    b.right <- ((R0.max - R0.min) * (q.mx + q.mn.right)) / ((q.mx - q.md) * (q.mn.right - q.md))
    a.right <- (-1 * b.right) / (q.mx + q.mn.right)
    c.right <- R0.min - a.right * q.md ** 2 - b.right * q.md
    
    fit1 = fit2 = fit3 =numeric(nrow(da_wea))
    
    # split the data into two sets (those >=q.md, and those <q.md)
    idx.left = which(spec.hum < q.md); idx.right = which(spec.hum >= q.md)
    
    # Full model:
    q1.left <- spec.hum[idx.left]; q1.right = spec.hum[idx.right]
    t1.left <- temp[idx.left]; t1.right = temp[idx.right]
    fit1[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit1[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    # Reduced model:
    q1 <- spec.hum
    q1[q1 < q.mn.cut] <- q.mn.cut; q1[q1 > q.mx] <- q.mx
    t1 <- temp; t1[t1 < Tmin] <- Tmin
    q1.left <- q1[idx.left]; q1.right = q1[idx.right]
    t1.left <- t1[idx.left]; t1.right = t1[idx.right]
    fit2[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit2[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    # Reduced temp model (keep humidity:
    q1 <- spec.hum
    # q1[q1 < q.mn.cut] <- q.mn.cut; q1[q1 > q.mx] <- q.mx
    t1 <- temp; t1[t1 < Tmin] <- Tmin
    q1.left <- q1[idx.left]; q1.right = q1[idx.right]
    t1.left <- t1[idx.left]; t1.right = t1[idx.right]
    fit3[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit3[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    # Store results in matrices:
    res.temp[, ix] <- fit1; res.temp.red2[, ix] <- fit2; res.temp.red1[, ix] <- fit3
    
  }
  if(smooth){
    res.temp = res.temp %>% apply(2, stats::filter, filter = rep(1/3,3), circular = T)
    res.temp.red2 = res.temp.red2 %>% apply(2, stats::filter, filter = rep(1/3,3), circular = T)
    res.temp.red1 = res.temp.red1 %>% apply(2, stats::filter, filter = rep(1/3,3), circular = T)
  }
  # Return results:
  return(list(res.temp, res.temp.red2, res.temp.red1))
}


# by state
fn_getRelR0.loc = function(da.t, ref.wk, Rwea_parm.bounds, smooth = T){
  # loc.t = 'UK'
  da = da.t %>% data.table()
  
  # num_ens = 500
  parm_Rwea = lhs(num_ens,rect = Rwea_parm.bounds) 
  colnames(parm_Rwea) = rownames(Rwea_parm.bounds)
  
  
  tmp = calc_R0wea(da_wea = da, parm_Rwea)
  relR0 = tmp[[2]] 
  
  # relR0 = relR0 / matrix(colMeans(relR0[ref.wk + (-1:1),]), nrow(relR0), num_ens, byrow=T) # relative to to spring
  # relative to the mean?
  relR0 = relR0 / matrix(colMeans(relR0[-53,,drop=F]), nrow(relR0), num_ens, byrow=T) # relative to to spring
  # smooth the curve
  relR0 = relR0 %>% apply(2, stats::filter, filter = rep(1/3,3), circular = T)
  
  relR0
}


# function to tranform the seasonality trend
fn_R0trans = function(tda, 
                      relR0min = .9, # lowest value for relR0
                      dur0 = 15, # duration/number of weeks relR0>1
                      pshift = 0 # small shift in peak timing
){
  # save the orignal
  tda0 = copy(tda)
  setnames(tda0, 'x', 'wk')
  tda = tda[x %in% 1:52] # %>% filter(x %in% 1:52)
  
  imax = which.max(tda$y)
  idx.high = which(tda$y > 1)
  # break it down to different segments
  # test whether there are multiple peak
  # (idx.high + 52) %% 52
  dif = idx.high[-length(idx.high)] - idx.high[-1]
  iseg = which(dif < -1)
  segs = NULL;
  if(length(iseg) > 0){
    for(i in 1:length(iseg)){
      if(i==1){
        segs[[i]] = idx.high[1]:idx.high[iseg[i]]
      } else {
        segs[[i]] = idx.high[iseg[i-1]+1]:(idx.high[iseg[i]]) # 
      }
    }
    # add the last one
    segs[[i+1]] = idx.high[iseg[i]+1]:tail(idx.high,1)
    # combine first and last if they are adjacent
    if(segs[[1]][1]==1 & tail(segs[[i+1]],1) == 52){
      segs[[1]] = c(segs[[1]], segs[[i+1]])
      segs[[i+1]] = NULL
    }
  } else {
    segs[[1]] = idx.high
  }
  
  fit.high.tran = NULL
  for(iseg in 1:length(segs)){
    idx.t = segs[[iseg]]
    # b/c it's part of the full year - need to account for the weeks not included?
    # imax = which.max(tda[idx.t]$y) + idx.t[1] - 1
    imax = which(tda$y == max(tda[idx.t]$y))
    ishift = (idx.t - imax)  
    i1 = which(ishift > 26); new1 = ishift[ishift > 26]  - 52
    i2 = which(ishift < -26); new2 = ishift[ishift < -26]  + 52
    ishift[i1] = new1; ishift[i2] = new2
    # ishift[ishift < -26] = ishift[ishift < -26]  + 52
    vtran = (length(idx.high) / dur0) %>% pmin(2) # up to 2, this include all peaks
    itran = ishift / vtran
    # shift the timing of the peak
    # imax.new = (imax + round(pshift / vtran)) %% 52
    imax.new = (imax + round(pshift)) %% 52
    # tda.high = data.table(x = ishift,y = tda$y[idx.high])
    # transform the x axis, then use the model to map to week
    tda.high = data.table(x = itran,y = tda$y[idx.t])
    fit.high = lm(y~poly(x, 5, raw = T), data = tda.high)
    fit.high.tran.t = data.table(x = round(itran) %>% unique)
    fit.high.tran.t$y = 0
    fit.high.tran.t$y = predict(fit.high, newdata = fit.high.tran.t)
    # convert back to the calendar week
    rishift = fit.high.tran.t$x
    # back to the center
    # rishift[rishift<0] = rishift[rishift<0] + imax 
    # rishift[rishift>0] = (rishift[rishift>0] + imax)  %% 52
    # don't overwrite it yet
    i1 = which(rishift<0); new1 = (rishift[rishift<0] + imax.new) %% 52
    i2 = which(rishift>0); new2 = (rishift[rishift>0] + imax.new)  %% 52
    rishift[rishift==0] = imax.new
    rishift[i1] = new1
    rishift[i2] = new2
    rishift[rishift==0] = 52
    
    # rishift
    fit.high.tran.t$wk = rishift
    
    fit.high.tran = rbind(fit.high.tran, fit.high.tran.t)
  }
  
  idx.low = which(tda$y <= 1)
  # break it down to different segments
  # (idx.low + 52) %% 52
  dif = idx.low[-length(idx.low)] - idx.low[-1]
  iseg = which(dif < -1)
  segs = NULL;
  if(length(iseg) > 0){
    
    for(i in 1:length(iseg)){
      if(i==1){
        segs[[i]] = idx.low[1]:idx.low[iseg[i]]
      } else {
        segs[[i]] = idx.low[iseg[i-1]+1]:(idx.low[iseg[i]]) # 
      }
    }
    # add the last one
    segs[[i+1]] = idx.low[iseg[i]+1]:tail(idx.low,1)
    # combine first and last if they are adjacent
    if(segs[[1]][1]==1 & tail(segs[[i+1]],1) == 52){
      segs[[1]] = c(segs[[1]], segs[[i+1]])
      segs[[i+1]] = NULL
    }
  } else {
    segs[[1]] = idx.low
  }
  
  fit.low.tran = NULL
  for(iseg in 1:length(segs)){
    idx.t = segs[[iseg]]
    
    # imin = which.min(tda[idx.t]$y) + idx.t[1] - 1
    imin = which(tda$y == min(tda[idx.t]$y))
    ishift = (idx.t - imin)  
    ishift[ishift > 26] = ishift[ishift > 26]  - 52
    ishift[ishift < -26] = ishift[ishift < -26]  + 52
    vtran = (52 - nrow(fit.high.tran))/length(idx.low) # up to 2, include all
    itran = ishift * vtran  # stretch out
    # shift the trough timing
    imin.new = imin + round(pshift) # round(pshift * vtran)
    # tda.high = data.table(x = ishift,y = tda$y[idx.high])
    # transform the x axis, then use the model to map to week
    tda.low = data.table(x = itran,y = tda$y[idx.t])
    fit.low = lm(y~poly(x, 5, raw = T), data = tda.low)
    # fit.low.tran = data.table(x = min(round(itran,0)): max(round(itran,0)) %>% unique)
    range.low = range(round(itran)%>% unique)
    fit.low.tran.t = data.table(x =  range.low[1] : range.low[2])
    fit.low.tran.t$y = 0
    fit.low.tran.t$y = predict(fit.low, newdata = fit.low.tran.t)
    # convert back to the calendar week
    rishift = fit.low.tran.t$x
    i1 = which(rishift<0); new1 = (rishift[rishift<0] + imin.new) %% 52
    i2 = which(rishift>0); new2 = (rishift[rishift>0] + imin.new) %% 52
    rishift[rishift==0] = imin.new
    rishift[i1] = new1
    rishift[i2] = new2
    rishift[rishift==0] = 52
    rishift
    fit.low.tran.t$wk = rishift
    
    fit.low.tran = rbind(fit.low.tran, fit.low.tran.t)
  }
  
  fit.low.tran$y.tran = 1 - (1-fit.low.tran$y) * pmin(1, min(fit.low.tran$y) / relR0min)
  
  fit.tran = rbind(fit.high.tran[,c('wk','y')], fit.low.tran[,c('wk','y.tran')],use.names=F) %>% .[order(wk)] # arrange(wk)
  # there may be overlaps
  fit.tran = fit.tran[, list(y = mean(y)), by = 'wk']
  fit = lm(y~poly(wk, 5, raw = T), data = fit.tran)
  fit.tran.full = data.table(wk = 1:52)
  fit.tran.full$y = 0
  fit.tran.full$y = predict(fit, newdata = fit.tran.full)
  
  fit.tran.full$y = fit.tran.full$y / mean(fit.tran.full$y)
  fit.tran.full$yma = fit.tran.full$y %>% stats::filter(., filter = rep(1/3,3), circular = T) %>% as.numeric
  
  # add week 53
  fit.tran.full = rbind(fit.tran.full, 
                        data.table(wk = 53, 
                                   y = mean(fit.tran.full$y[c(1,52)]),
                                   yma = mean(fit.tran.full$yma[c(1,52)])))
  fit.tran.full = merge(fit.tran.full, tda0[,c('wk','yo')], by = 'wk')
  
  fit.tran.full
}


# function to tranform the seasonality trend
fn_R0trans_allowtie = function(tda, 
                      relR0min = .9, # lowest value for relR0
                      dur0 = 15, # duration/number of weeks relR0>1
                      pshift = 0 # small shift in peak timing
){
  # save the orignal
  tda0 = copy(tda)
  setnames(tda0, 'x', 'wk')
  tda = tda[x %in% 1:52] # %>% filter(x %in% 1:52)
  
  y.max = tda$y %>% max # 8/30/22 to make sure the maximum is not too extreme
  y.raw.max = tda$yraw %>% max
  
  imax = which.max(tda$y)
  idx.high = which(tda$y > 1)
  # break it down to different segments
  # test whether there are multiple peak
  # (idx.high + 52) %% 52
  dif = idx.high[-length(idx.high)] - idx.high[-1]
  iseg = which(dif < -1)
  segs = NULL;
  if(length(iseg) > 0){
    for(i in 1:length(iseg)){
      if(i==1){
        segs[[i]] = idx.high[1]:idx.high[iseg[i]]
      } else {
        segs[[i]] = idx.high[iseg[i-1]+1]:(idx.high[iseg[i]]) # 
      }
    }
    # add the last one
    segs[[i+1]] = idx.high[iseg[i]+1]:tail(idx.high,1)
    # combine first and last if they are adjacent
    if(segs[[1]][1]==1 & tail(segs[[i+1]],1) == 52){
      segs[[1]] = c(segs[[1]], segs[[i+1]])
      segs[[i+1]] = NULL
    }
  } else {
    segs[[1]] = idx.high
  }
  
  fit.high.tran = NULL
  for(iseg in 1:length(segs)){
    idx.t = segs[[iseg]]
    # b/c it's part of the full year - need to account for the weeks not included?
    # imax = which.max(tda[idx.t]$y) + idx.t[1] - 1
    imax = which(tda$y == max(tda[idx.t]$y))
    
    # if the relR0 passed is capped (per T/AH), then there may be multiple maximum
    # order those by timing, and get the middle pt?
    # based the timing on the un-capped relR0
    if(length(imax) > 1)
      imax =  which(tda$yraw == max(tda[idx.t]$yraw))
    
    ishift = (idx.t - imax)  
    i1 = which(ishift > 26); new1 = ishift[ishift > 26]  - 52
    i2 = which(ishift < -26); new2 = ishift[ishift < -26]  + 52
    ishift[i1] = new1; ishift[i2] = new2
    # ishift[ishift < -26] = ishift[ishift < -26]  + 52
    vtran = (length(idx.high) / dur0) %>% pmin(2) # up to 2, this include all peaks
    itran = ishift / vtran
    # shift the timing of the peak
    # imax.new = (imax + round(pshift / vtran)) %% 52
    imax.new = (imax + round(pshift)) %% 52
    # tda.high = data.table(x = ishift,y = tda$y[idx.high])
    # transform the x axis, then use the model to map to week
    tda.high = data.table(x = itran,y = tda$y[idx.t])
    fit.high = lm(y~poly(x, 5, raw = T), data = tda.high)
    fit.high.tran.t = data.table(x = round(itran) %>% unique)
    fit.high.tran.t$y = 0
    fit.high.tran.t$y = predict(fit.high, newdata = fit.high.tran.t)
    # convert back to the calendar week
    rishift = fit.high.tran.t$x
    # back to the center
    # rishift[rishift<0] = rishift[rishift<0] + imax 
    # rishift[rishift>0] = (rishift[rishift>0] + imax)  %% 52
    # don't overwrite it yet
    i1 = which(rishift<0); new1 = (rishift[rishift<0] + imax.new) %% 52
    i2 = which(rishift>0); new2 = (rishift[rishift>0] + imax.new)  %% 52
    rishift[rishift==0] = imax.new
    rishift[i1] = new1
    rishift[i2] = new2
    rishift[rishift==0] = 52
    
    # rishift
    fit.high.tran.t$wk = rishift
    
    fit.high.tran = rbind(fit.high.tran, fit.high.tran.t)
  }
  
  idx.low = which(tda$y <= 1)
  # break it down to different segments
  # (idx.low + 52) %% 52
  dif = idx.low[-length(idx.low)] - idx.low[-1]
  iseg = which(dif < -1)
  segs = NULL;
  if(length(iseg) > 0){
    
    for(i in 1:length(iseg)){
      if(i==1){
        segs[[i]] = idx.low[1]:idx.low[iseg[i]]
      } else {
        segs[[i]] = idx.low[iseg[i-1]+1]:(idx.low[iseg[i]]) # 
      }
    }
    # add the last one
    segs[[i+1]] = idx.low[iseg[i]+1]:tail(idx.low,1)
    # combine first and last if they are adjacent
    if(segs[[1]][1]==1 & tail(segs[[i+1]],1) == 52){
      segs[[1]] = c(segs[[1]], segs[[i+1]])
      segs[[i+1]] = NULL
    }
  } else {
    segs[[1]] = idx.low
  }
  
  fit.low.tran = NULL
  for(iseg in 1:length(segs)){
    idx.t = segs[[iseg]]
    
    # imin = which.min(tda[idx.t]$y) + idx.t[1] - 1
    imin = which(tda$y == min(tda[idx.t]$y))
    ishift = (idx.t - imin)  
    ishift[ishift > 26] = ishift[ishift > 26]  - 52
    ishift[ishift < -26] = ishift[ishift < -26]  + 52
    vtran = (52 - nrow(fit.high.tran))/length(idx.low) # up to 2, include all
    itran = ishift * vtran  # stretch out
    # shift the trough timing
    imin.new = imin + round(pshift) # round(pshift * vtran)
    # tda.high = data.table(x = ishift,y = tda$y[idx.high])
    # transform the x axis, then use the model to map to week
    tda.low = data.table(x = itran,y = tda$y[idx.t])
    fit.low = lm(y~poly(x, 5, raw = T), data = tda.low)
    # fit.low.tran = data.table(x = min(round(itran,0)): max(round(itran,0)) %>% unique)
    range.low = range(round(itran)%>% unique)
    fit.low.tran.t = data.table(x =  range.low[1] : range.low[2])
    fit.low.tran.t$y = 0
    fit.low.tran.t$y = predict(fit.low, newdata = fit.low.tran.t)
    # convert back to the calendar week
    rishift = fit.low.tran.t$x
    i1 = which(rishift<0); new1 = (rishift[rishift<0] + imin.new) %% 52
    i2 = which(rishift>0); new2 = (rishift[rishift>0] + imin.new) %% 52
    rishift[rishift==0] = imin.new
    rishift[i1] = new1
    rishift[i2] = new2
    rishift[rishift==0] = 52
    rishift
    fit.low.tran.t$wk = rishift
    
    fit.low.tran = rbind(fit.low.tran, fit.low.tran.t)
  }
  
  fit.low.tran$y.tran = 1 - (1-fit.low.tran$y) * pmin(1, min(fit.low.tran$y) / relR0min)
  
  fit.tran = rbind(fit.high.tran[,c('wk','y')], fit.low.tran[,c('wk','y.tran')],use.names=F) %>% .[order(wk)] # arrange(wk)
  # there may be overlaps
  fit.tran = fit.tran[, list(y = mean(y)), by = 'wk']
  fit = lm(y~poly(wk, 5, raw = T), data = fit.tran)
  fit.tran.full = data.table(wk = 1:52)
  fit.tran.full$y = 0
  fit.tran.full$y = predict(fit, newdata = fit.tran.full)
  
  # 8/30/22 - cap y at the maximum of yraw to make sure it is not artificially increased?
  fit.tran.full = merge(fit.tran.full, fit.tran, by = 'wk', all.x = T, suffixes = c('','.before'))
  fit.tran.full = fit.tran.full %>% mutate(y = case_when(y > y.before ~ y.before, 
                                                                 T ~ y)) 
  fit.tran.full$y.before = NULL
  
  fit.tran.full$y = fit.tran.full$y / mean(fit.tran.full$y)
  fit.tran.full$yma = fit.tran.full$y %>% stats::filter(., filter = rep(1/3,3), circular = T) %>% as.numeric
  
  # add week 53
  fit.tran.full = rbind(fit.tran.full, 
                        data.table(wk = 53, 
                                   y = mean(fit.tran.full$y[c(1,52)]),
                                   yma = mean(fit.tran.full$yma[c(1,52)])))
  fit.tran.full = merge(fit.tran.full, tda0[,c('wk','yo')], by = 'wk')
  
  fit.tran.full
}

fn_getR0trans = function(relR0, num_ens, relR0uncapped = NULL, 
                         bound.relR0min = c(.8, 1.2), # bound.relR0min = c(.8, 1.1), # 
                         bound.dur0 = c(15,20), # bound.dur0 = c(15,25), #  # bound.dur0 = c(12,18), 
                         bound.pshift=c(-2, 2)){
  if(is.null(relR0uncapped))
    relR0uncapped = relR0
  tda = data.table(x = c(1:53), yo = c(relR0[1:53]), 
                   y = c(relR0[1:53]) %>% stats::filter(., filter = rep(1/3,3), circular = T) %>% as.numeric,
                   yraw = relR0uncapped %>% stats::filter(., filter = rep(1/3,3), circular = T) %>% as.numeric # to allow ties [used when there are multiple maxima], use moving average so it'd be more stable
                   )
  nhigh = which(tda[1:52]$y > 1) %>% length()
  relR0min = min(tda[1:52]$y)
  ymin1 = pmax(1, 1/min(relR0) * .7)
  ymin2 = 1/min(relR0)
  parm_vec = lhs(num_ens, rect = rbind(# bound.relR0min, # to adjust the minial
                                      # c(pmin(relR0min*1.1, bound.relR0min[1]), pmin(relR0min * 1.5, bound.relR0min[2])),
                                      c(relR0min, pmin(relR0min * 2, bound.relR0min[2])), # no shrinking to 50% shrinking
                                       # c(pmin(round(nhigh/2),bound.dur0[1]), pmin(round(nhigh*.7), bound.dur0[2])), # for the number of week with relR0 >1
                                      c(pmax(12, pmin(round(nhigh * .6),bound.dur0[1])), pmin(round(nhigh*.9), bound.dur0[2])), 
                                      bound.pshift)) # shift in peak/trough timing, likely earlier than the flu
  colnames(parm_vec) = c('relR0min','dur0','pshift')
  fit.tran = matrix(0, nrow = 53, ncol = num_ens)
  for(i in 1:num_ens){
    tmp = fn_R0trans_allowtie(tda, relR0min=parm_vec[i,'relR0min'], dur0=parm_vec[i,'dur0'], pshift = parm_vec[i,'pshift'])
    fit.tran[,i] = tmp$yma
  }
  fit.tran
}
