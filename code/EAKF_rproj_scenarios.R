## function to generate the forecasts
## use EAKF for model training first
## then generate 2 types of forecasts: 
# 1) asIs: assuming no new variants emerging; 2) newV: assuming continued emergence of new variants  

# note: 
# 1. The SEIRSVimmLoss was used here, it only included 2 doses of vaccinations modeled here. 
# For later phase, combined doses 1&2 as 1st dose, and booster as 2nd dose in the model
# 2. As noted in the study, to reduce uncertainty and for comparison of approaches, 
#   the retrospective forecasts used available estimates/data during the forecast period
#   If such estimates/data are not available, it will generate real-time forecasts

library("truncnorm"); library("tgp"); library('mvtnorm'); # for lhs
library("MASS"); # for multivariate normal distribution

# supporting functions
{
  Fn_checkDA<-function(xnew,bound.low,bound.up){
    b.low=bound.low;
    b.up=bound.up;
    n.var=nrow(xnew); n.ens=ncol(xnew);
    for(vi in 1:n.var){
      #  Corrects if <b.low
      ug=min(xnew[vi,]);
      if (ug<b.low[vi]){  
        for (jj in 1:n.ens){
          if (xnew[vi,jj]<b.low[vi]){
            # xnew[vi,jj]=b.low[vi];
            xnew[vi,jj]=pmax(b.low[vi],runif(1, min=pmax(b.low[vi],quantile(xnew[vi,],.25)), max=pmax(b.low[vi],quantile(xnew[vi,],.75)))); # biased high
          }
        }
      }
      ug=max(xnew[vi,]);
      if (ug>b.up[vi]){  
        for (jj in 1:n.ens){
          if (xnew[vi,jj]>b.up[vi]){
            # xnew[vi,jj]=b.up[vi];
            # apply to non-S variables
            if(grepl('S', names(b.up[vi]))){
              xnew[vi,jj]=b.up[vi];
            } else {
              xnew[vi,jj]=runif(1, min=min(b.up[vi]/2,quantile(xnew[vi,],.5)), max=min(quantile(xnew[vi,],.80),b.up[vi]));
            }
            
          }
        }
      }
    }
    xnew;
  }
  
  Fn_getR0_SEIR=function(PARMS){
    with(as.list(PARMS), {
      Ro = beta.mean * Tir.mean
      
      Ro
    })
  }
  
  Fn_getRt_SEIR=function(PARMS){
    with(as.list(PARMS), {
      Rt = beta.mean * Tir.mean * S/N
      
      Rt
    })
  }
  
  fn_getRelMob = function(rel.mob.t, p.mob.t){ # return the scaled moblity for adjusting tx
    (rel.mob.t * p.mob.t) %>% pmin(1.5) # make sure it's <=1
  }
  
  fn_getImmLoss = function(N, S.t, E.t, I.t, tmstep, Trs.t, ts.ImmLoss.t, 
                           cum.ImmLoss.t){
    # Trs.t = mean(state0['Trs',])
    # ts.ImmLoss[tt] = (N - Spost[tt] - Epost[tt] - Ipost[tt])/Trs.t * tmstep
    ts.ImmLoss.t = ts.ImmLoss.t %>% unlist
    ts.ImmLoss = append(ts.ImmLoss.t, 
                        (N - S.t - E.t - I.t)/Trs.t * tmstep
    )
    
    wk.immloss = pmin(2*52, round(Trs.t / 7 * .75, 0)) # go back half of the immunity period? but no more than 2 years
    
    t.end = length(ts.ImmLoss)
    cum.ImmLoss = sum(ts.ImmLoss[pmax(1, t.end-wk.immloss):t.end])  # wk.immloss: how far do we wanna go back?
    # print(paste('% cum.ImmLoss:', round(cum.ImmLoss / N * 100, 2)))
    
    return(list(ts.ImmLoss = ts.ImmLoss, cum.ImmLoss = cum.ImmLoss))
  }
  
  getSRDAbounds = function(date.t, # date for the current week
                           SRbounds0, # original bounds
                           DAbounds0, # original bounds
                           parm.bound_vec  # changes by period
                           # percSR.normal = .1,
                           # percSR.extra = .2
  ){
    
    # potential updates needed, alpha, ifr, beta
    SRbounds.t = SRbounds0
    DAbounds.t = DAbounds0
    
    parm.bound_vec$date.start = parm.bound_vec$date.start %>% as.Date()
    parm.bound_vec$date.end = parm.bound_vec$date.end %>% as.Date()
    date.t = date.t  %>% as.Date()
    percSR.t = NULL
    for(type.t in c('SR', 'DA')){
      for(parm.t in c('alpha', 'ifr', 'beta')){
        tmp = parm.bound_vec[parm == parm.t & type == type.t]
        
        i.t = which(tmp$date.start <= date.t & tmp$date.end >= date.t)
        eval(parse(text = paste(type.t,'bounds.t["',parm.t,'",]= unlist(tmp[',i.t,', c("lwr","upr"),with=F])', sep = '')))
        
        if(type.t == 'SR'){
          
          if(tmp[i.t]$SR.level == 'extra.first8wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 8*7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'extra.first7wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 7*7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'extra.first6wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 6*7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'extra.first5wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 5*7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'extra.first4wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 28, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'extra.first3wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 21, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'extra.first2wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 14, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'extra'))
          } else if(tmp[i.t]$SR.level == 'Extra.first4wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 28, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra'))
          } else if(tmp[i.t]$SR.level == 'Extra.first3wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 21, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra'))
          } else if(tmp[i.t]$SR.level == 'Extra.first2wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 14, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra'))
          } else if(tmp[i.t]$SR.level == 'Extra.first1wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra'))
          } else if(tmp[i.t]$SR.level == 'Extra1.first1wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra1'))
          } else if(tmp[i.t]$SR.level == 'Extra2.first1wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra2'))
          } else if(tmp[i.t]$SR.level == 'Extra3.first1wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'Extra3'))
          } else if (tmp[i.t]$SR.level == 'EXTRA.first1wk' & date.t %in% seq(tmp[i.t]$date.start, length.out = 7, by = 'day')){
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.extra))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'EXTRA'))
          } else {
            # percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = percSR.normal))
            percSR.t = rbind(percSR.t, data.table(parm = parm.t, percSR = 'normal'))
          }
        } # end SR
      }
    }
    
    return(list(SRbounds.t = SRbounds.t, DAbounds.t = DAbounds.t, percSR.t = percSR.t))
  }
  
  getVE = function(date.t, 
                   date.delta = date.delta, #  as.Date('2021/07/01'), 
                   date.omicron = date.omicron,
                   VE1wt = VE1wt,
                   VE2wt = VE2wt,  # higher b/c we are using mortality data too
                   VE1delta = VE1delta,
                   VE2delta = VE2delta,  # higher b/c we are using mortality data too
                   VE1omicron = VE1omicron,
                   VE2omicron = VE2omicron
  ){
    if(date.t < date.delta){
      return(list(VE1 = VE1wt, VE2 = VE2wt))
    } else if(date.t < date.omicron){
      return(list(VE1 = VE1delta, VE2 = VE2delta))
    } else {
      return(list(VE1 = VE1omicron, VE2 = VE2omicron))
    }
  }
  
  # get parameters needed to specific the probability dist of immune loss for those gained immunity from vaccination
  get_parmVimmLoss = function(date.t, 
                              date.delta = date.delta, #  as.Date('2021/07/01'), 
                              date.omicron = date.omicron,
                              # wt, alpha
                              tm.imm.wt = tm.imm.wt, # during of vaccine-incuded protection against infection
                              # tm.ini.imm.wt = tm.ini.imm.wt, # initial period with near 0 imm loss
                              p.imm.wane.max.wt = p.imm.wane.max.wt, # maximal level of immunity loss (=1 or lower)
                              k.wt = k.wt, # logistic function tuning parameter
                              
                              tm.imm.delta = tm.imm.delta, # during of vaccine-incuded protection against infection
                              # tm.ini.imm.delta = tm.ini.imm.delta, # initial period with near 0 imm loss
                              p.imm.wane.max.delta = p.imm.wane.max.delta, # maximal level of immunity loss (=1 or lower)
                              k.delta = k.delta, # logistic function tuning parameter
                              
                              tm.imm.omicron = tm.imm.omicron, # during of vaccine-incuded protection against infection
                              # tm.ini.imm.omicron = tm.ini.imm.omicron, # initial period with near 0 imm loss
                              p.imm.wane.max.omicron = p.imm.wane.max.omicron, # maximal level of immunity loss (=1 or lower)
                              k.omicron = k.omicron # logistic function tuning parameter
  ){
    if(date.t < date.delta){
      return(list(tm.imm = tm.imm.wt, p.imm.wane.max = p.imm.wane.max.wt, k = k.wt))  # tm.ini.imm = tm.ini.imm.wt, 
    } else if(date.t < date.omicron){
      return(list(tm.imm = tm.imm.delta, p.imm.wane.max = p.imm.wane.max.delta, k = k.delta))  #  tm.ini.imm = tm.ini.imm.delta,
    } else {
      return(list(tm.imm = tm.imm.omicron, p.imm.wane.max = p.imm.wane.max.omicron, k = k.omicron)) # tm.ini.imm = tm.ini.imm.omicron, 
    }
  }
  
  # check immune evasion this wave - the guage level of SR on S
  fn_cum.dS.t = function(date.t, 
                         stat.dS = stat.dS,
                         xpost = xpost,
                         dS.cut.mn = .85, # threshold, if >75% already, no more major SR on this one
                         dS.cut.ens = .6
  ){
    
    allowSRprobe = T; restrictSome = F; idx.no.update = NULL; ens_dImm.max = NULL
    
    parm.bound_vec$date.start = parm.bound_vec$date.start %>% as.Date()
    parm.bound_vec$date.end = parm.bound_vec$date.end %>% as.Date()
    date.t = date.t  %>% as.Date()
    
    tm_variant = parm.bound_vec[type == 'variant characterization']
    # exclude duplicate if the overall omicron period is included
    if(any(duplicated(tm_variant$date.start))){
      tm_variant = tm_variant[variant != 'Omicron']
    }
    iv = which(tm_variant$date.start <= date.t & tm_variant$date.end >= date.t)
    # immunity - imm evasion
    idx0 = which(as.Date(Week.starts) >=  as.Date(tm_variant[iv]$date.start)) %>% head(1)
    idx1 = which(as.Date(Week.starts) <=  as.Date(tm_variant[iv]$date.end)) %>% tail(1)
    idx.full.t = idx0: idx1
    
    if(! is.null(stat.dS)){
      dS.t = stat.dS[week %in% idx.full.t] # , c('week', 'S0.mean', 'dS.mean', 'dS.sd')
      # cum.dS.t = dS.t$dS.mean %>% sum()
      # cum.imm0.t = N - dS.t$S0.mean[1]
      # dImm.t = cum.dS.t / cum.imm0.t
      if(nrow(dS.t) > 0){
        dS.tot.mn = sum(dS.t$dS.mean)
        S0.t = dS.t$S0.mean[1] # susceptiblity when it started the adjustment
        wk.t = dS.t$week[1] - 1  # the week prior to the first adjustment
        
        S0.prior.strt = xpost['S1',,wk.t] 
        S0.prior.strt.mn = S0.prior.strt %>% mean
        
        # no further probing, if the mean already exceed the cutoff
        if(dS.tot.mn / (N - S0.prior.strt.mn) > dS.cut.mn){
          allowSRprobe = F
        } else {
          # look at if there are individual ens needing some restriction
          ens_dS.t = colSums(dS.t[,4+1:num_ens,with=F])
          
          perc.dImm.t = ens_dS.t / (N - S0.prior.strt)
          idx.no.update = which(perc.dImm.t > dS.cut.ens)
          
          # also recall the maximum residual dImm 
          # state0['S1', Sidx.t] = state0['S1', Sidx.t] + (N - state0['S1', Sidx.t]) * runif(length(Sidx.t), s.lwr.t, s.upr.t)
          ens_dImm.max = ((N - S0.prior.strt) * (.85 - perc.dImm.t)) %>% pmax(0) # don't go above 95%
          # idx = which(is.na(ens_dImm.max))
          ens_dImm.max[is.na(ens_dImm.max)] = 0
          
          if(length(idx.no.update)>0) 
            restrictSome = T
        }
        
      }
    }
    
    
    return(list(allowSRprobe = allowSRprobe, restrictSome = restrictSome, idx.no.update = idx.no.update, ens_dImm.max = ens_dImm.max))
  }
  
  fn_get.fcast.mob = function(fcast.wk.starts = fcast.wk.starts, 
                              weeks.fcast = weeks.fcast, da.mob.t = da.mob.t, 
                              num_wk.fcast= num_wk.fcast){
    fcast.mob.t = da.mob.t[as.Date(date) %in% fcast.wk.starts & data.type == 'business']$mob
    if(length(fcast.mob.t) == num_wk.fcast){
      fcast.mob = matrix(da.mob.t[as.Date(date) %in% fcast.wk.starts & data.type == 'business']$mob, nrow = num_wk.fcast, ncol = 1)
    } else if(length(fcast.mob.t) < num_wk.fcast){
      # incomplete data, use the latest for the missing week
      # but if the latest week is the summer, due to hot weather mob might be low
      # use the same calendar week instead
      mob_by_week = da.mob.t %>% dplyr::filter(date > as.Date('2020/4/1')) %>% # exclude pre-pandemic
        dcast(., week ~ year, value.var = 'mob')
      
      tmp.mob = fcast.mob.t %>% tail(2) %>% mean
      fcast.mob.full.t = data.table(week = weeks.fcast, mob = rep(-1, num_wk.fcast))
      fcast.mob.full.t = fcast.mob.full.t %>% left_join(x = ., y = mob_by_week, by = 'week') %>% data.table()
      # there could be no data at all
      if(length(fcast.mob.t) > 0){
        fcast.mob.full.t[1:length(fcast.mob.t)]$mob = fcast.mob.t
      } 
      fcast.mob.full.t[mob < 0]$mob = fcast.mob.full.t[mob < 0, -c(1,2)] %>% apply(., 1, max, na.rm = T) %>% pmax(., tmp.mob, na.rm = T)
      fcast.mob = matrix(fcast.mob.full.t$mob, nrow = num_wk.fcast, ncol = 1)
    }
    
    fcast.mob
    
  }
  
  fn_adj.parm4winter = function(parm.t = state0['alpha',], # parm to adjust
                                parm.name.t = 'alpha',
                                parm.low.cut = .15, parm.winter.amp = .25, parm.phase.shift = 10,
                                vdate.t = vdate.t){
    if(mean(parm.t) < parm.low.cut){
      parm.t = parm.t *  pmax(1, (1 + parm.winter.amp * cos(week(as.Date(vdate.t) + parm.phase.shift)/52 * 4*(base::pi)))) # 1.2
      print(paste('incr', parm.name.t, 'b/c winter'))
    }
    parm.t
  }
  
  # get adjustment for holiday tx
  fn_get.p.beta.summer.holiday = function(vdate.t, # this week
                                          wk.train.end, # week last estimate was made
                                          p.beta.summer.holiday = 1.2, # maximum incr
                                          p.beta.after.holiday = .9 # if the last estimate was during the holiday, then to adjust it down
  ){
    yr.this = format(as.Date(vdate.t),'%Y')
    summer.holidays_wide = seq(as.Date(paste0(yr.this, '/6/24')), as.Date(paste0(yr.this, '/7/21')), by = 'day')
    summer.holidays0a = seq(as.Date(paste0(yr.this, '/6/24')), length.out = 7, by = 'day') # estimates are often a bit delayed, so shift by 1 week
    summer.holidays0b = seq(as.Date(paste0(yr.this, '/7/1')), length.out = 7, by = 'day') # estimates are often a bit delayed, so shift by 1 week
    summer.holidays1 = seq(as.Date(paste0(yr.this, '/7/1')) + 7, length.out = 7, by = 'day') # estimates are often a bit delayed, so shift by 1 week
    summer.holidays2 = seq(as.Date(paste0(yr.this, '/7/8')) + 7, length.out = 7, by = 'day') # estimates are often a bit delayed, so shift by 1 week
    summer.holidays3 = seq(as.Date(paste0(yr.this, '/7/15')) + 7, length.out = 7, by = 'day')
    
    if(!as.Date(wk.train.end) %in% summer.holidays_wide){
      if(as.Date(vdate.t) %in% c(seq(as.Date(paste0(yr.this,'/6/24')), length.out = 7, by = 'day') # ,
                                 # seq(as.Date(paste0(yr.this,'/7/22')), length.out = 3, by = 'day')
      )){
        p.beta.summer.holiday.t = p.beta.summer.holiday - .15 # 1.1
      } else if(as.Date(vdate.t) %in% seq(as.Date(paste0(yr.this,'/7/1')), length.out = 14, by = 'day')){
        p.beta.summer.holiday.t = p.beta.summer.holiday # 1.2
      } else if(as.Date(vdate.t) %in% seq(as.Date(paste0(yr.this,'/7/15')), length.out = 7, by = 'day')){
        p.beta.summer.holiday.t = p.beta.summer.holiday - .05 # 1.15
      } else {
        p.beta.summer.holiday.t = 1
      }
    } else if(as.Date(wk.train.end) %in% summer.holidays0a){
      if(as.Date(vdate.t) %in% seq(as.Date(paste0(yr.this,'/7/1')), length.out = 14, by = 'day')){
        p.beta.summer.holiday.t = p.beta.summer.holiday  # 1.2
      } else if(as.Date(vdate.t) %in% seq(as.Date(paste0(yr.this,'/7/15')), length.out = 7, by = 'day')){
        p.beta.summer.holiday.t = p.beta.summer.holiday - .05 # 1.15
      } else {
        p.beta.summer.holiday.t = 1
      }
    } else if(as.Date(wk.train.end) %in% summer.holidays0b){
      if(as.Date(vdate.t) %in% seq(as.Date(paste0(yr.this,'/7/8')), length.out = 7, by = 'day')){
        p.beta.summer.holiday.t = p.beta.summer.holiday - .05 # 1.2
      } else if(as.Date(vdate.t) %in% seq(as.Date(paste0(yr.this,'/7/15')), length.out = 7, by = 'day')){
        p.beta.summer.holiday.t = p.beta.summer.holiday - .1 # 1.15
      } else {
        p.beta.summer.holiday.t = 1
      }
    } else if(as.Date(wk.train.end) %in% summer.holidays1){
      if(as.Date(vdate.t) %in% summer.holidays2){
        p.beta.summer.holiday.t = 1 # also high
      } else if(as.Date(vdate.t) %in% summer.holidays3){
        p.beta.summer.holiday.t = pmin(1, p.beta.after.holiday + .1)
      } else {
        p.beta.summer.holiday.t = pmin(1, p.beta.after.holiday + .05)
      }
      
    } else if(as.Date(wk.train.end) %in% summer.holidays2){
      if(as.Date(vdate.t) %in% summer.holidays3){
        p.beta.summer.holiday.t = pmin(1, p.beta.after.holiday + .1)
      } else {
        p.beta.summer.holiday.t = pmin(1, p.beta.after.holiday + .05) 
      }
    } else if(as.Date(wk.train.end) %in% summer.holidays3){
      p.beta.summer.holiday.t = pmin(1, p.beta.after.holiday + .1)
    } else {
      p.beta.summer.holiday.t = 1
    }
    
    p.beta.summer.holiday.t
  }
  
  # function to check if there is a rising new variant in the past few weeks?
  fn_checkNewV = function(DAT.VARIANT, # variant data
                          loc.t, # location
                          fcast.start, # week of fcast
                          exclBA.1 = T){
    da.variant = DAT.VARIANT %>% dplyr::filter(location == loc.t & as.Date(week) < (fcast.start - 7*2) & as.Date(week) >= (fcast.start - 7 * 6)) %>% data.table # assume a 2 week delay in variant data reporting
    # include more weeks for recording tm2dominant
    da.variant2 = DAT.VARIANT %>% dplyr::filter(location == loc.t & as.Date(week) < (fcast.start - 7*2) & as.Date(week) >= (fcast.start - 7 * 18)) %>% data.table # assume a 2 week delay in variant data reporting
    
    majorVs = da.variant[,!colnames(da.variant) %in% c('location','week','Omicron','Others'),with=F]
    perc.max = majorVs %>% apply(., 2, max)
    majorVs = majorVs[,perc.max > .1, with =F]
    # exclude BA.1, b/c it has been accounted for
    if(exclBA.1)
      majorVs = majorVs[,!grepl('BA.1', names(majorVs)), with = F]
    # corresponding variant with longer history
    majorVs2 = da.variant2[,names(majorVs), with=F]
    flag.newVrising = 0
    wk2newVreplace = 0 # also project wk the new variant will reach ~100%
    dur.newVreplace = 0 # record weeks it would take to replace from beginning
    v.growth.rate = 0
    newV = NULL
    QC = '' # for quality check
    if(ncol(majorVs) >= 1){
      for(j in 1:ncol(majorVs)){
        x = 1:nrow(majorVs)
        
        # linear model, assume exp growth -> overestimate
        # do not use to project replacement time
        y = log((majorVs[,j,with=F] %>% unlist) + 1e-5)
        v.gr.lm.model = lm(y~x)
        v.gr.t = v.gr.lm.model$coefficients['x'] %>% exp %>% round(2) # this is y1/y0, fold incr by week
        v.gr.model = stats::smooth.spline(x = x, y = majorVs[,j,with=F] %>% unlist)  # , all.knots = F
        
        if(names(majorVs)[j]=='Omicron_nonBA.1o2' & v.gr.t > 1.5) # this combines BA4/BA5, so would be higher than others
          v.gr.t = v.gr.t * ifelse(v.gr.t > 1.7, 3/4, .85) # assume 1/4 is BA.4
        
        if(v.gr.t >1.1 & tail(unlist(majorVs[,j,with=F]),1)<.8){ # 10% increase per week
          flag.newVrising = flag.newVrising + unname(v.gr.t)
          # project using the model to see when the new variant will reach 100%
          v.gr.pred = data.table(x=seq(nrow(majorVs)+1,length.out = 16, by =1), y = NA)
          v.gr.pred$y = predict(v.gr.model, x = v.gr.pred$x)$y # %>% exp
          wk2newVreplace.t = which(v.gr.pred$y > .9) %>% head(1) 
          
          if(length(wk2newVreplace.t) < 1)
            wk2newVreplace.t = nrow(v.gr.pred)
          
          # record weeks it would take to replace from beginning
          v.gr.pred2 = rbind(data.table(x = 1:nrow(majorVs2), y = majorVs2[,j,with=F] %>% unlist),
                             v.gr.pred)
          idx0 = which(v.gr.pred2$y > .05) %>% head(1)
          idx1 = which(v.gr.pred2$y > .75) %>% head(1)
          if(length(idx0) < 1)
            idx0 = 1
          if(length(idx1) < 1)
            idx1 = nrow(v.gr.pred2)
          dur.newVreplace.t = idx1 - idx0 + 1
          if(names(majorVs)[j]=='Omicron_nonBA.1o2' & dur.newVreplace.t < 8) # this combines BA4/BA5, so would be higher than others
            dur.newVreplace.t = dur.newVreplace.t * 1.25 # assume 20% longer
          QC.t = ''
          if(mean(tail(unlist(majorVs[,j,with=F]),3))<.1)
            QC.t = 'low perc'
          
          if(wk2newVreplace==0){
            wk2newVreplace = wk2newVreplace.t # first new v
            v.growth.rate = v.gr.t
            newV = names(majorVs)[j]
            dur.newVreplace = dur.newVreplace.t
            QC = QC.t
            
          } else {
            # if multiple variants
            wk2newVreplace = append(wk2newVreplace, wk2newVreplace.t)
            v.growth.rate = append(v.growth.rate, v.gr.t)
            newV = append(newV, names(majorVs)[j])
            dur.newVreplace = append(dur.newVreplace, dur.newVreplace.t)
            QC = append(QC, QC.t)
          }
        } 
        
        # allow multiple new variants co-emerging
        # and based on growth rate
      }
    }
    return(list(flag.newVrising = flag.newVrising, 
                newv.pred = data.table(newV = newV,
                                       v.growth.rate = v.growth.rate,
                                       wk2newVreplace = wk2newVreplace,
                                       dur.newVreplace = dur.newVreplace,
                                       QC = QC)))
  }
}


EAKF_rproj = function(epi.model=epi.model, num_ens=num_ens,inflat=1.03, fcast.deflat = .97,
                obs_i=obs_i, obs_vars_i=obs_vars_i, # case
                obs_d=obs_d, obs_vars_d=obs_vars_d,
                weeks=weeks,Week.starts=Week.starts,
                parm.bounds=parm.bounds, DAbounds=DAbounds, SRbounds=SRbounds, 
                parm.names = rownames(parm.bounds), rel.mob = rel.mob,
                state0=state0, state.names=rownames(state0),
                severity = severity,
                tm.ini=1, tmstep=7,
                newI.previous = NULL,
                parm.bound_vec = parm.bound_vec,
                # for the projection
                weeks.fcast, # week of the year to get seasonality 
                fcast.wk.starts,
                incl.sce.newV = T, # whether to include scenarios of new variant during projection
                save.cumIetc = F # whether to save addition model outputs on ens cumI, S, etc.
){
  # parm.bounds: prior bounds for the parms
  # x.prior: priors for all state variables, and model parameters
  # obs_i: observations
  # idx.obs: idx for obseverations to filter through
  
  # save initial condition passed in
  state00 = state0
  DAbounds00=DAbounds;
  SRbounds00=SRbounds;
  newI.previous00 = newI.previous;
  
  DAbounds.t = DAbounds;
  SRbounds.t = SRbounds;
  
  
  if(!exists('redn.priority')) redn.priority = 1
  if(!exists('exclExtreme')) exclExtreme = F # wether to exclude extreme values before the projection
  if(!exists('perc2excl')) perc2excl = .1 # default: exclude 10%, 5% lowest and 5% highest
  if(!exists('v.deflat')) v.deflat = rownames(state0) # set which variables to apply deflation
  
  if(is.null(dim(obs_i))){ # only 1 observation (no subgroup)
    num_times_i=length(obs_i);  
  } else {
    num_times_i = nrow(obs_i)
  }
  if(is.null(dim(obs_d))){ # only 1 observation (no subgroup)
    num_times_d=length(obs_d);  
  } else {
    num_times_d = nrow(obs_d)
  }
  num_var = nrow(state0)
  
  cumlike=NULL; # cummulative likelihood
  
  xprior=array(0,c(num_var,num_ens,num_times_i+1));
  xpost=array(0,c(num_var,num_ens,num_times_i));
  dimnames(xpost)[1] = list(state.names)
  
  # to estimate daily cases
  xprior.daily = matrix(0,(num_times_i+1)*tmstep, num_ens);
  xpost.daily = matrix(0,(num_times_i)*tmstep, num_ens);
  
  Eprior = Iprior = Sprior = numeric(num_times_i+1)
  Epost = Ipost = Spost = Itotpost = numeric(num_times_i)
  
  stat.dS = NULL  # to record the changes in S - immune evasion
  s.lwr.t = .3; s.upr.t = .8 # level of imm evasion, for SR on S
  # cntSR.S = 0  # need to be reset for subsequent imm evasive wave
  cntSR.Scut = 4 # no more than 4
  cntSR.S = cntSR.Scut  # set the initial value to the max so that the weeks can be specified seperately
  dS.cut.mn = .85; # threshold, if >85% already, no more major SR on this one
  dS.cut.ens = .85; 
  
  {
    dist_tm.to.detect = NULL
    for(ii in 1:num_ens){
      tmp = generation.time(dist_tm.to.detect.name,c(state0['Td.mean',ii],state0['Td.sd',ii]),truncate = tm.to.detect.max)
      dist_tm.to.detect=cbind(dist_tm.to.detect,tmp$GT[-1]); 
    }
    
    
    dist_tm.to.death = NULL  # time to death
    for(ii in 1:num_ens){
      # tmp = generation.time(dist_tm.to.death.name,c(state0['Td.mean',ii]+diff.dd,state0['Td.sd',ii]+diff.sd2),truncate = tm.to.death.max)
      # do not link it to Td
      tmp = generation.time(dist_tm.to.death.name,c(tm.to.outcome.mn['tm.to.death',ii], tm.to.outcome.sd['tm.to.death',ii]),truncate = tm.to.death.max)
      dist_tm.to.death=cbind(dist_tm.to.death,tmp$GT[-1]); 
    }
    if(tm.to.deathFromDiag){ # if the distribution of time to death is from diagnosis, not infectious
      # add time from infectious to diagnosis
      dist_tm.to.death = rbind(matrix(0,round(tm.to.diag,0),num_ens),dist_tm.to.death)
    }
  }
  
  # integrate forward 1 step to get the prior
  tm_strt = tm.ini+1; tm_end = tm.ini + tmstep;
  cur.wk = weeks[1];
  seed = seed
  vdate.t = Week.starts[1] %>% as.Date()
  
  severity.t = severity
  
  tmp = getVE(date.t = vdate.t,
              date.delta = date.delta, #  as.Date('2021/07/01'), 
               date.omicron = date.omicron,
               VE1wt = VE1wt,
               VE2wt = VE2wt,  # higher b/c we are using mortality data too
               VE1delta = VE1delta,
               VE2delta = VE2delta,  # higher b/c we are using mortality data too
               VE1omicron = VE1omicron,
               VE2omicron = VE2omicron
              )
  VE1 = tmp$VE1; VE2 = tmp$VE2
  
  newVacc.t = NULL # for tracking immune loss for those vaccinated
  cumVimmLoss.t = rep(0, num_ens) # to exclude changes in S due to this immune loss
  Rimm.t  = NULL # for tracking those recovered and remain immune
  ts_cumVimmLoss = NULL # tracking immune loss of vaccinees
  ts_cumIimmLoss = NULL # tracking immune loss of recoverees
  ts_R = NULL # tracking recoverees
  # get prior for week 1
  {
    if(!is.null(newI.previous) & vdate.t >= vax.start){
      tm.t = nrow(newI.previous)
      # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
      # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
      # 2/5/21 set t1 to 1 yr given the slow rollout
      t1 = 365
      cumI.t = colSums(newI.previous[pmax(1,tm.t-t1) : (tm.t),]) #  %>% apply(1, median) # excl last two weeks and get the median
      # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
      # higher infection rate for the priority group
      # tm.t = pmax(1, tm.t - t1 + 1) # re-aline timing with start of vac
      tm.imm = 365*2.5 # assume 3 yr immunity
      p.imm.wane.max = .8; k = .015  # 1/tm.imm  
      # since only the last year is included, should be:
      p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(pmin(t1, tm.t) + 60 - tm.imm/2))) # not all infected from day 1
      # earlier infections are likely to be in the high priority groups 
      p.imm = 1 *  p.imm.wane * redn.priority # assume 100% prior infection provide immunity, but wane over time
      # and multiple by % excluded if there is prior testing before vax: p.prior.test
      percSmax.t = 1 - cumI.t / N * p.imm
      # also compare to current susceptibility, in case of immune evasion that increases susceptiblity
      percSmax.t = pmax(percSmax.t, state0[paste0('S',1:num_gr),] / N)
      # no lower than 50%, in case of outliers
      # percSmax.t = pmax(percSmax.t, .5)
      # print(c('cohort %S:',round(summary(mean(percSmax.t)),2)), quote = F)
    } else {
      percSmax.t = 0
      # print('no vax yet')
    }
    
    beta_tt = state0['beta',] * fn_getRelMob(rel.mob[1,], state0['p.mob',])
    if(seasonality) {
      # beta_tt = state0['beta',] * seasonal.cycle[seasonal.cycle$week==cur.wk,]$relR0 
      beta_tt = beta_tt * relR0[cur.wk,] 
    } 
    
    if(epi.model == 'SEIRS'){
      
      simEpi=SEIRS(tm_strt, tm_end, tm_step=1, # 1 day time-step
                   tmstep = tmstep,
                   state0 = state0,
                   S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                   I0=state0[paste0('I',1:num_gr),], 
                   beta=beta_tt, 
                   Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                   seed=seed, stoch=stoch, 
                   severity = severity.t,
                   newI.previous = newI.previous,
                   dist_tm.to.detect = dist_tm.to.detect,
                   dist_tm.to.death = dist_tm.to.death) # include the delay reporting etc.
      
    } else if(epi.model == 'SEIRSV') {
      
      daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
      
      if(nrow(daVacc.t)<1){  # no data yet
        V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
      } else { # yes data
        
        daVacc.t$date = daVacc.t$date %>% as.Date
        
        # make sure it includes a full week
        dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
        daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
        daVacc.t[is.na(daVacc.t)] = 0
        V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
        V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
        
        # print('start vacc!')
        
      }
      simEpi=SEIRSV(tm_strt, tm_end, tm_step=1, # 1 day time-step
                    tmstep = tmstep,
                    state0 = state0,
                    S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                    I0=state0[paste0('I',1:num_gr),], 
                    beta=beta_tt, 
                    Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                    seed=seed, stoch=stoch, 
                    severity = severity.t,
                    newI.previous = newI.previous,
                    dist_tm.to.detect = dist_tm.to.detect,
                    dist_tm.to.death = dist_tm.to.death,
                    percSmax.t = percSmax.t,
                    V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                    # these are total number of vaccinees with unknown immunity
                    # but pre-ajust for time lag from vaccination to immune protection
                    VE1 = VE1, VE2=VE2 # Vaccine efficacy, need further adjustment by prior immunity 
      ) # include the delay reporting etc.
    } else if (epi.model == 'SEIRSVimmLoss'){
      
      daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
      
      if(nrow(daVacc.t)<1){  # no data yet
        V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
      } else { # yes data
        
        daVacc.t$date = daVacc.t$date %>% as.Date
        
        # make sure it includes a full week
        dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
        daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
        daVacc.t[is.na(daVacc.t)] = 0
        V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
        V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
        
        # print('start vacc!')
        
      }
      
      # get vaccine-induced immunity related parameter
      parmVimmLoss.t = get_parmVimmLoss(date.t = vdate.t, 
                                        date.delta = date.delta, #  as.Date('2021/07/01'), 
                                        date.omicron = date.omicron,
                                        # wt, alpha
                                        tm.imm.wt = tm.imm.wt, # during of vaccine-incuded protection against infection
                                        p.imm.wane.max.wt = p.imm.wane.max.wt, # maximal level of immunity loss (=1 or lower)
                                        k.wt = k.wt, # logistic function tuning parameter
                                        
                                        tm.imm.delta = tm.imm.delta, # during of vaccine-incuded protection against infection
                                        p.imm.wane.max.delta = p.imm.wane.max.delta, # maximal level of immunity loss (=1 or lower)
                                        k.delta = k.delta, # logistic function tuning parameter
                                        
                                        tm.imm.omicron = tm.imm.omicron, # during of vaccine-incuded protection against infection
                                        p.imm.wane.max.omicron = p.imm.wane.max.omicron, # maximal level of immunity loss (=1 or lower)
                                        k.omicron = k.omicron # logistic function tuning parameter
      )
      
      # number of people infected, recovered, and remain immune
      if(is.null(Rimm.t )){
        Rimm.t  = N - state0[paste0('S',1:num_gr),] - state0[paste0('E',1:num_gr),] - state0[paste0('I',1:num_gr),]
        Rimm.t[Rimm.t < 0] = 0
        # Rimm.t  = rep(0, num_ens)
      }
      
      
      simEpi=SEIRSVimmLoss(tm_strt, tm_end, tm_step=1, # 1 day time-step
                    tmstep = tmstep,
                    state0 = state0,
                    S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                    I0=state0[paste0('I',1:num_gr),], 
                    beta=beta_tt, 
                    Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                    seed=seed, stoch=stoch, 
                    severity = severity.t,
                    newI.previous = newI.previous,
                    dist_tm.to.detect = dist_tm.to.detect,
                    dist_tm.to.death = dist_tm.to.death,
                    percSmax.t = percSmax.t,
                    V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                    # these are total number of vaccinees with unknown immunity
                    # but pre-ajust for time lag from vaccination to immune protection
                    VE1 = VE1, VE2=VE2, # Vaccine efficacy, need further adjustment by prior immunity 
                    newVacc = newVacc.t, # prior vaccinated & protected from infection
                    parmVimmLoss = parmVimmLoss.t, # parameters related to the immunity period of vaccine, against infection
                    R0 = Rimm.t   # to check those recovered from infection and remaining immune to infection
      ) # include the delay reporting etc.
      
      # update hist newVacc.t
      newVacc.t = simEpi$newVacc
      cumVimmLoss.t = simEpi$cumVimmLoss
      cumIimmLoss.t = simEpi$cumIimmLoss #  %>% c()
      Rimm.t  = simEpi$R
      simEpi$newVacc = NULL
      simEpi$cumVimmLoss = NULL
      simEpi$R = NULL
      simEpi$cumIimmLoss = NULL
      
      ts_cumVimmLoss = rbind(ts_cumVimmLoss, cumVimmLoss.t)
      ts_cumIimmLoss = rbind(ts_cumIimmLoss, cumIimmLoss.t)
      ts_R = rbind(ts_R, Rimm.t)
    }
    
    # re-assemble to the same order as the prior: state0
    n.end = tm_end - tm_strt + 2
    state.new = NULL
    for(i in 1:(length(simEpi)-1)){
      tmp = simEpi[[i]][n.end,,drop=F]; 
      rownames(tmp)=gsub('cumI','newI',paste0(names(simEpi)[i],1:num_gr))
      state.new = rbind(state.new,tmp)
    }
    
    state.new = rbind(state.new, state0[parm.names,])
    state.new = state.new[rownames(state0),] # make sure the order is the same
    
    
    xprior[,,1]= state.new
    xprior.daily[1:tmstep,] = simEpi$daily.newItot
    # save prior for E and I
    Eprior[1] = state.new['E1',] %>% mean
    Iprior[1] = state.new['I1',] %>% mean
  }
  
  dx.t.newitot = dx.t.newiobs = dx.t.e = dx.t.i = numeric(num_times_i)
  #### Begin looping through observations
  for (tt in 1:num_times_i){ # num_times_i
    
    xmn=rowMeans(xprior[,,tt]); 
    xnew=inflat*(xprior[,,tt]-xmn%*%matrix(1,1,num_ens))+xmn%*%matrix(1,1,num_ens);
    row.names(xnew)=state.names
    for(H in idx.obs_i){  # LOOP THROUGH OBSERVATIONS
      ####  Get the variance of the ensemble
      ####  Define the mapping operator H
      ####  HH is the number of variables and parameters in state space
      ####  H is only the number of columns being observed
      io = H - min(idx.obs_i) + 1
      obs_var = obs_vars_i[tt,io];
      prior_var = pmax(var(xnew[H,]),1);
      # post_var = 1/(1/prior_var + 1/obs_var);
      post_var = prior_var*obs_var/(prior_var+obs_var);
      
      prior_mean = mean(xnew[H,]);
      post_mean = post_var*(prior_mean/prior_var + obs_i[tt,io]/obs_var);
      
      #### Compute alpha and adjust distribution to conform to posterior moments
      alp = sqrt(obs_var/(obs_var+prior_var));
      
      dy = post_mean + alp*(xnew[H,]-prior_mean)-xnew[H,];
      ###  Getting the covariance of the prior state space and
      ###  observations  (which could be part of state space, e.g. infections)
      rr=NULL;
      for (j in 1:num_var){
        C=cov(xnew[j,],xnew[H,])/prior_var;  # covariance/varance of x.obs
        
        if(donotUpdateS & vdate.t %in% tm_rednUpdateS){  # restrict updating of S - before the surge of immune evasive variant
          if(grepl('S',state.names[j])) C = C / 10
        }
        if(rednUpdateEI  & vdate.t %in% tm_rednUpdateEI){  # key period of time to restrict updating of E and I
          if((grepl('E',state.names[j]) | grepl('I',state.names[j])) &
             !grepl('new',state.names[j])){
            C = C / 10 # reduce the level of update by a factor of 10
          }
        }
        
        rr=append(rr,C);
      }
      dx=rr%*%t(dy);
      
      # record the level of adjustment
      dx.t = rowMeans(dx); names(dx.t) = state.names
      dx.t.newitot[tt] = dx.t['newItot1']
      dx.t.newiobs[tt] = dx.t['newIobs1']
      dx.t.e[tt] = dx.t['E1']
      dx.t.i[tt] = dx.t['I1']
      
      ###  Get the new ensemble and save prior and posterior
      xnew = xnew + dx;  # ADJUSTED NUM_OBS TIMES, ONCE BY EACH OBSERVATION
    }
    
    # filter using the mortality data
    if(tt %in% 1:num_times_d){
      for(H in idx.obs_d){  # LOOP THROUGH OBSERVATIONS
        ####  Get the variance of the ensemble
        ####  Define the mapping operator H
        ####  HH is the number of variables and parameters in state space
        ####  H is only the number of columns being observed
        io = H - min(idx.obs_d) + 1
        obs_var = obs_vars_d[tt,io];
        prior_var = pmax(var(xnew[H,]),1);
        # post_var = 1/(1/prior_var + 1/obs_var);
        post_var = prior_var*obs_var/(prior_var+obs_var);
        
        prior_mean = mean(xnew[H,]);
        post_mean = post_var*(prior_mean/prior_var + obs_d[tt,io]/obs_var);
        
        #### Compute alpha and adjust distribution to conform to posterior moments
        alp = sqrt(obs_var/(obs_var+prior_var));
        
        dy = post_mean + alp*(xnew[H,]-prior_mean)-xnew[H,];
        ###  Getting the covariance of the prior state space and
        ###  observations  (which could be part of state space, e.g. infections)
        rr=NULL;
        for (j in 1:num_var){
          C=cov(xnew[j,],xnew[H,])/prior_var;  # covariance/varance of x.obs
          
          if(donotUpdateS & vdate.t %in% tm_rednUpdateS){
            if(grepl('S',state.names[j])) C = C / 10
          }
          if(rednUpdateEI & vdate.t %in% tm_rednUpdateEI){ # & tt> end1stWave
            if((grepl('E',state.names[j]) | grepl('I',state.names[j])) &
               !grepl('new',state.names[j])){
              C = C / 10 # reduce the level of update by a factor of 10
            }
          }
          
          rr=append(rr,C);
        }
        dx=rr%*%t(dy);
        
        ###  Get the new ensemble and save prior and posterior
        xnew = xnew + dx;  # ADJUSTED NUM_OBS TIMES, ONCE BY EACH OBSERVATION
      }
    }
    
    #  Corrections to DA produced aphysicalities
    xnew=Fn_checkDA(xnew, bound.low = DAbounds.t[,1], bound.up =  DAbounds.t[,2]);
    row.names(xnew)=state.names
    xpost[,,tt]=xnew;
    state0 = xnew
    
    # save the state if this will be the start of the forecast
    if(as.Date(Week.starts[tt]) == (as.Date(fcast.start) - 7)){
      state0.fcast = state0
      # also record the current number ppl recovered and remain immune
      Rimm.t.fcast = Rimm.t
      newVacc.t.fcast = newVacc.t # tracking immune loss in vaccinees
    }
    
    # save post for E and I
    Epost[tt] = state0['E1',] %>% mean
    Ipost[tt] = state0['I1',] %>% mean
    Spost[tt] = state0['S1',] %>% mean
    Itotpost[tt] = state0['newItot1',] %>% mean
    
    # update daily estimates as well
    {
      
      # avoid dividing by 0 here!
      if(tmstep>1){
        f.adj = state0[idx.newItot,]/colSums(xprior.daily[1:tmstep+(tt-1)*tmstep,])
      } else {
        f.adj = state0[idx.newItot,]/(xprior.daily[1:tmstep+(tt-1)*tmstep,])
      }
      
      
      {
        i0=which(is.na(f.adj) | is.infinite(f.adj)); 
        xx=1:num_ens; 
        i.non0=xx[!(xx %in% i0)]
        if(length(i.non0)>0){
          f.adj[i0] = sample(f.adj[i.non0],size=length(i0),replace = T)
        } else {
          f.adj[i0] = 1 # all are 0
        }
        
      }
      f.adj_all = matrix(f.adj, tmstep, num_ens, byrow = T)
      # if obs=0 for the whole week, set all days to 0 - that is taken care of by default
      # so all the problem comes from obs!=0, but the prior says all days have 0 cases
      # in that case, distribute the cases evenly
      xpost.daily[1:tmstep+(tt-1)*tmstep,]=xprior.daily[(1:tmstep)+(tt-1)*tmstep,]*f.adj_all
      
    } # update daily estimates
    
    
    # compute the posterior of increase in S - immune evasion
    if(!(vdate.t %in% tm_rednUpdateS) & # this week S can be adjusted
       tt > 1 # for the omicron wave, start probing S earliest from 2nd week
      ){ 
      
      if(Spost[tt] - Spost[tt-1] + Itotpost[tt] - mean(cumVimmLoss.t) - mean(cumIimmLoss.t) > .5/100 * N  # greater than .5% of the population
         ){ # there is a large adjustment to S - suggesting immune evasion
        dS.t = xpost['S1',,tt] - xpost['S1',, tt-1] + xpost['newItot1',,tt] - cumVimmLoss.t - cumIimmLoss.t
        # not good b/c the eakf uses the mean for indiv variables
        # dS.t = mean(xpost['S1',,tt]) - mean(xpost['S1',, tt-1]) + mean(xpost['newItot1',,tt])
        
        xx1 = xpost['S1',,tt]; xx1mn = xx1 %>% mean; 
        xx2 = xpost['S1',, tt-1]; xx2mn = xx2 %>% mean; 
        xx3=xpost['newItot1',,tt]; xx3mn = xx3 %>% mean; 
        xx4 = cumVimmLoss.t; xx4mn = xx4 %>% mean; 
        # var1 = var(xx1); var2 = var(xx2);  var3=var(xx3)
        m.cov = cov(x=cbind(xx1, xx2, xx3, xx4))
        
        dS.sd = deltamethod(~ x1 - x2 + x3 - x4, mean = c(xx1mn, xx2mn, xx3mn, xx4mn), 
                            cov = m.cov
        ) %>% sqrt
        
        dS.mn = dS.t  %>% mean
        
        stat.dS = rbind(stat.dS, 
                        data.table(week = tt, S0.mean = Spost[tt-1], # susceptiblity prior to the change
                                   dS.mean = dS.mn,
                                   dS.sd = dS.sd,
                                   dS.t %>% t
                        ))
        
        # also update the Rimm.t box
        Rimm.t = Rimm.t - dS.mn #  - dS.t  Use the instead to reduce variance
        if(as.Date(Week.starts[tt]) == (as.Date(fcast.start) - 7)){
          # also record the current number ppl recovered and remain immune
          Rimm.t.fcast = Rimm.t
        }
      }
      
    }
    
    #  Integrate forward one time step
    tcurrent = tm.ini+tmstep*tt;
    tm_strt = tcurrent+1; tm_end = tcurrent + tmstep
    cur.wk = weeks[tt+1];
    if(is.na(cur.wk))
      cur.wk = weeks[tt];
    # print(paste('Week:',cur.wk),quote=F)
    vdate.t = Week.starts[pmin(tt+1,length(Week.starts))]
    
    newI.previous = xpost.daily[1:(tt*tmstep),] # these are at weekly level!
    
    tmp = getVE(date.t = vdate.t,
                date.delta = date.delta, #  as.Date('2021/07/01'), 
                 date.omicron = date.omicron,
                 VE1wt = VE1wt,
                 VE2wt = VE2wt,  # higher b/c we are using mortality data too
                 VE1delta = VE1delta,
                 VE2delta = VE2delta,  # higher b/c we are using mortality data too
                 VE1omicron = VE1omicron,
                 VE2omicron = VE2omicron
                )
    VE1 = tmp$VE1; VE2 = tmp$VE2
    
    # check the level of adjust for newItot to gauge how the filter has been working
    mn.pr.p = xprior[idx.newItot,,tt-1] %>% mean
    mn.po.p = xpost[idx.newItot,,tt-1] %>% mean
    mn.pr = xprior[idx.newItot,,tt] %>% mean
    mn.po = xpost[idx.newItot,,tt] %>% mean
    e.po = xpost[idx.e,,tt] %>% mean
    i.po = xpost[idx.i,,tt] %>% mean
    e.pr = xprior[idx.e,,tt] %>% mean
    i.pr = xprior[idx.i,,tt] %>% mean
    
    # SR - 
    if(doSR){
      
      # get SR/DA bounds
      tmp = getSRDAbounds(date.t = vdate.t, # date for the current week
                          SRbounds0 = SRbounds00, # original bounds
                          DAbounds0 = DAbounds00, # original bounds
                          parm.bound_vec = parm.bound_vec   # changes by period
      )
      SRbounds.t = tmp$SRbounds.t
      DAbounds.t = tmp$DAbounds.t
      percSR.t = tmp$percSR.t  # number of ensemble members to SR
      
      num_SR.t = round(num_ens*(SR.perc.local + SR.perc.full), 0)
      num_SR.local = round(num_ens*SR.perc.local, 0)
      num_SR.full = round(num_ens*SR.perc.full, 0)
      SR.idx.t = sample(1:num_ens, num_SR.t, replace = F)
      SR.idx.local.t = head(SR.idx.t, num_SR.local)
      SR.idx.full.t = tail(SR.idx.t, num_SR.full)
      
      sr.local.mean = rowMeans(state0[SR.var.local,]) # local SR bounds
      sr.local.lwr = apply(state0[SR.var.local,], 1, quantile, .1) # .2
      sr.local.upr = apply(state0[SR.var.local,], 1, quantile, .9) # .7
      sr.local.lwr = pmin(sr.local.mean * .7, sr.local.lwr) # .75
      sr.local.upr = pmax(sr.local.mean * 1.3, sr.local.upr) # 1.25
      sr.local.lwr = pmax(sr.local.lwr, DAbounds.t[SR.var.local,1]) # .75
      sr.local.upr = pmin(sr.local.upr, DAbounds.t[SR.var.local,2]) # 1.25
      
      state0[SR.var.local,SR.idx.local.t] = t(lhs(num_SR.local, rect =  cbind(sr.local.lwr, sr.local.upr)))
      
      state0[SR.var.full,SR.idx.full.t] = t(lhs(num_SR.full, rect = SRbounds.t[SR.var.full,]))
      
      # any parm need extra SR?
      # idxExtra = which(tolower(percSR.t$percSR) == 'extra')
      idxExtra = which(grepl('extra', tolower(percSR.t$percSR))==T)
      if(length(idxExtra)>0){
        for(idx.t in idxExtra){
          if(percSR.t[idx.t]$percSR == 'EXTRA'){
            num_SR.t = round(num_ens * SR.perc.EXTRA, 0)
          } else if (percSR.t[idx.t]$percSR == 'Extra') {
            num_SR.t = round(num_ens * SR.perc.Extra, 0)
          } else if (percSR.t[idx.t]$percSR == 'Extra1') {
            num_SR.t = round(num_ens * SR.perc.Extra1, 0)
          } else if (percSR.t[idx.t]$percSR == 'Extra2') {
            num_SR.t = round(num_ens * SR.perc.Extra2, 0)
          } else if (percSR.t[idx.t]$percSR == 'Extra3') {
            num_SR.t = round(num_ens * SR.perc.Extra3, 0)
          } else {
            num_SR.t = round(num_ens * SR.perc.extra, 0)
          }
          
          parm.t = percSR.t[idx.t]$parm
          SR.idx.t = sample(1:num_ens, num_SR.t, replace = F)
          state0[parm.t, SR.idx.t] = runif(num_SR.t, min = SRbounds.t[parm.t,1], max = SRbounds.t[parm.t, 2])
        }
      }
      
      # when the new variant is immune evasive, beta tend to overshot then decrease, as it takes several weeks to get S up
      # very small probing on S, when the variant is immune evasive?
      {
        # check if it's a new imm evasive wave so that cntSR.S be reset to 0
        if(vdate.t %in% tm_reset_cntSR.S){
          # print(paste(vdate.t, 'new wave, reset cntSR.S'), quote = F)
          cntSR.S = 0
        }
        if(!(vdate.t %in% tm_rednUpdateS) & # this week S can be adjusted
           tt > 1 # for omicron wave, start earliest from week 2
        ){ 
          
          if(Spost[tt] - Spost[tt-1] + Itotpost[tt] - mean(cumVimmLoss.t) - mean(cumIimmLoss.t) > .5/100 * N &  # greater than .5% of the population
             cntSR.S < cntSR.Scut   # 1/25/22 use lower threshold so the filter can act promptly (raise the bar for % below to control the level - balancing the two)
          ){ 
            # print(paste0(tt, ': SR on S'))
            cntSR.S = cntSR.S + 1
            # allow SR.perc.imm.t to change with the filter adjustment
            # SR.perc.imm.t = pmax(SR.perc.imm, pmin(SR.perc.imm * 4, 
            #                      (Spost[tt] - Spost[tt-1] + Itotpost[tt]) / (4/100 * N) * SR.perc.imm * .8^cntSR.S
            # ))
            
            # try a lower bound, as for some places, the above is not sensitive to changes enough
            # further subtracting mean(cumVimmLoss.t) and mean(cumIimmLoss.t) raise the bar, so should use lower denominator instead
            SR.perc.imm.t = pmax(SR.perc.imm, pmin(SR.perc.imm * 4, 
                                                   (Spost[tt] - Spost[tt-1] + Itotpost[tt] - mean(cumVimmLoss.t) - mean(cumIimmLoss.t)) / (1/100 * N) * SR.perc.imm * .8^cntSR.S  # (3/100 * N)
            ))
            
            # print(paste(tt, '% SRimm', SR.perc.imm.t), quote = F)
          
            num_SR.t = round(num_ens * SR.perc.imm.t, 0)
            Sidx = sample(1:num_ens, num_SR.t, replace = F)
            
            # 
            
            # check current cumulative level of dImm first - if this already very high, do not allow further large probing
            tmp = fn_cum.dS.t(date.t = vdate.t,
                              stat.dS = stat.dS,
                              xpost = xpost,
                                   dS.cut.mn = dS.cut.mn, # threshold, if >75% already, no more major SR on this one
                                   dS.cut.ens = dS.cut.ens)
            
            if(tmp$allowSRprobe){
              if(!tmp$restrictSome){
                state0['S1', Sidx] = state0['S1', Sidx] + (N - state0['S1', Sidx]) * runif(num_SR.t, s.lwr.t, s.upr.t)
              } else {
                Sidx.t = Sidx[! Sidx %in% tmp$idx.no.update]
                # state0['S1', Sidx.t] = state0['S1', Sidx.t] + (N - state0['S1', Sidx.t]) * runif(length(Sidx.t), s.lwr.t, s.upr.t)
                state0['S1', Sidx.t] = state0['S1', Sidx.t] + pmin((N - state0['S1', Sidx.t]) * runif(length(Sidx.t), s.lwr.t, s.upr.t),
                                                                   tmp$ens_dImm.max[Sidx.t])
                
                # print(paste('restrict SR on ', num_SR.t - length(Sidx.t), 'ens'))
              }
            }
            
            
          }
        }
        
        
      }
    }
    
    
    # get prior
    {
      if(!is.null(newI.previous) & vdate.t >= vax.start){
        tm.t = nrow(newI.previous)
        # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
        # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
        # 2/5/21 set t1 to 1 yr given the slow rollout
        t1 = 365
        cumI.t = colSums(newI.previous[pmax(1,tm.t-t1) : (tm.t),]) #  %>% apply(1, median) # excl last two weeks and get the median
        # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
        # higher infection rate for the priority group
        # tm.t = pmax(1, tm.t - t1 + 1) # re-aline timing with start of vac
        tm.imm = 365* 2.5 # assume 3 yr immunity
        p.imm.wane.max = .8; k = .015  # 1/tm.imm  
        # since only the last year is included, should be:
        p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(pmin(t1, tm.t) + 60 - tm.imm/2))) # not all infected from day 1
        # earlier infections are likely to be in the high priority groups 
        p.imm = 1 *  p.imm.wane * redn.priority # assume 100% prior infection provide immunity, but wane over time
        # and multiple by % excluded if there is prior testing before vax: p.prior.test
        percSmax.t = 1 - cumI.t / N * p.imm
        # also compare to current susceptibility, in case of immune evasion that increases susceptiblity
        percSmax.t = pmax(percSmax.t, state0[paste0('S',1:num_gr),] / N)
        # no lower than 50%, in case of outliers
        # percSmax.t = pmax(percSmax.t, .5)
        # print(c('cohort %S:',round(summary(mean(percSmax.t)),2)), quote = F)
      } else {
        percSmax.t = 0
        # print('no vax yet')
      }
      
      beta_tt = state0['beta',] * fn_getRelMob(rel.mob[pmin(tt+1, nrow(rel.mob)),],state0['p.mob',]) 
      if(seasonality) {
        beta_tt = beta_tt * relR0[cur.wk,] 
      } 
      
      severity.t = severity
      severity.t['death',] = state0['ifr',]
      
      # [2/23/21] check this if problem!
      dist_tm.to.detect = NULL
      for(ii in 1:num_ens){
        tmp = generation.time(dist_tm.to.detect.name,c(state0['Td.mean',ii],state0['Td.sd',ii]),truncate = tm.to.detect.max)
        dist_tm.to.detect=cbind(dist_tm.to.detect,tmp$GT[-1]); 
      }
      
      # [1/20/22] also update time from infection to death, b/c longer time-lag for omicron
      if(vdate.t >= (date.omicron+7) & vdate.t < (date.omicron+14)){
        print('update tm to death for omicron')
        
        tm.to.death.max = tm.to.death.max_omicron
        if(tm.to.deathFromDiag){
          tm.from.inf.to.death.max = tm.to.diag + tm.to.death.max
        } else {
          tm.from.inf.to.death.max = tm.to.death.max
        }
        
        dist_tm.to.death = NULL  # time to death
        for(ii in 1:num_ens){
          tmp = generation.time(dist_tm.to.death.name,c(tm.to.outcome.mn['tm.to.death',ii] + tm2death_adj_omicron, 
                                                        tm.to.outcome.sd['tm.to.death',ii] # not good to increase sd
          ),
          truncate = tm.to.death.max  # updated above
          )
          dist_tm.to.death=cbind(dist_tm.to.death,tmp$GT[-1]); 
        }
        if(tm.to.deathFromDiag){ # if the distribution of time to death is from diagnosis, not infectious
          # add time from infectious to diagnosis
          dist_tm.to.death = rbind(matrix(0,round(tm.to.diag,0),num_ens),dist_tm.to.death)
        }
      }
      
      if(epi.model == 'SEIRS'){
        
        simEpi=SEIRS(tm_strt, tm_end, tm_step=1, # 1 day time-step
                     tmstep = tmstep,
                     state0 = state0,
                     S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                     I0=state0[paste0('I',1:num_gr),], 
                     beta=beta_tt, 
                     Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                     seed=seed, stoch=stoch,
                     severity = severity.t,
                     newI.previous = newI.previous,
                     dist_tm.to.detect = dist_tm.to.detect,
                     dist_tm.to.death = dist_tm.to.death)
        
      } else if(epi.model == 'SEIRSV'){
        daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
        
        if(nrow(daVacc.t)<1){  # no data yet
          V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
        } else { # yes data
          
          daVacc.t$date = daVacc.t$date %>% as.Date
          
          # make sure it includes a full week
          dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
          daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
          daVacc.t[is.na(daVacc.t)] = 0
          V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
          V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
          
          # print('start vacc!')
          
        }
        simEpi=SEIRSV(tm_strt, tm_end, tm_step=1, # 1 day time-step
                      tmstep = tmstep,
                      state0 = state0,
                      S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                      I0=state0[paste0('I',1:num_gr),], 
                      beta=beta_tt, 
                      Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                      seed=seed, stoch=stoch, 
                      severity = severity.t,
                      newI.previous = newI.previous,
                      dist_tm.to.detect = dist_tm.to.detect,
                      dist_tm.to.death = dist_tm.to.death,
                      percSmax.t = percSmax.t,
                      V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                      # these are total number of vaccinees with unknown immunity
                      # but pre-ajust for time lag from vaccination to immune protection
                      VE1 = VE1, VE2=VE2 # Vaccine efficacy, need further adjustment by prior immunity 
        ) # include the delay reporting etc.
      } else if (epi.model == 'SEIRSVimmLoss'){
        
        daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
        
        if(nrow(daVacc.t)<1){  # no data yet
          V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
        } else { # yes data
          
          daVacc.t$date = daVacc.t$date %>% as.Date
          
          # make sure it includes a full week
          dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
          daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
          daVacc.t[is.na(daVacc.t)] = 0
          V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
          V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
          
          # print('start vacc!')
          
        }
        
        # get vaccine-induced immunity related parameter
        parmVimmLoss.t = get_parmVimmLoss(date.t = vdate.t, 
                                          date.delta = date.delta, #  as.Date('2021/07/01'), 
                                          date.omicron = date.omicron,
                                          # wt, alpha
                                          tm.imm.wt = tm.imm.wt, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.wt = p.imm.wane.max.wt, # maximal level of immunity loss (=1 or lower)
                                          k.wt = k.wt, # logistic function tuning parameter
                                          
                                          tm.imm.delta = tm.imm.delta, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.delta = p.imm.wane.max.delta, # maximal level of immunity loss (=1 or lower)
                                          k.delta = k.delta, # logistic function tuning parameter
                                          
                                          tm.imm.omicron = tm.imm.omicron, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.omicron = p.imm.wane.max.omicron, # maximal level of immunity loss (=1 or lower)
                                          k.omicron = k.omicron # logistic function tuning parameter
        )
        simEpi=SEIRSVimmLoss(tm_strt, tm_end, tm_step=1, # 1 day time-step
                             tmstep = tmstep,
                             state0 = state0,
                             S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                             I0=state0[paste0('I',1:num_gr),], 
                             beta=beta_tt, 
                             Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                             seed=seed, stoch=stoch, 
                             severity = severity.t,
                             newI.previous = newI.previous,
                             dist_tm.to.detect = dist_tm.to.detect,
                             dist_tm.to.death = dist_tm.to.death,
                             percSmax.t = percSmax.t,
                             V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                             # these are total number of vaccinees with unknown immunity
                             # but pre-ajust for time lag from vaccination to immune protection
                             VE1 = VE1, VE2=VE2, # Vaccine efficacy, need further adjustment by prior immunity 
                             newVacc = newVacc.t, # prior vaccinated & protected from infection
                             parmVimmLoss = parmVimmLoss.t, # parameters related to the immunity period of vaccine, against infection
                             R0 = Rimm.t   # to check those recovered from infection and remaining immune to infection
        ) # include the delay reporting etc.
        
        # update hist newVacc.t, if it's not just the prior
        if(tt < num_times_i){
          cumVimmLoss.t = simEpi$cumVimmLoss # for computing changes in S
          newVacc.t = simEpi$newVacc
          Rimm.t  = simEpi$R
          cumIimmLoss.t = simEpi$cumIimmLoss # %>% c # convert it to vector
          ts_cumVimmLoss = rbind(ts_cumVimmLoss, cumVimmLoss.t)
          ts_cumIimmLoss = rbind(ts_cumIimmLoss, cumIimmLoss.t)
          ts_R = rbind(ts_R, Rimm.t)
        }
        simEpi$newVacc = NULL
        simEpi$cumVimmLoss = NULL
        simEpi$R = NULL
        simEpi$cumIimmLoss = NULL

      }
      
      # re-assemble to the same order as the prior: state0
      n.end = tm_end - tm_strt + 2
      state.new = NULL
      for(i in 1:(length(simEpi)-1)){
        tmp = simEpi[[i]][n.end,,drop=F]; 
        rownames(tmp)=gsub('cumI','newI',paste0(names(simEpi)[i],1:num_gr))
        state.new = rbind(state.new,tmp)
      }
      
      state.new = rbind(state.new, state0[parm.names,])
      state.new = state.new[rownames(state0),] # make sure the order is the same
      
      
      xprior[,,tt+1]= state.new
      xprior.daily[1:tmstep+tt*tmstep, ]=simEpi$daily.newItot # we want the daily total new cases, without delay, without under-report
      
      # save prior for E and I
      Eprior[tt+1] = state.new['E1',] %>% mean
      Iprior[tt+1] = state.new['I1',] %>% mean
    } # end get prior
    
  } # end for-loop 
  
  
  # calculate the mean of ensemble
  xprior_mean=xpost_mean=xprior_sd=xpost_sd=matrix(0,num_times_i,num_var)
  for (tt in 1:num_times_i){
    xprior_mean[tt,]=apply(xprior[,,tt],1,mean, na.rm=T)
    xprior_sd[tt,]=apply(xprior[,,tt],1,sd, na.rm=T)
    xpost_mean[tt,]=apply(xpost[,,tt],1,mean, na.rm=T)
    xpost_sd[tt,]=apply(xpost[,,tt],1,sd, na.rm=T)
    
  }
  colnames(xprior_mean)=colnames(xpost_mean)=colnames(xpost_sd)=colnames(xprior_sd)=state.names
  xpost_mean = data.table(Week.start = Week.starts, xpost_mean)
  xprior_mean = data.table(Week.start = Week.starts, xprior_mean)
  
  
  
  # get stats instead
  dimnames(xpost)[1] = list(state.names)
  states_stats = NULL
  for(var in state.names){
    tmp = xpost[var,,]
    tmp = tmp %>% apply(2, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
    colnames(tmp) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
    
    tmp.mean = xpost[var,,] %>% apply(2, mean)
    states_stats = rbind(states_stats, 
                         data.table(Week.start = Week.starts, state = var, mean = tmp.mean, tmp))
  }
  
  # add imm loss
  tmpI = ts_cumIimmLoss %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
  tmpV = ts_cumVimmLoss %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
  immLoss_stats = rbind(data.table(Week.start = Week.starts, state = 'IimmLoss', mean = rowMeans(ts_cumIimmLoss), tmpI), 
                        data.table(Week.start = Week.starts, state = 'VimmLoss', mean = rowMeans(ts_cumVimmLoss), tmpV))
  colnames(immLoss_stats) = c('Week.start', 'state', 'mean', 'median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
  
  # compute Rt
  Rt_ens = matrix(0, num_times_i, num_ens)
  R0_ens = matrix(0, num_times_i, num_ens)
  Rtx_ens = matrix(0, num_times_i, num_ens)
  for(tt in 1:num_times_i){
    
    cur.wk = weeks[tt]
    
    for(ii in 1:num_ens){
      
      beta.mean = xpost['beta',ii,tt] * fn_getRelMob(rel.mob[tt],xpost['p.mob',ii,tt]) 
      # add seasonality? 
      if(seasonality) {
        # beta.mean =beta.mean * seasonal.cycle[seasonal.cycle$week==weeks[tt],]$relR0
        beta.mean = beta.mean * relR0[cur.wk,ii] # * rel.mob[tt]
      }
      PARMS = list(beta.mean = beta.mean,  Tir.mean= xpost['Tir',ii,tt], S = xpost['S1',ii,tt], N = N)
      Rt_ens[tt, ii] = Fn_getRt_SEIR(PARMS)
      
      PARMS = list(beta.mean = xpost['beta',ii,tt] * ifelse(seasonality, relR0[cur.wk,ii], 1),  Tir.mean= xpost['Tir',ii,tt], S = xpost['S1',ii,tt], N = N) # it'd be the same as Rtx if use this parm list
      R0_ens[tt, ii] = Fn_getR0_SEIR(PARMS) # # cp Rtx: additionally account for seasonality
      
      Rtx_ens[tt, ii] = xpost['beta',ii,tt] * xpost['Tir',ii,tt] # no seasonality
    } # ens
  }
  
  Rt_stats = Rt_ens %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
  colnames(Rt_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
  
  Rt_stats = data.table(Week.start = Week.starts, 
                        mean = Rt_ens %>% apply(1, mean),
                        sd = Rt_ens %>% apply(1, sd),
                        Rt_stats)
  
  R0_stats = R0_ens %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
  colnames(R0_stats) = c('median', 'iqr.lwr', 'iqr.upr','ci95.lwr','ci95.upr')
  
  R0_stats = data.table(Week.start = Week.starts, 
                        mean = R0_ens %>% apply(1, mean),
                        sd = R0_ens %>% apply(1, sd),
                        R0_stats)
  
  Rtx_stats = Rtx_ens %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
  colnames(Rtx_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
  
  Rtx_stats = data.table(Week.start = Week.starts, 
                         mean = Rtx_ens %>% apply(1, mean),
                         sd = Rtx_ens %>% apply(1, sd),
                         Rtx_stats)
  
  xpost.last=xpost[,,tt] # save all the ens members
  row.names(xpost.last)=state.names
  
  cumIperc = NULL # place holders
  cumIperc_stats = NULL
  Susceptibility = NULL
  Susceptibility_stats = NULL
  if(!save.cumIetc)
    xpost.last = NULL # also don't save the last xpost ens 
  if(save.cumIetc){ # don't need this, don't save unless specified to reduce file size
    # cumulative infection
    cumIperc = xpost['newItot1',,] %>% t %>% apply(2, cumsum)
    cumIperc = cumIperc / N * 100 # %
    cumIperc_stats = cumIperc %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
    colnames(cumIperc_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
    cumIperc_stats = data.table(Week.start = Week.starts, 
                                mean = cumIperc %>% apply(1, mean),
                                sd = cumIperc %>% apply(1, sd),
                                cumIperc_stats)
    # for susceptible
    Susceptibility = xpost['S1',,] %>% t 
    Susceptibility = Susceptibility / N * 100 # %
    Susceptibility_stats = Susceptibility %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
    colnames(Susceptibility_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
    Susceptibility_stats = data.table(Week.start = Week.starts, 
                                      mean = Susceptibility %>% apply(1, mean),
                                      sd = Susceptibility %>% apply(1, sd),
                                      Susceptibility_stats)
  }
  
  # do projection
  # integrate forward 1 step to get the prior
  # baseline projection
  p.beta.summer.holiday = 1.2 # assuming during the 7/4 holi transmission is 5% higher
  p.beta.after.holiday = .9 # assume estimates made during the holiday is higher and need to be adjusted
  print('do projection sce asIs...')
  {
    num_wk.fcast= length(weeks.fcast)
    fcast.inf=matrix(0,num_wk.fcast,num_ens)
    fcast.case=matrix(0,num_wk.fcast,num_ens)
    fcast.death=matrix(0,num_wk.fcast,num_ens)
    seed0 = seed
    nwk.train = Week.starts[!Week.starts %in% fcast.wk.starts]
    wk.train.end = nwk.train %>% tail(1) %>% as.Date
    newI.previous = newI.previous[1:(length(nwk.train)*7),]
    
    if(useMean){ # use the mean for projection
      means = state0.fcast %>% apply(1, mean) 
      sds = state0.fcast %>% apply(1, sd) 
      spreads = state0.fcast %>% apply(1, quantile, prob = c(.4, .5, .6)) %>% t
      pleft = spreads[,2] / spreads[,1] 
      pright = spreads[,3] / spreads[,2]
      tmp = cbind(left = pleft, right = pright, p = pleft / pright)
      
      means = cbind(means - sds * useMean.wd, means + sds * useMean.wd * tmp[,'p'])  # +/- x% of sd, adjust for asymmetry 
      state0.fcast.mean = lhs(num_ens, rect = means) %>% t
      rownames(state0.fcast.mean) = rownames(state0.fcast)
      
      state0 = state0.fcast.mean
      Rimm.t = Rimm.t.fcast
      newVacc.t = newVacc.t.fcast
      
    } else if(exclExtreme) { # exclude extreme values before projection, this could help with long term projection - less likely to deplete the susceptible
      # most importantly, susceptibility, beta, Tir, maybe Tei, pmob, alpha, ifr
      # sample as a whole to reserve the parameter connections
      state0.fcast.redn = state0.fcast
      Rimm.t.fcast.redn = Rimm.t.fcast
      extremes = state0.fcast %>% apply(1, quantile, prob = c(perc2excl/2, 1 - perc2excl/2)) %>% t
      idx = which(!(state0.fcast['S1',] > extremes['S1',1] & state0.fcast['S1',] < extremes['S1',2]) |
                    !(state0.fcast['beta',] > extremes['beta',1] & state0.fcast['beta',] < extremes['beta',2]) |
                   !(state0.fcast['Tir',] > extremes['Tir',1] & state0.fcast['Tir',] < extremes['Tir',2]) |
                   !(state0.fcast['alpha',] > extremes['alpha',1] & state0.fcast['alpha',] < extremes['alpha',2]) |
                   !(state0.fcast['p.mob',] > extremes['p.mob',1] & state0.fcast['p.mob',] < extremes['p.mob',2]) # |
                   # !(state0.fcast['ifr',] > extremes['ifr',1] & state0.fcast['ifr',] < extremes['ifr',2])
      )
      idx = which(!(state0.fcast['E1',] > extremes['E1',1] & state0.fcast['E1',] < extremes['E1',2]) |
                    !(state0.fcast['I1',] > extremes['I1',1] & state0.fcast['I1',] < extremes['I1',2])
      )
      tmp = state0.fcast[,idx]
      smp.all = 1:num_ens
      smp.set = smp.all[!smp.all %in% idx]
      idx.new = sample(smp.set, length(idx), replace = F)
      state0.fcast.redn[,idx] = state0.fcast.redn[,idx.new]
      Rimm.t.fcast.redn[,idx] = Rimm.t.fcast[,idx.new]
      state0 = state0.fcast.redn
      Rimm.t = Rimm.t.fcast.redn
      newVacc.t = newVacc.t.fcast
      
      print(paste('exclude',length(idx),'extreme ens'))
    } else { # use the last training estimate for projection
      state0 = state0.fcast
      Rimm.t = Rimm.t.fcast
      newVacc.t = newVacc.t.fcast # tracking immune loss in vaccinees
    }
    
    fcast.state0 = state0
    
    fcast.mob = fn_get.fcast.mob(fcast.wk.starts = fcast.wk.starts, 
                                weeks.fcast = weeks.fcast, da.mob.t = da.mob.t,
                                num_wk.fcast = num_wk.fcast)
    
     
    
    fda.vacc = da.vacc
    
    source(paste0(dir_code,'set_tm2event.R')) # reset the time to event variables
    
    for(iwk in 1:num_wk.fcast){
      cur.wk = weeks.fcast[iwk]
      tm_strt = 1; tm_end = tmstep;
      vdate.t = fcast.wk.starts[iwk];
      wk.train.t = which(Week.starts == fcast.wk.starts[iwk]) # to match the wk of training and get the posterior estimate of ifr and detection rate
      flag.higher.det.rate = F
      flag.higher.ifr = F
      if(length(wk.train.t)<1){ # if the historical estimates did not cover the entire projection period
        wk.train.t = length(Week.starts) # use the last estimates
        # there is no estimate for detection rate
        # allow higher detection rate if it's the winter & the last estimate was very low
        # from 11/15 to 1/15
        if((as.Date(vdate.t)) %>% format('%m') %>% as.numeric() %in% c(11, 12, 1) & as.Date(vdate.t) > as.Date('2022/10/1') &
            # also only when the last estimate was generated during summer/fall and likely low
           as.Date(tail(nwk.train,1)) %>% format('%m') %>% as.numeric %in% 5:10
           )
          flag.higher.det.rate = T
        # similarly for ifr
        if((as.Date(vdate.t)) %>% format('%m') %>% as.numeric() %in% c(12, 1, 2, 3) & as.Date(vdate.t) > as.Date('2022/10/1') &
           # also only when the last estimate was generated during summer/fall and likely low
           as.Date(tail(nwk.train,1)) %>% format('%m') %>% as.numeric %in% 6:11
        )
          flag.higher.ifr = T
      }
      
      # shrink the spread? only the state variables that change as the trajectories unfold, or for parms as well?
      if(fcast.deflat !=1){
        xmn=rowMeans(state0); 
        
        # only adjust the state variables
        mat.deflat = matrix(1, nrow = nrow(state0), ncol = ncol(state0))
        rownames(mat.deflat) = rownames(state0)
        mat.deflat[v.deflat,] = fcast.deflat # only adjust the state variables
        
        xnew=mat.deflat*(state0-xmn%*%matrix(1,1,num_ens))+xmn%*%matrix(1,1,num_ens); 
        
        state0 = xnew
      }
      
      # set seeding based on national nevel case data and local mobility
      seed.t = pmin(1, pmax(.1, case.us[date == as.Date(vdate.t)]$US / 7 * .005)) * ifelse(mean(fcast.mob[iwk,])>1, mean(fcast.mob[iwk,])*2, mean(fcast.mob[iwk,]))
      if(length(seed.t)<1) # no data
        seed.t = tail(rowMeans(fcast.mob), 1)
      # print(paste(vdate.t, round(seed.t,2)), quote = F)
      # get VE depending time period / circulating variant 
      tmp = getVE(date.t = vdate.t,
                  date.delta = date.delta, #  as.Date('2021/07/01'), 
                  date.omicron = date.omicron,
                  VE1wt = VE1wt,
                  VE2wt = VE2wt,  # higher b/c we are using mortality data too
                  VE1delta = VE1delta,
                  VE2delta = VE2delta,  # higher b/c we are using mortality data too
                  VE1omicron = VE1omicron,
                  VE2omicron = VE2omicron
      )
      VE1 = tmp$VE1; VE2 = tmp$VE2
      
      if(!is.null(newI.previous) & vdate.t >= vax.start){
        tm.t = nrow(newI.previous)
        # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
        # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
        # 2/5/21 set t1 to 1 yr given the slow rollout
        t1 = 365
        cumI.t = colSums(newI.previous[pmax(1,tm.t-t1) : (tm.t),]) #  %>% apply(1, median) # excl last two weeks and get the median
        # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
        # higher infection rate for the priority group
        # tm.t = pmax(1, tm.t - t1 + 1) # re-aline timing with start of vac
        tm.imm = 365* 2.5 # assume 3 yr immunity
        p.imm.wane.max = .8; k = .015  # 1/tm.imm  
        p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(tm.t + 60 - tm.imm/2))) # not all infected from day 1
        # earlier infections are likely to be in the high priority groups 
        p.imm = 1 *  p.imm.wane * redn.priority # assume 100% prior infection provide immunity, but wane over time
        # and multiple by % excluded if there is prior testing before vax: p.prior.test
        percSmax.t = 1 - cumI.t / N * p.imm
        # also compare to current susceptibility, in case of immune evasion that increases susceptiblity
        percSmax.t = pmax(percSmax.t, state0[paste0('S',1:num_gr),] / N)
        # no lower than 50%, in case of outliers
        percSmax.t = pmax(percSmax.t, .3)
        # print(c('cohort %S:',round(summary(mean(percSmax.t)),2)), quote = F)
      } else {
        percSmax.t = 0
        # print('no vax yet')
      }
      
      beta_tt = state0['beta',] * fn_getRelMob(fcast.mob[iwk,],state0['p.mob',]) 
      if(seasonality) {
        beta_tt = beta_tt * relR0[cur.wk,] 
      } 
      # add holiday eff
      if(as.Date(vdate.t) %in% seq(as.Date(paste0(format(as.Date(vdate.t), '%Y'),'/6/25')), length.out = 30, by = 'day') |
         as.Date(wk.train.end) %in% seq(as.Date(paste0(format(as.Date(vdate.t), '%Y'),'/7/1')), length.out = 30, by = 'day')
         )
        beta_tt = beta_tt * fn_get.p.beta.summer.holiday(vdate.t, wk.train.end, p.beta.summer.holiday, p.beta.after.holiday)
      
      # update severity
      severity.t['death',] = xpost['ifr',,wk.train.t]
      if(flag.higher.ifr)
        severity.t['death',] = fn_adj.parm4winter(parm.t = severity.t['death',], # parm to adjust
                                   parm.name.t = 'ifr',
                                   parm.low.cut = 1, # set it to a high level so it is invoked
                                   parm.winter.amp = .25, parm.phase.shift = -30,
                                   vdate.t = vdate.t)
      # update detection rate
      state0['alpha',] = xpost['alpha',,wk.train.t]
      # higher det rate in the winter
      if(flag.higher.det.rate)
        state0['alpha',] = fn_adj.parm4winter(parm.t = state0['alpha',], # parm to adjust
                                   parm.name.t = 'alpha',
                                   parm.low.cut = .14, parm.winter.amp = .25, parm.phase.shift = 10,
                                   vdate.t = vdate.t)
      
      
      # [2/23/21] check this if problem!
      dist_tm.to.detect = NULL
      for(ii in 1:num_ens){
        tmp = generation.time(dist_tm.to.detect.name,c(state0['Td.mean',ii],state0['Td.sd',ii]),truncate = tm.to.detect.max)
        dist_tm.to.detect=cbind(dist_tm.to.detect,tmp$GT[-1]); 
      }
      
      # [1/20/22] also update time from infection to death, b/c longer time-lag for omicron
      if((vdate.t >= (date.omicron+7) & vdate.t < (date.omicron+14)) |
         # or first week of fcast and it hasn't been updated
         (iwk == 1 & vdate.t >= (date.omicron+7))
         ){
        print('update tm to death for omicron')
        
        tm.to.death.max = tm.to.death.max_omicron
        if(tm.to.deathFromDiag){
          tm.from.inf.to.death.max = tm.to.diag + tm.to.death.max
        } else {
          tm.from.inf.to.death.max = tm.to.death.max
        }
        
        dist_tm.to.death = NULL  # time to death
        for(ii in 1:num_ens){
          tmp = generation.time(dist_tm.to.death.name,c(tm.to.outcome.mn['tm.to.death',ii] + tm2death_adj_omicron, 
                                                        tm.to.outcome.sd['tm.to.death',ii] # not good to increase sd
          ),
          truncate = tm.to.death.max # updated above
          )
          dist_tm.to.death=cbind(dist_tm.to.death,tmp$GT[-1]); 
        }
        if(tm.to.deathFromDiag){ # if the distribution of time to death is from diagnosis, not infectious
          # add time from infectious to diagnosis
          dist_tm.to.death = rbind(matrix(0,round(tm.to.diag,0),num_ens),dist_tm.to.death)
        }
      }
      
      if(epi.model == 'SEIRS'){
        
        simEpi=SEIRS(tm_strt, tm_end, tm_step=1, # 1 day time-step
                     tmstep = tmstep,
                     state0 = state0,
                     S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                     I0=state0[paste0('I',1:num_gr),], 
                     beta=beta_tt, 
                     Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                     seed=seed.t, stoch=stoch,
                     severity = severity.t,
                     newI.previous = newI.previous,
                     dist_tm.to.detect = dist_tm.to.detect,
                     dist_tm.to.death = dist_tm.to.death)
        
      } else if(epi.model == 'SEIRSV'){
        daVacc.t = fda.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
        
        if(nrow(daVacc.t)<1){  # no data yet
          V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
        } else { # yes data
          
          daVacc.t$date = daVacc.t$date %>% as.Date
          
          # make sure it includes a full week
          dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
          daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
          daVacc.t[is.na(daVacc.t)] = 0
          V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
          V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
          
          # print('start vacc!')
          
        }
        simEpi=SEIRSV(tm_strt, tm_end, tm_step=1, # 1 day time-step
                      tmstep = tmstep,
                      state0 = state0,
                      S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                      I0=state0[paste0('I',1:num_gr),], 
                      beta=beta_tt, 
                      Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                      seed=seed.t, stoch=stoch, 
                      severity = severity.t,
                      newI.previous = newI.previous,
                      dist_tm.to.detect = dist_tm.to.detect,
                      dist_tm.to.death = dist_tm.to.death,
                      percSmax.t = percSmax.t,
                      V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                      # these are total number of vaccinees with unknown immunity
                      # but pre-ajust for time lag from vaccination to immune protection
                      VE1 = VE1, VE2=VE2 # Vaccine efficacy, need further adjustment by prior immunity 
        ) # include the delay reporting etc.
      } else if (epi.model == 'SEIRSVimmLoss'){
        
        daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
        
        if(nrow(daVacc.t)<1){  # no data yet
          V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
        } else { # yes data
          
          daVacc.t$date = daVacc.t$date %>% as.Date
          
          # make sure it includes a full week
          dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
          daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
          daVacc.t[is.na(daVacc.t)] = 0
          V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
          V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
          
          # print('start vacc!')
          
        }
        
        # get vaccine-induced immunity related parameter
        parmVimmLoss.t = get_parmVimmLoss(date.t = as.Date(vdate.t), 
                                          date.delta = date.delta, #  as.Date('2021/07/01'), 
                                          date.omicron = date.omicron,
                                          # wt, alpha
                                          tm.imm.wt = tm.imm.wt, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.wt = p.imm.wane.max.wt, # maximal level of immunity loss (=1 or lower)
                                          k.wt = k.wt, # logistic function tuning parameter
                                          
                                          tm.imm.delta = tm.imm.delta, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.delta = p.imm.wane.max.delta, # maximal level of immunity loss (=1 or lower)
                                          k.delta = k.delta, # logistic function tuning parameter
                                          
                                          tm.imm.omicron = tm.imm.omicron, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.omicron = p.imm.wane.max.omicron, # maximal level of immunity loss (=1 or lower)
                                          k.omicron = k.omicron # logistic function tuning parameter
        )
        simEpi=SEIRSVimmLoss(tm_strt, tm_end, tm_step=1, # 1 day time-step
                             tmstep = tmstep,
                             state0 = state0,
                             S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                             I0=state0[paste0('I',1:num_gr),], 
                             beta=beta_tt, 
                             Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                             seed=seed.t, stoch=stoch, 
                             severity = severity.t,
                             newI.previous = newI.previous,
                             dist_tm.to.detect = dist_tm.to.detect,
                             dist_tm.to.death = dist_tm.to.death,
                             percSmax.t = percSmax.t,
                             V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                             # these are total number of vaccinees with unknown immunity
                             # but pre-ajust for time lag from vaccination to immune protection
                             VE1 = VE1, VE2=VE2, # Vaccine efficacy, need further adjustment by prior immunity 
                             newVacc = newVacc.t, # prior vaccinated & protected from infection
                             parmVimmLoss = parmVimmLoss.t, # parameters related to the immunity period of vaccine, against infection
                             R0 = Rimm.t   # to check those recovered from infection and remaining immune to infection
        ) # include the delay reporting etc.
        
        # update hist newVacc.t
        newVacc.t = simEpi$newVacc
        Rimm.t  = simEpi$R
        simEpi$newVacc = NULL
        simEpi$cumVimmLoss = NULL
        simEpi$R = NULL
        simEpi$cumIimmLoss = NULL
      }
      
      # re-assemble to the same order as the prior: state0
      n.end = tm_end - tm_strt + 2
      state.new = NULL
      for(i in 1:(length(simEpi)-1)){
        tmp = simEpi[[i]][n.end,,drop=F]; 
        rownames(tmp)=gsub('cumI','newI',paste0(names(simEpi)[i],1:num_gr))
        state.new = rbind(state.new,tmp)
      }
      
      # state.new = rbind(state.new, state0[parm.names,])
      state.new = rbind(state.new, fcast.state0[parm.names,]) # restore the parameters, so the deflation doesn't accumulate
      state.new = state.new[rownames(state0),] # make sure the order is the same
      
      state0 = state.new # update
      newi = simEpi$cumItot
      newi = newi[-1,] - newi[-nrow(newi),]
      newI.previous = rbind(newI.previous,newi) # include the most recent week
      fcast.inf[iwk,] = state.new[idx.newItot,]
      fcast.case[iwk,] = state.new[idx.obs_i,]
      fcast.death[iwk,] = state.new[idx.obs_d,]
      
      if(any(state0 < 0))
        break
      
    } # end week
    
    # get summary stats
    prob_vec = c(.5, .25, .75, .025, .975, .05, .95, .1, .9)
    prob.name_vec = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr','ci90.lwr','ci90.upr','ci80.lwr','ci80.upr')
    tmp.inf = fcast.inf %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.inf) = prob.name_vec
    tmp.case = fcast.case %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.case) = prob.name_vec
    tmp.death = fcast.death %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.death) = prob.name_vec
    
    fcast_stats = rbind(data.table(measure = 'Infections', Week.start=fcast.wk.starts,tmp.inf),
                        data.table(measure = 'Cases', Week.start=fcast.wk.starts,tmp.case),
                        data.table(measure = 'Deaths', Week.start=fcast.wk.starts,tmp.death)
    )
    
    # cumulative
    # get summary stats
    tmp.inf = fcast.inf %>% apply(2,cumsum) %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.inf) = prob.name_vec
    tmp.case = fcast.case %>% apply(2,cumsum) %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.case) = prob.name_vec
    tmp.death = fcast.death %>% apply(2,cumsum) %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.death) = prob.name_vec
    fcast_stats = rbind(fcast_stats, data.table(measure = 'Cumulative Infections', Week.start=fcast.wk.starts,tmp.inf),
                        data.table(measure = 'Cumulative Cases', Week.start=fcast.wk.starts,tmp.case),
                        data.table(measure = 'Cumulative Deaths', Week.start=fcast.wk.starts,tmp.death)
    )
    
    # get probability dist for computing log score
    fcastDist = fn_getProbDist(fcast.case, fcast.death, bins.case, bins.death)
    
    fcast_stats_base = fcast_stats
    fcastDist_base = fcastDist
    
    fcast_stats_base$scenario = 'asIs'
    fcastDist_base$scenario = 'asIs'
  } # end projection
  
  print('do projection sce newV...')
  if(incl.sce.newV){
    num_wk.fcast= length(weeks.fcast)
    fcast.inf=matrix(0,num_wk.fcast,num_ens)
    fcast.case=matrix(0,num_wk.fcast,num_ens)
    fcast.death=matrix(0,num_wk.fcast,num_ens)
    seed0 = seed
    nwk.train = Week.starts[!Week.starts %in% fcast.wk.starts]
    wk.train.end = nwk.train %>% tail(1) %>% as.Date
    newI.previous = newI.previous[1:(length(nwk.train)*7),]
    Sfcast = NULL
    
    # check if there is a recent major update due to a new virus already
    # based on major changes in S and maybe beta as well
    # wk.train.end = nwk.train %>% tail(1) %>% as.Date
    states.train = xpost_mean[1:length(nwk.train),]
    dS.train = states.train$S1[-1] - states.train$S1[-nrow(states.train)] + states.train$newItot1[-1]
    # only count increases
    i0 = which(dS.train > .05 * N) %>% head(1)  # first major increase
    imm0 = N - xpost_mean$S1[i0-1]  # cumulative infections / immunity
    dS.train[dS.train<0] = 0
    dtx.train = states.train$beta[-1] - states.train$beta[-nrow(states.train)]
    beta0 = states.train %>% filter(Week.start < as.Date('2020/9/1')) %>% .$beta %>% mean
    if(is.na(beta0))
      beta0 = states.train$beta %>% head(5) %>% mean
    # cum.dImm.train = cumsum(dS.train) / imm0 # relative to prior cumulative immunity, rather than N
    # recent.major.S.incr = any(tail(cum.dS.train, 4 * 3) > .3)
    # look at each 1 month sliding window
    flag.recent.major.dimm = F # whether there is a major immune erosion prior to the fcast
    width.train = ifelse(wk.train.end > as.Date('2022/2/28'), 6, 6) # 6 weeks or 4 weeks 
    cum.dimm.cut = ifelse(wk.train.end > as.Date('2022/2/28'), .25, .25)  
    hist.idx = length(dS.train) :1
    if(max(hist.idx) < width.train){
      width.train = max(hist.idx)
    }
    hist.idx = hist.idx[hist.idx >= width.train]
    # at least include all available hist data, if training period is short
    
    cum.dimm = NULL
    cum.dimm_vec = NULL
    if(length(hist.idx) > 0){
      for(i in hist.idx){
        # print(i)
        cum.dimm.t = (sum(dS.train[i:(i-width.train+1)])) / pmax((N - xpost_mean$S1[i-width.train+1]), sum(xpost_mean$newItot1[1:i])) # in case previous increase has increase S
        cum.dimm = c(cum.dimm, cum.dimm.t)
        cum.dimm_vec = rbind(cum.dimm_vec, data.table(wk = i, cum.dimm = cum.dimm.t)) # confusing time - record this to double check - looks correct
      }
      
      # only look at the recent 2 months
      flag.recent.major.dimm = any(head(cum.dimm, width.train) > cum.dimm.cut)
    }
    
    
    # tx change
    flag.recent.major.dtx = F # whether there is a major tx incr prior to the fcast
    width.train = ifelse(wk.train.end > date.omicron, 8, 10) # 6 weeks or 4 weeks 
    cum.dtx.cut = ifelse(wk.train.end > date.omicron, .3, .4)  
    hist.idx = length(dS.train) :1
    if(max(hist.idx) < width.train){
      width.train = max(hist.idx)
    }
    hist.idx = hist.idx[hist.idx >= width.train] 
    cum.dtx = NULL
    if(length(hist.idx) > 0){
      for(i in hist.idx){
        cum.tx.t = sum(dtx.train[i:(i-width.train+1)]) / beta0
        cum.dtx = c(cum.dtx, cum.tx.t)
      }
      # only look at the recent 2 months
      flag.recent.major.dtx = any(head(cum.dtx, width.train) > cum.dtx.cut)
    }
    flag.recent.major.newV = F
    # whether the recent changes due to a newV has been accounted for
    flag.recent.major.newV = flag.recent.major.dimm & flag.recent.major.dtx
    
    # 9/14/22 check if there is rising new variant based on AVAILABLE genomic data
    # assume a 2 week delay in data collection
    tmp.newv = fn_checkNewV(DAT.VARIANT, # variant data
                 loc.t, # location
                 fcast.start)
    flag.newVrising = tmp.newv$flag.newVrising
    newv.pred = tmp.newv$newv.pred
    # print(newv.pred)
    
    # 8/10/22 NEW VARIANT SCENARIO 1 - based on the real historical record before the forecast only
    # WAS THERE RECENT MASS WAVE WITH LARGE INFECTION?
    # not necessarily recent 
    # 2 months sliding window
    flag.mass.wave.hist = F # whether there is a massive wave recently (def: 1/3 infected w/i 2 months)
    cumi.mean = NULL
    cumi.mean_vec = NULL # for quality check
    # width = ifelse(variant.tag == 'Omicron', 60, 90) # 2 months
    # cumi.cut = ifelse(variant.tag == 'Omicron', 1/3, 1/4) 
    width = ifelse(wk.train.end > date.omicron, 60, 90) # 2 months
    cumi.cut = ifelse(wk.train.end > date.omicron, 1/3, 1/4) 
    # hist.idx = nrow(newI.previous) : (nrow(newI.previous) - width+1)
    hist.idx = nrow(newI.previous) :1
    if(max(hist.idx) < width){
      width = max(hist.idx)
    }
    hist.idx = hist.idx[hist.idx >= width]
    if(length(hist.idx) > 0){
      for(i in hist.idx){
        cumi.t = (colSums(newI.previous[i:(i-width+1),]) %>% mean) / N
        cumi.mean = c(cumi.mean, cumi.t)
        cumi.mean_vec = rbind(cumi.mean_vec, data.table(day = i, cumi.mean = cumi.t)) # looks correct
      }
      
      # only look at the recent 2 months
      flag.mass.wave.hist = any(head(cumi.mean, width) > cumi.cut)
    }
    
    cum.dS.fcast = 0 # keeping check how much S has been increased when assuming immune erosion
    cum.dS.fcast.cut = .2 # relative to population size, .25 IS TOO HIGH
    cum.dimm.recentNfcast = cum.dimm[1]  # include the recent 1 month
    cum.dimm.recentNfcast.cut = ifelse(wk.train.end > date.omicron, .6, .4) # max increase including recent ~1 month, and relative to the immunity, rather than N
    
    if(useMean){ # use the mean for projection
      means = state0.fcast %>% apply(1, mean) 
      sds = state0.fcast %>% apply(1, sd) 
      spreads = state0.fcast %>% apply(1, quantile, prob = c(.4, .5, .6)) %>% t
      pleft = spreads[,2] / spreads[,1] 
      pright = spreads[,3] / spreads[,2]
      tmp = cbind(left = pleft, right = pright, p = pleft / pright)
      
      means = cbind(means - sds * useMean.wd, means + sds * useMean.wd * tmp[,'p'])  # +/- x% of sd, adjust for asymmetry 
      state0.fcast.mean = lhs(num_ens, rect = means) %>% t
      rownames(state0.fcast.mean) = rownames(state0.fcast)
      
      state0 = state0.fcast.mean
      Rimm.t = Rimm.t.fcast
      newVacc.t = newVacc.t.fcast
      
    } else if(exclExtreme) { # exclude extreme values before projection, this could help with long term projection - less likely to deplete the susceptible
      # most importantly, susceptibility, beta, Tir, maybe Tei, pmob, alpha, ifr
      # sample as a whole to reserve the parameter connections
      state0.fcast.redn = state0.fcast
      Rimm.t.fcast.redn = Rimm.t.fcast
      extremes = state0.fcast %>% apply(1, quantile, prob = c(perc2excl/2, 1 - perc2excl/2)) %>% t
      idx = which(!(state0.fcast['S1',] > extremes['S1',1] & state0.fcast['S1',] < extremes['S1',2]) |
                    !(state0.fcast['beta',] > extremes['beta',1] & state0.fcast['beta',] < extremes['beta',2]) |
                    !(state0.fcast['Tir',] > extremes['Tir',1] & state0.fcast['Tir',] < extremes['Tir',2]) |
                    !(state0.fcast['alpha',] > extremes['alpha',1] & state0.fcast['alpha',] < extremes['alpha',2]) |
                    !(state0.fcast['p.mob',] > extremes['p.mob',1] & state0.fcast['p.mob',] < extremes['p.mob',2]) # |
                  # !(state0.fcast['ifr',] > extremes['ifr',1] & state0.fcast['ifr',] < extremes['ifr',2])
      )
      idx = which(!(state0.fcast['E1',] > extremes['E1',1] & state0.fcast['E1',] < extremes['E1',2]) |
                    !(state0.fcast['I1',] > extremes['I1',1] & state0.fcast['I1',] < extremes['I1',2])
      )
      tmp = state0.fcast[,idx]
      smp.all = 1:num_ens
      smp.set = smp.all[!smp.all %in% idx]
      idx.new = sample(smp.set, length(idx), replace = F)
      state0.fcast.redn[,idx] = state0.fcast.redn[,idx.new]
      Rimm.t.fcast.redn[,idx] = Rimm.t.fcast[,idx.new]
      state0 = state0.fcast.redn
      Rimm.t = Rimm.t.fcast.redn
      newVacc.t = newVacc.t.fcast
      
      print(paste('exclude',length(idx),'extreme ens'))
    } else { # use the last training estimate for projection
      state0 = state0.fcast
      Rimm.t = Rimm.t.fcast
      newVacc.t = newVacc.t.fcast # tracking immune loss in vaccinees
    }
    
    fcast.state0 = state0
    
    
    fcast.mob = fn_get.fcast.mob(fcast.wk.starts = fcast.wk.starts, 
                                 weeks.fcast = weeks.fcast, da.mob.t = da.mob.t,
                                 num_wk.fcast= num_wk.fcast)
    
    fda.vacc = da.vacc
    
    source(paste0(dir_code,'set_tm2event.R')) # reset the time to event variables
    
    for(iwk in 1:num_wk.fcast){
      cur.wk = weeks.fcast[iwk]
      tm_strt = 1; tm_end = tmstep;
      vdate.t = fcast.wk.starts[iwk];
      wk.train.t = which(Week.starts == fcast.wk.starts[iwk]) # to match the wk of training and get the posterior estimate of ifr and detection rate
      
      # check if for this prediction wk, it is detected that a new variant is rising and is still on its way to become predominant
      # allow adjustment up to 8 weeks (6 weeks after accounting for the lag) -> so it won't be over adjusted -> too much, reduce it to 6 (4 weeks)
      # as slow progressing new variants may be replaced by others before reaching high level
      flag.newVrising.t = flag.newVrising > 0 & iwk < (pmin(6, max(newv.pred$wk2newVreplace)) - ifelse(max(newv.pred$wk2newVreplace)>3, 1, 0)) # if more than 3 weeks, take off 1 week to account for the 2-wk delay a bit
      
      flag.higher.det.rate = F
      flag.higher.ifr = F
      if(length(wk.train.t)<1){ # if the historical estimates did not cover the entire projection period
        wk.train.t = length(Week.starts) # use the last estimates
        # there is no estimate for detection rate
        # allow higher detection rate if it's the winter & the last estimate was very low
        # from 11/15 to 1/15
        if((as.Date(vdate.t)) %>% format('%m') %>% as.numeric() %in% c(11, 12, 1) & as.Date(vdate.t) > as.Date('2022/10/1') &
           # also only when the last estimate was generated during summer/fall and likely low
           as.Date(tail(nwk.train,1)) %>% format('%m') %>% as.numeric %in% 5:10
        )
          flag.higher.det.rate = T
        # similarly for ifr
        if((as.Date(vdate.t)) %>% format('%m') %>% as.numeric() %in% c(12, 1, 2, 3) & as.Date(vdate.t) > as.Date('2022/10/1') &
           # also only when the last estimate was generated during summer/fall and likely low
           as.Date(tail(nwk.train,1)) %>% format('%m') %>% as.numeric %in% 6:11
        )
          flag.higher.ifr = T
      }
      
      # shrink the spread? only the state variables that change as the trajectories unfold, or for parms as well?
      if(fcast.deflat !=1){
        xmn=rowMeans(state0); 
        
        # only adjust the state variables
        mat.deflat = matrix(1, nrow = nrow(state0), ncol = ncol(state0))
        rownames(mat.deflat) = rownames(state0)
        mat.deflat[v.deflat,] = fcast.deflat # only adjust the state variables
        
        xnew=mat.deflat*(state0-xmn%*%matrix(1,1,num_ens))+xmn%*%matrix(1,1,num_ens); 
        
        state0 = xnew
      }
      
      # set seeding based on national nevel case data and local mobility
      seed.t = pmin(1, pmax(.1, case.us[date == as.Date(vdate.t)]$US / 7 * .005)) * ifelse(mean(fcast.mob[iwk,])>1, mean(fcast.mob[iwk,])*2, mean(fcast.mob[iwk,]))
      if(length(seed.t)<1) # no data
        seed.t = tail(rowMeans(fcast.mob), 1)
      # print(paste(vdate.t, round(seed.t,2)), quote = F)
      # tmp = getVE(date.t = vdate.t)  # get VE depending time period / circulating variant 
      tmp = getVE(date.t = vdate.t,
                  date.delta = date.delta, #  as.Date('2021/07/01'), 
                  date.omicron = date.omicron,
                  VE1wt = VE1wt,
                  VE2wt = VE2wt,  # higher b/c we are using mortality data too
                  VE1delta = VE1delta,
                  VE2delta = VE2delta,  # higher b/c we are using mortality data too
                  VE1omicron = VE1omicron,
                  VE2omicron = VE2omicron
      )
      VE1 = tmp$VE1; VE2 = tmp$VE2
      
      if(!is.null(newI.previous) & vdate.t >= vax.start){
        tm.t = nrow(newI.previous)
        # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
        # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
        # 2/5/21 set t1 to 1 yr given the slow rollout
        t1 = 365
        cumI.t = colSums(newI.previous[pmax(1,tm.t-t1) : (tm.t),]) #  %>% apply(1, median) # excl last two weeks and get the median
        # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
        # higher infection rate for the priority group
        # tm.t = pmax(1, tm.t - t1 + 1) # re-aline timing with start of vac
        tm.imm = 365* 2.5 # assume 3 yr immunity
        p.imm.wane.max = .8; k = .015  # 1/tm.imm  
        p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(tm.t + 60 - tm.imm/2))) # not all infected from day 1
        # earlier infections are likely to be in the high priority groups 
        p.imm = 1 *  p.imm.wane * redn.priority # assume 100% prior infection provide immunity, but wane over time
        # and multiple by % excluded if there is prior testing before vax: p.prior.test
        percSmax.t = 1 - cumI.t / N * p.imm
        # also compare to current susceptibility, in case of immune evasion that increases susceptiblity
        percSmax.t = pmax(percSmax.t, state0[paste0('S',1:num_gr),] / N)
        # no lower than 50%, in case of outliers
        percSmax.t = pmax(percSmax.t, .3)
        # print(c('cohort %S:',round(summary(mean(percSmax.t)),2)), quote = F)
      } else {
        percSmax.t = 0
        # print('no vax yet')
      }
      
      # 8/10/22 NEW VARIANT SCENARIO 1 - careful this could lead to a vicious cycle!
      # WAS THERE RECENT MASS WAVE WITH LARGE INFECTION?
      # not necessarily recent 
      # 2 months sliding window
      flag.mass.wave = F # whether there is a massive wave recently (def: 1/3 infected w/i 2 months)
      cumi.mean = NULL
      cumi.mean_vec = NULL # for quality check
      # width = ifelse(variant.tag == 'Omicron', 60, 90) # 2 months
      # cumi.cut = ifelse(variant.tag == 'Omicron', 1/3, 1/4) 
      width = ifelse(vdate.t > date.omicron, 60, 90) # 2 months
      cumi.cut = ifelse(vdate.t > date.omicron, 1/3, 1/4) 
      # hist.idx = nrow(newI.previous) : (nrow(newI.previous) - width+1)
      hist.idx = nrow(newI.previous) :1
      hist.idx = hist.idx[hist.idx >= width]
      if(length(hist.idx) > 0){
        for(i in hist.idx){
          cumi.t = (colSums(newI.previous[i:(i-width+1),]) %>% mean) / N
          cumi.mean = c(cumi.mean, cumi.t)
          cumi.mean_vec = rbind(cumi.mean_vec, data.table(day = i, cumi = cumi.t)) # checked
        }
        
        # only look at the recent 2 months
        flag.mass.wave = any(head(cumi.mean, width) > cumi.cut)
      }
      
      # look at only the forecast period
      flag.mass.wave.fcast = F # whether there is a massive wave recently (def: 1/3 infected w/i 2 months)
      cumi.mean = NULL
      cumi.mean_vec = NULL # for quality check
      # width = ifelse(variant.tag == 'Omicron', 60, 90) # 2 months
      # cumi.cut = ifelse(variant.tag == 'Omicron', 1/3, 1/4) 
      width = ifelse(vdate.t > date.omicron, 60, 90) # 2 months
      cumi.cut = ifelse(vdate.t > date.omicron, 1/3, 1/4) * .7
      # hist.idx = nrow(newI.previous) : (nrow(newI.previous) - width+1)
      # 10/18/22 bug
      # newI.previous.f = newI.previous[nrow(newI.previous) : (length(nwk.train) * 7),,drop = F] # from recent to older dates
      newI.previous.f = newI.previous[(length(nwk.train) * 7): nrow(newI.previous),,drop = F] # keep the same timeline, older on top, more recent in the bottom
      rownames(newI.previous.f) = (length(nwk.train) * 7) : nrow(newI.previous) # for quality check re time idx
      hist.idx = nrow(newI.previous.f) : 1
      hist.idx = hist.idx[hist.idx >= width]
      if(length(hist.idx) > 0){
        for(i in hist.idx){
          cumi.t = (colSums(newI.previous.f[i:(i-width+1),]) %>% median) / N
          cumi.mean = c(cumi.mean, cumi.t)
          cumi.mean_vec = rbind(cumi.mean_vec, data.table(day = i, cumi = cumi.t)) # checked
        }
        
        # only look at the recent 2 months
        flag.mass.wave.fcast = any(head(cumi.mean, width) > cumi.cut)
      }
      # 8/10/22 NEW VARIANT SCENARIO 2
      # IS SOUTHERN HEMESPHERE WINTER?
      # flag.sh.winter = (vdate.t+14) %>% month() %in% 6:8 # only do the first 2 months, but just in case monsoon season too
      # with a 2-wk delay
      # 8/18/22 - maybe don't need that delay, considering it might start in the late fall, rather than winter?
      flag.sh.winter = (vdate.t) %>% month() %in% 6:8 # only do the first 2 months, but just in case monsoon season too
      flag.nh.winter = ((vdate.t) + 15) %>% month() %in% c(12, 1:2)  # northern hemisphere winter, 11/15 - 2/15
      
      beta_tt = state0['beta',] * fn_getRelMob(fcast.mob[iwk,],state0['p.mob',]) 
      if(seasonality) {
        beta_tt = beta_tt * relR0[cur.wk,] 
      } 
      # add holiday eff
      if(as.Date(vdate.t) %in% seq(as.Date(paste0(format(as.Date(vdate.t), '%Y'),'/6/25')), length.out = 30, by = 'day') |
         as.Date(wk.train.end) %in% seq(as.Date(paste0(format(as.Date(vdate.t), '%Y'),'/7/1')), length.out = 30, by = 'day')
      )
        beta_tt = beta_tt * fn_get.p.beta.summer.holiday(vdate.t, wk.train.end, p.beta.summer.holiday, p.beta.after.holiday)
      
      # update severity
      severity.t['death',] = xpost['ifr',,wk.train.t]
      if(flag.higher.ifr){
        severity.t['death',] = fn_adj.parm4winter(parm.t = severity.t['death',], # parm to adjust
                                                  parm.name.t = 'ifr',
                                                  parm.low.cut = 1, # set it to a high level so it is invoked
                                                  parm.winter.amp = .25, parm.phase.shift = -30,
                                                  vdate.t = vdate.t)
      }
        
      # update detection rate
      state0['alpha',] = xpost['alpha',,wk.train.t]
      # higher det rate in the winter
      if(flag.higher.det.rate)
        state0['alpha',] = fn_adj.parm4winter(parm.t = state0['alpha',], # parm to adjust
                                              parm.name.t = 'alpha',
                                              parm.low.cut = .14, parm.winter.amp = .25, parm.phase.shift = 10,
                                              vdate.t = vdate.t)
      
      
      if((flag.recent.major.newV & iwk < 8) & (flag.mass.wave | flag.sh.winter))
        print(paste(vdate.t, 'newV not activated b/c recent major update'))
      
      print(paste('flag.mass.wave =',flag.mass.wave))
      
      if(((flag.mass.wave & flag.nh.winter) | # yes, there is a massive wave recently  & and it is northern hemmisphere (in season for the US) winter
          flag.sh.winter) & # yes, southern hemisphere winter
         ! (flag.recent.major.newV & iwk < 8) | # and the recent changes (w/in 2 months) haven't been update
         flag.newVrising.t # detected there is a new variant rising at time of the forecast, matched to this week
      ){ 
        # check susceptiblity level
        Smn.t = (state0['S1', ] %>% mean) / N
        
        
        # flag.mass.wave.hist
        # if it is triggered by a forecasted outbreak - need to be very careful, as it might be part of vicious cycle
        # should we restrict it to only based on real recorded recent major outbreak?
        frac.resmp.t = .5 # baseline resampling fraction to account for new variant
        dur.newVreplace.base.t = ifelse(as.Date(fcast.start) < as.Date('2021/5/1'), 13, # earlier days, not that much genomic surveillance
        ifelse(as.Date(fcast.start) < as.Date(date.omicron), 10+1, 8)) # during increased genomic surveillance, faster detection/reporting (reflected in replacement)
        # v.gr.base.t = ifelse(as.Date(fcast.start) < as.Date(date.omicron), 1.2, 1.4) # during increased genomic surveillance, faster detection/reporting (reflected in replacement)
        if(flag.mass.wave.fcast){  # it is triggered by a forecast major outbreak - no S incr
          s.incr.lwr = 0
          s.incr.upr = 0
          # print('no update on S b/c the model has projected a large wave')
        } else if(flag.newVrising.t){ # detected potential newV rising in recent week
          if(nrow(newv.pred)>1){
            # if more than 1, use the more aggressive one?
            newv.pred.t = newv.pred %>% dplyr::filter(wk2newVreplace == min(wk2newVreplace))
          } else {
            newv.pred.t = newv.pred
          }
          # print(newv.pred)
          # very fast replacement? let it scale with the speed?
          v.gr.t = newv.pred.t$v.growth.rate
          dur.newVreplace.t = newv.pred.t$dur.newVreplace
          QC.t = newv.pred.t$QC
          # if it started from a very very low case rate, 
          # est growth rate could be very high & may not reflect absolute changes
          # need some quality control
          if(mean(tail(obs_i[1:length(nwk.train),], 3)) / N * 100 < .08 | QC.t != ''){
            v.gr.t = pmin(v.gr.t, 1.5)
            dur.newVreplace.t = dur.newVreplace.t * 1.2
          }
          # conversely, if it starts at a very high case level
          # given the high background, reduce the timing a bit
          if(mean(tail(obs_i[1:length(nwk.train),], 3)) / N * 100 > .3){
            if(dur.newVreplace.t > 12)
              dur.newVreplace.t = dur.newVreplace.t * .8
          }
          frac.resmp.t = frac.resmp.t * # test28
            pmax(2/3, pmin(1.6, (v.gr.t / 1.2 * # scale with growth rate, assuming a baseline of 20% incr
                                   dur.newVreplace.base.t/dur.newVreplace.t)^.7)) # inversely scale with speed of replacement
          s.incr.lwr = 3/100 * N * 
            pmax(2/3, pmin(2.5, v.gr.t / 1.2 * # scale with growth rate, assuming a baseline of 20% incr
                             dur.newVreplace.base.t/dur.newVreplace.t)) # inversely scale with speed of replacement
          s.incr.upr = 9/100 * N * 
            pmax(2/3, pmin(2.5, v.gr.t / 1.2 * # scale with growth rate, assuming a baseline of 20% incr
                             dur.newVreplace.base.t/dur.newVreplace.t)) # inversely scale with speed of replacement
          # print(paste('incrS by', round(s.incr.lwr/N*100 * frac.resmp.t, 1),'-',round(s.incr.upr/N*100 * frac.resmp.t, 1),'%'))
        } else {
          
          s.incr.lwr = ifelse(flag.mass.wave.hist & flag.sh.winter, 3/100 * N, 2/100 * N) * ifelse((vdate.t) %>% month() %in% 8, .5, 1) # if both are true, higher increase
          # s.incr.upr = ifelse(flag.mass.wave.hist & flag.sh.winter, 6/100 * N, 4/100 * N) # smaller increase? - too low? test27, 28
          s.incr.upr = ifelse(flag.mass.wave.hist & flag.sh.winter, 9/100 * N, 6/100 * N) * ifelse((vdate.t) %>% month() %in% 8, .5, 1) # test 29
        }
        
        
        # is %S very low now? 
        # need some balance so that it won't increase S due to a vicious cycle with outbreaks reduce S and 
        # S keep getting increased out of nowhere continuously
        if(Smn.t < ifelse(flag.newVrising.t, .55, .4) & 
           # cum.dS.fcast < .2 & # also count how much S has been increased, if it's already a lot, don't do it
           cum.dS.fcast < max(cum.dS.fcast.cut, .015 * iwk) & # 9/14/22 allow it to scale with the number of week projected, assuming a 1.5% ave max per week
           cum.dimm.recentNfcast < cum.dimm.recentNfcast.cut # & relative to prior immunity
        ){
          # idx = sample(1:num_ens, size = round(num_ens * .5, 0)) # only do 50%
          # % resmp also based on whether there is a new variant
          idx = sample(1:num_ens, size = round(num_ens * frac.resmp.t, 0)) # only do 50%
          
          s0 = state0['S1',] %>% mean
          state0['S1', idx] = state0['S1', idx] + runif(length(idx), min = s.incr.lwr, max = s.incr.upr) 
          s1 = state0['S1',] %>% mean
          cum.dS.fcast = cum.dS.fcast + (s1 - s0) / N
          cum.dimm.recentNfcast = cum.dimm.recentNfcast + (s1 - s0) / pmax(N - states.train[nrow(states.train)]$S1, sum(states.train$newItot1))
          # increase by 2-6% for 20% of the ensemble
          # so, roughtly 4 * .5  = ~2% increase
          # print('incr S')
          # print(paste0('cum.dS.fcast = ', cum.dS.fcast))
          # print(paste0('cum.dimm.recentNfcast =',cum.dimm.recentNfcast))
        }
        
        # also widen beta a bit? 
        # only widen the larger half
        # idx = which(beta_tt > median(beta_tt) & beta_tt < quantile(beta_tt, .9))  # exclude the extreme - they are high already
        # beta_tt[idx] = 1.1 * (beta_tt[idx] - mean(beta_tt)) + mean(beta_tt)
        
        # 8/24/22 -for gauge the level of increase allowed to account for newV, do it before the seasonal component is added
        # 8/24/22 - allow larger increase if it is early stage where there is more room for increase
        # beta.base = xpost_mean %>% filter(Week.start < as.Date('2020/9/1')) %>% .$beta %>% mean
        # only do it if the major tx incr hasn't been accounted for
        if(! flag.recent.major.dtx){
          p.beta.newV = ifelse(mean(state0['beta',]) < 1.3 * beta0 & vdate.t < date.omicron, 1.3, 1.1)
          beta.newV.upr = ifelse(mean(state0['beta',]) < 1.3 * beta0, .95, .9)
          
          idx = which(beta_tt > mean(beta_tt) & beta_tt < quantile(beta_tt, beta.newV.upr))  # exclude the extreme - they are high already
          beta_tt[idx] = p.beta.newV * (beta_tt[idx] - mean(beta_tt)) + mean(beta_tt)
          
          # print(paste0('p.beta.newV = ', p.beta.newV, '; beta.newV.upr = ', beta.newV.upr))
        }
        
        # print(paste(vdate.t, 'newV activated'))
      }
      
      Sfcast = rbind(Sfcast, 
                     data.table(Week.start = vdate.t, Smn = round(mean(state0['S1',]) / N *100, 1)))
      
      # [2/23/21] check this if problem!
      dist_tm.to.detect = NULL
      for(ii in 1:num_ens){
        tmp = generation.time(dist_tm.to.detect.name,c(state0['Td.mean',ii],state0['Td.sd',ii]),truncate = tm.to.detect.max)
        dist_tm.to.detect=cbind(dist_tm.to.detect,tmp$GT[-1]); 
      }
      
      # [1/20/22] also update time from infection to death, b/c longer time-lag for omicron
      if((vdate.t >= (date.omicron+7) & vdate.t < (date.omicron+14)) |
         # first week of forecast and it hasn't been updated yet
         (iwk == 1 & vdate.t >= (date.omicron + 7))
      ){
        print('update tm to death for omicron')
        tm.to.death.max = tm.to.death.max_omicron
        if(tm.to.deathFromDiag){
          tm.from.inf.to.death.max = tm.to.diag + tm.to.death.max
        } else {
          tm.from.inf.to.death.max = tm.to.death.max
        }
        
        dist_tm.to.death = NULL  # time to death
        for(ii in 1:num_ens){
          tmp = generation.time(dist_tm.to.death.name,c(tm.to.outcome.mn['tm.to.death',ii] + tm2death_adj_omicron, 
                                                        tm.to.outcome.sd['tm.to.death',ii] # not good to increase sd
          ),
          truncate = tm.to.death.max # updated above
          )
          dist_tm.to.death=cbind(dist_tm.to.death,tmp$GT[-1]); 
        }
        if(tm.to.deathFromDiag){ # if the distribution of time to death is from diagnosis, not infectious
          # add time from infectious to diagnosis
          dist_tm.to.death = rbind(matrix(0,round(tm.to.diag,0),num_ens),dist_tm.to.death)
        }
      }
      
      if(epi.model == 'SEIRS'){
        
        simEpi=SEIRS(tm_strt, tm_end, tm_step=1, # 1 day time-step
                     tmstep = tmstep,
                     state0 = state0,
                     S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                     I0=state0[paste0('I',1:num_gr),], 
                     beta=beta_tt, 
                     Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                     seed=seed.t, stoch=stoch,
                     severity = severity.t,
                     newI.previous = newI.previous,
                     dist_tm.to.detect = dist_tm.to.detect,
                     dist_tm.to.death = dist_tm.to.death)
        
      } else if(epi.model == 'SEIRSV'){
        daVacc.t = fda.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
        
        if(nrow(daVacc.t)<1){  # no data yet
          V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
        } else { # yes data
          
          daVacc.t$date = daVacc.t$date %>% as.Date
          
          # make sure it includes a full week
          dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
          daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
          daVacc.t[is.na(daVacc.t)] = 0
          V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
          V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
          
          # print('start vacc!')
          
        }
        simEpi=SEIRSV(tm_strt, tm_end, tm_step=1, # 1 day time-step
                      tmstep = tmstep,
                      state0 = state0,
                      S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                      I0=state0[paste0('I',1:num_gr),], 
                      beta=beta_tt, 
                      Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                      seed=seed.t, stoch=stoch, 
                      severity = severity.t,
                      newI.previous = newI.previous,
                      dist_tm.to.detect = dist_tm.to.detect,
                      dist_tm.to.death = dist_tm.to.death,
                      percSmax.t = percSmax.t,
                      V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                      # these are total number of vaccinees with unknown immunity
                      # but pre-ajust for time lag from vaccination to immune protection
                      VE1 = VE1, VE2=VE2 # Vaccine efficacy, need further adjustment by prior immunity 
        ) # include the delay reporting etc.
      } else if (epi.model == 'SEIRSVimmLoss'){
        
        daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
        
        if(nrow(daVacc.t)<1){  # no data yet
          V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
        } else { # yes data
          
          daVacc.t$date = daVacc.t$date %>% as.Date
          
          # make sure it includes a full week
          dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
          daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
          daVacc.t[is.na(daVacc.t)] = 0
          V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
          V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
          
          # print('start vacc!')
          
        }
        
        # get vaccine-induced immunity related parameter
        parmVimmLoss.t = get_parmVimmLoss(date.t = as.Date(vdate.t), 
                                          date.delta = date.delta, #  as.Date('2021/07/01'), 
                                          date.omicron = date.omicron,
                                          # wt, alpha
                                          tm.imm.wt = tm.imm.wt, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.wt = p.imm.wane.max.wt, # maximal level of immunity loss (=1 or lower)
                                          k.wt = k.wt, # logistic function tuning parameter
                                          
                                          tm.imm.delta = tm.imm.delta, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.delta = p.imm.wane.max.delta, # maximal level of immunity loss (=1 or lower)
                                          k.delta = k.delta, # logistic function tuning parameter
                                          
                                          tm.imm.omicron = tm.imm.omicron, # during of vaccine-incuded protection against infection
                                          p.imm.wane.max.omicron = p.imm.wane.max.omicron, # maximal level of immunity loss (=1 or lower)
                                          k.omicron = k.omicron # logistic function tuning parameter
        )
        simEpi=SEIRSVimmLoss(tm_strt, tm_end, tm_step=1, # 1 day time-step
                             tmstep = tmstep,
                             state0 = state0,
                             S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                             I0=state0[paste0('I',1:num_gr),], 
                             beta=beta_tt, 
                             Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                             seed=seed.t, stoch=stoch, 
                             severity = severity.t,
                             newI.previous = newI.previous,
                             dist_tm.to.detect = dist_tm.to.detect,
                             dist_tm.to.death = dist_tm.to.death,
                             percSmax.t = percSmax.t,
                             V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                             # these are total number of vaccinees with unknown immunity
                             # but pre-ajust for time lag from vaccination to immune protection
                             VE1 = VE1, VE2=VE2, # Vaccine efficacy, need further adjustment by prior immunity 
                             newVacc = newVacc.t, # prior vaccinated & protected from infection
                             parmVimmLoss = parmVimmLoss.t, # parameters related to the immunity period of vaccine, against infection
                             R0 = Rimm.t   # to check those recovered from infection and remaining immune to infection
        ) # include the delay reporting etc.
        
        # update hist newVacc.t
        newVacc.t = simEpi$newVacc
        Rimm.t  = simEpi$R
        simEpi$newVacc = NULL
        simEpi$cumVimmLoss = NULL
        simEpi$R = NULL
        simEpi$cumIimmLoss = NULL
      }
      
      # re-assemble to the same order as the prior: state0
      n.end = tm_end - tm_strt + 2
      state.new = NULL
      for(i in 1:(length(simEpi)-1)){
        tmp = simEpi[[i]][n.end,,drop=F]; 
        rownames(tmp)=gsub('cumI','newI',paste0(names(simEpi)[i],1:num_gr))
        state.new = rbind(state.new,tmp)
      }
      
      # state.new = rbind(state.new, state0[parm.names,])
      state.new = rbind(state.new, fcast.state0[parm.names,]) # restore the parameters, so the deflation doesn't accumulate
      state.new = state.new[rownames(state0),] # make sure the order is the same
      
      state0 = state.new # update
      newi = simEpi$cumItot
      newi = newi[-1,] - newi[-nrow(newi),]
      newI.previous = rbind(newI.previous,newi) # include the most recent week
      fcast.inf[iwk,] = state.new[idx.newItot,]
      fcast.case[iwk,] = state.new[idx.obs_i,]
      fcast.death[iwk,] = state.new[idx.obs_d,]
      
      if(any(state0 < 0))
        break
      
    } # end week
    
    # get summary stats
    prob_vec = c(.5, .25, .75, .025, .975, .05, .95, .1, .9)
    prob.name_vec = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr','ci90.lwr','ci90.upr','ci80.lwr','ci80.upr')
    tmp.inf = fcast.inf %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.inf) = prob.name_vec
    tmp.case = fcast.case %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.case) = prob.name_vec
    tmp.death = fcast.death %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.death) = prob.name_vec
    
    fcast_stats = rbind(data.table(measure = 'Infections', Week.start=fcast.wk.starts,tmp.inf),
                        data.table(measure = 'Cases', Week.start=fcast.wk.starts,tmp.case),
                        data.table(measure = 'Deaths', Week.start=fcast.wk.starts,tmp.death)
    )
    
    # cumulative
    # get summary stats
    tmp.inf = fcast.inf %>% apply(2,cumsum) %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.inf) = prob.name_vec
    tmp.case = fcast.case %>% apply(2,cumsum) %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.case) = prob.name_vec
    tmp.death = fcast.death %>% apply(2,cumsum) %>% apply(1, quantile, prob = prob_vec) %>% t
    colnames(tmp.death) = prob.name_vec
    fcast_stats = rbind(fcast_stats, data.table(measure = 'Cumulative Infections', Week.start=fcast.wk.starts,tmp.inf),
                        data.table(measure = 'Cumulative Cases', Week.start=fcast.wk.starts,tmp.case),
                        data.table(measure = 'Cumulative Deaths', Week.start=fcast.wk.starts,tmp.death)
    )
    
    # get probability dist for computing log score
    fcastDist = fn_getProbDist(fcast.case, fcast.death, bins.case, bins.death)
    
    fcast_stats_newV = fcast_stats
    fcastDist_newV = fcastDist
    fcast_stats_newV$scenario = 'newV'
    fcastDist_newV$scenario = 'newV'
  } # end projection
  
  
  fcast_stats = rbind(fcast_stats_base, fcast_stats_newV)
  fcastDist = rbind(fcastDist_base, fcastDist_newV)
  rm(fcast_stats_base, fcast_stats_newV, fcastDist_base, fcastDist_newV)
  
  return(list(Rt_stats = Rt_stats, R0_stats = R0_stats, 
              Rtx_stats = Rtx_stats, 
              Rtx_ens = Rtx_ens,
              states_stats = states_stats, 
              xpost_mean = xpost_mean, xpost_sd = xpost_sd, 
              xprior_mean = xprior_mean, xprior_sd = xprior_sd,
              cumIperc_stats = cumIperc_stats,
              cumIperc_ens = cumIperc,
              Susceptibility_stats = Susceptibility_stats,
              Susceptibility_ens = Susceptibility,
              xpost.last=xpost.last,
              immLoss_stats = immLoss_stats,
              fcast_stats = fcast_stats, 
              fcastDist = fcastDist
              ))
  
}

