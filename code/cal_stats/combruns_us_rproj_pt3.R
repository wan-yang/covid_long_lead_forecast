# compare log score / point prediction accuracy
# 10/4/22 -
# the long-term n-wk ahead forecast using non.Omicron data extended to the Omicron period
# these should be excluded for evaluation
# criteria: truncated at end of Dec 2021

dir_data = './data/'
dir_code = './code/'
dir_res = paste0('./results/')
dir_out = paste0('./outputs/')

source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code,'fn_util.R'))


num_runs = 10
stoch = T
date.tag = ''
data.tag = 'nyt' # whether new york times data or hopkins data are use
tag.model = ''

vtags = c('non.Omicron', 'Omicron')

locs = c('California', 'Florida', 'Iowa', 'Massachusetts', 'Michigan', 
         'New York', 'Pennsylvania', 'Texas', 'Washington', 'Wyoming')

# read the datas key variants became predominant in different states
vdates = read.csv(paste0(dir_data, 'dates_variant_by_state.csv')) %>% data.table()
vdates.dom = read.csv(paste0(dir_data, 'tm.variant_by_state_smoothed.wkly.csv'))  %>% data.table()
vdates.dom$date.start = vdates.dom$date.start %>% as.Date()
vdates.dom$date.end = vdates.dom$date.end %>% as.Date()

if(T){
  res.eval00 = NULL
  deflat_vec = c(.9, .95, 1)
  for(deflat.t in deflat_vec){
    for(loc.t in locs){
      rm(res.eval)
      res.eval = NULL
      
      try(load(paste0(dir_res, 'res.eval', paste0('_df',deflat.t), '_', loc.t,date.tag,'.RData')))
      
      res.eval00 = rbind(res.eval00, res.eval)
      
    }
  }
  
  res.eval = res.eval00
  rm(res.eval00)
  
  deflat_vec = res.eval %>% .$fcast.deflat %>% unique %>% sort(., decreasing = T)
  sce_vec =  res.eval %>% .$scenario %>% unique
  sn_vec  =  res.eval %>% .$seasonality %>% unique %>% as.character()
  fcast.start.month_vec = res.eval$fcast.start.month %>% unique %>% sort
  fcast.start.week_vec = res.eval$fcast.start.week %>% unique %>% sort
}

sn_vec = c('no','y/fix','y/tf')

res.eval$fcast.start.month = NULL
res.eval$seasonality = factor(res.eval$seasonality, levels = c('No seasonality', 'Seasonality assumed', 'Transformed seasonality'),
                              labels = c('no','y/fix','y/tf'))
nrow(res.eval)
res.eval = res.eval %>% filter(measure %in% c('Cases','Deaths')) %>% # rowwise() %>% 
  mutate(week.t = fcast.start.week + as.numeric(gsub('wk ahead', '', target)) * 7) %>%
  filter(!(variant == 'non.Omicron' & (!is.na(week.t) & as.Date(week.t) > as.Date('2021/12/31')))) %>% 
  mutate(wave.t = NULL)
# filter(!(variant == 'non.Omicron' & as.Date(week.t) > as.Date('2021/12/31'))) # this will exclude all those week.t is na - the overall targets

# add wether it is respiratory virus season
resp.months = c(11:12, 1:4) # shift by half a month, so roughly 10/15 to 4/15
# if it is the long-term targets (peak week etc), count the ones covering at least 3 months of the resp.sn (as the forecast was for the next 6 months)
# that would be 7/15/ - 8/14, 8/15 - 9/14, 9/15-10/14, 10/15 - 11/14, 11/15 - 12/14, and 12/15 - 1/14
resp.months4sntargets = c(8:12, 1)

stat.tag = 'wilcoxon'


######################################################
# look at pt estimates
######################################################
t1 = Sys.time()
# COMPARE SCENARIOS
{
  acc.asIs.v.newV = res.eval %>% 
    melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
    dcast(., loc + variant + fcast.type + fcast.deflat + seasonality + measure + metric + fcast.start.week + target + variable ~ scenario) %>% 
    filter(variable == 'mean') %>%
    mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                            labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
    mutate(diff = `newV` - `asIs`) %>% 
    rename('opt' = `newV`,'base' = `asIs`) %>% mutate(cp = 'asIs.v.newV') %>% data.table()
  acc.asIs.v.newV$week.t = acc.asIs.v.newV$fcast.start.week + as.numeric(gsub('wk ahead', '', acc.asIs.v.newV$target)) * 7
  wave.t = acc.asIs.v.newV[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
  acc.asIs.v.newV = acc.asIs.v.newV %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
  acc.cp.sce = acc.asIs.v.newV; 
  rm(acc.asIs.v.newV)
  
  # add respiratory season indicator
  acc.cp.sce = acc.cp.sce %>% 
    mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                               target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                               !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
    mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))
  
  cp.sce_vec = c('asIs.v.newV')
  cp.sce_vec.lab = 'none vs new variants'
  stats = NULL; da.t = acc.cp.sce
  for(cp.t in cp.sce_vec){
    for(sn.t in rev(sn_vec)){
      for(deflat.t in rev(deflat_vec)){
        print(paste(cp.t, sn.t, deflat.t, collapse = ','))
        rm(stats.t, ttda.t)
        ttda.t = da.t %>% filter(cp == cp.t & seasonality == sn.t & fcast.deflat == deflat.t)
        stats.t = fn_getAllStats(tda.t = ttda.t, v.cp.t = 'opt', v.ref.t = 'base') %>%
          mutate(cp = cp.t, seasonality = sn.t, fcast.deflat = deflat.t)
        
        stats = rbind(stats, stats.t)
        
        time.t = Sys.time()
        print(paste('time since start: ', time.t - t1))
      }
    }
  }
  rm(da.t, stats.t, ttda.t)
  stats.acc.cp.sce = stats
  save(stats.acc.cp.sce, file = paste0(dir_res, 'stats.acc.cp.sce_', stat.tag, '.RData'))
  save(acc.cp.sce, file = paste0(dir_res, 'acc.cp.sce_', stat.tag, '.RData'))
  rm(stats.acc.cp.sce, acc.cp.sce)
  
}






