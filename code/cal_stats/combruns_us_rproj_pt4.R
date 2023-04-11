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
# COMPARE SEASONALITY
t1 = Sys.time()
cp.deflat_vec = paste0('deflat', c('1.v.095','1.v.09','095.v.09'))
cp.deflat_vec.lab = c('none vs .95', 'none vs .9', '.95 vs .9')
cp.sn_vec = c('no.v.yfix','no.v.ytf','yfix.v.ytf')
cp.sn_vec.lab = c('none vs fixed', 'none vs trans', 'fixed vs trans')

{
  acc.no.v.yfix = res.eval %>% filter(seasonality %in% c('no','y/fix')) %>%
    melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
    dcast(., loc + variant + fcast.type + fcast.deflat + measure + metric + fcast.start.week + target + scenario + variable ~ seasonality) %>% 
    filter(variable == 'mean') %>% 
    mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                            labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
    mutate(diff = `y/fix` - `no`) %>% 
    rename('opt' = `y/fix`,'base' = `no`) %>% mutate(cp = 'no.v.yfix') %>% data.table()
  acc.no.v.yfix$week.t = acc.no.v.yfix$fcast.start.week + as.numeric(gsub('wk ahead', '', acc.no.v.yfix$target)) * 7
  mainv.t= acc.no.v.yfix[, list(main.variant = fn_getVariant(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom)), by = c('loc', 'week.t')]
  acc.no.v.yfix = acc.no.v.yfix %>% left_join(x = ., y = mainv.t, by = c('loc', 'week.t'))  %>% data.table()
  acc.no.v.yfix$vlab = factor(acc.no.v.yfix$main.variant, 
                              levels = c('Alpha', 'Delta','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2"), 
                              labels = c('a', 'd', 'o', 'o', 'o', 'o'))
  wave.t = acc.no.v.yfix[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
  acc.no.v.yfix = acc.no.v.yfix %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
  
  acc.no.v.ytf = res.eval %>% filter(seasonality %in% c('no','y/tf')) %>%
    melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
    dcast(., loc + variant + fcast.type + fcast.deflat + measure + metric + fcast.start.week + target + scenario + variable ~ seasonality) %>%
    filter(variable == 'mean') %>% 
    mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                            labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
    mutate(diff = `y/tf` - `no`)  %>% 
    rename('opt' = `y/tf`,'base' = `no`) %>% mutate(cp = 'no.v.ytf') %>% data.table()
  acc.no.v.ytf$week.t = acc.no.v.ytf$fcast.start.week + as.numeric(gsub('wk ahead', '', acc.no.v.ytf$target)) * 7
  mainv.t= acc.no.v.ytf[, list(main.variant = fn_getVariant(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom)), by = c('loc', 'week.t')]
  acc.no.v.ytf = acc.no.v.ytf %>% left_join(x = ., y = mainv.t, by = c('loc', 'week.t')) %>% data.table()
  acc.no.v.ytf$vlab = factor(acc.no.v.ytf$main.variant, 
                             levels = c('Alpha', 'Delta','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2"), 
                             labels = c('a', 'd', 'o', 'o', 'o', 'o'))
  wave.t = acc.no.v.ytf[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
  acc.no.v.ytf = acc.no.v.ytf %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
  
  acc.yfix.v.ytf = res.eval %>% filter(seasonality %in% c('y/fix','y/tf')) %>%
    melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
    dcast(., loc + variant + fcast.type + fcast.deflat + measure + metric + fcast.start.week + target + scenario + variable ~ seasonality) %>%
    filter(variable == 'mean') %>% 
    mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                            labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
    mutate(diff = `y/tf` - `y/fix`) %>% 
    rename('opt' = `y/tf`,'base' = `y/fix`) %>% mutate(cp = 'yfix.v.ytf') %>% data.table()
  acc.yfix.v.ytf$week.t = acc.yfix.v.ytf$fcast.start.week + as.numeric(gsub('wk ahead', '', acc.yfix.v.ytf$target)) * 7
  mainv.t= acc.yfix.v.ytf[, list(main.variant = fn_getVariant(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom)), by = c('loc', 'week.t')]
  acc.yfix.v.ytf = acc.yfix.v.ytf %>% left_join(x = ., y = mainv.t, by = c('loc', 'week.t')) %>% data.table()
  acc.yfix.v.ytf$vlab = factor(acc.yfix.v.ytf$main.variant, 
                               levels = c('Alpha', 'Delta','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2"), 
                               labels = c('a', 'd', 'o', 'o', 'o', 'o'))
  wave.t = acc.yfix.v.ytf[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
  acc.yfix.v.ytf = acc.yfix.v.ytf %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
  
  acc.cp.sn = rbind(acc.no.v.yfix, acc.no.v.ytf, acc.yfix.v.ytf)
  rm(acc.no.v.yfix, acc.no.v.ytf, acc.yfix.v.ytf)
  # add if it is respiratory season
  acc.cp.sn = acc.cp.sn %>% 
    mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                               target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                               !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
    mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))
  

  stats = NULL; da.t = acc.cp.sn
  for(cp.t in cp.sn_vec){
    for(sce.t in sce_vec){
      for(deflat.t in deflat_vec){
        
        print(paste(cp.t, sce.t, deflat.t, collapse = ','))
        rm(stats.t, ttda.t)
        ttda.t = da.t %>% filter(cp == cp.t & scenario == sce.t & fcast.deflat == deflat.t)
        stats.t = fn_getAllStats(tda.t = ttda.t, v.cp.t = 'opt', v.ref.t = 'base') %>%
          mutate(cp = cp.t, scenario = sce.t, fcast.deflat = deflat.t)
        
        stats = rbind(stats, stats.t)
        
        time.t = Sys.time()
        print(paste('time since start: ', time.t - t1))
      }
    }
  }
  rm(da.t, stats.t, ttda.t)
  stats.acc.cp.sn = stats
  save(stats.acc.cp.sn, file = paste0(dir_res, 'stats.acc.cp.sn_', stat.tag, '.RData'))
  save(acc.cp.sn, file = paste0(dir_res, 'acc.cp.sn_', stat.tag, '.RData'))
  rm(stats.acc.cp.sn, acc.cp.sn)
}