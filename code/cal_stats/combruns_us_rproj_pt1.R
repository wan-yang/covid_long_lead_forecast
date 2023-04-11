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
      res.eval = NULL; 
      try(load(paste0(dir_res, 'res.eval', paste0('_df',deflat.t), '_', loc.t,date.tag,'.RData')))
      # try(load(paste0(dir_res, 'res.eval', paste0('_df',deflat.t), '_', loc.t,'.RData')))
      
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


# 10/4/22 -
# the long-term n-wk ahead forecast using non.Omicron data extended to the Omicron period
# these should be excluded for evaluation
# criteria: truncated at end of Dec 2021
# can only do this for the n-wk ahead forecast and not the overall targets

res.eval$fcast.start.month = NULL
res.eval$seasonality = factor(res.eval$seasonality, levels = c('No seasonality', 'Seasonality assumed', 'Transformed seasonality'),
                              labels = sn_vec)
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
{
  # COMPARE baseline (no deflation, no sn, no new variant) vs. opt (.9 deflat, trans sn, new v)
  cut0 = .1 # count those with abs rel error within 10% of the observed as accurate
  cut1 = .25 # count those with abs rel error within 25% 10% of the observed as accurate
  cut2 = .5 # count those with abs rel error within 50% 10% of the observed as accurate
  acc.base.v.opt = rbind(res.eval %>% filter(fcast.deflat == 1 & seasonality == 'no' & scenario == 'asIs') %>% mutate(tag = 'base'),
                         res.eval %>% filter(fcast.deflat == .9 & seasonality == 'y/tf' & scenario == 'newV') %>% mutate(tag = 'opt')) %>%
    melt(., id.vars =  c('tag', 'loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
    dcast(., loc + variant + measure + metric + fcast.start.week + target + variable ~ tag) %>% 
    filter(variable == 'mean') %>%
    mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                            labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
    mutate(diff = `opt` - `base`) %>% data.table()
  
  acc.base.v.opt$week.t = acc.base.v.opt$fcast.start.week + as.numeric(gsub('wk ahead', '', acc.base.v.opt$target)) * 7
  wave.t = acc.base.v.opt[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
  acc.base.v.opt = acc.base.v.opt %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
  # add if it is respiratory season
  acc.base.v.opt = acc.base.v.opt %>% 
    mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                               target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                               !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
    mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))
  
  setDT(acc.base.v.opt)
  
  stats.acc.base.v.opt = fn_getAllStats(tda.t = acc.base.v.opt, v.cp.t = 'opt', v.ref.t = 'base')
  save(acc.base.v.opt, file = paste0(dir_res, 'acc.base.v.opt_', stat.tag, '.RData'))
  save(stats.acc.base.v.opt, file = paste0(dir_res, 'stats.acc.base.v.opt_', stat.tag, '.RData'))
  
}

time.t = Sys.time()
print(paste('time since start: ', time.t - t1))

rm(res.eval)

######################################################
# look at scores 
######################################################
if(T){
  res.score00 = NULL
  deflat_vec = c(.9, .95, 1)
  for(deflat.t in deflat_vec){
    for(loc.t in locs){
      rm(res.score)
      res.score = NULL
      
      try(load(paste0(dir_res, 'res.score', paste0('_df',deflat.t), '_', loc.t,date.tag,'.RData')))
      
      res.score00 = rbind(res.score00, res.score)
    }
  }
  
  res.score = res.score00
  rm(res.score00)
  
  deflat_vec = res.score %>% .$fcast.deflat %>% unique %>% sort(., decreasing = T)
  sce_vec =  res.score %>% .$scenario %>% unique
  sn_vec  =  res.score %>% .$seasonality %>% unique %>% as.character()
  fcast.start.month_vec = res.score$fcast.start.month %>% unique %>% sort
  fcast.start.week_vec = res.score$fcast.start.week %>% unique %>% sort
}

sn_vec = c('no','y/fix','y/tf')

res.score$fcast.start.month = NULL
res.score$fcast.start.week = res.score$fcast.start.week %>% as.Date
res.score$target = factor(res.score$target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'))
# res.score$target2 = factor(res.score$target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
#                           labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))
res.score$fcast.deflat = factor(res.score$fcast.deflat)
res.score$seasonality = factor(res.score$seasonality, levels = c('No seasonality', 'Seasonality assumed', 'Transformed seasonality'),
                               labels = sn_vec)

# 10/4/22 -
# the long-term n-wk ahead forecast using non.Omicron data extended to the Omicron period
# these should be excluded for evaluation
# criteria: truncated at end of Dec 2021
# can only do this for the n-wk ahead forecast and not the overall targets
res.score = res.score %>% # filter(loc == 'California' & seasonality == 'no' & fcast.deflat==.9 & measure == 'case' & metric == 'logScore.loose' & scenario == 'asIs') %>% 
  # rowwise() %>% 
  mutate(week.t = fcast.start.week + as.numeric(gsub('wk ahead', '', target)) * 7) %>%
  filter(!(variant == 'non.Omicron' & (!is.na(week.t) & as.Date(week.t) > as.Date('2021/12/31')))) %>% 
  mutate(wave.t = NULL)


t1 = Sys.time()

## COMPARISON
# COMPARE baseline (no deflation, no sn, no new variant) vs. opt (.9 deflat, trans sn, new v)
tda.base.v.opt = rbind(res.score %>% filter(fcast.deflat == 1 & seasonality == 'no' & scenario == 'asIs') %>% mutate(tag = 'base'),
                       res.score %>% filter(fcast.deflat == .9 & seasonality == 'y/tf' & scenario == 'newV') %>% mutate(tag = 'opt')) %>%
  melt(., id.vars =  c('tag', 'loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + measure + metric + fcast.start.week + target + variable ~ tag) %>% 
  filter(variable == 'mean') %>%
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `opt` - `base`) %>% data.table()
tda.base.v.opt$week.t = tda.base.v.opt$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.base.v.opt$target)) * 7
wave.t = tda.base.v.opt[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.base.v.opt = tda.base.v.opt %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
# add if it is respiratory season
tda.base.v.opt = tda.base.v.opt %>% 
  mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                             target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                             !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))

setDT(tda.base.v.opt)

stats.base.v.opt = fn_getAllStats(tda.t = tda.base.v.opt, v.cp.t = 'opt', v.ref.t = 'base')
save(stats.base.v.opt, file = paste0(dir_res, 'stats.base.v.opt_', stat.tag, '.RData'))
save(tda.base.v.opt, file = paste0(dir_res, 'tda.base.v.opt_', stat.tag, '.RData'))

# COMPARE DEFLATE
tda.1.v.095 = res.score %>% filter(fcast.deflat %in% c(1,.95)) %>%
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + seasonality + measure + metric + fcast.start.week + target + scenario + variable ~ fcast.deflat) %>% 
  filter(variable == 'mean') %>%
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `0.95` - `1`) %>% 
  rename('opt' = `0.95`,'base' = `1`) %>% mutate(cp = 'deflat1.v.095') %>% data.table()
tda.1.v.095$week.t = tda.1.v.095$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.1.v.095$target)) * 7
wave.t = tda.1.v.095[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.1.v.095 = tda.1.v.095 %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()

tda.1.v.09 = res.score %>% filter(fcast.deflat %in% c(1,.9)) %>%
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + seasonality + measure + metric + fcast.start.week + target + scenario + variable ~ fcast.deflat) %>% 
  filter(variable == 'mean') %>%
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `0.9` - `1`) %>% 
  rename('opt' = `0.9`,'base' = `1`) %>% mutate(cp = 'deflat1.v.09') %>% data.table()
tda.1.v.09$week.t = tda.1.v.09$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.1.v.09$target)) * 7
wave.t = tda.1.v.09[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.1.v.09 = tda.1.v.09 %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()

tda.095.v.09 = res.score %>% filter(fcast.deflat %in% c(.95,.9)) %>%
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + seasonality + measure + metric + fcast.start.week + target + scenario + variable ~ fcast.deflat) %>% 
  filter(variable == 'mean') %>%
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `0.9` - `0.95`) %>% 
  rename('opt' = `0.9`,'base' = `0.95`) %>% mutate(cp = 'deflat095.v.09') %>% data.table()
tda.095.v.09$week.t = tda.095.v.09$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.095.v.09$target)) * 7
wave.t = tda.095.v.09[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.095.v.09 = tda.095.v.09 %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()

tda.cp.deflat = rbind(tda.1.v.095, tda.1.v.09, tda.095.v.09)
rm(tda.1.v.095, tda.1.v.09, tda.095.v.09)

tda.cp.deflat = tda.cp.deflat %>% 
  mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                             target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                             !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))


stats = NULL; da.t = tda.cp.deflat
cp.deflat_vec = paste0('deflat', c('1.v.095','1.v.09','095.v.09'))
cp.deflat_vec.lab = c('none vs .95', 'none vs .9', '.95 vs .9')
for(cp.t in cp.deflat_vec){
  for(sn.t in sn_vec){
    for(sce.t in sce_vec){
      print(paste(cp.t, sn.t, sce.t))
      rm(stats.t, ttda.t)
      ttda.t = da.t %>% filter(cp == cp.t & seasonality == sn.t & scenario == sce.t)
      stats.t = fn_getAllStats(tda.t = ttda.t, v.cp.t = 'opt', v.ref.t = 'base') %>%
        mutate(cp = cp.t, seasonality = sn.t, scenario = sce.t)
      
      stats = rbind(stats, stats.t)
      
      time.t = Sys.time()
      print(paste('time since start: ', time.t - t1))
      
    }
  }
}
rm(da.t, stats.t, ttda.t)
stats.cp.deflat = stats
save(stats.cp.deflat, file = paste0(dir_res, 'stats.cp.deflat_', stat.tag, '.RData'))
save(tda.cp.deflat, file = paste0(dir_res, 'tda.cp.deflat_', stat.tag, '.RData'))

# COMPARE SCENARIOS
tda.asIs.v.newV = res.score %>% 
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + fcast.deflat + seasonality + measure + metric + fcast.start.week + target + variable ~ scenario) %>% 
  filter(variable == 'mean') %>%
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `newV` - `asIs`) %>% 
  rename('opt' = `newV`,'base' = `asIs`) %>% mutate(cp = 'asIs.v.newV') %>% data.table()
tda.asIs.v.newV$week.t = tda.asIs.v.newV$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.asIs.v.newV$target)) * 7
wave.t = tda.asIs.v.newV[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.asIs.v.newV = tda.asIs.v.newV %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()
tda.cp.sce = tda.asIs.v.newV; 
rm(tda.asIs.v.newV)

# add respiratory season indicator
tda.cp.sce = tda.cp.sce %>% 
  mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                             target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                             !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))

cp.sce_vec = c('asIs.v.newV')
cp.sce_vec.lab = 'none vs new variants'
stats = NULL; da.t = tda.cp.sce
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
stats.cp.sce = stats
save(stats.cp.sce, file = paste0(dir_res, 'stats.cp.sce_', stat.tag, '.RData'))
save(tda.cp.sce, file = paste0(dir_res, 'tda.cp.sce_', stat.tag, '.RData'))

# COMPARE SEASONALITY
print('start summarizing scores...')
tda.no.v.yfix = res.score %>% filter(seasonality %in% c('no','y/fix')) %>%
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + fcast.deflat + measure + metric + fcast.start.week + target + scenario + variable ~ seasonality) %>% 
  filter(variable == 'mean') %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `y/fix` - `no`) %>% 
  rename('opt' = `y/fix`,'base' = `no`) %>% mutate(cp = 'no.v.yfix') %>% data.table()
tda.no.v.yfix$week.t = tda.no.v.yfix$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.no.v.yfix$target)) * 7
mainv.t= tda.no.v.yfix[, list(main.variant = fn_getVariant(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom)), by = c('loc', 'week.t')]
tda.no.v.yfix = tda.no.v.yfix %>% left_join(x = ., y = mainv.t, by = c('loc', 'week.t'))  %>% data.table()
tda.no.v.yfix$vlab = factor(tda.no.v.yfix$main.variant, 
                            levels = c('Alpha', 'Delta','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2"), 
                            labels = c('a', 'd', 'o', 'o', 'o', 'o'))
wave.t = tda.no.v.yfix[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.no.v.yfix = tda.no.v.yfix %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()

tda.no.v.ytf = res.score %>% filter(seasonality %in% c('no','y/tf')) %>%
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + fcast.deflat + measure + metric + fcast.start.week + target + scenario + variable ~ seasonality) %>%
  filter(variable == 'mean') %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `y/tf` - `no`)  %>% 
  rename('opt' = `y/tf`,'base' = `no`) %>% mutate(cp = 'no.v.ytf') %>% data.table()
tda.no.v.ytf$week.t = tda.no.v.ytf$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.no.v.ytf$target)) * 7
mainv.t= tda.no.v.ytf[, list(main.variant = fn_getVariant(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom)), by = c('loc', 'week.t')]
tda.no.v.ytf = tda.no.v.ytf %>% left_join(x = ., y = mainv.t, by = c('loc', 'week.t')) %>% data.table()
tda.no.v.ytf$vlab = factor(tda.no.v.ytf$main.variant, 
                           levels = c('Alpha', 'Delta','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2"), 
                           labels = c('a', 'd', 'o', 'o', 'o', 'o'))
wave.t = tda.no.v.ytf[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.no.v.ytf = tda.no.v.ytf %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()

tda.yfix.v.ytf = res.score %>% filter(seasonality %in% c('y/fix','y/tf')) %>%
  melt(., id.vars =  c('loc', 'variant', 'fcast.type', 'fcast.deflat', 'fcast.start.week','measure', 'metric', 'target', 'scenario', 'seasonality')) %>%
  dcast(., loc + variant + fcast.type + fcast.deflat + measure + metric + fcast.start.week + target + scenario + variable ~ seasonality) %>%
  filter(variable == 'mean') %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(diff = `y/tf` - `y/fix`) %>% 
  rename('opt' = `y/tf`,'base' = `y/fix`) %>% mutate(cp = 'yfix.v.ytf') %>% data.table()
tda.yfix.v.ytf$week.t = tda.yfix.v.ytf$fcast.start.week + as.numeric(gsub('wk ahead', '', tda.yfix.v.ytf$target)) * 7
mainv.t= tda.yfix.v.ytf[, list(main.variant = fn_getVariant(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom)), by = c('loc', 'week.t')]
tda.yfix.v.ytf = tda.yfix.v.ytf %>% left_join(x = ., y = mainv.t, by = c('loc', 'week.t')) %>% data.table()
tda.yfix.v.ytf$vlab = factor(tda.yfix.v.ytf$main.variant, 
                             levels = c('Alpha', 'Delta','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2"), 
                             labels = c('a', 'd', 'o', 'o', 'o', 'o'))
wave.t = tda.yfix.v.ytf[, list(wave = fn_getWave(week.t = week.t, loc.t = loc, vdates.dom = vdates.dom, variant = variant)), by = c('loc', 'variant', 'week.t')]
tda.yfix.v.ytf = tda.yfix.v.ytf %>% left_join(x = ., y = wave.t, by = c('loc', 'variant', 'week.t')) %>% data.table()

tda.cp.sn = rbind(tda.no.v.yfix, tda.no.v.ytf, tda.yfix.v.ytf)
rm(tda.no.v.yfix, tda.no.v.ytf, tda.yfix.v.ytf)
# add if it is respiratory season
tda.cp.sn = tda.cp.sn %>% 
  mutate(resp.sn = case_when(as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Respiratory season',
                             target %in% c('peak week','peak intensity','total') & (as.numeric(format(as.Date(fcast.start.week)+15, '%m')) %in% resp.months4sntargets) ~ 'Respiratory season',
                             !as.numeric(format(as.Date(week.t)+15, '%m')) %in% resp.months ~ 'Off season')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('Respiratory season','Off season')))


cp.sn_vec = c('no.v.yfix','no.v.ytf','yfix.v.ytf')
cp.sn_vec.lab = c('none vs fixed', 'none vs trans', 'fixed vs trans')
stats = NULL; da.t = tda.cp.sn
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
stats.cp.sn = stats
save(stats.cp.sn, file = paste0(dir_res, 'stats.cp.sn_', stat.tag, '.RData'))
save(tda.cp.sn, file = paste0(dir_res, 'tda.cp.sn_', stat.tag, '.RData'))

time.t = Sys.time()
print(paste('time since score start: ', time.t - t1))

