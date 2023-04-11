# to read and combine the arima runs (comparison model for the US long lead forecast study)
# 3/29/23

dir_data = './data/'
dir_code = './code/'

dir_res = './results/rproj_arima/'
dir_out = './outputs/'

if(! file.exists(dir_res))  dir.create(dir_res,recursive = T)


num_runs = 10
stoch = T

date.tag = ''
epida.tag = '_upto20221002' # timestamp for the epi data
data.tag = 'nyt' # whether new york times data or hopkins data are use

tag.model = ''

# scripts
source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code,'Fn_util.R'))
source(paste0(dir_code,'scorer.R'))



vtags = c('non.Omicron', 'Omicron')
truth00 = NULL
DAT.EPI00 = NULL
for(variant.tag in vtags){
  load(paste0(dir_data, 'truth_',data.tag, '_',variant.tag, epida.tag, '.RData'))
  truth$variant = variant.tag
  truth00 = rbind(truth00, truth)
  
  DAT.EPI = read.csv(paste0(dir_data, 'da_case_death_us_', data.tag, '_',variant.tag, epida.tag, '.csv')) %>% data.table()
  DAT.EPI$date = DAT.EPI$date %>% as.Date
  DAT.EPI$variant = variant.tag
  DAT.EPI00 = rbind(DAT.EPI00, DAT.EPI, fill=T)
  
  rm(truth, DAT.EPI)
}
truth = truth00; DAT.EPI = DAT.EPI00
rm(truth00, DAT.EPI00)

vax.start = as.Date('2020/12/14')
date.start = as.Date('2020/03/01') # start from 3/1/20
mob.type = 'business'


# read mobility data
DAT.MOB = read.csv(paste0(dir_data, 'da_mobility_us.csv')) %>% data.table()
DAT.MOB = DAT.MOB[data.type == mob.type]

# read vaccination data
DAT.VAC = read.csv(paste0(dir_data,'da_vx_perM_us_lagged_booster.csv')) %>% data.table()

# read seasonal trend
DAT.SN = read.csv(paste0(dir_data, 'est.sn.2000t2020_us.csv')) %>% data.table()


N = 1e6; # per 1 M
epi.model = 'SEIRSVimmLoss' # 'SEIRSV' # susceptible-exposed-infectious-recovered-susc
# for seeding: p.home.lwr - p.home.upr from home
stoch = T
# seasonality = T; # seasonality not good
start.date = '2020/3/1' %>% as.Date()

list.loc.names = c(state.name, "District of Columbia")
list.locs = gsub(' ', '', list.loc.names)
names(list.locs) = list.loc.names

locs = c('California', 'Florida', 'Iowa', 'Massachusetts', 'Michigan', 
         'New York', 'Pennsylvania', 'Texas', 'Washington', 'Wyoming')

# READ AND COMBINE RUNS BY THE ARIMAX MODELS - DONE LOCALLY DUE TO LARGE FILE SIZE
if(F){
  files = list.files(dir_res, full.names = T)
  files = files[grepl('.RData',files)]
  files = files[grepl('train.proj', files)]
  # files = files[grepl('Mean', files)] # transformed seasonal trend
  # compare performance across runs?
  length(files) # 15000 per state
  length(files)
  
  source(paste0(dir_code,'getCombinedRes_rproj_arima.R'))
}

# load the combined files
res.eval00 = NULL
res.score00 = NULL
for(loc.t in locs){
  try(load(paste0(dir_res, 'res.eval_', loc.t,date.tag,'.RData')))
  try(load(paste0(dir_res, 'res.score_', loc.t,date.tag,'.RData')))
  
  res.eval00 = rbind(res.eval00, res.eval)
  res.score00 = rbind(res.score00, res.score)
}

res.eval = res.eval00
res.score = res.score00
rm(res.eval00, res.score00)

# load the ensemble performance
load("./results/acc.base.v.opt_wilcoxon.RData")
load("./results/tda.base.v.opt_wilcoxon.RData")
tda.base.v.opt = tda.base.v.opt  %>%
  mutate(measure = factor(measure, levels = c('case','death'), labels = c('Cases','Deaths')))

met.score.t = 'logScore.loose'
met.pt.t = 'acc.loose1'

# the long-term n-wk ahead forecast using non.Omicron data extended to the Omicron period
# these should be excluded for evaluation
# criteria: truncated at end of Dec 2021
# can only do this for the n-wk ahead forecast and not the overall targets
res.eval = res.eval %>% filter(measure %in% c('Cases','Deaths')) %>% # rowwise() %>% 
  mutate(week.t = fcast.start.week + as.numeric(gsub('wk ahead', '', target)) * 7) %>%
  filter(!(variant == 'non.Omicron' & (!is.na(week.t) & as.Date(week.t) > as.Date('2021/12/31')))) %>% 
  mutate(wave.t = NULL)

res.score = res.score %>% # filter(loc == 'California' & seasonality == 'no' & fcast.deflat==.9 & measure == 'case' & metric == 'logScore.loose' & scenario == 'asIs') %>% 
  # rowwise() %>% 
  mutate(week.t = fcast.start.week + as.numeric(gsub('wk ahead', '', target)) * 7) %>%
  filter(!(variant == 'non.Omicron' & (!is.na(week.t) & as.Date(week.t) > as.Date('2021/12/31')))) %>% 
  mutate(wave.t = NULL)

# some arima model did not work for certain weeks at so have no forecast for many weeks
# check which models have the more complete set of forecasts
n.fcast.by.model = res.score %>% filter(target == 'peak week' & metric == met.score.t) %>%
  group_by(fcast.type) %>%
  summarise(n = n())
# the full model including vaccination data only work for ~1/3 of the time (634 vs 1588 of the base ARIMA model)
# arima, arimax.mob are most complete; a few were missing for arimax.sn, and arima.ms

models = n.fcast.by.model %>% filter(n > max(n) * .9) %>% .$fcast.type # missed no more than 10%
# check overall which arima models did the best 
res.score %>% filter(metric == met.score.t & fcast.type %in% models) %>%
  group_by(fcast.type) %>%
  summarise(mean.score = mean(mean)) %>%
  arrange(., -mean.score)

res.score %>% filter(metric == met.score.t & fcast.type %in% models) %>%
  group_by(fcast.type, measure) %>%
  summarise(mean.score =mean(mean)) %>%
  arrange(., measure, -mean.score)


res.eval %>% filter(metric == met.pt.t & fcast.type %in% models) %>%
  group_by(fcast.type) %>%
  summarise(mean.acc = mean(mean)) %>%
  arrange(., -mean.acc)

res.eval %>% filter(metric == met.pt.t & fcast.type %in% models) %>%
  group_by(fcast.type, measure) %>%
  summarise(mean.acc = mean(mean)) %>%
  arrange(., measure, -mean.acc)

# overall, it looks like the arimax.sn model did the best
# when tallied by measure, ranked 3rd for cases (arimax.mob model was the best for cases) and 1st for deaths per logScore
# ranked 1st for both cases and deaths per point prediction accuracy 

# => so use the arimax.sn model as the null model for comparison
res.score = res.score %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(measure = factor(measure, levels = c('case','death'), labels = c('Cases','Deaths')))

res.eval = res.eval %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) 

# cp the arima models
cp.score.arimas = rbind(res.score %>% filter(metric == met.score.t & fcast.type %in% models) %>%
                          group_by(fcast.type, measure) %>%
                          summarise(score = mean(mean) %>% round(., 2) %>% as.character()) %>%
                          mutate(fcast.type = factor(fcast.type, levels = c('arima', 'arimax.mob','arimax.sn','arimax.ms'),
                                                     labels = c('ARIMA','ARIMAX.MOB','ARIMAX.SN','ARIMAX.MOB.SN'))) %>%
                          dcast(., measure ~ fcast.type, value.var = 'score') %>%
                          mutate(target = 'all'),
                        res.score %>% filter(metric == met.score.t & fcast.type %in% models) %>%
  group_by(fcast.type, measure, target2) %>%
  summarise(score = mean(mean) %>% round(., 2) %>% as.character()) %>%
  rename(., target = target2) %>%
  mutate(fcast.type = factor(fcast.type, levels = c('arima', 'arimax.mob','arimax.sn','arimax.ms'),
                             labels = c('ARIMA','ARIMAX.MOB','ARIMAX.SN','ARIMAX.MOB.SN'))) %>%
  dcast(., measure + target ~ fcast.type, value.var = 'score')) %>%
  mutate(metric = 'Log score')

cp.acc.arimas = rbind(res.eval %>% filter(metric == met.pt.t & fcast.type %in% models) %>%
                          group_by(fcast.type, measure) %>%
                          summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                          mutate(fcast.type = factor(fcast.type, levels = c('arima', 'arimax.mob','arimax.sn','arimax.ms'),
                                                     labels = c('ARIMA','ARIMAX.MOB','ARIMAX.SN','ARIMAX.MOB.SN'))) %>%
                          dcast(., measure ~ fcast.type, value.var = 'acc') %>%
                          mutate(target = 'all'),
                        res.eval %>% filter(metric == met.pt.t & fcast.type %in% models) %>%
                          group_by(fcast.type, measure, target2) %>%
                          summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                          rename(., target = target2) %>%
                          mutate(fcast.type = factor(fcast.type, levels = c('arima', 'arimax.mob','arimax.sn','arimax.ms'),
                                                     labels = c('ARIMA','ARIMAX.MOB','ARIMAX.SN','ARIMAX.MOB.SN'))) %>%
                          dcast(., measure + target ~ fcast.type, value.var = 'acc')) %>%
  mutate(metric = 'Accuracy')

tab1 = rbind(cp.score.arimas, cp.acc.arimas) %>% 
  mutate(target = factor(target, levels = c('all', '1-8wk ahead', '9-16wk ahead', '17-26wk ahead','peak intensity', 'peak week', 'total')),
         metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>% 
  setcolorder(., c('target', 'metric', 'measure')) %>%
  arrange(., target, metric, measure)

m.best = 'arimax.sn'
score.arima = rbind(res.score %>% filter(metric == met.score.t & fcast.type == m.best) %>%
                      group_by(measure) %>%
                      summarise(score = mean(mean) %>% round(., 2) %>% as.character()) %>%
                      mutate(target = 'all'),
                    res.score %>% filter(metric == met.score.t & fcast.type == m.best) %>%
                      group_by(measure, target2) %>%
                      summarise(score = mean(mean) %>% round(., 2) %>% as.character()) %>%
                      rename(., target = target2)) %>%
  mutate(metric = 'Log score', method = m.best)

acc.arima = rbind(res.eval %>% filter(metric == met.pt.t & fcast.type == m.best) %>%
                    group_by( measure) %>%
                    summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                    mutate(target = 'all'),
                  res.eval %>% filter(metric == met.pt.t & fcast.type == m.best) %>%
                    group_by(measure, target2) %>%
                    summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                    rename(., target = target2)) %>%
  mutate(metric = 'Accuracy', method = m.best)
  

score.opt = rbind(tda.base.v.opt %>% filter(metric == met.score.t) %>%
                    group_by(measure) %>%
                    summarise(score = mean(opt) %>% round(.,2) %>% as.character()) %>%
                    mutate(target = 'all'), 
                  tda.base.v.opt %>% filter(metric == met.score.t) %>%
  group_by(measure, target2) %>%
  summarise(score = mean(opt) %>% round(.,2) %>% as.character()) %>%
  rename(., target = target2)) %>%
  mutate(metric = 'Log score', method = 'Best-performing')

score.base = rbind(tda.base.v.opt %>% filter(metric == met.score.t) %>%
                    group_by(measure) %>%
                    summarise(score = mean(base) %>% round(.,2) %>% as.character()) %>%
                    mutate(target = 'all'), 
                  tda.base.v.opt %>% filter(metric == met.score.t) %>%
                    group_by(measure, target2) %>%
                    summarise(score = mean(base) %>% round(.,2) %>% as.character()) %>%
                    rename(., target = target2)) %>%
  mutate(metric = 'Log score', method = 'Base')

acc.opt = rbind(acc.base.v.opt %>% filter(metric == met.pt.t) %>%
                  group_by(measure) %>%
                  summarise(acc = (mean(opt, na.rm =T) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                  mutate(target = 'all'), 
                acc.base.v.opt %>% filter(metric == met.pt.t) %>%
  group_by(measure, target2) %>%
  summarise(acc = (mean(opt, na.rm =T) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
  rename(., target = target2)) %>%
  mutate(metric = 'Accuracy', method = 'Best-performing')

acc.base = rbind(acc.base.v.opt %>% filter(metric == met.pt.t) %>%
                   group_by(measure) %>%
                   summarise(acc = (mean(base, na.rm =T) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                   mutate(target = 'all'), 
                 acc.base.v.opt %>% filter(metric == met.pt.t) %>%
                   group_by(measure, target2) %>%
                   summarise(acc = (mean(base, na.rm =T) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                   rename(., target = target2)) %>%
  mutate(metric = 'Accuracy', method = 'Base')

tab3 = rbind(rbind(score.arima, score.base, score.opt) %>%
  dcast(., metric + measure + target ~ method, value.var = 'score'),
  rbind(acc.arima, acc.base, acc.opt) %>%
    dcast(., metric + measure + target ~ method, value.var = 'acc')) %>%
  mutate(target = factor(target, levels = c('all', '1-8wk ahead', '9-16wk ahead', '17-26wk ahead','peak intensity', 'peak week', 'total')),
       metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  rename(., ARIMAX.SN = arimax.sn) %>%
  setcolorder(., c('metric','measure','target','ARIMAX.SN', 'Base','Best-performing')) %>%
  arrange(., target, metric, measure)
  

tab3v2 = merge(rbind(score.arima, score.base, score.opt) %>%
                 dcast(., measure + target ~ metric + method, value.var = 'score'),
               rbind(acc.arima, acc.base, acc.opt) %>%
                 dcast(., measure + target ~ metric + method, value.var = 'acc'), by = c('target','measure')) %>%
  mutate(target = factor(target, levels = c('all', '1-8wk ahead', '9-16wk ahead', '17-26wk ahead','peak intensity', 'peak week', 'total'))
  ) %>%
  arrange(., target, measure)

write.csv(tab1, paste0(dir_out,'TableS5_performance_arima_models.csv'), row.names = F)
write.csv(tab3v2, paste0(dir_out, 'TableS6_cp_arima_base_opt.csv'), row.names = F)

