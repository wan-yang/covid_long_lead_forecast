# to evaluate the last forecast (real-time) generated for the long-lead forecast study
# 3/31/23

# note re data issue:
# NYT data used for the forecast appeared off for more recent months for some states, likely due to infrequent/irregular reporting period
# as such, use CDC data for evaluation instead, which look more reasonable
# however, for Washington, mortality data from CDC prior to Dec 2022 look off (timing shifted by ~3 months) and NYT data look ok
# so for WA mortality, use NYT data for evaluation instead
# for the 3/30/23 CDC data release, FL case and deaths for the week of 3/26/23 were recorded as 0, due to data issues
# thus, need to exclude the last week for FL

dir_data = './data/'
dir_code = './code/'
dir_res = paste0('./results/last.fcast/')
dir_out = paste0('./outputs/')

if(!file.exists(dir_out)) dir.create(dir_out)

num_runs = 10
stoch = T
date.tag = ''
data.tag = 'cdc' # use cdc data instead, 
# b/c NYT data were not updated regularly and with the cumulative tallies, the weekly tallies became irregular
tag.model = ''

ff.tag = '20221002' # last forecast made in the long-lead forecast study

# scripts
source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code,'Fn_util.R'))
source(paste0(dir_code,'scorer.R'))

vtags = c('non.Omicron', 'Omicron')
# however, for Washington, mortality data from CDC look off (timing shifted by ~3 months) and NYT data look ok
# so for WA mortality, use NYT data for evaluation instead
truth00 = NULL
DAT.EPI00 = NULL
for(variant.tag in vtags){
  load(paste0(dir_data, 'truth_','mixed', '_',variant.tag,'.RData'))  # generated using combined cdc/nyt data correcting for potential error
  truth$variant = variant.tag
  truth00 = rbind(truth00, truth)
  
  DAT.EPI_cdc = read.csv(paste0(dir_data, 'da_case_death_us_', 'cdc', '_',variant.tag,'.csv')) %>% data.table()
  DAT.EPI_cdc$date = DAT.EPI_cdc$date %>% as.Date
  DAT.EPI_nyt = read.csv(paste0(dir_data, 'da_case_death_us_', 'nyt' , '_',variant.tag,'.csv')) %>% data.table()
  DAT.EPI_nyt$date = DAT.EPI_nyt$date %>% as.Date
  
  date.t = as.Date('2022/12/1') # data after this week are about the same from both cdc and nyt
  DAT.EPI = rbind(DAT.EPI_cdc %>% filter(!(state == 'Washington' & data.type == 'death' & date <= date.t)),
                  DAT.EPI_nyt %>% filter((state == 'Washington' & data.type == 'death' & date <= date.t))) %>%
    mutate(date.1st = NULL)
  DAT.EPI$variant = variant.tag
  
  
  DAT.EPI00 = rbind(DAT.EPI00, DAT.EPI)
  
  rm(truth, DAT.EPI)
}

truth = truth00; DAT.EPI = DAT.EPI00
rm(truth00, DAT.EPI00)

vax.start = as.Date('2020/12/14')
date.start = as.Date('2020/03/01') # start from 3/1/20
start.date = '2020/3/1' %>% as.Date()

N = 1e6; # per 1 M
epi.model = 'SEIRSVimmLoss' # 'SEIRSV' # susceptible-exposed-infectious-recovered-susc
# for seeding: p.home.lwr - p.home.upr from home
stoch = T
# seasonality = T; # seasonality not good


fcast.deflat_vec = seq(.9, 1, by = .05)

locs = c('California', 'Florida', 'Iowa', 'Massachusetts', 'Michigan', 
         'New York', 'Pennsylvania', 'Texas', 'Washington', 'Wyoming')

# READ AND COMBINE RUNS BY THE ARIMAX MODELS - DONE LOCALLY DUE TO LARGE FILE SIZE
if(F){
  files = list.files(dir_res, full.names = T)
  files = files[grepl('.RData',files)]
  files = files[grepl('train.proj', files)]
  files = files[grepl(ff.tag, files)] # get the last fcast
  
  source(paste0(dir_code,'getCombinedRes_rproj.R'))
  
}

da.case = DAT.EPI %>% filter(date >= as.Date(ff.tag, '%Y%m%d') & state %in% locs & data.type == 'case' & variant == 'Omicron') %>%
  dcast(., date ~ state, value.var = 'value')


locs = c("California","Florida","Iowa","Massachusetts", "Michigan","New York","Pennsylvania","Texas","Washington","Wyoming")
res.eval00 = NULL
res.score00 = NULL
deflat_vec = .9 # c(.9, .95, 1)
for(loc.t in locs){
  rm(res.eval, res.score)
  res.eval = NULL
  res.score = NULL
  
  try(load(paste0(dir_res, 'res.eval_', loc.t,date.tag,'.RData')))
  try(load(paste0(dir_res, 'res.score_', loc.t,date.tag,'.RData')))
  
  res.eval00 = rbind(res.eval00, res.eval)
  res.score00 = rbind(res.score00, res.score)
  
}
res.eval = res.eval00
res.score = res.score00
rm(res.eval00, res.score00)

measures_vec = c('Infections','Cases','Deaths', "Cumulative Infections", "Cumulative Cases", "Cumulative Deaths")
met.score.t = 'logScore.loose'
met.pt.t = 'acc.loose1'

res.score = res.score %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) %>%
  mutate(measure = factor(measure, levels = c('case','death'), labels = c('Cases','Deaths')))

res.eval = res.eval %>% 
  mutate(target2 = factor(target, levels = c(paste0(1:26, 'wk ahead'), 'peak intensity', 'peak week', 'total'),
                          labels = c(rep('1-8wk ahead', 8), rep('9-16wk ahead', 8), rep('17-26wk ahead',10),'peak intensity', 'peak week', 'total'))) 

df.t = .9; sce.t = 'newV'; sn.t = 'Transformed seasonality'
res.score.t = res.score %>% filter(fcast.deflat == df.t & scenario == sce.t & seasonality == sn.t)
res.eval.t = res.eval %>% filter(fcast.deflat == df.t & scenario == sce.t & seasonality == sn.t & measure %in% c('Cases','Deaths'))

# cp the arima models
score.last = rbind(res.score.t %>% filter(metric == met.score.t) %>%
                        group_by(measure) %>%
                        summarise(score = mean(mean) %>% round(.,2) %>% as.character()) %>%
                        mutate(target = 'all', loc = 'All'), 
                      res.score.t %>% filter(metric == met.score.t) %>%
                        group_by(loc, measure) %>%
                        summarise(score = mean(mean) %>% round(.,2) %>% as.character()) %>%
                        mutate(target = 'all'), 
                      res.score.t %>% filter(metric == met.score.t) %>%
                        group_by(measure, target2) %>%
                        summarise(score = mean(mean) %>% round(.,2) %>% as.character()) %>%
                        rename(., target = target2) %>%
                        mutate(loc = 'All'),
                      res.score.t %>% filter(metric == met.score.t) %>%
                        group_by(loc, measure, target2) %>%
                        summarise(score = mean(mean) %>% round(.,2) %>% as.character()) %>%
                        rename(., target = target2)) %>%
  rename(., State = loc) %>%
  mutate(metric = 'Log score', target = factor(target, levels = c('all', '1-8wk ahead','9-16wk ahead','17-26wk ahead','peak intensity','peak week', 'total')))

tda = score.last %>% filter(target == 'all') %>% 
  dcast(.,State ~ measure, value.var = 'score')

tda2 = score.last %>% filter(State == 'All') %>% 
  dcast(., target ~ measure, value.var = 'score')

acc.last = rbind(res.eval.t %>% filter(metric == met.pt.t) %>%
                     group_by(measure) %>%
                     summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                     mutate(target = 'all', loc = 'All'), 
                   res.eval.t %>% filter(metric == met.pt.t) %>%
                     group_by(loc, measure) %>%
                     summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                     mutate(target = 'all'), 
                   res.eval.t %>% filter(metric == met.pt.t) %>%
                     group_by(measure, target2) %>%
                     summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                     rename(., target = target2) %>%
                     mutate(loc = 'All'),
                   res.eval.t %>% filter(metric == met.pt.t) %>%
                     group_by(loc, measure, target2) %>%
                     summarise(acc = (mean(mean) * 100) %>% round(.,0) %>% paste0(., '%')) %>%
                     rename(., target = target2)) %>%
  rename(., State = loc) %>%
  mutate(metric = 'Accuracy', target = factor(target, levels = c('all', '1-8wk ahead','9-16wk ahead','17-26wk ahead','peak intensity','peak week', 'total')))

# merge(score.last %>% dplyr::select(-metric), acc.last %>% dplyr::select(-metric), by = c('measure','State','target')) 

tab1 = merge(score.last %>% dcast(., State + target + metric ~ measure, value.var = 'score') %>% dplyr::select(-metric), 
             acc.last %>% dcast(., State + target + metric ~ measure, value.var = 'acc') %>% dplyr::select(-metric),
             by = c('State','target'), suffix = c('.score','.acc')) %>%
  mutate(target = factor(target, levels = c('all', '1-8wk ahead', '9-16wk ahead', '17-26wk ahead','peak intensity', 'peak week', 'total'))) %>%
  arrange(., State, target)

write.csv(tab1, paste0(dir_out, 'TableS7_stats_last.fcast.csv'), row.names = F)



