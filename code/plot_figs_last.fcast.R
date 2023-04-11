# generate Figs 8 and 9 for the ms 
# 3/23/23 adding observed cases/deaths for comparison with the forecast

# note re data issue:
# NYT data used for the forecast appeared off for more recent months for some states, likely due to infrequent/irregular reporting period
# as such, use CDC data for evaluation instead, which look more reasonable
# however, for Washington, mortality data from CDC prior to Dec 2022 look off (timing shifted by ~3 months) and NYT data look ok
# so for WA mortality, use NYT data for evaluation instead
# for the 3/30/23 CDC data release, FL case and deaths for the week of 3/26/23 were recorded as 0, due to data issues
# thus, need to exclude the last week for FL

dir_data = './data/'
dir_code = './code/'
dir_res = './results/'
dir_out = './outputs/'

num_runs = 10
stoch = T

date.tag = ''
data.tag = 'mixed' # whether new york times data or hopkins data are use
tag.model = ''

# scripts
source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code,'getPlot.R'))
source(paste0(dir_code,'fn_util.R'))

# LABELS/ETC
vtags = c('non.Omicron', 'Omicron')
score.met.t = 'logScore.loose'
pt.met.t = 'acc.loose1' # for pt est
# combine all targets
state_vec = c('all', "California","Florida","Iowa","Massachusetts", "Michigan","New York","Pennsylvania","Texas","Washington","Wyoming")
state_lab = c('All', "California","Florida","Iowa","Massachusetts", "Michigan","New York","Pennsylvania","Texas","Washington","Wyoming")
variant_vec = c('non.Omicron', 'Omicron', 'all')
variant_lab = c('Pre-Omicron', 'Omicron', 'all')
sn_vec = c('no', 'y/fix', 'y/tf')
sn_lab = c('No seasonality', 'Fixed seasonality', 'Transformed seasonality')
sce_vec = c('asIs', 'newV')
sce_lab = c('No new variants', 'New variants')
measure_vec = c('case','death')
measure_lab = c('Cases', 'Deaths')
cp.df_vec = paste0('deflat', c('1.v.095', '1.v.09', '095.v.09'))
cp.df_lab = c('0.95 vs none', '0.9 vs none', '0.9 vs 0.95')
target2_vec = c('1-8wk ahead','9-16wk ahead','17-26wk ahead','peak week','peak intensity', 'total')
target2all_vec = c('all',target2_vec)
target2all_vec.lab = c('all combined',target2_vec)
cp.sce_vec = c('asIs.v.newV')
cp.sce_lab = 'none vs new variants'
cp.sce_lab = 'new variants vs baseline'
wave2_vec = c('wave2', 'Alpha', 'Delta', 'Pre-Omicron', 'Omicron','Omicron_BA.1', "Omicron_BA.2","Omicron_BA.2.12.1", "Omicron_nonBA.1o2")
wave2_lab = c('2nd wave', 'Alpha', 'Delta', 'Pre-Omicron','Omicron','Omicron BA.1',rep('After BA.1', 3))
cp.sn_vec = c('no.v.yfix','no.v.ytf','yfix.v.ytf')
cp.sn_lab = c('None vs fixed seasonality', 'None vs transformed seasonality', 'Fixed vs transformed seasonality')
cp.sn_lab = c('Fixed vs no seasonality', 'Transformed vs no seasonality', 'Transformed vs fixed seasonality')
measures_vec = c('Infections','Cases','Deaths', "Cumulative Infections", "Cumulative Cases", "Cumulative Deaths")

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
  DAT.EPI$date = DAT.EPI$date %>% as.Date
  
  DAT.EPI00 = rbind(DAT.EPI00, DAT.EPI)
  
  rm(truth, DAT.EPI)
}

truth = truth00; DAT.EPI = DAT.EPI00
rm(truth00, DAT.EPI00)


# look at the last forecast
start.date = date.start = '2020/3/1' %>% as.Date()
locs = c("California","Florida","Iowa","Massachusetts", "Michigan","New York","Pennsylvania","Texas","Washington","Wyoming")
res.train00 = NULL
res.proj00 = NULL
deflat_vec = .9 # c(.9, .95, 1)
for(deflat.t in deflat_vec){
  for(loc.t in locs){
    rm(res.train, res.proj)
    res.train = NULL
    res.proj = NULL
    
    try(load(paste0(dir_res, 'res.train', paste0('_df',deflat.t), '_', loc.t,date.tag,'.RData')))
    try(load(paste0(dir_res, 'res.proj', paste0('_df',deflat.t), '_', loc.t,date.tag,'.RData')))
    
    res.train00 = rbind(res.train00, res.train)
    res.proj00 = rbind(res.proj00, res.proj)
    
  }
}

res.train = res.train00
res.proj = res.proj00
rm(res.train00, res.proj00)

fcast.start.week_vec = res.proj$fcast.start.week %>% unique %>% sort
fcast.last = fcast.start.week_vec %>% max

res.proj$measure %>% unique

stat.tag = 'wilcoxon'
load(file = paste0(dir_res, 'stats.acc.base.v.opt_', stat.tag, '.RData'))  # Get the best-performing setting

# separate the ones with higher accuracy 
da.t = stats.acc.base.v.opt %>% mutate(sig = case_when(as.numeric(pvalue.wrs) < .05 ~ '*', 
                                                       T ~ '')) %>%  
  mutate(resp.sn = resp.sn %>% as.character()) %>% 
  replace_na(list(resp.sn = 'All', wave = 'All')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(# cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
    # measure = factor(measure, levels = measure_vec, labels = measure_lab),
    # scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
    # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
    target = factor(target, levels = c('all', target2_vec)),
    variant = factor(variant, levels = c('all', variant_vec), labels = c('All', variant_lab)),
    State = factor(loc, levels = state_vec, labels = state_lab))  %>% 
  # rowwise() %>% 
  # mutate(out = paste0((as.numeric(diff.median.wrs)*100) %>% round(0),
  #                    '% (',(as.numeric(diff.median.wrs.ci95lwr)*100) %>% round(0),'%, ',
  #                    (as.numeric(diff.median.wrs.ci95upr)*100) %>% round(0),'%) ', sig)) %>%
  mutate(base = (as.numeric(mean.base) * 100),
         opt = (as.numeric(mean.opt) * 100)) 
# resp season
tab2 = da.t %>% dplyr::filter(State != 'All' & resp.sn == 'Respiratory season' & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & variant == 'All') %>%
  dplyr::select(State, resp.sn, measure, target, opt) %>%
  dcast(., State ~ measure, value.var = 'opt') %>%
  rowwise() %>%
  mutate(mean.acc = mean(c(Cases, Deaths))) %>%
  dplyr::arrange(., by = -mean.acc)
locs.order = tab2$State
locs.high = head(locs.order, 5)
locs.low = tail(locs.order, 5) # tab2 %>% filter(mean.acc < median(tab2$mean.acc)) %>% .$State %>% as.character()

# summarize peak timing, peak intensity and cumulative totals
sce.t = 'newV'; sn.t = 'Transformed seasonality'
res.proj.t = res.proj %>% filter(seasonality == sn.t & scenario == sce.t &
                                   measure %in% c('Infections','Cases','Deaths') &
                                   fcast.start.week == fcast.last) %>%
            dcast(., loc + measure + Week.start ~ variable, value.var = 'value') %>%
  mutate(measure = factor(measure, levels = c('Infections','Cases','Deaths'))) %>%
  mutate(loc = factor(loc, levels = locs.order))

wks.fcast = res.proj.t$Week.start %>% unique()
res.train.t = res.train[seasonality == sn.t & variant == 'Omicron' &
                          fcast.start.week == fcast.last &
                          state %in% c('infection','case','death') &
                          !Week.start %in% wks.fcast] %>%
  mutate(measure = factor(state, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths')),
         variant = NULL) %>%
  dcast(., loc + measure + Week.start ~ variable, value.var = 'value')

wks.train = res.train.t$Week.start %>% unique() %>% sort
da.t = DAT.EPI[variant == 'Omicron'] %>% dcast(., state + date + year + week ~ data.type, value.var = 'value') %>% 
        melt(id.var = c('state', 'date','year','week')) %>% 
        data.table::copy() %>%
        setnames(c('state','date', 'variable', 'value'), c('loc','Week.start','measure','obs')) %>%
        mutate(measure = factor(measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths')) )

res.train.t = merge(res.train.t, da.t, by = c('loc', 'Week.start', 'measure'), all.x = T) %>%
  filter(Week.start >= fcast.last - 3 * 31) %>%
  mutate(loc = factor(loc, levels = locs.order))

# add observation
res.proj.t = merge(res.proj.t, da.t, by = c('loc', 'Week.start', 'measure'), all.x = T) %>%
  # filter(Week.start >= fcast.last - 3 * 31) %>%
  mutate(loc = factor(loc, levels = locs.order))

dates.t = seq(min(res.train.t$Week.start), max(res.proj.t$Week.start)+14, by = '4 weeks')

locs.t = locs.high
p1 = ggplot(res.train.t %>% filter(loc %in% locs.t)) +
  geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
  geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
  geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
  geom_line(data=res.proj.t %>% filter(loc %in% locs.t), aes(x = Week.start, y = median), color = 'red') +  # no ctrl
  geom_ribbon(data=res.proj.t %>% filter(loc %in% locs.t), aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'red', alpha = .1) +
  geom_ribbon(data=res.proj.t %>% filter(loc %in% locs.t), aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
  geom_vline(data = res.train.t %>% filter(loc %in% locs.t), aes(xintercept = max(res.train.t$Week.start)), linetype = 'dashed')+
  geom_point(data = res.train.t %>% filter(loc %in% locs.t), mapping = aes(x = Week.start, y=obs)) + 
  geom_point(data = res.proj.t %>% filter(loc %in% locs.t), mapping = aes(x = Week.start, y=obs)) + 
  facet_rep_wrap(~ loc +  measure, scales = 'free', repeat.tick.labels = T, ncol = 3, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
  labs(x = 'Week Start', y = 'Weekly number per 1 million people', title = '') +
  scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
  theme_minimal() +  theme.t6l # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
locs.t = locs.low
p2 = ggplot(res.train.t %>% filter(loc %in% locs.t)) +
  geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
  geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
  geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
  geom_line(data=res.proj.t %>% filter(loc %in% locs.t), aes(x = Week.start, y = median), color = 'red') +  # no ctrl
  geom_ribbon(data=res.proj.t %>% filter(loc %in% locs.t), aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'red', alpha = .1) +
  geom_ribbon(data=res.proj.t %>% filter(loc %in% locs.t), aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
  geom_vline(data = res.train.t %>% filter(loc %in% locs.t), aes(xintercept = max(res.train.t$Week.start)), linetype = 'dashed')+
  geom_point(data = res.train.t %>% filter(loc %in% locs.t), mapping = aes(x = Week.start, y=obs)) + 
  geom_point(data = res.proj.t %>% filter(loc %in% locs.t), mapping = aes(x = Week.start, y=obs)) + 
  facet_rep_wrap(~ loc +  measure, scales = 'free', repeat.tick.labels = T, ncol = 3, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
  labs(x = 'Week Start', y = 'Weekly number per 1 million people', title = '') +
  scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
  theme_minimal() +  theme.t6l # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))

pdf(paste0(dir_out,'Fig8_last.fcast_',data.tag,'.pdf'), width = 12, height = 8)
grid.arrange(p1, p2, ncol = 2, nrow = 1)
dev.off()



truth.t = truth %>% filter(wk.fcast == fcast.last & variant == 'Omicron')
da.t = DAT.EPI[variant == 'Omicron'] %>% dcast(., state + date + year + week ~ data.type, value.var = 'value') %>% 
  dplyr::filter(date >= min(res.proj.t$Week.start) & date <= max(res.proj.t$Week.start) &
                  state %in% unique(res.proj.t$loc)) %>% 
  melt(id.var = c('state', 'date','year','week')) %>% 
  data.table::copy() %>%
  setnames(c('state','date', 'variable', 'value'), c('loc','Week.start','measure','obs')) %>%
  group_by(loc, measure) %>% 
  summarise(obs = sum(obs)) %>% 
  mutate(obs = case_when(measure %in% c('infection') ~ obs / 1e6 * 100, # "Cumulative Cases",
                                measure %in% c('case') ~ obs / 1e6 * 100,
                                T ~ obs)) %>%
  mutate(measure = factor(measure, levels = c('infection','case','death'), labels = c('Infections (% population)', 'Cases (% population)', 'Deaths (per 1 million)')) )

tda = res.proj %>% filter(seasonality == sn.t & scenario == sce.t &
                            measure %in% c("Cumulative Infections", "Cumulative Cases", "Cumulative Deaths") &
                            fcast.start.week == fcast.last & Week.start == max(Week.start)) %>%
  mutate(value.show = case_when(measure %in% c('Cumulative Infections') ~ value / 1e6 * 100, # "Cumulative Cases",
                                measure %in% c('Cumulative Cases') ~ value / 1e6 * 100,
                                T ~ value)) %>%
  dcast(., loc + measure ~ variable, value.var = 'value.show') %>%
  mutate(measure = factor(measure, levels = c("Cumulative Infections", "Cumulative Cases", "Cumulative Deaths"),
                          labels = c('Infections (% population)', 'Cases (% population)', 'Deaths (per 1 million)'))) %>%
  mutate(loc = factor(loc, levels = locs.order)) # based on historical accuracy

# tda %>% filter(measure == 'Infections (% population)') %>% arrange(., by = -median) 
# tda %>% filter(measure == 'Deaths (per 1 million)') %>% arrange(., by = -median) 

p = ggplot(tda, aes(x = loc)) +
  geom_boxplot(aes(
    lower = iqr.lwr, 
    upper = iqr.upr, 
    middle = median, 
    ymin = ci95.lwr, 
    ymax = ci95.upr,
    group = loc),
    stat = "identity"
  ) + 
  geom_point(data = da.t, aes(x = loc, y = obs), color = 'red', pch = '*', size = 5) +
  facet_rep_wrap(~ measure, ncol = 4, scales = 'free_y',repeat.tick.labels = 'all') + 
  labs(x = '', y = 'Cumulative total') +
  theme_minimal() + theme.t6g

pdf(paste0(dir_out,'Fig9_last.fcast_totals_',data.tag,'.pdf'), width = 6, height = 3)
print(p)
dev.off()


