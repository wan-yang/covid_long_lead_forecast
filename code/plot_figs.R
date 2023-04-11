# generate plots for the ms

# do not load the full data, otherwise too big and slow
# note: use local environment for plotting
# ggplot()'s aes seems to look for variables in the data frame first, 
# and if not found, then in the global environment. 
# so if local variable is preferred, don't put it in the "aes()" if those variables are not in the data frame

dir_data = './data/'
dir_code = './code/'
dir_res = './results/'
dir_out = './outputs/'

if(! file.exists(dir_out))  
  dir.create(dir_out,recursive = T)


dtag = '_upto20221002' # time stamp for data used for the study
num_runs = 10
stoch = T

date.tag = ''

data.tag = 'nyt' # whether new york times data or hopkins data are use

tag.model = ''


# scripts
source(paste0(dir_code, 'loadPackages.R'))
source(paste0(dir_code,'getPlot.R'))
source(paste0(dir_code,'fn_util.R'))


# load pairwise diff
load(paste0(dir_res, "tda.cp.sn_wilcoxon.RData"))
load(paste0(dir_res, "tda.cp.sce_wilcoxon.RData"))
load(paste0(dir_res, "tda.cp.deflat_wilcoxon.RData"))
load(paste0(dir_res, "tda.base.v.opt_wilcoxon.RData"))

# load pt estimate evaluations
load(paste0(dir_res, "acc.cp.sn_wilcoxon.RData"))
load(paste0(dir_res, "acc.cp.sce_wilcoxon.RData"))
load(paste0(dir_res, "acc.cp.deflat_wilcoxon.RData"))
load(paste0(dir_res, "acc.base.v.opt_wilcoxon.RData"))

# load stats
load(paste0(dir_res, "stats.cp.sn_wilcoxon.RData"))
load(paste0(dir_res, "stats.cp.sce_wilcoxon.RData"))
load(paste0(dir_res, "stats.cp.deflat_wilcoxon.RData"))
load(paste0(dir_res, "stats.base.v.opt_wilcoxon.RData"))

load(paste0(dir_res, "stats.acc.cp.sn_wilcoxon.RData"))
load(paste0(dir_res, "stats.acc.cp.sce_wilcoxon.RData"))
load(paste0(dir_res, "stats.acc.cp.deflat_wilcoxon.RData"))
load(paste0(dir_res, "stats.acc.base.v.opt_wilcoxon.RData"))

# LABELS/ETC
vtags = c('non.Omicron', 'Omicron')
measures_vec = c('Infections','Cases','Deaths')

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

sz = 2.5

# read in data
truth00 = NULL
DAT.EPI00 = NULL
for(variant.tag in vtags){
  load(paste0(dir_data, 'truth_',data.tag, '_',variant.tag, dtag,'.RData'))
  truth$variant = variant.tag
  truth00 = rbind(truth00, truth)
  
  DAT.EPI = read.csv(paste0(dir_data, 'da_case_death_us_', data.tag, '_',variant.tag, dtag,'.csv')) %>% data.table()
  DAT.EPI$date = DAT.EPI$date %>% as.Date
  DAT.EPI$variant = variant.tag
  DAT.EPI00 = rbind(DAT.EPI00, DAT.EPI, fill=T)
  
  rm(truth, DAT.EPI)
}
truth = truth00; DAT.EPI = DAT.EPI00
rm(truth00, DAT.EPI00)



########################################################################
# FIG 1
# map of the 10 states
library(usmap)
# plot_usmap(regions = "state")
library(ggplot2)
date.start = as.Date('2020/3/1')
date.end = as.Date('2022/9/30') #  max(DAT.EPI$date)
N = 1e6  # per 1 million people

da = read.csv(paste0(dir_data, 'da_case_death_us_nyt',dtag,'.csv')) %>% data.table()
da$date = da$date %>%  as.Date
da$value = da$value %>% as.numeric()

state_total = da %>% filter(date >= date.start & date <= date.end) %>% dplyr::group_by(state, data.type) %>%
  dplyr::summarise(total = sum(value) / N * 100) %>% 
  # setnames(., 'state', 'region') %>%
  dcast(., state ~ data.type, value.var = 'total') %>% 
  dplyr::filter(state %in% locs)

pp1 = plot_usmap(regions = "state", data = state_total, values = 'case') +
  scale_fill_distiller(palette = "RdPu", direction = 1, name = "Incidence\nrate (%)", na.value = alpha('grey', .1)) +
  labs(x = '', y = '', title = paste0('(A) Cumulative incidence rate, ', format(date.start, '%b %Y'), ' - ',format(date.end, '%b %Y'))) + 
  theme_bw() + theme.t8a + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())

pp2 = plot_usmap(regions = "state", data = state_total, values = 'death') +
  scale_fill_distiller(palette = "RdPu", direction = 1, name = "Mortality\nrate (%)", na.value = alpha('grey', .1)) +
  labs(x = '', y = '', title = paste0('(B) Cumulative mortality rate, ', format(date.start, '%b %Y'), ' - ',format(date.end, '%b %Y'))) + 
  theme_bw() + theme.t8a + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())


tda = da %>% filter(date >= date.start & date <= date.end & state %in% locs & data.type == 'case')
pp3 = ggplot(tda, 
             aes(x = date, y = value, color = state)) +
  geom_line() +
  scale_x_date(breaks = seq(min(tda$date), max(tda$date), by = '3 month'),
               labels = format(seq(min(tda$date), max(tda$date), by = '3 month'),'%m/%d/%Y')) +
  labs(x = '', y = 'Weekly number of cases per 1 million people', color = 'State', title = '(C) Reported weekly number of COVID-19 cases') +
  theme_minimal() +  theme.t8b 

tda = da %>% filter(date >= date.start & date <= date.end & state %in% locs & data.type == 'death')
pp4 = ggplot(tda, 
             aes(x = date, y = value, color = state)) +
  geom_line() +
  scale_x_date(breaks = seq(min(tda$date), max(tda$date), by = '3 month'),
               labels = format(seq(min(tda$date), max(tda$date), by = '3 month'),'%m/%d/%Y')) +
  labs(x = '', y = 'Weekly number of cases per 1 million people', color = 'State', title = '(D) Reported weekly number of COVID-19 deaths') +
  theme_minimal() +  theme.t8b 
pdf(paste0(dir_out, 'Fig1_map.pdf'), width = 7.5, height = 6)
grid.arrange(pp1, pp2, pp3, pp4, nrow = 2, ncol = 2)
dev.off()
########################################################################


########################################################################
# Fig 2 - DONE LOCALLY, B/C FILES TOO LARGE & NOT UPLOADED
start.date = date.start = '2020/3/1' %>% as.Date()
loc.ex = 'New York'
res.train00 = NULL
res.proj00 = NULL
deflat_vec = c(.9, .95, 1)
for(deflat.t in deflat_vec){
  for(loc.t in loc.ex){
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

{
  da.full.t = DAT.EPI[state == loc.ex] %>% dcast(., variant + date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= date.start]
  
  loc.t = loc.ex; ifcast.ex = 1
  fcast.start.week.t = fcast.start.week_vec[ifcast.ex]
  sce.t = 'asIs'
  mea.t = 'Cases'
  seasonality.t = 'Seasonality assumed'
  # res.proj.t = res.proj[loc == gsub(' ','', loc.t) & fcast.start.month == fcast.start.month.t]
  res.proj.t = res.proj[loc == loc.t & fcast.start.week == fcast.start.week.t] #  & fcast.deflat == deflat.t
  
  wks.fcast = res.proj.t$Week.start %>% unique()
  res.train.t = res.train[loc == loc.t & # loc == gsub(' ','', loc.t) & 
                            seasonality == seasonality.t &
                            # fcast.deflat == deflat.t &
                            fcast.start.week == fcast.start.week.t & 
                            state %in% c('infection','case','death') &
                            !Week.start %in% wks.fcast]
  # duplicated(res.train.t) %>% any
  # res.train.t98 %>% filter(fcast.deflat == .9 & seasonality == 'No seasonality' & state == 'death' & variable == 'mean' & Week.start == as.Date('2021-11-14'))
  
  res.train.t = dcast(res.train.t, variant + loc + fcast.deflat + seasonality + state + Week.start ~ variable, value.var = 'value')
  res.train.t$measure = factor(res.train.t$state, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  wks.train = res.train.t$Week.start %>% unique() %>% sort
  da.t = da.full.t[date <= max(wks.train)] %>% melt(id.var = c('variant', 'date','year','week')) %>% 
    setnames(c('date', 'variable', 'value'), c('Week.start','measure','obs'))
  da.t$measure = factor(da.t$measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  res.train.t = merge(res.train.t, da.t, by = c('variant', 'Week.start', 'measure'), all.x = T)
  
  da.t = da.full.t[date %in% wks.fcast] %>% melt(id.var = c('variant', 'date','year','week')) %>% 
    setnames(c('date', 'variable', 'value'), c('Week.start','measure','obs'))
  da.t$measure = factor(da.t$measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  res.proj.t = dcast(res.proj.t, variant + scenario + loc + fcast.deflat +  seasonality + measure + Week.start ~ variable, value.var = 'value')
  res.proj.t = merge(res.proj.t, da.t, by = c('variant', 'Week.start', 'measure'), all.x = T)
  
  res.train.t$fcast.deflat = factor(res.train.t$fcast.deflat, levels = c(1, .95, .9), labels = c('No deflation', 'Deflation = 0.95', 'Deflation = 0.9'))
  res.proj.t$fcast.deflat = factor(res.proj.t$fcast.deflat, levels = c(1, .95, .9), labels = c('No deflation', 'Deflation = 0.95', 'Deflation = 0.9'))
  
  pp1 = getPlotProj9580deflat1sn(train.t = res.train.t[measure %in% measures_vec & seasonality == seasonality.t], 
                                    proj.t = res.proj.t %>% filter(seasonality == seasonality.t & measure %in% measures_vec & scenario == sce.t), #  & Week.start < as.Date('2021/12/1')
                                    mea.t = mea.t,
                                 ptitle.t = paste0('(A) Example forecasts comparing deflation settings\n(', loc.t, '; ', mea.t,'; no new variants; fixed seasonality; ', 
                                                   fcast.start.week.t, ')'),
                                 ylab.t = '',
                                 theme.t = theme.t6
                                            )
}

{
  da.full.t = DAT.EPI[state == loc.ex] %>% dcast(., variant + date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= date.start]
  
  loc.t = loc.ex; ifcast.ex = 75; # ifcast.ex = 74
  fcast.start.week.t = fcast.start.week_vec[ifcast.ex]
  mea.t = 'Cases'
  seasonality.t = 'Seasonality assumed'
  deflat.t = .9
  # res.proj.t = res.proj[loc == gsub(' ','', loc.t) & fcast.start.month == fcast.start.month.t]
  res.proj.t = res.proj[seasonality==seasonality.t & fcast.deflat == deflat.t & loc == loc.t & fcast.start.week == fcast.start.week.t] #  & fcast.deflat == deflat.t
  
  wks.fcast = res.proj.t$Week.start %>% unique()
  res.train.t = res.train[loc == loc.t & # loc == gsub(' ','', loc.t) & 
                            seasonality == seasonality.t &
                            fcast.deflat == deflat.t &
                            fcast.start.week == fcast.start.week.t & 
                            state %in% c('infection','case','death') &
                            !Week.start %in% wks.fcast]
  # duplicated(res.train.t) %>% any
  # res.train.t98 %>% filter(fcast.deflat == .9 & seasonality == 'No seasonality' & state == 'death' & variable == 'mean' & Week.start == as.Date('2021-11-14'))
  
  res.train.t = dcast(res.train.t, variant + loc + fcast.deflat + seasonality + state + Week.start ~ variable, value.var = 'value')
  res.train.t$measure = factor(res.train.t$state, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  wks.train = res.train.t$Week.start %>% unique() %>% sort
  da.t = da.full.t[date <= max(wks.train)] %>% melt(id.var = c('variant', 'date','year','week')) %>% 
    setnames(c('date', 'variable', 'value'), c('Week.start','measure','obs'))
  da.t$measure = factor(da.t$measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  res.train.t = merge(res.train.t, da.t, by = c('variant', 'Week.start', 'measure'), all.x = T)
  
  da.t = da.full.t[date %in% wks.fcast] %>% melt(id.var = c('variant', 'date','year','week')) %>% 
    setnames(c('date', 'variable', 'value'), c('Week.start','measure','obs'))
  da.t$measure = factor(da.t$measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  res.proj.t = dcast(res.proj.t, variant + scenario + loc + fcast.deflat +  seasonality + measure + Week.start ~ variable, value.var = 'value')
  res.proj.t = merge(res.proj.t, da.t, by = c('variant', 'Week.start', 'measure'), all.x = T)
  
  res.proj.t$scenario = factor(res.proj.t$scenario, levels = sce_vec, labels = sce_lab)
  
  pp2= getPlotProj9580Sce1sn(train.t = res.train.t[measure %in% measures_vec], 
                                 proj.t = res.proj.t %>% filter(measure %in% measures_vec), #  & Week.start < as.Date('2021/12/1')
                                 mea.t = mea.t,
                                 ptitle.t = paste0('(B) Example forecasts comparing new variants settings\n(', loc.t, '; ', mea.t,'; deflation = .9; fixed seasonality; ', 
                                                   fcast.start.week.t, ')'),
                             theme.t = theme.t6
  )
}

{
  da.full.t = DAT.EPI[state == loc.ex] %>% dcast(., variant + date + year + week ~ data.type, value.var = 'value')
  da.full.t$date = da.full.t$date %>% as.Date()
  da.full.t = da.full.t[date >= date.start]
  
  loc.t = loc.ex; ifcast.ex = 1
  fcast.start.week.t = fcast.start.week_vec[ifcast.ex]
  mea.t = 'Cases'
  sce.t = 'newV'
  deflat.t = .9
  # res.proj.t = res.proj[loc == gsub(' ','', loc.t) & fcast.start.month == fcast.start.month.t]
  res.proj.t = res.proj[scenario==sce.t & fcast.deflat == deflat.t & loc == loc.t & fcast.start.week == fcast.start.week.t] #  & fcast.deflat == deflat.t
  
  wks.fcast = res.proj.t$Week.start %>% unique()
  res.train.t = res.train[loc == loc.t & # loc == gsub(' ','', loc.t) & 
                            fcast.deflat == deflat.t &
                            fcast.start.week == fcast.start.week.t & 
                            state %in% c('infection','case','death') &
                            !Week.start %in% wks.fcast]
  # duplicated(res.train.t) %>% any
  # res.train.t98 %>% filter(fcast.deflat == .9 & seasonality == 'No seasonality' & state == 'death' & variable == 'mean' & Week.start == as.Date('2021-11-14'))
  
  res.train.t = dcast(res.train.t, variant + loc + fcast.deflat + seasonality + state + Week.start ~ variable, value.var = 'value')
  res.train.t$measure = factor(res.train.t$state, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  wks.train = res.train.t$Week.start %>% unique() %>% sort
  da.t = da.full.t[date <= max(wks.train)] %>% melt(id.var = c('variant', 'date','year','week')) %>% 
    setnames(c('date', 'variable', 'value'), c('Week.start','measure','obs'))
  da.t$measure = factor(da.t$measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  res.train.t = merge(res.train.t, da.t, by = c('variant', 'Week.start', 'measure'), all.x = T)
  
  da.t = da.full.t[date %in% wks.fcast] %>% melt(id.var = c('variant', 'date','year','week')) %>% 
    setnames(c('date', 'variable', 'value'), c('Week.start','measure','obs'))
  da.t$measure = factor(da.t$measure, levels = c('infection','case','death'), labels = c('Infections','Cases','Deaths'))
  res.proj.t = dcast(res.proj.t, variant + scenario + loc + fcast.deflat +  seasonality + measure + Week.start ~ variable, value.var = 'value')
  res.proj.t = merge(res.proj.t, da.t, by = c('variant', 'Week.start', 'measure'), all.x = T)
  
  res.train.t$seasonality = factor(res.train.t$seasonality, levels = c('No seasonality', 'Seasonality assumed','Transformed seasonality'), labels = sn_lab)
  res.proj.t$seasonality = factor(res.proj.t$seasonality, levels = c('No seasonality', 'Seasonality assumed','Transformed seasonality'), labels = sn_lab)
  
  pp3= getPlotProj9580Sns(train.t = res.train.t[measure %in% measures_vec], 
                             proj.t = res.proj.t %>% filter(measure %in% measures_vec), #  & Week.start < as.Date('2021/12/1')
                             mea.t = mea.t,
                             ptitle.t = paste0('(C) Example forecasts comparing seasonality settings\n(', loc.t, '; ', mea.t,'; deflation = .9; new variants; ', 
                                               fcast.start.week.t, ')'),
                          ylab.t = '',
                          theme.t = theme.t6
  )
}



pdf(paste0(dir_out, 'Fig2_example_3strategies.pdf'), width = 7, height = 7)
grid.arrange(pp1, pp2, pp3, ncol = 1, nrow = 3)
dev.off()
########################################################################


########################################################################
# FIG 3
# DEFLAT: THIS OPTIMIZATION SHOULD BE FOR ALL SEASONALITY, SCENARIO SETTINGS
# 1. plot comparison of deflation settings, by targets
da.t = tda.cp.deflat %>%  filter(metric == score.met.t) %>% # %>% filter(scenario == 'newV') # & seasonality == 'no'
  mutate(cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         target2 = factor(target2, levels = target2_vec),
         variant = factor(variant, levels = variant_vec, labels = variant_lab))
da.acc.t = stats.acc.cp.deflat %>% filter(metric == pt.met.t & measure %in% c('Cases','Deaths')) %>% 
  mutate(cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
         # measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec))

# heatmap
stats.t = stats.cp.deflat %>% filter(metric == score.met.t) %>% 
  mutate(cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec, 
                         labels = target2all_vec.lab))
stats.acc.t = stats.acc.cp.deflat %>% filter(metric == pt.met.t & measure %in% measure_lab) %>% 
  mutate(cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
         # measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec, 
                         labels = target2all_vec.lab))
# all combined
tda = stats.t %>% filter(target == 'all combined'  & variant == 'all' & loc != 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  mutate(value = diff.median.wrs %>% as.numeric(), f1 = measure, f2 = cp, f3 = paste0(scenario,', ', seasonality)) %>%
  mutate(f3 = factor(f3, levels = c(paste0('No new variants, ',sn_lab), paste0('New variants, ',sn_lab)))) #,
# labels = c(paste0('No new variants\n',sn_lab), paste0('New variants\n',sn_lab))))
pp1 = getPlotHeatmap.allTargets(tda, title.t = paste0('(', LETTERS[1], ') Impact on log score (all targets combined)'),
                                fill.name.t = 'Difference\nin log score',
                                ncol.t = 3, theme.tt = theme.t6k)
tda = stats.acc.t %>% filter(target == 'all combined'  & variant == 'all' & loc != 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  mutate(value = (diff.median.wrs %>% as.numeric()) * 100, f1 = measure, f2 = cp, f3 = paste0(scenario,', ', seasonality)) %>%
  mutate(f3 = factor(f3, levels = c(paste0('No new variants, ',sn_lab), paste0('New variants, ',sn_lab)))) #,
# labels = c(paste0('No new variants\n',sn_lab), paste0('New variants\n',sn_lab))))
pp2 = getPlotHeatmap.allTargets(tda, title.t = paste0('(', LETTERS[2], ') Impact on point prediction accuracy (all targets combined)'),
                                fill.name.t = 'Difference\nin accuracy (%)',
                                ncol.t = 3, theme.tt = theme.t6k)
p = ggarrange(pp1, pp2, ncol=1, nrow=2, common.legend = F, legend="right", heights = c(1, 1))
pdf(paste0(dir_out,'Fig3_cp.df.pdf'), width = 7, height = 6)
print(p)
dev.off()
########################################################################



########################################################################
# FIG 4
# SCENARIOS: SHOULD IMPROVE FOR ALL SEASONALITY SETTING
# compare scenario
df.opt = .9
stats.t = stats.cp.sce %>% filter(fcast.deflat == df.opt & metric == score.met.t) %>% 
  mutate(cp = factor(cp, levels = cp.sce_vec, labels = cp.sce_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         # scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec, 
                         labels = target2all_vec.lab))
stats.acc.t = stats.acc.cp.sce %>% filter(fcast.deflat == df.opt & metric == pt.met.t & measure %in% measure_lab) %>% 
  mutate(cp = factor(cp, levels = cp.sce_vec, labels = cp.sce_lab), 
         # measure = factor(measure, levels = measure_vec, labels = measure_lab),
         # scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec, 
                         labels = target2all_vec.lab))

# all targets, by wave
tda = rbind(stats.t %>% filter(variant == 'Pre-Omicron' & loc != 'all' & target == 'all combined' & !is.na(wave) & is.na(resp.sn)),
            stats.t %>% filter(variant == 'Omicron' & loc != 'all' & target == 'all combined' & is.na(wave) & is.na(resp.sn)) %>%
              mutate(wave = 'Omicron')) %>% 
  # filter(target == 'all combined')  %>% 
  mutate(value = diff.median.wrs %>% as.numeric(), 
         wave = factor(wave, levels = c('wave2','Alpha','Delta', 'Omicron'), 
                       labels =  c('2nd wave','Alpha','Delta', 'Omicron')),
         f1 = measure, f2 = seasonality, f3 = wave)
pp1 = getPlotHeatmap.allTargets(tda, title.t = paste0('(', LETTERS[1], ') Impact on log score'),
                                fill.name.t = 'Difference\nin log score',
                                ncol.t = 3, theme.tt = theme.t6k)

tda = rbind(stats.acc.t %>% filter(variant == 'Pre-Omicron' & loc != 'all' & target == 'all combined' & !is.na(wave) & is.na(resp.sn)),
            stats.acc.t %>% filter(variant == 'Omicron' & loc != 'all' & target == 'all combined' & is.na(wave) & is.na(resp.sn)) %>%
              mutate(wave = 'Omicron')) %>% 
  # filter(target == 'all combined')  %>% 
  mutate(value = (diff.median.wrs %>% as.numeric()) * 100, 
         wave = factor(wave, levels = c('wave2','Alpha','Delta', 'Omicron'), 
                       labels =  c('2nd wave','Alpha','Delta', 'Omicron')),
         f1 = measure, f2 = seasonality, f3 = wave)
pp2 = getPlotHeatmap.allTargets(tda, title.t = paste0('(', LETTERS[2], ') Impact on point prediction accuracy'),
                                fill.name.t = 'Difference\nin accuracy (%)',
                                ncol.t = 3, theme.tt = theme.t6k)
p = ggarrange(pp1, pp2, ncol=1, nrow=2, common.legend = F, legend="right")
pdf(paste0(dir_out,'Fig4_df',df.opt,'_cp.sce.bywave.pdf'), width = 8, height = 6)
print(p)
dev.off()

# check stats
stats.t %>% filter(variant == 'Omicron' & target == 'all combined' & loc == 'all' & 
                     measure == 'Deaths' & is.na(wave) & is.na(resp.sn))

stats.t %>% filter(variant == 'Omicron' & target == 'all combined' & loc == 'all' & 
                     measure == 'Cases' & is.na(wave) & is.na(resp.sn))

out1s = stats.t %>% filter(variant == 'Omicron' & target == 'all combined' & loc == 'all' & 
                            is.na(wave) & is.na(resp.sn))
out1a = stats.acc.t %>% filter(variant == 'Omicron' & target == 'all combined' & loc == 'all' & 
                             is.na(wave) & is.na(resp.sn))
out1s %>% filter(measure == 'Cases') %>% .$perc.diff %>% as.numeric() %>% range
out1s %>% filter(measure == 'Deaths') %>% .$perc.diff %>% as.numeric() %>% range

out1a %>% filter(measure == 'Cases') %>% .$perc.diff %>% as.numeric() %>% range
out1a %>% filter(measure == 'Deaths') %>% .$perc.diff %>% as.numeric() %>% range

################################################################################################
# FIG 5
# SEASONALITY
df.opt = .9
sce.opt = 'newV'

stats.t = stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & metric == score.met.t) %>% 
  mutate(cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec, 
                         labels = target2all_vec.lab))
stats.acc.t = stats.acc.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt  & metric == pt.met.t & measure %in% measure_lab) %>% 
  mutate(cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
         # measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         variant = factor(variant, levels = variant_vec, labels = variant_lab),
         target = factor(target, levels = target2all_vec, 
                         labels = target2all_vec.lab))

tda = stats.t %>% filter(variant == 'all' & loc != 'all' & is.na(wave) & !is.na(resp.sn)) %>% 
  mutate(value = diff.median.wrs %>% as.numeric(), f1 = measure, f2 = resp.sn, f3 = cp)
pp1 = getPlotHeatmap(tda, title.t = paste0('(', LETTERS[1], ') Impact on log score'),
                     fill.name.t = 'Difference\nin log score',
                     ncol.t = 6, theme.tt = theme.t6k)


tda = stats.acc.t %>% filter(variant == 'all' & loc != 'all' & is.na(wave) & !is.na(resp.sn)) %>% 
  mutate(value = (diff.median.wrs %>% as.numeric())*100, f1 = measure, f2 = resp.sn, f3 = cp)
# labels = c(paste0('No new variants\n',sn_lab), paste0('New variants\n',sn_lab))))
pp2 = getPlotHeatmap(tda, title.t = paste0('(', LETTERS[2], ') Impact on point prediction accuracy'),
                     fill.name.t = 'Difference\nin accuracy (%)',
                     ncol.t = 6, theme.tt = theme.t6k)
p = ggarrange(pp1, pp2, ncol=1, nrow=2, common.legend = F, legend="right", heights = c(1, 1))
pdf(paste0(dir_out,'Fig5_df',df.opt,'_',sce.opt,'_cp.sn.pdf'), width = 12.5, height = 8)
print(p)
dev.off()
################################################################################################


################################################################################################
# FIG 6
# PLOT THE OPTIMIZED VERSON THE BASELINE
sn.opt = 'y/tf'
df.base = 1; sce.base = 'asIs'; sn.base = 'no'
tda.t = tda.base.v.opt %>% filter(metric == score.met.t) %>%
  mutate(measure = factor(measure, levels = measure_vec, labels = measure_lab),
         target2 = factor(target2, levels = target2_vec),
         variant = factor(variant, levels = variant_vec, labels = variant_lab)) %>%
  mutate(wave2 = case_when(is.na(wave) & variant == 'Pre-Omicron' ~ 'Pre-Omicron',
                           is.na(wave) & variant == 'Omicron' ~ 'Omicron',
                           T ~ wave)) %>%
  mutate(wave2 = factor(wave2, levels = wave2_vec, labels = wave2_lab))

agg.tag = '_byMeasure_byVariantObyRespSn_byTarget'
# first plot it by variant
ylim.t = tda.t %>% group_by(target2, loc, measure, variant) %>% 
  summarise(count = n(), 
            ymin = quantile(diff, .05), ymax = quantile(diff, .95),
            ymin2 = quantile(diff, .25) * 1.5, ymax2 = quantile(diff, .75) * 1.5) %>% ungroup()
ymin.t = ylim.t %>% dplyr:: select(ymin, ymin2) %>% unlist %>% min
ymax.t = ylim.t %>% dplyr:: select(ymax, ymax2) %>% unlist %>% max
count.t = ylim.t %>% group_by(target2, variant,  measure) %>%
  summarise(n.min = min(count), n.max = max(count)) %>%
  rowwise() %>%
  mutate(lab = paste0(n.min, '-', n.max))%>% ungroup()

pp1 = ggplot(tda.t, environment = environment()) +
  geom_boxplot(aes(x = target2, y = diff, color = loc), size = .3, outlier.shape = NA) + 
  geom_text(aes(x = target2,label = lab), y = ymax.t, size = sz, data = count.t %>% filter(target2 == '1-8wk ahead'), vjust = 1, hjust = .4) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .8, size = sz, data = count.t %>% filter(target2 == '9-16wk ahead'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .6, size = sz, data = count.t %>% filter(target2 == '17-26wk ahead'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t, size = sz, data = count.t %>% filter(target2 == 'peak week'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .8, size = sz, data = count.t %>% filter(target2 == 'peak intensity'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .6, size = sz, data = count.t %>% filter(target2 == 'total'), vjust = 1, hjust = .6) +
  geom_hline(yintercept = 0, color = 'black', size = .08) +
  facet_rep_wrap(~ measure + variant, ncol = 5, # , scales = 'free_y'
                 repeat.tick.labels = T) + # , labeller = label_wrap_gen(multi_line=FALSE)
  labs(x = '', y = 'Diffence in log score', color = 'State') +
  ggtitle(paste0('(', LETTERS[1], ') Baseline vs optimized setting, by period')) + 
  guides(color = guide_legend(override.aes = list(size = 0.3))) +
  coord_cartesian(ylim = c(ymin.t, ymax.t)) +
  theme_minimal() +theme.t6c

# plot by resp sn
ylim.t = tda.t %>% group_by(target2, loc, measure, resp.sn) %>% 
  summarise(count = n(), 
            ymin = quantile(diff, .05), ymax = quantile(diff, .95),
            ymin2 = quantile(diff, .25) * 1.5, ymax2 = quantile(diff, .75) * 1.5) %>% ungroup()
ymin.t = ylim.t %>% dplyr:: select(ymin, ymin2) %>% unlist %>% min
ymax.t = ylim.t %>% dplyr:: select(ymax, ymax2) %>% unlist %>% max
count.t = ylim.t %>% group_by(target2, resp.sn,  measure) %>%
  summarise(n.min = min(count), n.max = max(count)) %>%
  rowwise() %>%
  mutate(lab = paste0(n.min, '-', n.max))%>% ungroup()

pp2 = ggplot(tda.t, environment = environment()) +
  geom_boxplot(aes(x = target2, y = diff, color = loc), size = .3, outlier.shape = NA) + 
  geom_text(aes(x = target2,label = lab), y = ymax.t, size = sz, data = count.t %>% filter(target2 == '1-8wk ahead'), vjust = 1, hjust = .4) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .8, size = sz, data = count.t %>% filter(target2 == '9-16wk ahead'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .6, size = sz, data = count.t %>% filter(target2 == '17-26wk ahead'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t, size = sz, data = count.t %>% filter(target2 == 'peak week'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .8, size = sz, data = count.t %>% filter(target2 == 'peak intensity'), vjust = 1) +
  geom_text(aes(x = target2,label = lab), y = ymax.t * .6, size = sz, data = count.t %>% filter(target2 == 'total'), vjust = 1, hjust = .6) +
  geom_hline(yintercept = 0, color = 'black', size = .08) +
  facet_rep_wrap(~ measure + resp.sn, ncol = 5, # , scales = 'free_y'
                 repeat.tick.labels = T) + # , labeller = label_wrap_gen(multi_line=FALSE)
  labs(x = '', y = 'Diffence in log score', color = 'State') +
  ggtitle(paste0('(', LETTERS[2], ') Baseline vs optimized setting, by virus season')) + 
  guides(color = guide_legend(override.aes = list(size = 0.3))) +
  coord_cartesian(ylim = c(ymin.t, ymax.t)) +
  theme_minimal() +theme.t6c

p = ggarrange(pp1, pp2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
pdf(paste0(dir_out,'Fig6_base.v.opt', agg.tag,'.pdf'), width = 8.5, height = 5)
print(p)
dev.off() 
################################################################################################


################################################################################################
# FIG 7
# plot the pt estimates for accuracy
agg.tag = '_acc'
# first plot it for all forecasts
tda.t = stats.acc.base.v.opt %>% filter(measure %in% c('Cases','Deaths') & loc != 'all' & target != 'all' & variant=='all' & metric == pt.met.t & is.na(wave) & is.na(resp.sn)) %>%
  mutate(base = as.numeric(mean.base) * 100, best = as.numeric(mean.opt) * 100, 
         mean.base = NULL, mean.opt = NULL,
         target = factor(target, levels = target2_vec))
tda.t2 = tda.t %>% dplyr::select(loc, variant, measure, target, base, best) %>% 
  melt(., id.var = c('loc', 'variant', 'measure', 'target')) %>%
  setnames(., 'variable', 'model')
ncol.t = tda.t %>% .$target %>% unique %>% length
pp1 = ggplot(tda.t, aes(x = loc, ymin = base, ymax = best), environment = environment()) +
  geom_linerange(aes(color = loc),size = .3, show.legend = F) + # position = position_dodge(width = 0.2), 
  facet_rep_wrap(~ measure + target, ncol = ncol.t, # , scales = 'free_y'
                 repeat.tick.labels = T) + # , labeller = label_wrap_gen(multi_line=FALSE)
  labs(x = '', y = 'Accuracy (%)', color = 'State') +
  ggtitle(paste0('(', LETTERS[1], ') All')) + 
  # coord_cartesian(ylim = c(ymin.t, ymax.t)) +
  theme_minimal() +theme.t6g

pp1 = pp1 + geom_point(data = tda.t2, mapping = aes(x = loc, y = value, shape = model), size = 1, inherit.aes = F)

# plot for resp sn only
tda.t = stats.acc.base.v.opt %>% filter(measure %in% c('Cases','Deaths') & loc != 'all' & target != 'all' &
                                       resp.sn == 'Respiratory season' & variant == 'all' & metric == pt.met.t & is.na(wave)) %>%
  mutate(base = as.numeric(mean.base) * 100, best = as.numeric(mean.opt) * 100, 
         mean.base = NULL, mean.opt = NULL,
         target = factor(target, levels = target2_vec))
tda.t2 = tda.t %>% dplyr::select(loc, variant, measure, target, base, best) %>% 
  melt(., id.var = c('loc', 'variant', 'measure', 'target')) %>%
  setnames(., 'variable', 'model')
pp2 = ggplot(tda.t, aes(x = loc, ymin = base, ymax = best), environment = environment()) +
  geom_linerange(aes(color = loc),size = .3, show.legend = F) + # position = position_dodge(width = 0.2), 
  facet_rep_wrap(~ measure + target, ncol = ncol.t, # , scales = 'free_y'
                 repeat.tick.labels = T) + # , labeller = label_wrap_gen(multi_line=FALSE)
  labs(x = '', y = 'Accuracy (%)', color = 'State') +
  ggtitle(paste0('(', LETTERS[2], ') Respiratory season only')) + 
  # coord_cartesian(ylim = c(ymin.t, ymax.t)) +
  theme_minimal() + theme.t6g

pp2 = pp2 + geom_point(data = tda.t2, mapping = aes(x = loc, y = value, shape = model), size = 1, inherit.aes = F) 


p = ggarrange(pp1, pp2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
pdf(paste0(dir_out,'Fig7_base.v.opt', agg.tag,'.pdf'), width = 9, height = 7)
print(p)
dev.off() 
################################################################################################

# FIG 8 & 9: use script "plot_figs_last.fcast.R"


################################################################################################
# output summary tables
# deflation
da.t = stats.cp.deflat %>% mutate(sig = case_when(as.numeric(pvalue.wrs) < .05 ~ '*', 
                                                   T ~ '')) %>%  
  mutate(resp.sn = resp.sn %>% as.character()) %>% 
  replace_na(list(resp.sn = 'All', wave = 'All')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
    measure = factor(measure, levels = measure_vec, labels = measure_lab),
    scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
    seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
    target = factor(target, levels = c('all', target2_vec)),
    variant = factor(variant, levels = c('all', variant_vec), labels = c('All', variant_lab)),
    State = factor(loc, levels = state_vec, labels = state_lab))  %>% 
  rowwise() %>% 
  mutate(out = paste(median.diff, sig)) %>%
  rowwise() %>% 
  mutate(out = paste0(perc.diff,'%, ', out))

tab0 = da.t  %>% filter(cp == '0.9 vs none' & metric == score.met.t & target == 'all' & variant == 'All'  & resp.sn == 'All') %>%
  dcast(., scenario + measure + State ~ seasonality, value.var = 'out')
tab1 = da.t %>% filter(cp == '0.9 vs none' & scenario == 'No new variants' & metric == score.met.t & target == 'all' & variant == 'All'  & resp.sn == 'All') %>%
  dcast(., measure + State ~ seasonality, value.var = 'out')
tab2 = da.t %>% filter(cp == '0.9 vs none' & scenario == 'New variants' & metric == score.met.t & target == 'all' & variant == 'All'  & resp.sn == 'All') %>%
  dcast(., measure + State ~ seasonality, value.var = 'out')

sheets = list('all' = tab0, "NoNewVariants" =tab1, "NewVariants" = tab2) 
File.name = paste0(dir_out,'Tab_stats_df1.v.09.xlsx')
wb <- openxlsx:: createWorkbook()
for(i in 1:length(sheets)){
  sheet.name = names(sheets)[i]
  sheet.cont = sheets[[i]]
  addWorksheet(wb, sheet.name)
  writeData(wb, i, sheet.cont, colNames = T)
  setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
}
openxlsx:: saveWorkbook(wb, File.name, overwrite = T)

# sce:
da.t = stats.cp.sce %>% mutate(sig = case_when(as.numeric(pvalue.wrs) < .05 ~ '*', 
                                                  T ~ '')) %>%  
  mutate(resp.sn = resp.sn %>% as.character()) %>% 
  replace_na(list(resp.sn = 'All', wave = 'All')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(cp = factor(cp, levels = cp.sce_vec, labels = cp.sce_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         # scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         target = factor(target, levels = c('all', target2_vec)),
         variant = factor(variant, levels = c('all', variant_vec), labels = c('All', variant_lab)),
         State = factor(loc, levels = state_vec, labels = state_lab))  %>% 
  rowwise() %>% 
  mutate(out = paste(median.diff, sig)) %>%
  rowwise() %>% 
  mutate(out = paste0(perc.diff,'%, ', out))

# tmp = tda.cp.sce %>% filter(seasonality == "y/tf" & fcast.deflat == df.opt & metric == score.met.t)
# tmp = stats.cp.sce %>% filter(seasonality == "y/tf" & fcast.deflat == 0.9)
tab1 = da.t %>% filter(fcast.deflat == df.opt & metric == score.met.t & target == 'all' & variant == 'All'& resp.sn == 'All') %>%
  dcast(., measure + State ~ seasonality, value.var = 'out')

sheets = list("CompareScenarios" =tab1) 
File.name = paste0(dir_out,'Tab_stats_sce.xlsx')
wb <- openxlsx:: createWorkbook()
for(i in 1:length(sheets)){
  sheet.name = names(sheets)[i]
  sheet.cont = sheets[[i]]
  addWorksheet(wb, sheet.name)
  writeData(wb, i, sheet.cont, colNames = T)
  setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
}
openxlsx:: saveWorkbook(wb, File.name, overwrite = T)

# seasonality
# all states combined
out.cp.sn = rbind(stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave)) %>% 
                    mutate(metric = 'Log score',
                           measure = factor(measure, levels = measure_vec, labels = measure_lab), 
                           cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                           scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                           # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                           perc.diff = paste0(perc.diff, '%')),
                  # acc
                  stats.acc.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave)) %>% 
                    mutate(metric = 'Accuracy',
                           cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                           scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                           # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                           perc.diff = paste0(perc.diff, '%')) %>%
                  mutate(median.diff = paste0(round(as.numeric(diff.median.wrs)*100, 2), '% (',round(as.numeric(diff.median.wrs.ci95lwr)*100, 2),'%, ',round(as.numeric(diff.median.wrs.ci95upr)*100, 2),'%)'))
                  ) %>%
  mutate(metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  mutate(season = resp.sn %>% as.character(), resp.sn = NULL) %>% 
  replace_na(list(season = 'All')) %>% 
  mutate(season = factor(season, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(sig = case_when(diff.median.wrs < .05 ~ '*', T ~ '')) %>%
  mutate(out = paste0(perc.diff, ', ', median.diff, sig))  %>%
  # mutate(out = paste0(perc.diff, ', ', round(as.numeric(diff.median.wrs)*100, 2), '% (',round(as.numeric(diff.median.wrs.ci95lwr)*100, 2),'%, ',round(as.numeric(diff.median.wrs.ci95upr)*100, 2),'%)', sig)) %>%
  dcast(season + metric + measure ~ cp, value.var = 'out') %>%
  mutate(wave = 'All')
out.cp.sn.offsn.bywave = rbind(rbind(stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & metric == score.met.t & target == 'all' & loc == 'all' & 
                                                              variant == 'non.Omicron' & wave == 'wave2' & resp.sn == 'Off season') %>% 
                                       mutate(time = 'pre-VOC off season'), 
                                     stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & metric == score.met.t & target == 'all' & loc == 'all' & 
                                                              (variant == 'Omicron' | wave %in% c('Alpha', 'Delta')) & resp.sn == 'Off season') %>% 
                                       mutate(time = 'VOC off season') %>% 
                                       mutate(wave = wave %>% as.character()) %>%
                                       replace_na(list(wave = 'Omicron'))) %>% 
                                 mutate(metric = 'Log score',
                                        measure = factor(measure, levels = measure_vec, labels = measure_lab), 
                                        cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                                        scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                                        # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                                        perc.diff = paste0(perc.diff, '%')),
                               # acc
                               rbind(stats.acc.cp.sn %>% filter(measure %in% c('Cases','Deaths') & fcast.deflat == df.opt & scenario == sce.opt & metric == pt.met.t & target == 'all' & loc == 'all' & 
                                                                  variant == 'non.Omicron' & wave == 'wave2' & resp.sn == 'Off season') %>% 
                                       mutate(time = 'pre-VOC off season'), 
                                     stats.acc.cp.sn %>% filter(measure %in% c('Cases','Deaths') & fcast.deflat == df.opt & scenario == sce.opt & metric == pt.met.t & target == 'all' & loc == 'all' & 
                                                                  (variant == 'Omicron' | wave %in% c('Alpha', 'Delta')) & resp.sn == 'Off season') %>% 
                                       mutate(time = 'VOC off season') %>% 
                                       mutate(wave = wave %>% as.character()) %>%
                                       replace_na(list(wave = 'Omicron'))) %>% 
                                 mutate(metric = 'Accuracy',
                                        # measure = factor(measure, levels = measure_vec, labels = measure_lab), 
                                        cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                                        scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                                        # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                                        perc.diff = paste0(perc.diff, '%'))  %>%
                               mutate(median.diff = paste0(round(as.numeric(diff.median.wrs)*100, 2), '% (',round(as.numeric(diff.median.wrs.ci95lwr)*100, 2),'%, ',round(as.numeric(diff.median.wrs.ci95upr)*100, 2),'%)'))
                              ) %>%
  mutate(wave = factor(wave, levels = wave2_vec, labels = wave2_lab),
         metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>% 
  mutate(sig = case_when(diff.median.wrs < .05 ~ '*', T ~ '')) %>%
  mutate(out = paste0(perc.diff, ', ', median.diff, sig))  %>%
  dcast(wave + metric + measure  ~ cp, value.var = 'out') %>% 
  mutate(season = 'Off season')

tab.cp.sn = rbind(out.cp.sn, out.cp.sn.offsn.bywave) %>%
  mutate(wave = factor(wave, levels = c('All', "2nd wave","Alpha", 'Delta',"Omicron")),
         season = factor(season, levels = c('All', 'Respiratory season', 'Off season'))) %>%
  setcolorder(., neworder = c('wave','season','metric','measure', "None vs fixed seasonality","None vs transformed seasonality", "Fixed vs transformed seasonality")) %>%
  arrange(wave, season)

tab.cp.sn.byloc.no.v.tf = stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == 'newV' & cp == "no.v.ytf" &
                                metric == score.met.t & target == 'all' & variant == 'all') %>% 
  mutate(sig = case_when(as.numeric(pvalue.wrs) < .05 ~ '*', 
                                               T ~ '')) %>%  
  mutate(resp.sn = resp.sn %>% as.character()) %>% 
  replace_na(list(resp.sn = 'All', wave = 'All')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         target = factor(target, levels = c('all', target2_vec)),
         variant = factor(variant, levels = c('all', variant_vec), labels = c('All', variant_lab)),
         State = factor(loc, levels = state_vec, labels = state_lab))  %>% 
  mutate(out = paste(median.diff, sig)) %>%
  mutate(out = paste0(perc.diff,'%, ', out)) %>% 
  mutate(metric = 'Log score') %>%
  dcast(., metric + measure + State ~ resp.sn, value.var = 'out')  %>%
  rbind(., 
        stats.acc.cp.sn %>% filter(fcast.deflat == df.opt & scenario == 'newV' & cp == "no.v.ytf" &
                                     metric == pt.met.t & target == 'all' & variant == 'all') %>% 
          mutate(sig = case_when(as.numeric(pvalue.wrs) < .05 ~ '*', 
                                               T ~ '')) %>%  
          mutate(resp.sn = resp.sn %>% as.character()) %>% 
          replace_na(list(resp.sn = 'All', wave = 'All')) %>% 
          mutate(resp.sn = factor(resp.sn, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
          mutate(cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                 # measure = factor(measure, levels = measure_vec, labels = measure_lab),
                 scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
                 # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
                 target = factor(target, levels = c('all', target2_vec)),
                 variant = factor(variant, levels = c('all', variant_vec), labels = c('All', variant_lab)),
                 State = factor(loc, levels = state_vec, labels = state_lab))  %>% 
          rowwise() %>% 
          mutate(out = paste0(round(as.numeric(diff.median.wrs)*100, 2), '% (',round(as.numeric(diff.median.wrs.ci95lwr)*100, 2),'%, ',round(as.numeric(diff.median.wrs.ci95upr)*100, 2),'%)', sig)) %>%
          mutate(out = paste0(perc.diff,'%, ', out)) %>% 
          mutate(metric = 'Accuracy')  %>%
        dcast(., metric + measure + State ~ resp.sn, value.var = 'out') 
        ) %>%
  mutate(metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  setcolorder(., neworder = c('State','metric','measure')) %>%
  arrange(State, metric, measure)

tab.cp.sn.by.loc = rbind(stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & metric == score.met.t & target == 'all' & loc != 'all' & variant == 'all' & is.na(wave)) %>% 
               mutate(metric = 'Log score',
                      measure = factor(measure, levels = measure_vec, labels = measure_lab), 
                      cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                      scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                      # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                      perc.diff = paste0(perc.diff, '%')),
             # acc
             stats.acc.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc != 'all' & variant == 'all' & is.na(wave)) %>% 
               mutate(metric = 'Accuracy',
                      cp = factor(cp, levels = cp.sn_vec, labels = cp.sn_lab), 
                      scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                      # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                      perc.diff = paste0(perc.diff, '%')) %>%
               mutate(median.diff = paste0(round(as.numeric(diff.median.wrs)*100, 2), '% (',round(as.numeric(diff.median.wrs.ci95lwr)*100, 2),'%, ',round(as.numeric(diff.median.wrs.ci95upr)*100, 2),'%)'))
) %>%
  mutate(metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  mutate(season = resp.sn %>% as.character(), resp.sn = NULL) %>% 
  replace_na(list(season = 'All')) %>% 
  mutate(season = factor(season, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(sig = case_when(diff.median.wrs < .05 ~ '*', T ~ '')) %>%
  mutate(out = paste0(perc.diff, ', ', median.diff, sig))  %>%
  # mutate(out = paste0(perc.diff, ', ', round(as.numeric(diff.median.wrs)*100, 2), '% (',round(as.numeric(diff.median.wrs.ci95lwr)*100, 2),'%, ',round(as.numeric(diff.median.wrs.ci95upr)*100, 2),'%)', sig)) %>%
  dcast(loc + season + metric + measure ~ cp, value.var = 'out')


sheets = list("All locations" =tab.cp.sn, "By location" = tab.cp.sn.by.loc) 
File.name = paste0(dir_out,'Tab_S3_S4_stats_sn.xlsx')
wb <- openxlsx:: createWorkbook()
for(i in 1:length(sheets)){
  sheet.name = names(sheets)[i]
  sheet.cont = sheets[[i]]
  addWorksheet(wb, sheet.name)
  writeData(wb, i, sheet.cont, colNames = T)
  setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
}
openxlsx:: saveWorkbook(wb, File.name, overwrite = T)


# baseline vs. best-performing
da.t = stats.base.v.opt %>% mutate(sig = case_when(as.numeric(pvalue.wrs) < .05 ~ '*', 
                                                   T ~ '')) %>%  
  mutate(resp.sn = resp.sn %>% as.character()) %>% 
  replace_na(list(resp.sn = 'All', wave = 'All')) %>% 
  mutate(resp.sn = factor(resp.sn, levels = c('All', 'Respiratory season', 'Off season'))) %>% 
  mutate(# cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
         measure = factor(measure, levels = measure_vec, labels = measure_lab),
         # scenario = factor(scenario, levels = sce_vec, labels = sce_lab),
         # seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab),
         target = factor(target, levels = c('all', target2_vec)),
         variant = factor(variant, levels = c('all', variant_vec), labels = c('All', variant_lab)),
         State = factor(loc, levels = state_vec, labels = state_lab))  %>% 
  rowwise() %>% 
  mutate(out = paste(median.diff, sig)) %>%
  rowwise() %>% 
  mutate(out = paste0(perc.diff,'%, ', out))
tab1 = da.t %>% filter(metric == score.met.t & target == 'all' & variant == 'All') %>%
  dcast(., State + measure  ~ resp.sn, value.var = 'out')

tab2 = da.t %>% filter(metric == score.met.t & target == 'all' & resp.sn == 'All' & wave == 'All') %>%
  dcast(., State + measure ~ variant, value.var = 'out')

tab0 = tab2 %>% left_join(x = ., y = tab1 %>% dplyr::select(-All), by = c('measure', 'State'))

sheets = list('Tabe 1' = tab0, "byRespSn" =tab1, "byPeriod" = tab2) 
File.name = paste0(dir_out,'Table1_stats_base.v.opt.xlsx')
wb <- openxlsx:: createWorkbook()
for(i in 1:length(sheets)){
  sheet.name = names(sheets)[i]
  sheet.cont = sheets[[i]]
  addWorksheet(wb, sheet.name)
  writeData(wb, i, sheet.cont, colNames = T)
  setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
}
openxlsx:: saveWorkbook(wb, File.name, overwrite = T)

# accuracy
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
  mutate(base = (as.numeric(mean.base) * 100) %>% round(0) %>% paste0('%'),
         opt = (as.numeric(mean.opt) * 100) %>% round(0) %>% paste0('%')) %>%
  mutate(out = paste(base,'v', opt, sig))

tab1 = da.t %>% filter(resp.sn == 'All' & measure %in% c('Cases','Deaths') & metric == pt.met.t & target != 'all' & variant == 'All') %>%
  dplyr::select(State, resp.sn, measure, target, out) %>%
  dcast(., State + measure ~ target, value.var = 'out')

tab2 = da.t %>% filter(resp.sn == 'Respiratory season' & measure %in% c('Cases','Deaths') & metric == pt.met.t & target != 'all' & variant == 'All') %>%
  dplyr::select(State, resp.sn, measure, target, out) %>%
  dcast(., State + measure ~ target, value.var = 'out')

tab3 = da.t %>% filter(resp.sn == 'Off season' & measure %in% c('Cases','Deaths') & metric == pt.met.t & target != 'all' & variant == 'All') %>%
  dplyr::select(State, resp.sn, measure, target, out) %>%
  dcast(., State + measure ~ target, value.var = 'out')

sheets = list("Table 2" =tab1, "Respiratory season" = tab2, "Off season" = tab3) 
File.name = paste0(dir_out,'Table2_acc_base.v.opt.xlsx')
wb <- openxlsx:: createWorkbook()
for(i in 1:length(sheets)){
  sheet.name = names(sheets)[i]
  sheet.cont = sheets[[i]]
  addWorksheet(wb, sheet.name)
  writeData(wb, i, sheet.cont, colNames = T)
  setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
}
openxlsx:: saveWorkbook(wb, File.name, overwrite = T)


# check stats
out1s = stats.cp.deflat %>% filter(cp == 'deflat1.v.095' & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out1s

out2s = stats.cp.deflat %>% filter(cp == 'deflat1.v.09' & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out2s

out3s = stats.cp.deflat %>% filter(cp == 'deflat095.v.09' & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out3s

out1a = stats.acc.cp.deflat %>% filter(measure %in% c('Cases','Deaths') & cp == 'deflat1.v.095' & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out1a

out2a = stats.acc.cp.deflat %>% filter(measure %in% c('Cases','Deaths') & cp == 'deflat1.v.09' & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out2a


out3a = stats.acc.cp.deflat %>% filter(measure %in% c('Cases','Deaths') & cp == 'deflat1.v.095' & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out3a

out.cp.df = rbind(stats.cp.deflat %>% filter(metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  mutate(metric = 'Log score',
         measure = factor(measure, levels = measure_vec, labels = measure_lab), 
         cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
         scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
         seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
         perc.diff = paste0(perc.diff, '%')),
  # acc
  stats.acc.cp.deflat %>% filter(measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
    mutate(metric = 'Accuracy',
           cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
           scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
           seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
           perc.diff = paste0(perc.diff, '%'))) %>%
  mutate(metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  dcast(metric + measure + scenario + seasonality ~ cp, value.var = 'perc.diff')

out.cp.sce = rbind(stats.cp.sce %>% filter(fcast.deflat == df.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
                    mutate(metric = 'Log score',
                           measure = factor(measure, levels = measure_vec, labels = measure_lab), 
                           # cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
                           # scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                           seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                           perc.diff = paste0(perc.diff, '%')),
                  # acc
                  stats.acc.cp.sce %>% filter(fcast.deflat == df.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
                    mutate(metric = 'Accuracy',
                           # cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
                           # scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                           seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                           perc.diff = paste0(perc.diff, '%'))) %>%
  mutate(metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  dcast(metric + measure ~ seasonality, value.var = 'perc.diff')


out.cp.sce.bywave = rbind(rbind(stats.cp.sce %>% filter(fcast.deflat == df.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'non.Omicron' & !is.na(wave) & is.na(resp.sn)),
                                stats.cp.sce %>% filter(fcast.deflat == df.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'Omicron' & is.na(wave) & is.na(resp.sn)) %>%
                                  mutate(wave = 'Omicron')) %>% 
                     mutate(metric = 'Log score',
                            wave = factor(wave, levels = wave2_vec, labels = wave2_lab),
                            measure = factor(measure, levels = measure_vec, labels = measure_lab), 
                            # cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
                            # scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                            seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                            perc.diff = paste0(perc.diff, '%')),
                   # acc
                   rbind(stats.acc.cp.sce %>% filter(fcast.deflat == df.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'non.Omicron' & !is.na(wave) & is.na(resp.sn)),
                         stats.acc.cp.sce %>% filter(fcast.deflat == df.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'Omicron' & is.na(wave) & is.na(resp.sn)) %>%
                           mutate(wave = 'Omicron')) %>% 
                     mutate(metric = 'Accuracy',
                            wave = factor(wave, levels = wave2_vec, labels = wave2_lab),
                            # cp = factor(cp, levels = cp.df_vec, labels = cp.df_lab), 
                            # scenario = factor(scenario, levels = sce_vec, labels = sce_lab), 
                            seasonality = factor(seasonality, levels = sn_vec, labels = sn_lab), 
                            perc.diff = paste0(perc.diff, '%'))) %>%
  mutate(metric = factor(metric, levels = c('Log score', 'Accuracy'))) %>%
  dcast(metric + measure + seasonality ~ wave, value.var = 'perc.diff')


sheets = list("Table S1" = out.cp.df, "Table S2" = out.cp.sce.bywave) 
File.name = paste0(dir_out,'Table_S1_S2.xlsx')
wb <- openxlsx:: createWorkbook()
for(i in 1:length(sheets)){
  sheet.name = names(sheets)[i]
  sheet.cont = sheets[[i]]
  addWorksheet(wb, sheet.name)
  writeData(wb, i, sheet.cont, colNames = T)
  setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
}
openxlsx:: saveWorkbook(wb, File.name, overwrite = T)


  

out4s = stats.cp.sce %>% filter(fcast.deflat == df.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out4s
out4a = stats.acc.cp.sce %>% filter(fcast.deflat == df.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out4a

stats.cp.sce %>% filter(measure == 'death' & fcast.deflat == df.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'Omicron' & is.na(resp.sn)) 

out5s = stats.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out5s
out5a = stats.acc.cp.sn %>% filter(fcast.deflat == df.opt & scenario == sce.opt & measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) %>%
  summarise(perc.diff.min = min(perc.diff %>% as.numeric), perc.diff.max = max(perc.diff %>% as.numeric))
out5a

# combine all
out6s = stats.base.v.opt %>% filter(metric == score.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) 
out6s

out6a = stats.acc.base.v.opt %>% filter(measure %in% c('Cases','Deaths') & metric == pt.met.t & target == 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure) 
out6a

out7a = stats.acc.base.v.opt %>% filter(measure %in% c('Cases','Deaths') & metric == pt.met.t & target != 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  group_by(measure, target) 
out7a

tmp = stats.base.v.opt %>% filter(metric == score.met.t  & target != 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
  dplyr::select(loc, variant, wave, resp.sn, measure, target, n.smp) %>% mutate(measure = factor(measure, levels = measure_vec, labels = measure_lab))  %>%
  full_join(., y = stats.acc.base.v.opt %>% filter(measure %in% c('Cases','Deaths') & metric == pt.met.t & target != 'all' & loc == 'all' & variant == 'all' & is.na(wave) & is.na(resp.sn)) %>% 
              dplyr::select(loc, variant, wave, resp.sn, measure, target, n.smp),
              by = c('loc', 'variant', 'wave', 'resp.sn', 'measure', 'target'),
            suffix = c('.s','.a'))


# SUPPLEMENTAL TABLES
# table S8 parameter bounds
{
  locs = c('California', 'Florida', 'Iowa', 'Massachusetts', 'Michigan', 
           'New York', 'Pennsylvania', 'Texas', 'Washington', 'Wyoming')
  PRAMBOUNDS = read_excel(paste0(dir_code, 'parm.bounds_pop.den_lowerbeta.xlsx'), sheet = 1) %>% data.table() # adjust beta per pop density slightly 
  PRAMBOUNDS$date.start = PRAMBOUNDS$date.start %>% as.Date()
  PRAMBOUNDS$date.end = PRAMBOUNDS$date.end %>% as.Date()
  PRAMBOUNDS$lwr = PRAMBOUNDS$lwr %>% as.numeric()
  PRAMBOUNDS$upr = PRAMBOUNDS$upr %>% as.numeric()
  parms = PRAMBOUNDS %>% filter(location %in% locs)
  
  # estimated from VE data
  est_parmVimmLoss = read.csv(paste0(dir_code, 'est_parmVimmLoss.csv')) %>% data.table()
  
  # summer2021 and Delta wave may have multiple settings depending on situations
  parms_edit = parms
  for(loc.t in locs){
    for(parm.t in c('alpha','beta')){
      for(pp in c('summer2021','Delta')){
        tmp = parms %>% filter(location == loc.t & type == 'SR' & parm == parm.t & period == pp) %>% arrange(., date.start)
        if(nrow(tmp) > 1){
          if(parm.t == 'alpha'){
            tmp$period = paste0(pp, 1:nrow(tmp))
          } else if(parm.t == 'beta') {
            # identify superspread
            tmp = tmp %>% mutate(period = case_when(lwr == max(lwr) & upr == max(upr) ~ paste0(period,'_holiday'),
                                                    T ~ period))
            tmp_h = tmp %>% filter(period == paste0(pp,'_holiday'))
            tmp_nh = tmp %>% filter(period != paste0(pp,'_holiday')) %>% arrange(., date.start) 
            tmp_nh = tmp_nh %>%
              mutate(period = paste0(period,1:nrow(tmp_nh)))
            tmp = rbind(tmp_h, tmp_nh)
          }
          
          
          parms_edit = rbind(parms_edit %>% filter(!(location == loc.t & type == 'SR' & parm == parm.t & period == pp)),
                             tmp)
        }
        
      }
    }
  }
  
  
  p.mob_bounds = c(.5, 1.5); # scaling for mobility
  Td.mean_bounds =  c(5,8) # mean Td: reporting delay
  Td.mean_SRbounds = c(5, 7)
  Td.sd_bounds = c(1,3) # Td, sd: reporting delay sd 
  
  # imm_bounds = c(2, 3) * 365
  if(variant.tag == 'non.Omicron'){
    imm_bounds = c(2, 3) * 365
  } else if (variant.tag == 'Omicron'){
    imm_bounds = c(1, 3) * 365
  }
  
  iniSomicron = c(.5, .9)
  seed_max = (da.t$case[1] + .1) * 100 /20 
  
  VE1wt = .85
  VE2wt = .95  # higher b/c we are using mortality data too
  
  VE1delta = .5
  VE2delta = .8  # higher b/c we are using mortality data too
  
  VE1omicron = .1 # .5  # 1st dose
  VE2omicron = .7  # 2nd and 3rd dose
  init0 = rbind(data.table(Parameter = 'Initial susceptible', Symbol = 'S(t=0)', Range = paste0('All locations: non-Omicron period, U[99%, 100%] population; Omicron period, U[50%, 90%] population'), Note = 'n/a'),
                data.table(Parameter = 'Initial exposed', Symbol = 'E(t=0)', Range = paste0('All locations: U[5, 50] * no. cases during 1st week'), Note = 'n/a'),
                data.table(Parameter = 'Initial infectious', Symbol = 'I(t=0)', Range = paste0('All locations: U[5,50] * no. cases during 1st week'), Note = 'n/a'),
                data.table(Parameter = 'Infectious period', Symbol = 'D', Range = paste0('All locations: U[2, 5] days'), Note = 'n/a'),
                data.table(Parameter = 'Latency period', Symbol = 'Z', Range = paste0('All locations: U[2, 5] days'), Note = 'n/a'),
                data.table(Parameter = 'Duration of immunity (from prior infection)', Symbol = 'L', Range = paste0('All locations: non-Omicron period, U[2, 3] years; Omicron period: U[1, 3] years'), Note = 'n/a'),
                data.table(Parameter = 'Time-to-detection, mean', Symbol = 'Td, mean', Range = paste0('All locations: U[5, 8] days'), Note = 'n/a'),
                data.table(Parameter = 'Time-to-detection, sd', Symbol = 'Td, sd', Range = paste0('All locations: U[1, 3] days'), Note = 'To allow variation in time to diagnosis/reporting'),
                data.table(Parameter = 'Scaling of NPI effectiveness', Symbol = 'e', Range = paste0('All locations: U[0.5, 1.5]'), Note = 'Around 1, with a large bound to be flexible'),
                data.table(Parameter = 'Vaccine efficacy (VE)', Symbol = 'n/a', Range = 'All locations: before Delta, VE1=85%, VE2 = 95%; Delta, VE1 = 50%, VE2 = 80%; Omicron, VE1 = 10%, VE2 (combined 2nd and 3rd doses) = 70%', Note = 'Used higher VE values, as the observations included both cases/infections and deaths; i.e., here VE is for both infections and mortality'),
                data.table(Parameter = 'VE waning', Symbol = 'rho', 
                           Range = paste0('rho(t) = 1/(1+exp(-k * (t - tm.imm/2); for wildtype: k = ', round(est_parmVimmLoss[variant=='wt']$k, 3), '; tm.imm = ', round(est_parmVimmLoss[variant=='wt']$tm.imm, 0), 
                                          '; for Delta: k = ', round(est_parmVimmLoss[variant=='delta']$k, 3), '; tm.imm = ', round(est_parmVimmLoss[variant=='delta']$tm.imm, 0),
                                          '; for Omicron: k = ', round(est_parmVimmLoss[variant=='omicron']$k, 3), '; tm.imm = ', round(est_parmVimmLoss[variant=='omicron']$tm.imm, 0)),
                           Note = 'Parameter in the logistic function fitted based on data from UKHSA')
  )
  init1 = parms_edit %>% filter(type == 'initialization' & parm!='beta') %>% dplyr::select(parm, lwr, upr) %>% unique %>%
    mutate(Range = paste0('all locations: U [', lwr, ', ', upr,']'))
  init2 = parms_edit %>% filter(type == 'initialization' & parm =='beta') %>% dplyr::select(location, parm, lwr, upr) %>% unique %>%
    mutate(Range = paste0('U [', lwr, ', ', upr,']')) %>% group_by(Range) %>%
    summarise(location_all = paste(location,collapse = ', ')) %>%
    mutate(Range = paste0(location_all,': ', Range))
  init = rbind(init1 %>% dplyr::select(parm, Range), data.table(parm = 'beta', Range = paste(init2$Range, collapse = '; '))) %>%
    setnames(.,'parm','Symbol') %>%
    mutate(Parameter = factor(Symbol, levels = c('beta','alpha','ifr'), labels = c('Transmission rate', 'Infection-detection rate','Infection-fatality risk')),
           Note = 'n/a') %>% 
    rbind(init0, .)
  init$Type = 'Initialization'
  init$Parameter = factor(init$Parameter, levels = c('Initial susceptible', 'Initial exposed','Initial infectious','Transmission rate',
                                                     'Infectious period','Latency period','Duration of immunity (from prior infection)',
                                                     'Infection-detection rate','Time-to-detection, mean','Time-to-detection, sd',
                                                     'Infection-fatality risk', 'Scaling of NPI effectiveness', 'Vaccine efficacy (VE)', 'VE waning'))
  init = init %>% arrange(., by = 'Parameter')
  # SR Range, for alpha
  parm_sr = NULL
  for(parm.t in c('beta', 'alpha', 'ifr')){
    parm.sr.t  = parms_edit %>% filter(type == 'SR' & parm == parm.t & period != 'start') %>% dplyr::select(date.start, date.end, parm, lwr, upr, period, location) %>% unique #  %>%
    # mutate(period = paste0(date.start, ' to ', date.end, ' (', period, ')'))
    # mutate(period = case_when(grepl('Omicron', period) ~ paste0(date.start, ' to ', date.end, ' (Omicron)'),
    #                         T ~ paste0(date.start, ' to ', date.end, ' (', period, ')')))
    
    tmp_periods = parm.sr.t   %>% .$period %>% unique() %>% sort
    sr.t = NULL
    for(pp in tmp_periods){
      tmp = parm.sr.t %>% filter(parm == parm.t & period == pp) %>% dplyr::select(date.start, date.end, lwr, upr, location) %>%
        mutate(Range = paste0('U [', lwr, ', ', upr,']')) %>% group_by(Range)
      if((tmp$location %>% unique() %>% length()) == length(locs) & length(unique(tmp$Range))==1){
        # same for all
        date_Range = paste(min(tmp$date.start), 'to', max(tmp$date.end))
        tmp = data.table(parm = parm.t, Range = paste0(pp, ' (', date_Range, '): all locations, ', tmp$Range[1]))
      } else {
        # multiple settings
        date_Range = paste(min(tmp$date.start), 'to', max(tmp$date.end))
        date_Range = paste(min(tmp$date.start), 'to', max(tmp$date.end))
        tmp = tmp %>% summarise(location_all = paste(location,collapse = ', ')) %>%
          mutate(Range = paste0(location_all,', ', Range))
        tmp =  data.table(parm = parm.t, Range = paste0(pp, ' (', date_Range, '): ', paste(tmp$Range, collapse = '; ')))
      }
      sr.t = rbind(sr.t, tmp)
    }
    parm_sr = rbind(parm_sr, data.table(Type = 'SR', parm = parm.t, Range = paste(sr.t$Range,collapse = '; \n')))
  }
  parm_bounds = rbind(init, 
                      parm_sr %>% setnames(.,'parm','Symbol') %>%
                        mutate(Parameter = factor(Symbol, levels = c('beta','alpha','ifr'), labels = c('Transmission rate', 'Infection-detection rate','Infection-fatality risk')),
                               Note = 'n/a')
  ) %>% setcolorder(., neworder = c('Type','Parameter','Symbol','Range','Note'))
  write_xlsx(parm_bounds, paste0(dir_out,'Table_S8.xlsx'))
  sheets = list("Table S8" = parm_bounds) 
  File.name = paste0(dir_out,'Table_S8_parm_bounds.xlsx')
  wb <- openxlsx:: createWorkbook()
  for(i in 1:length(sheets)){
    sheet.name = names(sheets)[i]
    sheet.cont = sheets[[i]]
    addWorksheet(wb, sheet.name)
    writeData(wb, i, sheet.cont, colNames = T)
    setColWidths(wb, sheet = i, cols = 1:ncol(sheet.cont), widths = 'auto')
  }
  openxlsx:: saveWorkbook(wb, File.name, overwrite = T)
}

# SUPPLEMENTAL FIGURES (some may not be included b/c need to load extra model outputs that are large in size)
# plot example of seasonality trends
source(paste0(dir_code,'get_relR0.R'))
# read seasonal trend
DAT.SN = read.csv(paste0(dir_data, 'est.sn.2000t2020_us.csv')) %>% data.table()
DAT.SNraw = read.csv(paste0(dir_data, 'est.raw.sn.2000t2020_us.csv')) %>% data.table()
# best parameter settings for the transformed seasonality
sn.tran.parms = read.csv(paste0(dir_code,'best_sn.tran.parms.csv')) %>% data.table()
num_ens = 100

locs = c('California', 'Florida', 'Iowa', 'Massachusetts', 'Michigan', 
         'New York', 'Pennsylvania', 'Texas', 'Washington', 'Wyoming')
pdf(paste0(dir_out,'FigS1_sn.trends.pdf'), width = 7, height = 10)
par(mfrow = c(5, 2), mar = c(2.5, 2.5, 1.5, .5), cex = .8, cex.lab = .95, cex.axis = .9, mgp = c(1.2, .3, 0), tck = -.02)
for(loc.t in locs){  
  # get the best parm settings for the transformed seasonality
  sn.tran.parm = sn.tran.parms %>% dplyr::filter(loc == loc.t)
  sn.relR0min.adj = sn.tran.parm$sn.relR0min.adj
  sn.dur0 = sn.tran.parm$sn.dur0;
  sn.pshift = sn.tran.parm$sn.pshift;
  
  
  # seasonality
  relR0 = DAT.SN[state == loc.t] %>% .[order(week)]
  relR0 = matrix(relR0$value, nrow = nrow(relR0), ncol = num_ens)
  
  # transformed seasonal trend
  relR0tran = fn_getR0trans(relR0=DAT.SN[state == loc.t] %>% .[order(week)] %>% .$value, # 8/10/22 used the capped relR0 so that it won't be too extreme
                            num_ens = num_ens,
                            relR0uncapped = DAT.SNraw[state == loc.t] %>% .[order(week)] %>% .$value,
                            bound.relR0min = sn.relR0min.adj+c(.8, 1.2),
                            bound.dur0 = sn.dur0 + c(-4, 4),
                            bound.pshift = sn.pshift + c(-2, 2))
  # relR0tran = fn_getR0trans(relR0=DAT.SN[state == loc.t] %>% .[order(week)] %>% .$value, num_ens)
  
  matplot(relR0tran, type = 'l', ylim = c(min(c(relR0, relR0tran)), max(c(relR0, relR0tran))), 
          xlab = 'Week of the year', ylab = expression('Relative infection risk, b'[t]),
          col = 'grey', lwd = .5)
  lines(rowMeans(relR0tran), col = 'black', lwd =1.5)
  lines(relR0[,1], col = 'blue', lwd = 1.5)
  abline(h=1, lwd = .5, col = 'red')
  legend('top', legend = c('Fixed','Transformed (ensemble)', 'Transformed (mean)'),
         col = c('blue','grey','black'), lty = 1, bty = 'n', cex = .8)
  mtext(loc.t, side = 3, outer = F, line = .1, adj = 0, cex = .8, font = 2)
}
dev.off()



