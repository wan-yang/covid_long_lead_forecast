theme.t = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=8,angle = 30),
                plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=8), legend.text=element_text(size=8),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                legend.key.size = unit(.2, 'cm'), #change legend key size
                legend.key.height = unit(.5, 'cm'), #change legend key height
                legend.key.width = unit(.2, 'cm')) #change legend key width)
theme.t2 = theme(plot.title = element_text(v=0, size = 10, margin=margin(2,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=8,angle = 45, hjust = .9),
                 plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=8), legend.text=element_text(size=8),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.2, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)
theme.t3 = theme(plot.title = element_text(v=0, size = 10, margin=margin(2,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=8,angle = 0, hjust = .9),
                 plot.margin=unit(c(c(.3, 1, -.7, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=8), legend.text=element_text(size=8),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.2, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t4 = theme(plot.title = element_text(v=0, size = 10, margin=margin(2,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=8,angle = 90, hjust = 1),
                 plot.margin=unit(c(c(.3, 1, -.7, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=8), legend.text=element_text(size=8),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.2, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t5 = theme(plot.title = element_text(v=0, size = 10, margin=margin(2,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=8,angle = 0, hjust = .5),
                 plot.margin=unit(c(c(.3, 1, 0, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=8), legend.text=element_text(size=8),
                 legend.margin=margin(0,0,0,3),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.2, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)
theme.t6 = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=8,angle = 30),
                plot.margin=unit(c(c(.3, 1, -1, .5)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=8), legend.text=element_text(size=8),
                legend.margin=margin(0,0,0,3),
                legend.box.margin=margin(-10,-10,-10,-10),
                legend.key.size = unit(.2, 'cm'), #change legend key size
                legend.key.height = unit(.5, 'cm'), #change legend key height
                legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6b = theme(plot.title = element_text(v=0, size = 12, margin=margin(0,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 11, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =11, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=11, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=11,angle = 30),
                 plot.margin=unit(c(c(.3, 1, -1, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=11), legend.text=element_text(size=11),
                 legend.margin=margin(0,0,0,3),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.5, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6c = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=8,angle = 30),
                 plot.margin=unit(c(c(.3, 1, -1.5, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=8), legend.text=element_text(size=8),
                 legend.margin=margin(0,12,0,0),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.5, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)
theme.t6d = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 15),
                  plot.margin=unit(c(c(.3, 1, -1, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,12,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6e = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                  axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 0),
                  plot.margin=unit(c(c(.3, 1, 0, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,12,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6f = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(-.3, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 25, hjust = .9),
                  plot.margin=unit(c(c(.3, 1, -1, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,12,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6g = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(-.3, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 40, hjust = .9),
                  plot.margin=unit(c(c(.3, 1, -.8, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,12,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6h = theme(plot.title = element_text(v=0, size = 9, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(-.3, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 8, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 25, hjust = .9),
                  plot.margin=unit(c(c(.3, 1, -1, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,12,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6i = theme(plot.title = element_text(v=0, size = 9, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(0, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 8, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =8, margin=margin(0,0,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0,0,0)), 
                  axis.text.x = element_text(size=8, margin=margin(-5,0,0,0),angle = 45, hjust = 1, vjust = 1),
                  plot.margin=unit(c(c(.3, 1, 0, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8, margin=margin(0,0,1,1)), legend.text=element_text(size=8, margin=margin(0,0,0,0)),
                  legend.margin=margin(1,12,0,3),# top, right, bottom, left
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)


theme.t6k = theme(plot.title = element_text(v=0, size = 9, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(0, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 8, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =8, margin=margin(0,0,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0,0,0)), 
                  axis.text.x = element_text(size=8, margin=margin(-5,0,0,0),angle = 45, hjust = 1, vjust = 1),
                  plot.margin=unit(c(c(.3, 1, 0, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8, margin=margin(0,0,1,1)), legend.text=element_text(size=8, margin=margin(0,0,0,0)),
                  legend.margin=margin(1,12,0,3),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t6l = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(0, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(5,0,1,0)),
                  axis.title = element_text(size =8, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 40, hjust = .9),
                  plot.margin=unit(c(c(0, 1, 0.1, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,0,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t7 = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                 axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                 axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=8,angle = 45, hjust = .9),
                 plot.margin=unit(c(c(.3, 1, -1, -.5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=8), legend.text=element_text(size=8),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-10,-10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.5, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t8a = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                  axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 30),
                  plot.margin=unit(c(c(0, 1.3, -3, -.3)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,0,0,5),# top, right, bottom, left
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)
theme.t8b = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  panel.spacing.y = unit(-.3, "lines"), # panel.margin = unit(.2, "lines"),
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(0,0,0,0)),
                  axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 40, hjust = .9),
                  plot.margin=unit(c(c(.5, 1.3, 0, .7)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,0,0,0),
                  legend.box.margin=margin(0,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.5, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)


getPlotHeatmap = function(tda, title.t, fill.name.t,
                         ncol.t = 6, theme.tt = theme.t){
  # only show results with at least 5 samples
  tda = tda %>% filter(as.numeric(n.smp) >=10) %>% data.table()
  
  if(max(tda$value, na.rm = T) > 100 & nrow(tda[value > 100]) >= 2){
    # use gradient
    
    cut.pt = 100
    n.neg = 2
    n.pos1 = 4
    n.pos2 = 1
    tda$value2 = ''
    tda[value > cut.pt]$value2 = cut_interval(tda[value >cut.pt]$value, n.pos2) # %>% as.character()
    tda[value >= 0 & value <=cut.pt]$value2 = cut_interval(tda[value >= 0 & value <=cut.pt]$value, n.pos1) # %>% as.character()
    tda[value < 0]$value2 = cut_interval(tda[value < 0]$value, n.neg) # %>% as.character()
    
    breaks.t = c(ggplot2:::breaks(c(min(tda[value < 0]$value), 0),"width",n = n.neg),
      ggplot2:::breaks(c(0,max(tda[value >= 0 & value <=cut.pt]$value)),"width",n = n.pos1),
      ggplot2:::breaks(c(cut.pt,max(tda[value >cut.pt]$value)),"width",n = n.pos2)
    ) %>% unique %>% sort
    
    cols.levels = c(levels(cut_interval(tda[value < 0]$value, n.neg)), 
                    levels(cut_interval(tda[value >= 0 & value <=cut.pt]$value, n.pos1)),
                    levels(cut_interval(tda[value >cut.pt]$value, n.pos2))) %>% rev
    
    tda = tda %>% mutate(value2 = factor(value2, levels = cols.levels))
    
    # use gradient
    # steps.pos <- c("white", 'gold', 'orange', 'red4')
    # steps.neg <- c('blue', 'white')
    # pal.pos=colorRampPalette(steps.pos,space="rgb")
    # pal.neg=colorRampPalette(steps.neg,space="rgb")
    # cols.pos <- pal.pos(n.pos1)
    # cols.neg <- pal.neg(8)
    # cols.t = c(tail(cols.neg, 2), cols.pos, 'brown') 
    
    
    cols.pos = brewer.pal(c(n.pos1*2+n.pos2), "YlOrRd") 
    cols.pos = c(cols.pos[seq(1, c(n.pos1)*2, by = 2)], tail(cols.pos, n.pos2)) # c(head(cols.pos, n.pos1), )
    cols.t = c(head(brewer.pal(c(n.pos1+n.pos2), "Blues"), n.neg) %>% rev, cols.pos)  %>% rev
    
    p =  ggplot(tda, aes(x = loc, y = target, fill = value2)) +
      geom_tile() + 
      facet_wrap(~ f1 + f2 + f3, ncol = ncol.t) +
      # scale_fill_distiller(palette = "RdPu", direction = direction.t) +
      scale_fill_manual(values = cols.t, na.value = alpha('grey', 50)) +
      ggtitle(title.t) + labs(x = '', y = '', fill = fill.name.t) + 
      theme_minimal() + theme.tt
  } else if(min(tda$value, na.rm = T) < -100 & nrow(tda[value < -100]) >= 2){
    cut.pt = -100
    n.pos = 4
    n.neg1 = 4
    n.neg2 = 1
    tda$value2 = ''
    tda[value < cut.pt]$value2 = cut_interval(tda[value <cut.pt]$value, n.neg2) # %>% as.character()
    tda[value >= cut.pt & value <0]$value2 = cut_interval(tda[value >= cut.pt & value <0]$value, n.neg1) # %>% as.character()
    tda[value >= 0]$value2 = cut_interval(tda[value >= 0]$value, n.pos) # %>% as.character()
    
    
    cols.levels = c(levels(cut_interval(tda[value < cut.pt]$value, n.neg2)), 
                    levels(cut_interval(tda[value >= cut.pt & value <0]$value, n.neg1)),
                    levels(cut_interval(tda[value >= 0]$value, n.pos))) %>% rev
    
    tda = tda %>% mutate(value2 = factor(value2, levels = cols.levels))
    
    cols.neg = brewer.pal(c(n.neg1*2+n.neg2), "Blues") 
    cols.neg = c(tail(cols.neg, n.neg2), cols.neg[seq(1, c(n.neg1)*2, by = 2)] %>% rev)  # c(head(cols.pos, n.pos1), )
    cols.t = c(cols.neg, head(brewer.pal(c(n.neg1+n.neg2), "YlOrRd"), n.pos)) %>% rev
    
    p =  ggplot(tda, aes(x = loc, y = target, fill = value2)) +
      geom_tile() + 
      facet_wrap(~ f1 + f2 + f3, ncol = ncol.t) +
      # scale_fill_distiller(palette = "RdPu", direction = direction.t) +
      scale_fill_manual(values = cols.t, na.value = alpha('grey', 50)) +
      ggtitle(title.t) + labs(x = '', y = '', fill = fill.name.t) + 
      theme_minimal() + theme.tt
  } else {
    p = ggplot(tda, aes(x = loc, y = target, fill = value)) +
      geom_tile() + 
      facet_wrap(~ f1 + f2 + f3, ncol = ncol.t) +
      # scale_fill_distiller(palette = "RdPu", direction = direction.t) +
      # scale_fill_gradient2(low = muted('blue'), high = muted('red')) +
      scale_fill_gradient2(low = 'blue', high = 'red', na.value = alpha('grey', 50)) +
      ggtitle(title.t) + labs(x = '', y = '', fill = fill.name.t) + 
      theme_minimal() + theme.tt
  }
  
  
  
  p
}

getPlotHeatmap.allTargets = function(tda, title.t, fill.name.t,
                          ncol.t = 6, theme.tt = theme.t){
  # only show results with at least 5 samples
  tda = tda %>% filter(as.numeric(n.smp) >=10) %>% data.table()
  
  if(max(tda$value, na.rm = T) > 100 & nrow(tda[value > 100]) >= 2){
    # use gradient
    
    cut.pt = 100
    n.neg = 2
    n.pos1 = 4
    n.pos2 = 1
    tda$value2 = ''
    tda[value > cut.pt]$value2 = cut_interval(tda[value >cut.pt]$value, n.pos2) # %>% as.character()
    tda[value >= 0 & value <=cut.pt]$value2 = cut_interval(tda[value >= 0 & value <=cut.pt]$value, n.pos1) # %>% as.character()
    tda[value < 0]$value2 = cut_interval(tda[value < 0]$value, n.neg) # %>% as.character()
    
    breaks.t = c(ggplot2:::breaks(c(min(tda[value < 0]$value), 0),"width",n = n.neg),
                 ggplot2:::breaks(c(0,max(tda[value >= 0 & value <=cut.pt]$value)),"width",n = n.pos1),
                 ggplot2:::breaks(c(cut.pt,max(tda[value >cut.pt]$value)),"width",n = n.pos2)
    ) %>% unique %>% sort
    
    cols.levels = c(levels(cut_interval(tda[value < 0]$value, n.neg)), 
                    levels(cut_interval(tda[value >= 0 & value <=cut.pt]$value, n.pos1)),
                    levels(cut_interval(tda[value >cut.pt]$value, n.pos2))) %>% rev
    
    tda = tda %>% mutate(value2 = factor(value2, levels = cols.levels))
    
    cols.pos = brewer.pal(c(n.pos1*2+n.pos2), "YlOrRd") 
    cols.pos = c(cols.pos[seq(1, c(n.pos1)*2, by = 2)], tail(cols.pos, n.pos2)) # c(head(cols.pos, n.pos1), )
    cols.t = c(head(brewer.pal(c(n.pos1+n.pos2), "Blues"), n.neg) %>% rev, cols.pos)  %>% rev
    
    p =  ggplot(tda, aes(x = loc, y = f3, fill = value2)) +
      geom_tile() + 
      facet_wrap(~ f1 + f2, ncol = ncol.t) +
      # scale_fill_distiller(palette = "RdPu", direction = direction.t) +
      scale_fill_manual(values = cols.t, na.value = alpha('grey', 50)) +
      ggtitle(title.t) + labs(x = '', y = '', fill = fill.name.t) + 
      theme_minimal() + theme.tt
  } else if(min(tda$value, na.rm = T) < -100 & nrow(tda[value < - 100]) >= 2){
    cut.pt = -100
    n.pos = 4
    n.neg1 = 4
    n.neg2 = 1
    tda$value2 = ''
    tda[value < cut.pt]$value2 = cut_interval(tda[value <cut.pt]$value, n.neg2) # %>% as.character()
    tda[value >= cut.pt & value <0]$value2 = cut_interval(tda[value >= cut.pt & value <0]$value, n.neg1) # %>% as.character()
    tda[value >= 0]$value2 = cut_interval(tda[value >= 0]$value, n.pos) # %>% as.character()
    
    
    cols.levels = c(levels(cut_interval(tda[value < cut.pt]$value, n.neg2)), 
                    levels(cut_interval(tda[value >= cut.pt & value <0]$value, n.neg1)),
                    levels(cut_interval(tda[value >= 0]$value, n.pos))) %>% rev
    
    tda = tda %>% mutate(value2 = factor(value2, levels = cols.levels))
    
    cols.neg = brewer.pal(c(n.neg1*2+n.neg2), "Blues") 
    cols.neg = c(tail(cols.neg, n.neg2), cols.neg[seq(1, c(n.neg1)*2, by = 2)] %>% rev)  # c(head(cols.pos, n.pos1), )
    cols.t = c(cols.neg, head(brewer.pal(c(n.neg1+n.neg2), "YlOrRd"), n.pos)) %>% rev
    
    p =  ggplot(tda, aes(x = loc, y = f3, fill = value2)) +
      geom_tile() + 
      facet_wrap(~ f1 + f2, ncol = ncol.t) +
      # scale_fill_distiller(palette = "RdPu", direction = direction.t) +
      scale_fill_manual(values = cols.t, na.value = alpha('grey', 50)) +
      ggtitle(title.t) + labs(x = '', y = '', fill = fill.name.t) + 
      theme_minimal() + theme.tt
  } else {
    p = ggplot(tda, aes(x = loc, y = f3, fill = value)) +
      geom_tile() + 
      facet_wrap(~ f1 + f2, ncol = ncol.t) +
      # scale_fill_distiller(palette = "RdPu", direction = direction.t) +
      # scale_fill_gradient2(low = muted('blue'), high = muted('red')) +
      scale_fill_gradient2(low = 'blue', high = 'red', na.value = alpha('grey', 50)) +
      # scale_x_discrete(expand = c(0, 0)) +
      # scale_y_discrete(expand = c(0, 0)) +
      ggtitle(title.t) + labs(x = '', y = '', fill = fill.name.t) + 
      theme_minimal() + theme.tt
  }
 
  
  p
}

getPlot = function(tda){
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    # geom_point(mapping = aes(x = Week.start, y=obs)) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_rep_wrap(~loc + state, scales = 'free_y', repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),
                 labels = format(seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() +  theme.t # theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlot_cpObs = function(tda){
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    geom_point(mapping = aes(x = Week.start, y=obs)) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_rep_wrap(~loc + state, scales = 'free_y', repeat.tick.labels = T, ncol = 2, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),
                 labels = format(seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() + theme.t # theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlot_cpObs2 = function(tda){
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    geom_point(mapping = aes(x = Week.start, y=obs)) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_wrap(~loc + state, scales = 'free_y',  ncol = 2, labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),
                 labels = format(seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() + theme.t # theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlot_wVEsce = function(tda){
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    geom_point(mapping = aes(x = Week.start, y=obs)) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_rep_wrap(~loc + state + VEpriorInf, scales = 'free_y', repeat.tick.labels = T, ncol = 6) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),
                 labels = format(seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() + theme.t # theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlotStates_wVEsce = function(res, var = parm.names){
  p = ggplot(res[state %in% var]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    facet_rep_wrap(~ loc  + VEpriorInf + state, scales = 'free_y', repeat.tick.labels = T, ncol=length(var), labeller = label_wrap_gen(multi_line=FALSE)) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(res$Week.start), max(res$Week.start), by = 'month'),
                 labels = format(seq(min(res$Week.start), max(res$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlotStates = function(res, var = parm.names){
  p = ggplot(res[state %in% var]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    facet_rep_wrap(~ loc + state, scales = 'free_y', repeat.tick.labels = T, ncol=length(var)) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(res$Week.start), max(res$Week.start), by = 'month'),
                 labels = format(seq(min(res$Week.start), max(res$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlotProj = function(train.t, proj.t){
  
  p = ggplot(train.t) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t, aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t, aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t, mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t, mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure, scales = 'free_y', repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),
                 labels = format(seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),'%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}

getPlotProj9595 = function(train.t, proj.t){
  
  p = ggplot(train.t) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t, aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t, aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t, mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t, mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure, scales = 'free_y', repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),
                 labels = format(seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),'%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}

getPlotProj80 = function(train.t, proj.t){
  
  p = ggplot(train.t) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t, aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t, aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t, mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t, mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure, scales = 'free_y', repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),
                 labels = format(seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),'%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}

getPlotProj1 = function(train.t, proj.t){

  p = ggplot(train.t) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t, aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t, aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t, mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t, mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure + seasonality, scales = 'free_y', repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),
                 labels = format(seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),'%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}

# diff deflation factor
getPlotProj9580deflat = function(train.t, proj.t, mea.t, ptitle.t){
  
  dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
  if(length(dates.t) > 20)
    dates.t = dates.t[seq(1,length(dates.t), by = 8)]
  
  ncol.t = proj.t$fcast.deflat %>% unique %>% length()
  
  # mea.t = 'Cases'
  p = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ seasonality +  fcast.deflat, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR, 80% CI)', title = ptitle.t) +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  
  p
}

getPlotProj9580deflat1sn = function(train.t, proj.t, mea.t, ptitle.t, ylab.t = 'Estimate/Forecast (median, IQR, 80% CI)', theme.t){
  
  dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
  if(length(dates.t) > 20)
    dates.t = dates.t[seq(1,length(dates.t), by = 8)]
  
  ncol.t = proj.t$fcast.deflat %>% unique %>% length()
  
  # mea.t = 'Cases'
  p = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ fcast.deflat, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = '', y = ylab.t, title = ptitle.t) +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  
  p
}


getPlotProj9580Sce1sn = function(train.t, proj.t, mea.t, ptitle.t, ylab.t = 'Estimate/Forecast (median, IQR, 80% CI)', theme.t){
  
  dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
  if(length(dates.t) > 20)
    dates.t = dates.t[seq(1,length(dates.t), by = 8)]
  
  ncol.t = proj.t$scenario %>% unique %>% length()
  
  # mea.t = 'Cases'
  p = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ scenario, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = '', y = ylab.t, title = ptitle.t) +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}


getPlotProj9580Sns = function(train.t, proj.t, mea.t, ptitle.t, ylab.t = 'Estimate/Forecast (median, IQR, 80% CI)', theme.t){
  
  dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
  if(length(dates.t) > 20)
    dates.t = dates.t[seq(1,length(dates.t), by = 8)]
  
  ncol.t = proj.t$seasonality %>% unique %>% length()
  
  # mea.t = 'Cases'
  p = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ seasonality, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = '', y = ylab.t, title = ptitle.t) +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}

getPlotProj9580_2sn = function(train.t, proj.t, mea.t, sn0 = 'No seasonality', sn1 = 'Transformed seasonality', ptitle.t, theme.t = theme.t6, col.train = 'blue', dates.t=NULL, ymax = NULL, ymin = 0){
  
 
  if(is.null(dates.t)){
    dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
    if(length(dates.t) > 20)
      dates.t = dates.t[seq(1,length(dates.t), by = 8)]
  } # otherwise, use fixed time frame
   
  
  # ncol.t = proj.t$fcast.deflat %>% unique %>% length()
  col0 = 'blue'; col1 = 'red'; 
  # mea.t = 'Cases'
  if(! is.null(ymax)){
    p = ggplot(train.t[measure==mea.t & seasonality == sn0]) +
      geom_line(aes(x = Week.start, y = median), color = col.train) +  # no ctrl
      geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = col.train, alpha = .1) +
      geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col.train, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & seasonality == sn0], aes(x = Week.start, y = median), color = col0) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn0], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col0, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn0], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col0, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & seasonality == sn1], aes(x = Week.start, y = median), color = col1) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn1], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col1, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn1], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col1, alpha = .2) +
      geom_vline(data = train.t[measure==mea.t & seasonality == sn0], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
      geom_point(data = train.t[measure==mea.t & seasonality == sn0], mapping = aes(x = Week.start, y=obs)) + 
      geom_point(data = proj.t[measure==mea.t & seasonality == sn0], mapping = aes(x = Week.start, y=obs)) + 
      # facet_rep_wrap(~ seasonality +  fcast.deflat, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
      labs(x = '', y = 'Estimate/Forecast (median, IQR, 80% CI)', title = ptitle.t) +
      scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d'), limits = c(as.Date(min(dates.t)), as.Date(max(dates.t)))) +
      # ylim(c(ymin.t, ymax.t)) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  } else {
    p = ggplot(train.t[measure==mea.t & seasonality == sn0]) +
      geom_line(aes(x = Week.start, y = median), color = col.train) +  # no ctrl
      geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = col.train, alpha = .1) +
      geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col.train, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & seasonality == sn0], aes(x = Week.start, y = median), color = col0) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn0], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col0, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn0], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col0, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & seasonality == sn1], aes(x = Week.start, y = median), color = col1) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn1], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col1, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & seasonality == sn1], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col1, alpha = .2) +
      geom_vline(data = train.t[measure==mea.t & seasonality == sn0], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
      geom_point(data = train.t[measure==mea.t & seasonality == sn0], mapping = aes(x = Week.start, y=obs)) + 
      geom_point(data = proj.t[measure==mea.t & seasonality == sn0], mapping = aes(x = Week.start, y=obs)) + 
      # facet_rep_wrap(~ seasonality +  fcast.deflat, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
      labs(x = '', y = 'Estimate/Forecast (median, IQR, 80% CI)', title = ptitle.t) +
      scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d'), limits = c(as.Date(min(dates.t)), as.Date(max(dates.t)))) +
      theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  }
  
  p
}

getPlotProj9580_2sce = function(train.t, proj.t, mea.t, sce0 = 'asIs', sce1 = 'newI', ptitle.t, theme.t = theme.t6, col.train = 'blue', dates.t=NULL, ymax = NULL, ymin = 0){
  
  
  if(is.null(dates.t)){
    dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
    if(length(dates.t) > 20)
      dates.t = dates.t[seq(1,length(dates.t), by = 8)]
  } # otherwise, use fixed time frame
  
  
  # ncol.t = proj.t$fcast.deflat %>% unique %>% length()
  col0 = 'blue'; col1 = 'red'; 
  # mea.t = 'Cases'
  if(! is.null(ymax)){
    p = ggplot(train.t[measure==mea.t]) +
      geom_line(aes(x = Week.start, y = median), color = col.train) +  # no ctrl
      geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = col.train, alpha = .1) +
      geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col.train, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & scenario == sce0], aes(x = Week.start, y = median), color = col0) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce0], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col0, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce0], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col0, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & scenario == sce1], aes(x = Week.start, y = median), color = col1) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce1], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col1, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce1], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col1, alpha = .2) +
      geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
      geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
      geom_point(data = proj.t[measure==mea.t & scenario == sce0], mapping = aes(x = Week.start, y=obs)) + 
      # facet_rep_wrap(~ seasonality +  fcast.deflat, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
      labs(x = '', y = 'Estimate/Forecast (median, IQR, 80% CI)', title = ptitle.t) +
      scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d'), limits = c(as.Date(min(dates.t)), as.Date(max(dates.t)))) +
      # ylim(c(ymin.t, ymax.t)) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  } else {
    p = ggplot(train.t[measure==mea.t]) +
      geom_line(aes(x = Week.start, y = median), color = col.train) +  # no ctrl
      geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = col.train, alpha = .1) +
      geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col.train, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & scenario == sce0], aes(x = Week.start, y = median), color = col0) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce0], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col0, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce0], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col0, alpha = .2) +
      geom_line(data=proj.t[measure==mea.t & scenario == sce1], aes(x = Week.start, y = median), color = col1) +  # no ctrl
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce1], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = col1, alpha = .1) +
      geom_ribbon(data=proj.t[measure==mea.t & scenario == sce1], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = col1, alpha = .2) +
      geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
      geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
      geom_point(data = proj.t[measure==mea.t & scenario == sce0], mapping = aes(x = Week.start, y=obs)) + 
      # facet_rep_wrap(~ seasonality +  fcast.deflat, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
      labs(x = '', y = 'Estimate/Forecast (median, IQR, 80% CI)', title = ptitle.t) +
      scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d'), limits = c(as.Date(min(dates.t)), as.Date(max(dates.t)))) +
      theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  }
  
  p
}


getPlotProj9580 = function(train.t, proj.t, ptitle.t){
  
  dates.t = c(train.t$Week.start, proj.t$Week.start) %>% unique %>% as.Date %>% sort
  if(length(dates.t) > 20)
    dates.t = dates.t[seq(1,length(dates.t), by = 4)]
  
  ncol.t = 4 # proj.t$seasonality %>% unique %>% length()
  
  mea.t = 'Cases'
  p1 = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure + seasonality+ fcast.type, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  mea.t = 'Deaths'
  p2 = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure  + seasonality+ fcast.type, scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t, labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  mea.t = 'Infections'
  p3 = ggplot(train.t[measure==mea.t]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t[measure==mea.t], aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = ci80.lwr, ymax = ci80.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t[measure==mea.t], aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t[measure==mea.t], aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    geom_point(data = proj.t[measure==mea.t], mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ measure + seasonality + fcast.type , scales = 'fixed', repeat.tick.labels = T, ncol = ncol.t,labeller = label_wrap_gen(multi_line=FALSE,width=60)) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = dates.t, labels = format(dates.t,'%Y/%m/%d')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p = grid.arrange(top = textGrob(ptitle.t, gp=gpar(fontsize=11,font=1)),
    grobs = list(p1, p2, p3),
    layout_matrix = rbind(c(1, 1),
                          c(2, 2),
                          c(3, 3))
  )
  p
}

getPlotMultiV = function(tda, y.lab = 'Estimate'){
  
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = v.median), color = 'blue', size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    facet_rep_wrap(~ variant, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median, IQR)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() +  theme.t #  + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlotMultiVoverlay = function(tda, title.t, y.lab = 'Estimate',withObs = F, obs.t = NULL){
 
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) + ggtitle(title.t) +
    geom_line(aes(x = Week.start, y = v.median, color = variant), size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = variant), alpha = .3) +
    facet_rep_wrap(~ measure, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median, IQR)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 20))
  
  if(withObs & is.null(obs.t$variant)){
    p = p + geom_point(data = obs.t[date %in% dates.t], aes(x = date, y = value))
  } else if (withObs & !is.null(obs.t$variant)){
    p = p + geom_point(data = obs.t[date %in% dates.t], aes(x = date, y = value, color = variant, shape = variant))
  } 
  p
}

getPlotMultiVoverlayMedian = function(tda, title.t, y.lab = 'Estimate'){
  
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) + ggtitle(title.t) +
    geom_line(aes(x = Week.start, y = v.median, color = variant), size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    # geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = variant), alpha = .3) +
    facet_rep_wrap(~ measure, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() +  theme.t # + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 20))
  
  p
}


getPlotMultiVage.overlay = function(tda, title.t, y.lab = 'Estimate'){
  
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) + ggtitle(title.t) +
    geom_line(aes(x = Week.start, y = v.median, color = age.grp), size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = age.grp), alpha = .15) +
    facet_rep_wrap(~ variant, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median, IQR)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal()+  theme.t #  + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 20))
  
  p
}

