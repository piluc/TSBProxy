library(ggplot2)
library(scales)
library(gridExtra)
library(tools)
library(plyr)
library(xtable)
library(stringi)

compute_statistics = function(data,
                              measurevar,
                              groupvars,
                              conf.interval = 0.95){
  data = na.omit(data[order(data[[sequence_column]]),])
  data_c = ddply(.data = data, groupvars, .fun = function(xx, col) {
    c(N = length(xx[[col]]), mean = mean(xx[[col]]), sd = sd(xx[[col]]))
  }, measurevar)
  data_c = rename(data_c, c(mean = measurevar))
  data_c$se = data_c$sd/sqrt(data_c$N)
  conf_interval_multiplier = qt(conf.interval/2 + 0.5, data_c$N - 1)
  data_c$ci = data_c$se * conf_interval_multiplier 
  
  print(data_c)
  return(data_c)
}

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == 'guide-box')
  step3 <- step1$grobs[[step2]]
  return(step3)
}

plot2d = function(table,
                  meansx,
                  meansy,
                  x_value,
                  y_value,
                  sequence_column,
                  sequence_column_levels,
                  colors,
                  shapes,
                  plotwidth = 15,
                  plotheight = 5,
                  point_size = 4.5,
                  stroke_size = 1.5,
                  font_base_size = 15,
                  errorbar_size = 1.25
) {
  
  table[[sequence_column]] =  factor(table[[sequence_column]], levels = sequence_column_levels)
  table[['network']] = factor(table[['network']], levels = networks)
  meansx[[sequence_column]] = factor(meansx[[sequence_column]], levels = sequence_column_levels)
  meansy[[sequence_column]] = factor(meansy[[sequence_column]], levels = sequence_column_levels)
  table = na.omit(table[order(table[[sequence_column]]),])
  meansx = na.omit(meansx[order(meansx[[sequence_column]]),])
  meansy = na.omit(meansy[order(meansy[[sequence_column]]),])
  
  min_x = min(min(table[[x_value]]), min(meansx[[x_value]] - meansx[['sd']])) 
  max_x = max(max(table[[x_value]]), max(meansx[[x_value]] + meansx[['sd']])) 
  min_y = min(min(table[[y_value]]), min(meansy[[y_value]] - meansy[['sd']]))
  max_y = max(max(table[[y_value]]), max(meansy[[y_value]] + meansy[['sd']]))
  
  theme_set(theme_bw(base_size = font_base_size))
  
  scatter_p = ggplot(table,
                     aes(x=get(x_value),
                         y=get(y_value),
                         shape=factor(network),
                         colour=factor(get(sequence_column)))) +
    geom_point(size = point_size, stroke = stroke_size) +
    labs(x = 'weighted kendall\'s tau', y = 'time ratio', fill='proxy') +
    scale_x_continuous(limits=c(min_x, max_x)) +
    scale_y_continuous(limits=c(min_y, max_y)) +
    scale_shape_manual(values=shapes) +
    scale_colour_manual(values = colors) +
    theme_bw(base_size = font_base_size) +
    theme(legend.title = element_blank(),
          legend.position= 'bottom',
          plot.margin=unit(c(0, 0, 0, 0.5), 'lines'),
          legend.text = element_text(size=15)) 

  shared_legend = extract_legend(scatter_p)
  scatter_p = scatter_p + theme(legend.position = 'none')
  
  theme0 = function(...) theme( legend.position = 'none',
                                panel.background = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.ticks = element_blank(),
                                axis.text.x = element_blank(),
                                axis.text.y = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.length = unit(0, 'null'),
                                panel.border=element_rect(color = NA),
                                ...)
  
  x_p = ggplot(meansx, aes_string(x = x_value, y = 0,
                                  color=sequence_column)) +
    geom_point(size = point_size, stroke = stroke_size) +
    scale_x_continuous(limits=c(min_x, max_x)) +
    scale_colour_manual(values = colors) +
    theme_bw() +
    theme0(plot.margin = unit(c(0.25, 0., 0.25, 4.5), 'lines'))

  y_p = ggplot(meansy, aes_string(x=y_value, y=0,
                                  color=sequence_column)) +
    geom_point(size = point_size, stroke = stroke_size) +
    scale_x_continuous(limits=c(min_y, max_y)) +
    scale_colour_manual(values = colors) +
    theme_bw() +
    coord_flip()  +
    theme0(plot.margin = unit(c(0, 0, 2.5, 0 ), 'lines'))
  
  pdf(NULL)
  plot = grid.arrange(
    arrangeGrob(x_p, ncol = 2, widths = c(3, 0.15)),
    arrangeGrob(scatter_p, y_p, ncol = 2, widths = c(3, 0.15)),
    shared_legend,
    heights = c(0.3, 3, 0.87))

  ggsave(paste('experiment_1/plot-', x_value, '-', y_value, '.pdf', sep = ''), 
    plot=plot, width = plotwidth, height = plotheight, limitsize = FALSE)
  system("rm Rplots.pdf")
}

latex_format = function(table_proxies, table_onbras){
  
  time_ratios =  c('time_ratio_prefix', 'time_ratio_egostb', 'time_ratio_egoprefix', 'time_ratio_ptd', 'time_ratio_onbra', 'time_ratio_onbra_std')
  spearmans = c('spearman_prefix', 'spearman_egostb', 'spearman_egoprefix', 'spearman_ptd', 'spearman_onbra', 'spearman_onbra_std')
  ktaus = c('ktau_prefix', 'ktau_egostb', 'ktau_egoprefix', 'ktau_ptd', 'ktau_onbra', 'ktau_onbra_std')
  wktaus = c('wktau_prefix', 'wktau_egostb', 'wktau_egoprefix', 'wktau_ptd', 'wktau_onbra', 'wktau_onbra_std')
  
  column_names1 = c('Network', '$t(\\textsc{TempBrandes})$', time_ratios, wktaus)
  table1 = data.frame(matrix(ncol = 14, nrow = 0))
  colnames(table1) = column_names1
  
  column_names2 = c('Network', spearmans, ktaus)
  table2 = data.frame(matrix(ncol = 13, nrow = 0))
  colnames(table2) = column_names2
  
  for(i in 1:nrow(table_proxies)) {
    row = table_proxies[i, ]
    row_o = table_onbras[i, ]
    print(row['Network'])
    print(row_o['Network'])
    stopifnot(row_o['Network'] == row['Network'])
    
    exact_time = row['time_exact']
    
    new_row1 = data.frame(matrix(ncol=14, nrow=0))
    new_row1 = rbind(new_row1, 
            c(paste('\\texttt{', row['Network'], '}', sep=''), exact_time, 
                row['time_prefix']/exact_time, row['time_egotsb']/exact_time, row['time_egoprefix']/exact_time, row['time_ptd']/exact_time, row_o['avg_time_onbra_twice']/exact_time, row_o['std_time_onbra_twice']/exact_time,
                row['wktau_prefix'], row['wktau_egotsb'], row['wktau_egoprefix'], row['wktau_ptd'], row_o['avg_wktau_onbra_twice'], row_o['std_wktau_onbra_twice']))

    new_row2 = data.frame(matrix(ncol=13, nrow=0))
    new_row2 = rbind(new_row2, 
            c(paste('\\texttt{', row['Network'], '}', sep=''), 
                row['spearman_prefix'], row['spearman_egotsb'], row['spearman_egoprefix'], row['spearman_ptd'], row_o['avg_spearman_onbra_twice'], row_o['std_spearman_onbra_twice'],
                row['ktau_prefix'], row['ktau_egotsb'], row['ktau_egoprefix'], row['ktau_ptd'], row_o['avg_ktau_onbra_twice'], row_o['std_ktau_onbra_twice']))
                
    colnames(new_row1) = column_names1
    colnames(new_row2) = column_names2
    
    table1 = rbind(table1, new_row1)
    table2 = rbind(table2, new_row2)
    
  }
  
  
  table1['$t(\\textsc{TempBrandes})$'] = sapply(sanitize.numbers(format(table1['$t(\\textsc{TempBrandes})$'], digits=1)[[1]], type = "latex"), function(x) paste("$", x, "$", sep=""))
  
  for(cname in c(time_ratios)){
    table1[cname] = sanitize.numbers(format(table1[cname], scientific = TRUE, digits=2)[[1]], type = "latex", math.style.exponents = TRUE)
  }
  for(cname in c(wktaus)){
    table1[cname] = sapply(sanitize.numbers(format(table1[cname], digits=2)[[1]], type = "latex"), function(x) if(stri_sub(x, -1, -1) == 'N'){x}else{paste("$", x, "$", sep="")})
  }
  for(cname in c(spearmans, ktaus)){
    table2[cname] = sapply(sanitize.numbers(format(table2[cname], digits=2)[[1]], type = "latex"), function(x) if(stri_sub(x, -1, -1) == 'N'){x}else{paste("$", x, "$", sep="")})
  }
  
  for (i in 1:nrow(table1)) {
    row = table1[i, ]

    if(stri_sub(row['time_ratio_onbra'], -1, -1) != 'N'){
      row['time_ratio_onbra'] = paste(stri_sub(row['time_ratio_onbra'], 1 , -2), stri_sub(row['time_ratio_onbra_std'], 2 , -1), sep=" \\pm ")
    }
    if(stri_sub(row['wktau_onbra'], -1, -1) != 'N'){
      row['wktau_onbra'] = paste(stri_sub(row['wktau_onbra'], 1 , -2), stri_sub(row['wktau_onbra_std'], 2 , -1), sep=" \\pm ")
    }
    table1[i, ] = row
  }
  
  for (i in 1:nrow(table2)) {
    row = table2[i, ]

    if(stri_sub(row['spearman_onbra'], -1, -1) != 'N'){
      row['spearman_onbra'] = paste(stri_sub(row['spearman_onbra'], 1 , -2), stri_sub(row['spearman_onbra_std'], 2 , -1), sep=" \\pm ")
    }
    if(stri_sub(row['ktau_onbra'], -1, -1) != 'N'){
      row['ktau_onbra'] = paste(stri_sub(row['ktau_onbra'], 1 , -2), stri_sub(row['ktau_onbra_std'], 2 , -1), sep=" \\pm ")
    }
    table2[i, ] = row
  }
  
  table1[, "time_ratio_onbra_std"] = NULL
  table1[, "wktau_onbra_std"] = NULL
  table2[, "spearman_onbra_std"] = NULL
  table2[, "ktau_onbra_std"] = NULL
  
  table1 = data.frame(lapply(table1, function(x) {gsub("times", "cdot", x)}))
  
  print(colnames(table1))
  table11 = table1[, c('Network', 'X.t..textsc.TempBrandes...', 'time_ratio_prefix', 'time_ratio_egostb', 'time_ratio_egoprefix', 'time_ratio_ptd', 'time_ratio_onbra')]
  table12 = table1[, c('Network', 'wktau_prefix', 'wktau_egostb', 'wktau_egoprefix', 'wktau_ptd', 'wktau_onbra')]
  
  colnames(table11) = c('Network', '$t(\\textsc{TempBrandes})$', 
                    '\\textsc{Prefix}', '\\textsc{EgoSTB}', '\\textsc{EgoPrefix}', '\\textsc{PTD}', '\\textsc{Onbra}')
  colnames(table12) = c('Network',
                    '\\textsc{Prefix}', '\\textsc{EgoSTB}', '\\textsc{EgoPrefix}', '\\textsc{PTD}', '\\textsc{Onbra}')
  colnames(table2) = c('Network',
                    '\\textsc{Prefix}', '\\textsc{EgoSTB}', '\\textsc{EgoPrefix}', '\\textsc{PTD}', '\\textsc{Onbra}', 
                    '\\textsc{Prefix}', '\\textsc{EgoSTB}', '\\textsc{EgoPrefix}', '\\textsc{PTD}', '\\textsc{Onbra}')
  rename_pairs = list( c("venice", "Venice"), 
                    c("college_msg", "College msg"), 
                    c("email_eu", "Email EU"), 
                    c("bordeaux", "Bordeaux"),  
                    c("infectious", "Infectious"), 
                    c("topology", "Topology"),
                    c("wiki_elections", "Wiki elections"), 
                    c("facebook_wall", "Facebook wall"), 
                    c("digg_reply", "Digg reply"))
  for(p in rename_pairs){ 
    table11['Network'] = unlist(lapply(table11['Network'], gsub, pattern = p[1], replacement = p[2], fixed = TRUE))
    table12['Network'] = unlist(lapply(table12['Network'], gsub, pattern = p[1], replacement = p[2], fixed = TRUE))
    table2['Network'] = unlist(lapply(table2['Network'], gsub, pattern = p[1], replacement = p[2], fixed = TRUE))

  }  
  write.table(table11, file="latex_tables/experiment11.txt", sep="    &    ", quote=FALSE, col.names = TRUE, row.names = FALSE, eol='\\\\\\hline\n')
  write.table(table12, file="latex_tables/experiment12.txt", sep="    &    ", quote=FALSE, col.names = TRUE, row.names = FALSE, eol='\\\\\\hline\n')
  write.table(table2, file="latex_tables/experiment13.txt", sep="    &    ", quote=FALSE, col.names = TRUE, row.names = FALSE, eol='\\\\\\hline\n')
}

reformat_proxies = function(table){
  cat("\nProxies input table:\n")
  print(table)
  column_names = c('network', 'algorithm', 'running_time', 'time_ratio', 'wktau')
  new_table = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(new_table) = column_names
  
  for(i in 1:nrow(table)) {
    row = table[i, ]
    exact_time = row['time_exact']
    
    exact_row = c(row['Network'], 'exact', exact_time, 1, 1)
    names(exact_row) = column_names
    new_table = rbind(new_table, exact_row)
    
    for(proxy in proxies){
      new_row = c(row['Network'], proxy, row[paste('time_', proxy, sep='')], 
      row[paste('time_', proxy, sep='')] / exact_time, 
      row[paste('wktau_', proxy, sep='')])
      names(new_row) = column_names
      new_table = rbind(new_table, new_row)
    }
  }
  cat("\nReformatted proxies table:\n")
  print(new_table)
}


reformat_onbras = function(table){
  cat("\nOnbra input table:\n")
  print(table)
  column_names = c('network', 'algorithm', 'running_time', 'time_ratio', 'wktau')
  new_table = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(new_table) = column_names
  
  for(i in 1:nrow(table)) {
    row = table[i, ]
    exact_time = row['time_exact']
    
    exact_row = c(row['Network'], 'exact', exact_time, 1, 1)

    for(onbra in onbras){
      new_row = c(row['Network'], onbra, row[paste('avg_time_', onbra, sep='')], 
      row[paste('avg_time_', onbra, sep='')] / exact_time, 
      row[paste('avg_wktau_', onbra, sep='')])
      names(new_row) = column_names
      new_table = rbind(new_table, new_row)
    }
  }
  cat("\nReformatted onbras table:\n")
  print(new_table)
}

cr_plot = function(table,
                  network,
                  x_value,
                  y_value,
                  sequence_column,
                  sequence_column_levels,
                  colors,
                  shapes,
                  nu_rows,
                  plotwidth = 11,
                  plotheight = 5.5,
                  point_size = 3.5,
                  stroke_size = 1.5,
                  font_base_size = 25,
                  errorbar_width = 2,
                  errorbar_size = 1.25
) {
    
  data = table
  data[[sequence_column]] = factor(data[[sequence_column]], levels = sequence_column_levels)
  data = na.omit(data[order(data[[sequence_column]]),])
  
  theme_set(theme_bw(base_size = font_base_size))
  
  plot = ggplot(data, aes_string(x=x_value, y=y_value, shape=sequence_column, color=sequence_column))
  plot = plot + geom_point(size = point_size, stroke = stroke_size)
  plot = plot + scale_shape_manual(values=shapes)
  plot = plot + scale_colour_manual(values = colors)
  plot = plot + labs(x = "time ratio", y = "weighted kendall\'s tau", fill="Algorithm")
  plot = plot + theme_bw(base_size = font_base_size) +
              theme(legend.title = element_blank(),
                    legend.position= "bottom",
                    legend.text = element_text(size=20))

  
  for(i in seq.int(1, 3)){
    comp = sequence_column_levels[i]
    plot = plot + geom_hline(yintercept = table[table$algorithm==comp, y_value], color = colors[i], size= 1)
    plot = plot + geom_vline(xintercept = table[table$algorithm==comp, x_value], color = colors[i], size= 1)
  }

  # 
  plot_f_name = paste("experiment_2/", network, "-", x_value, "-", y_value, ".pdf", sep = "")
  ggsave(plot_f_name, width = plotwidth, height = plotheight, limitsize = FALSE)
}

cr_last_plot = function(table,
                  x_value,
                  y_value,
                  filename,
                  plotwidth = 11,
                  plotheight = 5.5,
                  point_size = 3.5,
                  stroke_size = 1.5,
                  font_base_size = 20,
                  errorbar_width = 2,
                  errorbar_size = 1.25
) {
    
  data = table
  data[["network"]] = factor(data[["network"]])
  
  theme_set(theme_bw(base_size = font_base_size))
  plot = ggplot(data, aes_string(x=x_value, y=y_value, group=1))
  plot = plot + geom_point(size = point_size, stroke = stroke_size) 
  plot = plot + scale_x_discrete(breaks=c("1","2","3", "4", "5", "6"), labels=c("Slashdot reply", " Mathoverflow", "Wiki talk", "Email Enron", "Askubuntu", "Superuser"))#labels=c("1" = "x", "2" = "x", "3" = "y", "4" = "z", "5" = 't', "6" = 'u'))
  plot = plot + geom_line()
  plot = plot + labs(x = "Network", y = "time ratio", fill="Algorithm")
  plot = plot + theme_bw(base_size = font_base_size) +
              theme(legend.title = element_blank(),
                    legend.position= "bottom",
                    legend.text = element_text(size=20))

  ggsave(filename, width = plotwidth, height = plotheight, limitsize = FALSE)
}



## Experiment 1
proxies = c('prefix', 'ptd')
onbras = c('onbra_twice')
networks = c('venice', 'college_msg', 'email_eu', 'bordeaux', 'infectious', 'SMS', 'topology', 'wiki_elections', 'facebook_wall', 'digg_reply')

# read
proxies_table = read.table('./proxies.csv', header = TRUE, sep=';')
onbras_table = read.table('./onbras.csv', header = TRUE, sep=';')

dir.create('latex_tables')

latex_format(proxies_table, onbras_table)

proxies_table = reformat_proxies(proxies_table)
onbras_table = reformat_onbras(onbras_table)
onbras_table = onbras_table[onbras_table[, 'algorithm']=='onbra_twice', ]
onbras_table['algorithm'] = unlist(lapply(onbras_table['algorithm'], gsub, pattern = 'onbra_twice', replacement = 'onbra', fixed = TRUE))

table = rbind(proxies_table, onbras_table)
cat("\nComplete reformatted table:\n")
print(table)

dir.create('experiment_1')
# plot
x_value = 'wktau'
y_value = 'time_ratio'
sequence_column = c('algorithm')
sequence_column_levels = c(proxies, 'onbra')

colors = c('#b8860b', '#a30000', '#000057')
shapes = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

cat(paste("\nStatistics for ", x_value, ":\n", sep=""))
meansx = compute_statistics(table, x_value, sequence_column)
cat(paste("\nStatistics for ", y_value, ":\n", sep=""))
meansy = compute_statistics(table, y_value, sequence_column)
plot2d(table, meansx, meansy, x_value, y_value, sequence_column, sequence_column_levels, colors, shapes)



## Experiment 2
networks = c('venice', 'college_msg', 'email_eu', 'bordeaux', 'infectious', 'SMS', 'topology', 'wiki_elections', 'facebook_wall', 'digg_reply')
rfolder = 'onbra_evolution'
folders = c('01_venice', '02_college_msg', '03_email_eu', '04_bordeaux', '06_infectious', '07_SMS', '08_topology', '09_wiki_elections', '10_facebook_wall', '11_digg_reply')

dir.create('experiment_2')
colors = c('#b8860b', '#a30000', '#2F5233', '#000057')
shapes = c(0, 0, 0, 0)

for (f in folders){
  n = stri_sub(f, 4, -1)
  print(n)
  otable = read.table(paste(rfolder, f, 'evolution.txt', sep='/'), header=TRUE, sep=':')
  otable = cbind(otable, algorithm='onbra', network=n)

  # otable = cbind(otable, time_ratio=)
  others = table[table[,'network']==n,]
  others = cbind(others, num_samples=0)

  cnames = c('network', 'algorithm', 'running_time', 'wktau', 'num_samples')
  others = others[, cnames]
  colnames(others) = c('network', 'algorithm', 'time', 'wtau', 'num_samples')
  otable = otable[, c('network', 'algorithm', 'time', 'wtau', 'num_samples')]

  tab = rbind(others, otable)
  tab = cbind(tab, time_ratio = tab[, 'time']/tab[tab[,'algorithm']=='exact',][,'time'])

  x_value = 'time_ratio'
  y_value = 'wtau'

  sequence_column_levels = c('prefix', 'ptd', 'exact', 'onbra')

  cr_plot(tab, f, x_value, y_value, sequence_column, sequence_column_levels, colors, shapes, 1)


}


## Experiment 3
dir.create('experiment_3')
in_table = read.table('prac_vs_theo.txt', header=TRUE)
table1 = data.frame(matrix(ncol = 0, nrow = 6))
table2 = data.frame(matrix(ncol = 0, nrow = 6))
table1 = cbind(table1, time_ratio=in_table[, 'nMlogm'], tvp='theory', network = c(1, 2, 3, 4, 5, 6))
table2 = cbind(table2, time_ratio=in_table[, 'tPrefix.tPTD'], tvp='practice', network = c(1, 2, 3, 4, 5, 6))

cr_last_plot(table1, 'network', 'time_ratio', 'experiment_3/theory.pdf')
cr_last_plot(table2, 'network', 'time_ratio', 'experiment_3/practice.pdf')
