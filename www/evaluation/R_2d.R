library(ggplot2)
library(scales)
library(gridExtra)
library(tools)
library(plyr)

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
  data_c$ci = data_c$se * conf_interval_multiplier # round(x, digits = 3)
  data_file_name = paste('./data/data-', sub('[.]', '_', file_path_sans_ext(file_name)),
                         '-', measurevar, '.csv', sep = '')
  
  write.csv(data_c, data_file_name, row.names = FALSE)
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
                  plotwidth = 12,
                  plotheight = 5,
                  point_size = 3.5,
                  stroke_size = 1.5,
                  font_base_size = 16,
                  errorbar_size = 1.25
) {
  
  table[[sequence_column]] =  factor(table[[sequence_column]], levels = sequence_column_levels)
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
    labs(x = 'spearman correlation', y = 'running time', fill='proxy') +
    scale_x_continuous(limits=c(min_x, max_x)) +
    scale_y_continuous(limits=c(min_y, max_y)) +
    scale_shape_manual(values=shapes) +
    scale_color_brewer(palette='Dark2') +
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
    geom_point(size = point_size, stroke = stroke_size,
                position = position_jitter(h = 0.4, w = 0.0, seed = 123)) +
    scale_x_continuous(limits=c(min_x, max_x)) +
    scale_color_brewer(palette='Dark2') +
    theme_bw() +
    geom_errorbar(aes(xmin=.data[[x_value]] - .data[['sd']],
                      xmax=.data[[x_value]] + .data[['sd']]),
                  size = errorbar_size,
                  alpha=0.5, 
                  position = position_jitter(h = 0.4, w = 0.0, seed = 123)) + 
    theme0(plot.margin = unit(c(0.25, 0., 0.25, 4.5), 'lines'))

  y_p = ggplot(meansy, aes_string(x=y_value, y=0,
                                  color=sequence_column)) +
    geom_point(size = point_size, stroke = stroke_size,
                position = position_jitter(h = 0.4, w = 0.0, seed = 123)) +
    scale_x_continuous(limits=c(min_y, max_y)) +
    scale_color_brewer(palette='Dark2') +
    theme_bw() +
    geom_errorbar(aes(xmin=meansy[[y_value]] - meansy[['sd']],
                      xmax=meansy[[y_value]] + meansy[['sd']]),
                  size = errorbar_size,
                  alpha=0.5, 
                  position = position_jitter(h = 0.4, w = 0.0, seed = 123)) + 
    coord_flip()  +
    theme0(plot.margin = unit(c(0, 0.25, 2.5, 0.25), 'lines'))
  
  pdf(NULL)
  plot = grid.arrange(
    arrangeGrob(x_p, ncol = 2, widths = c(3, 0.15)),
    arrangeGrob(scatter_p, y_p, ncol = 2, widths = c(3, 0.15)),
    shared_legend,
    heights = c(0.3, 3, 0.87))

  plot_f_name = paste('./plots/plot-', sub('[.]', '_', file_path_sans_ext(file_name)),
                      '-', x_value, '-', y_value, '.pdf', sep = '')
  ggsave(plot_f_name, plot=plot, width = plotwidth, height = plotheight, limitsize = FALSE)
  system("rm Rplots.pdf")
}

reformat = function(table){
  cat("\nInput table:\n")
  print(table)
  column_names = c('network', 'algorithm', 'running_time', 'time_ratio', 'correlation')#, 'h_10', 'h_25', 'h_50', 'h_10.')
  new_table = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(new_table) = column_names
  
  for(i in 1:nrow(table)) {
    row = table[i, ]
    exact_time = row[paste('time_', 'exact', sep='')]
    
    exact_row = c(row['Network'], 'exact', exact_time, 1, 1)#, 10, 25, 50, NA)
    names(exact_row) = column_names
    new_table = rbind(new_table, exact_row)
    
    for(proxy in proxies){
      new_row = c(row['Network'], proxy, row[paste('time_', proxy, sep='')], 
      row[paste('time_', proxy, sep='')] / exact_time, 
      row[paste('spearman_', proxy, sep='')])
      # , row[paste('h_10_', proxy, sep='')], 
      # row[paste('h_25_', proxy, sep='')], row[paste('h_50_', proxy, sep='')], 
      # row[paste('h_10._', proxy, sep='')])
      names(new_row) = column_names
      new_table = rbind(new_table, new_row)
    }
    for(proxy in onbras){
      new_row = c(row['Network'], proxy, row[paste('avg_time_', proxy, sep='')], 
      row[paste('avg_time_', proxy, sep='')] / exact_time, 
      row[paste('avg_spearman_', proxy, sep='')])
      # , row[paste('h_10_', proxy, sep='')], 
      # row[paste('h_25_', proxy, sep='')], row[paste('h_50_', proxy, sep='')], 
      # row[paste('h_10._', proxy, sep='')])
      names(new_row) = column_names
      new_table = rbind(new_table, new_row)
    }
  }
  cat("\nReformatted table:\n")
  print(new_table)
}


experiments = list(c('./', 'results.csv'))
proxies = c('egotsb', 'egoprefix', 'prefix', 'ptd')
onbras = c('onbra_half', 'onbra_equal', 'onbra_twice')

for(experiment in experiments){
  folder_name = experiment[1]
  file_name = experiment[2]

  # read
  table = read.table(paste(folder_name, file_name, sep=''), header = TRUE, sep=':')
  table = reformat(table)

  # plot
  x_value = 'correlation'
  y_value = 'running_time'
  sequence_column = c('algorithm')
  sequence_column_levels = c(proxies, onbras)

  colors = c('#b8860b', '#a30000', '#000057', '#C77CFF')
  shapes = c(1, 2, 3, 4, 5, 6, 7, 8)

  dir.create('data')
  dir.create('plots')

  cat(paste("\nStatistics for ", x_value, ":\n", sep=""))
  meansx = compute_statistics(table, x_value, sequence_column)
  cat(paste("\nStatistics for ", y_value, ":\n", sep=""))
  meansy = compute_statistics(table, y_value, sequence_column)
  plot2d(table, meansx, meansy, x_value, y_value, sequence_column, sequence_column_levels, colors, shapes)
}
