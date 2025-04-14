library(ggstats)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rlang)

#define a function to visualize population statistics

plot_population_stats<-function(data, column=NULL, group=NULL, stat_filter=NULL,
                                total_pop=F, bins=50, group_means=T, 
                                stats_in_label=T, panel_labels="AUTO",
                                x_label=NULL, margin=unit(c(1,1.5,0,0), 'cm')){
  
  if(!missing(column)){
  
  
  column<-enquo(column)
  col_name <- quo_name(column)

  col_data <- data[[col_name]]
  is_factor_column<-is.factor(col_data)
  
  
  if (!is_factor_column){
    pop_mean <- mean(col_data)
  }else{
    group_means<-F

  }

  if (!missing(group) & quo_name(enquo(group))!="NULL"){
    
    group<-enquo(group)

    
    group_name<-quo_name(group)
    
    if(!is_factor_column){
      gp <- ggplot(data, aes(x=!!column, fill=!!group, after_stat(ndensity)))
    }else{
      gp <- ggplot(data, aes(x=!!column, fill=!!group))
    }
    
  }else{
    group_name<-"NULL"
    if(!is_factor_column){
      gp <- ggplot(data, aes(x=!!column, after_stat(ndensity)))
    }else{
      gp <- ggplot(data, aes(x=!!column))
    }
  }
  
  
  if(!is_factor_column){#add step plot, fill with histogram(s), and a vertical line for the pop. mean
    
    if(total_pop){
      gp<-gp+geom_step(stat="bin", bins=bins, direction="mid",
                       aes(x=!!column, fill=NULL),
                       linetype=6,
                       linewidth=0.8,
                       show.legend = c(fill=F, color=T))
    }
      gp<-gp+geom_histogram(color=NA, alpha=.3, bins=bins,
                     position = "identity")
    if (group_name!="NULL"){
      gp<-gp+geom_step(stat="bin", bins=bins, direction="mid", show.legend = F)
      }
    if (group_means){#add mean plot
      if(total_pop){
      gp<-gp+geom_vline(aes(xintercept=pop_mean, fill=NULL),
                        linetype="dotted", 
                        linewidth=1.5, 
                        show.legend = F)
      }
      if (group_name!="NULL"){#add corresponding group plots if necessary
        group_stats <- as_tibble(tapply(col_data, data[[group_name]], mean),
                                 rownames = group_name)
        gp<-gp+geom_vline(data=group_stats, aes(xintercept=value, color=!!group),
                     linetype="dotted", linewidth=1.5, show.legend = F)
      }
    }
    gp<-gp+labs(y="density (norm.)",
                x=ifelse(stats_in_label,
                         sprintf("%s (mean=%0.2g, std-dev=%0.2g)", col_name, pop_mean, sd(col_data)), ifelse(is.null(x_label),col_name,x_label)
                         )
                )
  }else{
    #add bar graph
    gp<-gp+geom_bar( width=0.5, alpha=.3, linewidth=0.8)+
      geom_text(aes(label = scales::percent(after_stat(prop), accuracy = 1)),
                stat = "prop",
                position = position_stack(.5)
      )+
      labs(y="count")
    
  }
  gp+theme(legend.position="top",text = element_text(size = 14))+theme(plot.margin=margin)
  }else{
    plts<-lapply(lapply(colnames(data), sym),
                 function(c) plot_population_stats(data, column=!!c, group=!!sym(quo_name(enquo(group))),#this seems like an extreme solution but idk how else to make things work
                                                   total_pop=total_pop, group_means=group_means, bins=bins,
                                                   stats_in_label=stats_in_label, margin=margin ) 
                )
    sq<-ceiling(sqrt(length(colnames(data))))
    ggarrange(plotlist=plts, ncol=sq, nrow=sq, common.legend=T, labels=panel_labels, align='hv')
  }
}



geom_cross_line<-function(x=0.0,y=0.0,...){
  list(geom_vline(aes(xintercept=x),...),geom_hline(aes(yintercept=y),...))
}

average_geom_point<-function(data, post_process=NULL, fill="", ...){
  
  data<-summarise_if(data,is.numeric, mean)
  if(!is.null(post_process)){
    data<-post_process(data)
  }
  
  if(is.function(fill)){
    fill<-fill(data)
  }
  
  list(geom_point(data=data, aes(fill=fill),...),
       guides(fill=guide_legend(title="averaged parameters")))
  
}

geom_smooth_group_and_mean<-function(group,linetypes=c("dotted","solid"), linewidths=c(0.5, 1), mean_color="black",...){
  list(geom_smooth(aes(group=!!enquo(group)),linetype=linetypes[[1]], linewidth=linewidths[[1]],...),
       geom_smooth(aes(group=NA),linetype=linetypes[[2]], linewidth=linewidths[[2]], color=mean_color,...))
  
}

captioned<-function(gp, text, horz_space=6, vert_space=1, heights=c(10,1)){
  ggarrange(gp, ggparagraph(text=text )+
                theme(plot.margin=unit(c(t = vert_space, r = horz_space, b = 1.5, l = horz_space),"lines")), 
            nrow=2, ncol=1, heights=heights, align='h')
}


filtered_stats_in_legend<-function(data, cond, group=NULL, scale=scale_fill_discrete, success_string=NULL){
  
  if(!missing(group)){
    if(is.null(success_string)){
      success_string<-"%s (%0.1f%% successful)"
    }
    group<-enquo(group)
    group_name<-quo_name(group)
    
    filtered_stats <- filtered_fraction(data, cond=!!enquo(cond), group=!!group)
    labels <- sprintf(success_string, filtered_stats[[group_name]], 100*filtered_stats$fraction)
  }else{
    
    labels<-sprintf(success_string, "", 100*filtered_fraction(data, cond=!!enquo(cond)))
  }
  
  scale(labels=labels)
}

filtered_group_stats_in_legend<-function(data, cond, group=NULL, scale=scale_fill_discrete, success_string=NULL){
  
  if(!missing(group)){
    if(is.null(success_string)){
      success_string<-"%s (%0.1f%% successful)"
    }
    group<-enquo(group)
    group_name<-quo_name(group)
    
    filtered_stats <- filtered_percent(data, cond=!!enquo(cond), group=!!group)
    labels <- sprintf(success_string, filtered_stats[[group_name]], filtered_stats$percent)
  }else{
    
    labels<-sprintf(success_string, "", filtered_percent(data, cond=!!enquo(cond)))
  }
  
  scale(labels=labels)
}


filtered_fraction<-function(data,cond,group=NULL){
  
  
  
  if(!missing(group)){
    
    group<-enquo(group)
    group_name<-quo_name(group)
    grouped<-group_by(data, !!group)
    
    filtered<-filter(grouped, !!enquo(cond), .preserve=T)
    
    out<-data.frame(fraction=count(filtered)$n/count(grouped)$n)
    out[[group_name]]<-factor(levels(factor(data[[group_name]])))
    out
    
  }else{
    (count(filter(data, !!enquo(cond)))$n/count(data)$n)[[1]]
  }
  
}

filtered_percent<-function(data,cond,group=NULL){
  if(!missing(group))
    fraction_stats<-filtered_fraction(data,!!enquo(cond),group=!!enquo(group))
  else{
    fraction_stats<-filtered_fraction(data,!!enquo(cond))
  }
  mutate(fraction_stats, percent=100*fraction, .keep="unused")
}




sub_sample <- function(data, n, group=NULL){

  if(!missing(group)){
    bind_rows(lapply(split(data,data[[as_string(enexpr(group))]]), function(data) sub_sample(data, n) ))
  }
  else{
    data[sample(nrow(data), n),]
  }

}



bind_with<-function(data,x, .id=NULL){
  n<-nrow(data)
  out <- vector("list", n)
  for(i in 1:n){
    out[[i]]<-eval_tidy(enquo(x), data[i,])
  }
  bind_rows(out,.id=.id)
}