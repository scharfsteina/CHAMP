library(igraph)
library(igraphdata)
library(tidyverse)
library(ggthemes)
library(ggplot2)
library(reshape2)
champ <- function(network, min = 0, weighted = T, plot = T) {
  
  source("research/iterate_gamma.R")
  source("research/unique_partitions.R")
  source("research/best_partitions.R")
  source("research/partition_level_champ.R")
  source("research/find_min.R")
  
  if (weighted == F) {
    E(network)$weight <- 1
  }
  
  partitions <- unique_partitions(iterate_gamma(network))
  best_partitions <- best_partitions(c=partitions[[1]], p=partitions[[2]],min) 
  
  
  a <- list()
  p <- list()
  gammas <- seq(0,2,0.025)
  modularity <- array(NA, dim = c(length(gammas),
                                  length(best_partitions[[2]])))
  k <- 0
  for (p in best_partitions[[2]]) {
    k <- k+1
    ret <- partition_level_champ(p) # MAKE THIS MORE EFFICIENT
    a[k] <- ret$a_partition
    p[k] <- ret$p_partition
  }
  
  # Create descending lists
  ordered_a <- order(unlist(a), decreasing = TRUE) # correct descending order of a_hats
  a_descending <- rep(NA,length(ordered_a))
  p_descending <- rep(NA,length(ordered_a))
  partitions_descending <- rep(NA,length(ordered_a))
  for (i in ordered_a) {
    p_descending[i] <- p[ordered_a[i]]
    a_descending[i] <- a[ordered_a[i]]
    partitions_descending[i] <- best_partitions[[2]][ordered_a[i]]
  }
  
  gamma_matrix <- matrix(data = NA, nrow = length(ordered_a), ncol = length(ordered_a))
  for (i in 1:length(a_descending)) {
    for (j in 1:length(a_descending)) {
      if (!is.null(a_descending[[i]]) & !is.null(p_descending[[i]])) {
        gamma_matrix[i,j] <- (a_descending[[i]]-a_descending[[j]])/(p_descending[[i]]-p_descending[[j]])
      }
    }
  }
  
  minimum <- find_min(gamma_matrix)
  best_gammas = minimum$mins
  corresponding_partitions = minimum$partitions
  best_modularities <- c()
  for (g in 1:length(best_gammas)) {
    partition_index <- corresponding_partitions[g]
    best_modularities[g] <- a[[partition_index]]-(best_gammas[g]*p[[partition_index]])
  }
  
  for (k in 1:length(best_partitions[[2]])) {
    for (g in 1:length(gammas)) {
      modularity[g,k] <- a[[k]]-(gammas[g]*p[[k]])
    }
  }
  
  # Plot data frames
  
  if (plot == T) {
    all <- data.frame(x = gammas)
    for (i in 1:dim(modularity)[2]) {
      all <- cbind(all, i=modularity[,i])
    }
    colnames(all) <- c('x',paste("", 1:dim(modularity)[2], sep = ""))
    all <- melt(all, id = 'x')
    colnames(all)[2] <- "partition_num"
    best <- data.frame(best_gammas, best_modularities, corresponding_partitions)
    
    segments <- data.frame(x1 = c(0,best_gammas), 
                           y1 = c(modularity[1,1], best_modularities), 
                           x2 = c(best_gammas,NA), 
                           y2 = c(best_modularities,NA),
                           partitions = c(corresponding_partitions,0))
    
    ggplot() +
      geom_line(data = all,
                mapping = aes(x = x, 
                              y = value,
                              group = partition_num), 
                show.legend = F,
                color = all$partition_num,
                na.rm = T) +
      geom_point(data = best,
                 mapping = aes(x = best_gammas, 
                               y = best_modularities)) +
      geom_segment(data = best,
                   mapping = aes(x = best_gammas,
                                 xend = best_gammas,
                                 y = best_modularities,
                                 yend = -Inf),
                   linetype = "dashed") +
      scale_y_continuous(breaks = seq(0,max(segments$y1),length = 10), 
                         labels = round(seq(0,max(segments$y1),length = 10)),
                         expand = c(0,0),
                         limits = c(0,max(segments$y1))) +
      scale_x_continuous(breaks = c(segments$x1,2), 
                         labels = c(round(segments$x1,2),2),
                         limits = c(0,2),
                         expand = c(0,0)) +
      labs(x = expression(paste("Resolution Parameter (", gamma,")")),
           y = "Modularity",
           title = "CHAMP") + # ??
      theme_few() +
      theme(axis.text = element_text(size = 8))
    
    ggplot(segments) + 
      geom_segment(aes(x=x1, 
                       y=y1, 
                       xend = x2, 
                       yend = y2),
                   color = segments$partitions,
                   na.rm = T)+
      geom_segment(aes(x=x1,
                       y=y1,
                       xend=x1,
                       yend=-Inf), linetype = "dashed") +
      guides(color = "none") +
      labs(x = expression(paste("Resolution Parameter (", gamma,")")),
           y = "Modularity",
           title = "CHAMP")+ # ??
      scale_y_continuous(breaks = seq(0,max(segments$y1),length = 10), 
                         labels = round(seq(0,max(segments$y1),length = 10)),
                         expand = c(0,0),
                         limits = c(0,max(segments$y1))) +
      scale_x_continuous(breaks = c(segments$x1,2), 
                         labels = c(round(segments$x1,2),2),
                         limits = c(0,2),
                         expand = c(0,0)) +
      geom_point(aes(x = x1, y = y1)) +
      theme_few() +
      theme(axis.text = element_text(size = 8))
    
    ggplot() +
      geom_line(data = all,
                mapping = aes(x = x, 
                              y = value, 
                              group = partition_num), 
                show.legend = F,
                color = all$partition_num,
                alpha = .3,
                na.rm = T) +
      geom_segment(data = segments, 
                   mapping = aes(x=x1, 
                                 y=y1, 
                                 xend = x2, 
                                 yend = y2),
                   color = "#63666A",
                   size = 1.5,
                   na.rm = T) +
      geom_segment(data = best,
                   mapping = aes(x = best_gammas,
                                 xend = best_gammas,
                                 y = best_modularities,
                                 yend = -Inf),
                   linetype = "dashed",
                   color = "black",
                   na.rm = T) +
      geom_point(data = best,
                 mapping = aes(x = best_gammas, 
                               y = best_modularities),
                 color = "black",
                 na.rm = T) +
      labs(x = expression(paste("Resolution Parameter (", gamma,")")),
           y = "Modularity",
           title = "CHAMP")+
      scale_y_continuous(breaks = seq(0,max(segments$y1),length = 10), 
                         labels = round(seq(0,max(segments$y1),length = 10)),
                         expand = c(0,0),
                         limits = c(0,max(segments$y1))) +
      scale_x_continuous(breaks = c(segments$x1,2), 
                         labels = c(round(segments$x1,2),2),
                         limits = c(0,2),
                         expand = c(0,0)) +
      theme_few() +
      theme(axis.text = element_text(size = 8))
  }
  
  # Ordering line segments by length
  partition_summary <- data.frame(matrix(ncol = 5, nrow = length(segments[,1])-1))
  colnames(partition_summary) <- c("segment_length",
                                   "starting_gamma",
                                   "ending_gamma",
                                   "gamma_range",
                                   "partition_num")
  
  for (x in 1:lengths(partition_summary)[1]) {
    partition_summary[x,"segment_length"] <- sqrt((segments[x, "x1"]-segments[x, "x2"])**2+
                                                    (segments[x, "y1"]-segments[x, "y2"])**2)
    partition_summary[x,"starting_gamma"] <- segments[x,"x1"]
    partition_summary[x,"ending_gamma"] <- segments[x,"x2"]
    partition_summary[x,"gamma_range"] <- abs(segments[x,"x1"]-segments[x,"x2"])
    partition_summary[x,"partition_num"] <- segments[x,"partitions"]
    partition_summary[x,"num_clusters"] <- best_partitions[[2]][[x]]$nb_clusters
  }
  
  partition_summary[order(-partition_summary$gamma_range),]
  
  return(partition_summary)

}



