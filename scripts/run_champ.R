# Parameters: 
# - network: a network you want to run community detection and champ analysis on
# - weighted: defaulted at false, performs network analysis on the unweighted network
# - min: defaulted at 0, minimum number of gamma values for which a partition is created
# - name: name of network (eg. "Karate")
# - plot: if you want to plot the partitions
# - n_plots: top # of partitions you want to plot 
run_champ <- function(network, weighted = F, min = 0, name, plot = T, n_plots = 2) {
  
  # Import relevant functions
  source("scripts/iterate_gamma.R")
  source("scripts/unique_partitions.R")
  source("scripts/best_partitions.R")
  source("scripts/partition_level_champ.R")
  source("scripts/find_min.R")
  
  # Load in libraries
  library(igraph)
  library(igraphdata)
  library(tidyverse)
  library(ggthemes)
  library(ggplot2)
  library(reshape2)
  
  # Turn network into an unweighted network
  if (weighted==F) {
    E(network)$weight <- 1
  }
  
  # Find the unique partitions
  partitions <- unique_partitions(iterate_gamma(network))
  
  # Extract the unique partitions that were created from a minimum number of gamma values
  # Equivalent to partitions if min is set to 0
  best_partitions <- best_partitions(partitions$count, partitions$partitions, min) 
  
  # CHAMP method:
  # Get A-hat and P-hat for each partition, gives you a intercept and slope. 
  # Vary gamma to get a line. Do this for each partition and you get a series of lines.
  a <- list(); p <- list(); gammas <- seq(0,2,0.025)
  modularity <- array(NA, dim = c(length(gammas), length(best_partitions$partitions)))
  k <- 0
  for (partition in best_partitions$partitions) {
    k <- k+1
    ret <- partition_level_champ(network, partition)
    a[k] <- ret$a_partition
    p[k] <- ret$p_partition
  }
  
  # Put a and p in descending order to construct a matrix of gamma values
  ordered_a <- order(unlist(a), unlist(p), decreasing = TRUE)
  a_descending <- rep(NA,length(ordered_a)); p_descending <- rep(NA,length(ordered_a))
  partitions_descending <- rep(NA,length(ordered_a))
  # Sort by a and then p (be careful if there are ties)
  for (i in ordered_a) {
    p_descending[i] <- p[ordered_a[i]]
    a_descending[i] <- a[ordered_a[i]] 
    partitions_descending[i] <- best_partitions$partitions[ordered_a[i]]
  }
  report <- ""
  gamma_matrix <- matrix(data = NA, nrow = length(ordered_a), ncol = length(ordered_a))
  prev <- NULL
  for (i in 1:length(a_descending)) {
    for (j in 1:length(a_descending)) {
      gamma_matrix[i,j] <- (a_descending[[i]]-a_descending[[j]])/(p_descending[[i]]-p_descending[[j]])
    }
    if (!is.null(prev)) {
      if (all(mapply(identical, gamma_matrix[i,], prev))) {
        report <- str_c(report,"Partition ", 
                        ordered_a[i-1], 
                        " and Partition ", 
                        ordered_a[i],
                        " have the same modularity line. ")
      }
    }
    prev <- gamma_matrix[i,]
  }
  
  # Find the upper envelope of gamma values
  min <- find_min(gamma_matrix)
  best_gammas <- min$mins
  corresponding_partitions <- min$partitions
  best_modularities <- c()
  for (g in 1:length(best_gammas)) {
    partition_index <- corresponding_partitions[g]
    best_modularities[g] <- a[[partition_index]]-(best_gammas[g]*p[[partition_index]])
  }
  
  for (k in 1:length(best_partitions$partitions)) {
    for (g in 1:length(gammas)) {
      modularity[g,k] <- a[[k]]-(gammas[g]*p[[k]])
    }
  }
  
  # Plot data frames
  all <- data.frame(x = gammas)
  for (i in 1:dim(modularity)[2]) {
    all <- cbind(all, i=modularity[,i])
  }
  colnames(all) <- c('x',paste("", 1:dim(modularity)[2], sep = ""))
  all <- melt(all, id = 'x')
  colnames(all)[2] <- "partition_num"
  last <- corresponding_partitions[length(corresponding_partitions)]
  best_gammas <- c(0,best_gammas,2)
  best_modularities <- c(modularity[1,1], best_modularities, a_descending[[last]]-(2*p_descending[[last]]))
  last <- length(best_modularities)
  best <- data.frame(best_gammas, best_modularities)
  segments <- data.frame(x1 = best_gammas[1:last-1], 
                         y1 = best_modularities[1:last-1], 
                         x2 = best_gammas[2:last], 
                         y2 = best_modularities[2:last],
                         partitions = corresponding_partitions)
  
  w <- ifelse(weighted, "Weighted", "Unweighted")
  title <- str_c("Champ on", w, name, "Network", sep = " ")
  
  # Champ visualizations
  
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
    geom_text(data = segments,
              aes(x = (x1+x2)/2, 
                  y = (y1+y2)/2, 
                  label = partitions),
              color = segments$partitions,
              vjust = -.5) +
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
         title = title) +
    theme_few() +
    theme(axis.text = element_text(size = 8))
  
  ggsave("figures/figure1.png", width = 12)
  
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
    geom_text(data = segments,
              aes(x = (x1+x2)/2, 
                  y = (y1+y2)/2, 
                  label = partitions),
              color = segments$partitions,
              vjust = -.5) +
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
         title = title)+
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
  
  ggsave("figures/figure2.png", width = 12)
  
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
    geom_text(aes(x = (x1+x2)/2, 
                  y = (y1+y2)/2, 
                  label = partitions),
              color = segments$partitions,
              vjust = -.5) +
    guides(color = "none") +
    labs(x = expression(paste("Resolution Parameter (", gamma,")")),
         y = "Modularity",
         title = title)+
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

  
  ggsave("figures/figure3.png", width = 12)
  
  # Partition summary
  # Ordering line segments by length
  partition_summary <- data.frame(matrix(ncol = 5, nrow = length(segments[,1])-1))
  colnames(partition_summary) <- c("segment_length", "starting_gamma", "ending_gamma", "gamma_range", "partition_num")
  
  for (x in 1:lengths(partition_summary)[1]) {
    partition_summary[x,"segment_length"] <- sqrt((segments[x, "x1"]-segments[x, "x2"])**2+(segments[x, "y1"]-segments[x, "y2"])**2)
    partition_summary[x,"starting_gamma"] <- segments[x,"x1"]
    partition_summary[x,"ending_gamma"] <- segments[x,"x2"]
    partition_summary[x,"gamma_range"] <- abs(segments[x,"x1"]-segments[x,"x2"])
    partition_summary[x,"partition_num"] <- segments[x,"partitions"]
    partition_summary[x,"num_clusters"] <- best_partitions$partitions[[x]]$nb_clusters
  }
  
  partition_summary <- partition_summary[order(-partition_summary$gamma_range),]
  print(partition_summary)
  if (report != "") {
    print(report)
  }
  
  # bump ideanet channel with github
  
  # Top n_plots partitions
  top <- partition_summary[1:n_plots,]
  top_partitions <- top[,"partition_num"]
  top_clusters <- top[,"num_clusters"]
  for (i in 1:n_plots) {
    partitions <- best_partitions$partitions
    p <- top_partitions[[i]]
    png(str_c("figures/","partition",top_partitions[i],".png"),
        width = 680, 
        height = 680, 
        units = "px",
        res = 100)
    plot(partitions[[p]], network, main = str_c(name,
                                                "Partition",
                                                top_partitions[i],
                                                "with",
                                                top_clusters[i],
                                                "clusters",
                                                sep = " "))
    dev.off()
  }
  return(list("partitions" = best_partitions$partitions,
              "summary" = partition_summary,
              "warning" = report))
}
