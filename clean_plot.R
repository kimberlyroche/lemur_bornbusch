library(ggplot2)

rm(list = ls())

plot_network <- function(input_file, which_sign = NULL, show_plot = FALSE) {

  # ============================================================================================================
  #   parse data
  # ============================================================================================================
  
  correlators <- tryCatch({
    read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    # empty file
    data.frame(feature_1 = NA, feature_2 = NA, interaction_strength = NA)
  })

  # grab the unique features in this table, e.g.
  # feature_1    feature_2    interaction_strength
  #         A            B    ...
  #         A            C    ...
  #         B            C    ...
  #
  # gives c(A, B, C)
  items <- sort(unique(c(correlators$feature_1, correlators$feature_2)))
  n_items <- length(items)
  n_interactions <- nrow(correlators)
  
  # ============================================================================================================
  #   assign plot node positions (points around a circle of radius `radius`)
  # ============================================================================================================
  
  radius <- 10 # this doesn't really matter
  
  nodes <- data.frame(item = character(n_items), x = numeric(n_items), y = numeric(n_items),
                      stringsAsFactors = FALSE)
  for(i in 1:n_items) {
    x <- radius * cos(2 * pi * i / n_items)
    y <- radius * sin(2 * pi * i / n_items)
    nodes$item[i] <- items[i]
    nodes$x[i] <- -x
    nodes$y[i] <- y
  }
  nodes$label <- 1:n_items

  # find the (vertically) topmost element, we'll call this no. 1
  topmost_idx <- which(nodes$y == max(nodes$y))
  if(topmost_idx > 1) {
    nodes$label <- c(seq(from = (n_items-topmost_idx+2), to = n_items), seq(from = 1,
                                                                            to = (n_items - topmost_idx + 1)))
  }

  # ============================================================================================================
  #   now that we have positions, label the interaction edges
  # ============================================================================================================
  
  edges <- data.frame(x1 = numeric(n_interactions), y1 = numeric(n_interactions),
                      x2 = numeric(n_interactions), y2 = numeric(n_interactions),
                      weight = numeric(n_interactions), sign = numeric(n_interactions))
  for(i in 1:n_interactions) {
    # string matching could be a problem here; in which case we'd just want unique numeric indices for
    # "features"
    edges$x1[i] <- nodes[nodes$item == correlators[i,]$feature_1,]$x
    edges$y1[i] <- nodes[nodes$item == correlators[i,]$feature_1,]$y
    edges$x2[i] <- nodes[nodes$item == correlators[i,]$feature_2,]$x
    edges$y2[i] <- nodes[nodes$item == correlators[i,]$feature_2,]$y
    edges$weight[i] <- abs(correlators[i,]$interaction_strength)*8 # the scalar here is just because this is the easiest way to adjust
    # edge thickness in the plot
    edges$sign[i] <- sign(correlators[i,]$interaction_strength)
  }

  # ============================================================================================================
  #   PLOT WITH LABELS
  # ============================================================================================================

  colors <- c("#E6194B", # green
              "#3CB44B") # red
  if(!is.null(which_sign)) {
    if(which_sign == "negative") {
      edges <- edges[edges$sign < 0,]
      colors <- colors[1]
    } else {
      edges <- edges[edges$sign > 0,]
      colors <- colors[2]
    }
  }
  
  p <- ggplot() +
    geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2, color = as.factor(sign), size = weight), alpha = 0.66) +
    scale_color_manual(values = colors) +
    scale_size_identity() + # use the width specified by `weight`
    geom_point(data = nodes, aes(x = x, y = y), size = 7) +
    geom_text(data = nodes, aes(x = x, y = y, label = label), size = 4, color = "#FFFFFF") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none") +
    xlim(-radius*1.2, radius*4) +
    ylim(-radius*1.2, radius*1.2)
  
  # crudely work out some spacing
  legend_tl_x <- radius
  legend_tl_y <- max(nodes$y)
  item_space_x <- radius*0.08 # this is the horizontal distance from bullet to feature name in the legend
  item_space_y <- radius*0.08 # this is the vertical distance between bullets in the legend

  reorder.order <- order(nodes$label)
  legend_df <- data.frame(item = nodes$item[reorder.order], label = nodes$label[reorder.order])
  legend_df$x_label <- radius*1.5
  legend_df$y_label <- legend_tl_y - seq(from = 0, to = n_items*item_space_y, length.out = n_items)
  legend_df$x_item <- legend_df$x_label + item_space_x
  legend_df$y_item <- legend_df$y_label

  # manually plot legend
  p <- p +
    geom_point(data = legend_df, aes(x = x_label, y = y_label), size = 7) + 
    geom_text(data = legend_df, aes(x = x_label, y = y_label, label = label), size = 4, color = "#FFFFFF") +
    geom_text(data = legend_df, aes(x = x_item, y = y_item, label = item), size = 4, hjust = 0)
  if(show_plot) {
    show(p)
  }
  ggsave(paste0(which_sign,"_interactions.png"), p, units="in", dpi=100, height=10, width=20)
}

plot_network("network_input.txt", which_sign = NULL, show_plot = TRUE)























