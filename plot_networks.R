library(ggplot2)

evaluate <- "16S"

# ============================================================================================================
#   PARSE 16S CORRELATIONS
# ============================================================================================================

CON <- tryCatch({
  read.table(paste0("correlations_",evaluate,"_CON.txt"), sep="\t", header=FALSE)
}, error = function(e) {
  # empty file
  data.frame(V1=NA, V2=NA, V3=NA)
})
ABX <- tryCatch({
  read.table(paste0("correlations_",evaluate,"_ABX.txt"), sep="\t", header=FALSE)
}, error = function(e) {
  # empty file
  data.frame(V1=NA, V2=NA, V3=NA)
})
ABXFT <- tryCatch({
  read.table(paste0("correlations_",evaluate,"_ABXFT.txt"), sep="\t", header=FALSE)
}, error = function(e) {
  # empty file
  data.frame(V1=NA, V2=NA, V3=NA)
})

correlators <- cbind(as.data.frame(CON), treatment="CON")
correlators <- rbind(correlators, cbind(as.data.frame(ABX), treatment="ABX"))
correlators <- rbind(correlators, cbind(as.data.frame(ABXFT), treatment="ABXFT"))
colnames(correlators) <- c("taxon1", "taxon1type", "tax1name", "taxon2", "tax2type", "tax2name", "value", "treatment")
correlators$taxon1 <- as.character(correlators$taxon1)
correlators$taxon2 <- as.character(correlators$taxon2)
correlators$treatment <- as.character(correlators$treatment)

correlators <- correlators[complete.cases(correlators),]

items <- sort(unique(c(correlators$taxon1, correlators$taxon2)))

# ============================================================================================================
#   ASSIGN PLOT POSITIONS (POINTS AROUND A CIRCLE OF RADIUS r)
# ============================================================================================================

# item_positions <- data.frame(item=c(), x=c(), y=c())
n <- length(items)
r <- 15
x_scale <- 0.4
y_scale <- 0.4
position_list <- list(x = c(), y = c())
for(i in 1:n) {
  x <- r * cos(2 * pi * i / n)
  y <- r * sin(2 * pi * i / n)
  # add item label positions; these need to be oriented such that they'll not overlap the points -- a pain!
  # x_label <- x
  # y_label <- y
  # label_len <- str_length(items[i])
  # if(y == r) {
  #   y_label <- y_label + r*y_scale
  # } else if(y == -r) {
  #   y_label <- y_label - r*y_scale
  # } else {
  #   if(x <= 0) {
  #     x_label <- x_label - label_len*x_scale
  #   } else {
  #     x_label <- x_label + label_len*x_scale
  #   }
  # }
  #item_positions <- rbind(item_positions, data.frame(item=items[i], x=x, y=y, x_label=x_label, y_label=y_label))
  position_list$x[i] <- -x
  position_list$y[i] <- y
}
# find top element
idx <- which(position_list$y == max(position_list$y))
rotate_order <- 1:n
if(idx > 1) {
  rotate_order <- c(idx:n, 1:(idx-1))
}

item_positions <- data.frame(item = items, x = position_list$x[rotate_order], y = position_list$y[rotate_order])

# ============================================================================================================
#   UPDATE CORRELATION EDGES NOW THAT WE HAVE POSITIONS
# ============================================================================================================

path_positions <- data.frame(x1=c(), y1=c(), x2=c(), y2=c(), weight=c(), sign=c())
for(i in 1:nrow(correlators)) {
  entry <- correlators[i,]
  path_positions <- rbind(path_positions, data.frame(x1=item_positions[item_positions == entry$taxon1,]$x,
                                                     y1=item_positions[item_positions == entry$taxon1,]$y,
                                                     x2=item_positions[item_positions == entry$taxon2,]$x,
                                                     y2=item_positions[item_positions == entry$taxon2,]$y,
                                                     weight=abs(entry$value)*8,
                                                     sign=sign(entry$value),
                                                     treatment=entry$treatment))
}
path_positions$sign <- as.factor(path_positions$sign)

item_positions$index <- 1:nrow(item_positions)

# ============================================================================================================
#   PLOT WITH LABELS
# ============================================================================================================

plot_network <- function(treatment, which_sign = "positive") {
  if(which_sign == "positive") {
    path_positions_filtered <- path_positions[as.numeric(path_positions$sign) == 2,]
    path_positions_filtered <- path_positions_filtered[path_positions_filtered$treatment == treatment,]
    color <- "#3CB44B"
  } else {
    path_positions_filtered <- path_positions[as.numeric(path_positions$sign) == 1,]
    path_positions_filtered <- path_positions_filtered[path_positions_filtered$treatment == treatment,]
    color <- "#E6194B"
  }
  p <- ggplot() +
    geom_segment(data=path_positions_filtered, aes(x=x1, y=y1, xend=x2, yend=y2, color=sign, size=weight), alpha=0.66) +
    scale_color_manual(values=color) +
    scale_size_identity() + # use the width specified by `weight`
    geom_point(data=item_positions, aes(x=x, y=y), size=5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
    #geom_text(data=item_positions, aes(x=x_label, y=y_label, label=item), size=4) +
    geom_text(data=item_positions, aes(x=x, y=y, label=index), size=3, color="#FFFFFF") +
    xlim(-r*3, r*3) +
    ylim(-r*3, r*3)
  p
  ggsave(paste0("images/",evaluate,"_correlation_network_",treatment,"_",which_sign,".png"), p, units="in", dpi=150, height=10, width=10)
  
  # work through the labels here, which are faulty
  if(FALSE) {
    sink(paste0("images/",treatment,"_network_labels.txt"))
    for(i in 1:nrow(item_positions)) {
      if(i < 10) {
        cat(" ")
      }
      item <- item_positions$item[i]
      correlators_subset <- correlators[correlators$treatment == treatment,]
      label <- correlators_subset[correlators_subset$taxon1 == item,]$tax1name
      if(is.null(label)) {
        label <- correlators[correlators$taxon2 == item_positions$item[i],]$tax2name
      }
      cat(paste0(i,"  ",as.character(label),"\n"))
    }
    sink()
  }
}

plot_network("CON", which_sign = "positive")
plot_network("CON", which_sign = "negative")
plot_network("ABX", which_sign = "positive")
plot_network("ABX", which_sign = "negative")
plot_network("ABXFT", which_sign = "positive")
plot_network("ABXFT", which_sign = "negative")

























