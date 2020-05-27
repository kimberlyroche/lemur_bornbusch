rm(list = ls())

library(tidyverse)
library(driver)
library(stray)
library(Rcpp)
library(abind)

# uncollapsing/posterior sampling of Lambda, Sigma; lifted directly from stray::basset implementations
sourceCpp("sampling.cpp")

# =====================================================================================================
#   EXPLORATORY FNS.
# =====================================================================================================

# plot the change in (1) the mean log abundance and (2) very strict interquartile mean to get a sense of
# how large the low-variance cohort is
explore_mean_stability <- function(data, lower.quantile = 0.25, upper.quantile = 0.75) {
  # I'm combining years though; I think that's ok, though the days may only roughtly line up
  df <- NULL
  for(treatment in c("CON", "ABX", "ABXFT")) {
    for(year in c(1,2)) {
      use_samples <- data$`16S`[[paste0("year",year)]]$metadata$Treatment == treatment
      use_md <- data$`16S`[[paste0("year",year)]]$metadata[use_samples,c("Description","Day")]
      use_counts <- data$`16S`[[paste0("year",year)]]$filtered[,use_samples]
      log_counts <- log(use_counts + 0.5)
      df.temp <- data.frame(logmean = apply(log_counts, 2, mean), treatment = treatment, timepoint = use_md$Day)
      if(is.null(df)) {
        df <- df.temp
      } else {
        df <- rbind(df, df.temp)
      }
    }
  }
  p <- ggplot(df) +
    geom_point(aes(x = timepoint, y = logmean, color = treatment))
  ggsave("images/change_in_log_mean.png", p, dpi = 100, units = "in", height = 6, width = 10)
  
  # this is somewhat mitigated by the IQLR but you have to choose very tight quartiles!
  df <- NULL
  for(treatment in c("CON", "ABX", "ABXFT")) {
    for(year in c(1,2)) {
      use_samples <- data$`16S`[[paste0("year",year)]]$metadata$Treatment == treatment
      use_md <- data$`16S`[[paste0("year",year)]]$metadata[use_samples,c("Description","Day")]
      use_counts <- data$`16S`[[paste0("year",year)]]$filtered[,use_samples]
      # evaluate the CLR variance (variance relative to the mean for each taxon)
      clr.counts <- t(clr(t(use_counts) + 0.5)) # taxa x samples
      clr.var <- apply(clr.counts, 1, var)
      # chose the most typical
      qs <- quantile(clr.var, probs = c(lower.quantile, upper.quantile))
      include_taxa <- clr.var > qs[1] & clr.var < qs[2]
      log_counts <- log(use_counts[include_taxa,] + 0.5)
      df.temp <- data.frame(logmean = apply(log_counts, 2, mean), treatment = treatment, timepoint = use_md$Day)
      if(is.null(df)) {
        df <- df.temp
      } else {
        df <- rbind(df, df.temp)
      }
    }
  }
  p <- ggplot(df) +
    geom_point(aes(x = timepoint, y = logmean, color = treatment))
  ggsave("images/change_in_IQL_mean.png", p, dpi = 100, units = "in", height = 6, width = 10)
}

# =====================================================================================================
#   DATA PARSING
# =====================================================================================================

# parse out the taxonomic labels from either 16S or ITS count table
get_counts <- function(df) {
  max_tax_id <- max(which(sapply(df[1,], is.numeric) == FALSE))
  return(df[,(max_tax_id+1):ncol(df)])
}

get_exclude_sample_list <- function(data_type) {
  if(data_type == "16S") {
    exclude_samples <- list(year1=c("LCAXC1","LCAXC2"),
                            year2=c("CON1", "CON3", "LCAX038",
                                    "LCAX050", "LCAX071", "LCAX073", "LCAX114", "LCAX130",
                                    "LCAX137", "LCAX141", "LCAX142", "LCAX149", "LCAX152",
                                    "LCAX157", "LCAX158", "LCAX173"))
  } else {
    exclude_samples <- list(year1=c("CON1", "CON2", "CON5", "LCAXC1"),
                            year2=c("CON3", "LCAX038",
                                    "LCAX050", "LCAX071", "LCAX073", "LCAX114", "LCAX130",
                                    "LCAX137", "LCAX141", "LCAX142", "LCAX149", "LCAX152",
                                    "LCAX157", "LCAX158", "LCAX173"))
  }
  return(exclude_samples)
}

# parse sequence count data AND consolidate ASVs between years
# arguments:
#   data_type: "16S" or "ITS"
# returns: a pair of data.frames of taxa (rows) x phylogeny (columns 1:7) + samples (columns 8:K)
read_sequence_counts <- function(data_type = "16S") {
  # excluded samples
  exclude_samples <- get_exclude_sample_list(data_type)

  if(data_type == "16S") {
    sequences.y1 <- read.csv(paste0("data/Bornbusch_16S_Y1_fecal_count-table.tsv"), header=TRUE, sep='\t')
  } else {
    sequences.y1 <- read.csv(paste0("data/Bornbusch_ITS_Y1_fecal_count-table.tsv"), header=TRUE, sep='\t')
  }
  
  rownames(sequences.y1) <- sequences.y1$OTU.ID
  tax.y1 <- as.list(sequences.y1$TAX)
  names(tax.y1) <- rownames(sequences.y1)
  sequences.y1 <- sequences.y1[,!(colnames(sequences.y1) %in% c("OTU.ID", "TAX", "CONF", exclude_samples$year1))]
  
  if(data_type == "16S") {
    sequences.y2 <- read.csv(paste0("data/Bornbusch_16S_Y2_fecal_count-table.tsv"), header=TRUE, sep='\t')
  } else {
    sequences.y2 <- read.csv(paste0("data/Bornbusch_ITS_Y2_fecal_count-table.tsv"), header=TRUE, sep='\t')
  }
  
  rownames(sequences.y2) <- sequences.y2$OTU.ID
  tax.y2 <- as.list(sequences.y2$TAX)
  names(tax.y2) <- rownames(sequences.y2)
  sequences.y2 <- sequences.y2[,!(colnames(sequences.y2) %in% c("OTU.ID", "TAX", "CONF", exclude_samples$year2))]
  
  # consolidate sequence variant ID's between years
  # this is time-consuming, so save the results when complete
  if(data_type == "16S" & file.exists("data/tax_collapsed_y1_16S.rds") & file.exists("data/tax_collapsed_y1_16S.rds")) {
    tax.all <- readRDS("data/tax_all_16S.rds")
    tax.collapsed.y1 <- readRDS("data/tax_collapsed_y1_16S.rds")
    tax.collapsed.y2 <- readRDS("data/tax_collapsed_y2_16S.rds")
  } else if(data_type == "ITS" & file.exists("data/tax_collapsed_y1_ITS.rds") & file.exists("data/tax_collapsed_y1_ITS.rds")) {
    tax.all <- readRDS("data/tax_all_ITS.rds")
    tax.collapsed.y1 <- readRDS("data/tax_collapsed_y1_ITS.rds")
    tax.collapsed.y2 <- readRDS("data/tax_collapsed_y2_ITS.rds")
  } else {
    # coerce all sequence variant ID's to match
    addend.y1 <- rownames(sequences.y2)[!(rownames(sequences.y2) %in% rownames(sequences.y1))]
    addend.mat.y1 <- as.data.frame(matrix(0, length(addend.y1), ncol(sequences.y1)))
    rownames(addend.mat.y1) <- addend.y1
    colnames(addend.mat.y1) <- colnames(sequences.y1)
    sequences.y1 <- rbind(sequences.y1, addend.mat.y1)
    
    OTU_list.y1 <- rownames(sequences.y1)
    sample_list.y1 <- colnames(sequences.y1)
    
    addend.y2 <- rownames(sequences.y1)[!(rownames(sequences.y1) %in% rownames(sequences.y2))]
    addend.mat.y2 <- as.data.frame(matrix(0, length(addend.y2), ncol(sequences.y2)))
    rownames(addend.mat.y2) <- addend.y2
    colnames(addend.mat.y2) <- colnames(sequences.y2)
    sequences.y2 <- rbind(sequences.y2, addend.mat.y2)
    
    OTU_list.y2 <- rownames(sequences.y2)
    sample_list.y2 <- colnames(sequences.y2)
    
    # get OTU's in both data.frames in the same order
    OTU_resort.y1 <- order(OTU_list.y1)
    OTU_list.y1 <- OTU_list.y1[OTU_resort.y1]
    sequences.y1 <- sequences.y1[OTU_resort.y1,]
    
    OTU_resort.y2 <- order(OTU_list.y2)
    OTU_list.y2 <- OTU_list.y2[OTU_resort.y2]
    sequences.y2 <- sequences.y2[OTU_resort.y2,]
    
    counts.y1 <- driver::gather_array(sequences.y1, "value", "OTU", "sample")
    counts.y1 <- cbind(counts.y1, year=1)
    counts.y1$OTU <- as.factor(counts.y1$OTU)
    levels(counts.y1$OTU) <- OTU_list.y1
    counts.y1$sample <- as.factor(counts.y1$sample)
    levels(counts.y1$sample) <- sample_list.y1
    
    counts.y2 <- driver::gather_array(sequences.y2, "value", "OTU", "sample")
    counts.y2 <- cbind(counts.y2, year=2)
    counts.y2$OTU <- as.factor(counts.y2$OTU)
    levels(counts.y2$OTU) <- OTU_list.y2
    counts.y2$sample <- as.factor(counts.y2$sample)
    levels(counts.y2$sample) <- sample_list.y2
    
    counts.all <- rbind(counts.y1, counts.y2)
    
    OTU_list.all <- rownames(sequences.y1)
    
    # merge taxonomic lookups!
    keys <- unique(c(names(tax.y1), names(tax.y2)))
    tax.all <- list()
    for(k in 1:length(keys)) {
      OTU_ID <- keys[k]
      OTU_str <- tax.y1[[OTU_ID]]
      if(is.null(OTU_str)) {
        OTU_str <- tax.y2[[OTU_ID]]
      }
      if(is.null(OTU_str)) {
        cat("Not found in either list?!\n")
      }
      tax.all[[OTU_ID]] <- OTU_str
    }
    
    tax.df <- cbind(counts.all, data.frame(kingdom=NA, phylum=NA, class=NA, order=NA, family=NA, genus=NA, species=NA))
    #for(i in 1:length(OTU_list.all)) {
    for(i in 1:30) {
      if(i %% 10 == 0) {
        cat("Parsing OTU:",i,"\n")
      }
      OTU_ID <- OTU_list.all[i]
      if(data_type == "16S") {
        str_pieces <- strsplit(as.character(tax.all[[OTU_ID]]), "(;?)D_\\d__", perl=TRUE)[[1]]
      } else {
        str_pieces <- strsplit(as.character(tax.all[[OTU_ID]]), "(;?).__", perl=TRUE)[[1]]
      }
      str_pieces <- str_pieces[2:length(str_pieces)]
      if(length(str_pieces) < 7) {
        str_pieces <- c(str_pieces, rep(NA, 7-length(str_pieces)))
      }
      # the order seems to get fucked up when I assign these as a vector (!!!)
      assign_list <- which(tax.df$OTU == OTU_ID)
      tax.df[assign_list,]$kingdom <- str_pieces[1]
      tax.df[assign_list,]$phylum <- str_pieces[2]
      tax.df[assign_list,]$class <- str_pieces[3]
      tax.df[assign_list,]$order <- str_pieces[4]
      tax.df[assign_list,]$family <- str_pieces[5]
      tax.df[assign_list,]$genus <- str_pieces[6]
      tax.df[assign_list,]$species <- str_pieces[7]
    }
    
    tax.df <- tax.df[!is.na(tax.df$kingdom),]
    
    tax.collapsed <- tax.df %>%
      group_by(kingdom, phylum, class, order, family, genus, species, sample) %>%
      summarise(agg_val = sum(value)) %>%
      spread(sample, agg_val)
    
    # omit a few ; this is very few counts (23 / > 100 million)
    tax.collapsed <- tax.collapsed[which(tax.collapsed$phylum != "Unassigned"),]
    if(data_type == "16S") {
      tax.collapsed <- tax.collapsed[which(tax.collapsed$order != "Chloroplast"),]
    }
    
    # separate the years back out
    tax.collapsed.y1 <- as.data.frame(tax.collapsed)
    tax.collapsed.y1 <- tax.collapsed.y1[,colnames(tax.collapsed.y1) %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species", sample_list.y1)]
    tax.collapsed.y2 <- as.data.frame(tax.collapsed)
    tax.collapsed.y2 <- tax.collapsed.y2[,colnames(tax.collapsed.y2) %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species", sample_list.y2)]
    
    append_str <- data_type
    saveRDS(tax.collapsed, paste0("data/tax_all_",append_str,".rds"))
    saveRDS(tax.collapsed.y1, paste0("data/tax_collapsed_y1_",append_str,".rds"))
    saveRDS(tax.collapsed.y2, paste0("data/tax_collapsed_y2_",append_str,".rds"))
  }

  # print stats  
  n_taxa <- nrow(tax.collapsed.y1)
  counts.all <- cbind(get_counts(tax.collapsed.y1), get_counts(tax.collapsed.y2))
  p_zero <- sum(counts.all == 0)/(nrow(counts.all) * ncol(counts.all))
  cat(paste0("Agglomerated data has ",n_taxa," taxa and is ",round(p_zero*100),"% zeros!\n"))
  return(list(year1 = tax.collapsed.y1, year2 = tax.collapsed.y2))
}

# read metadata per year
# returns: a pair of data.frames of samples (rows) x attributes (columns)
read_metadata <- function(data_type = "16S") {
  metadata <- list(year1 = NULL, year2 = NULL)
  exclude_samples <- get_exclude_sample_list(data_type)
  for(year in c(1, 2)) {
    year_label <- paste0("Y",year)
    if(data_type == "16S") {
      metadata.y <- read.csv(paste0("data/16S_",year_label,"_fecal_metadata.txt"), header=TRUE, sep='\t', stringsAsFactors = FALSE)
    } else {
      metadata.y <- read.csv(paste0("data/ITS_",year_label,"_fecal_metadata.txt"), header=TRUE, sep='\t', stringsAsFactors = FALSE)
    }
    
    metadata.y$SampleID <- as.character(metadata.y$SampleID)
    metadata.y <- metadata.y[!(metadata.y$SampleID %in% exclude_samples[[paste0("year",year)]]),]

    # render periods as "Pre", "During", "Post"
    metadata.y$Period <- as.character(metadata.y$Period)
    metadata.y$Period <- unname(sapply(metadata.y$Period, function(x) {
      str_replace(x, "\\d+", "")
    }))
    
    # alter year labels: "Y1" -> 1, "Y2" -> 2
    metadata.y$Year <- unname(sapply(metadata.y$Year, function(x) {
      str_replace(x, "Y", "")
    }))
    metadata.y$Year <- as.numeric(metadata.y$Year)
    
    # standardize names
    metadata.y$Animal <- as.character(metadata.y$Animal)
    metadata.y$Animal <- sapply(metadata.y$Animal, function(x) sub(".*Nikos.*", "Nikos", x, perl=TRUE))
    metadata.y$Animal <- sapply(metadata.y$Animal, function(x) sub(".*Jones.*", "Jones", x, perl=TRUE))
    # metadata.y$Animal <- as.factor(metadata.y$Animal)
    # Nikos received antibiotics, move him from "CON" to "ABX" treatment groups
    # this is indicated in Experimental_Groups.xlsx
    # metadata.y[metadata.y$Animal == "Nikos",]$Treatment <- "ABX"
    # fix this awful space that prevents string matching
    # new_levels <- levels(metadata.y$Animal)
    # new_levels[which(new_levels == "Randy ")] <- "Randy"
    # levels(metadata.y$Animal) <- new_levels
    
    # fix Randy_
    metadata.y$Animal <- unname(sapply(metadata.y$Animal, function(x) {
      str_replace(x, "Randy ", "Randy")
    }))
    
    metadata[[paste0("year",year)]] <- metadata.y
  }
  return(metadata)
}

# =====================================================================================================
#   DATA AGGLOMERATION AND FILTERING
# =====================================================================================================

get_tax_labels <- function() {
  return(c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
}

filter_data <- function(raw_data, metadata, level = "genus", filter_percent = NULL) {
  # NOTE: things undefined at the specified taxonomic level will be filtered into "Other"
  tax_pieces <- get_tax_labels()
  tax_pieces <- tax_pieces[1:which(tax_pieces == level)]
  tax_pieces <- c(tax_pieces, "experiment")

  # look for all-zero-count taxa to remove first
  all_zero_taxa.y1 <- apply(as.matrix(get_counts(raw_data$year1)), 1, function(x) sum(x) == 0)
  all_zero_taxa.y2 <- apply(as.matrix(get_counts(raw_data$year2)), 1, function(x) sum(x) == 0)
  all_zero_taxa <- all_zero_taxa.y1 & all_zero_taxa.y2
  raw_data$year1 <- raw_data$year1[!all_zero_taxa,]
  raw_data$year2 <- raw_data$year2[!all_zero_taxa,]
  
  # print stats
  n_taxa <- nrow(raw_data$year1)
  counts.all <- cbind(get_counts(raw_data$year1), get_counts(raw_data$year2))
  p_zero <- sum(counts.all == 0)/(nrow(counts.all) * ncol(counts.all))
  cat(paste0("After filtering out all-zero taxa, data has ",n_taxa," taxa and is ",round(p_zero*100),"% zeros!\n"))
  
  # filter out taxa that are present < k times in the entire data set
  # I'm seeing lots of singletons or near singletons
  # arguably we can't measure correlation on these anyway, so we'll exclude them
  # setting k = 5 (minimum appearance count) retains over 99% of total counts!
  filter_rule <- function(x) { sum(x != 0) < 5 }
  vrare_taxa.y1 <- apply(as.matrix(get_counts(raw_data$year1)), 1, filter_rule)
  vrare_taxa.y2 <- apply(as.matrix(get_counts(raw_data$year2)), 1, filter_rule)
  vrare_taxa <- vrare_taxa.y1 | vrare_taxa.y2
  raw_data$year1 <- raw_data$year1[!vrare_taxa,]
  raw_data$year2 <- raw_data$year2[!vrare_taxa,]

  # print stats
  n_taxa <- nrow(raw_data$year1)
  counts.all <- cbind(get_counts(raw_data$year1), get_counts(raw_data$year2))
  p_zero <- sum(counts.all == 0)/(nrow(counts.all) * ncol(counts.all))
  cat(paste0("After filtering out taxa with very low presence, data has ",n_taxa," taxa and is ",round(p_zero*100),"% zeros!\n"))
    
  # exclude (into "other") anything not resolved to the desired level
  # we'll incorporate these counts into Other later
  unresolved <- is.na(raw_data$year1[,c(level)]) | is.na(raw_data$year2[,c(level)])
  
  # package up the keepers
  min_sample_idx.y1 <- min(which(sapply(raw_data$year1[1,], is.numeric)))
  max_sample_idx.y1 <- max(which(sapply(raw_data$year1[1,], is.numeric)))
  min_sample_idx.y2 <- min(which(sapply(raw_data$year2[1,], is.numeric)))
  max_sample_idx.y2 <- max(which(sapply(raw_data$year2[1,], is.numeric)))
  collapsed_unresolved.y1 <- colSums(raw_data$year1[unresolved,min_sample_idx.y1:max_sample_idx.y1])
  collapsed_unresolved.y2 <- colSums(raw_data$year2[unresolved,min_sample_idx.y2:max_sample_idx.y2])
  resolved <- list(year1 = raw_data$year1,
                   year2 = raw_data$year2)
  resolved$year1 <- resolved$year1[!unresolved,]
  resolved$year2 <- resolved$year2[!unresolved,]
  
  # print stats
  n_taxa <- nrow(resolved$year1)
  counts.all <- cbind(get_counts(resolved$year1), get_counts(resolved$year1))
  p_zero <- sum(counts.all == 0)/(nrow(counts.all) * ncol(counts.all))
  cat(paste0("After filtering out unresolved taxa, data has ",n_taxa," taxa and is ",round(p_zero*100),"% zeros!\n"))
  
  agg <- function(data, min_idx, max_idx, tax_pieces) {
    suppressWarnings(as.data.frame(data %>%
                    gather(experiment, value, min_idx:max_idx) %>%
                    group_by(!!!syms(tax_pieces)) %>%
                    summarise(agg_val = sum(value)) %>%
                    spread(experiment, agg_val)))
  }
  collected <- list(year1 = agg(resolved$year1, min_sample_idx.y1, max_sample_idx.y1, tax_pieces),
                       year2 = agg(resolved$year2, min_sample_idx.y2, max_sample_idx.y2, tax_pieces))
  all_collected <- cbind(collected$year1, get_counts(collected$year2))
  
  # now, looking across both years, exclude any taxa that are below a minimum total abundance
  # this and the first (near-singleton) filtering were implemented at different stages and should probably
  # be combined at this point
  counts <- get_counts(all_collected) # data.frame of taxa (rows) x samples (columns)
  level_idx <- length(tax_pieces) - 1 # translate genus=6, species=7, etc.
  
  filter_rule <- function(df, filter_percent) { (apply(df, 1, sum)/sum(df)) > (filter_percent/100) }
  retain_idx <- filter_rule(counts, filter_percent)
  retained_tax <- all_collected[retain_idx, 1:level_idx]
  # collapse to create "other" category and tack it on as the last row in the data.frame
  all_filtered <- counts[retain_idx,]
  all_filtered[nrow(all_filtered)+1,] <- apply(counts[!retain_idx,], 2, sum) + c(collapsed_unresolved.y1, collapsed_unresolved.y2)
  retained_tax[nrow(retained_tax)+1,] <- "Other"
  
  # print stats
  n_taxa <- nrow(all_filtered)
  p_zero <- sum(all_filtered == 0)/(nrow(all_filtered) * ncol(all_filtered))
  cat(paste0("After filtering very low abundance taxa into 'other', data has ",n_taxa," taxa and is ",round(p_zero*100),"% zeros!\n"))
  
  return(list(counts = all_filtered, level = level_idx, tax = retained_tax))
}

# we'll build this data structure as a list of lists
#
# all_data
#  |
#  --- tax
#  |
#  --- 16S
#  |    |
#  |    --- year1
#  |    |    |
#  |    |    --- filtered
#  |    |    |
#  |    |    --- metadata
#  |    |
#  |    --- year2
#  |         |
#  |         --- filtered
#  |         |
#  |         --- metadata
#  |
#  --- ITS
#       |
#       --- year1
#       |    |
#       |    --- filtered
#       |    |
#       |    --- metadata
#       |
#       --- year2
#            |
#            --- filtered
#            |
#            --- metadata
#
#   filtered -- taxa (rows) x {phylogeny, samples} (columns) as kingdom ... species LCAX001 LCAX002 ... LCAX104
#   metadata -- samples (rows) x annotations (columns) e.g. SampleID, Animal, Year, Date, Day, Treatment
#
# if data_type == NULL, pull 16S and ITS
get_data <- function(data_type = NULL, remove_baseline_samples = TRUE) {
  all_data <- list()
  data_types <- data_type
  if(is.null(data_type)) {
    data_types <- c("16S", "ITS")
  }
  for(data_type in data_types) {
    raw_data <- read_sequence_counts(data_type = data_type)
    metadata <- read_metadata(data_type = data_type)
    tax_level <- "genus"
    tax_labels <- get_tax_labels()
    tax_labels <- tax_labels[1:which(tax_labels == tax_level)]
    filtered_data <- filter_data(raw_data, metadata, level = tax_level, filter_percent = 0.01)
    
    # split the data out into separate years again
    all_data[[data_type]] <- list(tax = filtered_data$tax,
                                  year1 = list(filtered = filtered_data$counts[,colnames(filtered_data$counts) %in% metadata$year1$SampleID],
                                           metadata = metadata$year1),
                                  year2 = list(filtered = filtered_data$counts[,colnames(filtered_data$counts) %in%  metadata$year2$SampleID],
                                           metadata = metadata$year2))
    # str(all_data, max.level = 3)
    
    # test this: if in a treatment group, exclude before and after samples to see if signal increases
    if(remove_baseline_samples) {
      tag_baseline <- function(data) {
        include_samples <- !str_detect(data$metadata$Period, "Pre")
        # chop "Pre" samples out of metadata
        data$metadata <- data$metadata[include_samples,]
        # used filtered metadata sample IDs/descriptions
        data$filtered <- data$filtered[,colnames(data$filtered) %in% c(tax_labels, as.character(data$metadata$Description))]
        return(data)
      }
      all_data[[data_type]]$year1 <- tag_baseline(all_data[[data_type]]$year1)
      all_data[[data_type]]$year2 <- tag_baseline(all_data[[data_type]]$year2)
    }
  }
  return(all_data)
}

# match samples between 16S, ITS
match_sample_IDs_across_data_types <- function(all_data) {
  for(year in c(1,2)) {
    year_label <- paste0("year",year)
    shared_sample_IDs <- intersect(c(all_data$`16S`[[year_label]]$metadata$SampleID), c(all_data$ITS[[year_label]]$metadata$SampleID))
    for(data_type in c("16S", "ITS")) {
      year_label <- paste0("year",year)
      all_data[[data_type]][[year_label]]$metadata <- all_data[[data_type]][[year_label]]$metadata[all_data[[data_type]][[year_label]]$metadata$SampleID %in% shared_sample_IDs,]
      all_data[[data_type]][[year_label]]$filtered <- all_data[[data_type]][[year_label]]$filtered[,colnames(all_data[[data_type]][[year_label]]$filtered) %in% shared_sample_IDs]
    }
  }
  return(all_data)
}

# subsets the data (across years) given a list of sample IDs
subselect_from_sampleIDs <- function(data, target_sampleIDs) {
  for(year_label in c("year1","year2")) {
    for(data_type in names(data)) {
      data[[data_type]][[year_label]]$filtered <- data[[data_type]][[year_label]]$filtered[,target_sampleIDs[[year_label]]]
      if(ncol(data[[data_type]][[year_label]]$filtered) > 0) {
        data[[data_type]][[year_label]]$metadata <- suppressMessages(left_join(data.frame(SampleID = target_sampleIDs[[year_label]], stringsAsFactors = FALSE),
                                                                               data[[data_type]][[year_label]]$metadata))
      } else {
        # otherwise, empty selection; no metadata
        data[[data_type]][[year_label]]$metadata <- NULL
      }
    }
  }
  return(data)
}

# pulls metadata for a given label_type across years
# year : 1 or 2
get_label_type <- function(all_data, label_type, year = NULL) {
  arb_data_type <- names(all_data)[[1]]
  if(!is.null(year)) {
    labels <- all_data[[arb_data_type]][[paste0("year", year)]]$metadata[[label_type]]
  } else {
    # grab both years
    labels <- c()
    for(year in c(1,2)) {
      year_label <- paste0("year", year)
      labels <- c(labels, all_data[[arb_data_type]][[year_label]]$metadata[[label_type]])
    }
  }
  return(labels)
}

# given an individual animal, year, and reference animal ordering, calculates the X values
# for that animals observation days, e.g.
# animal "A" observed at 1, 3, 7, 8 will return bumped_X = 1, 3, 7, 8 and
#                                               predicted bumped_X = 1:8
# animal "B" observed at 1, 2, 5, 7 will return X = 1001, 1002, 1005, 1007
#                                               predicted bumped_X = 1001:1007
# X and bumped_X effectively serve as a mapping from Days in the metadata for a given
# individual to observed days as the appear in the design matrix
map_index <- function(data, animal, year, animal_order, bump = 1000, predict = FALSE) {
  day_labels <- get_label_type(data, "Day", year = year)
  animal_labels <- get_label_type(data, "Animal", year = year)
  animal_idx <- which(animal_order == animal)
  if(predict) {
    # interpolate spacing
    base_X <- min(day_labels[animal_labels == animal]):max(day_labels[animal_labels == animal])
  } else {
    # leave spaces
    base_X <- day_labels[animal_labels == animal]
  }
  year_offset <- 0
  if(year == 2) {
    # get year bump/offset; this is almost exact
    year_offset <- length(unique(get_label_type(data, "Animal", year = 1)))*bump
  }
  bumped_X <- base_X + bump*(animal_idx - 1) + year_offset
  return(list(X = base_X, bumped_X = bumped_X))
}

# reorder such that individual animal's observations are in blocks; this format is just
# easier to deal with down the road
reorder_samples <- function(data) {
  sample_IDs <- list()
  for(year in c(1,2)) {
    arb_data_type <- names(data)[[1]]
    day_labels <- get_label_type(data, "Day", year = year)
    animal_labels <- get_label_type(data, "Animal", year = year)
    sample_labels <- get_label_type(data, "SampleID", year = year)
    unique_animals <- unique(animal_labels)
    year_label <- paste0("year",year)
    sample_IDs[[year_label]] <- c()
    for(animal in unique_animals) {
      samples_for_animal_year <- sample_labels[animal_labels == animal]
      # there samples are not always in order of time
      days_for_subset <- unname(sapply(samples_for_animal_year, function(x) {
        day_labels[sample_labels == x]
      }))
      reordering <- order(days_for_subset)
      sample_IDs[[year_label]] <- c(sample_IDs[[year_label]], samples_for_animal_year)
    }
  }
  data <- subselect_from_sampleIDs(data, sample_IDs)
  return(data)
}

pull_treatment_data <- function(data, treatment) {
  sample_labels <- get_label_type(data, "SampleID")
  treatment_labels <- get_label_type(data, "Treatment")
  year_labels <- get_label_type(data, "Year")
  sample_IDs <- list()
  sample_IDs[["year1"]] <- sample_labels[treatment_labels == treatment & year_labels == 1]
  sample_IDs[["year2"]] <- sample_labels[treatment_labels == treatment & year_labels == 2]
  treatment_data <- subselect_from_sampleIDs(data, sample_IDs)
}

# note: this returns design matrices of characters; time will need to be converted to numeric later
build_design_matrix <- function(data, animal_order) {
  X_time <- c()
  X_time_predict <- c()
  X_animal_labels <- c()
  X_animal_labels_predict <- c()
  X_year <- c()
  X_year_predict <- c()
  for(year in c(1,2)) {
    for(animal in animal_order) {
      if(animal %in% unique(get_label_type(data, "Animal", year = year))) {
        # if an animal has samples in this treatment condition in this year
        cat("Adding animal",animal,"in year",year,"\n")
        X_animal <- map_index(data, animal, year, animal_order, predict = FALSE)$bumped_X
        X_animal_labels <- c(X_animal_labels, rep(animal, length(X_animal)))
        X_time <- c(X_time, X_animal)
        X_animal_predict <- map_index(data, animal, year, animal_order, predict = TRUE)$bumped_X
        X_animal_labels_predict <- c(X_animal_labels_predict, rep(animal, length(X_animal_predict)))
        X_time_predict <- c(X_time_predict, X_animal_predict)
        X_year <- c(X_year, rep(year, length(X_animal)))
        X_year_predict <- c(X_year_predict, rep(year, length(X_animal_predict)))
      }    
    }
  }
  X <- rbind(as.numeric(X_time), X_animal_labels, X_year)
  X_predict <- rbind(X_time_predict, X_animal_labels_predict, X_year_predict)
  return(list(X = X, X_predict = X_predict))
}

# multinomial-Dirichlet resample the count data set
# called by resample_dataset
# data set returned are CLR values with dimensions: {taxa (16S or ITS)} x {samples (Y1 + Y2)}
resample_data_type <- function(data, data_type) {
  counts <- cbind(data[[data_type]]$year1$filtered, data[[data_type]]$year2$filtered)
  clr.resampled <- matrix(NA, nrow(counts), ncol(counts))
  for(j in 1:ncol(counts)) {
    augmented_counts <- unlist(counts[,j] + 0.5)
    clr.resampled[,j] <- c(LaplacesDemon::rdirichlet(1, augmented_counts))
  }
  clr.resampled <- t(driver::clr(t(clr.resampled)))
  return(clr.resampled)
}

get_ntaxa <- function(data) {
  D.16S <- 0
  D.ITS <- 0
  if("16S" %in% names(data)) {
    if("year1" %in% names(data$`16S`)) {
      D.16S <- nrow(data$`16S`[["year1"]]$filtered)
    } else if("year1" %in% names(data$`16S`)) {
      D.16S <- nrow(data$`16S`[["year2"]]$filtered)
    }
  }
  if("ITS" %in% names(data)) {
    if("year1" %in% names(data$`16S`)) {
      D.ITS <- nrow(data$ITS[["year1"]]$filtered)
    } else if("year1" %in% names(data$`16S`)) {
      D.ITS <- nrow(data$ITS[["year2"]]$filtered)
    }
  }
  return(list(D.16S = D.16S, D.ITS = D.ITS))
}

resample_dataset <- function(data) {
  clr.resampled <- NULL
  if("16S" %in% names(data)) {
    clr.resampled <- resample_data_type(data, data_type = "16S")
  }
  if("ITS" %in% names(data)) {
    clr.resampled.ITS <- resample_data_type(data, data_type = "ITS")
    if(is.null(clr.resampled)) {
      clr.resampled <- clr.resampled.ITS
    } else {
      clr.resampled <- rbind(clr.resampled, clr.resampled.ITS)
    }
  }
  return(list(clr.resampled = clr.resampled))
}

fit_model <- function(data, X) {
  # resample
  D.taxa <- get_ntaxa(data)
  D.16S <- D.taxa$D.16S
  D.ITS <- D.taxa$D.ITS
  resample_results <- resample_dataset(data)
  clr.resampled <- resample_results$clr.resampled
  
  # hyperparameters
  dc <- 0.01 # desired minimum correlation
  dd_se <- 7
  rho_se <<- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
  g_sigma <<- 1.5
  
  upsilon <- (D.16S + D.ITS + 2) + 10 # low certainty
  Xi <- diag(D.16S + D.ITS)
  
  # make and assign baselines
  clr.baselines <- list() # named list of individual CLR mean abundances
  for(animal in unique(X[2,])) {
    replace_idx <- which(X[2,] == animal)
    clr.baselines[[animal]] <- rowMeans(clr.resampled[,replace_idx])
  }
  
  # clr.baselines must be in the global workspace
  Theta.train <- function(X) {
    Theta <- matrix(NA, length(clr.baselines[[1]]), ncol(X))
    for(j in 1:ncol(X)) {
      Theta[,j] <- clr.baselines[[X[2,j]]]
    }
    return(Theta)
  }
  
  Gamma <- function(X) {
    X_time <- as.numeric(X[1,])
    dim(X_time) <- c(1, length(X_time))
    SE(X_time, sigma=g_sigma, rho=rho_se, jitter=1e-10) # time only
  }
  
  posterior_samples <- uncollapse_simple(clr.resampled,
                                         diag(ncol(X)), # per fit_basset.R: 53
                                         Theta.train(X),
                                         Gamma(X),
                                         Xi,
                                         upsilon)
  
  return(list(posterior_samples = posterior_samples,
              clr.baselines = clr.baselines,
              clr.resampled = clr.resampled,
              Gamma = Gamma))
  
}

pull_predictions <- function(X, X_predict, clr.baselines, Lambda, Sigma, Gamma) {
  n_samples <- 1
  D <- dim(Sigma)[1]
  Theta.test <- function(X) {
    # X is now cbind(X, X_predict)
    Theta <- matrix(NA, length(clr.baselines[[1]]), ncol(X))
    for(j in 1:ncol(X)) {
      Theta[,j] <- clr.baselines[[X[2,j]]]
    }
    return(Theta)
  }
  nnew <- ncol(X_predict)
  # Set up Function Evaluation
  obs <- c(rep(TRUE, ncol(X)), rep(FALSE, nnew)) 
  Theta_realized <- Theta.test(cbind(X, X_predict))
  Gamma_realized <- Gamma(cbind(X, X_predict))
  # Predict Lambda
  Gamma_oo <- Gamma_realized[obs, obs, drop=F]
  Gamma_ou <- Gamma_realized[obs, !obs, drop=F]
  Gamma_uu <- Gamma_realized[!obs, !obs, drop=F]
  Gamma_ooIou <- solve(Gamma_oo, Gamma_ou)
  Gamma_schur <- Gamma_uu - t(Gamma_ou) %*% Gamma_ooIou 
  U_Gamma_schur <- chol(Gamma_schur)
  Theta_o <- Theta_realized[,obs, drop=F]
  Theta_u <- Theta_realized[,!obs, drop=F]
  Lambda_u <- array(0, dim=c(D, nnew, n_samples))
  # function for prediction - sample one Lambda_u
  lu <- function(Lambda_o, Sigma){
    Z <- matrix(rnorm(nrow(Gamma_uu)*D), D, nrow(Gamma_uu))
    Theta_u + (Lambda_o-Theta_o)%*%Gamma_ooIou + t(chol(Sigma))%*%Z%*%U_Gamma_schur
  }
  # Fill in and predict
  for(i in 1:n_samples) {
    Lambda_u[,,i] <- lu(Lambda[,,i], Sigma[,,i])
  }
  Eta <- array(0, dim=dim(Lambda_u))
  zEta <- array(rnorm((D)*nnew*n_samples), dim = dim(Eta))
  for(i in 1:n_samples) {
    Eta[,,i] <- Lambda_u[,,i] + t(chol(Sigma[,,i]))%*%zEta[,,i]
  }
  return(Eta)
}

# =====================================================================================================
#   VISUALIZATION
# =====================================================================================================

# input is a data.frame in the format
#
#     taxon sample resample     value
# 2       2      1        1 -4.180696
# 24      2      2        1 -4.293836
# 46      2      3        1 -5.359160
# 68      2      4        1 -4.601737
# 90      2      5        1 -7.509395
# 112     2      6        1 -5.492191
#
# 'sample' indexes the experimental sample
# 'resample' indexes the posterior sample
get_quantiles <- function(df) {
  out_df <- df %>%
    group_by(sample) %>%
    summarise(p2.5 = quantile(value, prob=0.025),
              p25 = quantile(value, prob=0.25),
              mean = mean(value),
              p75 = quantile(value, prob=0.75),
              p97.5 = quantile(value, prob=0.975)) %>%
    ungroup()
  return(out_df)  
}

get_ribbon_plot <- function(prediction_quantiles, data_points = NULL, title = NULL) {
  p <- ggplot() +
    geom_ribbon(data=prediction_quantiles, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(data=prediction_quantiles, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    geom_line(data=prediction_quantiles, aes(x=sample, y=mean), color="blue") +
    xlab("day") +
    ylab(paste0("CLR(abundance)")) +
    theme_minimal()
  if(!is.null(data_points)) {
    p <- p +
      geom_point(data=data_points, aes(x=sample, y=value))
  }
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

plot_predictive <- function(predictions, data, X_predict, plot_animal, plot_year, plot_taxon) {
  # get desired (dim 2) indices in the predictions for this animal and year
  idx_in_predictX <- which(X_predict[2,] == plot_animal & as.numeric(X_predict[3,]) == plot_year)
  
  # subset to desired indices and get quantiles
  subset_predictions <- predictions[,idx_in_predictX,]
  plot_df <- gather_array(subset_predictions, "value", "taxon", "sample", "resample")
  prediction_quantiles <- get_quantiles(plot_df[plot_df$taxon == plot_taxon,])

  # (laboriously) get something as near a grouth truth as we can manage
  # note: this will be more accurate for highly abundant taxa than lowly abundant ones
  sample_IDs <- get_label_type(data, "SampleID", year = plot_year)
  animals <- get_label_type(data, "Animal", year = plot_year)
  days <- get_label_type(data, "Day", year = plot_year)
  # "near truth" will be the empirical CLR values of these samples
  neartruth_days <- days[animals == plot_animal]
  # the 2nd argument to subselect_from_sampleIDs must be a list named by year
  neartruth_sample_IDs <- list(c(sample_IDs[animals == plot_animal]))
  names(neartruth_sample_IDs) <- paste0("year",plot_year)
  neartruth_counts <- subselect_from_sampleIDs(data, neartruth_sample_IDs)
  
  # CLR transform and select focal taxon
  D.taxa <- get_ntaxa(data)
  year_label <- paste0("year",plot_year)
  if(plot_taxon <= D.taxa$D.16S) {
    neartruth_subset <- neartruth_counts$`16S`[[year_label]]$filtered
  } else if(plot_taxon <= D.taxa$D.16S + D.taxa$D.ITS) {
    neartruth_subset <- neartruth_counts$ITS[[year_label]]$filtered
  } else {
    stop("Can't plot index larger than combined number of bacterial and fungal taxa!\n")
  }
  clr.counts <- t(driver::clr(t(neartruth_subset) + 0.5))
  clr.neartruth <- clr.counts[plot_taxon,] 
  neartruth_df <- data.frame(sample = neartruth_days, value = clr.neartruth)
  # plot
  get_ribbon_plot(prediction_quantiles, data_points = neartruth_df)
}

# =====================================================================================================
#   MAIN LOGIC WORKFLOW
# =====================================================================================================

# NOTES:
# (1) need to remove Nikos and baseline samples
# (2) match_sample_IDs_across_data_types should only be called if data_type = NULL

all_data <- get_data(data_type = NULL)
all_data <- match_sample_IDs_across_data_types(all_data)
all_data <- reorder_samples(all_data)

canonical_animal_order <- unique(get_label_type(all_data, "Animal"))

treatment_data <- pull_treatment_data(all_data, treatment = "CON")

design_matrices <- build_design_matrix(treatment_data, canonical_animal_order)
X <- design_matrices$X
X_predict <- design_matrices$X_predict

resampling_iterations <- 100

# right now fit_model and pull_predictions pull one sample each
all_predictions <- NULL
for(r_it in 1:resampling_iterations) {
  cat("Iteration:",r_it,"\n")
  # results
  #  |
  #  --- posterior_samples
  #  |    |
  #  |    --- Lambda
  #  |    |
  #  |    --- Sigma
  #  |
  #  --- clr.baselines
  #  |
  #  --- clr.resampled (taxa x treatment-specific samples)
  #  |
  #  --- Gamma
  results <- fit_model(treatment_data, X)
  # predictions is (taxa x interpolated samples x 1)
  # awkwardly, fit_model must have been run before this to make sure g_sigma
  # and rho_se are added to the GLOBAL workspace
  predictions <- pull_predictions(X, X_predict, results$clr.baselines,
                                  results$posterior_samples$Lambda,
                                  results$posterior_samples$Sigma,
                                  results$posterior_samples$Gamma)
  if(is.null(all_predictions)) {
    all_predictions <- predictions
  } else {
    all_predictions <- abind(all_predictions, predictions)
  }
}

plot_predictive(all_predictions, treatment_data, X_predict, plot_animal = "Onyx", plot_year = 1, plot_taxon = 3)

# still to do:
# (1) collecting high-confidence interactions (copy, paste)

  













