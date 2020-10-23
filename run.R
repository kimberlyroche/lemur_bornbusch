rm(list = ls())

library(tidyverse)
library(driver)
library(stray)
library(Rcpp)
library(abind)
library(grid)
library(gridExtra)

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
  ggsave(file.path("output","change_in_log_mean.png"), p, dpi = 100, units = "in", height = 6, width = 10)
  
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
  ggsave(file.path("images","change_in_IQL_mean.png"), p, dpi = 100, units = "in", height = 6, width = 10)
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
    for(i in 1:length(OTU_list.all)) {
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
      summarise(agg_val = sum(value), .groups = 'keep') %>%
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
                    summarise(agg_val = sum(value), .groups = 'keep') %>%
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
get_data <- function(data_type = NULL, post_treatment_only = TRUE) {
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
                                  year2 = list(filtered = filtered_data$counts[,colnames(filtered_data$counts) %in% metadata$year2$SampleID],
                                           metadata = metadata$year2))
    # str(all_data, max.level = 3)
    
    # added later: remove baseline samples from treatment groups so as not to "dilute" signal
    tag_samples <- function(data, label) {
      include_samples <- !str_detect(data$metadata$Period, label)
      # chop "Pre" samples out of metadata
      data$metadata <- data$metadata[include_samples,]
      # used filtered metadata sample IDs/descriptions
      data$filtered <- data$filtered[,colnames(data$filtered) %in% c(tax_labels, as.character(data$metadata$Description))]
      return(data)
    }

    cat("Removing PRE-treatment samples...\n")
    all_data[[data_type]]$year1 <- tag_samples(all_data[[data_type]]$year1, label = "Pre")
    all_data[[data_type]]$year2 <- tag_samples(all_data[[data_type]]$year2, label = "Pre")

    if(post_treatment_only) {
      cat("Removing DURING-treatment samples...\n")
      all_data[[data_type]]$year1 <- tag_samples(all_data[[data_type]]$year1, label = "During")
      all_data[[data_type]]$year2 <- tag_samples(all_data[[data_type]]$year2, label = "During")
    }
    
    # added later: remove Nikos' samples
    tag_Nikos <- function(data) {
      include_samples <- !str_detect(data$metadata$Animal, "Nikos")
      data$metadata <- data$metadata[include_samples,]
      # used filtered metadata sample IDs/descriptions
      data$filtered <- data$filtered[,colnames(data$filtered) %in% c(tax_labels, as.character(data$metadata$Description))]
      return(data)
    }
    
    all_data[[data_type]]$year1 <- tag_Nikos(all_data[[data_type]]$year1)
    all_data[[data_type]]$year2 <- tag_Nikos(all_data[[data_type]]$year2)
    
  }
  if("16S" %in% data_types & "ITS" %in% data_types) {
    # take intersection of 16S and ITS samples
    all_data <- match_sample_IDs_across_data_types(all_data)
  }
  # make sure samples are blocked by individual (this make model fitting and visualization easier)
  all_data <- reorder_samples(all_data)
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
get_label_type <- function(data, label_type, year = NULL) {
  arb_data_type <- names(data)[[1]]
  if(!is.null(year)) {
    labels <- data[[arb_data_type]][[paste0("year", year)]]$metadata[[label_type]]
  } else {
    # grab both years
    labels <- c()
    for(year in c(1,2)) {
      year_label <- paste0("year", year)
      labels <- c(labels, data[[arb_data_type]][[year_label]]$metadata[[label_type]])
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
      
      reordering <- order(days_for_subset) # there are a handful of samples out of order; fix this
      samples_for_animal_year <- samples_for_animal_year[reordering]
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

# =====================================================================================================
#   MODEL FITTING
# =====================================================================================================

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
    if("year1" %in% names(data$`ITS`)) {
      D.ITS <- nrow(data$ITS[["year1"]]$filtered)
    } else if("year1" %in% names(data$`ITS`)) {
      D.ITS <- nrow(data$ITS[["year2"]]$filtered)
    }
  }
  return(list(D.16S = D.16S, D.ITS = D.ITS))
}

# multinomial-Dirichlet resample the count data set
# called by resample_dataset
# data set returned are CLR values with dimensions: {taxa (16S or ITS)} x {samples (Y1 + Y2)}
resample_data_type <- function(data, data_type, resample_iterations = 100) {
  counts <- cbind(data[[data_type]]$year1$filtered, data[[data_type]]$year2$filtered)
  clr.resampled <- array(NA, dim = c(nrow(counts), ncol(counts), resample_iterations))
  for(i in 1:resample_iterations) {
    for(j in 1:ncol(counts)) {
      augmented_counts <- unlist(counts[,j] + 0.5)
      clr.resampled[,j,i] <- c(LaplacesDemon::rdirichlet(1, augmented_counts))
    }
    clr.resampled[,,i] <- t(driver::clr(t(clr.resampled[,,i])))
  }
  return(clr.resampled)
}

resample_dataset <- function(data, resample_iterations = 100) {
  clr.resampled <- NULL
  if("16S" %in% names(data)) {
    clr.resampled <- resample_data_type(data, data_type = "16S", resample_iterations = resample_iterations)
  }
  if("ITS" %in% names(data)) {
    clr.resampled.ITS <- resample_data_type(data, data_type = "ITS", resample_iterations = resample_iterations)
    if(is.null(clr.resampled)) {
      clr.resampled <- clr.resampled.ITS
    } else {
      clr.resampled <- abind(clr.resampled, clr.resampled.ITS, along = 1)
    }
  }
  return(clr.resampled)
}

fit_GP <- function(data, X, iterations = 100) {
  # resample
  D.taxa <- get_ntaxa(data)
  D.16S <- D.taxa$D.16S
  D.ITS <- D.taxa$D.ITS
  clr.resampled <- resample_dataset(data, resample_iterations = iterations)
  
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
    # calculate baselines from one sample
    clr.baselines[[animal]] <- rowMeans(clr.resampled[,replace_idx,1,drop=F])
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
  D <- dim(Sigma)[1]
  n_samples <- dim(Sigma)[3]
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

fit_model <- function(data, animal_order, resample_iterations = 100, get_predictions = FALSE) {
  design_matrices <- build_design_matrix(data, animal_order)
  # fit
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
  fit <- fit_GP(data, design_matrices$X, iterations = resample_iterations)
  correlations <- fit$posterior_samples$Sigma
  for(i in 1:resample_iterations) {
    correlations[,,i] <- cov2cor(correlations[,,i])
  }
  if(get_predictions) {
    predictions <- pull_predictions(design_matrices$X,
                                   design_matrices$X_predict,
                                   fit$clr.baselines,
                                   fit$posterior_samples$Lambda,
                                   fit$posterior_samples$Sigma,
                                   fit$Gamma)
  }
  return(list(correlations = correlations, predictions = predictions,
              X = design_matrices$X, X_predict = design_matrices$X_predict))
}

get_taxon_label <- function(data, taxon_idx, data_type) {
  if(data_type == "16S") {
    tax_obj <- data$`16S`$tax
  } else if(data_type == "ITS") {
    tax_obj <- data$ITS$tax
  } else {
    D.taxa <- get_ntaxa(data)
    D.16S <- D.taxa$D.16S
    D.ITS <- D.taxa$D.ITS
    if(taxon_idx <= D.16S) {
      tax_obj <- data$`16S`$tax
    } else {
      tax_obj <- data$ITS$tax
      taxon_idx <- taxon_idx - D.16S
    }
  }
  
  label_level <- ncol(tax_obj)
  phylo_names <- colnames(tax_obj)
  label_content <- NA
  while(is.na(label_content) | str_detect(label_content, "uncultured") | str_detect(label_content, "gut metagenome") | str_detect(label_content, "unidentified")) {
    label_content <- tax_obj[taxon_idx,label_level]
    label <- paste(phylo_names[label_level], label_content)
    label_level <- label_level - 1
  }
  return(label)
}

get_high_conf_associations <- function(data, correlations, threshold = 0.5) {
  resample_iterations <- dim(correlations)[3]
  # NA out half the correlation matrix so we can dump it from the data.frame later and have
  # only unique associations
  # apply doesn't work here for whatever reason
  for(resample_iteration in 1:resample_iterations) {
    temp <- correlations[,,resample_iteration]
    temp[upper.tri(temp, diag=TRUE)] <- NA
    correlations[,,resample_iteration] <- temp
  }
  correlations_df <- driver::gather_array(correlations, "value", "taxon1", "taxon2", "sample")
  # remove redundant entries and diagonal (tagged as NA)
  correlations_df <- correlations_df[complete.cases(correlations_df),]
  # define 95% posterior inverval value
  cat("Calculating correlation quantiles for",resample_iterations,"iterations...\n")
  correlations_quantiles <- correlations_df %>%
    group_by(taxon1, taxon2) %>%
    summarise(p2.5 = quantile(value, prob=0.025),
              mean = mean(value),
              p97.5 = quantile(value, prob=0.975),
              .groups = 'keep') %>%
    ungroup()
  # "high confidence" indices have a 95% interval that doesn't span zero (i.e. they're consistent with respect to interaction sign)
  high_conf_idx <- which(sign(correlations_quantiles$p2.5) == sign(correlations_quantiles$p97.5))
  
  D.taxa <- get_ntaxa(data)
  D.16S <- D.taxa$D.16S
  D.ITS <- D.taxa$D.ITS
  
  high_conf_correlators <- correlations_quantiles[high_conf_idx,]
  high_conf_df <- data.frame(idx1 = c(), type1 = c(), label1 = c(),
                             idx2 = c(), type2 = c(), label2 = c(),
                             mean_correlation = c())
  for(i in 1:nrow(high_conf_correlators)) {
    if(abs(high_conf_correlators$mean[i]) >= threshold) {
      tax1_idx <- high_conf_correlators[i,]$taxon1
      tax2_idx <- high_conf_correlators[i,]$taxon2
      if(D.16S == 0) {
        tax1_type <- "ITS"
        tax2_type <- "ITS"
      } else if(D.ITS == 0) {
        tax1_type <- "16S"
        tax2_type <- "16S"
      } else {
        if(tax1_idx <= D.16S) {
          tax1_type <- "16S"
        } else {
          tax1_type <- "ITS"
          tax1_idx <- tax1_idx - D.16S
        }
        if(tax2_idx <= D.16S) {
          tax2_type <- "16S"
        } else {
          tax2_type <- "ITS"
          tax2_idx <- tax2_idx - D.16S
        }
      }
      tax1_label <- get_taxon_label(data, tax1_idx, data_type = tax1_type)
      tax2_label <- get_taxon_label(data, tax2_idx, data_type = tax2_type)
      high_conf_df <- rbind(high_conf_df, data.frame(idx1 = tax1_idx, type1 = tax1_type, label1 = tax1_label,
                                                     idx2 = tax2_idx, type2 = tax2_type, label2 = tax2_label,
                                                     mean_correlation = high_conf_correlators[i,]$mean))
    }    
  }
  return(high_conf_df)
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
              p97.5 = quantile(value, prob=0.975),
              .groups = 'keep') %>%
    ungroup()
  return(out_df)  
}

get_ribbon_plot <- function(prediction_quantiles, data_points = NULL, ylabel = NULL, title = NULL) {
  p <- ggplot() +
    geom_ribbon(data=prediction_quantiles, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(data=prediction_quantiles, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    geom_line(data=prediction_quantiles, aes(x=sample, y=mean), color="blue") +
    xlab("observation day")
  if(is.null(ylabel)) {
    p <- p + 
      ylab(paste0("CLR(abundance)"))
  } else {
    p <- p + 
      ylab(ylabel)
  }
  p <- p +
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

plot_predictive <- function(predictions, data, X_predict, plot_animal, plot_year, plot_taxon, ylabel = NULL) {
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
    adj_taxon_idx <- plot_taxon
  } else if(plot_taxon > D.taxa$D.16S & plot_taxon <= D.taxa$D.16S + D.taxa$D.ITS) {
    adj_taxon_idx <- plot_taxon - D.taxa$D.16S
    neartruth_subset <- neartruth_counts$ITS[[year_label]]$filtered
  } else {
    stop("Can't plot index larger than combined number of bacterial and fungal taxa!\n")
  }
  clr.counts <- t(driver::clr(t(neartruth_subset) + 0.5))
  clr.neartruth <- clr.counts[adj_taxon_idx,] 
  neartruth_df <- data.frame(sample = neartruth_days, value = clr.neartruth)
  
  # get desired (dim 2) indices in the predictions for this animal and year
  idx_in_predictX <- which(X_predict[2,] == plot_animal & as.numeric(X_predict[3,]) == plot_year)
  if(is.null(idx_in_predictX)) {
    stop("Animal has no samples in this year!\n")
  }
  
  # subset to desired indices and get quantiles
  subset_predictions <- predictions[,idx_in_predictX,]
  plot_df <- gather_array(subset_predictions, "value", "taxon", "sample", "resample")
  plot_df$sample <- plot_df$sample + min(neartruth_days) - 1 # e.g. start at 3 instead of 1 if the first
                                                             #      sample day is day 3
  
  prediction_quantiles <- get_quantiles(plot_df[plot_df$taxon == plot_taxon,])
  
  # pass to plotting fn.
  return(get_ribbon_plot(prediction_quantiles, data_points = neartruth_df, ylabel = ylabel))
}

# plot random animal x year x taxon
diagnostic_plot <- function(data, model_output, data_type = NULL) {
  # diagnostic; fails for "Nikos" with one sample in Y2
  total_taxa <- get_ntaxa(data)
  D.16S <- total_taxa$D.16S
  D.ITS <- total_taxa$D.ITS
  D.all <- D.16S + D.ITS
  y_idx <- sample(1:2)[1]
  ua <- unique(get_label_type(data, "Animal", year = y_idx))
  animal <- "Nikos"
  while(animal == "Nikos") {
    a_idx <- sample(1:length(ua))[1]
    animal <- ua[a_idx]
  }
  if(is.null(data_type)) {
    t_idx <- sample(1:D.all)[1]
  } else if(data_type == "16S") {
    t_idx <- sample(1:D.16S)[1]
  } else if(data_type == "ITS") {
    t_idx <- sample(1:D.ITS)[1] + D.16S
  } else {
    stop("Not an allowable data type!")
  }
  cat("Plotting taxon",t_idx,"in animal",animal,"in year",y_idx,"\n")
  p <- plot_predictive(model_output$predictions, data, model_output$X_predict,
                  plot_animal = animal,
                  plot_year = y_idx,
                  plot_taxon = t_idx)
  show(p)
}

plot_network <- function(data_type = "16S", conditions = c("CON", "ABX", "ABXFT"), which_sign = NULL, show_plot = FALSE) {

  # ============================================================================================================
  #   parse data
  # ============================================================================================================
  
  filenames <- list("during-post" = file.path("output","network_input_during-post.txt"),
                    "post" = file.path("output","network_input_post.txt"))
  setups <- names(filenames)
  correlators <- NULL
  for(setup in setups) {
    input_file <- filenames[[setup]]
    setup_correlators <- tryCatch({
      read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
      # empty file
      data.frame(feature_1 = NA, feature_2 = NA, interaction_strength = NA)
    })
    setup_correlators$setup <- setup
    if(is.null(correlators)) {
      correlators <- setup_correlators
    } else {
      correlators <- rbind(correlators, setup_correlators)
    }
  }
  
  correlators <- correlators[correlators$feature_1_type == data_type & correlators$feature_2_type == data_type,]
  correlators <- correlators[correlators$condition %in% conditions,]
  
  if(nrow(correlators) == 0) {
    return(NULL)
  }
  
  View(correlators)
  
  # For testing
  correlators <- correlators[c(2:4, 70:72),]

  # grab the unique features in this table, e.g.
  # feature_1    feature_2    interaction_strength
  #         A            B    ...
  #         A            C    ...
  #         B            C    ...
  #
  # gives c(A, B, C)
  items <- sort(unique(c(correlators$feature_1, correlators$feature_2)))
  n_items <- length(items)

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
  topmost_idx <- which(nodes$y == max(nodes$y))[1]
  if(topmost_idx > 1) {
    nodes$label <- c(seq(from = (n_items-topmost_idx+2), to = n_items), seq(from = 1,
                                                                            to = (n_items - topmost_idx + 1)))
  }

  for(condition in conditions) {
    
    # ============================================================================================================
    #   now that we have positions, label the interaction edges
    # ============================================================================================================
    for(setup in setups) { 
      n_interactions_cond <- nrow(correlators[correlators$condition == condition & correlators$setup == setup,])
      correlators_cond <- correlators[correlators$condition == condition & correlators$setup == setup,]
      edges <- data.frame(x1 = numeric(n_interactions_cond), y1 = numeric(n_interactions_cond),
                          x2 = numeric(n_interactions_cond), y2 = numeric(n_interactions_cond),
                          weight = numeric(n_interactions_cond), sign = character(n_interactions_cond))
      for(i in 1:n_interactions_cond) {
        # string matching could be a problem here; in which case we'd just want unique numeric indices for
        # "features"
        edges$x1[i] <- nodes[nodes$item == correlators_cond[i,]$feature_1,]$x
        edges$y1[i] <- nodes[nodes$item == correlators_cond[i,]$feature_1,]$y
        edges$x2[i] <- nodes[nodes$item == correlators_cond[i,]$feature_2,]$x
        edges$y2[i] <- nodes[nodes$item == correlators_cond[i,]$feature_2,]$y
        edges$weight[i] <- abs(correlators_cond[i,]$interaction_strength)*11 # fat line weight
        #edges$weight[i] <- abs(correlators_cond[i,]$interaction_strength)*5 # moderate line weight
        # edge thickness in the plot
        # edges$sign[i] <- sign(correlators_cond[i,]$interaction_strength)
        if(correlators_cond[i,]$interaction_strength >= 0) {
          edges$sign[i] <- "positive"
        } else {
          edges$sign[i] <- "negative"
        }
      }
      dim(edges)
      
      # ============================================================================================================
      #   PLOT WITH LABELS
      # ============================================================================================================
    
      colors <- c("#E6194B", # red
                  "#3CB44B") # green
      # if(!is.null(which_sign)) {
      #   if(which_sign == "negative") {
      #     edges <- edges[edges$sign < 0,]
      #     colors <- colors[1]
      #   } else {
      #     edges <- edges[edges$sign > 0,]
      #     colors <- colors[2]
      #   }
      # }
      
      # assign factor to signs
      edges$sign <- factor(edges$sign, levels = c("negative", "positive"))
      
      p <- ggplot() +
        # geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2, color = sign, size = weight)) +
        geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2, color = as.factor(sign), size = weight), alpha = 0.66) +
        # scale_color_manual(values = colors) +
        scale_color_manual(values = c("positive" = "#3CB44B", "negative" = "#E6194B")) +
        scale_size_identity() + # use the width specified by `weight`
        geom_point(data = nodes, aes(x = x, y = y), size = 9) +
        geom_text(data = nodes, aes(x = x, y = y, label = label), size = 4, color = "#FFFFFF") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        theme(legend.position = "none")
      
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
    
      margin <- 1.05
      # xlimits <- c(min(c(edges$x1, edges$x2))*margin, max(c(edges$x1, edges$x2, legend_df$x_item*2))*margin)
      # ylimits <- c(min(c(edges$y1, edges$y2, legend_df$y_item))*margin, max(c(edges$y1, edges$y2, legend_df$y_item))*margin)
      xlimits <- c(min(nodes$x)*margin, max(c(nodes$x, legend_df$x_item*2))*margin)
      ylimits <- c(min(c(nodes$x, legend_df$y_item))*margin, max(c(nodes$x, legend_df$y_item))*margin)

      # manually plot legend
      p <- p +
        geom_point(data = legend_df, aes(x = x_label, y = y_label), size = 7) + 
        geom_text(data = legend_df, aes(x = x_label, y = y_label, label = label), size = 4, color = "#FFFFFF") +
        geom_text(data = legend_df, aes(x = x_item, y = y_item, label = item), size = 4, hjust = 0) +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        xlim(xlimits) +
        ylim(ylimits)
      if(show_plot) {
        show(p)
      }
      sign_append <- ""
      if(!is.null(which_sign)) {
        sign_append <- paste0("_",which_sign)
      }
      for(extension in c(".eps",".png")) {
        ggsave(file.path("output",paste0(data_type,"_",condition,"_interactions",sign_append,"_",setup,extension)),
               p, units="in", dpi=150, height = (abs(ylimits[1]) + abs(ylimits[2]))/4, width = (abs(xlimits[1]) + abs(xlimits[2]))/4)
      }
    }
  }
}

# =====================================================================================================
#   MAIN LOGIC WORKFLOW
# =====================================================================================================

# improvements:
# (1) switch CLR to IQLR?
# (2) could calculate baselines from overall resampling in fit_GP()
#
# for Brianna, Zack -- how can we generalize this such that sub-item YEARS don't need to be present?!
#
# basically no correlators show up in the ITS data unless you relax the correlation threshold to 0.4
# and these are not convincing

threshold <- 0.5 # if threshold is zero, we'll evaluate thresholding: visualizing correlators within range
                 # of apparent posterior correlation (e.g. [0.4, 0.5))

use_extra_viz <- FALSE # extra visualization
                      # (1) correlation matrices per condition
                      # (2) plot trajectories ("squiggle plots") for a few selected pairs of taxa

plot_all_taxa <- FALSE

for(post_treatment_only in c(TRUE, FALSE)) {
  for(data_type in c("16S")) {
    if(data_type == "BOTH") {
      all_data <- get_data(data_type = NULL, post_treatment_only = post_treatment_only)
    } else {
      all_data <- get_data(data_type = data_type, post_treatment_only = post_treatment_only)
    }
    
    conditions <- c("CON", "ABX", "ABXFT")
    for(treatment in conditions) {
      cat("Evaluating",data_type,"x",treatment,"\n")
      treatment_data <- pull_treatment_data(all_data, treatment = treatment)
      # set the order of individuals in the combined series
      canonical_animal_order <- unique(get_label_type(treatment_data, "Animal"))
      
      model_output <- fit_model(treatment_data, canonical_animal_order, resample_iterations = 500, get_predictions = TRUE)
      # sanity check the fit via: diagnostic_plot(treatment_data, model_output, data_type = NULL) # plot a random feature
      
      # saveRDS(apply(model_output$correlations, c(1,2), mean), file = paste0(treatment, "_Sigma.rds"))
    
      if(data_type == "16S" & use_extra_viz) {
        # plot correlation matrices
        # we'll look at the average Sigma
        mean_Sigma_corr <- apply(model_output$correlations, c(1,2), mean)
        mean_Sigma_long <- driver::gather_array(mean_Sigma_corr, "value", "row", "col")
        p <- ggplot(mean_Sigma_long, aes(row, col, fill=value)) + 
          geom_tile() +
          scale_fill_gradient2(low = "darkblue", high = "darkred", name = "correlation")
        for(extension in c(".png")) {
          ggsave(file.path("output",paste0(data_type,"_",treatment,"_correlmat",extension)),
                 p, dpi = 100, units = "in", height = 8, width = 10)
        }
        # plot selected trajectories
        if(treatment == "CON") {
          idx1 <- 157
          idx2 <- 127
        }
        if(treatment == "ABX") {
          idx1 <- 7
          idx2 <- 26
        }
        if(treatment == "ABXFT") {
          idx1 <- 110
          idx2 <- 26
        }
        label1 <- get_taxon_label(data = treatment_data, taxon_idx = idx1, data_type = data_type)
        label2 <- get_taxon_label(data = treatment_data, taxon_idx = idx2, data_type = data_type)
        plots <- list()
        plot_pair <- sample(1:ncol(model_output$X_predict))[1]
        plot_animal <- unname(model_output$X_predict[2,plot_pair])
        plot_year <- as.numeric(unname(model_output$X_predict[3,plot_pair]))
        p1 <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
                              plot_animal = plot_animal, plot_year = plot_year, plot_taxon = idx1)
        p1 <- p1 +
          ylab(label1) +
          theme(plot.title = element_text(size = 12))
        plots[[length(plots)+1]] <- p1
        p2 <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
                              plot_animal = plot_animal, plot_year = plot_year, plot_taxon = idx2)
        p2 <- p2 +
          ylab(label2) +
          theme(plot.title = element_text(size = 12))
        plots[[length(plots)+1]] <- p2
        p <- grid.arrange(grobs = plots, ncol = 1, top = textGrob(paste0(plot_animal," (Y",plot_year,")")))
        ggsave(file.path("output",paste0(data_type,"_",treatment,"_strongcorrelator.png")),
               p, dpi = 100, units = "in", height = 6, width = 8)
      }
      
      if(threshold > 0) {
        # evaluate high confidence correlators at the specified threshold mean (high confidence) correlation
        hcc <- get_high_conf_associations(treatment_data, model_output$correlations, threshold = threshold)
        cat("Found",nrow(hcc),"high confidence interactors\n")
        if(nrow(hcc) > 0) {
          write_tsv(x = hcc, path = file.path("output",paste0(data_type,"_",treatment,"_highconftable.tsv")))
          
          # # combinations of individual x year
          # indiv_year_combos <- data.frame(individual = model_output$X[2,], year = model_output$X[3,])
          # indiv_year_combos <- as.data.frame(indiv_year_combos %>%
          #   group_by(individual, year) %>%
          #   tally())
          # 
          # plots <- list()
          # #for(i in 1:nrow(indiv_year_combos)) {
          # for(i in 1:2) {
          #   cat(i,"\n")
          #   p1 <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
          #                         plot_animal = indiv_year_combos[i,]$individual, plot_year = indiv_year_combos[i,]$year,
          #                         plot_taxon = hcc[2,]$idx1, ylabel = paste0(hcc[2,]$label1, "\n", indiv_year_combos[i,]$individual))
          #   p2 <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
          #                         plot_animal = indiv_year_combos[i,]$individual, plot_year = indiv_year_combos[i,]$year,
          #                         plot_taxon = hcc[2,]$idx2, ylabel = paste0(hcc[2,]$label2, "\n", indiv_year_combos[i,]$individual))
          #   plots[[length(plots) + 1]] <- p1
          #   plots[[length(plots) + 1]] <- p2
          # }
          # do.call("grid.arrange", c(plots, ncol = 2))
  
          if(post_treatment_only) {
            network_file <- file.path("output","network_input_post.txt")
          } else {
            network_file <- file.path("output","network_input_during-post.txt")
          }
          if(!file.exists(network_file)) {
            write("feature_1\tfeature_2\tinteraction_strength\tfeature_1_type\tfeature_2_type\tcondition", file = network_file)
          }
          for(i in 1:nrow(hcc)) {
            write(paste0(hcc$label1[i],"\t",hcc$label2[i],"\t",hcc$mean_correlation[i],"\t",
                         hcc$type1[i],"\t",hcc$type2[i],"\t",treatment),
                  file = network_file,
                  append = TRUE)
          }
        }
      } else if(data_type != "BOTH") {
        # This is no longer in use but I'm loath to throw it away.
        correlations <- model_output$correlations
        resample_iterations <- dim(correlations)[3]
        for(resample_iteration in 1:resample_iterations) {
          temp <- correlations[,,resample_iteration]
          temp[upper.tri(temp, diag=TRUE)] <- NA
          correlations[,,resample_iteration] <- temp
        }
        correlations_df <- driver::gather_array(correlations, "value", "taxon1", "taxon2", "sample")
        # remove redundant entries and diagonal (tagged as NA)
        correlations_df <- correlations_df[complete.cases(correlations_df),]
        # this takes ~30 sec. to run
        median_post_corr <- as.data.frame(correlations_df %>%
          group_by(taxon1, taxon2) %>%
          summarise(median = median(value),
                    lower95 = quantile(value, probs = c(0.025)),
                    upper95 = quantile(value, probs = c(0.975)),
                    .groups = 'keep'))
        bounds <- seq(0.3, 0.9, by = 0.1)
        for(bound in bounds) {
          cat("Evaluating pairs with correlation",bound,"to",(bound + 0.1),"\n")
          criterion1 <- median_post_corr$median > bound
          criterion2 <- median_post_corr$median < (bound + 0.1)
          qualified_pairs <- which(criterion1 & criterion2)
          cat("\tFound",length(qualified_pairs),"pairs...\n")
          if(length(qualified_pairs) > 0) {
            idx <- qualified_pairs[1:min(3,length(qualified_pairs))]
            for(i in 1:length(idx)) {
              r_idx <- idx[i]
              idx1 <- median_post_corr$taxon1[r_idx]
              idx2 <- median_post_corr$taxon2[r_idx]
              # get labels for these taxa
              label1 <- get_taxon_label(data = treatment_data, taxon_idx = idx1, data_type = data_type)
              label2 <- get_taxon_label(data = treatment_data, taxon_idx = idx2, data_type = data_type)
  
              plots <- list()
              pairs_used <- c()
              for(j in 1:min(4, length(unique(apply(model_output$X_predict[2:3,], 2, function(x) paste0(x[1], x[2])))))) {
                # plot k instances of this pair
                new_pair <- FALSE
                while(!new_pair) {
                  plot_pair <- sample(1:ncol(model_output$X_predict))[1]
                  plot_animal <- unname(model_output$X_predict[2,plot_pair])
                  plot_year <- as.numeric(unname(model_output$X_predict[3,plot_pair]))
                  slug <- paste0(plot_animal,plot_year)
                  if(!(slug %in% pairs_used)) {
                    new_pair <- TRUE
                  }
                }
                pairs_used <- c(pairs_used, slug)
  
                p1 <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
                                      plot_animal = plot_animal, plot_year = plot_year, plot_taxon = idx1)
                p1 <- p1 +
                  ggtitle(paste0(plot_animal," (Y",plot_year,")")) +
                  theme(plot.title = element_text(size = 12))
                plots[[length(plots)+1]] <- p1
                p2 <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
                                      plot_animal = plot_animal, plot_year = plot_year, plot_taxon = idx2)
                p2 <- p2 +
                  ggtitle(paste0(plot_animal," (Y",plot_year,")")) +
                  theme(plot.title = element_text(size = 12))
                plots[[length(plots)+1]] <- p2
              }
  
              p <- grid.arrange(grobs = plots, ncol = 2, top = textGrob(paste0(label1," x ",label2,"\n",
                                                                               "median correlation: ",round(median_post_corr$median[r_idx], 2),"\n",
                                                                               "95% interval: (",round(median_post_corr$lower95[r_idx], 2),", ",
                                                                                                round(median_post_corr$upper95[r_idx], 2),")")))
              ggsave(file.path("output","diagnostic",paste0(data_type,"_",treatment,"_interval_",bound,"_",i,".png")),
                     p, dpi = 100, units = "in", height = 10, width = 10)
            }
          }
        }
      }
  
      if(plot_all_taxa) {
        # plot logratio change in all taxa, so we can see which were killed off by the antibiotics
    
        # combinations of individual x year
        indiv_year_combos <- data.frame(individual = model_output$X[2,], year = model_output$X[3,])
        indiv_year_combos <- as.data.frame(indiv_year_combos %>%
          group_by(individual, year) %>%
          tally())
        
        limit <- 10
        plot_height <- 12
        if(treatment %in% c("ABX", "ABXFT")) {
          limit <- nrow(treatment_data$`16S`$year1$filtered)
          plot_height <- 8
        }
        for(tax_idx in 1:limit) {
          cat("Plotting taxon",tax_idx,"\n")
          plots <- list()
          tax_label <- get_taxon_label(treatment_data, taxon_idx = tax_idx, data_type = data_type)
          for(i in 1:nrow(indiv_year_combos)) {
            p <- plot_predictive(model_output$predictions, treatment_data, model_output$X_predict,
                                  plot_animal = indiv_year_combos[i,]$individual, plot_year = indiv_year_combos[i,]$year,
                                  plot_taxon = tax_idx, ylabel = paste0(indiv_year_combos[i,]$individual, " (Year ",indiv_year_combos[i,]$year,")"))
            plots[[length(plots) + 1]] <- p
          }
          p <- grid.arrange(grobs = plots, top = tax_label, ncol = 2)
          ggsave(file.path("output", treatment, paste0("taxon_",tax_idx,".png")), plot = p, units = "in",
                 dpi = 150, height = plot_height, width = 6)
        }
        
        # taxa are really different! validate by looking at direct data?
        # md <- treatment_data$`16S`$year2$metadata
        # counts <- treatment_data$`16S`$year2$filtered
        # clr.counts <- as.matrix(counts)
        # clr.counts <- t(clr(t(clr.counts) + 0.5))
        # selected_samples <- md[md$Animal == "Onyx" & md$Year == 2,]$Description
        # selected_days <- md[md$Animal == "Onyx" & md$Year == 2,]$Day
        # plot(selected_days, clr.counts[1,colnames(counts) %in% selected_samples])
      }
    }
  }
}

if(data_type != "BOTH") {
  plot_network(data_type = "16S", conditions = c("CON", "ABX", "ABXFT"), which_sign = NULL, show_plot = FALSE)
}





