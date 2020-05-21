library(tidyverse)
  
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
      metadata.y <- read.csv(paste0("data/16S_",year_label,"_fecal_metadata.txt"), header=TRUE, sep='\t')
    } else {
      metadata.y <- read.csv(paste0("data/ITS_",year_label,"_fecal_metadata.txt"), header=TRUE, sep='\t')
    }
    
    metadata.y$SampleID <- as.character(metadata.y$SampleID)
    metadata.y <- metadata.y[!(metadata.y$SampleID %in% exclude_samples[[paste0("year",year)]]),]

    # standardize names
    metadata.y$Animal <- as.character(metadata.y$Animal)
    metadata.y$Animal <- sapply(metadata.y$Animal, function(x) sub(".*Nikos.*", "Nikos", x, perl=TRUE))
    metadata.y$Animal <- sapply(metadata.y$Animal, function(x) sub(".*Jones.*", "Jones", x, perl=TRUE))
    metadata.y$Animal <- as.factor(metadata.y$Animal)
    # Nikos received antibiotics, move him from "CON" to "ABX" treatment groups
    # this is indicated in Experimental_Groups.xlsx
    metadata.y[metadata.y$Animal == "Nikos",]$Treatment <- "ABX"
    # fix this awful space that prevents string matching
    new_levels <- levels(metadata.y$Animal)
    new_levels[which(new_levels == "Randy ")] <- "Randy"
    levels(metadata.y$Animal) <- new_levels
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

# so far treatment is only 
all_data <- get_data(data_type = "16S")

# plot changes in the mean log abundance over the samples
# how affected by antibiotic treatment is this in 16S? in ITS?
df <- NULL
for(treatment in c("CON", "ABX", "ABXFT")) {
  for(year in c(1, 2)) {
    use_samples <- all_data$`16S`[[paste0("year",year)]]$metadata$Treatment == treatment
    use_md <- all_data$`16S`[[paste0("year",year)]]$metadata[use_samples,c("Description","Day")]
    use_counts <- all_data$`16S`[[paste0("year",year)]]$filtered[,use_samples]
    log_counts <- log(use_counts + 0.5)
    df.temp <- data.frame(logmean = apply(log_counts, 2, mean), treatment = treatment, timepoint = use_md$Day)
    if(is.null(df)) {
      df <- df.temp
    } else {
      df <- rbind(df, df.temp)
    }
  }
}
ggplot(df) +
  geom_point(aes(x = timepoint, y = logmean, color = treatment))

# this is somewhat mitigated by the IQLR but you have to choose very tight quartiles!
year <- 1
treatment <- "ABX"
use_samples <- all_data$`16S`[[paste0("year",year)]]$metadata$Treatment == treatment
use_md <- all_data$`16S`[[paste0("year",year)]]$metadata[use_samples,c("Description","Day")]
use_counts <- all_data$`16S`[[paste0("year",year)]]$filtered[,use_samples]
clr.counts <- t(clr(t(use_counts) + 0.5)) # taxa x samples
clr.var <- apply(clr.counts, 1, var)
qs <- quantile(clr.var, probs = c(0.45, 0.55))
include_taxa <- clr.var > qs[1] & clr.var < qs[2]
log_counts <- log(use_counts[include_taxa,] + 0.5)
df <- data.frame(logmean = apply(log_counts, 2, mean), timepoint = use_md$Day)
ggplot(df) +
  geom_point(aes(x = timepoint, y = logmean, color = treatment))


str(all_data, max.level = 3)


# =====================================================================================================
#   NEXT -- TBD!!!
# =====================================================================================================

# go through each sample individually and resample counts
resample_counts <- function(all_data) {
  data_types <- names(all_data)
  for(data_type in data_types) {
    counts <- NULL
    if(is.null(counts)) {
      counts <- cbind(all_data[[data_type]]$year1$filtered, all_data[[data_type]]$year2$filtered)
    }
  }
  
  if(evaluate == "16S" | evaluate == "ITS") {
    Y <- as.matrix(cbind(counts.y1, counts.y2))
    eta <- t(driver::clr(t(Y) + 0.5)) # we'll use this to calculate overall CLR means
    # (for now)rm(list)
    N <- ncol(Y)
    resample_it <- 200
    eta.resampled <- array(NA, dim=c(D, N, resample_it))
    for(j in 1:resample_it) {
      # sample posterior of a multinomial-Dirichlet
      resampled_count_proportions <- Y
      for(samp in 1:ncol(Y)) {
        alpha.j <- unlist(resampled_count_proportions[,samp] + 0.5)
        resampled_count_proportions[,samp] <- t(LaplacesDemon::rdirichlet(1, alpha.j))
      }
      eta.resampled[,,j] <- t(driver::clr(t(resampled_count_proportions)))
    }
  } else if(evaluate == "both") {
    Y <- as.matrix(cbind(counts.y1, counts.y2))
    D <- D.16S + D.ITS
    eta <- rbind(t(driver::clr(t(Y[1:D.16S,]) + 0.5)),
                 t(driver::clr(t(Y[(D.16S+1):D,]) + 0.5)))
    N <- ncol(Y)
    resample_it <- 200
    eta.resampled <- array(NA, dim=c(D, N, resample_it))
    for(j in 1:resample_it) {
      # sample posterior of a multinomial-Dirichlet
      resampled_count_proportions <- Y
      for(samp in 1:ncol(Y)) {
        alpha.j.16S <- unlist(resampled_count_proportions[1:D.16S,samp] + 0.5)
        resampled_count_proportions[1:D.16S,samp] <- t(LaplacesDemon::rdirichlet(1, alpha.j.16S))
        alpha.j.ITS <- unlist(resampled_count_proportions[(D.16S+1):D,samp] + 0.5)
        resampled_count_proportions[(D.16S+1):D,samp] <- t(LaplacesDemon::rdirichlet(1, alpha.j.ITS))
      }
      eta.resampled[,,j] <- rbind(t(driver::clr(t(resampled_count_proportions[1:D.16S,]))),
                                  t(driver::clr(t(resampled_count_proportions[(D.16S+1):D,]))))
    }
  }
}



