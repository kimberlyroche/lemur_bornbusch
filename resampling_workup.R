rm(list=ls())

# these can probably be pared down
library(driver)
library(ggplot2)
#library(vegan)
library(tidyverse)
#library(RColorBrewer)
#library(stringr)
library(stray) # for SE function, overkill
library(driver)
#library(ggraph)
#library(igraph)
library(Rcpp)
library(LaplacesDemon)
library(gridExtra)

setwd("C:/Users/kim/Documents/lemur_bornbusch/")

# pre-treatment:
# (1) Renamed 16S_*_fecal_metdata.txt to 16S_*_fecal_metadata.txt
# (2) removed pound-sign from headers in *_Y1_fecal_metadata.txt, *_Y2_fecal_metadata.txt,
#     Bornbusch_*_Y1_fecal_count-table.tsv, and Bornbusch_*_Y2_fecal_count-table.tsv
# (3) Added header labels TAX and CONF to last two columns of Bornbusch_*_Y1_fecal_count-table.tsv
#     and Bornbusch_*_Y2_fecal_count-table.tsv
# (4) Removed extra empty lines from end of ITS_*_fecal_metadata.txt
# (5) In ITS_*_fecal_metadata.txt changed "Sample ID" header to "SampleID"

# ============================================================================================================
#   HIGH-LEVEL VARIABLES
# ============================================================================================================

# agglomeration levels
agg_16S_level <- "genus"
agg_ITS_level <- "genus"

# define covariates and hyperparameters
treatment <- "ABX" # "CON", "ABX", "ABXFT"

# do we need to subset to matched samples for 16S and ITS?
evaluate <- "both" # "both", "16S", "ITS"

# ============================================================================================================
#   FUNCTIONS
# ============================================================================================================

read_sequence_counts <- function(which_data="16S") {
  if(which_data == "16S") {
    sequences.y1 <- read.csv(paste0("data/Bornbusch_16S_Y1_fecal_count-table.tsv"), header=TRUE, sep='\t')
  } else {
    sequences.y1 <- read.csv(paste0("data/Bornbusch_ITS_Y1_fecal_count-table.tsv"), header=TRUE, sep='\t')
  }
  rownames(sequences.y1) <- sequences.y1$OTU.ID
  tax.y1 <- as.list(sequences.y1$TAX)
  names(tax.y1) <- rownames(sequences.y1)
  sequences.y1 <- sequences.y1[,!(colnames(sequences.y1) %in% c("OTU.ID", "TAX", "CONF", exclude_samples$year1))]
  
  if(which_data == "16S") {
    sequences.y2 <- read.csv(paste0("data/Bornbusch_16S_Y2_fecal_count-table.tsv"), header=TRUE, sep='\t')
  } else {
    sequences.y2 <- read.csv(paste0("data/Bornbusch_ITS_Y2_fecal_count-table.tsv"), header=TRUE, sep='\t')
  }
  rownames(sequences.y2) <- sequences.y2$OTU.ID
  tax.y2 <- as.list(sequences.y2$TAX)
  names(tax.y2) <- rownames(sequences.y2)
  sequences.y2 <- sequences.y2[,!(colnames(sequences.y2) %in% c("OTU.ID", "TAX", "CONF", exclude_samples$year2))]
  
  if(which_data == "16S" & file.exists("data/tax_collapsed_y1_16S.rds") & file.exists("data/tax_collapsed_y1_16S.rds")) {
    tax.all <- readRDS("data/tax_all_16S.rds")
    tax.collapsed.y1 <- readRDS("data/tax_collapsed_y1_16S.rds")
    tax.collapsed.y2 <- readRDS("data/tax_collapsed_y2_16S.rds")
  } else if(which_data == "ITS" & file.exists("data/tax_collapsed_y1_ITS.rds") & file.exists("data/tax_collapsed_y1_ITS.rds")) {
    tax.all <- readRDS("data/tax_all_ITS.rds")
    tax.collapsed.y1 <- readRDS("data/tax_collapsed_y1_ITS.rds")
    tax.collapsed.y2 <- readRDS("data/tax_collapsed_y2_ITS.rds")
  } else {
    # coerce all OTU id's to match by 
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
      if(which_data == "16S") {
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
    if(which_data == "16S") {
      tax.collapsed <- tax.collapsed[which(tax.collapsed$order != "Chloroplast"),]
    }
    
    # separate the years back out
    tax.collapsed.y1 <- as.data.frame(tax.collapsed)
    tax.collapsed.y1 <- tax.collapsed.y1[,colnames(tax.collapsed.y1) %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species", sample_list.y1)]
    tax.collapsed.y2 <- as.data.frame(tax.collapsed)
    tax.collapsed.y2 <- tax.collapsed.y2[,colnames(tax.collapsed.y2) %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species", sample_list.y2)]
    
    append_str <- which_data
    saveRDS(tax.collapsed, paste0("data/tax_all_",append_str,".rds"))
    saveRDS(tax.collapsed.y1, paste0("data/tax_collapsed_y1_",append_str,".rds"))
    saveRDS(tax.collapsed.y2, paste0("data/tax_collapsed_y2_",append_str,".rds"))
  }
  
  # at the genus level this is about 74% zeros
  return(list(tax=tax.all, raw_year1=tax.collapsed.y1, raw_year2=tax.collapsed.y2))
}

read_metadata <- function(year=1, which_data="16S") {
  if(year == 1) {
    year_label <- "Y1"
  }
  if(year == 2) {
    year_label <- "Y2"
  }
  
  if(which_data == "16S") {
    metadata <- read.csv(paste0("data/16S_",year_label,"_fecal_metadata.txt"), header=TRUE, sep='\t')
  } else {
    metadata <- read.csv(paste0("data/ITS_",year_label,"_fecal_metadata.txt"), header=TRUE, sep='\t')
  }
  
  metadata$SampleID <- as.character(metadata$SampleID)
  
  if(year == 1) {
    #metadata <- metadata[!(metadata$SampleID %in% c("LCAXC1","LCAXC2")),]
    metadata <- metadata[!(metadata$SampleID %in% exclude_samples$year1),]
  }
  if(year == 2) {
    #metadata <- metadata[!(metadata$SampleID %in% c("CON1", "CON3", "LCAX038",
    #                                                "LCAX050", "LCAX071", "LCAX073", "LCAX114", "LCAX130",
    #                                                "LCAX137", "LCAX141", "LCAX142", "LCAX149", "LCAX152",
    #                                                "LCAX157", "LCAX158", "LCAX173")),]
    metadata <- metadata[!(metadata$SampleID %in% exclude_samples$year2),]
  }
  
  # standardize names
  metadata$Animal <- as.character(metadata$Animal)
  metadata$Animal <- sapply(metadata$Animal, function(x) sub(".*Nikos.*", "Nikos", x, perl=TRUE))
  metadata$Animal <- sapply(metadata$Animal, function(x) sub(".*Jones.*", "Jones", x, perl=TRUE))
  metadata$Animal <- as.factor(metadata$Animal)
  # Nikos received antibiotics, move him from "CON" to "ABX" treatment groups
  # this is indicated in Experimental_Groups.xlsx
  metadata[metadata$Animal == "Nikos",]$Treatment <- "ABX"
  # fix this awful space that prevents string matching
  new_levels <- levels(metadata$Animal)
  new_levels[which(new_levels == "Randy ")] <- "Randy"
  levels(metadata$Animal) <- new_levels
  return(metadata)
}

# parse out the taxonomic labels from either 16S or ITS count table
get_counts <- function(df) {
  max_tax_id <- max(which(sapply(df[1,], is.numeric) == FALSE))
  return(df[,(max_tax_id+1):ncol(df)])
}

# input is counts of a single taxon over n samples
filter_low_taxa_sub <- function(x, min_abundance=10, min_representation=0.2) {
  sum(x >= min_abundance)/length(x) >= min_representation
}

# filter out taxa below a threshold count in some percentage of samples
filter_low_taxa <- function(data.y1, data.y2, verbose=FALSE) {
  counts.y1 <- get_counts(data.y1)
  counts.y2 <- get_counts(data.y2)
  counts.all <- cbind(counts.y1, counts.y2)
  retain_idx <- apply(counts.all, 1, filter_low_taxa_sub)
  filtered.y1 <- data.y1[retain_idx,]
  filtered.y2 <- data.y2[retain_idx,]
  filtered.y1 <- rbind(filtered.y1,
                       c(list(kingdom=NA,phylum=NA,class=NA,order=NA,family=NA,genus=NA,species=NA),
                         apply(counts.y1[!retain_idx,], 2, sum)))
  filtered.y2 <- rbind(filtered.y2,
                       c(list(kingdom=NA,phylum=NA,class=NA,order=NA,family=NA,genus=NA,species=NA),
                         apply(counts.y2[!retain_idx,], 2, sum)))
  if(verbose) {
    counts.y1 <- get_counts(filtered.y1)
    counts.y2 <- get_counts(filtered.y2)
    cat("Year 1:\n")
    cat("\tPercent zeros:",round(sum(counts.y1 == 0)/(nrow(counts.y1)*ncol(counts.y1))*100),"% zeros\n")
    cat("\tTaxa x samples:",nrow(counts),"x",ncol(counts),"\n")
  }
  return(list(year1=filtered.y1, year2=filtered.y2))
}

# filter out taxa < 1 total counts
filter_low_taxa_alt <- function(df, percent_of_whole=1) {
  (apply(df, 1, sum)/sum(df)) > (percent_of_whole/100)
}

agglomerate_data <- function(data.y1, data.y2, level="genus", filter_percent=NULL) {
  # ignore variable names -- these can be sequence counts of bacterial or fungal sequence variants
  # note: things undefined at the specified taxonomic level will be filtered into "Other"
  tax_pieces <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_pieces <- tax_pieces[1:which(tax_pieces == level)]
  tax_pieces <- c(tax_pieces, "experiment")
  min_sample_idx <- min(which(sapply(data.y1$raw_data[1,], is.numeric)))
  max_sample_idx <- max(which(sapply(data.y1$raw_data[1,], is.numeric)))
  # first exclude anything not resolved to this level
  # we'll incorporate these counts into Other later
  unresolved <- is.na(data.y1$raw_data[,c(level)]) | is.na(data.y2$raw_data[,c(level)])
  collapsed_unresolved.y1 <- colSums(data.y1$raw_data[unresolved,min_sample_idx:max_sample_idx])
  resolved.y1 <- data.y1
  resolved.y1$raw_data <- resolved.y1$raw_data[!unresolved,]
  resolved.y1$metadata <- resolved.y1$metadata[!unresolved,]
  bacteria.agg.y1 <- as.data.frame(resolved.y1$raw_data %>%
                                     gather(experiment, value, min_sample_idx:max_sample_idx) %>%
                                     group_by(!!!syms(tax_pieces)) %>%
                                     summarise(agg_val = sum(value)) %>%
                                     spread(experiment, agg_val))
  min_sample_idx <- min(which(sapply(data.y2$raw_data[1,], is.numeric)))
  max_sample_idx <- max(which(sapply(data.y2$raw_data[1,], is.numeric)))
  collapsed_unresolved.y2 <- colSums(data.y2$raw_data[unresolved,min_sample_idx:max_sample_idx])
  resolved.y2 <- data.y2
  resolved.y2$raw_data <- resolved.y2$raw_data[!unresolved,]
  resolved.y2$metadata <- resolved.y2$metadata[!unresolved,]
  bacteria.agg.y2 <- as.data.frame(resolved.y2$raw_data %>%
                                     gather(experiment, value, min_sample_idx:max_sample_idx) %>%
                                     group_by(!!!syms(tax_pieces)) %>%
                                     summarise(agg_val = sum(value)) %>%
                                     spread(experiment, agg_val))
  bacteria.agg <- cbind(bacteria.agg.y1,
                        bacteria.agg.y2[,(max(which(sapply(bacteria.agg.y2[1,], is.numeric) == FALSE)) + 1):ncol(bacteria.agg.y2)])
  counts <- get_counts(bacteria.agg)
  level_idx <- length(tax_pieces)-1 # genus=6, species=7, etc.
  # filter on percent abundance if necessary
  # roll unresolved counts into last row ("Other")
  if(is.null(filter_percent)) {
    retain_idx <- 1:nrow(counts)
    retained_tax <- bacteria.agg[retain_idx,1:level_idx]
    bacteria.chopped <- counts[retain_idx,]
    bacteria.chopped <- rbind(bacteria.chopped, c(collapsed_unresolved.y1, collapsed_unresolved.y2))
    retained_tax <- rbind(retained_tax, rep("Other", level_idx))
  } else {
    retain_idx <- filter_low_taxa_alt(counts, percent_of_whole=filter_percent)
    retained_tax <- bacteria.agg[retain_idx,1:level_idx]
    # collapse "Other"; if there are D taxa, retained tax will be D-1 items long
    #   the last item (or row in the counts table) is the assumed "Other" group
    bacteria.chopped <- counts[retain_idx,]
    bacteria.chopped[nrow(bacteria.chopped)+1,] <- apply(counts[!retain_idx,], 2, sum) +
      c(collapsed_unresolved.y1, collapsed_unresolved.y2)
    retained_tax[nrow(retained_tax)+1,] <- "Other"
  }
  list(counts=bacteria.chopped, level=level_idx, tax=retained_tax)
  
  return(list(counts=bacteria.chopped, level=level_idx, tax=retained_tax))
}

# ============================================================================================================
#   DATA PARSING
# ============================================================================================================

# we'll build this data structure (list of lists)

# all_data
#  |
#  --- 16S
#  |    |
#  |    --- year1
#  |    |    |
#  |    |    --- filtered
#  |    |    |
#  |    |    --- raw_data
#  |    |    |
#  |    |    --- metadata
#  |    |
#  |    --- year2
#  |         |
#  |         --- filtered
#  |         |
#  |         --- raw_data
#  |         |
#  |         --- metadata
#  |
#  --- ITS
#       |
#       --- year1
#       |    |
#       |    --- filtered
#       |    |
#       |    --- raw_data
#       |    |
#       |    --- metadata
#       |
#       --- year2
#            |
#            --- filtered
#            |
#            --- raw_data
#            |
#            --- metadata
#
#   filtered -- taxa (rows) x {phylogeny, samples} (columns) as kingdom ... species LCAX001 LCAX002 ... LCAX104
#   raw_data -- as above but uncollapsed
#   metadata -- samples (rows) x {SampleID, BarcodeSequence, LinkerPrimerSequence, BarcodePlate, Well,
#                                  Animal, Year, Date, Period, Session, Day, Treatment}

all_data <- list()

if(evaluate == "16S" | evaluate == "both") {
  
  exclude_samples <- list(year1=c("LCAXC1","LCAXC2"),
                          year2=c("CON1", "CON3", "LCAX038",
                                  "LCAX050", "LCAX071", "LCAX073", "LCAX114", "LCAX130",
                                  "LCAX137", "LCAX141", "LCAX142", "LCAX149", "LCAX152",
                                  "LCAX157", "LCAX158", "LCAX173"))
  
  # returned as a tibble with columns: kingdom phylum class order family genus species LCAX001 LCAX002 LCAX003 LCAX004 ...
  year.all <- read_sequence_counts(which_data="16S")
  tax <- year.all$tax
  year1 <- list(raw_data=year.all$raw_year1)
  year2 <- list(raw_data=year.all$raw_year2)
  year1$metadata <- read_metadata(year=1, which_data="16S")
  year2$metadata <- read_metadata(year=2, which_data="16S")
  
  # agglomerate and filter out by 1% percent abundance
  filtered.all <- agglomerate_data(year1, year2, level=agg_16S_level, filter_percent=1)
  year1$filtered <- cbind(filtered.all$tax, filtered.all$counts[,colnames(filtered.all$counts) %in% year1$metadata$SampleID])
  year2$filtered <- cbind(filtered.all$tax, filtered.all$counts[,colnames(filtered.all$counts) %in% year2$metadata$SampleID])
  
  all_data[["16S"]] <- list(year1=year1, year2=year2)

  # print pre-/post-filtering zero percent
  # zero.y1.before <- sum(get_counts(year1$raw_data) == 0)/(nrow(get_counts(year1$raw_data))*ncol(get_counts(year1$raw_data)))
  # zero.y2.before <- sum(get_counts(year2$raw_data) == 0)/(nrow(get_counts(year2$raw_data))*ncol(get_counts(year2$raw_data)))
  # zero.y1.after <- sum(get_counts(year1$filtered) == 0)/(nrow(get_counts(year1$filtered))*ncol(get_counts(year1$filtered)))
  # zero.y2.after <- sum(get_counts(year2$filtered) == 0)/(nrow(get_counts(year2$filtered))*ncol(get_counts(year2$filtered)))
  # cat("Year 1 zeros:",round(zero.y1.before*100, 2),"% ->",round(zero.y1.after*100, 2),"%\n")
  # cat("Year 2 zeros:",round(zero.y2.before*100, 2),"% ->",round(zero.y2.after*100, 2),"%\n")
}

if(evaluate == "ITS" | evaluate == "both") {
  # get ITS data
  
  exclude_samples <- list(year1=c("CON1", "CON2", "CON5", "LCAXC1"),
                          year2=c("CON3", "LCAX038",
                                  "LCAX050", "LCAX071", "LCAX073", "LCAX114", "LCAX130",
                                  "LCAX137", "LCAX141", "LCAX142", "LCAX149", "LCAX152",
                                  "LCAX157", "LCAX158", "LCAX173"))
  
  # returned as a tibble with columns: kingdom phylum class order family genus species LCAX001 LCAX002 LCAX003 LCAX004 ...
  year.all <- read_sequence_counts(which_data="ITS")
  tax <- year.all$tax
  year1 <- list(raw_data=year.all$raw_year1)
  year2 <- list(raw_data=year.all$raw_year2)
  year1$metadata <- read_metadata(year=1, which_data="ITS")
  year2$metadata <- read_metadata(year=2, which_data="ITS")
  
  # filter at a few levels and look at sparsity in ITS data
  # NOTE: really these should be CLR abundances etc. but I doubt the distributions would look much different
  if(FALSE) {
    # combine individuals from a given treatment across years and assess the normality/lack of this distribution
    for(local_level in c("family", "order", "class", "phylum")) {
      # agglomerate and filter out by 1% percent abundance
      filtered.all <- agglomerate_data(year1, year2, level=local_level, filter_percent=1)
      year1$filtered <- cbind(filtered.all$tax,
                              filtered.all$counts[,colnames(filtered.all$counts) %in% year1$metadata$SampleID])
      year2$filtered <- cbind(filtered.all$tax,
                              filtered.all$counts[,colnames(filtered.all$counts) %in% year2$metadata$SampleID])
      # temporarily add to the larger data structure
      all_data[["ITS"]] <- list(year1=year1, year2=year2)
      
      for(treatment in c("CON", "ABX", "ABXFT")) {
        cat("Plotting histogram of",treatment,"at level",local_level,"...\n")
        sampleIDs.y1 <- all_data[["ITS"]]$year1$metadata[year1$metadata$Treatment == treatment,]$SampleID
        sampleIDs.y2 <- all_data[["ITS"]]$year2$metadata[year2$metadata$Treatment == treatment,]$SampleID
        counts.y1 <- get_counts(all_data[["ITS"]]$year1$filtered)[,sampleIDs.y1]
        counts.y2 <- get_counts(all_data[["ITS"]]$year2$filtered)[,sampleIDs.y2]
        counts.all_years <- as.matrix(cbind(counts.y1, counts.y2))
        tax_labels <- filtered.all$tax[[local_level]]
        df <- data.frame(x=c(), tax_label=c())
        for(i in sample(1:nrow(counts.all_years))) {
          df <- rbind(df, data.frame(x=log(counts.all_years[i,] + 0.5), tax_label=paste0(local_level, " ", tax_labels[i])))
        }
        ncol <- 8
        p <- ggplot(df, aes(x)) +
          geom_histogram(bins=10) +
          facet_wrap(. ~ tax_label, ncol=ncol) +
          xlab("log(abundance)") +
          theme(strip.text.x = element_text(size = 6))
        img_width <- 10
        img_height <- 10
        if(local_level == "phylum") {
          img_height <- 3
        } else if(local_level == "class") {
          img_height <- 3
        } else if(local_level == "order") {
          img_height <- 4
        } else if(local_level == "family") {
          img_height <- 5
        }
        ggsave(paste0("images/ITS_",local_level,"_",treatment,"_logcount_distros.png"), p, units="in", dpi=150, height=img_height, width=img_width)
      }
    }
  }
  
  # moving on with final ITS filtering
  filtered.all <- agglomerate_data(year1, year2, level=agg_ITS_level, filter_percent=1)
  year1$filtered <- cbind(filtered.all$tax,
                          filtered.all$counts[,colnames(filtered.all$counts) %in% year1$metadata$SampleID])
  year2$filtered <- cbind(filtered.all$tax,
                          filtered.all$counts[,colnames(filtered.all$counts) %in% year2$metadata$SampleID])
  
  all_data[["ITS"]] <- list(year1=year1, year2=year2)
}

str(all_data, max.level=3)
# List of 2
# $ 16S:List of 2
# ..$ year1:List of 3
# .. ..$ raw_data:'data.frame':	953 obs. of  165 variables:   (all taxa x kingdom::species + sample ID)
#   .. ..$ metadata:'data.frame':	158 obs. of  13 variables:
#   .. ..$ filtered:'data.frame':	21 obs. of  164 variables:  (filtered taxa x kingdom::genus + sample ID)
#   ..$ year2:List of 3
# .. ..$ raw_data:'data.frame':	953 obs. of  238 variables:
#   .. ..$ metadata:'data.frame':	231 obs. of  15 variables:
#   .. ..$ filtered:'data.frame':	21 obs. of  237 variables:
#   $ ITS:List of 2
# ..$ year1:List of 3
# .. ..$ raw_data:'data.frame':	1452 obs. of  116 variables:
#   .. ..$ metadata:'data.frame':	109 obs. of  15 variables:
#   .. ..$ filtered:'data.frame':	17 obs. of  115 variables:
#   ..$ year2:List of 3
# .. ..$ raw_data:'data.frame':	1452 obs. of  168 variables:
#   .. ..$ metadata:'data.frame':	161 obs. of  15 variables:
#   .. ..$ filtered:'data.frame':	17 obs. of  167 variables:

# print pre-/post-filtering zero percent; this data set is much sparser!
zero.y1.before <- sum(get_counts(year1$raw_data) == 0)/(nrow(get_counts(year1$raw_data))*ncol(get_counts(year1$raw_data)))
zero.y2.before <- sum(get_counts(year2$raw_data) == 0)/(nrow(get_counts(year2$raw_data))*ncol(get_counts(year2$raw_data)))
zero.y1.after <- sum(get_counts(year1$filtered) == 0)/(nrow(get_counts(year1$filtered))*ncol(get_counts(year1$filtered)))
zero.y2.after <- sum(get_counts(year2$filtered) == 0)/(nrow(get_counts(year2$filtered))*ncol(get_counts(year2$filtered)))
cat("Year 1 zeros:",round(zero.y1.before*100, 2),"% ->",round(zero.y1.after*100, 2),"%\n")
cat("Year 2 zeros:",round(zero.y2.before*100, 2),"% ->",round(zero.y2.after*100, 2),"%\n")

if(FALSE & evaluate == "both") {
  # how ab-normal are these? pretty bad...
  df <- data.frame(x=c(log(as.matrix(get_counts(all_data$`16S`$year1$filtered))+0.5)), type="16S year1")
  df <- rbind(df, data.frame(x=c(log(as.matrix(get_counts(all_data$`16S`$year2$filtered))+0.5)), type="16S year2"))
  df <- rbind(df, data.frame(x=c(log(as.matrix(get_counts(all_data$`ITS`$year1$filtered))+0.5)), type="ITS year1"))
  df <- rbind(df, data.frame(x=c(log(as.matrix(get_counts(all_data$`ITS`$year2$filtered))+0.5)), type="ITS year2"))
  
  p <- ggplot(df, aes(x)) +
    geom_density() +
    facet_wrap(. ~ type, nrow=1)
  ggsave("images/log_densities.png", p, units="in", dpi=150, height=3, width=12)
}

# delete these to keep ourselves honeset and prevent reusing the wrong year's content
rm(year.all)
rm(tax)
rm(year1)
rm(year2)
rm(filtered.all)
rm(exclude_samples)

# ============================================================================================================
#   PRELIMINARY LABELS AND STUFF
# ============================================================================================================

# get the sample IDs we want and animal names assoc. with them
# these track with sample indices for each year
if(evaluate == "16S" | evaluate == "ITS") {
  sampleIDs.y1 <- all_data[[evaluate]]$year1$metadata[all_data[[evaluate]]$year1$metadata$Treatment == treatment,]$SampleID
  sampleIDs.y2 <- all_data[[evaluate]]$year2$metadata[all_data[[evaluate]]$year2$metadata$Treatment == treatment,]$SampleID
  animals.y1 <- all_data[[evaluate]]$year1$metadata[all_data[[evaluate]]$year1$metadata$Treatment == treatment,]$Animal
  animals.y2 <- all_data[[evaluate]]$year2$metadata[all_data[[evaluate]]$year2$metadata$Treatment == treatment,]$Animal
  days.y1 <- all_data[[evaluate]]$year1$metadata[all_data[[evaluate]]$year1$metadata$Treatment == treatment,]$Day
  days.y2 <- all_data[[evaluate]]$year2$metadata[all_data[[evaluate]]$year2$metadata$Treatment == treatment,]$Day
  counts.y1 <- get_counts(all_data[[evaluate]]$year1$filtered)[,sampleIDs.y1]
  counts.y2 <- get_counts(all_data[[evaluate]]$year2$filtered)[,sampleIDs.y2]
  D <- nrow(counts.y1)
} else if(evaluate == "both") {
  # for now intersect ITS, 16S samples but we'll have to re-evaluate this later
  sampleIDs.y1.16S <- all_data$`16S`$year1$metadata[all_data$`16S`$year1$metadata$Treatment == treatment,]$SampleID
  sampleIDs.y2.16S <- all_data$`16S`$year2$metadata[all_data$`16S`$year2$metadata$Treatment == treatment,]$SampleID
  sampleIDs.y1.ITS <- all_data$`ITS`$year1$metadata[all_data$`ITS`$year1$metadata$Treatment == treatment,]$SampleID
  sampleIDs.y2.ITS <- all_data$`ITS`$year2$metadata[all_data$`ITS`$year2$metadata$Treatment == treatment,]$SampleID

  # filter these to respective 16S/ITS sample ID sets
  subsetted.md.16S.y1 <- all_data$`16S`$year1$metadata[all_data$`16S`$year1$metadata$Treatment == treatment & all_data$`16S`$year1$metadata$SampleID %in% sampleIDs.y1.16S,]
  subsetted.md.16S.y2 <- all_data$`16S`$year2$metadata[all_data$`16S`$year2$metadata$Treatment == treatment & all_data$`16S`$year2$metadata$SampleID %in% sampleIDs.y2.16S,]
  subsetted.md.ITS.y1 <- all_data$`ITS`$year1$metadata[all_data$`ITS`$year1$metadata$Treatment == treatment & all_data$`ITS`$year1$metadata$SampleID %in% sampleIDs.y1.ITS,]
  subsetted.md.ITS.y2 <- all_data$`ITS`$year2$metadata[all_data$`ITS`$year2$metadata$Treatment == treatment & all_data$`ITS`$year2$metadata$SampleID %in% sampleIDs.y2.ITS,]
  
  # there's pretty substantial missingness between 16S and ITS samples, let's visualize this
  days.y1.16S <- subsetted.md.16S.y1$Day
  days.y2.16S <- subsetted.md.16S.y2$Day
  days.y1.ITS <- subsetted.md.ITS.y1$Day
  days.y2.ITS <- subsetted.md.ITS.y2$Day

  animals.y1.16S <- subsetted.md.16S.y1$Animal
  animals.y2.16S <- subsetted.md.16S.y2$Animal
  animals.y1.ITS <- subsetted.md.ITS.y1$Animal
  animals.y2.ITS <- subsetted.md.ITS.y2$Animal
  
  # define factor labels for labels "Teres 16S" (etc.) ahead of time so we can control their order in levels()
  # this will make for better plotting
  ht_levels <- sort(unique(c(paste(animals.y1.16S, "16S"),
                      paste(animals.y2.16S, "16S"),
                      paste(animals.y1.ITS, "ITS"),
                      paste(animals.y2.ITS, "ITS"))), decreasing=TRUE)
  df <- data.frame(day=days.y1.16S, hosttype=sapply(paste(animals.y1.16S, "16S"), function(x) factor(x, ht_levels)), assay="16S", sampleID=sampleIDs.y1.16S)
  df <- rbind(df, data.frame(day=days.y2.16S, hosttype=sapply(paste(animals.y2.16S, "16S"), function(x) factor(x, ht_levels)), assay="16S", sampleID=sampleIDs.y2.16S))
  df <- rbind(df, data.frame(day=days.y1.ITS, hosttype=sapply(paste(animals.y1.ITS, "ITS"), function(x) factor(x, ht_levels)), assay="ITS", sampleID=sampleIDs.y1.ITS))
  df <- rbind(df, data.frame(day=days.y2.ITS, hosttype=sapply(paste(animals.y2.ITS, "ITS"), function(x) factor(x, ht_levels)), assay="ITS", sampleID=sampleIDs.y2.ITS))
  df$alpha <- -1
  
  # crude but for now models matched 16S/ITS samples only; have to think about how to incorporate updates otherwise
  sampleIDs.y1 <- intersect(sampleIDs.y1.16S, sampleIDs.y1.ITS)
  sampleIDs.y2 <- intersect(sampleIDs.y2.16S, sampleIDs.y2.ITS)
  
  df[df$sampleID %in% c(sampleIDs.y1, sampleIDs.y2),]$alpha <- 1.0
  
  p <- ggplot(df, aes(x=day, y=hosttype, color=assay, alpha=alpha)) +
    geom_point(size=3) +
    scale_alpha(range = c(0.33, 1.0), guide="none") # hide legend
  ggsave(paste0("images/16S_ITS_sample_pairing_",treatment,".png"), plot=p, units="in", dpi=150, height=8, width=15)
  
  # re-filter to these matched sampleIDs only (same for both 16S and ITS)
  subsetted.md.y1 <- all_data$`16S`$year1$metadata[all_data$`16S`$year1$metadata$Treatment == treatment & all_data$`16S`$year1$metadata$SampleID %in% sampleIDs.y1,]
  subsetted.md.y2 <- all_data$`16S`$year2$metadata[all_data$`16S`$year2$metadata$Treatment == treatment & all_data$`16S`$year2$metadata$SampleID %in% sampleIDs.y2,]
  animals.y1 <- subsetted.md.y1$Animal
  animals.y2 <- subsetted.md.y2$Animal
  days.y1 <- subsetted.md.y1$Day
  days.y2 <- subsetted.md.y2$Day
  D.16S <- nrow(all_data$`16S`$year1$filtered) # same between years
  D.ITS <- nrow(all_data$`ITS`$year1$filtered)
  counts.y1 <- rbind(get_counts(all_data$`16S`$year1$filtered)[,sampleIDs.y1],
                     get_counts(all_data$`ITS`$year1$filtered)[,sampleIDs.y1])
  counts.y2 <- rbind(get_counts(all_data$`16S`$year2$filtered)[,sampleIDs.y2],
                     get_counts(all_data$`ITS`$year2$filtered)[,sampleIDs.y2])
}

# re-order all for clarity (1) by animal
new_order.y1 <- order(animals.y1)
new_order.y2 <- order(animals.y2)
sampleIDs.y1 <- sampleIDs.y1[new_order.y1]
sampleIDs.y2 <- sampleIDs.y2[new_order.y2]
animals.y1 <- animals.y1[new_order.y1]
animals.y2 <- animals.y2[new_order.y2]
days.y1 <- days.y1[new_order.y1]
days.y2 <- days.y2[new_order.y2]

nsamp.y1 <- length(sampleIDs.y1)
nsamp.y2 <- length(sampleIDs.y2)

# ============================================================================================================
#   COUNT RESAMPLING
# ============================================================================================================

if(evaluate == "16S" | evaluate == "ITS") {
  Y <- as.matrix(cbind(counts.y1, counts.y2))
  eta <- t(driver::clr(t(Y) + 0.5)) # we'll use this to calculate overall CLR means
                                    # (for now)
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

# ============================================================================================================
#   SET UP COVARIATES, HYPERPARAMETERS FOR A GIVEN CONDITION
# ============================================================================================================

for(animal in unique(animals.y1)) {
  subset_idx <- which(animals.y1 == animal)
  new_order <- order(days.y1[subset_idx])
  sampleIDs.y1[subset_idx] <- sampleIDs.y1[subset_idx][new_order]
  animals.y1[subset_idx] <- animals.y1[subset_idx][new_order]
  days.y1[subset_idx] <- days.y1[subset_idx][new_order]
}
for(animal in unique(animals.y2)) {
  subset_idx <- which(animals.y2 == animal)
  new_order <- order(days.y2[subset_idx])
  sampleIDs.y2[subset_idx] <- sampleIDs.y2[subset_idx][new_order]
  animals.y2[subset_idx] <- animals.y2[subset_idx][new_order]
  days.y2[subset_idx] <- days.y2[subset_idx][new_order]
}

# get animal in year1 and year2
unique.animals.y1 <- unique(as.character(animals.y1))
unique.animals.y2 <- unique(as.character(animals.y2))

# set up kernel over samples
# parameter settings are utterly fucking heuristic rn
dc <- 0.01 # desired minimum correlation
dd_se <- 7
rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
g_sigma <- 1.5

Gamma <- function(X) {
  SE(X[1,,drop=F], sigma=g_sigma, rho=rho_se, jitter=1e-10) # time only
}

X <- matrix(NA, 2, N)
X_predict <- NULL
new_sample_order <- c()
bump <- 1000 # hack; days to separate each individual's series by
# such that they're effectively independent
offset <- 0
idx.animals.y1 <- list()
days.animals.y1 <- list()
clr_means.y1 <- list()
for(a in 1:length(unique.animals.y1)) {
  animal <- as.character(unique.animals.y1)[a]
  bump.animal <- (a-1)*bump
  idx.animal.y1 <- which(animals.y1 == animal)
  # NOTE: THIS SHOULD PROBABLY BE MODIFIED TO USE THE CLR MEAN OF EACH RESAMPLE (TBD)
  if(evaluate == "16S" | evaluate == "ITS") {
    clr_means.y1[[animal]] <- rowMeans(eta[,idx.animal.y1])
  } else {
    clr_means.y1[[animal]] <- c(rowMeans(eta[1:D.16S,idx.animal.y1]),
                                    rowMeans(eta[(D.16S+1):D,idx.animal.y1]))
  }
  days.animal.y1 <- days.y1[idx.animal.y1]
  days.animal.y1 <- days.animal.y1 - min(days.animal.y1) # 0 at baseline
  days.animal.y1 <- days.animal.y1 + bump.animal
  days.animals.y1[[animal]] <- days.animal.y1
  X[1,(offset + 1):(offset + length(days.animal.y1))] <- days.animal.y1
  X[2,(offset + 1):(offset + length(days.animal.y1))] <- a
  idx.animals.y1[[animal]] <- (offset + 1):(offset + length(days.animal.y1))
  offset <- offset + length(days.animal.y1)
  X_predict.sub <- matrix(0, 2, max(days.animal.y1)-min(days.animal.y1)+1)
  X_predict.sub[1,] <- min(days.animal.y1):max(days.animal.y1)
  X_predict.sub[2,] <- a
  if(is.null(X_predict)) {
    X_predict <- X_predict.sub
  } else {
    X_predict <- cbind(X_predict, X_predict.sub)
  }
}
nsamp.predict.y1 <- ncol(X_predict)
idx.animals.y2 <- list()
days.animals.y2 <- list()
clr_means.y2 <- list()
for(aa in 1:length(unique.animals.y2)) {
  animal <- as.character(unique.animals.y2)[aa]
  bump.animal <- ((a + aa)-1)*bump
  idx.animal.y2 <- which(animals.y2 == animal)
  if(evaluate == "16S" | evaluate == "ITS") {
    clr_means.y2[[animal]] <- rowMeans(eta[,(idx.animal.y2+nsamp.y1)])
  } else {
    clr_means.y2[[animal]] <- c(rowMeans(eta[1:D.16S,(idx.animal.y2+nsamp.y1)]),
                                rowMeans(eta[(D.16S+1):D,(idx.animal.y2+nsamp.y1)]))
  }
  days.animal.y2 <- days.y2[idx.animal.y2]
  days.animal.y2 <- days.animal.y2 - min(days.animal.y2) # 0 at baseline
  days.animal.y2 <- days.animal.y2 + bump.animal
  days.animals.y2[[animal]] <- days.animal.y2
  X[1,(offset + 1):(offset + length(days.animal.y2))] <- days.animal.y2
  X[2,(offset + 1):(offset + length(days.animal.y2))] <- a + aa
  idx.animals.y2[[animal]] <- (offset + 1):(offset + length(days.animal.y2))
  offset <- offset + length(days.animal.y2)
  X_predict.sub <- matrix(0, 2, max(days.animal.y2)-min(days.animal.y2)+1)
  X_predict.sub[1,] <- min(days.animal.y2):max(days.animal.y2)
  X_predict.sub[2,] <- (a + aa)
  if(is.null(X_predict)) {
    X_predict <- X_predict.sub
  } else {
    X_predict <- cbind(X_predict, X_predict.sub)
  }
}

# per-animal named list of indices within Y, eta
# idx.animals.y1
# idx.animals.y2

# per-animal named list of day offsets within Y, eta (separated by 1000-day bumps)
# days.animals.y1
# days.animals.y2

# per-animal named list of D-vectors of average CLR values
# clr_means.y1
# clr_means.y2

upsilon <- D-1+10 # lesser certainty
Xi <- diag(D)

ids.y1 <- X[2,1:nsamp.y1]
ids.y2 <- X[2,(nsamp.y1+1):N]
ids.y2 <- ids.y2 - max(ids.y1)

Theta.assignments <- c(paste0(unique.animals.y1[ids.y1],"Y1"), paste0(unique.animals.y2[ids.y2],"Y2"))
clr_means <- list()
for(animal in unique.animals.y1[ids.y1]) {
  clr_means[[paste0(animal,"Y1")]] <- clr_means.y1[[animal]]
}
for(animal in unique.animals.y2[ids.y2]) {
  clr_means[[paste0(animal,"Y2")]] <- clr_means.y2[[animal]]
}

Theta <- function(X) {
  # need these things as global:
  #  Theta.assignments: per-sample label e.g. "AracusY1"
  #  clr_means: named list (labeled as "AracusY1" etc.) of mean clr
  temp <- matrix(NA, D, ncol(X))
  for(ta in 1:length(Theta.assignments)) {
    a <- Theta.assignments[ta]
    temp[,ta] <- clr_means[[a]]
  }
  return(temp)
}

# ============================================================================================================
#   GET POSTERIOR ESTIMATES FOR LAMBDA, SIGMA FROM BAYESIAN MULTIVARIATE LINEAR REGRESSION FORMS
# ============================================================================================================

# these calculations are lifted directly from their forms in stray::pibble
# implementing them in R (especially the Inv-Wishart sampling is just too unstable)

# uncollapsing/posterior sampling of Lambda, Sigma
sourceCpp("sampling.cpp")

# GP prediction
# is the CLR going to cause problems here?
GPpredict <- function(X, X_predict, Lambda, Sigma, D, iter) {
  
  # what should the format of X_predict be?
  nnew <- ncol(X_predict)
  
  # Set up Function Evaluation
  obs <- c(rep(TRUE, ncol(X)), rep(FALSE, nnew)) 
  Theta_realized <- Theta(cbind(X, X_predict))
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
  
  #Lambda_u <- array(0, dim=c(D-1, nnew, iter))
  Lambda_u <- array(0, dim=c(D, nnew, iter))
  
  # function for prediction - sample one Lambda_u
  lu <- function(Lambda_o, Sigma){
    Z <- matrix(rnorm(nrow(Gamma_uu)*D), D, nrow(Gamma_uu))
    Theta_u + (Lambda_o-Theta_o)%*%Gamma_ooIou + t(chol(Sigma))%*%Z%*%U_Gamma_schur
  }
  
  # Fill in and predict
  for (i in 1:iter){
    Lambda_u[,,i] <- lu(Lambda[,,i], Sigma[,,i])
  }
  
  Eta <- array(0, dim=dim(Lambda_u))
  zEta <- array(rnorm((D)*nnew*iter), dim = dim(Eta))
  for (i in 1:iter){
    Eta[,,i] <- Lambda_u[,,i] + t(chol(Sigma[,,i]))%*%zEta[,,i]
  }
  return(Eta)
}

X_realized <- diag(ncol(X))
posterior <- uncollapse_simple(eta.resampled,
                         X_realized, # per fit_basset.R: 53
                         Theta(X),
                         Gamma(X),
                         Xi,
                         upsilon)

# png(paste0("C:/Users/kim/Desktop/model2_",treatment,".png"))
# image(apply(posterior$Sigma, c(1,2), mean))
# dev.off()

ids.y1 <- X_predict[2,1:nsamp.predict.y1]
ids.y2 <- X_predict[2,(nsamp.predict.y1+1):ncol(X_predict)]
ids.y2 <- ids.y2 - max(ids.y1)
Theta.assignments.append <- c(paste0(unique.animals.y1[ids.y1],"Y1"), paste0(unique.animals.y2[ids.y2],"Y2"))
Theta.assignments <- c(Theta.assignments, Theta.assignments.append)

predictions <- GPpredict(X, X_predict, posterior$Lambda, posterior$Sigma, D, resample_it)

# ============================================================================================================
#   QUICK VISUALIZATION, SANITY CHECK THE ESTIMATES
# ============================================================================================================

# look at the average Sigma
mean_Sigma <- apply(posterior$Sigma, c(1,2), mean)
mean_Sigma_corr <- cov2cor(mean_Sigma)
mean_Sigma_long <- driver::gather_array(mean_Sigma_corr, "value", "row", "col")

p <- ggplot(mean_Sigma_long, aes(row, col, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low="darkblue", high="darkred", name="correlation")
ggsave(paste0("images/paired_correlation_",evaluate,"_",treatment,".png"), plot=p, units="in", dpi=150, height=6, width=7.5)

# diagnostic stuff; I think this is just looking at random fits
if(FALSE) {
  if(evaluate == "both") {
    lr_coords <- c(sample(1:D.16S)[1:2], sample(1:D.ITS)[1:2])
    for(it in 1:4) {
      # choose a logratio_coordinate
      lr_coord <- lr_coords[it]
      
      # choose an individual ID
      indiv.y1 <- sample(unique(ids.y1))[1]
      indiv.y2 <- sample(unique(ids.y2))[1]
      
      if(it < 3) {
        # 16S
        label_level <- max(which(sapply(all_data$`16S`$year1$filtered[lr_coord,], is.character)))
        label <- paste(names(all_data$`16S`$year1$filtered[label_level]), all_data$`16S`$year1$filtered[lr_coord,label_level])
        cat("Examining 16S log ratio coordinate",lr_coord,"(",label,") in individuals",unique.animals.y1[[indiv.y1]],"and",unique.animals.y2[[indiv.y2]],"...\n")
      } else {
        # ITS
        label_level <- max(which(sapply(all_data$`ITS`$year1$filtered[lr_coord,], is.character)))
        label <- paste(names(all_data$`ITS`$year1$filtered[label_level]), all_data$`ITS`$year1$filtered[lr_coord,label_level])
        cat("Examining ITS log ratio coordinate",lr_coord,"(",label,") in individuals",unique.animals.y1[[indiv.y1]],"and",unique.animals.y2[[indiv.y2]],"...\n")
      }
      
      # get truth and predictions for YEAR 1
      ground.eta.y1 <- gather_array(eta[,idx.animals.y1[[indiv.y1]]], value, "taxon", "sample")
      ground.eta.y1 <- ground.eta.y1[ground.eta.y1$taxon == lr_coord,]
      ground.eta.y1$sample <- plyr::mapvalues(ground.eta.y1$sample, 
                                              from=ground.eta.y1$sample, 
                                              to=days.y1[idx.animals.y1[[indiv.y1]]])
      
      # get truth and predictions for YEAR 2
      ground.eta.y2 <- gather_array(eta[,idx.animals.y2[[indiv.y2]]], value, "taxon", "sample")
      ground.eta.y2 <- ground.eta.y2[ground.eta.y2$taxon == lr_coord,]
      ground.eta.y2$sample <- plyr::mapvalues(ground.eta.y2$sample, 
                                              from=ground.eta.y2$sample, 
                                              to=days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]) # map 1 2 3 ... to 4001 4002 4003 ...
    
      # get predictions for YEAR 1
      predicted.eta.y1 <- gather_array(predictions[,which(ids.y1 == indiv.y1),], "value", "taxon", "sample", "resample")
      predicted.eta.y1 <- predicted.eta.y1[predicted.eta.y1$taxon == lr_coord,]
      predicted.eta.y1$sample <- plyr::mapvalues(predicted.eta.y1$sample, 
                                                 from=min(predicted.eta.y1$sample):max(predicted.eta.y1$sample), 
                                                 to=min(days.y1[idx.animals.y1[[indiv.y1]]]):max(days.y1[idx.animals.y1[[indiv.y1]]]))
      
      # get predictions for YEAR 2
      predicted.eta.y2 <- gather_array(predictions[,length(ids.y1) + which(ids.y2 == indiv.y2),], "value", "taxon", "sample", "resample")
      predicted.eta.y2 <- predicted.eta.y2[predicted.eta.y2$taxon == lr_coord,]
      predicted.eta.y2$sample <- plyr::mapvalues(predicted.eta.y2$sample, 
                                                 from=min(predicted.eta.y2$sample):max(predicted.eta.y2$sample), 
                                                 to=min(days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]):max(days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]))
      
      sample_quantiles.y1 <- predicted.eta.y1 %>%
        group_by(sample) %>%
        summarise(p2.5 = quantile(value, prob=0.025),
                  p5 = quantile(value, prob=0.05),
                  p10 = quantile(value, prob=0.1),
                  p25 = quantile(value, prob=0.25),
                  p50 = quantile(value, prob=0.5),
                  mean = mean(value),
                  p75 = quantile(value, prob=0.75),
                  p90 = quantile(value, prob=0.9),
                  p95 = quantile(value, prob=0.95),
                  p97.5 = quantile(value, prob=0.975)) %>%
        ungroup()
      
      sample_quantiles.y2 <- predicted.eta.y2 %>%
        group_by(sample) %>%
        summarise(p2.5 = quantile(value, prob=0.025),
                  p5 = quantile(value, prob=0.05),
                  p10 = quantile(value, prob=0.1),
                  p25 = quantile(value, prob=0.25),
                  p50 = quantile(value, prob=0.5),
                  mean = mean(value),
                  p75 = quantile(value, prob=0.75),
                  p90 = quantile(value, prob=0.9),
                  p95 = quantile(value, prob=0.95),
                  p97.5 = quantile(value, prob=0.975)) %>%
        ungroup()
      
      p.y1 <- ggplot() +
        geom_ribbon(data=sample_quantiles.y1, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
        geom_ribbon(data=sample_quantiles.y1, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
        geom_line(data=sample_quantiles.y1, aes(x=sample, y=mean), color="blue") +
        xlab("day") +
        ylab(paste0("CLR(abundance)")) +
        #theme(axis.title=element_text(size=12)) +
        ggtitle(paste0("Treatment = ",treatment,", Individual = ",unique.animals.y1[[indiv.y1]],", Year = 1, Taxon = ",label))
      p.y1 <- p.y1 + 
        geom_point(data=ground.eta.y1, aes(x=sample, y=value))
      p.y2 <- ggplot() +
        geom_ribbon(data=sample_quantiles.y2, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
        geom_ribbon(data=sample_quantiles.y2, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
        geom_line(data=sample_quantiles.y2, aes(x=sample, y=mean), color="blue") +
        xlab("day") +
        ylab(paste0("CLR(abundance)")) +
        #theme(axis.title=element_text(size=12)) +
        ggtitle(paste0("Treatment = ",treatment,", Individual = ",unique.animals.y2[[indiv.y2]],", Year = 2, Taxon = ",label))
      p.y2 <- p.y2 + 
        geom_point(data=ground.eta.y2, aes(x=sample, y=value))
      
      g <- grid.arrange(p.y1, p.y2, nrow=2)
      if(it < 3) {
        ggsave(paste0("images/paired_sanity_",treatment,"_16S_",lr_coord,".png"), plot=g, units="in", dpi=150, height=6, width=12)
      } else {
        ggsave(paste0("images/paired_sanity_",treatment,"_ITS_",lr_coord,".png"), plot=g, units="in", dpi=150, height=6, width=12)
      }
    }
  }
}

# ============================================================================================================
#   FILTER FOR HIGH-CONFIDENCE CORRELATORS
# ============================================================================================================

# convert sampled covariance matrices to correlation
correlations <- posterior$Sigma
# apply doesn't work here because ???
for(i in 1:resample_it) {
  temp <- correlations[,,i]
  temp <- cov2cor(temp)
  temp[upper.tri(temp, diag=TRUE)] <- NA
  correlations[,,i] <- temp
}
# long-ify
correlations_long <- driver::gather_array(correlations, "value", "taxon1", "taxon2", "sample")
# remove redundant entries and diagonal (tagged as NA)
correlations_long <- correlations_long[complete.cases(correlations_long),]
# define 95% posterior inverval value
correlations_quantiles <- correlations_long %>%
  group_by(taxon1, taxon2) %>%
  summarise(p2.5 = quantile(value, prob=0.025),
            mean = mean(value),
            p97.5 = quantile(value, prob=0.975)) %>%
  ungroup()
head(correlations_quantiles)

high_conf_idx <- which(sign(correlations_quantiles$p2.5) == sign(correlations_quantiles$p97.5))

corr_subset <- correlations_quantiles[high_conf_idx,]
if(evaluate == "both") {
  corr_subset <- corr_subset[((corr_subset$taxon1 <= D.16S & corr_subset$taxon2 > D.16S) | (corr_subset$taxon1 > D.16S & corr_subset$taxon2 <= D.16S)),]
}
cat("Average association strength:",round(mean(abs(corr_subset$mean)),3),"\n")
cat("Max association strength:",round(max(abs(corr_subset$mean)),3),"\n")

get_16S_label <- function(idx) {
  label_level <- max(which(sapply(all_data$`16S`$year1$filtered[idx,], is.character)))
  label_content <- NA
  while(is.na(label_content)) {
    label_content <- all_data$`16S`$year1$filtered[idx,label_level]
    label <- paste(names(all_data$`16S`$year1$filtered[label_level]), label_content)
    label_level <- label_level - 1
  }
  return(label)
}

get_ITS_label <- function(idx) {
  label_level <- max(which(sapply(all_data$`ITS`$year1$filtered[idx,], is.character)))
  label_content <- NA
  while(is.na(label_content)) {
    label_content <- all_data$`ITS`$year1$filtered[idx,label_level]
    label <- paste(names(all_data$`ITS`$year1$filtered[label_level]), label_content)
    label_level <- label_level - 1
  }
  return(label)
}

do_plot <- TRUE
fileout <- paste0("correlations_",evaluate,"_",treatment,".txt")
unlink(fileout)

# output these by name
for(hcf in high_conf_idx) {
  if(evaluate == "16S" | evaluate == "ITS") {
    tax1.idx <- correlations_quantiles[hcf,]$taxon1
    tax2.idx <- correlations_quantiles[hcf,]$taxon2
    if(evaluate == "16S") {
      tax1.label <- get_16S_label(tax1.idx)
      tax2.label <- get_16S_label(tax2.idx)
    } else {
      tax1.label <- get_ITS_label(tax1.idx)
      tax2.label <- get_ITS_label(tax2.idx)
    }
    corr.value <- correlations_quantiles[hcf,]$mean
    if(abs(corr.value) > 0.4) {
      # write to stdout
      cat("Labels:",hcf,",",tax1.label,",",tax2.label,", corr=",corr.value,"\n")
      # write to file
      write(paste0(tax1.label,"\t",tax2.label,"\t",round(corr.value,3)), file=fileout, append=TRUE)
    }
  } else {
    # evaluate both
    tax1.idx <- correlations_quantiles[hcf,]$taxon1
    tax2.idx <- correlations_quantiles[hcf,]$taxon2
    tax1.idx.offset <- tax1.idx
    tax2.idx.offset <- tax2.idx
    if(tax1.idx > D.16S) {
      tax1.idx <- tax1.idx - D.16S
      tax1.type <- "ITS"
      tax1.label <- get_ITS_label(tax1.idx)
    } else {
      tax1.type <- "16S"
      tax1.label <- get_16S_label(tax1.idx)
    }
    if(tax2.idx > D.16S) {
      tax2.idx <- tax2.idx - D.16S
      tax2.type <- "ITS"
      tax2.label <- get_ITS_label(tax2.idx)
    } else {
      tax2.type <- "16S"
      tax2.label <- get_16S_label(tax2.idx)
    }
    
    # get mean correlation
    corr.value <- correlations_quantiles[hcf,]$mean
    if(abs(corr.value) > 0.5) {
      # write to stdout
      cat("Labels:",hcf,",",tax1.label,"(",tax1.type,"),",tax2.label,"(",tax2.type,"), corr=",corr.value,"\n")
      # write to file
      write(paste0(tax1.label," (",tax1.type,")\t",tax2.label," (",tax2.type,")\t",round(corr.value,3)), file=fileout, append=TRUE)
      
      # sanity check plots
      # REUSE FROM ABOVE BUT WE'LL FIX THIS LATER
      
      # indiv.y1 <- sample(unique(ids.y1))[1]
      # indiv.y2 <- sample(unique(ids.y2))[1]
      
      #if(do_plot & (tax1.type != tax2.type)) {
      if(do_plot) {
        # choose a random individual
        for(indiv.y1 in sample(unique(ids.y1))[1]) {
          for(indiv.y2 in sample(unique(ids.y2))[1]) {
            cat("\tPlotting",tax1.type,"coord",tax1.label,"and",tax2.type,"coord",tax2.label,"in",unique.animals.y1[[indiv.y1]],"&",unique.animals.y2[[indiv.y2]],"\n")
          
            # get truth for TAX 1, YEAR 1
            ground.eta.tax1.y1 <- gather_array(eta[,idx.animals.y1[[indiv.y1]]], value, "taxon", "sample")
            ground.eta.tax1.y1 <- ground.eta.tax1.y1[ground.eta.tax1.y1$taxon == tax1.idx.offset,]
            ground.eta.tax1.y1$sample <- plyr::mapvalues(ground.eta.tax1.y1$sample, 
                                                    from=ground.eta.tax1.y1$sample, 
                                                    to=days.y1[idx.animals.y1[[indiv.y1]]])
            
            # get truth for TAX 1, YEAR 2
            ground.eta.tax1.y2 <- gather_array(eta[,idx.animals.y2[[indiv.y2]]], value, "taxon", "sample")
            ground.eta.tax1.y2 <- ground.eta.tax1.y2[ground.eta.tax1.y2$taxon == tax1.idx.offset,]
            ground.eta.tax1.y2$sample <- plyr::mapvalues(ground.eta.tax1.y2$sample, 
                                                    from=ground.eta.tax1.y2$sample, 
                                                    to=days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]) # map 1 2 3 ... to 4001 4002 4003 ...
            
            # get truth for TAX 2, YEAR 1
            ground.eta.tax2.y1 <- gather_array(eta[,idx.animals.y1[[indiv.y1]]], value, "taxon", "sample")
            ground.eta.tax2.y1 <- ground.eta.tax2.y1[ground.eta.tax2.y1$taxon == tax2.idx.offset,]
            ground.eta.tax2.y1$sample <- plyr::mapvalues(ground.eta.tax2.y1$sample, 
                                                   from=ground.eta.tax2.y1$sample, 
                                                   to=days.y1[idx.animals.y1[[indiv.y1]]])
            
            # get truth for TAX 2, YEAR 2
            ground.eta.tax2.y2 <- gather_array(eta[,idx.animals.y2[[indiv.y2]]], value, "taxon", "sample")
            ground.eta.tax2.y2 <- ground.eta.tax2.y2[ground.eta.tax2.y2$taxon == tax2.idx.offset,]
            ground.eta.tax2.y2$sample <- plyr::mapvalues(ground.eta.tax2.y2$sample, 
                                                         from=ground.eta.tax2.y2$sample, 
                                                         to=days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]) # map 1 2 3 ... to 4001 4002 4003 ...
            
            # get predictions for TAX 1, YEAR 1
            predicted.eta.tax1.y1 <- gather_array(predictions[,which(ids.y1 == indiv.y1),], "value", "taxon", "sample", "resample")
            predicted.eta.tax1.y1 <- predicted.eta.tax1.y1[predicted.eta.tax1.y1$taxon == tax1.idx.offset,]
            predicted.eta.tax1.y1$sample <- plyr::mapvalues(predicted.eta.tax1.y1$sample, 
                                                       from=min(predicted.eta.tax1.y1$sample):max(predicted.eta.tax1.y1$sample), 
                                                       to=min(days.y1[idx.animals.y1[[indiv.y1]]]):max(days.y1[idx.animals.y1[[indiv.y1]]]))
            
            # get predictions for TAX 1, YEAR 2
            predicted.eta.tax1.y2 <- gather_array(predictions[,length(ids.y1) + which(ids.y2 == indiv.y2),], "value", "taxon", "sample", "resample")
            predicted.eta.tax1.y2 <- predicted.eta.tax1.y2[predicted.eta.tax1.y2$taxon == tax1.idx.offset,]
            predicted.eta.tax1.y2$sample <- plyr::mapvalues(predicted.eta.tax1.y2$sample, 
                                                       from=min(predicted.eta.tax1.y2$sample):max(predicted.eta.tax1.y2$sample), 
                                                       to=min(days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]):max(days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]))
            
            # get predictions for TAX 2, YEAR 1
            predicted.eta.tax2.y1 <- gather_array(predictions[,which(ids.y1 == indiv.y1),], "value", "taxon", "sample", "resample")
            predicted.eta.tax2.y1 <- predicted.eta.tax2.y1[predicted.eta.tax2.y1$taxon == tax2.idx.offset,]
            predicted.eta.tax2.y1$sample <- plyr::mapvalues(predicted.eta.tax2.y1$sample, 
                                                      from=min(predicted.eta.tax2.y1$sample):max(predicted.eta.tax2.y1$sample), 
                                                      to=min(days.y1[idx.animals.y1[[indiv.y1]]]):max(days.y1[idx.animals.y1[[indiv.y1]]]))
            
            # get predictions for TAX 2, YEAR 2
            predicted.eta.tax2.y2 <- gather_array(predictions[,length(ids.y1) + which(ids.y2 == indiv.y2),], "value", "taxon", "sample", "resample")
            predicted.eta.tax2.y2 <- predicted.eta.tax2.y2[predicted.eta.tax2.y2$taxon == tax2.idx.offset,]
            predicted.eta.tax2.y2$sample <- plyr::mapvalues(predicted.eta.tax2.y2$sample, 
                                                            from=min(predicted.eta.tax2.y2$sample):max(predicted.eta.tax2.y2$sample), 
                                                            to=min(days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]):max(days.y2[idx.animals.y2[[indiv.y2]]-nsamp.y1]))
            
            
            sample_quantiles.tax1.y1 <- predicted.eta.tax1.y1 %>%
              group_by(sample) %>%
              summarise(p2.5 = quantile(value, prob=0.025),
                        p5 = quantile(value, prob=0.05),
                        p10 = quantile(value, prob=0.1),
                        p25 = quantile(value, prob=0.25),
                        p50 = quantile(value, prob=0.5),
                        mean = mean(value),
                        p75 = quantile(value, prob=0.75),
                        p90 = quantile(value, prob=0.9),
                        p95 = quantile(value, prob=0.95),
                        p97.5 = quantile(value, prob=0.975)) %>%
              ungroup()
        
            sample_quantiles.tax1.y2 <- predicted.eta.tax1.y2 %>%
              group_by(sample) %>%
              summarise(p2.5 = quantile(value, prob=0.025),
                        p5 = quantile(value, prob=0.05),
                        p10 = quantile(value, prob=0.1),
                        p25 = quantile(value, prob=0.25),
                        p50 = quantile(value, prob=0.5),
                        mean = mean(value),
                        p75 = quantile(value, prob=0.75),
                        p90 = quantile(value, prob=0.9),
                        p95 = quantile(value, prob=0.95),
                        p97.5 = quantile(value, prob=0.975)) %>%
              ungroup()
            
            sample_quantiles.tax2.y1 <- predicted.eta.tax2.y1 %>%
              group_by(sample) %>%
              summarise(p2.5 = quantile(value, prob=0.025),
                        p5 = quantile(value, prob=0.05),
                        p10 = quantile(value, prob=0.1),
                        p25 = quantile(value, prob=0.25),
                        p50 = quantile(value, prob=0.5),
                        mean = mean(value),
                        p75 = quantile(value, prob=0.75),
                        p90 = quantile(value, prob=0.9),
                        p95 = quantile(value, prob=0.95),
                        p97.5 = quantile(value, prob=0.975)) %>%
              ungroup()
            
            sample_quantiles.tax2.y2 <- predicted.eta.tax2.y2 %>%
              group_by(sample) %>%
              summarise(p2.5 = quantile(value, prob=0.025),
                        p5 = quantile(value, prob=0.05),
                        p10 = quantile(value, prob=0.1),
                        p25 = quantile(value, prob=0.25),
                        p50 = quantile(value, prob=0.5),
                        mean = mean(value),
                        p75 = quantile(value, prob=0.75),
                        p90 = quantile(value, prob=0.9),
                        p95 = quantile(value, prob=0.95),
                        p97.5 = quantile(value, prob=0.975)) %>%
              ungroup()
            
            p.t1.y1 <- ggplot() +
              geom_ribbon(data=sample_quantiles.tax1.y1, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
              geom_ribbon(data=sample_quantiles.tax1.y1, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
              geom_line(data=sample_quantiles.tax1.y1, aes(x=sample, y=mean), color="blue") +
              xlab("day") +
              ylab(paste0("CLR(abundance)")) +
              ggtitle(paste0("Treatment = ",treatment,", Individual = ",unique.animals.y1[[indiv.y1]],", Year = 1, Taxon = ",tax1.label," (",tax1.type,")")) +
              theme(axis.title=element_text(size=10),
                    plot.title=element_text(size=10))
            p.t1.y1 <- p.t1.y1 + 
              geom_point(data=ground.eta.tax1.y1, aes(x=sample, y=value))
            
            p.t2.y1 <- ggplot() +
              geom_ribbon(data=sample_quantiles.tax2.y1, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
              geom_ribbon(data=sample_quantiles.tax2.y1, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
              geom_line(data=sample_quantiles.tax2.y1, aes(x=sample, y=mean), color="blue") +
              xlab("day") +
              ylab(paste0("CLR(abundance)")) +
              ggtitle(paste0("Treatment = ",treatment,", Individual = ",unique.animals.y1[[indiv.y1]],", Year = 1, Taxon = ",tax2.label," (",tax2.type,")")) +
              theme(axis.title=element_text(size=10),
                    plot.title=element_text(size=10))
            p.t2.y1 <- p.t2.y1 + 
              geom_point(data=ground.eta.tax2.y1, aes(x=sample, y=value))
            
            p.t1.y2 <- ggplot() +
              geom_ribbon(data=sample_quantiles.tax1.y2, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
              geom_ribbon(data=sample_quantiles.tax1.y2, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
              geom_line(data=sample_quantiles.tax1.y2, aes(x=sample, y=mean), color="blue") +
              xlab("day") +
              ylab(paste0("CLR(abundance)")) +
              ggtitle(paste0("Treatment = ",treatment,", Individual = ",unique.animals.y2[[indiv.y2]],", Year = 2, Taxon = ",tax1.label," (",tax1.type,")")) +
              theme(axis.title=element_text(size=10),
                    plot.title=element_text(size=10))
            p.t1.y2 <- p.t1.y2 + 
              geom_point(data=ground.eta.tax1.y2, aes(x=sample, y=value))
            
            p.t2.y2 <- ggplot() +
              geom_ribbon(data=sample_quantiles.tax2.y2, aes(x=sample, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
              geom_ribbon(data=sample_quantiles.tax2.y2, aes(x=sample, ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
              geom_line(data=sample_quantiles.tax2.y2, aes(x=sample, y=mean), color="blue") +
              xlab("day") +
              ylab(paste0("CLR(abundance)")) +
              ggtitle(paste0("Treatment = ",treatment,", Individual = ",unique.animals.y2[[indiv.y2]],", Year = 2, Taxon = ",tax2.label," (",tax2.type,")")) +
              theme(axis.title=element_text(size=10),
                    plot.title=element_text(size=10))
            p.t2.y2 <- p.t2.y2 + 
              geom_point(data=ground.eta.tax2.y2, aes(x=sample, y=value))
            
            g <- grid.arrange(p.t1.y1, p.t1.y2,
                              p.t2.y1, p.t2.y2, nrow=2)
            ggsave(paste0("images/paired_highconf_",hcf,"_",indiv.y1,"x",indiv.y2,"_",treatment,".png"), plot=g, units="in", dpi=150, height=3, width=12)
          }
        }
      }
    }
  }
}



