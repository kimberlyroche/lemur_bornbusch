library(driver)
library(ggplot2)
library(vegan)
library(tidyverse)
library(RColorBrewer)
library(stringr)
library(stray)
library(driver)
library(ggraph)
library(igraph)

# pre-treatment:
# (1) Renamed 16S_*_fecal_metdata.txt to 16S_*_fecal_metadata.txt
# (2) removed pound-sign from headers in *_Y1_fecal_metadata.txt, *_Y2_fecal_metadata.txt,
#     Bornbusch_*_Y1_fecal_count-table.tsv, and Bornbusch_*_Y2_fecal_count-table.tsv
# (3) Added header labels TAX and CONF to last two columns of Bornbusch_*_Y1_fecal_count-table.tsv
#     and Bornbusch_*_Y2_fecal_count-table.tsv
# (4) Removed extra empty lines from end of ITS_*_fecal_metadata.txt
# (5) In ITS_*_fecal_metadata.txt changed "Sample ID" header to "SampleID"

# ============================================================================================================
#   FUNCTIONS
# ============================================================================================================

get_sample_metadata <- function(sample_ID, attribute, metadata) {
  # sample_ID e.g. LCAX001
  metadata[metadata$SampleID == sample_ID,attribute]
}

get_taxon_identity_from_str <- function(OTU_ID, tax) {
  str_pieces <- strsplit(as.character(tax[[OTU_ID]]), "(;?)D_\\d__", perl=TRUE)[[1]]
  str_pieces <- str_pieces[2:length(str_pieces)]
  if(length(str_pieces) == 7) {
    return(paste0("genus ",str_pieces[6]," :: species ",str_pieces[7]))
  }
  if(length(str_pieces) == 6) {
    return(paste0("family ",str_pieces[5]," :: genus ",str_pieces[6]))
  }
  if(length(str_pieces) == 5) {
    return(paste0("order ",str_pieces[4]," :: family ",str_pieces[5]))
  }
  if(length(str_pieces) == 4) {
    return(paste0("class ",str_pieces[3]," :: order ",str_pieces[4]))
  }
  if(length(str_pieces) == 3) {
    return(paste0("phylum ",str_pieces[2]," :: class ",str_pieces[3]))
  }
  if(length(str_pieces) == 2) {
    return(paste0("kingdom ",str_pieces[1]," :: phylum ",str_pieces[2]))
  }
  if(length(str_pieces) == 1) {
    return(paste0("kingdom ",str_pieces[1]))
  }
  return("NA")
}

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

read_metadata <- function(year=1) {
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

get_counts <- function(df) {
  max_tax_id <- max(which(sapply(df[1,], is.numeric) == FALSE))
  return(df[,(max_tax_id+1):ncol(df)])
}

# input is counts of a single taxon over n samples
filter_low_taxa_sub <- function(x, min_abundance=10, min_representation=0.2) {
  sum(x >= min_abundance)/length(x) >= min_representation
}

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
  tax_pieces <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_pieces <- tax_pieces[1:which(tax_pieces == level)]
  tax_pieces <- c(tax_pieces, "experiment")
  min_sample_idx <- min(which(sapply(data.y1$raw_data[1,], is.numeric)))
  max_sample_idx <- max(which(sapply(data.y1$raw_data[1,], is.numeric)))
  bacteria.agg.y1 <- as.data.frame(data.y1$raw_data %>%
                                     gather(experiment, value, min_sample_idx:max_sample_idx) %>%
                                     group_by(!!!syms(tax_pieces)) %>%
                                     summarise(agg_val = sum(value)) %>%
                                     spread(experiment, agg_val))
  min_sample_idx <- min(which(sapply(data.y2$raw_data[1,], is.numeric)))
  max_sample_idx <- max(which(sapply(data.y2$raw_data[1,], is.numeric)))
  bacteria.agg.y2 <- as.data.frame(data.y2$raw_data %>%
                                     gather(experiment, value, min_sample_idx:max_sample_idx) %>%
                                     group_by(!!!syms(tax_pieces)) %>%
                                     summarise(agg_val = sum(value)) %>%
                                     spread(experiment, agg_val))
  bacteria.agg <- cbind(bacteria.agg.y1,
                        bacteria.agg.y2[,(max(which(sapply(bacteria.agg.y2[1,], is.numeric) == FALSE)) + 1):ncol(bacteria.agg.y2)])
  counts <- get_counts(bacteria.agg)
  level_idx <- length(tax_pieces)-1 # genus=6, species=7, etc.
  if(is.null(filter_percent)) {
    retain_idx <- 1:nrow(counts)
    retained_tax <- bacteria.agg[retain_idx,1:level_idx]
    bacteria.chopped <- counts[retain_idx,]
  } else {
    retain_idx <- filter_low_taxa_alt(counts, percent_of_whole=filter_percent)
    retained_tax <- bacteria.agg[retain_idx,1:level_idx]
    # collapse "Other"; if there are D taxa, retained tax will be D-1 items long
    #   the last item (or row in the counts table) is the assumed "Other" group
    bacteria.chopped <- counts[retain_idx,]
    bacteria.chopped[nrow(bacteria.chopped)+1,] <- apply(counts[!retain_idx,], 2, sum)
    retained_tax[nrow(retained_tax)+1,] <- "Other"
  }
  return(list(counts=bacteria.chopped, level=level_idx, tax=retained_tax))
}

plot_logratios <- function(data, tax, metadata, animal_list, tax_list, save_slug) {
  plot_df_full <- NULL
  first_samples <- metadata %>%
    group_by(Animal) %>%
    filter(Day == min(Day)) %>% 
    filter(1:n() == 1) %>%
    select(c(Animal, SampleID))
  for(lemur in animal_list) {
    l_samples <- metadata[metadata$Animal == lemur,]$SampleID
    l_first <- as.character(first_samples[first_samples$Animal == lemur,]$SampleID)
    l_counts <- matrix(NA, nrow(data), length(l_samples))
    for(i in 1:length(l_samples)) {
      l_counts[,i] <- data[,as.character(l_samples[i])]
    }
    l_clr <- driver::clr(t(l_counts) + 0.5)
    
    df <- gather_array(l_clr, value, sample, taxon)
    df$taxon <- as.factor(df$taxon)
    levels(df$taxon) <- c(apply(tax[,3:5], 1, function(x) paste(x, collapse=" \n ")), "Other")
    
    plot_df <- df[df$taxon %in% levels(df$taxon)[tax_idx],]
    plot_df$animal <- lemur
    if(is.null(plot_df_full)) {
      plot_df_full <- plot_df
    } else {
      plot_df_full <- rbind(plot_df_full, plot_df)
    }
  }
  
  p <- ggplot(plot_df_full, aes(x=sample, y=value)) +
    geom_point(size=2) +
    facet_grid(animal ~ taxon) +
    geom_smooth(aes(x=sample, y=value), method="loess")
  ggsave(paste0("images/",which_data,"/lr_",save_slug,".png"), plot=p, dpi=100, height=6, width=15, units="in")
}

# lifted from Rules of Life code

default_ALR_prior <- function(D, log_var_scale=1) {
  upsilon <- D-1+10 # lesser certainty
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  Xi <- GG%*%(diag(D)*log_var_scale)%*%t(GG) # take diag as covariance over log abundances
  Xi <- Xi*(upsilon-D-1)
  return(list(upsilon=upsilon, Xi=Xi))
}

# ============================================================================================================
#   DATA PARSING
# ============================================================================================================

which_data <- "16S"
#which_data <- "ITS"

# global vars
dist_metric <- "Aitchison" # BC isn't so distinct but maybe count normalization needs to happen
save_images <- TRUE

if(which_data == "16S") {
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

# returned as a tibble with columns: kingdom phylum class order family genus species LCAX001 LCAX002 LCAX003 LCAX004 ...
year.all <- read_sequence_counts(which_data=which_data)
tax <- year.all$tax
year1 <- list(raw_data=year.all$raw_year1)
year2 <- list(raw_data=year.all$raw_year2)
year1$metadata <- read_metadata(year=1)
year2$metadata <- read_metadata(year=2)

# don't agglomerate, filter by minimum common abundance
# filtered.all <- filter_low_taxa(year1$raw_data, year2$raw_data)
# year1$filtered <- filtered.all$year1
# year2$filtered <- filtered.all$year2

# agglomerate, filter by percent abundance
level <- "genus"
filtered.all <- agglomerate_data(year1, year2, level=level, filter_percent=1)
year1$filtered <- cbind(filtered.all$tax, filtered.all$counts[,colnames(filtered.all$counts) %in% year1$metadata$SampleID])
year2$filtered <- cbind(filtered.all$tax, filtered.all$counts[,colnames(filtered.all$counts) %in% year2$metadata$SampleID])

zero.y1.before <- sum(get_counts(year1$raw_data) == 0)/(nrow(get_counts(year1$raw_data))*ncol(get_counts(year1$raw_data)))
zero.y2.before <- sum(get_counts(year2$raw_data) == 0)/(nrow(get_counts(year2$raw_data))*ncol(get_counts(year2$raw_data)))
zero.y1.after <- sum(get_counts(year1$filtered) == 0)/(nrow(get_counts(year1$filtered))*ncol(get_counts(year1$filtered)))
zero.y2.after <- sum(get_counts(year2$filtered) == 0)/(nrow(get_counts(year2$filtered))*ncol(get_counts(year2$filtered)))
cat("Year 1 zeros:",round(zero.y1.before*100, 2),"% ->",round(zero.y1.after*100, 2),"%\n")
cat("Year 2 zeros:",round(zero.y2.before*100, 2),"% ->",round(zero.y2.after*100, 2),"%\n")

# ============================================================================================================
#   ORDINATION PLOTS
# ============================================================================================================

year1$clr <- t(driver::clr(t(get_counts(year1$filtered)) + 0.5)) # returned as taxa x samples
year2$clr <- t(driver::clr(t(get_counts(year2$filtered)) + 0.5))

if(dist_metric == "Aitchison") {
  dist.y1 <- dist(t(year1$clr))
  dist.y2 <- dist(t(year2$clr))
}

if(dist_metric == "BC") {
  counts.y1 <- get_counts(year1$filtered)
  counts.y2 <- get_counts(year2$filtered)
  counts.all <- cbind(counts.y1, counts.y2)
  raremax.all <- min(apply(counts.all, 2, sum))
  year1.sampled <- suppressWarnings(rrarefy(counts.y1, raremax.all))
  year2.sampled <- suppressWarnings(rrarefy(counts.y2, raremax.all))
  dist.y1 <- vegdist(t(year1.sampled), method="bray")
  dist.y2 <- vegdist(t(year2.sampled), method="bray")
}

y1.embed <- cmdscale(dist.y1)
y2.embed <- cmdscale(dist.y2)

# get time (days) from baseline
#samples <- colnames(bacteria.filtered)[8:ncol(bacteria.filtered)]
#sample_times <- lapply(samples, function(x) strptime(get_sample_metadata(x, "Date"), format="%d-%b-%y"))
#sample_times_ISO <- sapply(sample_times, function(x) as.numeric(as.POSIXct(x)))
#baseline_idx <- which(sample_times_ISO == min(sample_times_ISO))
#baseline_time <- sample_times[[baseline_idx]]
#sample_times <- sapply(sample_times, function(x) difftime(x, baseline_time, units="days"))

df.y1 <- data.frame(x=y1.embed[,1], y=y1.embed[,2], 
                    animal=year1$metadata$Animal,
                    treatment=year1$metadata$Treatment)
df.y2 <- data.frame(x=y2.embed[,1], y=y2.embed[,2], 
                    animal=year2$metadata$Animal,
                    treatment=year2$metadata$Treatment)

p <- ggplot(df.y1) +
  geom_point(aes(x=x, y=y, col=treatment), size=3) +
  xlab("PCoA 1") +
  ylab("PCoA 2")
show(p)
if(save_images) {
  ggsave(paste0("images/",which_data,"/ordination_by_treatment_Y1.png"), plot=p, dpi=100, width=8, height=6, units="in")
}

p <- ggplot(df.y2) +
  geom_point(aes(x=x, y=y, col=treatment), size=3) +
  xlab("PCoA 1") +
  ylab("PCoA 2")
show(p)
if(save_images) {
  ggsave(paste0("images/",which_data,"/ordination_by_treatment_Y2.png"), plot=p, dpi=100, width=8, height=6, units="in")
}

p <- ggplot(df.y1) +
  geom_point(aes(x=x, y=y, col=animal), size=3) +
  xlab("PCoA 1") +
  ylab("PCoA 2")
show(p)
if(save_images) {
  ggsave(paste0("images/",which_data,"/ordination_by_animal_Y1.png"), plot=p, dpi=100, width=8, height=6, units="in")
}

p <- ggplot(df.y2) +
  geom_point(aes(x=x, y=y, col=animal), size=3) +
  xlab("PCoA 1") +
  ylab("PCoA 2")
show(p)
if(save_images) {
  ggsave(paste0("images/",which_data,"/ordination_by_animal_Y2.png"), plot=p, dpi=100, width=8, height=6, units="in")
}

# ============================================================================================================
#   ALPHA-DIVERSITY PLOTS
# ============================================================================================================

df.alpha <- data.frame(time=c(), shannon=c(), animal=c())
counts.y1 <- get_counts(year1$filtered)
for(sample in colnames(counts.y1)) {
  sample_data <- counts.y1[,sample]
  time_offset <- get_sample_metadata(sample, "Day", year1$metadata)
  animal_label <- as.character(get_sample_metadata(sample, "Animal", year1$metadata))
  treatment_label <- as.character(get_sample_metadata(sample, "Treatment", year1$metadata))
  df.alpha <- rbind(df.alpha, data.frame(time=time_offset, shannon=vegan::diversity(sample_data),
                                         animal=animal_label, treatment=treatment_label, year=1))
}
counts.y2 <- get_counts(year2$filtered)
for(sample in colnames(counts.y2)) {
  sample_data <- counts.y2[,sample]
  time_offset <- get_sample_metadata(sample, "Day", year2$metadata)
  animal_label <- as.character(get_sample_metadata(sample, "Animal", year2$metadata))
  treatment_label <- as.character(get_sample_metadata(sample, "Treatment", year2$metadata))
  df.alpha <- rbind(df.alpha, data.frame(time=time_offset, shannon=vegan::diversity(sample_data),
                                         animal=animal_label, treatment=treatment_label, year=2))
}

p <- ggplot(df.alpha[df.alpha$year == 1,]) +
  geom_line(aes(x=time, y=shannon, group=animal, color=treatment)) +
  geom_point(aes(x=time, y=shannon, group=animal, color=treatment))
show(p)
if(save_images) {
  ggsave(paste0("images/",which_data,"/alphadiversity_time_Y1.png"), plot=p, dpi=100, width=8, height=5, units="in")
}

p <- ggplot(df.alpha[df.alpha$year == 2,]) +
  geom_line(aes(x=time, y=shannon, group=animal, color=treatment)) +
  geom_point(aes(x=time, y=shannon, group=animal, color=treatment))
show(p)
if(save_images) {
  ggsave(paste0("images/",which_data,"/alphadiversity_time_Y2.png"), plot=p, dpi=100, width=8, height=5, units="in")
}

# ============================================================================================================
#   PROPORTIONAL PLOTS
# ============================================================================================================

agg_levels <- c("family")

for(agg_level in agg_levels) {
  agg_obj <- agglomerate_data(year1, year2, level=agg_level, filter_percent=1)
  
  counts.chopped.y1 <- agg_obj$counts[,colnames(agg_obj$counts) %in% year1$metadata$SampleID]
  prop.y1 <- apply(counts.chopped.y1, 2, function(x) { x / sum(x) })

  df.y1 <- data.frame(index=c(), day=c(), taxon=c(), value=c(), animal=c())
  for(sample_idx in 1:ncol(counts.chopped.y1)) {
    sample <- colnames(counts.chopped.y1)[sample_idx]
    sample_data <- prop.y1[,sample]
    time_offset <- get_sample_metadata(sample, "Day", year1$metadata)
    animal_label <- get_sample_metadata(sample, "Animal", year1$metadata)
    for(taxon_idx in 1:length(sample_data)) {
      taxon_label_readable <- NULL
      for(t in agg_obj$level:1) {
        taxon_label_readable <- agg_obj$tax[taxon_idx,t]
        if(!is.na(taxon_label_readable)) {
          break
        }
      }
      taxon_value <- sample_data[taxon_idx]
      df.y1 <- rbind(df.y1, data.frame(index=sample_idx,
                                       day=time_offset,
                                       taxon=taxon_label_readable,
                                       value=taxon_value,
                                       animal=animal_label))
    }
  }
  
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df.y1$taxon)))
  coul[which(unique(df.y1$taxon) == "Other")] <- "#AAAAAA"
  
  plot_animals <- unique(year1$metadata$Animal)
  for(plot_animal in plot_animals) {
    animal_subset <- df.y1[df.y1$animal==plot_animal,]
    animal_subset$index <- as.factor(animal_subset$index)
    levels(animal_subset$index) <- 1:length(animal_subset$index)
    p <- ggplot(animal_subset, aes(x=index, y=value, fill=taxon)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values=coul)
    #show(p)
    if(save_images) {
      if(which_data == "16S") {
        im_height <- 4
      } else {
        im_height <- 6
      }
      ggsave(paste0("images/",which_data,"/timecourse_",agg_level,"_",plot_animal,"_Y1.png"), plot=p, dpi=100, width=10, height=im_height, units="in")
    }
  }
  
  counts.chopped.y2 <- agg_obj$counts[,colnames(agg_obj$counts) %in% year2$metadata$SampleID]
  prop.y2 <- apply(counts.chopped.y2, 2, function(x) { x / sum(x) })
  
  df.y2 <- data.frame(index=c(), day=c(), taxon=c(), value=c(), animal=c())
  for(sample_idx in 1:ncol(counts.chopped.y2)) {
    sample <- colnames(counts.chopped.y2)[sample_idx]
    sample_data <- prop.y2[,sample]
    time_offset <- get_sample_metadata(sample, "Day", year2$metadata)
    animal_label <- get_sample_metadata(sample, "Animal", year2$metadata)
    for(taxon_idx in 1:length(sample_data)) {
      taxon_label_readable <- NULL
      for(t in agg_obj$level:1) {
        taxon_label_readable <- agg_obj$tax[taxon_idx,t]
        if(!is.na(taxon_label_readable)) {
          break
        }
      }
      taxon_value <- sample_data[taxon_idx]
      df.y2 <- rbind(df.y2, data.frame(index=sample_idx, day=time_offset, taxon=taxon_label_readable, value=taxon_value, animal=animal_label))
    }
  }
  
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df.y2$taxon)))
  coul[which(unique(df.y2$taxon) == "Other")] <- "#AAAAAA"
  
  plot_animals <- unique(year2$metadata$Animal)
  for(plot_animal in plot_animals) {
    animal_subset <- df.y2[df.y2$animal==plot_animal,]
    animal_subset$index <- as.factor(animal_subset$index)
    levels(animal_subset$index) <- 1:length(animal_subset$index)
    p <- ggplot(animal_subset, aes(x=index, y=value, fill=taxon)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values=coul)
    #show(p)
    if(save_images) {
      if(which_data == "16S") {
        im_height <- 4
      } else {
        im_height <- 6
      }
      ggsave(paste0("images/",which_data,"/timecourse_",agg_level,"_",plot_animal,"_Y2.png"), plot=p, dpi=100, width=10, height=4, units="in")
    }
  }
}

# ============================================================================================================
#   PLOT SOME TIME EVOLVING LOG RELATIVE ABUNDANCES
# ============================================================================================================

animals <- unique(year1$metadata$Animal)
CON_animals.y1 <- unique(year1$metadata[year1$metadata$Treatment == "CON",]$Animal)
ABX_animals.y1 <- unique(year1$metadata[year1$metadata$Treatment == "ABX",]$Animal)
ABXFT_animals.y1 <- unique(year1$metadata[year1$metadata$Treatment == "ABXFT",]$Animal)

agg_obj <- agglomerate_data(year1, year2, level="family", filter_percent=0.05)
#bacteria.chopped.y1 <- agg_obj$counts[,colnames(agg_obj$counts) %in% year1$metadata$SampleID]

if(which_data == "16S") {
  idx <- c(39, 50, 52, 104, 117, 120, 174, 225)
} else {
  idx <- sample(1:nrow(agg_obj$tax))[1:8]
}
tax_idx <- which(rownames(agg_obj$tax) %in% idx)

plot_logratios(agg_obj$counts, agg_obj$tax, year1$metadata, CON_animals.y1, tax_idx, "CON_Y1")
plot_logratios(agg_obj$counts, agg_obj$tax, year1$metadata, ABX_animals.y1, tax_idx, "ABX_Y1")
plot_logratios(agg_obj$counts, agg_obj$tax, year1$metadata, ABXFT_animals.y1, tax_idx, "ABXFT_Y1")

# ============================================================================================================
#   FIT BASSET TO EACH EXPERIMENT INDIVIDUALLY
# ============================================================================================================

# collapse a bit for computational ease
if(!exists("filtered.all")) {
  filter_percent <- 1
  if(which_data == "ITS") {
    # MUCH sparser
    filter_percent <- 0.25
  }
  filtered.all <- agglomerate_data(year1, year2, level="genus", filter_percent=filter_percent)
}

# so the average correlation between samples is ~= 0
if(FALSE) {
  for(animal in unique(year1$metadata$Animal)) {
    ilr.samples <- driver::ilr(t(filtered.all$counts[,year1$metadata$SampleID]) + 0.5)
    animal.samples <- ilr.samples[year1$metadata$Animal == animal,]
    animal.scaled <- scale(t(animal.samples), center=T, scale=F) # K samples x D-1 taxa
    animal.cov <- cov(t(animal.scaled))
    animal.corr <- cov2cor(animal.cov)
    image(animal.corr)
    cat("Avg. for",animal,":",mean(animal.corr[upper.tri(animal.corr, diag=F)]),"\n")
  }
}

treatment <- "CON"
#treatment <- "ABX"
#treatment <- "ABXFT"

# get the sample IDs we want and animal names assoc. with them
sampleIDs.y1 <- year1$metadata[year1$metadata$Treatment == treatment,]$SampleID
sampleIDs.y2 <- year2$metadata[year2$metadata$Treatment == treatment,]$SampleID
animals.y1 <- year1$metadata[year1$metadata$Treatment == treatment,]$Animal
animals.y2 <- year2$metadata[year2$metadata$Treatment == treatment,]$Animal
days.y1 <- year1$metadata[year1$metadata$Treatment == treatment,]$Day
days.y2 <- year2$metadata[year2$metadata$Treatment == treatment,]$Day

# re-order all for clarity (1) by animal
new_order.y1 <- order(animals.y1)
new_order.y2 <- order(animals.y2)
sampleIDs.y1 <- sampleIDs.y1[new_order.y1]
sampleIDs.y2 <- sampleIDs.y2[new_order.y2]
animals.y1 <- animals.y1[new_order.y1]
animals.y2 <- animals.y2[new_order.y2]
days.y1 <- days.y1[new_order.y1]
days.y2 <- days.y2[new_order.y2]

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

nsamp.y1 <- length(sampleIDs.y1)

control.y1 <- filtered.all$counts[,sampleIDs.y1]
cat("Year 1:\n")
cat("\tPercent zeros:",round(sum(control.y1 == 0)/(nrow(control.y1)*ncol(control.y1))*100),"% zeros\n")
cat("\tTaxa x samples:",nrow(control.y1),"x",ncol(control.y1),"\n")

control.y2 <- filtered.all$counts[,sampleIDs.y2]
cat("Year 2:\n")
cat("\tPercent zeros:",round(sum(control.y2 == 0)/(nrow(control.y2)*ncol(control.y2))*100),"% zeros\n")
cat("\tTaxa x samples:",nrow(control.y2),"x",ncol(control.y2),"\n")

Y <- as.matrix(cbind(control.y1, control.y2))
D <- nrow(Y)

alr_ys <- driver::alr((t(Y) + 0.5)) # returned as: N samples x D-1 taxa

clr_ys <- driver::clr((t(Y) + 0.5))

# build kernels for basset
N <- ncol(Y)

unique.animals.y1 <- unique(as.character(animals.y1))
unique.animals.y2 <- unique(as.character(animals.y2))

unique.animals.y1
unique.animals.y2

X <- matrix(NA, 2, N)
X_predict <- NULL
new_sample_order <- c()
bump <- 1000 # hack; days to separate each individual's series by
            # such that they're effectively independent
offset <- 0
idx.animals.y1 <- list()
days.animals.y1 <- list()
alr_means.y1 <- list()
for(a in 1:length(unique.animals.y1)) {
  animal <- as.character(unique.animals.y1)[a]
  bump.animal <- (a-1)*bump
  idx.animal.y1 <- which(animals.y1 == animal)
  alr_means.y1[[animal]] <- colMeans(alr_ys[idx.animal.y1,])
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
alr_means.y2 <- list()
for(aa in 1:length(unique.animals.y2)) {
  animal <- as.character(unique.animals.y2)[aa]
  bump.animal <- ((a + aa)-1)*bump
  idx.animal.y2 <- which(animals.y2 == animal)
  alr_means.y2[[animal]] <- colMeans(alr_ys[(idx.animal.y2+nsamp.y1),])
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

idx.animals.y1
idx.animals.y2

days.animals.y1
days.animals.y2

alr_means.y1
alr_means.y2

dc <- 0.01 # desired minimum correlation
if(which_data == "16S") {
  dd_se <- 10
  rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
  g_sigma <- 1
} else {
  dd_se <- 10
  rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
  g_sigma <- 2
}
if(treatment == "ABX" | treatment == "ABXFT") {
  if(which_data == "16S") {
    dd_se <- 10
    rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
    g_sigma <- 1.5
  } else {
    dd_se <- 10
    rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
    g_sigma <- 3
  }
}

Gamma <- function(X) {
  SE(X[1,,drop=F], sigma=g_sigma, rho=rho_se, jitter=1e-10) # time only
}

prior_obj <- default_ALR_prior(D)

ids.y1 <- X[2,1:nsamp.y1]
ids.y2 <- X[2,(nsamp.y1+1):N]
ids.y2 <- ids.y2 - max(ids.y1)
Theta.assignments <- c(paste0(unique.animals.y1[ids.y1],"Y1"), paste0(unique.animals.y2[ids.y2],"Y2"))
alr_means <- list()
for(animal in unique.animals.y1[ids.y1]) {
  alr_means[[paste0(animal,"Y1")]] <- alr_means.y1[[animal]]
}
for(animal in unique.animals.y2[ids.y2]) {
  alr_means[[paste0(animal,"Y2")]] <- alr_means.y2[[animal]]
}

Theta <- function(X) {
  # need these things as global:
  #  Theta.assignments: per-sample label e.g. "AracusY1"
  #  alr_means: named list (labeled as "AracusY1" etc.) of mean alr
  temp <- matrix(NA, D-1, ncol(X))
  for(ta in 1:length(Theta.assignments)) {
    a <- Theta.assignments[ta]
    temp[,ta] <- alr_means[[a]]
  }
  return(temp)
}

# fit on all
fit <- stray::basset(Y, X, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi, n_samples=100)

# need the prediction year assignments

ids.y1 <- X_predict[2,1:nsamp.predict.y1]
ids.y2 <- X_predict[2,(nsamp.predict.y1+1):ncol(X_predict)]
ids.y2 <- ids.y2 - max(ids.y1)
Theta.assignments.append <- c(paste0(unique.animals.y1[ids.y1],"Y1"), paste0(unique.animals.y2[ids.y2],"Y2"))
Theta.assignments <- c(Theta.assignments, Theta.assignments.append)

fit.clr <- to_clr(fit)

# png(paste0("C:/Users/kim/Desktop/model1_",treatment,".png"))
# image(apply(fit.clr$Sigma, c(1,2), mean))
# dev.off()

Sigma_mean_samples <- fit.clr$Sigma
for(i in 1:dim(Sigma_mean_samples)[3]) {
  temp <- Sigma_mean_samples[,,i]
  temp <- cov2cor(temp)
  temp[upper.tri(temp, diag=T)] <- NA
  Sigma_mean_samples[,,i] <- temp
}

Sigma_samples <- gather_array(Sigma_mean_samples, "corr_value", "row", "col", "sample")
Sigma_samples <- Sigma_samples[complete.cases(Sigma_samples),]

# ============================================================================================================
#   RENDER STRONGEST CORRELATORS AS A GRAPH
# ============================================================================================================

post_quantiles <- Sigma_samples %>%
  unite(combo, c("row", "col")) %>%
  group_by(combo) %>%
  summarise(min = min(corr_value),
            corr_value = mean(corr_value),
            max = max(corr_value)) %>%
  ungroup()
post_quantiles <- as.data.frame(post_quantiles)

relations <- post_quantiles[sign(post_quantiles$min) == sign(post_quantiles$max),c("combo", "corr_value")]
relations <- relations %>% separate(combo, into=c("from", "to"))
relations <- relations %>%
  mutate(sign=as.factor(sign(corr_value))) %>%
  mutate(magnitude=abs(corr_value))

# filter to relations greater than some strength
relations <- relations[relations$magnitude >= 0.3,]
# filter to top k edges
#relations <- relations[order(relations$magnitude, decreasing=T)[1:50],]

nodes <- data.frame(name=c(), label=c())
for(i in 1:D) {
  if(level == "family") {
    nodes <- rbind(nodes, data.frame(name=i, label=filtered.all$tax[i,5]))
  }
  if(level == "genus") {
    label <- filtered.all$tax[i,6]
    if(is.na(label)) {
      label <- paste0(filtered.all$tax[i,5],"/NA")
    }
    nodes <- rbind(nodes, data.frame(name=i, label=label))
    
  }
}

# Cytoscape graphs
if(FALSE) {
  rp <- relations
  rp$from <- as.numeric(rp$from)
  rp$to <- as.numeric(rp$to)
  rp$sign <- as.numeric(rp$sign)
  for(nl in unique(c(rp$from, rp$to))){
    rpl_idx <- which(rp$from == nl)
    rp[rpl_idx,]$from <- rep(as.character(nodes[nodes$name == nl,]$label), length(rpl_idx))
    rpl_idx <- which(rp$to == nl)
    rp[rpl_idx,]$to <- rep(as.character(nodes[nodes$name == nl,]$label), length(rpl_idx))
  }
  rp$corr_value <- round(rp$corr_value, 2)
  rp <- rp %>%
    mutate(order = order(from))
  write.table(rp, file=paste0("cytoscape_input_",level,"_",treatment,".tsv"),
              quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
}

# R graphs
g <- graph_from_data_frame(relations, directed=FALSE, vertices=nodes)

# layouts: 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl'
p <- ggraph(g, layout='circle') + 
  geom_edge_link(aes(color=sign, width=magnitude)) +
  #geom_node_point(color="black", size=3, stroke=1) +
  scale_edge_color_manual(values = c("red", "black")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
show(p)
ggsave(paste0("images/",which_data,"/graph_",level,"_",treatment,".png"),
       plot=p,units="in",dpi=300,height=6,width=7)
p <- p +
  geom_node_text(aes(label=label), check_overlap=TRUE)
show(p)
ggsave(paste0("images/",which_data,"/graph_",level,"_",treatment,"_labels.png"),
       plot=p,units="in",dpi=300,height=6,width=7)

# ============================================================================================================
#   GET RANKED LIST OF STRONG CORRELATORS
# ============================================================================================================

predicted <- predict(fit.clr, X_predict, iter=100) # predicts samples from the posterior (default = 2000)

sample_quantiles <- Sigma_samples %>%
  group_by(row, col) %>%
  summarise(p2.5 = quantile(corr_value, prob=0.025),
            p5 = quantile(corr_value, prob=0.05),
            p10 = quantile(corr_value, prob=0.1),
            p25 = quantile(corr_value, prob=0.25),
            p50 = quantile(corr_value, prob=0.5),
            mean = mean(corr_value),
            p75 = quantile(corr_value, prob=0.75),
            p90 = quantile(corr_value, prob=0.9),
            p95 = quantile(corr_value, prob=0.95),
            p97.5 = quantile(corr_value, prob=0.975)) %>%
  ungroup()

head(as.data.frame(sample_quantiles))
temp <- as.data.frame(sample_quantiles)
temp2 <- which((temp$p2.5 > 0 & temp$p97.5 > 0) | (temp$p2.5 < 0 & temp$p97.5 < 0), arr.ind=T)
temp3 <- temp[temp2,]
temp4 <- temp3[order(abs(temp3$mean),decreasing=T),]

readable_tax <- function(tax_pieces) {
  viable <- which(!is.na(tax_pieces))
  upper <- max(viable)
  lower <- max(1, upper-1)
  if(lower == 1) { lower_str <- paste0("Kingdom ",tax_pieces[lower]) }
  if(lower == 2) { lower_str <- paste0("Phylum ",tax_pieces[lower]) }
  if(lower == 3) { lower_str <- paste0("Class ",tax_pieces[lower]) }
  if(lower == 4) { lower_str <- paste0("Order ",tax_pieces[lower]) }
  if(lower == 5) { lower_str <- paste0("Family ",tax_pieces[lower]) }
  if(lower == 6) { lower_str <- paste0("Genus ",tax_pieces[lower]) }
  if(lower == 7) { lower_str <- paste0("Species ",tax_pieces[lower]) }
  if(upper == 1) { upper_str <- paste0("Kingdom ",tax_pieces[upper]) }
  if(upper == 2) { upper_str <- paste0("Phylum ",tax_pieces[upper]) }
  if(upper == 3) { upper_str <- paste0("Class ",tax_pieces[upper]) }
  if(upper == 4) { upper_str <- paste0("Order ",tax_pieces[upper]) }
  if(upper == 5) { upper_str <- paste0("Family ",tax_pieces[upper]) }
  if(upper == 6) { upper_str <- paste0("Genus ",tax_pieces[upper]) }
  if(upper == 7) { upper_str <- paste0("Species ",tax_pieces[upper]) }
  return(paste0(lower_str, "/", upper_str))
}

lines <- c()
max_cor <- c()
max_val <- -Inf
min_cor <- c()
min_val <- Inf
for(sample_idx in 1:20) {
  row <- temp4[sample_idx,]$row
  col <- temp4[sample_idx,]$col
  t1.name <- readable_tax(filtered.all$tax[row,])
  cat(t1.name,"\n")
  t2.name <- readable_tax(filtered.all$tax[col,])
  cat(t2.name,"\n")
  cor_val <- round(temp[temp$row==row & temp$col==col,]$mean, 2)
  if(cor_val > max_val) {
    max_val <- cor_val
    max_cor <- c(row, col)
  }
  if(cor_val < min_val) {
    min_val <- cor_val
    min_cor <- c(row, col)
  }
  cat("\t",cor_val,"\n")
  lines <- c(lines, c(paste0(t1.name,"\t",t2.name,"\t",cor_val)))
}
fileConn <- file(paste0("lab_mtg/rankedlist_",treatment,"_",which_data,".txt"))
writeLines(lines, fileConn)
close(fileConn)

Sigma_set <- fit.clr$Sigma
for(i in 1:100) {
  Sigma_set[,,i] <- cov2cor(Sigma_set[,,i])
}
meanSigma <- apply(Sigma_set, c(1,2), mean)
df <- gather_array(meanSigma, value, row, col)
p <- ggplot(df, aes(row, col, fill=value)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-1,1),
                       colours=c("navyblue", "white", "darkred"))
show(p)
ggsave(paste0("images/",which_data,"/corrmat_",treatment,".png"),plot=p,units="in",dpi=100,height=6,width=7)

# plot a few samples of Sigma's posterior
if(treatment == "CON") {
  sample_idx <- sample(1:100)[1:5]
  for(s in sample_idx) {
    df <- gather_array(cov2cor(fit.clr$Sigma[,,s]), value, row, col)
    p <- ggplot(df, aes(row, col, fill=value)) + 
      geom_tile() +
      scale_fill_gradientn(limits = c(-1,1),
                           colours=c("navyblue", "white", "darkred"))
    #show(p)
    ggsave(paste0("images/",which_data,"/corrmat_",treatment,"_sample",s,".png"),plot=p,units="in",dpi=100,height=6,width=7)
  }
}

log_ratios <- list(max_cor[1], max_cor[2], min_cor[1], min_cor[2])
names(log_ratios) <- c(paste0("(+) corr 1\n",filtered.all$tax[max_cor[1],4],"\n",filtered.all$tax[max_cor[1],5]),
                       paste0("(+) corr 2\n",filtered.all$tax[max_cor[2],4],"\n",filtered.all$tax[max_cor[2],5]),
                       paste0("(-) corr 1\n",filtered.all$tax[min_cor[1],4],"\n",filtered.all$tax[min_cor[1],5]),
                       paste0("(-) corr 2\n",filtered.all$tax[min_cor[2],4],"\n",filtered.all$tax[min_cor[2],5]))

post_quantiles_all <- NULL
df_true_all <- NULL
for(animal.idx in unique(X_predict[2,])) {
  year <- 1
  animal.idx.adj <- animal.idx
  if(animal.idx > max(ids.y1)) {
    year <- 2
    animal.idx.adj <- animal.idx - max(ids.y1)
  }
  if(year == 1) {
    animal <- unique.animals.y1[animal.idx]
  } else {
    animal <- unique.animals.y2[animal.idx.adj]
  }
  cat("Plotting predictions for",animal,"in year",year,"\n")
  plot.idx.predict <- which(X_predict[2,] == animal.idx)
  if(year == 1) {
    plot.idx.true <- idx.animals.y1[[animal]]
  } else {
    plot.idx.true <- idx.animals.y2[[animal]]
  }

  for(lr in names(log_ratios)) {
    lr_idx <- log_ratios[[lr]]
    posterior_samples <- gather_array(predicted[lr_idx,plot.idx.predict,], "LR_value", "observation", "sample_no")
    
    post_quantiles <- posterior_samples %>%
      group_by(observation) %>%
      summarise(p2.5 = quantile(LR_value, prob=0.025),
                p5 = quantile(LR_value, prob=0.05),
                p10 = quantile(LR_value, prob=0.1),
                p25 = quantile(LR_value, prob=0.25),
                p50 = quantile(LR_value, prob=0.5),
                mean = mean(LR_value),
                p75 = quantile(LR_value, prob=0.75),
                p90 = quantile(LR_value, prob=0.9),
                p95 = quantile(LR_value, prob=0.95),
                p97.5 = quantile(LR_value, prob=0.975)) %>%
      ungroup()
  
    animal_label <- paste0(animal," (Y",year,")")
    
    if(is.null(post_quantiles_all)) {
      post_quantiles_all <- as.data.frame(post_quantiles)
      post_quantiles_all <- cbind(post_quantiles_all, animal=animal_label, lr=lr)
    } else {
      post_quantiles_appended <- as.data.frame(post_quantiles)
      post_quantiles_appended <- cbind(post_quantiles_appended, animal=animal_label, lr=lr)
      post_quantiles_all <- rbind(post_quantiles_all, post_quantiles_appended)
    }
    
    if(year == 1) {
      x_span <- days.y1[idx.animals.y1[[animal]]]
      x_span <- x_span - min(x_span) + 1
    } else {
      x_span <- days.y2[idx.animals.y2[[animal]] - nsamp.y1]
      x_span <- x_span - min(x_span) + 1
    }
    #mean_eta <- apply(fit$Eta, c(1,2), mean)
    mean_eta <- apply(fit.clr$Eta, c(1,2), mean)
    df_true <- data.frame(x=x_span, y=mean_eta[lr_idx,plot.idx.true])
    
    if(is.null(df_true_all)) {
      df_true_all <- cbind(df_true, animal=animal_label, lr=lr)
    } else {
      df_true_appended <- cbind(df_true, animal=animal_label, lr=lr)
      df_true_all <- rbind(df_true_all, df_true_appended)
    }
  }
}

p <- ggplot(post_quantiles_all, aes(x=observation, y=mean)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
  geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
  geom_line(color="blue") +
  geom_point(data=df_true_all, aes(x=x, y=y), alpha=0.5) +
  facet_grid(lr ~ animal, scales="free_y")
show(p)
im_width <- 15
if(treatment == "ABX") {
  im_width <- 12
}
if(treatment == "ABXFT") {
  im_width <- 10
}
ggsave(paste0("images/",which_data,"/lr_",treatment,".png"),plot=p,units="in",dpi=100,height=6,width=im_width)




