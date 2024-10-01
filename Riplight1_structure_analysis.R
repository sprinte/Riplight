#### RIPLIGHT1 Analysis of community data ####
#### Last updated: August 16, 2024

#### USER SETTINGS ####

save_output <- FALSE
skip_stats <- TRUE

# determine in- and output directories
dir_in <- paste("D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Data/Structure & Abundance") # where is all your data
dir_out <- paste("D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Results/Structure & Abundance") # where the results and plots are saved
dir_stats <- "Stats"
dir_tables <- "Tables"

# set variables
parameter <- "Fungi" # either "Bacteria", "Fungi", "Algae", or "Diatoms"
project <- "RipLight1" # to name output files

# organism labeling
heterotrophs <- c("bacteria", "fungi")
autotrophs <- c("algae", "diatoms", "cyanobacteria", "green_algae")

# include terrestrial samples? (in plots only)
include.terrestrial <- TRUE # if TRUE, terrestrial samples (only available for heterotrophs) will be included in the NMDS and alpha diversity plots (but not in stats!), if FALSE, they will be dropped

# include initial samples? (for betadiversity analyses)
include.initial <- TRUE

# paths
parameter <- tolower(parameter)

# file names
if (parameter %in% heterotrophs) { # File names for bacteria and fungi
  # sample_names<-as.character(187:222)
  metafile <- "Meta_Sequencing.txt"
  otufile <- paste("ASV_", str_to_title(parameter), "_", project, ".txt", sep = "") # for bacteria and fungi
  taxonomyfile <- paste("Taxonomy_cleaned_", str_to_title(parameter), "_", project, ".txt", sep = "") # for bacteria and fungi
  abundancefile <- "qPCR.txt"
} else if (parameter %in% autotrophs) { # File names for diatoms
  sample_names <- as.character(1:24)
  metafile <- "Meta_algae.txt"
  otufile <- "Morphotypes_no_algae.txt"
  taxonomyfile <- "Taxonomy_algae.txt"
}
envirfile <- "Riplight1_Enviro.txt"

# file names output
filename_vectors <- paste0(paste("beta_vectors", parameter, sep = "_"), ".rds")
filename_nmds <- paste0(paste("beta_nmds", parameter, sep = "_"), ".rds")
filename_venn_full <- paste0(paste("venn_full", parameter, sep = "_"), ".rds")
filename_venn <- paste0(paste("venn", parameter, sep = "_"), ".rds")
filename_dbrda <- paste0(paste("dbRDA", parameter, sep = "_"), ".rds")

## VARIABLE SETTINGS
# variable 1: sediment type
var1 <- "sediment"
levels_var1 <- c("mixed" = "Mixed", "aquatic" = "Aquatic", "terrestrial" = "Terrestrial") # how factor 1 (sediment type) shall be renamed
var1_selected <- c("Aquatic", "Mixed", "Terrestrial") # subset that will be displayed
# if (parameter %in% autotrophs) {
#   var1_selected <- c("Aquatic", "Mixed") # subset that will be displayed
# } else if (parameter %in% heterotrophs) {
#   var1_selected <- c("Aquatic", "Mixed", "Terrestrial") # subset that will be displayed
# }

# variable 2: transport state
var2 <- "transport"
levels_var2 <- c("initial" = "Initial", "migrating" = "Migrating", "stationary" = "Stationary") # how factor 2 (transport regime) shall be renamed
levels_var2a <- c("initial" = "Initial", "final" = "Final")
var2_selected <- c("Initial", "Migrating", "Stationary") # subset that will be displayed

# algae relabelling
algae_levels <- c("diatoma" = "diatoms", "cyanobacteria" = "cyanobacteria", "green algae" = "green_algae")

## CALCULATION SETTINGS
# use_scaled<-TRUE # decide if abundance data should be scaled

# taxa settings (for bacteria and fungi)
if (parameter == "bacteria") {
  taxranks <- c("domain", "phylum", "class", "order", "family", "genus")
  taxadepth <- "class" # set depth of classification
  threshold <- 0.01 # OTU's less than ... % will be put into "Other ..." (limited legend display in plots)
} else if (parameter == "fungi") {
  taxranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  taxadepth <- "class" # set depth of classification
  threshold <- 0.01 # OTU's less than ... % will be put into "Other ..." (limited legend display in plots)
} else if (parameter %in% autotrophs) {
  taxadepth <- "morphotype"
}

#### END OF SETTINGS ###

#### LIBRARIES ####
# function that checks if a specific package is installed (if not, it will be installed), and loads the package afterwards
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need)
    lapply(need, require, character.only = TRUE)
  }
}

# packages needed
packages <- c(
  "AICcPermanova",
  "broom",
  "car",
  "cowplot",
  "data.table",
  "dplyr",
  "effectsize",
  "emmeans",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggsci",
  "grid",
  "gridExtra",
  "lattice",
  "lme4",
  "lmerTest",
  "lmPerm",
  "patternplot",
  "phyloseq",
  "png",
  "RColorBrewer",
  "scales",
  "stringr",
  "tidyverse",
  "vegan"
)

using(packages)

#### FUNCTIONS ####
paste3 <- function(..., sep = ", ") {
  L <- list(...)
  L <- lapply(L, function(x) {
    x[is.na(x)] <- ""
    x
  })
  ret <- gsub(
    paste0("(^", sep, "|", sep, "$)"), "",
    gsub(
      paste0(sep, sep), sep,
      do.call(paste, c(L, list(sep = sep)))
    )
  )
  is.na(ret) <- ret == ""
  ret
}

se <- function(x) sqrt(var(x, na.rm = T) / length(x))

give.n <- function(y) {
  return(data.frame(y = ifelse(quantile(y, 0.75) >= 0, quantile(y, 0.75) * 1.1, quantile(y, 0.75) * 0.9), label = paste("n =", length(y))))
}

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

annotation_compass <- function(label, position = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"), fontsize = 20, padding = grid::unit(c(0.5, 0.5), "line"), ...) {
  position <- match.arg(position)
  x <- switch(position,
              N = 0.5,
              NE = 1,
              E = 1,
              SE = 1,
              S = 0.5,
              SW = 0,
              W = 0,
              NW = 0
  )
  y <- switch(position,
              N = 1,
              NE = 1,
              E = 0.5,
              SE = 0,
              S = 0,
              SW = 0,
              W = 0.5,
              NW = 1
  )
  hjust <- switch(position,
                  N = 0.5,
                  NE = 1,
                  E = 1,
                  SE = 1,
                  S = 0.5,
                  SW = 0,
                  W = 0,
                  NW = 0
  )
  vjust <- switch(position,
                  N = 1,
                  NE = 1,
                  E = 0.5,
                  SE = 0,
                  S = 0,
                  SW = 0,
                  W = 0.5,
                  NW = 1
  )
  f1 <- switch(position,
               N = 0,
               NE = -1,
               E = -1,
               SE = -1,
               S = 0,
               SW = 1,
               W = 1,
               NW = 1
  )
  f2 <- switch(position,
               N = -1,
               NE = -1,
               E = 0,
               SE = 1,
               S = 1,
               SW = 1,
               W = 0,
               NW = -1
  )
  annotation_custom(grid::textGrob(label,
                                   gp = gpar(col = "black", fontsize = fontsize),
                                   x = grid::unit(x, "npc") + f1 * padding[1],
                                   y = grid::unit(y, "npc") + f2 * padding[2],
                                   hjust = hjust, vjust = vjust, ...
  ))
}

tax_zoom <- function(physeq, zoomtaxa = NULL, zoomlevel = NULL, taxdepth = "class", taxranks = rank_names(physeq), return.means = FALSE) {
  
  options(scipen=999) # disable scientific notation
  
  meta<-sample_data(physeq)
  taxonomy<-tax_table(physeq)
  allcols<-c("sample_ID", names(meta))
  factors<-names(meta)[-which(names(meta) %in% "replicate")]
  grouping.factors<-lapply(factors, as.symbol)
  
  taxranks <- taxranks[apply(taxonomy,2, function(x)!all(is.na(x)))]
  
  taxalevels_depth<-vector()
  for (i in seq_along(taxranks)) {
    if (taxranks[i]!=taxdepth) {
      taxalevels_depth[i]<-taxranks[i]
    } else if (taxranks[i]==taxdepth) {
      taxalevels_depth[i]<-taxdepth
      break
    }
  }
  
  subset_taxa2<-function(physeq, subset_level=zoom_level, subset_taxa=zoomtaxa) {
    
    physeq_otu<-otu_table(physeq)
    physeq_taxonomy<-tax_table(physeq)
    physeq_sample<-sample_data(physeq)
    
    new_taxonomy<-physeq_taxonomy %>%
      as.data.frame() %>%
      filter(!!as.name(subset_level)==subset_taxa) %>%
      as.matrix() %>%
      tax_table()
    
    sub_otu<-row.names(new_taxonomy)
    new_otu<-physeq_otu %>%
      as.data.frame() %>% 
      select(all_of(sub_otu)) %>%
      otu_table(taxa_are_rows=FALSE)
    
    new_physeq<-phyloseq(new_otu, new_taxonomy, physeq_sample)
    
    return(new_physeq)
  }
  
  if (!is.null(zoomtaxa) & !is.null(zoomlevel)) {
    physeq_zoom<-subset_taxa2(physeq, subset_level = zoomlevel, subset_taxa = zoomtaxa)
  } else {
    physeq_zoom<-physeq
  }
  
  # calculate relative abundances based on OTU's found
  otu_counts_zoom<-vector("list", length=length(taxalevels_depth))
  otu_rel_zoom<-vector("list", length=length(taxalevels_depth))
  relab_zoom<-vector("list", length=length(taxalevels_depth))
  absab_zoom<-vector("list", length=length(taxalevels_depth))
  names(otu_counts_zoom)<-taxalevels_depth
  names(otu_rel_zoom)<-taxalevels_depth
  names(relab_zoom)<-taxalevels_depth
  names(absab_zoom)<-taxalevels_depth
  
  for (i in seq_along(taxalevels_depth)) {
    otu_counts_zoom[[i]]<-tax_glom(physeq_zoom, taxrank=taxalevels_depth[i])
    otu_rel_zoom[[i]]<-transform_sample_counts(otu_counts_zoom[[i]], function(x) x/sum(x))
    
    abscounts_zoom<-as.data.frame(t(otu_table(otu_counts_zoom[[i]])))
    relcounts_zoom<-as.data.frame(t(otu_table(otu_rel_zoom[[i]])))
    
    taxnames_zoom<-tax_table(otu_rel_zoom[[i]])
    taxnames_zoom<-taxnames_zoom[which(row.names(taxnames_zoom) %in% row.names(relcounts_zoom)), ]
    taxnames_zoom<-t(as.data.frame(t(taxnames_zoom)))
    
    relab_zoom[[i]]<-merge(taxnames_zoom, relcounts_zoom, by="row.names")
    #names(relab_zoom[[i]])[1]<-"ASV_ID"
    relab_zoom[[i]]<-relab_zoom[[i]][, -1]
    
    absab_zoom[[i]]<-merge(taxnames_zoom, abscounts_zoom, by="row.names")
    #names(absab_zoom[[i]])[1]<-"ASV_ID"
    absab_zoom[[i]]<-absab_zoom[[i]][, -1]
  }
  
  relab_zoom<-lapply(relab_zoom, function(x) x[do.call(order, x[, taxalevels_depth]), ]) # order alphabetically after tax ranks
  relab_zoom<-lapply(relab_zoom, function(x) x[, colSums(is.na(x)) != nrow(x)]) # remove taxrank columns that are not needed (= entirely NA)
  absab_zoom<-lapply(absab_zoom, function(x) x[do.call(order, x[, taxalevels_depth]), ]) # order alphabetically after tax ranks
  absab_zoom<-lapply(absab_zoom, function(x) x[, colSums(is.na(x)) != nrow(x)]) # remove taxrank columns that are not needed (= entirely NA)
  
  relab_long<-vector("list", length=length(taxalevels_depth))
  names(relab_long)<-taxalevels_depth
  absab_long<-vector("list", length=length(taxalevels_depth))
  names(absab_long)<-taxalevels_depth
  
  relab_single<-relab_long
  relab_single_big<-relab_single
  relab_single_unclass<-relab_single
  relab_single_other<-relab_single
  relab_cleaned_single<-relab_single
  relab_cleaned_mean<-relab_single
  
  # relab_mean<-relab_long
  # relab_mean_big<-relab_mean
  # relab_mean_unclass<-relab_mean
  # relab_mean_other<-relab_mean
  # relab_cleaned<-relab_mean
  
  for (i in seq_along(taxalevels_depth)) {
    relab_long[[i]]<-reshape2::melt(relab_zoom[[i]], id.vars=taxalevels_depth[1:i], variable.name = "sample_ID", value.name = "rel_abundance")
    relab_long[[i]]<-merge(cbind(sample_ID=row.names(meta), meta), relab_long[[i]], by="sample_ID")
    
    absab_long[[i]]<-reshape2::melt(absab_zoom[[i]], id.vars=taxalevels_depth[1:i], variable.name = "sample_ID", value.name = "abs_abundance")
    absab_long[[i]]<-merge(cbind(sample_ID=row.names(meta), meta), absab_long[[i]], by="sample_ID")
    
    dots<-lapply(taxalevels_depth[1:i], as.symbol)
    dots2<-lapply(taxalevels_depth[1:i-1], as.symbol)
    
    relab_long[[i]]<-relab_long[[i]] %>%
      left_join(absab_long[[i]], by=(c(allcols, taxalevels_depth[1:i])))
    # ## option 1: apply threshold individually (taxa A will only be displayed in samples where it is above threshold, in all others it will be put into the "Other" category)
    # relab_mean[[i]]<-relab_long[[i]] %>%
    #   group_by(!!as.name(var1), .dots=dots) %>%
    #   summarize(mean_rel_ab=mean(rel_abundance)) %>%
    #   mutate(flag = ifelse(mean_rel_ab >= threshold, 1, 0)) %>%
    #   ungroup() %>%
    #   mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
    #   as.data.frame()
    # ##
    
    ## option 2: apply threshold to all (if taxa A in one sample is above threshold, it will be displayed for all other samples, too, although there it might be below threshold)
    # relab_mean[[i]]<-relab_long[[i]] %>%
    #   group_by(.dots=c(grouping.factors, dots)) %>%
    #   summarize(mean_rel_ab=mean(rel_abundance), mean_abs_ab=mean(abs_abundance)) %>%
    #   ungroup() %>%
    #   group_by(.dots=dots) %>%
    #   mutate(flag = ifelse(max(mean_rel_ab) >= threshold, 1, 0)) %>%
    #   ungroup() %>%
    #   mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
    #   as.data.frame()
    
    relab_single[[i]]<-relab_long[[i]] %>%
      group_by(.dots=c(grouping.factors, dots)) %>%
      mutate(mean_rel_ab=mean(rel_abundance)) %>%
      ungroup() %>%
      group_by(.dots=dots) %>%
      mutate(flag = ifelse(max(mean_rel_ab) >= threshold, 1, 0)) %>%
      select(-mean_rel_ab) %>%
      ungroup() %>%
      mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
      as.data.frame()
    ##
    
    # relab_mean_big[[i]]<-relab_mean[[i]] %>% # only taxa above threshold, plus no unclassified taxa
    #   filter(flag==1) %>%
    #   #filter(rowSums(. == "unclassified")==0) %>%
    #   filter(apply(., 1, function(r) !any(grepl("*unclassified*", r)))) %>%
    #   select(-flag) %>%
    #   droplevels() %>%
    #   as.data.frame()
    
    relab_single_big[[i]]<-relab_single[[i]] %>% # only taxa above threshold, plus no unclassified taxa
      filter(flag==1) %>%
      #filter(rowSums(. == "unclassified")==0) %>%
      filter(apply(., 1, function(r) !any(grepl("*unclassified*", r)))) %>%
      select(-flag) %>%
      droplevels() %>%
      as.data.frame()
    
    # relab_mean_unclass[[i]]<-relab_mean[[i]] %>% # only unclassified taxa
    #   select(-flag) %>%
    #   #filter(apply(., 1, function(r) any(r %in% "unclassified"))) 
    #   filter(apply(., 1, function(r) any(grepl("*unclassified*", r))))
    
    relab_single_unclass[[i]]<-relab_single[[i]] %>% # only unclassified taxa
      select(-flag) %>%
      #filter(apply(., 1, function(r) any(r %in% "unclassified"))) 
      filter(apply(., 1, function(r) any(grepl("*unclassified*", r))))
    
    # if (nrow(relab_mean_unclass[[i]])!=0) {
    #   relab_mean_unclass[[i]]<-relab_mean_unclass[[i]] %>%
    #     mutate(!!(taxalevels_depth[i]):="Unclassified") %>%
    #     mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
    #     #group_by(.dots=c(grouping.factors, zoomlevel)) %>%
    #     #summarize(sums=sum(mean_rel_ab)) %>%
    #     droplevels() %>%
    #     as.data.frame()
    # }
    if (nrow(relab_single_unclass[[i]])!=0) {
      relab_single_unclass[[i]]<-relab_single_unclass[[i]] %>%
        mutate(!!(taxalevels_depth[i]):="Unclassified") %>%
        mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
        #group_by(.dots=c(grouping.factors, zoomlevel)) %>%
        #summarize(sums=sum(mean_rel_ab)) %>%
        droplevels() %>%
        as.data.frame()
    }
    
    if (i == 1) {
      # relab_mean_other[[i]]<-relab_mean[[i]] %>% # only taxa below threshold, summed to "Other"
      #   filter(flag==0) %>%
      #   #filter(rowSums(. == "unclassified")==0) %>%
      #   filter(apply(., 1, function(r) !any(grepl("*unclassified*", r)))) %>%
      #   group_by(.dots=grouping.factors) %>%
      #   summarize("mean_rel_ab"=sum(mean_rel_ab)) %>%
      #   mutate(!!(taxalevels_depth[i]):="Other") %>%
      #   mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
      #   droplevels() %>%
      #   as.data.frame()
      relab_single_other[[i]]<-relab_single[[i]] %>% # only taxa below threshold, summed to "Other"
        filter(flag==0) %>%
        #filter(rowSums(. == "unclassified")==0) %>%
        filter(apply(., 1, function(r) !any(grepl("*unclassified*", r)))) %>%
        group_by(.dots=c(grouping.factors, "replicate")) %>%
        summarize("sum_rel_ab"=sum(rel_abundance), "sum_abs_ab"=sum(abs_abundance)) %>%
        mutate(!!(taxalevels_depth[i]):="Other") %>%
        mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
        droplevels() %>%
        as.data.frame()
    } else {
      # relab_mean_other[[i]]<-relab_mean[[i]] %>% # only taxa below threshold, summed to "Other"
      #   filter(flag==0) %>%
      #   #filter(rowSums(. == "unclassified")==0) %>%
      #   filter(apply(., 1, function(r) !any(grepl("*unclassified*", r)))) %>%
      #   group_by(.dots=c(grouping.factors, dots2)) %>%
      #   summarize("mean_rel_ab"=sum(mean_rel_ab)) %>%
      #   mutate(!!(taxalevels_depth[i]):="Other") %>%
      #   mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
      #   droplevels() %>%
      #   as.data.frame()
      relab_single_other[[i]]<-relab_single[[i]] %>% # only taxa below threshold, summed to "Other"
        filter(flag==0) %>%
        #filter(rowSums(. == "unclassified")==0) %>%
        filter(apply(., 1, function(r) !any(grepl("*unclassified*", r)))) %>%
        group_by(.dots=c(grouping.factors, "replicate", dots2)) %>%
        summarize("sum_rel_ab"=sum(rel_abundance), "sum_abs_ab"=sum(abs_abundance)) %>%
        mutate(!!(taxalevels_depth[i]):="Other") %>%
        mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
        rename(rel_abundance=sum_rel_ab, abs_abundance=sum_abs_ab) %>%
        droplevels() %>%
        as.data.frame()
    }
    
    # combine to new single datafram with relative abundances
    #relab_cleaned[[i]]<-bind_rows(relab_mean_big[[i]], relab_mean_unclass[[i]], relab_mean_other[[i]])
    
    relab_cleaned_single[[i]]<-bind_rows(relab_single_big[[i]], relab_single_unclass[[i]], relab_single_other[[i]])
    
    relab_cleaned_mean[[i]]<-relab_cleaned_single[[i]] %>% 
      group_by(.dots=c(grouping.factors, dots)) %>% 
      summarize(mean_rel_ab=mean(rel_abundance), mean_abs_ab=mean(abs_abundance)) %>%
      as.data.frame()
  }
  
  options(scipen=0) # back to scientific notation
  abundance_means<-relab_cleaned_mean[[i]] %>% droplevels()
  abundance_single<-relab_cleaned_single[[i]] %>% droplevels()
  
  abundance_taxa_mean<-abundance_means %>%
    rename(absolute_mean=mean_abs_ab) %>%
    mutate(percent_mean=mean_rel_ab*100) %>%
    select(-mean_rel_ab) %>%
    as.data.frame()
  
  abundance_taxa_single<-abundance_single %>%
    rename(absolute_abundance=abs_abundance) %>%
    mutate(percent_abundance=rel_abundance*100) %>%
    select(-rel_abundance) %>%
    as.data.frame()
  
  # taxonomy_selected<-taxonomy %>%
  #   as.data.frame() %>%
  #   group_by(.dots=dots) %>%
  #   summarize()
  # 
  # abundance_taxa_compact<-abundance_taxa %>%
  #   group_by(.dots=grouping.factors, class) %>%
  #   summarize(sum_percent_mean=sum(percent_mean))  %>%
  #   mutate("percent_mean"=sum_percent_mean) %>%
  #   select(-sum_percent_mean) %>%
  #   arrange(!!as.name(factors), !!!taxdepth) %>%
  #   left_join(taxonomy_selected, by=taxdepth) %>%
  #   select(!!as.name(factors), !!!taxalevels_depth, percent_mean) %>%
  #   as.data.frame()
  # 
  # for (i in rev(seq_along(taxalevels_depth))) {
  #   abundance_taxa_compact<-abundance_taxa_compact %>%
  #     mutate(!!as.name(taxalevels_depth[i]):=ifelse(is.na(!!as.name(taxalevels_depth[i])), !!as.name(taxdepth), !!as.name(taxalevels_depth[i]))) %>%
  #     as.data.frame()
  # }
  
  #abundance_taxa<-abundance_taxa_compact
  if (return.means==TRUE) {
    abundance_taxa<-abundance_taxa_mean
  } else {
    abundance_taxa<-abundance_taxa_single
  }
  
  abundance_taxa<-abundance_taxa %>%
    mutate_if(is.factor, as.character) %>%
    replace(., sapply(abundance_taxa, function(.) grepl("*unclassified*", .)), "Unclassified") %>%
    mutate_if(is.character, as.factor) %>%
    as.data.frame()
  
  return(abundance_taxa)
}

recode <- dplyr::recode
select <- dplyr::select
filter <- dplyr::filter

#### DATA IMPORT ####
setwd(dir_in)

## METADATA
meta <- read.table(metafile, sep = "\t", header = T, row.names = 1)
# rename factor levels
meta[, which(names(meta) %in% var1)] <- dplyr::recode(meta[, which(names(meta) %in% var1)], !!!levels_var1)
meta[, which(names(meta) %in% var1)] <- factor(meta[, which(names(meta) %in% var1)], levels = levels_var1)
meta[, which(names(meta) %in% var2)] <- recode(meta[, which(names(meta) %in% var2)], !!!levels_var2)
meta[, which(names(meta) %in% var2)] <- factor(meta[, which(names(meta) %in% var2)], levels = levels_var2)
# choose only selected rows
meta <- meta[which(meta[, which(names(meta) %in% var1)] %in% var1_selected), ]
meta <- meta[which(meta[, which(names(meta) %in% var2)] %in% var2_selected), ]
meta <- meta %>%
  mutate(logclass = ifelse(transport == "Initial", "Initial", "Final")) %>%
  as.data.frame()

# meta<-meta$sediment[which(meta$logclass=="Initial")]
meta <- droplevels(meta)
meta <- meta[order(as.numeric(row.names(meta))), ] # rows in numerical order

meta_all <- rownames_to_column(meta, var = "id")

# exclude Terrestrial from stats
rn <- meta_all$id[meta_all$sediment == "Terrestrial"]
meta_stats <- meta_all %>%
  filter(!(id %in% rn)) %>%
  droplevels() %>%
  as.data.frame()
row.names(meta_stats) <- meta_stats$id

if (include.initial == FALSE) {
  rn <- meta_all$id[meta_all$transport == "Initial"]
  meta_all <- meta_all %>%
    filter(!(id %in% rn)) %>%
    droplevels() %>%
    as.data.frame()
  row.names(meta_all) <- meta_all$id
  
  rn <- meta_stats$id[meta_stats$transport == "Initial"]
  meta_stats <- meta_stats %>%
    filter(!(id %in% rn)) %>%
    droplevels() %>%
    as.data.frame()
  row.names(meta_stats) <- meta_stats$id
}

## OTU DATA (or algae morphotypes)

otu <- read.table(otufile, sep = "\t", dec = ",", header = T, row.names = 1, check.names = F) # OTU table

# filter out algae that are not diatoms
# if (parameter == "diatoms") {
#   not_diatoms <- c("Merismopedia", "Pseudoanabaena", "Scenedesmus")
#   otu <- otu %>% 
#     rownames_to_column() %>%
#     filter(!(rowname %in% not_diatoms)) %>% 
#     column_to_rownames("rowname")
# }

# all OTUs
otu_all <- otu[, names(otu) %in% row.names(meta_all)]
otu_all <- otu_all[, order(as.numeric(names(otu_all)))] # columns in numerical order

# selected OTUs (with/without terrestrial/initial)
otu_stats <- otu[, names(otu) %in% row.names(meta_stats)]
otu_stats <- otu_stats[, order(as.numeric(names(otu_stats)))] # columns in numerical order

# trim meta in case of missing samples in otu file
meta_all <- meta_all[row.names(meta_all) %in% names(otu_all), ]
meta_stats <- meta_stats[row.names(meta_stats) %in% names(otu_stats), ]

## TAXONOMY DATA
if (parameter %in% heterotrophs) {
  taxonomy <- read.table(taxonomyfile, sep = "\t", header = T, row.names = 1, check.names = F) # Taxa table
  # names(taxonomy)[which(names(taxonomy) %in% "domain")]<-"kingdom" # for equal data handling bacteria/fungi
} else if (parameter %in% autotrophs) {
  taxonomy <- read.table(taxonomyfile, sep = "\t", header = T, check.names = F) # Taxa table

  taxonomy$group <- recode(taxonomy$group, !!!algae_levels)
  taxonomy$group <- factor(taxonomy$group, levels = algae_levels)

  taxonomy_all <- taxonomy

  if (parameter != "algae") {
    taxonomy <- taxonomy[taxonomy$group == parameter, ]
    
    otu_all$taxa <- row.names(otu_all)
    otu_all <- taxonomy %>%
      left_join(otu_all, by = "taxa") %>%
      as.data.frame()
    rownames(otu_all) <- otu_all$taxa
    otu_all <- otu_all %>%
      select(-group, -taxa) %>%
      as.data.frame()
    
    otu_stats$taxa <- row.names(otu_stats)
    otu_stats <- taxonomy %>%
      left_join(otu_stats, by = "taxa") %>%
      as.data.frame()
    rownames(otu_stats) <- otu_stats$taxa
    otu_stats <- otu_stats %>%
      select(-group, -taxa) %>%
      as.data.frame()
  }
}

## ABUNDANCE
if (parameter %in% heterotrophs) {
  abun <- read.table(abundancefile, sep = "\t", header = T, check.names = F)
  # rename factor levels
  abun[, which(names(abun) %in% var1)] <- recode(abun[, which(names(abun) %in% var1)], !!!levels_var1)
  abun[, which(names(abun) %in% var1)] <- factor(abun[, which(names(abun) %in% var1)], levels = levels_var1)
  abun[, which(names(abun) %in% var2)] <- recode(abun[, which(names(abun) %in% var2)], !!!levels_var2)
  abun[, which(names(abun) %in% var2)] <- factor(abun[, which(names(abun) %in% var2)], levels = levels_var2)

  abun <- abun %>%
    mutate(logclass = ifelse(transport == "Initial", "Initial", "Final")) %>%
    as.data.frame()

  abun_all <- rownames_to_column(abun, var = "id")
  
  # exclude terrestrial from stats
  rn <- abun_all$id[abun_all$sediment == "Terrestrial"]
  abun_stats <- abun_all %>%
    filter(!(id %in% rn)) %>%
    droplevels() %>%
    as.data.frame()
  row.names(abun_stats) <- abun_stats$id

  if (include.initial == FALSE) {
    rn <- abun_all$id[abun_all$transport == "Initial"]
    abun_all <- abun_all %>%
      filter(!(id %in% rn)) %>%
      droplevels() %>%
      as.data.frame()
    row.names(abun_all) <- abun_all$id
    
    rn <- abun_stats$id[abun_stats$transport == "Initial"]
    abun_stats <- abun_stats %>%
      filter(!(id %in% rn)) %>%
      droplevels() %>%
      as.data.frame()
    row.names(abun_stats) <- abun_stats$id
  }
}

## ENVIRONMENTAL DATA
envir <- read.table(envirfile, header = T, dec = ".", sep = "\t")

# select only samples that have OTU data (i.e. match with meta file)
envir_all <- envir %>%
  right_join(meta_all, c("transport", "sediment", "replicate")) %>%
  relocate(id) %>%
  arrange(as.numeric(id)) %>% 
  select(-GPP) # deselect GPP (focus only on NCP)

#### DATA PREPARATION ####
# create phyloseq objects
# meta
meta_s <- sample_data(meta_all)
sample_names(meta_s) <- row.names(meta_all)

meta_s_stats <- sample_data(meta_stats)
sample_names(meta_s_stats) <- row.names(meta_stats)

# OTU
otu_t <- as.data.frame(t(otu_all))
OTU <- otu_table(otu_t, taxa_are_rows = F)

otu_t_stats <- as.data.frame(t(otu_stats))
OTU_stats <- otu_table(otu_t_stats, taxa_are_rows = F)

# find "empty" OTU's and remove from file
empty_otu <- colSums(otu_t)[colSums(otu_t) == 0]
if (length(empty_otu) > 0) {
  otu_t <- otu_t[-which(names(otu_t) %in% names(empty_otu))]
}
empty_otu <- colSums(otu_t_stats)[colSums(otu_t_stats) == 0]
if (length(empty_otu) > 0) {
  otu_t_stats <- otu_t_stats[-which(names(otu_t_stats) %in% names(empty_otu))]
}

# OTU's found (without rarefaction)
if (all(row.names(otu_t) == row.names(meta_all)) == TRUE) {
  meta2 <- cbind(meta_all, "OTUs" = rowSums(otu_t > 0), "individuals" = rowSums(otu_t))
} else {
  warning("Data.frames not merged: Row names of otu_t and meta_all do not match!")
}
if (all(row.names(otu_t_stats) == row.names(meta_stats)) == TRUE) {
  meta2_stats <- cbind(meta_stats, "OTUs" = rowSums(otu_t_stats > 0), "individuals" = rowSums(otu_t_stats))
} else {
  warning("Data.frames not merged: Row names of otu_t_stats and meta_stats do not match!")
}

## TAXONOMY DATA
if (parameter %in% heterotrophs) {
  taxonomy <- taxonomy[do.call("order", taxonomy[, which(names(taxonomy) %in% taxranks)]), ]
  taxonomy <- taxonomy[apply(taxonomy, 1, function(x) !any(grepl("Chloroplast", x))), ] # remove Chloroplasts (optional)
  taxonomy <- taxonomy[apply(taxonomy, 1, function(x) !any(grepl("Cyanobacteria", x))), ] # remove Cyanobacteria (optional)
  # create phyloseq object
  taxa <- tax_table(as.matrix(taxonomy))
} else if (parameter %in% autotrophs) {
  taxa <- tax_table(as.matrix(taxonomy))
}

## ENVIRONMENTAL DATA
# rename factor levels
envir_all[, which(names(envir_all) %in% var1)] <- recode(envir_all[, which(names(envir_all) %in% var1)], !!!levels_var1)
envir_all[, which(names(envir_all) %in% var1)] <- factor(envir_all[, which(names(envir_all) %in% var1)], levels = levels_var1)
envir_all[, which(names(envir_all) %in% var2)] <- recode(envir_all[, which(names(envir_all) %in% var2)], !!!levels_var2)
envir_all[, which(names(envir_all) %in% var2)] <- factor(envir_all[, which(names(envir_all) %in% var2)], levels = levels_var2)

# choose only selected rows
envir_all <- envir_all[which(envir_all[, which(names(envir_all) %in% var1)] %in% var1_selected), ]
envir_all <- envir_all[which(envir_all[, which(names(envir_all) %in% var2)] %in% var2_selected), ]
envir_all <- droplevels(envir_all)
envir_all <- envir_all %>%
  mutate(logclass = ifelse(transport == "Initial", "Initial", "Final")) %>%
  relocate(id, transport, sediment, replicate, logclass)

rn <- envir_all$id[envir_all$sediment == "Terrestrial"]
envir_stats <- envir_all %>%
  filter(!(id %in% rn)) %>%
  droplevels() %>%
  as.data.frame()
row.names(envir_stats) <- envir_stats$id

if (include.initial == FALSE) {
  rn <- envir_all$id[envir_all$transport == "Initial"]
  envir_all <- envir_all %>%
    filter(!(id %in% rn)) %>%
    droplevels() %>%
    as.data.frame()
  row.names(envir_all) <- envir_all$id
  
  rn <- envir_stats$id[envir_stats$transport == "Initial"]
  envir_stats <- envir_stats %>%
    filter(!(id %in% rn)) %>%
    droplevels() %>%
    as.data.frame()
  row.names(envir_stats) <- envir_stats$id
}

# get average initial and final metabolism rates
cols <- names(envir_all)[-c(1:5)]
envir_init <- envir_all %>%
  filter(logclass == "Initial") %>%
  group_by(logclass, transport, sediment) %>%
  summarize(across(all_of(cols), ~ mean(.x, na.rm = TRUE))) %>%
  select(-logclass) %>%
  mutate(across(all_of(cols), ~ ifelse(is.nan(.), NA, .))) %>%
  as.data.frame()

envir_final <- envir_all %>%
  filter(logclass == "Final") %>%
  group_by(logclass, transport, sediment) %>%
  summarize(across(all_of(cols), ~ mean(.x, na.rm = TRUE))) %>%
  select(-logclass) %>%
  mutate(across(all_of(cols), ~ ifelse(is.nan(.), NA, .))) %>%
  as.data.frame()

envir_means <- envir_final %>%
  rbind.data.frame(envir_init)

# create envir file to be fitted to NMDS
envirfit_init <- envir_all %>%
  filter(logclass == "Initial")

envirfit_final <- envir_all %>%
  filter(logclass == "Final")

## MERGE META, OTU & TAXA DATA to phyloseq object
if (parameter %in% heterotrophs) {
  seqdata <- phyloseq(OTU, taxa, meta_s)
  seqdata_stats <- phyloseq(OTU_stats, taxa, meta_s_stats)
} else if (parameter %in% autotrophs) {
  seqdata <- phyloseq(OTU, meta_s)
  seqdata_stats <- phyloseq(OTU_stats, meta_s_stats)
}

#### DATA ANALYSIS ####
#### (1) ABSOLUTE ABUNDANCE ####
if (parameter %in% heterotrophs) {
  abundance_group <- abun_all %>%
    pivot_longer(!c(id, transport, sediment, replicate, logclass), names_to = "variable", values_to = "value") %>%
    as.data.frame()
} else if (parameter %in% autotrophs) {
  otu_all$taxa <- rownames(otu_all)
  taxotu <- taxonomy %>%
    left_join(otu_all, by = "taxa") %>%
    pivot_longer(!c(group, taxa), names_to = "id", values_to = "value") %>%
    as.data.frame()

  abun <- meta_all %>%
    left_join(taxotu, by = "id") %>%
    as.data.frame()
  
  abundance_group <- abun %>%
    group_by(id, transport, sediment, replicate, group, logclass) %>%
    summarise(total = sum(value)) %>%
    ungroup() %>%
    mutate(variable = group, value = total) %>%
    select(-group, -total) %>%
    as.data.frame()
}

# calculate total abundance of organisms
abundance_tot <- abundance_group %>%
  group_by(id, transport, sediment, replicate, logclass) %>%
  summarise(total = sum(value)) %>%
  ungroup() %>%
  mutate(variable = "total", value = total) %>%
  select(-total) %>%
  as.data.frame()

abundance_group <- abundance_group %>%
  rbind.data.frame(abundance_tot) %>%
  arrange(transport, sediment, replicate, logclass, variable) %>%
  group_by(id, sediment, transport, replicate, logclass) %>%
  mutate(percent = ifelse(variable != "total", value[variable != "total"] / sum(value[variable == "total"]) * 100, 100)) %>%
  as.data.frame()

abundance_group_long <- abundance_group %>%
  pivot_longer(!c(id, transport, sediment, replicate, logclass, variable), names_to = "var2", values_to = "value") %>%
  mutate(var1 = paste(variable, var2, sep = "_")) %>%
  select(-variable, -var2) %>%
  mutate(variable = var1) %>%
  select(-var1) %>%
  select(id, transport, sediment, replicate, logclass, variable, value) %>%
  as.data.frame()

# mean abundance of groups (heterotrophs: bacteria/fungi, autotrophs: cyanobacteria/diatoms/green algae)
abundance_group_means <- abundance_group %>%
  group_by(transport, sediment, variable, logclass) %>%
  summarize(mean_value = mean(value, na.rm = T), sd_value = sd(value, na.rm = T)) %>%
  ungroup() %>%
  group_by(sediment, transport, logclass) %>%
  mutate(percent_mean = ifelse(variable != "total", mean_value[variable != "total"] / sum(mean_value[variable != "total"]) * 100, 100)) %>%
  mutate(percent_sd = sd_value / sum(mean_value) * 100) %>%
  as.data.frame()

abundance_group_means_long <- abundance_group_means %>%
  pivot_longer(!c(transport, sediment, variable, logclass), names_to = "var2", values_to = "values") %>%
  mutate(var2 = gsub("_value", "", var2)) %>%
  mutate(var2 = gsub("_mean", "", var2)) %>%
  filter(var2 != "sd", var2 != "percent_sd") %>%
  mutate(var1 = paste(variable, var2, sep = "_")) %>%
  select(-variable, -var2) %>%
  mutate(variable = var1) %>%
  mutate(values = format(values, scientific = FALSE)) %>%
  select(-var1) %>%
  select(transport, sediment, logclass, variable, values) %>%
  as.data.frame()

#### (2) RELATIVE ABUNDANCE ####
# File preparation (like pivot table)
if (parameter %in% heterotrophs) {
  options(scipen = 999) # disable scientific notation

  abundance_taxa <- tax_zoom(seqdata, taxdepth = "class", zoomtaxa = NULL, zoomlevel = NULL, return.means = TRUE)
  abundance_taxa_single <- tax_zoom(seqdata, taxdepth = "class", zoomtaxa = NULL, zoomlevel = NULL, return.means = FALSE)
  if (parameter == "bacteria") {
    # zoom in specific taxa
    taxa_zoom <- "Proteobacteria"
    level_zoom <- "phylum"
    abundance_taxa_zoom <- tax_zoom(seqdata, zoomtaxa = taxa_zoom, zoomlevel = level_zoom, taxdepth = "order", return.means = TRUE)
    abundance_taxa_zoom_single <- tax_zoom(seqdata, zoomtaxa = taxa_zoom, zoomlevel = level_zoom, taxdepth = "order", return.means = FALSE)
  } else {
    abundance_taxa_zoom <- NULL
    abundance_taxa_zoom_single <- NULL
  }
   
  abundance_taxa_single %>%
    #filter(domain != "Archaea") %>%
    filter(transport == "Initial", sediment == "Terrestrial") %>%
    group_by(transport, sediment, phylum, class) %>%
    summarize(n(), mean_percent_abundance = mean(percent_abundance, na.rm = TRUE), sd_percent_abundance = sd(percent_abundance, na.rm = TRUE)) %>%
    arrange(-mean_percent_abundance) %>%
    View()
  
  # # OLD CODE
  # taxalevels_depth <- vector()
  # for (i in seq_along(taxranks)) {
  #   if (taxranks[i] != taxadepth) {
  #     taxalevels_depth[i] <- taxranks[i]
  #   } else if (taxranks[i] == taxadepth) {
  #     taxalevels_depth[i] <- taxadepth
  #     break
  #   }
  # }
  # 
  # # calculate relative abundances based on OTU's found
  # otu_counts <- vector("list", length = length(taxalevels_depth))
  # otu_rel <- vector("list", length = length(taxalevels_depth))
  # relab <- vector("list", length = length(taxalevels_depth))
  # names(otu_counts) <- taxalevels_depth
  # names(otu_rel) <- taxalevels_depth
  # names(relab) <- taxalevels_depth
  # for (i in seq_along(taxalevels_depth)) {
  #   otu_counts[[i]] <- tax_glom(seqdata, taxrank = taxalevels_depth[i])
  #   otu_rel[[i]] <- transform_sample_counts(otu_counts[[i]], function(x) x / sum(x))
  #   relcounts <- as.data.frame(t(otu_table(otu_rel[[i]])))
  #   taxnames <- tax_table(otu_rel[[i]])
  #   taxnames <- taxnames[which(row.names(taxnames) %in% row.names(relcounts)), ]
  #   taxnames <- t(as.data.frame(t(taxnames)))
  #   relab[[i]] <- merge(taxnames, relcounts, by = "row.names")
  #   relab[[i]] <- relab[[i]][, -1]
  # }
  # 
  # relab <- lapply(relab, function(x) x[do.call(order, x[, taxalevels_depth]), ]) # order alphabetically after tax ranks
  # relab <- lapply(relab, function(x) x[, colSums(is.na(x)) != nrow(x)]) # remove taxrank columns that are not needed (= entirely NA)
  # 
  # relab_long <- vector("list", length = length(taxalevels_depth))
  # names(relab_long) <- taxalevels_depth
  # relab_mean <- relab_long
  # relab_mean_big <- relab_mean
  # relab_mean_unclass <- relab_mean
  # relab_mean_other <- relab_mean
  # relab_cleaned <- relab_mean
  # for (i in seq_along(taxalevels_depth)) {
  #   relab_long[[i]] <- reshape2::melt(relab[[i]], id.vars = taxalevels_depth[1:i], variable.name = "sample_ID", value.name = "rel_abundance")
  #   relab_long[[i]] <- merge(cbind(sample_ID = row.names(meta), meta), relab_long[[i]], by = "sample_ID")
  #   dots <- lapply(taxalevels_depth[1:i], as.symbol)
  #   dots2 <- lapply(taxalevels_depth[1:i - 1], as.symbol)
  # 
  #   # ## option 1: apply threshold individually (taxa A will only be displayed in samples where it is above threshold, in all others it will be put into the "Other" category)
  #   # relab_mean[[i]]<-relab_long[[i]] %>%
  #   #   group_by(!!as.name(var1), !!as.name(var2), .dots=dots) %>%
  #   #   summarize(mean_rel_ab=mean(rel_abundance),
  #   #             sd_rel_ab = sd(rel_abundance)) %>%
  #   #   mutate(flag = ifelse(mean_rel_ab >= threshold, 1, 0)) %>%
  #   #   ungroup() %>%
  #   #   mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
  #   #   as.data.frame()
  #   # ##
  # 
  #   ## option 2: apply threshold to all (if taxa A in one sample is above threshold, it will be displayed for all other samples, too, although there it might be below threshold)
  #   relab_mean[[i]] <- relab_long[[i]] %>%
  #     group_by(!!as.name(var1), !!as.name(var2), .dots = dots) %>%
  #     summarize(mean_rel_ab = mean(rel_abundance),
  #               sd_rel_ab = sd(rel_abundance)) %>%
  #     ungroup() %>%
  #     group_by(.dots = dots) %>%
  #     mutate(flag = ifelse(max(mean_rel_ab) >= threshold, 1, 0)) %>%
  #     ungroup() %>%
  #     mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
  #     as.data.frame()
  #   ##
  # 
  #   relab_mean_big[[i]] <- relab_mean[[i]] %>% # only taxa above threshold, plus no unclassified taxa
  #     filter(flag == 1) %>%
  #     filter(!!as.name(taxalevels_depth[i]) != "unclassified") %>%
  #     #filter(rowSums(. == "unclassified") == 0) %>%
  #     select(-flag) %>%
  #     droplevels() %>%
  #     as.data.frame()
  # 
  #   relab_mean_unclass[[i]] <- relab_mean[[i]] %>% # only unclassified taxa
  #     select(-flag) %>%
  #     filter(apply(., 1, function(r) any(r %in% "unclassified")))
  # 
  #   if (nrow(relab_mean_unclass[[i]]) != 0) {
  #     relab_mean_unclass[[i]] <- relab_mean_unclass[[i]] %>%
  #       mutate(!!(taxalevels_depth[i]) := "Unclassified") %>%
  #       mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
  #       droplevels() %>%
  #       as.data.frame()
  #   }
  # 
  #   if (i == 1) {
  #     relab_mean_other[[i]] <- relab_mean[[i]] %>% # only taxa below threshold, summed to "Other"
  #       filter(flag == 0) %>%
  #       filter(!!as.name(taxalevels_depth[i]) != "unclassified") %>%
  #       #filter(rowSums(. == "unclassified") == 0) %>%
  #       group_by(!!as.name(var1), !!as.name(var2)) %>%
  #       summarize("mean_rel_ab" = sum(mean_rel_ab),
  #                 "sd_rel_ab" = NA) %>%
  #       mutate(!!(taxalevels_depth[i]) := "Other") %>%
  #       mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
  #       droplevels() %>%
  #       as.data.frame()
  #   } else {
  #     relab_mean_other[[i]] <- relab_mean[[i]] %>% # only taxa below threshold, summed to "Other"
  #       filter(flag == 0) %>%
  #       filter(!!as.name(taxalevels_depth[i]) != "unclassified") %>%
  #       #filter(rowSums(. == "unclassified") == 0) %>%
  #       group_by(!!as.name(var1), !!as.name(var2), .dots = dots2) %>%
  #       summarize("mean_rel_ab" = sum(mean_rel_ab),
  #                 "sd_rel_ab" = NA) %>%
  #       mutate(!!(taxalevels_depth[i]) := "Other") %>%
  #       mutate_at(vars(!!(taxalevels_depth[1:i])), list(factor)) %>% # convert to factor
  #       droplevels() %>%
  #       as.data.frame()
  #   }
  # 
  #   # combine to new single datafram with relative abundances
  #   # levels(relab_mean_big[[i]][[taxranks[1]]])<-c(levels(relab_mean_big[[i]][[taxranks[1]]]), "Other", "Unclassified")
  # 
  #   relab_cleaned[[i]] <- bind_rows(relab_mean_big[[i]], relab_mean_unclass[[i]], relab_mean_other[[i]])
  #   # relab_cleaned[[i]][[taxranks[1]]]<-as.factor(relab_cleaned[[i]][[taxranks[1]]]) # needed?
  # 
  #   # # not needed anymore?
  #   # if (i > 1) {
  #   #   for (j in 2:i) {
  #   #     k<-taxalevels_depth[j]
  #   #     l<-taxalevels_depth[j-1]
  #   #     levels(relab_cleaned[[i]][[k]])<-unique(c(levels(relab_cleaned[[i]][[k]]), "Other", "Unclassified"))
  #   #     relab_cleaned[[i]][[k]][which(is.na(relab_cleaned[[i]][k]))]<-relab_cleaned[[i]][[l]][which(is.na(relab_cleaned[[i]][k]))]
  #   #   }
  #   # }
  # }
  # 
  # options(scipen = 0) # back to scientific notation
  # abundance_means <- relab_cleaned[[i]]
  # abundance_means <- droplevels(abundance_means)
  # names(abundance_means)[length(names(abundance_means))-1] <- "mean_value"
  # names(abundance_means)[length(names(abundance_means))] <- "sd_value"
  # 
  # abundance_taxa <- abundance_means %>%
  #   mutate(percent_mean = mean_value * 100,
  #          percent_sd = sd_value * 100) %>%
  #   select(-mean_value, -sd_value) %>%
  #   as.data.frame()
  # 
  # taxonomy_selected <- taxonomy %>%
  #   group_by(.dots = dots) %>%
  #   summarize()
  # 
  # abundance_taxa_compact <- abundance_taxa %>%
  #   group_by(sediment, transport, class) %>%
  #   summarize(sum_percent_mean = sum(percent_mean)) %>%
  #   mutate("percent_mean" = sum_percent_mean) %>%
  #   select(-sum_percent_mean) %>%
  #   arrange(sediment, transport, !!!taxadepth) %>%
  #   left_join(taxonomy_selected, by = taxadepth) %>%
  #   select(sediment, transport, !!!taxalevels_depth, percent_mean) %>%
  #   as.data.frame()
  # 
  # for (i in rev(seq_along(taxalevels_depth))) {
  #   abundance_taxa_compact <- abundance_taxa_compact %>%
  #     mutate(!!as.name(taxalevels_depth[i]) := ifelse(is.na(!!as.name(taxalevels_depth[i])), !!as.name(taxadepth), !!as.name(taxalevels_depth[i]))) %>%
  #     as.data.frame()
  # }
  # 
  # abundance_taxa <- abundance_taxa_compact
  
  ### manual percentage calculation for paper report
  options(scipen = 999)

  relab_long_class <- relab_long[["class"]]
  relab_overall <- relab_long_class %>%
    group_by(sediment, transport, phylum, class) %>%
    summarize(mean_rel_ab = mean(rel_abundance),
              sd_rel_ab = sd(rel_abundance),
              n_bac = n()) %>%
    ungroup() %>%
    arrange(sediment, transport, -mean_rel_ab) %>%
    group_by(phylum, class) %>%
    mutate(flag = ifelse(max(mean_rel_ab) >= threshold, 1, 0)) %>%
    ungroup() %>%
    filter(flag == 1) %>%
    as.data.frame()

  relab_paper <- relab_long_class %>%
    group_by(sediment, transport, domain, phylum, class) %>%
    summarize(mean_rel_ab = mean(rel_abundance),
              sd_rel_ab = sd(rel_abundance)) %>%
    ungroup() %>%
    group_by(domain, phylum, class) %>%
    mutate(flag = ifelse(max(mean_rel_ab) >= threshold, 1, 0)) %>%
    ungroup() %>%
    filter(flag == 1) %>%
    as.data.frame()
  
} else if (parameter %in% autotrophs) {
  taxalevels_depth <- "morphotype"
  abundance <- cbind("morphotype" = row.names(otu_all), otu_all)
  
  # convert to long format for gg-plotting
  abundance_long <- reshape2::melt(abundance, id.vars = "morphotype", variable.name = "sample_ID", value.name = "abundance")
  abundance_long <- merge(cbind(sample_ID = row.names(meta_all), meta_all), abundance_long, by = "sample_ID")
  abundance_long <- abundance_long[order(as.numeric(as.character(abundance_long$sample_ID))), ]
  abundance_long$abundance <- as.numeric(abundance_long$abundance)
  
  # if(parameter == "algae") {
  #   taxalevels_depth <- "algae_class"
  #   abundance_long <- abundance_long %>%
  #     mutate(algae_class = case_when(morphotype %in% c("Merismopedia", "Pseudoanabaena") ~ "cyanobacteria",
  #                                    morphotype == "Scenedesmus" ~ "green_algae",
  #                                    TRUE ~ "diatoms"))
  # } else {
  #   
  # }
  # calculate means from replicates
  abundance_means <- aggregate(abundance_long$abundance, by = as.list(abundance_long[, which(names(abundance_long) %in% c(var1, var2, taxalevels_depth))]), function(x) round(mean(x, na.rm = T), 3))
  abundance_sd <- aggregate(abundance_long$abundance, by = as.list(abundance_long[, which(names(abundance_long) %in% c(var1, var2, taxalevels_depth))]), function(x) round(sd(x, na.rm = T), 3))
  
  names(abundance_means)[length(names(abundance_means))] <- "mean_value"
  names(abundance_sd)[length(names(abundance_sd))] <- "sd_value"
  abundance_means <- abundance_means %>%
    left_join(abundance_sd, by = c(var1, var2, taxalevels_depth))
  
  abundance_taxa <- abundance_means %>%
    group_by(transport, sediment) %>%
    mutate(percent_mean = mean_value / sum(mean_value, na.rm = T) * 100,
           percent_sd = sd_value / sum(mean_value, na.rm = T) * 100) %>%
    select(-mean_value, -sd_value) %>%
    ungroup() %>%
    as.data.frame()
  
  abundance_taxa_single <- abundance_long %>%
    rename(absolute_abundance = abundance) %>%
    group_by(transport, sediment, replicate) %>%
    mutate(percent_abundance = absolute_abundance / sum(absolute_abundance, na.rm = T) * 100) %>%
    as.data.frame()
  
  abundance_rel <- abundance_long %>%
    group_by(sample_ID, id, transport, sediment, replicate) %>%
    mutate(rel_abun = abundance / sum(abundance, na.rm = T) * 100) %>%
    arrange(transport, sediment, -rel_abun) %>%
    ungroup()
  
  abundance_paper2 <- abundance_rel %>%
    filter(transport != "Initial") %>%
    group_by(transport, morphotype) %>%
    #group_by(transport, morphotype) %>%
    summarize(mean_rel_ab = mean(rel_abun, na.rm = T), sd_rel_ab = sd(rel_abun, na.rm = T), n_count = n()) %>%
    arrange(-mean_rel_ab) %>%
    ungroup()
  
  abundance_taxa_zoom <- NULL
  abundance_taxa_zoom_single <- NULL
}

abundance_taxa <- abundance_taxa %>%
  mutate(logclass = ifelse(transport == "Initial", "Initial", "Final"),
         parameter = parameter) %>%
  relocate(parameter, transport, sediment, logclass) %>%
  as.data.frame()

abundance_taxa_single <- abundance_taxa_single %>%
  mutate(logclass = ifelse(transport == "Initial", "Initial", "Final"),
         parameter = parameter) %>%
  relocate(parameter, transport, sediment, logclass) %>%
  select(-sample_ID, -id) %>%
  as.data.frame()

#### (3) ALPHA DIVERSITY ####
# select alpha diversity indices
alpha_indices <- c("Observed", "Shannon")

## CALCULATION OF ALPHA INDICESE
if (parameter %in% heterotrophs) {
  # Rarefy for even sequencing depth and calculate indices
  rafdata <- rarefy_even_depth(seqdata, rngseed = T, trimOTUs = T)
  rarefaction_plot <- rarecurve(otu_t, step = 1000, label = FALSE)
  
  alphadiv <- estimate_richness(rafdata)
  row.names(alphadiv) <- gsub("X", "", row.names(alphadiv))
  if (all(row.names(alphadiv) == row.names(meta_all)) == FALSE) {
    warning("Data.frames not merged: Row names of alphadiv and meta_all do not match!")
  }
  alphadiv <- cbind(meta_all, alphadiv) # combine design data with indices
  alphadiv <- alphadiv[, which(names(alphadiv) %in% c(names(meta_all), alpha_indices))]
} else if (parameter %in% autotrophs) {
  # calculate diversity indices
  H <- diversity(otu_t) # Shannon
  simpson <- diversity(otu_t, "simpson") # Simpson
  inv.simpson <- diversity(otu_t, "inv") # Inverse Simpson
  S <- specnumber(otu_t) # Observed (=species richness) (rowSums(BCI > 0) does the same...)
  J <- H / log(S) ## Pielou's evenness

  # combine to single df
  #alphadiv <- data.frame(meta_all, "Observed" = S, "Shannon" = H, "Simpson" = simpson, "InvSimpson" = inv.simpson, "P_evenness" = J)
  alphadiv <- data.frame(meta_all, "Observed" = S, "Shannon" = H)
}

alphadiv <- alphadiv %>%
  mutate(logclass = ifelse(transport == "Initial", "Initial", "Final"),
         parameter = parameter) %>%
  relocate(parameter, id, transport, sediment, replicate, logclass) %>%
  as.data.frame()

alphadiv_long <- alphadiv %>%
  pivot_longer(-c(parameter, id, transport, sediment, logclass, replicate), names_to = "variable", values_to = "value") %>%
  as.data.frame()

alpha_indices_calculated <- names(alphadiv)[which(!(names(alphadiv) %in% c("parameter", "id", var1, var2, "replicate", "logclass")))]

alphadiv_mean <- aggregate(x = alphadiv[, which(names(alphadiv) %in% alpha_indices_calculated)], by = as.list(alphadiv[names(alphadiv) %in% c(var1, var2, "logclass")]), mean)
# alphadiv_sd<-aggregate(x=alphadiv[ ,which(names(alphadiv) %in% alpha_indices_calculated)], by=as.list(alphadiv[names(alphadiv) %in% c(var1, var2)]), FUN=function(x) sd(x)/sqrt(length(x))) # standard error
alphadiv_sd <- aggregate(x = alphadiv[, which(names(alphadiv) %in% alpha_indices_calculated)], by = as.list(alphadiv[names(alphadiv) %in% c(var1, var2, "logclass")]), FUN = function(x) sd(x)) # standard deviation

alphadiv_mean_long <- alphadiv_mean %>%
  pivot_longer(!c(transport, sediment, logclass), names_to = "variable", values_to = "mean") %>%
  as.data.frame()
alphadiv_sd_long <- alphadiv_sd %>%
  pivot_longer(!c(transport, sediment, logclass), names_to = "variable", values_to = "sd") %>%
  as.data.frame()
alphadiv_means <- left_join(alphadiv_mean_long, alphadiv_sd_long, by = c("transport", "sediment", "logclass", "variable"))

if (parameter == "algae") {
  param <- "total"
} else {
  param <- parameter
}

if (skip_stats == FALSE) {
  #### pANOVA for alpha diversity indices and abundance of groups ####
  
  abundance_group_stats <- abundance_group_long %>%
    filter(sediment != "Terrestrial") %>% # stats without terrestrial data
    filter(., grepl(param, variable)) %>%
    mutate(variable = gsub(pattern = param, "abundance", variable)) %>%
    as.data.frame()
  
  aovdat <- alphadiv_long %>%
    filter(sediment != "Terrestrial") %>% # stats without terrestrial data
    select(-parameter) %>%
    rbind.data.frame(abundance_group_stats) %>%
    arrange(transport, sediment, logclass, replicate, variable) %>%
    droplevels() %>%
    mutate(variable = factor(variable), logclass = factor(logclass)) %>%
    as.data.frame()
  aovdat$logclass <- factor(aovdat$logclass, levels = rev(levels(aovdat$logclass)))
  
  # set contrasts
  options(contrasts = c("contr.sum", "contr.poly")) # required for Anova Type III SS
  
  # create lists for results
  template_list <- vector("list", length = length(levels(aovdat$variable)))
  aovdatlist <- template_list
  aovlist_inter <- template_list
  aovlist_add <- template_list
  emmlist_inter <- template_list
  emmlist_add <- template_list
  hedgeslist_inter <- template_list
  hedgeslist_add <- template_list
  bestmodellist <- template_list
  
  template_names <- levels(aovdat$variable)
  names(aovdatlist) <- template_names
  names(aovlist_inter) <- template_names
  names(aovlist_add) <- template_names
  names(emmlist_inter) <- template_names
  names(emmlist_add) <- template_names
  names(hedgeslist_inter) <- template_names
  names(hedgeslist_add) <- template_names
  names(bestmodellist) <- template_names
  
  # set individual contrasts (only contrasts of interests)
  IniMix <- c(1, 0, 0, 0, 0, 0)
  IniAqu <- c(0, 1, 0, 0, 0, 0)
  MigMix <- c(0, 0, 1, 0, 0, 0)
  MigAqu <- c(0, 0, 0, 1, 0, 0)
  StatMix <- c(0, 0, 0, 0, 1, 0)
  StatAqu <- c(0, 0, 0, 0, 0, 1)
  
  custom_contrasts <- list(
    "Mixed Initial - Aquatic Initial" = IniMix - IniAqu,
    "Mixed Initial - Mixed Migrating" = IniMix - MigMix,
    "Mixed Initial - Mixed Stationary" = IniMix - StatMix,
    "Aquatic Initial - Aquatic Migrating" = IniAqu - MigAqu,
    "Aquatic Initial - Aquatic Stationary" = IniAqu - StatAqu,
    "Mixed Migrating - Mixed Stationary" = MigMix - StatMix,
    "Mixed Migrating - Aquatic Migrating" = MigMix - MigAqu,
    "Mixed Stationary - Aquatic Stationary" = StatMix - StatAqu,
    "Aquatic Migrating - Aquatic Stationary" = MigAqu - StatAqu
  )
  
  mylevs <- levels(aovdat$variable)[-1]
  
  # calculate permANOVA and post-hoc test (selected pairwise comparisons) # here with aovp (lmPerm), can alternatively be done with aovperm (permuco)
  for (m in mylevs) {
    aovdat_sub <- aovdat %>%
      filter(variable == m) %>%
      droplevels() %>%
      mutate(treatment = paste(sediment, transport, sep = "-")) %>%
      mutate(treatment = as.factor(treatment))
    aovdatlist[[m]] <- aovdat_sub
  
    aovdat_sub_summary <- aovdat_sub %>%
      group_by(sediment, transport, treatment, variable, logclass) %>%
      filter(!is.na(value)) %>%
      summarize(replicates = n()) %>%
      relocate(sediment, transport, treatment, variable, logclass, replicates) %>%
      as.data.table()
  
    levelcomb <- combn(aovdat_sub_summary$treatment, 2)
    levelcombnames <- paste(levelcomb[1, ], levelcomb[2, ], sep = " - ")
    replsum <- combn(aovdat_sub_summary$replicates, 2, FUN = sum)
    combinations <- data.frame(contrast = levelcombnames, replicates = replsum)
  
    # pANOVA
    set.seed(111)
    model_add <- aovp(value ~ sediment + transport, aovdat_sub, perm = "Exact", seqs = F, maxIter = 99999)
    summary(model_add)
    set.seed(111)
    model_inter <- aovp(value ~ sediment * transport, aovdat_sub, perm = "Exact", seqs = F, maxIter = 99999)
    summary(model_inter)
    emmip(model_inter, sediment ~ transport) # plot linear prediction
  
    AIC_model_inter <- AIC(model_inter)
    AIC_model_add <- AIC(model_add)
    if (AIC_model_inter < AIC_model_add & abs(AIC_model_inter - AIC_model_add) > 2) {
      print("Best model: Interaction")
      bestmodellist[[m]] <- "interaction"
    } else if (AIC_model_inter > AIC_model_add | abs(AIC_model_inter - AIC_model_add) < 2) {
      print("Best model: Addition")
      bestmodellist[[m]]<- "addition"
    }
  
    df_inter <- as.data.frame(summary(model_inter)[[1]])
    df_inter$factor <- row.names(df_inter)
    df_inter <- df_inter %>%
      mutate(
        variable = m
      ) %>%
      relocate(variable, factor)
    row.names(df_inter) <- 1:nrow(df_inter)
    aovlist_inter[[m]] <- df_inter
  
    df_add <- as.data.frame(summary(model_add)[[1]])
    df_add$factor <- row.names(df_add)
    df_add <- df_add %>%
      mutate(
        variable = m
      ) %>%
      relocate(variable, factor)
    row.names(df_add) <- 1:nrow(df_add)
    aovlist_add[[m]] <- df_add
  
    # emmeans
    emm_inter <- emmeans(model_inter, specs = pairwise ~ sediment * transport) # calculate estimated marginal means
    emm_transport <- emmeans(model_add, specs = pairwise ~ transport) # calculate estimated marginal means
    emm_sediment <- emmeans(model_add, specs = pairwise ~ sediment) # calculate estimated marginal means
  
    # main contrasts (p-values)
    emm.contrasts_inter <- contrast(emm_inter, # custom contrasts (only selected pairwise comparisons)
      method = custom_contrasts,
      # method = "pairwise",
      adjust = "mvt"
    ) %>% as.data.frame()
  
    emm.contrasts_transport <- contrast(emm_transport, # custom contrasts (only selected pairwise comparisons)
      method = "pairwise",
      adjust = "mvt"
    ) %>% as.data.frame()
  
    emm.contrasts_sediment <- contrast(
      emm_sediment, # custom contrasts (only selected pairwise comparisons)
      method = "pairwise",
      adjust = "mvt"
    ) %>%
      as.data.frame()
  
    # calculate Hedge's g with package "effectsize"
    # addition model
    treatment_contrasts <- emm.contrasts_inter$contrast
    aovdat_aquatic <- aovdat_sub %>%
      filter(sediment == "Aquatic")
    aovdat_mixed <- aovdat_sub %>%
      filter(sediment == "Mixed")
    aovdat_migrating <- aovdat_sub %>%
      filter(transport == "Migrating")
    aovdat_stationary <- aovdat_sub %>%
      filter(transport == "Stationary")
    aovdat_initial <- aovdat_sub %>%
      filter(transport == "Initial")
    
    if (all(aovdat_aquatic$value == 0) | (parameter == "fungi" & m %in% c("Observed", "Shannon"))) {
      hedges_aquatic_MigStat <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_aquatic_IniStat <- hedges_aquatic_MigStat
      hedges_aquatic_IniMig <- hedges_aquatic_MigStat
    } else {
      hedges_aquatic_MigStat <- hedges_g(value ~ transport, data = aovdat_aquatic %>% filter(transport != "Initial"))
      hedges_aquatic_IniStat <- hedges_g(value ~ transport, data = aovdat_aquatic %>% filter(transport != "Migrating"))
      hedges_aquatic_IniMig <- hedges_g(value ~ transport, data = aovdat_aquatic %>% filter(transport != "Stationary"))
    }
    
    if (all(aovdat_mixed$value == 0)) {
      hedges_mixed_MigStat <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_mixed_IniStat <- hedges_mixed_MigStat
      hedges_mixed_IniMig <- hedges_mixed_MigStat
    } else {
      hedges_mixed_MigStat <- hedges_g(value ~ transport, data = aovdat_mixed %>% filter(transport != "Initial"))
      hedges_mixed_IniStat <- hedges_g(value ~ transport, data = aovdat_mixed %>% filter(transport != "Migrating"))
      hedges_mixed_IniMig <- hedges_g(value ~ transport, data = aovdat_mixed %>% filter(transport != "Stationary"))
    }
    
    if (all(aovdat_migrating$value == 0)) {
      hedges_migrating <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_migrating <- hedges_g(value ~ sediment, data = aovdat_migrating)
    }
    
    if (all(aovdat_stationary$value == 0)) {
      hedges_stationary <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_stationary <- hedges_g(value ~ sediment, data = aovdat_stationary)
    }
    
    if (all(aovdat_initial$value == 0)) {
      hedges_initial <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_initial <- hedges_g(value ~ sediment, data = aovdat_initial)
    }
    
    df_hedges_inter <- data.frame(
      "contrast" = treatment_contrasts,
      "hedges_g" = c(
        hedges_initial$Hedges_g,
        hedges_mixed_IniMig$Hedges_g,
        hedges_mixed_IniStat$Hedges_g,
        hedges_aquatic_IniMig$Hedges_g,
        hedges_aquatic_IniStat$Hedges_g,
        hedges_mixed_MigStat$Hedges_g,
        hedges_migrating$Hedges_g,
        hedges_stationary$Hedges_g,
        hedges_aquatic_MigStat$Hedges_g
      ),
      "CI_lower" = c(
        hedges_initial$CI_low,
        hedges_mixed_IniMig$CI_low,
        hedges_mixed_IniStat$CI_low,
        hedges_aquatic_IniMig$CI_low,
        hedges_aquatic_IniStat$CI_low,
        hedges_mixed_MigStat$CI_low,
        hedges_migrating$CI_low,
        hedges_stationary$CI_low,
        hedges_aquatic_MigStat$CI_low
      ),
      "CI_upper" = c(
        hedges_initial$CI_high,
        hedges_mixed_IniMig$CI_high,
        hedges_mixed_IniStat$CI_high,
        hedges_aquatic_IniMig$CI_high,
        hedges_aquatic_IniStat$CI_high,
        hedges_mixed_MigStat$CI_high,
        hedges_migrating$CI_high,
        hedges_stationary$CI_high,
        hedges_aquatic_MigStat$CI_high
      )
    )
    df_hedges_inter <- df_hedges_inter %>%
      mutate("interpretation" = interpret_hedges_g(df_hedges_inter$hedges_g, rules = "gignac2016")) %>%
      left_join(emm.contrasts_inter, by = "contrast") %>%
      mutate(signif = case_when(p.value < 0.05 ~ TRUE, TRUE ~ FALSE))
    
    # transport model
    treatment_contrasts_transport <- emm.contrasts_transport$contrast
    
    # sediment model
    treatment_contrasts_sediment <- emm.contrasts_sediment$contrast
    
    if (all(aovdat_sub$value == 0)) {
      hedges_transport_MigStat <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_transport_IniMig <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_transport_IniStat <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_sediment <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_transport_MigStat <- hedges_g(value ~ transport, data = aovdat_sub %>% filter(transport != "Initial"))
      hedges_transport_IniMig <- hedges_g(value ~ transport, data = aovdat_sub %>% filter(transport != "Migrating"))
      hedges_transport_IniStat <- hedges_g(value ~ transport, data = aovdat_sub %>% filter(transport != "Stationary"))
      hedges_sediment <- hedges_g(value ~ sediment, data = aovdat_sub)
    }
    
    df_hedges_transport <- data.frame(
      "contrast" = treatment_contrasts_transport,
      "hedges_g" = c(hedges_transport_IniMig$Hedges_g, 
                     hedges_transport_IniStat$Hedges_g, 
                     hedges_transport_MigStat$Hedges_g),
      "CI_lower" = c(hedges_transport_IniMig$CI_low, 
                     hedges_transport_IniStat$CI_low, 
                     hedges_transport_MigStat$CI_low),
      "CI_upper" = c(hedges_transport_IniMig$CI_high, 
                     hedges_transport_IniStat$CI_high, 
                     hedges_transport_MigStat$CI_high)
    )
    
    df_hedges_sediment <- data.frame(
      "contrast" = treatment_contrasts_sediment,
      "hedges_g" = hedges_sediment$Hedges_g,
      "CI_lower" = hedges_sediment$CI_low,
      "CI_upper" = hedges_sediment$CI_high
    )
    
    df_hedges_transport <- df_hedges_transport %>%
      mutate("interpretation" = interpret_hedges_g(df_hedges_transport$hedges_g, rules = "gignac2016")) %>%
      left_join(emm.contrasts_transport, by = "contrast") %>%
      mutate(signif = case_when(p.value < 0.05 ~ TRUE, TRUE ~ FALSE))
    
    df_hedges_sediment <- df_hedges_sediment %>%
      mutate("interpretation" = interpret_hedges_g(df_hedges_sediment$hedges_g, rules = "gignac2016")) %>%
      left_join(emm.contrasts_sediment, by = "contrast") %>%
      mutate(signif = case_when(p.value < 0.05 ~ TRUE, TRUE ~ FALSE))
    
    # combine both factor of addition model
    emm.contrasts_add <- emm.contrasts_transport %>%
      rbind(emm.contrasts_sediment)
    df_hedges_add <- df_hedges_transport %>%
      rbind(df_hedges_sediment)
    
    # write to output list
    emmlist_inter[[m]] <- emm.contrasts_inter
    hedgeslist_inter[[m]] <- df_hedges_inter
    emmlist_add[[m]] <- emm.contrasts_add
    hedgeslist_add[[m]] <- df_hedges_add
    
    # add variable
    emmlist_inter[[m]]$variable <- m
    hedgeslist_inter[[m]]$variable <- m
    emmlist_add[[m]]$variable <- m
    hedgeslist_add[[m]]$variable <- m
  }
  
  # combine to single dt
  aovtable <- rbindlist(aovlist_add)
  bestmodeltable <- unlist(bestmodellist)
  dt_bestmodel <- data.table(parameter = names(bestmodeltable), best_model = bestmodeltable)
  
  emm_inter <- rbindlist(emmlist_inter)
  emm_add <- rbindlist(emmlist_add)
  
  hedges_inter <- rbindlist(hedgeslist_inter) %>%
    mutate(variable = factor(variable, levels = mylevs)) %>%
    mutate(contrast = factor(contrast))
  hedges_add <- rbindlist(hedgeslist_add) %>%
    mutate(variable = factor(variable, levels = mylevs)) %>%
    mutate(contrast = factor(contrast)) %>%
    mutate(contrast = recode_factor(contrast,
                                    "Mixed - Aquatic" = "sediment"
    ))
}

#### (4) BETA-DIVERSITY WITHOUT TERRESTRIAL SAMPLES ####
otu_t_in_use <- otu_t_stats
meta_in_use <- meta_stats
envir_in_use <- envir_stats

# transform & standardize data -> Hellinger transformation (asymmetrical transformation to reduce problems with double zero data.. best transformation for density data)
otu_t_notrans <- otu_t_in_use # without transformation
otu_t_hell <- decostand(otu_t_in_use, method = "hellinger") # with transformation
#otu_t_scaled <- otu_table(otu_t_hell, taxa_are_rows = FALSE)

# choose according to Betadisper if to use Hellinger transformed data
if (parameter %in% autotrophs | parameter == "fungi") {
  otudata <- otu_t_hell
  #otudata <- otu_t_scaled
} else if (parameter == "bacteria") {
  otudata <- otu_t_notrans
}

# trim otu's to dataset without Chloroplast, unclassified etc.
if (parameter %in% heterotrophs) {
  otudata <- otudata[, colnames(otudata) %in% row.names(taxonomy)]
}

# remove zero columns (ASVs that do not occur in any sample)
ncol_before <- ncol(otudata)
otudata <- otudata[, colSums(otudata, na.rm = TRUE) != 0]
ncol_after <- ncol(otudata)
print(paste(ncol_before - ncol_after, "ASVs removed"))

if (skip_stats == FALSE) {
  # factors to test
  transport <- as.factor(meta_in_use$transport)
  names(transport) <- row.names(meta_in_use)
  sediment <- as.factor(meta_in_use$sediment)
  names(sediment) <- row.names(meta_in_use)
  
  # create distance matrix
  Fdis <- vegdist(otudata, method = "bray", na.rm = TRUE) # generate a similarity distance matrix
  
  ## betadisper (= multivariate version of Levene's test for homogeneity of variances)
  # sediment
  Fbeta_sed <- betadisper(Fdis, sediment)
  set.seed(111)
  permutest(Fbeta_sed, permutations = 999)
  TukeyHSD(Fbeta_sed)
  boxplot(Fbeta_sed)
  
  # transport
  Fbeta_trans <- betadisper(Fdis, transport)
  set.seed(111)
  permutest(Fbeta_trans, permutations = 999)
  TukeyHSD(Fbeta_trans)
  boxplot(Fbeta_trans)
  
  # CAUTION: order of factors matters!! (Type 2 SS)
  # overall test
  set.seed(111)
  permanova_otu_all <- adonis2(otudata ~ transport * sediment, data = meta_in_use, permutations = 999, method = "bray", by = NULL)
  # test by terms
  set.seed(111)
  permanova_otu_interact <- adonis2(otudata ~ transport * sediment, data = meta_in_use, permutations = 999, method = "bray")
  set.seed(111)
  permanova_otu_add <- adonis2(otudata ~ transport + sediment, data = meta_in_use, permutations = 999, method = "bray")
  
  # calculate AIC corrected for small sample size (AICc)
  permanova_interact_AIC <- AICc_permanova2(permanova_otu_interact)
  permanova_interact_AIC
  permanova_add_AIC <- AICc_permanova2(permanova_otu_add)
  permanova_add_AIC
  
  df_permanova_otu_interact <- permanova_otu_interact %>% 
    as.data.frame() %>%
    rownames_to_column(var = "variable")
  df_permanova_otu_add <- permanova_otu_add %>% 
    as.data.frame() %>%
    rownames_to_column(var = "variable")

  # NMDS: Non metric multidimensional scaling
  # set dissimilarity index
  if (parameter %in% heterotrophs) {
    disind <- "bray" # Bray-Curtis
  } else if (parameter %in% autotrophs) {
    # disind <- "kul" # Kulczynski
    disind <- "bray"
  }
  
  otu_nmds <- metaMDS(otudata, k = 2, distance = disind, autotransform = FALSE, trace = TRUE)
  
  # calling enviromental parameters as vectors, it shows the correlation of vectors to the B-diversity
  env_selected <- envir_in_use %>%
    select(CR:SRP)
  row.names(env_selected) <- row.names(envir_in_use)
  
  # fit environmental vectors to nMDS
  env_selected_metamatch <- env_selected %>%
    filter(row.names(.) %in% row.names(meta_in_use))
  set.seed(111)
  ef <- envfit(otu_nmds, env_selected_metamatch, perm = 999, na.rm = TRUE) # env factors fit
  ef
  
  #### BIOENV ####
  
  # without initial samples
  initial_ids <- meta_stats$id[meta_stats$transport == "Initial"]
  meta_stats_no_init <- meta_stats %>%
    filter(!(id %in% initial_ids))
  
  # choose environmental variables for bioenv procedure
  if (parameter %in% autotrophs) {
    env_bioenv <- env_selected %>%
      select(CR, NCP, DOC, NH4.N, NOx.N, SRP, bacteria, fungi)
  } else if (parameter %in% heterotrophs) {
    env_bioenv <- env_selected %>%
      select(CR, NCP, DOC, NH4.N, NOx.N, SRP, cyanobacteria, diatoms)
  }
  
  env_bioenv_metamatch <- env_bioenv %>%
    filter(row.names(.) %in% row.names(meta_stats))
  env_bioenv_metamatch_no_init <- env_bioenv %>%
    filter(row.names(.) %in% row.names(meta_stats_no_init))
  
  # standardize data
  env_bioenv_scaled <- as.data.frame(scale(env_bioenv))
  env_bioenv_scaled_metamatch <- as.data.frame(scale(env_bioenv_metamatch))
  env_bioenv_scaled_metamatch_no_init <- as.data.frame(scale(env_bioenv_metamatch_no_init))
  
  env_bioenv_scaled_metamatch_no_init_nona <- env_bioenv_scaled_metamatch_no_init %>%
    drop_na()
  otudata_no_init_nona <- otudata %>%
    filter(row.names(.) %in% row.names(env_bioenv_scaled_metamatch_no_init_nona)) %>%
    filter(row.names(.) %in% row.names(meta_stats_no_init))
  remaining_samples <- row.names(otudata_no_init_nona)
  
  # code transport and layer as dummy variables
  transport_no_init_nona <- transport[names(transport) %in% remaining_samples] %>% droplevels()
  sediment_no_init_nona <- sediment[names(sediment) %in% remaining_samples] %>% droplevels()
  
  dummy_transport <- scale(as.numeric(transport_no_init_nona == "Migrating"))
  dummy_sediment <- scale(as.numeric(sediment_no_init_nona == "Mixed"))
  
  env_bioenv_scaled_metamatch_no_init_nona_compl <- cbind("Transport" = dummy_transport, 
                                                          "Sediment" = dummy_sediment, 
                                                          env_bioenv_scaled_metamatch_no_init_nona)
  
  # select for most dominant environmental variables
  # with transport/sediment
  sol_compl <- bioenv(otudata_no_init_nona, env_bioenv_scaled_metamatch_no_init_nona_compl,
                      index = "bray", upto = ncol(env_bioenv_scaled_metamatch_no_init_nona_compl), trace = FALSE, partial = NULL,
                      method = "spearman",
                      metric = c("euclidean", "mahalanobis", "manhattan", "gower")
  )
  summary(sol_compl)
  
  sol_compl_summary <- summary(sol_compl)
  dominant_env_compl <- unlist(str_split(sol_compl_summary$variables[which.max(sol_compl_summary$correlation)], pattern = " "))
  dominant_env_compl
  dominant_env_compl_scaled <- env_bioenv_scaled_metamatch_no_init_nona_compl[colnames(env_bioenv_scaled_metamatch_no_init_nona_compl) %in% dominant_env_compl]
  
  # mantel test
  Fdis.hel_nona <- vegdist(otudata_no_init_nona, method = "bray", na.rm = TRUE) # use bray
  
  set.seed(111)
  res_mantel <- mantel(Fdis.hel_nona, vegdist(dominant_env_compl_scaled, method = "euc"), method = "spearman", permutations = 9999)
  mantel_corr <- res_mantel$statistic
  mantel_corr
  mantel_p <- res_mantel$signif
  mantel_p
  
  # without transport/sediment
  sol <- bioenv(otudata_no_init_nona, env_bioenv_scaled_metamatch_no_init_nona,
                index = "bray", upto = ncol(env_bioenv_scaled_metamatch_no_init_nona), trace = FALSE, partial = NULL,
                method = "spearman",
                metric = c("euclidean", "mahalanobis", "manhattan", "gower")
  )
  summary(sol)
  
  sol_summary <- summary(sol)
  dominant_env <- unlist(str_split(sol_summary$variables[which.max(sol_summary$correlation)], pattern = " "))
  dominant_env
  dominant_env_scaled <- env_bioenv_scaled_metamatch_no_init_nona[colnames(env_bioenv_scaled_metamatch_no_init_nona) %in% dominant_env]
  
  #### dbRDA ####
  
  #dominant_env_to_use <- dominant_env_scaled
  dominant_env_to_use <- dominant_env_compl_scaled
  
  if (parameter == "bacteria") {
    dbRDA <- dbrda(otudata_no_init_nona ~ Transport + Sediment + SRP + cyanobacteria + diatoms, dominant_env_to_use, dist = "bray")
  } else if (parameter == "fungi") {
    dbRDA <- dbrda(otudata_no_init_nona ~ Transport + DOC, dominant_env_to_use, dist = "bray")
  } else if (parameter == "algae") {
    dbRDA <- dbrda(otudata_no_init_nona ~ Sediment + NCP + NH4.N + SRP, dominant_env_to_use, dist = "bray")
  } else if (parameter == "diatoms") {
    dbRDA <- dbrda(otudata_no_init_nona ~ Transport + NH4.N + SRP, dominant_env_to_use, dist = "bray")
  }
  
  summary(dbRDA)
  set.seed(111)
  anova.cca(dbRDA) # overall test of the significant of the analysis
  set.seed(111)
  anova.cca(dbRDA, by = "axis", perm.max = 500) # test axes for significance
  set.seed(111)
  anova.cca(dbRDA, by = "terms")
}

#### (5) VENN DIAGRAM PREPARATION ####

# split count data by treatment
# id_init_aqua <- meta_stats %>%
#   mutate(id = row.names(.)) %>%
#   filter(transport == "Initial", sediment == "Aquatic") %>%
#   pull(id)
# id_init_mix <- meta_stats %>%
#   mutate(id = row.names(.)) %>%
#   filter(transport == "Initial", sediment == "Mixed") %>%
#   pull(id)
# 
# id_mig_aqua <- meta_stats %>%
#   mutate(id = row.names(.)) %>%
#   filter(transport == "Migrating", sediment == "Aquatic") %>%
#   pull(id)
# id_mig_mix <- meta_stats %>%
#   mutate(id = row.names(.)) %>%
#   filter(transport == "Migrating", sediment == "Mixed") %>%
#   pull(id)
# 
# id_stat_aqua <- meta_stats %>%
#   mutate(id = row.names(.)) %>%
#   filter(transport == "Stationary", sediment == "Aquatic") %>%
#   pull(id)
# id_stat_mix <- meta_stats %>%
#   mutate(id = row.names(.)) %>%
#   filter(transport == "Stationary", sediment == "Mixed") %>%
#   pull(id)
id_init_aqua <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Initial", sediment == "Aquatic") %>%
  pull(id)
id_init_mix <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Initial", sediment == "Mixed") %>%
  pull(id)
id_init_terr <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Initial", sediment == "Terrestrial") %>%
  pull(id)

id_mig_aqua <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Migrating", sediment == "Aquatic") %>%
  pull(id)
id_mig_mix <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Migrating", sediment == "Mixed") %>%
  pull(id)

id_stat_aqua <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Stationary", sediment == "Aquatic") %>%
  pull(id)
id_stat_mix <- meta_all %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "Stationary", sediment == "Mixed") %>%
  pull(id)

# list containing ASV in each treatment
ASV_init_aqua <- otu[c(id_init_aqua)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_init_mix <- otu[c(id_init_mix)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_init_terr <- otu[c(id_init_terr)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_mig_aqua <- otu[c(id_mig_aqua)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_mig_mix <- otu[c(id_mig_mix)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_stat_aqua <- otu[c(id_stat_aqua)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_stat_mix <- otu[c(id_stat_mix)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

# venn_list_full <- list("Init-Aqua" = ASV_init_aqua,
#                        "Init-Mix" = ASV_init_mix,
#                        "Mig-Aqua" = ASV_mig_aqua,
#                        "Mig-Mix" = ASV_mig_mix,
#                        "Stat-Aqua" = ASV_stat_aqua,
#                        "Stat-Mix" = ASV_stat_mix)
venn_list_full <- list("Init-Terr" = ASV_init_terr,
                       "Init-Mix" = ASV_init_mix,
                       "Init-Aqua" = ASV_init_aqua,
                       "Mig-Mix" = ASV_mig_mix,
                       "Mig-Aqua" = ASV_mig_aqua,
                       "Stat-Mix" = ASV_stat_mix,
                       "Stat-Aqua" = ASV_stat_aqua)

# venn_list <- list("Mig-Aqua" = ASV_mig_aqua, 
#                   "Mig-Mix" = ASV_mig_mix,
#                   "Stat-Aqua" = ASV_stat_aqua, 
#                   "Stat-Mix" = ASV_stat_mix)
venn_list <- list("Init-Terr" = ASV_init_terr, 
                  "Mig-Mix" = ASV_mig_mix,
                  "Mig-Aqua" = ASV_mig_aqua,
                  "Stat-Mix" = ASV_stat_mix,
                  "Stat-Aqua" = ASV_stat_aqua)


#### SAVING OUTPUT ####

if (save_output == TRUE) {
  
  ## TABLES
  # absolute abundance
  if ((parameter == autotrophs[1]) | (parameter == heterotrophs[1])) {
    parameter_new <- ifelse(parameter == autotrophs[1], "autotrophs", "heterotrophs")
    write.table(abundance_group_means, file.path(dir_out, dir_tables, paste0(paste3("Abundance_means", parameter_new, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
    write.table(abundance_group, file.path(dir_out, dir_tables, paste0(paste3("Abundance_single", parameter_new, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
  }
  
  # relative abundance
  write.table(abundance_taxa, file.path(dir_out, dir_tables, paste0(paste3("Relative_abundance_means", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
  write.table(abundance_taxa_single, file.path(dir_out, dir_tables, paste0(paste3("Relative_abundance_single", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
  if (!is.null(abundance_taxa_zoom)) {
    write.table(abundance_taxa_zoom, file.path(dir_out, dir_tables, paste0(paste3("Relative_abundance_means_zoom", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
    write.table(abundance_taxa_zoom_single, file.path(dir_out, dir_tables, paste0(paste3("Relative_abundance_single_zoom", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
  }
  
  # alpha diversity
  write.table(alphadiv_means, file.path(dir_out, dir_tables, paste0(paste3("Alphadiversity_means", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
  write.table(alphadiv, file.path(dir_out, dir_tables, paste0(paste3("Alphadiversity_single", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)

  # venn lists (as RDS objects)
  saveRDS(venn_list_full, file.path(dir_out, dir_tables, filename_venn_full))
  saveRDS(venn_list, file.path(dir_out, dir_tables, filename_venn))
  
  ## STATS
  if (skip_stats == FALSE) {
    
    # beta diversity (as RDS objects)
    write.table(meta_in_use, file.path(dir_out, dir_tables, paste0(paste3("Meta_betadiversity", parameter, project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
    saveRDS(ef, file.path(dir_out, dir_tables, filename_vectors))
    saveRDS(otu_nmds, file.path(dir_out, dir_tables, filename_nmds))
    
    # dbRDA (as RDS objects)
    saveRDS(dbRDA, file.path(dir_out, dir_tables, filename_dbrda))
    
    # alpha diversity and abundance (pANOVA)
    write.table(aovtable, file.path(dir_out, dir_stats, paste0("alpha_abundance_aovp_inter_", parameter, ".txt")), sep = "\t", row.names = FALSE)
    write.table(hedges_inter, file.path(dir_out, dir_stats, paste0("alpha_abundance_emmeans_hedges_inter_", parameter, ".txt")), sep = "\t", row.names = FALSE)
    write.table(hedges_add, file.path(dir_out, dir_stats, paste0("alpha_abundance_emmeans_hedges_add_", parameter, ".txt")), sep = "\t", row.names = FALSE)
    
    # beta diversity (PERMANOVA)
    write.table(df_permanova_otu_add, file.path(dir_out, dir_stats, paste0("beta_permanova_add_", parameter, ".txt")), sep = "\t", row.names = FALSE)
    write.table(df_permanova_otu_interact, file.path(dir_out, dir_stats, paste0("beta_permanova_interact_", parameter, ".txt")), sep = "\t", row.names = FALSE)
  }
}
