#### RIPLIGHT1 Plotting of community data ####
#### Last updated: August 16, 2024

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
  "broman",
  "broom",
  "car",
  "colorspace",
  "ComplexUpset",
  "cowplot",
  "data.table",
  "dplyr",
  "emmeans",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggsci",
  "ggVennDiagram",
  "grid",
  "gridExtra",
  "here",
  "lattice",
  "lme4",
  "lmerTest",
  "lmPerm",
  "magick",
  "patchwork",
  "patternplot",
  "phyloseq",
  "png",
  "RColorBrewer",
  "scales",
  "stringr",
  "tidyverse",
  "vegan",
  "VennDetail"
)

using(packages)

#### USER SETTINGS ####

# determine in- and output directories
dir_in <- "D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Results" # where is all your data
# dir_in <- "D:/RipLight1/Results"
dir_out <- "D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Results/Plots" # where the plots are saved
# dir_out <- "D:/Riplight1/Results/Plots"
dir_metab <- "Metabolism"
dir_struc <- "Structure & Abundance"
dir_tables <- "Tables"
dir_stats <- "Stats"

# set variables
parameters <- c("Algae", "Diatoms", "Bacteria", "Fungi")
parameterclass <- c("Autotrophs", "Heterotrophs")
parameters <- tolower(parameters)
parameterclass <- tolower(parameterclass)
project <- "RipLight1" # to name output files

# organism labeling
heterotrophs <- c("bacteria", "fungi")
autotrophs <- c("algae", "diatoms")

# include terrestrial samples?
# include.terrestrial <- FALSE # if TRUE, terrestrial samples (only available for heterotrophs) will be included in the NMDS and alpha diversity stats, if FALSE, they will be dropped

# filenames
# metabolism
file_metab <- "Metabolism_stats_single_RipLight1.txt"
filepath_metab <- here(dir_in, dir_metab, dir_tables, file_metab)
file_metab_means <- "Metabolism_stats_means_RipLight1.txt"
filepath_metab_means <- here(dir_in, dir_metab, dir_tables, file_metab_means)
file_metab_effsize <- "Metabolism_stats_effsize_add_RipLight1.txt"
filepath_metab_effsize <- here(dir_in, dir_metab, dir_stats, file_metab_effsize)
file_metab_all <- "Metabolism_stats_single_all_RipLight1.txt"
filepath_metab_all <- here(dir_in, dir_metab, dir_tables, file_metab_all)

# file patterns
# abundance
filepattern_abun <- "Abundance_single"
filepattern_abun_means <- "Abundance_means"

# alpha diversity
filepattern_alpha <- "Alphadiversity_single"
filepattern_alpha_means <- "Alphadiversity_means"

# relative abundance
filepattern_relabun <- "Relative_abundance_means"

# beta diversity
filepattern_nmds <- "beta_nmds"
filepattern_meta <- "Meta_betadiversity"
filepattern_dbrda <- "dbRDA"
filepattern_vectors <- "beta_vectors"
filepattern_venn_full <- "venn_full"
filepattern_venn <- "venn"

# effect size
filepattern_alpha_effsize <- "alpha_abundance_emmeans_hedges_add"

# colours & shapes
col_aquatic <- brewer.pal(6, "Blues")
# col_aquatic <- brewer.pal(6, "YlGnBu")
col_mixed <- brewer.pal(6, "Oranges")
# col_mixed <- brewer.pal(6, "YlOrBr")
col_terrestrial <- brewer.pal(7, "BrBG")[1]
# myfill <- c("white", col_aquatic[5], col_aquatic[4], "white", col_mixed[5], col_mixed[4])
myfill <- c("white", "white", col_mixed[5], col_mixed[4], "white", col_aquatic[5], col_aquatic[4])
# mycolour <- c(col_aquatic[5], "black", "black", col_mixed[5], "black", "black")
mycolour <- c(col_terrestrial, col_mixed[5], "black", "black", col_aquatic[5], "black", "black")
mycolour2 <- c(col_mixed[5], col_mixed[5], col_mixed[4], col_aquatic[5], col_aquatic[5], col_aquatic[4])

# without initial samples
myfill_no_init <- myfill[-c(1, 2, 5)]
mycolour_no_init <- mycolour[-c(1, 2, 5)]

myshape <- c(21, 24, 24, 21, 24, 24) # c(21, 24, 23) # define shapes

# # mycolour_barplots <- c(col_aquatic[4], col_mixed[4], col_aquatic[5], col_mixed[5], col_aquatic[5], col_mixed[5])
# mycolour_barplots <- c("white", "white", col_aquatic[5], col_mixed[5], col_aquatic[4], col_mixed[4])
# # myfill_barplots <- c(col_aquatic[4], col_mixed[4], col_aquatic[5], col_mixed[5], col_aquatic[5], col_mixed[5])
# myfill_barplots <- c("black", "black", col_aquatic[5], col_mixed[5], col_aquatic[4], col_mixed[4])

# mycolour_reduced <- c(col_aquatic[5], col_mixed[5])
# myfill_reduced <- mycolour_reduced

legend_key <- "rect"

# output
format <- ".png"
dpi <- 300

# output filenames
# Fig. 1
filepattern_metab_plot <- "Metabolism"
# Fig. 2
filepattern_abundance_plot <- "Abundance"
filepattern_relabundance_plot <- "Relative_abundance"
filepattern_abundance_relab_plot <- "Abundance_&_Rel_abundance"
# Fig. 3
filepattern_richness_plot <- "Richness"
filepattern_richness_legend <- "Richness_legend"
filepattern_venn_plot <- "Venn"
filepattern_richness_venn_plot <- "Richness_&_Venn"
filepattern_richness_upset_plot <- "Richness_&_Upset"

# Fig. 4
filepattern_betadiversity_plot <- "nMDS"
filepattern_betadiversity_legend <- "nMDS_legend"
# Fig. 5
filepattern_effectsize_plot <- "effectsizes"
# Fig. Sxx
filepattern_shannon_plot <- "Shannon"
filepattern_venn_plot2 <- "Venn2"
filepattern_metab_dyn_plot <- "Metabolism_dynamics"
filepattern_dbRDA_plot <- "dbRDA"
filepattern_dbRDA_legend <- "dbRDA_legend"

# variable 1: sediment type
var1 <- "sediment"
levels_var1 <- c("Terrestrial" = "Terrestrial", "Mixed" = "Mixed", "Aquatic" = "Aquatic") # how factor 1 (sediment type) shall be renamed

# variable 2: transport state
var2 <- "transport"
levels_var2 <- c("Initial" = "Initial", "Migrating" = "Migrating ripple", "Stationary" = "Stationary") # how factor 2 (transport regime) shall be renamed

# for plotting
levels_treatment <- c(
  "Terrestrial-Initial",
  "Mixed-Initial",
  "Mixed-Migrating",
  "Mixed-Stationary",
  "Aquatic-Initial",
  "Aquatic-Migrating",
  "Aquatic-Stationary"
)

levels_treatment_no_init <- c(
  "Mixed-Migrating",
  "Mixed-Stationary",
  "Aquatic-Migrating",
  "Aquatic-Stationary"
)
# levels_treatment_no_init_short <- c("Aqua-Mig", "Aqua-Stat", "Mix-Mig", "Mix-Stat")
levels_treatment_no_init_short <- c("Terr", "Mix-M", "Mix-S", "Aqua-M", "Aqua-S")

# levels
autotroph_levels <- c("cyanobacteria" = "Cyanobacteria", "green_algae" = "Green algae", "diatoms" = "Diatoms", "total" = "Total")
heterotroph_levels <- c("bacteria" = "16S", "fungi" = "ITS")

#### END OF SETTINGS ###

#### FUNCTIONS ####
round_any <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

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

# ggplot theme (from RipDiv)
mytheme <- function(facet_labs = TRUE, x_labs = TRUE, no_space = FALSE) {
  finish_theme <- theme(
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10, unit = "pt"),
    plot.title = element_text(size = rel(1.5), face = "bold", hjust = 1),
    axis.text.x = element_text(size = 16, angle = 0), # , hjust = 0.8, vjust = 0.8),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks.x.bottom = element_line(linewidth = 0.5),
    axis.ticks.y.left = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_blank(), # line(colour="grey", size = 0.3),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.background = element_rect(fill = "white", colour = "white"),
    strip.background.x = element_rect(fill = "white"),
    strip.background.y = element_rect(fill = "white"),
    strip.text = element_text(size = 18),
    # legend.key.width = unit(0.1, "cm"), # 0.3
    # legend.key.height = unit(0.1, "cm"), # 0.3
    # legend.key = element_blank(),
    legend.key = element_rect(colour = "black"),
    title = element_text(size = 10, face = "plain")
  )

  if (facet_labs == FALSE) {
    finish_theme <- finish_theme +
      theme(strip.text = element_blank())
  }

  if (x_labs == FALSE) {
    finish_theme <- finish_theme +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x.bottom = element_blank()
      )
  }

  if (no_space == TRUE) {
    finish_theme <- finish_theme +
      theme(
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.grid = element_blank()
      )
  }

  return(finish_theme)
}

theme_no_space <- theme(
  plot.margin = unit(c(0, 0, 0, 0), "null"),
  panel.grid = element_blank()
)

## CREATE NICER BREAKS
nicerbreaks <- function(low = NULL, high = NULL, rng = NULL, length.out = 8, round_to_mult = NULL) {
  ## FUNCTION ARGUMENTS
  # low:                 [numeric] lower boundary
  # high:                [numeric] higher bounday
  # rng:                 [vector]  instead of low and high, specify vector of both (rng=c(low, high))
  # length.out:          [numeric] length of break vector that will be returned
  # round_to_mult:       [numeric, 1<=x<=10] if specified, will round to the nearest multiple of this value*10^x, with the exponent x depending on the data (so round_to_mult can be 1,2,3,4,5,6,7,8,9,10)

  # helper functions
  round_any <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }

  getdig <- function(x) nchar(sub("\\.[0-9]+", "", abs(x)))

  findexp <- function(x) {
    findexp.single <- function(y) {
      if (!(is.numeric(y))) {
        stop("x must be numeric")
      }

      # if (any(y==0)) {
      #   stop("Math Error: x cannot be 0")
      # }

      if (y > 0 & y < 1) {
        myfun <- function(z) -(as.numeric(sub("[0-9]+(?:\\.[0-9]+)?e-", "", format(z, scientific = TRUE))))
      } else if (y < 0 & y > (-1)) {
        myfun <- function(z) -(as.numeric(sub("-[0-9]+(?:\\.[0-9]+)?e-", "", format(z, scientific = TRUE))))
      } else if (y >= 1) {
        myfun <- function(z) as.numeric(sub("[0-9]+(?:\\.[0-9]+)?e+", "", format(z, scientific = TRUE)))
      } else if (y <= (-1)) {
        myfun <- function(z) as.numeric(sub("-[0-9]+(?:\\.[0-9]+)?e+", "", format(z, scientific = TRUE)))
      } else if (y == 0) {
        myfun <- function(z) 0 * z
      }

      ret <- myfun(y)

      return(ret)
    }

    if (length(x) == 1) {
      ret <- findexp.single(x)
    } else {
      ret <- sapply(x, findexp.single)
    }

    return(ret)
  }

  # check input arguments
  if (all(is.null(c(rng, low, high)))) {
    stop("Function argument is missing")
  }

  if (!is.null(rng) & (!is.null(low) | !is.null(high))) {
    warning("Too many arguments supplied, rng will be used")
  } else if (is.null(rng) & (is.null(low) | is.null(high))) {
    stop("Both low and high argument have to be specified if no rng is supplied")
  }

  if (!is.null(rng)) {
    if (length(rng) != 2) {
      stop("Rng has to be a vector with two values only")
    }
  } else {
    if (length(c(high, low)) != 2) {
      stop("Low and high need to be one-dimensional")
    }
  }

  # main part
  if (is.null(rng)) {
    orig <- c(low, high)
  } else {
    orig <- rng
  }

  if (max(orig) < 1) {
    myexp <- min(findexp(orig))
    new <- orig * 10^(-myexp)
  } else {
    new <- orig
  }

  ## CREATE DEFAULT BINS
  accur <- max(10^(getdig(new[1]) - 1), 10^(getdig(new[2]) - 1))
  vallims <- c(round_any(new[1], accur, f = floor), round_any(new[2], accur, f = ceiling))
  # vallims<-c(round_any(new[1], 10^(getdig(new[1])-1), f=floor), round_any(new[2], 10^(getdig(new[2])-1), f=ceiling))
  bins <- seq(vallims[1], vallims[2], length = length.out)

  # bins<-seq(new[1], new[2], length=length.out)

  if (max(orig) < 1) {
    bins <- bins * 10^(myexp)
  }

  ## CHANGE BINS
  lower <- min(bins)
  upper <- max(bins)
  # increment<-(upper-lower)/(length(bins)-1)

  # if(increment<round_to_mult) {
  #   round_to_mult<-NULL
  # }

  if (!is.null(round_to_mult)) {
    if (min(findexp(bins)) < 0) {
      newaccuracy <- round_to_mult * 10^-(max(abs(findexp(bins))))
    } else {
      newaccuracy <- round_to_mult * 10^(max(abs(findexp(bins))) - 1)
    }

    # newlower<-round_any(lower, accuracy = newaccuracy, f=trunc)
    # newupper<-round_any(upper, accuracy = newaccuracy, f=round)
    # bins<-seq(newlower, newupper, by=newaccuracy)
    newincrement <- min(diff(unique(round_any(bins, newaccuracy))))

    bins <- seq(lower, upper, by = newincrement) # better solution
    if (min(bins) > min(new)) {
      bins <- c(min(bins) - newincrement, bins)
    }
    if (max(bins) < max(new)) {
      bins <- c(bins, max(bins) + newincrement)
    }

    # bins<-sort(unique(as.numeric(as.character(c(bins, lower, upper)))))
  }

  # trim to adequate breaks (only one bin higher than max of range)

  minind <- ifelse(length(which(bins <= min(orig))) == 0, 1, max(which(bins <= min(orig))))
  maxind <- ifelse(length(which(bins >= max(orig))) == 0, length(bins), min(which(bins >= max(orig))))
  bins <- bins[minind:maxind]

  return(bins)
}

## FIND EXPONENT
findexp <- function(x) {
  findexp.single <- function(y) {
    if (!(is.numeric(y))) {
      stop("x must be numeric")
    }

    # if (any(y==0)) {
    #   stop("Math Error: x cannot be 0")
    # }

    if (y > 0 & y < 1) {
      myfun <- function(z) -(as.numeric(sub("[0-9]+(?:\\.[0-9]+)?e-", "", format(z, scientific = TRUE))))
    } else if (y < 0 & y > (-1)) {
      myfun <- function(z) -(as.numeric(sub("-[0-9]+(?:\\.[0-9]+)?e-", "", format(z, scientific = TRUE))))
    } else if (y >= 1) {
      myfun <- function(z) as.numeric(sub("[0-9]+(?:\\.[0-9]+)?e+", "", format(z, scientific = TRUE)))
    } else if (y <= (-1)) {
      myfun <- function(z) as.numeric(sub("-[0-9]+(?:\\.[0-9]+)?e+", "", format(z, scientific = TRUE)))
    } else if (y == 0) {
      myfun <- function(z) 0 * z
    }

    ret <- myfun(y)

    return(ret)
  }

  if (length(x) == 1) {
    ret <- findexp.single(x)
  } else {
    ret <- sapply(x, findexp.single)
  }

  return(ret)
}

#### DATA IMPORT ####

## METABOLISM

metab <- read.table(file = filepath_metab, sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
metabmeans <- read.table(file = filepath_metab_means, sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
metab_all <- read.table(file = filepath_metab_all, sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
metab_effsize <- read.table(file = filepath_metab_effsize, sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)

metab <- metab %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment_no_init)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  droplevels()

metabmeans <- metabmeans %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment_no_init)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  droplevels()

metab_all <- metab_all %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment_no_init)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  droplevels()

metab_effsize <- metab_effsize %>%
  mutate(logclass = str_to_title(logclass)) %>%
  mutate(variable = paste(logclass, variable, sep = " ")) %>%
  select(-logclass)

## STRUCTURE

# create lists for data
template_list <- vector("list", length = length(parameterclass))
names(template_list) <- parameterclass
abundancelist <- template_list
abundancemeanlist <- template_list

template_list <- vector("list", length = length(parameters))
names(template_list) <- parameters
alphalist <- template_list
alphameanlist <- template_list
nmdslist <- template_list
metalist <- template_list
vectorslist <- template_list
dbrdalist <- template_list
relablist <- template_list
vennfulllist <- template_list
vennlist <- template_list
effsizelist <- template_list

for (i_param in seq_along(parameterclass)) {
  current_param <- parameterclass[i_param]

  # complete filenames
  file_abun <- paste0(paste(filepattern_abun, current_param, project, sep = "_"), ".txt")
  file_abun_means <- paste0(paste(filepattern_abun_means, current_param, project, sep = "_"), ".txt")

  # import data
  abundancelist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_tables, file_abun), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
  abundancemeanlist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_tables, file_abun_means), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
}

# turn lists into data.frames
abun <- bind_rows(abundancelist, .id = "parameterclass")
abunmean <- bind_rows(abundancemeanlist, .id = "parameterclass")

for (i_param in seq_along(parameters)) {
  current_param <- parameters[i_param]

  # complete filenames
  file_alpha <- paste0(paste(filepattern_alpha, current_param, project, sep = "_"), ".txt")
  file_alpha_means <- paste0(paste(filepattern_alpha_means, current_param, project, sep = "_"), ".txt")
  file_relabun <- paste0(paste(filepattern_relabun, current_param, project, sep = "_"), ".txt")
  file_nmds <- paste0(paste(filepattern_nmds, current_param, sep = "_"), ".rds")
  file_meta <- paste0(paste(filepattern_meta, current_param, project, sep = "_"), ".txt")
  file_vectors <- paste0(paste(filepattern_vectors, current_param, sep = "_"), ".rds")
  file_dbrda <- paste0(paste(filepattern_dbrda, current_param, sep = "_"), ".rds")
  file_venn_full <- paste0(paste(filepattern_venn_full, current_param, sep = "_"), ".rds")
  file_venn <- paste0(paste(filepattern_venn, current_param, sep = "_"), ".rds")
  file_effsize <- paste0(paste(filepattern_alpha_effsize, current_param, sep = "_"), ".txt")

  # import data
  alphalist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_tables, file_alpha), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
  alphameanlist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_tables, file_alpha_means), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
  relablist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_tables, file_relabun), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
  nmdslist[[i_param]] <- readRDS(file = here(dir_in, dir_struc, dir_tables, file_nmds))
  metalist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_tables, file_meta), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
  vectorslist[[i_param]] <- readRDS(file = here(dir_in, dir_struc, dir_tables, file_vectors))
  dbrdalist[[i_param]] <- readRDS(file = here(dir_in, dir_struc, dir_tables, file_dbrda))
  vennfulllist[[i_param]] <- readRDS(file = here(dir_in, dir_struc, dir_tables, file_venn_full))
  vennlist[[i_param]] <- readRDS(file = here(dir_in, dir_struc, dir_tables, file_venn))
  effsizelist[[i_param]] <- read.table(file = here(dir_in, dir_struc, dir_stats, file_effsize), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
}

# turn lists into data.frames
meta <- bind_rows(metalist, .id = "parameter")
alpha <- bind_rows(alphalist, .id = "parameter")
alphamean <- bind_rows(alphameanlist, .id = "parameter")
relab <- bind_rows(relablist, .id = "parameter")
effsize <- bind_rows(effsizelist, .id = "parameter") %>%
  mutate(variable = str_to_title(fct_recode(variable, "abundance" = "abundance_value", "Richness" = "Observed"))) %>%
  mutate(variable = paste(str_to_title(parameter), variable, sep = " ")) %>%
  mutate(variable = factor(variable)) %>%
  mutate(variable = str_replace_all(variable, "Diatoms", "Diatom")) %>%
  select(-parameter)

# rename factor levels
abun <- abun %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  filter(!(parameterclass == "heterotrophs" & variable == "total")) %>%
  droplevels()

abunmean <- abunmean %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  filter(!(parameterclass == "heterotrophs" & variable == "total")) %>%
  droplevels()

alpha <- alpha %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  droplevels()

alphamean <- alphamean %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  droplevels()

relab <- relab %>%
  mutate(treatment = paste(sediment, transport, sep = "-")) %>%
  mutate(treatment = factor(treatment, levels = levels_treatment)) %>%
  mutate(sediment = recode(sediment, !!!levels_var1)) %>%
  mutate(sediment = factor(sediment, levels = levels_var1)) %>%
  mutate(transport = recode(transport, !!!levels_var2)) %>%
  mutate(transport = factor(transport, levels = levels_var2)) %>%
  droplevels()

#### (0) METABOLISM ALL ####

# plotdat <- metab_all %>%
#   filter(variable == "CR") %>%
#   droplevels() %>%
#   mutate(resp = -resp) %>% # display CR as -CR
#   mutate(logclass = factor(logclass, levels = c(paste("Day", 1:12)))) %>%
#   as.data.frame()
#
# y.min <- 0
# y.max <- 1
# br <- nicerbreaks(y.min, y.max, round_to_mult = 0.5)
#
# ytitle <- bquote("-CR" ~ "[" * mu * gO[2] ~ h^{
#   -1
# } ~ gDW^
#   {
#     -1
#   } * "]")
#
# # create plot
# cr_plot_all <- ggplot(data = plotdat, aes(x = as.factor(cycle), y = resp, fill = treatment)) +
#   geom_boxplot(outlier.shape = NA, key_glyph = legend_key) +
#   scale_fill_manual(name = "Sediment", values = myfill_no_init) +
#   # stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, alpha = 0.9, colour = "black", fill = "red") +
#   geom_point(
#     shape = 21, size = 3, fill = "grey20", colour = "gray60",
#     position = position_jitter(width = 0.2, height = 0)
#   ) +
#   geom_hline(yintercept = 0, colour = "grey", linetype = "solid", size = 0.2) +
#   # geom_vline(xintercept = 2.5, linetype = "dotted", size = 0.2) +
#   # scale_x_discrete(labels = rep(levels_transport[c(1:2)], 2)) +
#   scale_y_continuous(breaks = br, limits = c(min(br), max(br))) +
#   labs(title = NULL, fill = NULL, x = NULL, y = ytitle) +
#   # draw_label("Sed.type: 0.043; Transp: 0.725", x = 0, y = br[length(br)], hjust = -0.3, vjust = 0.6) +
#   mytheme() +
#   facet_grid(. ~ sediment) +
#   # theme(legend.key = element_rect(colour="black")) +
#   theme(legend.position = "none")
# cr_plot_all

#### (1) METABOLISM INITIAL/FINAL/DYNAMICS ####

## CR
cr_plot <- vector("list", length = length(levels(metab$logclass)))
names(cr_plot) <- levels(metab$logclass)

for (i_var in seq_along(levels(metab$logclass))) {
  current_var <- levels(metab$logclass)[i_var]

  plotdat <- metab %>%
    filter(variable == "CR", logclass == current_var) %>%
    droplevels() %>%
    mutate(resp = -resp) %>% # display CR as -CR
    as.data.frame()

  if (current_var == "dynamics") {
    plotdat <- plotdat %>%
      mutate(resp = resp / (10^(-3)))

    y.min <- -0.5
    y.max <- 3
    br <- nicerbreaks(y.min, y.max, round_to_mult = 1)

    ytitle <- bquote(.(str_to_title(current_var)) ~ "CR" ~ "[" * 10^{
      -3
    } ~ mu * gO[2] ~ h^{
      -2
    } ~ gDW^
      {
        -1
      } * "]")
  } else {
    y.min <- 0
    y.max <- 1
    br <- nicerbreaks(y.min, y.max, round_to_mult = 0.5)

    ytitle <- bquote(.(str_to_title(current_var)) ~ "CR" ~ "[" * mu * gO[2] ~ h^{
      -1
    } ~ gDW^
      {
        -1
      } * "]")
  }

  facet_labs <- TRUE
  x_labs <- ifelse(current_var == "dynamics", TRUE, FALSE)
  no_space <- TRUE

  # create plot
  cr_plot[[i_var]] <- ggplot(data = plotdat, aes(x = transport, y = resp, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, key_glyph = legend_key, size = 1) +
    scale_fill_manual(name = "Sediment", values = myfill_no_init) +
    # stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, alpha = 0.9, colour = "black", fill = "red") +
    geom_point(
      shape = 21, size = 3, fill = "grey20", colour = "gray60",
      position = position_jitter(width = 0.2, height = 0)
    ) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "solid", size = 0.6) +
    # geom_vline(xintercept = 2.5, linetype = "dotted", size = 0.2) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    scale_y_continuous(breaks = br, limits = c(min(br), max(br))) +
    labs(title = NULL, fill = NULL, x = NULL, y = ytitle) +
    # draw_label("Sed.type: 0.043; Transp: 0.725", x = 0, y = br[length(br)], hjust = -0.3, vjust = 0.6) +
    facet_grid(. ~ sediment) +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    theme(legend.position = "none")
  # cr_plot[[i_var]] <- ggdraw(cr_plot[[i_var]]) +
  #   draw_plot_label(x = 0.13, y = 0.9, label = str_to_title(current_var), hjust = 0, size = 16, fontface = "plain")
}

## NCP
ncp_plot <- vector("list", length = length(levels(metab$logclass)))
names(ncp_plot) <- levels(metab$logclass)

for (i_var in seq_along(levels(metab$logclass))) {
  current_var <- levels(metab$logclass)[i_var]

  plotdat <- metab %>%
    filter(variable == "NCP", logclass == current_var) %>%
    droplevels() %>%
    as.data.frame()

  if (current_var == "dynamics") {
    plotdat <- plotdat %>%
      mutate(resp = resp / (10^(-3)))

    y.min <- -0.5
    y.max <- 20
    br <- nicerbreaks(y.min, y.max, round_to_mult = 5)

    ytitle <- bquote(.(str_to_title(current_var)) ~ "NCP" ~ "[" * 10^{
      -3
    } ~ mu * gO[2] ~ h^{
      -2
    } ~ gDW^
      {
        -1
      } * "]")
  } else {
    y.min <- -0.5
    y.max <- ifelse(current_var == "initial", 0.2, 6)
    br <- nicerbreaks(y.min, y.max, round_to_mult = ifelse(current_var == "initial", 0.1, 1))

    ytitle <- bquote(.(str_to_title(current_var)) ~ "NCP" ~ "[" * mu * gO[2] ~ h^{
      -1
    } ~ gDW^
      {
        -1
      } * "]")
  }

  facet_labs <- ifelse(current_var == "dynamics", TRUE, FALSE)
  x_labs <- TRUE
  no_space <- TRUE

  # create plot
  ncp_plot[[i_var]] <- ggplot(data = plotdat, aes(x = transport, y = resp, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, key_glyph = legend_key, size = 1) +
    scale_fill_manual(name = "Sediment", values = myfill_no_init) +
    # stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, alpha = 0.9, colour = "black", fill = "red") +
    geom_point(
      shape = 21, size = 3, fill = "grey20", colour = "gray60",
      position = position_jitter(width = 0.2, height = 0)
    ) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "solid", size = 0.6) +
    # geom_vline(xintercept = 2.5, linetype = "dotted", size = 0.2) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    scale_y_continuous(breaks = br, limits = c(min(br), max(br))) +
    labs(title = NULL, fill = NULL, x = NULL, y = ytitle) +
    # draw_label("Sed.type: 0.043; Transp: 0.725", x = 0, y = br[length(br)], hjust = -0.3, vjust = 0.6) +
    facet_grid(. ~ sediment) +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    theme(legend.position = "none")
}

## GPP
gpp_plot <- vector("list", length = length(levels(metab$logclass)))
names(gpp_plot) <- levels(metab$logclass)

for (i_var in seq_along(levels(metab$logclass))) {
  current_var <- levels(metab$logclass)[i_var]

  plotdat <- metab %>%
    filter(variable == "GPP", logclass == current_var) %>%
    droplevels() %>%
    as.data.frame()

  if (current_var == "dynamics") {
    plotdat <- plotdat %>%
      mutate(resp = resp / (10^(-3)))

    y.min <- 0
    y.max <- 25
    br <- nicerbreaks(y.min, y.max, round_to_mult = 5)

    ytitle <- bquote(.(str_to_title(current_var)) ~ "GPP" ~ "[" * 10^{
      -3
    } ~ mu * gO[2] ~ h^{
      -2
    } ~ gDW^
      {
        -1
      } * "]")
  } else {
    y.min <- 0
    y.max <- ifelse(current_var == "initial", 0.3, 7)
    br <- nicerbreaks(y.min, y.max, round_to_mult = ifelse(current_var == "initial", 5, 1))

    ytitle <- bquote(.(str_to_title(current_var)) ~ "GPP" ~ "[" * mu * gO[2] ~ h^{
      -1
    } ~ gDW^
      {
        -1
      } * "]")
  }

  facet_labs <- FALSE
  x_labs <- TRUE
  no_space <- TRUE

  # create plot
  gpp_plot[[i_var]] <- ggplot(data = plotdat, aes(x = transport, y = resp, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, key_glyph = legend_key, size = 1) +
    scale_fill_manual(name = "Sediment", values = myfill_no_init) +
    # stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, alpha = 0.9, colour = "black", fill = "red") +
    geom_point(
      shape = 21, size = 3, fill = "grey20", colour = "grey60",
      position = position_jitter(width = 0.2, height = 0)
    ) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "solid", size = 0.6) +
    # geom_vline(xintercept = 2.5, linetype = "dotted", size = 0.2) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    scale_y_continuous(breaks = br, limits = c(min(br), max(br))) +
    labs(title = NULL, fill = NULL, x = NULL, y = ytitle) +
    # draw_label("Sed.type: 0.043; Transp: 0.725", x = 0, y = br[length(br)], hjust = -0.3, vjust = 0.6) +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    facet_grid(. ~ sediment) +
    theme(legend.position = "none")
}

## arrange Fig. 1
metab_legend <- get_legend(cr_plot[["initial"]] + theme(
  legend.position = "right",
  legend.box.margin = margin(0, 0, 0, 12), # create some space to the left of the legend
  legend.title = element_blank()
))

metab_plot <- (cr_plot[["initial"]] + cr_plot[["final"]]) / (ncp_plot[["initial"]] + ncp_plot[["final"]]) & theme_no_space
metab_plot <- metab_plot +
  plot_annotation(tag_levels = LETTERS[1:4]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
metab_plot

## arrange dynamics (Supplem. Material)
metab_plot_dyn <- (cr_plot[["dynamics"]] + ncp_plot[["dynamics"]]) & theme_no_space
metab_plot_dyn <- metab_plot_dyn +
  plot_annotation(tag_levels = LETTERS[1:2]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
metab_plot_dyn

#### (2A) ABSOLUTE ABUNDANCE ####

abun_params <- c("bacteria", "fungi", "total", "cyanobacteria", "diatoms")
abun_plot <- vector("list", length = length(abun_params))
names(abun_plot) <- abun_params

for (i_var in seq_along(abun_params)) {
  current_var <- abun_params[i_var]

  facet_labs <- ifelse(current_var %in% c("bacteria", "cyanobacteria"), TRUE, FALSE)
  x_labs <- ifelse(current_var %in% c("diatoms", "cyanobacteria"), TRUE, FALSE)
  no_space <- TRUE

  is_autotroph <- ifelse(current_var %in% c("total", "cyanobacteria", "diatoms"), TRUE, FALSE)

  # plotdat <- abunmean %>%
  #   filter(variable == current_var) %>%
  #   droplevels() %>%
  #   as.data.frame()
  #
  # myexponent <- max(findexp(plotdat$mean_value))
  #
  # plotdat <- plotdat %>%
  #   mutate(mean_value = mean_value / 10^myexponent, sd_value = sd_value / 10^myexponent) %>%
  #   mutate(lower = mean_value, upper = mean_value + sd_value) %>%
  #   as.data.frame()

  plotdat <- abun %>%
    filter(variable == current_var) %>%
    droplevels() %>%
    as.data.frame()
  mycolour_abun <- mycolour

  # add fake (invisible) data for equal spacing in autotroph dataset (total algae and diatoms)
  if (is_autotroph == TRUE) {
    plotdat <- plotdat %>%
      add_row(parameterclass = "autotrophs", transport = "Initial", sediment = "Terrestrial", logclass = "Initial", value = mean(plotdat$value), treatment = "Terrestrial-Initial") %>%
      mutate(
        sediment = factor(sediment, levels = levels_var1),
        transport = factor(transport, levels = levels_var2),
        treatment = factor(treatment, levels = levels_treatment)
      )

    mycolour_abun[1] <- "white"
  }

  current_var <- ifelse(current_var == "total", "total algae", current_var)

  # myexponent <- max(findexp(plotdat$value))
  # plotdat <- plotdat %>%
  #   mutate(value = value / 10^myexponent) %>%
  #   as.data.frame()
  #
  # ytitle <- bquote(.(str_to_title(current_var)) ~ "[" * 10^{
  #   .(myexponent)
  # } ~ cells ~ gDW^
  #   {
  #     -1
  #   } * "]")
  #
  # y.min <- round_any(min(plotdat$value), 0.1)
  # y.max <- round_any(max(plotdat$value), 0.1)
  # br <- nicerbreaks(y.min, y.max, round_to_mult = 2)
  ytitle <- bquote(.(str_to_title(current_var)) ~ "[" * cells ~ gDW^
    {
      -1
    } * "]")

  plot_tag <- LETTERS[c(T, F)][i_var]

  # create plot
  current_plot <- ggplot(data = plotdat, aes(x = transport, y = value, colour = treatment, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, key_glyph = legend_key, size = 1) +
    scale_fill_manual(name = "Sediment", values = myfill) +
    scale_colour_manual(name = "Sediment", values = mycolour_abun) +
    # stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, alpha = 0.9, colour = "black", fill = "red") +
    geom_point(
      shape = 21, size = 3, colour = "grey60", fill = "grey20",
      position = position_jitter(width = 0.2, height = 0)
    ) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "solid", size = 0.6) +
    # geom_vline(xintercept = 2.5, linetype = "dotted", size = 0.2) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    # scale_y_continuous(breaks = br, limits = c(min(br), max(br))) +
    labs(title = NULL, fill = NULL, x = NULL, y = ytitle) +
    # labs(tag = plot_tag, title = NULL, fill = NULL, x = NULL, y = ytitle) +
    # draw_label("Sed.type: 0.043; Transp: 0.725", x = 0, y = br[length(br)], hjust = -0.3, vjust = 0.6) +
    facet_grid(. ~ sediment, space = "free_x", scales = "free_x") +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    theme(legend.position = "none")

  if (is_autotroph == TRUE) {
    current_plot <- current_plot +
      geom_point(
        data = plotdat %>% filter(sediment == "Terrestrial"),
        shape = 21, size = 20, colour = "white", fill = "white"
      )
  }

  abun_plot[[i_var]] <- current_plot
}

## arrange Fig. 2 - left
# abun_legend <- get_legend(abun_plot[[1]] + theme(
#   legend.position = "right",
#   legend.box.margin = margin(0, 0, 0, 12), # create some space to the left of the legend
#   legend.title = element_blank()
# ))

abun_plot_composed <- (abun_plot[["bacteria"]] / abun_plot[["fungi"]] / abun_plot[["total"]] + abun_plot[["diatoms"]])
abun_plot_composed <- abun_plot_composed +
  # plot_annotation(tag_levels = LETTERS[1:4]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
abun_plot_composed

# cyanobacterial abundance (Suppl. Material)
abun_plot_cyanobacteria <- abun_plot[["cyanobacteria"]] +
  theme(plot.tag = element_text(size = 20, face = "bold"))
abun_plot_cyanobacteria

#### (2B) RELATIVE ABUNDANCE ####

relab_params <- c("bacteria", "fungi", "total", "diatoms")
relab_plot <- vector("list", length = length(relab_params))
names(relab_plot) <- relab_params
relab_plot_legend <- relab_plot

for (i_var in seq_along(relab_params)) {
  current_var <- relab_params[i_var]

  facet_labs <- ifelse(current_var == "bacteria", TRUE, FALSE)
  x_labs <- ifelse(current_var == "diatoms", TRUE, FALSE)
  no_space <- TRUE

  is_autotroph <- ifelse(current_var %in% c("total", "diatoms"), TRUE, FALSE)

  ytitle <- "Relative abundance [%]"

  if (current_var == "total") {
    plotdat <- abunmean %>%
      filter(parameterclass == "autotrophs") %>%
      filter(variable != "total") %>%
      mutate(variable = factor(variable, levels = c("cyanobacteria", "diatoms", "green_algae"))) %>%
      mutate(variable = recode(variable, !!!c("cyanobacteria" = "Cyanobacteria", "diatoms" = "Diatoms", "green_algae" = "Green Algae"))) %>%
      droplevels() %>%
      as.data.frame()

    set.seed(2) # for repeatability
    fill.arg <- "variable"
  } else if (current_var == "diatoms") {
    plotdat <- relab %>%
      filter(parameter == current_var) %>%
      select(-(domain:kingdom)) %>%
      droplevels() %>%
      as.data.frame()

    set.seed(1) # for repeatability
    fill.arg <- "morphotype"
  } else if (current_var == "bacteria") {
    plotdat <- relab %>%
      filter(parameter == current_var) %>%
      select(-morphotype, -kingdom) %>%
      droplevels() %>%
      as.data.frame()

    set.seed(6) # for repeatability
    fill.arg <- "class"
    # reorder "Other" and "Unclassified"
    ord.ind <- which(levels(plotdat$class) %in% "Other")
    uncl.ind <- which(levels(plotdat$class) %in% "Unclassified")
    no.ind <- which(!(levels(plotdat$class) %in% c("Other", "Unclassified")))
    plotdat$class <- factor(plotdat$class, levels = levels(plotdat$class)[c(no.ind, uncl.ind, ord.ind)])
  } else if (current_var == "fungi") {
    plotdat <- relab %>%
      filter(parameter == current_var) %>%
      select(-morphotype, -domain) %>%
      droplevels() %>%
      as.data.frame()

    set.seed(1) # for repeatability
    fill.arg <- "class"
    # reorder "Other" and "Unclassified"
    ord.ind <- which(levels(plotdat$class) %in% "Other")
    uncl.ind <- which(levels(plotdat$class) %in% "Unclassified")
    no.ind <- which(!(levels(plotdat$class) %in% c("Other", "Unclassified")))
    plotdat$class <- factor(plotdat$class, levels = levels(plotdat$class)[c(no.ind, uncl.ind, ord.ind)])
  }

  # create colours palette
  ncols <- length(unique(plotdat[[fill.arg]]))
  colorpalette <- colorRampPalette(pal_locuszoom()(7))(ncols + 1)
  colorpalette <- colorpalette[-length(colorpalette)] # remove grey colour
  distinct.colours <- sample(colorpalette, ncols)
  if (current_var == "total") {
    distinct.colours <- distinct.colours[c(2, 1, 3)]
  }

  # add fake (invisible) data for equal spacing in autotroph dataset (total algae and diatoms)
  if (current_var == "diatoms") {
    plotdat <- plotdat %>%
      add_row(transport = "Initial", sediment = "Terrestrial", logclass = "Initial", morphotype = "Achnantes", percent_mean = 100, treatment = "Terrestrial-Initial") %>%
      mutate(
        sediment = factor(sediment, levels = levels_var1),
        transport = factor(transport, levels = levels_var2),
        treatment = factor(treatment, levels = levels_treatment)
      )
  } else if (current_var == "total") {
    plotdat <- plotdat %>%
      add_row(transport = "Initial", sediment = "Terrestrial", logclass = "Initial", variable = "Diatoms", percent_mean = 100, treatment = "Terrestrial-Initial") %>%
      mutate(
        sediment = factor(sediment, levels = levels_var1),
        transport = factor(transport, levels = levels_var2),
        treatment = factor(treatment, levels = levels_treatment)
      )
  }

  # create plot
  if (current_var %in% c("bacteria", "fungi")) {
    relab_plot[[i_var]] <- ggplot(plotdat, aes(x = transport, y = percent_mean, fill = class))
  } else if (current_var == "diatoms") {
    relab_plot[[i_var]] <- ggplot(plotdat, aes(x = transport, y = percent_mean, fill = morphotype))
  } else if (current_var == "total") {
    relab_plot[[i_var]] <- ggplot(plotdat, aes(x = transport, y = percent_mean, fill = variable))
  }

  current_var <- ifelse(current_var == "total", "total algae", current_var)
  plot_tag <- LETTERS[c(F, T)][i_var]

  current_plot <- relab_plot[[i_var]] +
    geom_bar(stat = "identity", size = 0.3, width = 0.5) +
    labs(x = NULL, y = ytitle, fill = NULL) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 105), expand = c(0, NA)) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    scale_fill_manual(values = distinct.colours) +
    guides(fill = guide_legend(ncol = 1)) +
    # labs(tag = plot_tag) +
    facet_grid(. ~ sediment, scale = "free_x", space = "free_x") +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    # theme(legend.position = "right") +
    theme(legend.justification = c(0, 0.5))

  relab_plot_legend[[i_var]] <- get_legend(
    current_plot + theme(legend.box.margin = margin(0, 0, 0, 12)) # create some space to the left of the legend
  )

  if (is_autotroph == TRUE) {
    current_plot <- current_plot +
      geom_bar(
        data = plotdat %>% filter(sediment == "Terrestrial"),
        stat = "identity", size = 0.3, width = 0.5,
        colour = "white", fill = "white"
      )
  }

  relab_plot[[i_var]] <- current_plot
}

## arrange Fig. 2 - right
relab_plot_composed <- (relab_plot[["bacteria"]] / relab_plot[["fungi"]] / relab_plot[["total"]] + relab_plot[["diatoms"]])
relab_plot_composed <- relab_plot_composed +
  # plot_annotation(tag_levels = LETTERS[1:4]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
relab_plot_composed

## arrange complete Fig. 2 - left & right
abun_relab_plot_composed <- ((abun_plot[["bacteria"]] + relab_plot[["bacteria"]]) /
  (abun_plot[["fungi"]] + relab_plot[["fungi"]]) /
  (abun_plot[["total"]] + relab_plot[["total"]]) /
  (abun_plot[["diatoms"]] + relab_plot[["diatoms"]])) & theme_no_space

abun_relab_plot_composed <- abun_relab_plot_composed +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 12)
  )
abun_relab_plot_composed

#### (3A) ALPHA DIVERSITY: OBSERVED RICHNESS (RAREFIED DATA) AND SHANNON INDEX ####

alpha_params <- c("bacteria", "fungi", "algae", "diatoms")
richness_plot <- vector("list", length = length(alpha_params))
names(richness_plot) <- alpha_params
shannon_plot <- richness_plot

for (i_param in seq_along(alpha_params)) {
  current_param <- alpha_params[i_param]

  # facet_labs <- ifelse(current_param == "bacteria", TRUE, FALSE)
  facet_labs <- TRUE
  # x_labs <- ifelse(current_param == "algae", TRUE, FALSE)
  x_labs <- TRUE
  no_space <- TRUE

  is_autotroph <- ifelse(current_param %in% c("algae", "diatoms"), TRUE, FALSE)

  if (current_param == "bacteria") {
    richness_breaks <- seq(0, 600, 100)
    shannon_breaks <- seq(0, 6, 1)
  } else if (current_param == "fungi") {
    richness_breaks <- seq(0, 300, 50)
    shannon_breaks <- seq(0, 4.0, 0.5)
  } else if (current_param == "algae") {
    richness_breaks <- seq(0, 20, 2)
    shannon_breaks <- seq(0, 2.5, 0.5)
  } else if (current_param == "diatoms") {
    richness_breaks <- seq(0, 18, 2)
    shannon_breaks <- seq(0, 2.0, 0.5)
  }

  plotdat <- alpha %>%
    filter(parameter == current_param)

  mycolour_abun <- mycolour

  # add fake (invisible) data for equal spacing in autotroph dataset (total algae and diatoms)
  if (is_autotroph == TRUE) {
    plotdat <- plotdat %>%
      add_row(parameter = "algae", transport = "Initial", sediment = "Terrestrial", logclass = "Initial", Observed = mean(plotdat$Observed), Shannon = mean(plotdat$Shannon), treatment = "Terrestrial-Initial") %>%
      mutate(
        sediment = factor(sediment, levels = levels_var1),
        transport = factor(transport, levels = levels_var2),
        treatment = factor(treatment, levels = levels_treatment)
      )

    mycolour_abun[1] <- "white"
  }

  current_param <- ifelse(current_param == "total", "total algae", current_param)
  y_label <- str_to_title(current_param)

  current_plot <- ggplot(plotdat, aes(x = transport, y = Observed, colour = treatment, fill = treatment)) +
    stat_summary(geom = "bar", fun = mean, position = position_dodge(), size = 1) +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.5, size = 1) +
    # geom_point(aes(group = treatment, fill = sediment), shape = 21, size = 3, colour = "gray60",
    #            position = position_jitter(width = 0.2, height = 0)
    # ) +
    scale_fill_manual(name = "Sediment", values = myfill) +
    scale_colour_manual(name = "Sediment", values = mycolour_abun) +
    geom_point(
      shape = 21, size = 3, fill = "grey20", colour = "gray60",
      position = position_jitter(width = 0.2, height = 0)
    ) +
    # geom_hline(yintercept = 0, colour = "grey60", linetype = "solid", size = 0.6) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    labs(x = NULL, y = paste0("Observed richness ", "(", y_label, ")"), fill = NULL, title = NULL) +
    scale_y_continuous(
      breaks = richness_breaks, limits = c(min(richness_breaks), max(richness_breaks)),
      expand = c(0, NA)
    ) +
    facet_grid(. ~ sediment, scales = "free_x", space = "free_x") +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    theme(legend.position = "none")

  if (is_autotroph == TRUE) {
    current_plot <- current_plot +
      geom_point(
        data = plotdat %>% filter(sediment == "Terrestrial"),
        shape = 21, size = 20, colour = "white", fill = "white"
      )
  }
  richness_plot[[i_param]] <- current_plot

  current_plot <- ggplot(plotdat, aes(x = transport, y = Shannon, colour = treatment, fill = treatment)) +
    stat_summary(geom = "bar", fun = mean, position = position_dodge(), size = 1) +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.5, size = 1) +
    # geom_point(aes(group = treatment, fill = sediment), shape = 21, size = 3, colour = "gray60",
    #            position = position_jitter(width = 0.2, height = 0)
    # ) +
    scale_fill_manual(name = "Sediment", values = myfill) +
    scale_colour_manual(name = "Sediment", values = mycolour_abun) +
    geom_point(
      shape = 21, size = 3, fill = "grey20", colour = "gray60",
      position = position_jitter(width = 0.2, height = 0)
    ) +
    # geom_hline(yintercept = 0, colour = "grey60", linetype = "solid", size = 0.6) +
    scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n")) +
    labs(x = NULL, y = paste0("Shannon index ", "(", y_label, ")"), fill = NULL, title = NULL) +
    scale_y_continuous(
      breaks = shannon_breaks, limits = c(min(shannon_breaks), max(shannon_breaks)),
      expand = c(0, NA)
    ) +
    facet_grid(. ~ sediment, scales = "free_x", space = "free_x") +
    mytheme(facet_labs = facet_labs, x_labs = x_labs, no_space = no_space) +
    theme(legend.position = "none")

  if (is_autotroph == TRUE) {
    current_plot <- current_plot +
      geom_point(
        data = plotdat %>% filter(sediment == "Terrestrial"),
        shape = 21, size = 20, colour = "white", fill = "white"
      )
  }
  shannon_plot[[i_param]] <- current_plot
}

# arrange Fig. 3 - left
# abun_legend <- get_legend(abun_plot[[1]] + theme(
#   legend.position = "right",
#   legend.box.margin = margin(0, 0, 0, 12), # create some space to the left of the legend
#   legend.title = element_blank()
# ))

# richness_plot_composed <- (richness_plot[["bacteria"]] / richness_plot[["fungi"]] / richness_plot[["algae"]] + richness_plot[["diatoms"]])
richness_plot_composed <- (richness_plot[["bacteria"]] / richness_plot[["fungi"]] / richness_plot[["algae"]])
richness_plot_composed <- richness_plot_composed +
  # plot_annotation(tag_levels = LETTERS[1:4]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
richness_plot_composed

shannon_plot_composed <- (shannon_plot[["bacteria"]] / shannon_plot[["fungi"]] / shannon_plot[["algae"]])
shannon_plot_composed <- shannon_plot_composed +
  # plot_annotation(tag_levels = LETTERS[1:4]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
shannon_plot_composed

#### (3B) VENN DIAGRAM/UPSET PLOTS ####

venn_params <- names(vennlist)
venn_plot <- vector("list", length = length(venn_params))
names(venn_plot) <- venn_params
venn_plot2 <- venn_plot
upset_plot <- venn_plot
upset_list <- lapply(vennlist, UpSetR::fromList)

for (i_param in seq_along(venn_params)) {
  current_param <- venn_params[i_param]

  is_autotroph <- ifelse(current_param %in% c("algae", "diatoms"), TRUE, FALSE)

  venn_list <- vennlist[[current_param]]
  # venn_list <- venn_list[c(1, 3, 2, 4)]
  venn_list <- venn_list[c(1, 2, 4, 3, 5)]
  names(venn_list) <- levels_treatment_no_init_short # [-1]
  names(venn_list)
  upset_cols <- c(col_terrestrial, myfill_no_init)

  if (is_autotroph == TRUE) {
    venn_list <- venn_list[-1]
    upset_cols <- upset_cols[-1]
  }
  dt_upset <- UpSetR::fromList(venn_list)
  
  unique_ASVs <- vector("list", length = length(venn_list))
  names(unique_ASVs) <- names(venn_list)
  shared_ASVs <- unique_ASVs
  unique_ASVs_perc <- unique_ASVs
  shared_ASVs_perc <- unique_ASVs
  total_ASVs <- nrow(dt_upset)
  
  for (i_treat in seq_along(names(venn_list))) {
    current_treatment <- names(venn_list)[i_treat]
    dt_upset_unique <- dt_upset %>%
      filter(get(current_treatment) == 1 & rowSums(.) == 1)
    unique_ASVs[[i_treat]] <- nrow(dt_upset_unique)
    unique_ASVs_perc[[i_treat]] <- paste0(round(nrow(dt_upset_unique)/total_ASVs * 100, 1), "%")
  }
  dt_upset_shared <- dt_upset %>%
    filter(rowSums(.) == 5)
  shared_ASVs <- nrow(dt_upset_shared)
  shared_ASVs_perc <- paste0(round(nrow(dt_upset_shared)/total_ASVs * 100, 2), "%")
  shared_ASVs
  unique_ASVs_perc
  
  #min_percent <- ifelse(current_param == "bacteria", 1.5, 0.5) # minimum intersection size (percent of total ASVs)
  min_percent <- 1.5
  
  if (is_autotroph == TRUE) {
    upset_plot[[i_param]] <- ggdraw(ComplexUpset::upset(dt_upset, # ggrdaw is needed to "freeze" the plot, otherwise aes() uses non-standard evaluation (ggplot in loop is a bad idea)
      names(venn_list),
      name = paste0(str_to_title(current_param), "\n", "(", total_ASVs, " morphotypes detected)"),
      sort_sets = FALSE,
      # sort_intersections = FALSE,
      height_ratio = c(0.6, 0.4),
      queries = list(
        upset_query(set = "Mix-M", fill = upset_cols[1]),
        upset_query(set = "Mix-S", fill = upset_cols[2]),
        upset_query(set = "Aqua-M", fill = upset_cols[3]),
        upset_query(set = "Aqua-S", fill = upset_cols[4])
      ),
      width_ratio = 0.2,
      stripes = upset_stripes(colors = "white"),
      set_sizes = upset_set_size() +
        ylab("Set size"),
      themes = upset_modify_themes(
        list(
          "intersections_matrix" = theme(text = element_text(size = 20)),
          "overall_sizes" = theme(axis.text.x = element_text(angle = 0))
        )
      ),
      base_annotations = list(
        # with manual aes specification:
        #"Intersection size #" = intersection_size(text_mapping = aes(label = !!get_size_mode("exclusive_intersection"))),
        "Intersection size" = intersection_size(text_mapping = aes(label = paste0(round(!!get_size_mode("exclusive_intersection") / total_ASVs * 100, 1), "%")))
      )
    ))
  } else {
    upset_plot[[i_param]] <- ggdraw(ComplexUpset::upset(dt_upset,
      names(venn_list),
      name = paste0(str_to_title(current_param), "\n", "(", total_ASVs, " ASVs detected)"),
      sort_sets = FALSE,
      min_size = floor(min_percent/100*total_ASVs),
      # sort_intersections = FALSE,
      height_ratio = c(0.6, 0.4),
      queries = list(
        upset_query(set = "Terr", fill = upset_cols[1]),
        upset_query(set = "Mix-M", fill = upset_cols[2]),
        upset_query(set = "Mix-S", fill = upset_cols[3]),
        upset_query(set = "Aqua-M", fill = upset_cols[4]),
        upset_query(set = "Aqua-S", fill = upset_cols[5])
      ),
      width_ratio = 0.2,
      stripes = upset_stripes(colors = "white"),
      set_sizes = upset_set_size() +
        ylab("Set size"),
      themes = upset_modify_themes(
        list(
          "intersections_matrix" = theme(text = element_text(size = 20)),
          "overall_sizes" = theme(axis.text.x = element_text(angle = 0))
        )
      ),
      base_annotations = list(
        # with manual aes specification:
        #"Intersection size #" = intersection_size(text_mapping = aes(label = !!get_size_mode("exclusive_intersection"))),
        "Intersection size" = intersection_size(text_mapping = aes(label = paste0(round(!!get_size_mode("exclusive_intersection") / total_ASVs * 100, 1), "%")))
      )
    ))
  }


  if (is_autotroph == FALSE) {
    venn_list <- venn_list[-1]
  }
  # df_venn <- result(venndetail(venn_list))
  # unique_items <- df_venn %>%
  #   filter(Subset %in% c("Aqua-Mig", "Aqua-Stat", "Mix-Mix", "Mix-Stat")) %>%
  #   # filter(Subset %in% c("Terr-Init", "Mig-Aqua", "Stat-Aqua", "Mig-Mix", "Stat-Mix")) %>%
  #   group_by(Subset) %>%
  #   summarize(allD = paste(Detail, collapse = ", "))
  # unique_items

  # set breaks and colours
  if (current_param == "bacteria") {
    venn_breaks <- seq(0, 1200, 200)
  } else if (current_param == "fungi") {
    venn_breaks <- seq(0, 150, 25)
  } else if (current_param == "algae") {
    venn_breaks <- seq(0, 25, 5)
  } else if (current_param == "diatoms") {
    venn_breaks <- seq(0, 25, 5)
  }
  venn_nbr <- length(venn_breaks) - 1
  venn_palette <- "Heat"
  venn_colours <- rev(sequential_hcl(venn_palette, n = venn_nbr))

  venn_plot[[i_param]] <- ggVennDiagram(venn_list,
    category.names = levels_treatment_no_init_short[-1],
    label_alpha = 0,
    label_size = 6,
    set_size = 6,
    set_colour = myfill
  ) + # ,
    # category.names = c("Mig-S","Mig-U","Sta-S", "Sta-U")) +
    scale_fill_stepsn(
      breaks = venn_breaks,
      limits = c(min(venn_breaks), max(venn_breaks)),
      colours = venn_colours
    ) +
    # scale_fill_distiller(palette = "YlOrBr", direction = 1) +
    # scale_fill_continuous(palette = sequential_hcl(n = 6, palette = "Grays")) +
    labs(fill = NULL) +
    scale_x_continuous(expand = expansion(mult = .1)) +
    scale_y_continuous(expand = expansion(mult = .1)) +
    # guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 14))
  venn_plot[[i_param]]

  # venn <- Venn(venn_list)
  # venn_data <- process_data(venn)
  # ggplot() +
  #   # change mapping of color filling
  #   geom_polygon(aes(X, Y, fill = id, group = id),
  #                data = venn_regionedge(venn_data),
  #                show.legend = FALSE) +
  #   # adjust edge size and color
  #   geom_path(aes(X, Y, color = id, group = id),
  #             data = venn_setedge(data),
  #             linewidth = 3,
  #             show.legend = FALSE) +
  #   # show set label in bold
  #   geom_text(aes(X, Y, label = name),
  #             fontface = "bold",
  #             data = venn_setlabel(data)) +
  #   # add a alternative region name
  #   geom_label(aes(X, Y, label = id),
  #              data = venn_regionlabel(data),
  #              alpha = 0.5) +
  #   coord_equal() +
  #   theme_void()

  myfill_venn <- myfill_no_init # c(col_terrestrial, myfill_no_init)
  venn_plot2[[i_param]] <- ggvenn::ggvenn(
    venn_list,
    fill_color = myfill_venn,
    fill_alpha = 0.5,
    stroke_color = NA,
    # stroke_color = "black",
    stroke_size = 0.8, set_name_size = 6,
    text_size = 6, text_color = "black"
  ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 0, r = -4, b = 0, l = -4, unit = "cm"))
  venn_plot2[[i_param]]
}

# venn_plot2_composed <- (venn_plot2[["bacteria"]] / venn_plot2[["fungi"]] / venn_plot2[["algae"]] + venn_plot2[["diatoms"]])
venn_plot2_composed <- (venn_plot2[["bacteria"]] / venn_plot2[["fungi"]] / venn_plot2[["algae"]])
venn_plot2_composed <- venn_plot2_composed +
  # plot_annotation(tag_levels = LETTERS[1:4]) &
  theme(plot.tag = element_text(size = 20, face = "bold"))
venn_plot2_composed

## arrange complete Fig. 3 - left & right
mywidths <- c(1, 1)
# richness_venn_plot_composed <- ((richness_plot[["bacteria"]] + venn_plot2[["bacteria"]] + plot_layout(widths = mywidths)) /
#   (richness_plot[["fungi"]] + venn_plot2[["fungi"]] + plot_layout(widths = mywidths)) /
#   (richness_plot[["algae"]] + venn_plot2[["algae"]] + plot_layout(widths = mywidths)) /
#   (richness_plot[["diatoms"]] + venn_plot2[["diatoms"]] + plot_layout(widths = mywidths))) & theme_no_space
richness_venn_plot_composed <- ((richness_plot[["bacteria"]] + venn_plot2[["bacteria"]] + plot_layout(widths = mywidths)) /
  (richness_plot[["fungi"]] + venn_plot2[["fungi"]] + plot_layout(widths = mywidths)) /
  (richness_plot[["algae"]] + venn_plot2[["algae"]] + plot_layout(widths = mywidths))) & theme_no_space

richness_venn_plot_composed <- richness_venn_plot_composed +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 12)
  )
richness_venn_plot_composed

## alternative Fig. 3 - richness & upset plot
upset_plot_bacteria <- wrap_elements(upset_plot[["bacteria"]])
upset_plot_fungi <- wrap_elements(upset_plot[["fungi"]])
upset_plot_algae <- wrap_elements(upset_plot[["algae"]])

myheights <- c(0.8, 1)
wrap_richness_upset_bacteria <- wrap_elements((richness_plot[["bacteria"]] / upset_plot_bacteria + plot_annotation(tag_levels = list(c('A', 'B'))) &
                                                 theme(
                                                   plot.tag = element_text(size = 20, face = "bold"),
                                                   legend.text = element_text(size = 12)
                                                 )) + plot_layout(heights = myheights))
wrap_richness_upset_fungi <- wrap_elements((richness_plot[["fungi"]] / upset_plot_fungi + plot_annotation(tag_levels = list(c('C', 'D'))) &
                                              theme(
                                                plot.tag = element_text(size = 20, face = "bold"),
                                                legend.text = element_text(size = 12)
                                              )) + plot_layout(heights = myheights))
wrap_richness_upset_algae <- wrap_elements((richness_plot[["algae"]] / upset_plot_algae + plot_annotation(tag_levels = list(c('E', 'F'))) &
                                              theme(
                                                plot.tag = element_text(size = 20, face = "bold"),
                                                legend.text = element_text(size = 12)
                                              )) + plot_layout(heights = myheights))

richness_upset_composed <- (wrap_richness_upset_bacteria + plot_spacer() + wrap_richness_upset_fungi + plot_layout(widths = c(0.49, 0.02, 0.49))) / (plot_spacer() + wrap_richness_upset_algae + plot_spacer() + plot_layout(widths = c(0.255, 0.49, 0.255)))

# richness_upset_plot_composed <- richness_upset_composed +
#   plot_annotation(tag_levels = "A") &
#   theme(
#     plot.tag = element_text(size = 20, face = "bold"),
#     legend.text = element_text(size = 12)
#   )
richness_upset_plot_composed <- richness_upset_composed
richness_upset_plot_composed

#### (4) BETA DIVERSITY ####

beta_params <- names(nmdslist)
beta_plot <- vector("list", length = length(beta_params))
names(beta_plot) <- beta_params

for (i_param in seq_along(beta_params)) {
  current_param <- beta_params[i_param]

  meta_in_use <- metalist[[current_param]] %>%
    rename(site = id)
  otu_nmds <- nmdslist[[current_param]]
  ef <- vectorslist[[current_param]]

  ef_df <- as.data.frame(scores(ef, display = "vectors"))
  ef_df <- cbind(ef_df, "R2" = (ef$vectors)$r, "p" = (ef$vectors)$pvals)
  names(ef_df) <- c("NMDS1", "NMDS2", "R2", "p")

  if (current_param %in% autotrophs) {
    ef_df_signif <- ef_df %>%
      rownames_to_column() %>%
      filter(rowname %in% c("CR", "GPP", "NCP", "DOC", "NH4.N", "NOx.N", "SRP", "bacteria", "fungi")) %>%
      column_to_rownames("rowname") # %>%
    # filter(p < 0.05)
  } else if (current_param %in% heterotrophs) {
    ef_df_signif <- ef_df %>%
      rownames_to_column() %>%
      filter(rowname %in% c("CR", "GPP", "NCP", "DOC", "NH4.N", "NOx.N", "SRP", "cyanobacteria", "diatoms")) %>%
      column_to_rownames("rowname") # %>%
    # filter(p < 0.05)
  }

  vector_lengths <- sqrt(ef_df_signif$NMDS1^2 + ef_df_signif$NMDS2^2)

  otu_nmds_ndim <- otu_nmds$ndim
  otu_nmds_stress <- otu_nmds$stress
  stressplot(otu_nmds)

  # NMDS scores for samples ("sites")
  nmds_site_scores <- as.data.frame(scores(otu_nmds, "sites")) # using the scores function from vegan to extract the site scores and convert to a data.frame
  nmds_site_scores$site <- rownames(nmds_site_scores) # create a column of site names, from the rownames of data.scores
  nmds_site_scores <- merge(nmds_site_scores, meta_in_use, by = "site")
  nmds_site_scores$treatment <- paste(nmds_site_scores$sediment, nmds_site_scores$transport, sep = "-")
  nmds_site_scores$treatment <- factor(nmds_site_scores$treatment) # , levels = c("migrating-superficial", "migrating-underlying", "stationary-superficial", "stationary-underlying"))
  nmds_site_scores <- nmds_site_scores %>%
    arrange(as.numeric(site))

  # NMDS scores for species ("species")
  nmds_species_scores <- as.data.frame(scores(otu_nmds, "species")) # using the scores function from vegan to extract the species scores and convert to a data.frame
  nmds_species_scores$species <- rownames(nmds_species_scores)

  # scores for environmental vectors ("vectors")
  ordiArrowMul(ef)
  multiplier <- case_when(
    current_param == "bacteria" ~ 1,
    current_param == "fungi" ~ 1,
    current_param == "diatoms" ~ 1,
    current_param == "algae" ~ 1,
    TRUE ~ 1
  )

  # nmds_envir_scores <- as.data.frame(scores(ef, "vectors")) * multiplier
  nmds_envir_scores <- ef_df_signif[c("NMDS1", "NMDS2")] * multiplier
  nmds_envir_scores <- nmds_envir_scores %>%
    rownames_to_column() %>%
    mutate(rowname = str_replace(rowname, pattern = "\\.", replacement = "-")) %>%
    mutate(rowname = ifelse(rowname %in% c(parameters, "cyanobacteria"), str_to_title(rowname), rowname))

  # initiate new plot
  plotdat <- nmds_site_scores
  specdat <- nmds_species_scores

  plot.new()
  ord_nmds <- invisible(ordiellipse(otu_nmds, plotdat[["treatment"]], choices = c(1, 2), display = "sites", kind = "ehull", conf = 0.95, label = T))

  # find treatments with less than 2 observations
  major_treatments <- plotdat %>%
    group_by(transport, sediment, treatment) %>%
    summarize(n_obs = n()) %>%
    filter(n_obs >= 3) %>%
    droplevels() %>%
    pull(treatment)
  plotdat_major <- plotdat %>%
    filter(treatment %in% major_treatments) %>%
    droplevels()

  df_ell_nmds <- data.frame()
  for (g in levels(plotdat_major$treatment)) {
    df_ell_nmds <- rbind(df_ell_nmds, cbind(as.data.frame(with(
      plotdat_major[plotdat_major[["treatment"]] == g, ],
      veganCovEllipse(ord_nmds[[g]]$cov, ord_nmds[[g]]$center, ord_nmds[[g]]$scale)
    )), treatment = g))
  }
  df_ell_nmds <- df_ell_nmds %>%
    separate_wider_delim(treatment, names = c("transport", "sediment"), delim = "-", cols_remove = FALSE)

  names(mycolour) <- levels(plotdat$treatment)
  names(myfill) <- levels(plotdat$treatment)

  y_label <- str_to_title(current_param)

  beta_plot[[i_param]] <- ggplot() +
    geom_path(data = df_ell_nmds, aes(x = NMDS1, y = NMDS2, colour = treatment, group = treatment), linetype = 1, show.legend = FALSE) +
    geom_point(data = plotdat, aes(x = NMDS1, y = NMDS2, colour = treatment, shape = treatment, fill = treatment), size = 3.5, stroke = 1.5) + # add the point markers
    scale_y_continuous(breaks = breaks_pretty(n = 10)) +
    scale_x_continuous(breaks = breaks_pretty(n = 10)) +
    scale_colour_manual(values = mycolour2) +
    scale_fill_manual(values = myfill) +
    scale_shape_manual(values = myshape, na.translate = FALSE) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    labs(title = NULL, fill = NULL, shape = NULL) +
    guides(
      fill = guide_legend(override.aes = list(fill = myfill, shape = myshape, colour = mycolour2)),
      colour = "none", shape = "none"
    ) +
    annotation_compass(fontsize = 16, label = paste0("Stress = ", round(otu_nmds_stress, 3), ", ", "k = ", otu_nmds_ndim), position = "SE") +
    # annotation_compass(fontsize = 18, label = "p = 0.15", position = "SE") +
    annotation_compass(label = y_label, position = "NW", fontsize = 20) +
    geom_segment(data = nmds_envir_scores, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = "grey30", show.legend = F, alpha = 0.8) +
    geom_text_repel(data = nmds_envir_scores, max.overlaps = 50, point.padding = NA, nudge_x = 0.02, nudge_y = 0.02, aes(x = NMDS1, y = NMDS2, label = rowname), parse = TRUE, colour = "black", size = 4, show.legend = F) +
    # geom_text_repel(data = nmds_envir_scores, max.overlaps = 50, point.padding = NA, nudge_x = 0.02, nudge_y = 0.02, aes(x = NMDS1, y = NMDS2, label = str_replace(row.names(nmds_envir_scores), pattern = "\\.", replacement = "-")), parse = TRUE, colour = "black", size = 4, show.legend = F) +
    mytheme(no_space = TRUE) +
    theme(legend.key = element_rect(colour = "white", fill = NA)) +
    # theme(legend.position = "none", panel.background = element_rect(colour = "black")) +
    theme(aspect.ratio = 1)
  beta_plot[[i_param]]

  if (i_param == 1) {
    beta_plot_legend <- get_legend(beta_plot[[i_param]] + theme(
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0), # create some space to the left of the legend
      legend.title = element_blank()
    ))
    filepath_legend_beta <- here(dir_out, paste0(filepattern_betadiversity_legend, format))
    ggsave(filepath_legend_beta, beta_plot_legend, width = 10, height = 10, units = "cm", dpi = dpi)
  }

  beta_plot[[i_param]] <- beta_plot[[i_param]] +
    theme(legend.position = "none")
}

# read in legend as image
legend_beta <- image_read(filepath_legend_beta) %>%
  image_trim()

# compose plot
beta_plot_composed <- (beta_plot[["bacteria"]] + beta_plot[["fungi"]]) / (beta_plot[["algae"]] + beta_plot[["diatoms"]]) & theme_no_space

# tag plot
beta_plot_composed <- beta_plot_composed +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

# add legend as image on plot
beta_plot_composed <- ggdraw() +
  draw_plot(beta_plot_composed) +
  draw_image(legend_beta, x = 0.32, y = 0.6, height = 0.6, width = 0.2, scale = 0.65)

#### (5) dbRDA ####

dbRDA_plot <- vector("list", length = length(parameters))
names(dbRDA_plot) <- parameters

for (i_param in seq_along(parameters)) {
  current_param <- parameters[i_param]
  meta_in_use <- meta %>%
    filter(parameter == current_param) %>%
    droplevels()

  y_label <- str_to_title(current_param)

  dbRDA <- dbrdalist[[current_param]]

  dbRDA_summary <- summary(dbRDA)
  dbRDA1_explainedvar_fitted <- round(dbRDA_summary$concont$importance[2, "dbRDA1"] * 100, 1)
  dbRDA2_explainedvar_fitted <- round(dbRDA_summary$concont$importance[2, "dbRDA2"] * 100, 1)
  dbRDA1_explainedvar_total <- round(dbRDA_summary$cont$importance[2, "dbRDA1"] * 100, 1)
  dbRDA2_explainedvar_total <- round(dbRDA_summary$cont$importance[2, "dbRDA2"] * 100, 1)

  dbRDA_sites <- as.data.frame(scores(dbRDA)$sites)
  dbRDA_env <- as.data.frame(scores(dbRDA)$biplot) %>%
    rownames_to_column()

  if (current_param == "algae") {
    dbRDA_env$rowname <- c("Sediment", "NCP", "NH4-N", "SRP")
  } else if (current_param == "diatoms") {
    dbRDA_env$rowname <- c("Transport", "NH4-N", "SRP")
  } else if (current_param == "bacteria") {
    dbRDA_env$rowname <- c("Transport", "Sediment", "SRP", "Cyanobacteria", "Diatoms")
  } else if (current_param == "fungi") {
    dbRDA_env$rowname <- c("Transport", "DOC")
  }

  # NMDS scores for samples ("sites")
  dbRDA_sites$id <- rownames(dbRDA_sites) # create a column of site names, from the rownames of data.scores
  dbRDA_sites <- merge(dbRDA_sites, meta_in_use, by = "id")
  dbRDA_sites$transport <- factor(dbRDA_sites$transport) # , levels = c("migrating-superficial", "migrating-underlying", "stationary-superficial", "stationary-underlying"))
  dbRDA_sites$sediment <- factor(dbRDA_sites$sediment) # , levels = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree"))
  # dbRDA_sites$treatment <- paste(dbRDA_sites$transport, dbRDA_sites$sediment, sep = "-")
  dbRDA_sites$treatment <- paste(dbRDA_sites$sediment, dbRDA_sites$transport, sep = "-")
  dbRDA_sites$treatment <- factor(dbRDA_sites$treatment, levels = levels_treatment_no_init)
  dbRDA_sites

  range(dbRDA_sites$dbRDA1)
  range(dbRDA_sites$dbRDA2)

  if (current_param == "algae") {
    x_lims <- c(-0.9, 1.6)
    y_lims <- c(-1.2, 1.1)
  } else if (current_param == "diatoms") {
    x_lims <- c(-1.7, 0.9)
    y_lims <- c(-1.3, 1.3)
  } else if (current_param == "bacteria") {
    x_lims <- c(-1.0, 1.0)
    y_lims <- c(-1.1, 0.9)
  } else if (current_param == "fungi") {
    x_lims <- c(-1.2, 1.1)
    y_lims <- c(-1.3, 2.0)
  }

  dbRDA_plot[[i_param]] <- ggplot(dbRDA_sites, aes(x = dbRDA1, y = dbRDA2)) +
    geom_point(data = dbRDA_sites, aes(x = dbRDA1, y = dbRDA2, colour = treatment, fill = treatment), size = 3.5, stroke = 1.5, shape = 24) + # add the point markers
    # geom_text(data = dbRDA_sites, aes(x = dbRDA1, y = dbRDA2, label = id), size = 5) + # add the point markers
    scale_y_continuous(breaks = seq(-5, 5, by = 0.5), limits = y_lims) +
    scale_x_continuous(breaks = seq(-5, 5, by = 0.5), limits = x_lims) +
    scale_colour_manual(values = myfill_no_init) +
    scale_fill_manual(values = myfill_no_init) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    labs(
      title = NULL, fill = NULL, shape = NULL,
      x = paste("dbRDA1\n", paste0("(", dbRDA1_explainedvar_fitted, "% of fitted, ", dbRDA1_explainedvar_total, "% of total variation)")),
      y = paste("dbRDA2\n", paste0("(", dbRDA2_explainedvar_fitted, "% of fitted, ", dbRDA2_explainedvar_total, "% of total variation)"))
    ) +
    guides(
      fill = guide_legend(override.aes = list(fill = myfill_no_init, colour = myfill_no_init, shape = 24)),
      colour = "none"
    ) +
    # annotation_compass(fontsize = 16, label = paste0("Stress = ", round(com_nmds_stress, 3), ", ", "k = ", com_nmds_ndim), position = "SW") +
    # annotation_compass(fontsize = 18, label = "p = 0.15", position = "SE") +
    annotation_compass(label = y_label, position = "NW", fontsize = 20) +
    geom_segment(
      data = dbRDA_env, aes(x = 0, xend = dbRDA1, y = 0, yend = dbRDA2),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = "black", show.legend = F
    ) +
    geom_text_repel(data = dbRDA_env, seed = 1, point.padding = NA, nudge_x = 0.05, nudge_y = 0.1, aes(x = dbRDA1, y = dbRDA2, label = rowname), colour = "black", size = 6, show.legend = F) +
    coord_fixed() +
    mytheme(no_space = TRUE) +
    theme(legend.key = element_rect(colour = "white", fill = NA)) +
    theme(aspect.ratio = 1)
  dbRDA_plot[[i_param]]

  if (i_param == 1) {
    dbRDA_plot_legend <- get_legend(dbRDA_plot[[i_param]] + theme(
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0), # create some space to the left of the legend
      legend.title = element_blank()
    ))
    filepath_legend_dbrda <- here(dir_out, paste0(filepattern_dbRDA_legend, format))
    ggsave(filepath_legend_dbrda, dbRDA_plot_legend, width = 10, height = 10, units = "cm", dpi = dpi)
  }

  dbRDA_plot[[i_param]] <- dbRDA_plot[[i_param]] +
    theme(legend.position = "none")
}

# read in legend as image
legend_dbrda <- image_read(filepath_legend_dbrda) %>%
  image_trim()

# compose plot
# dbRDA_plot_composed <- ((dbRDA_plot[["bacteria"]] + theme_no_space) + (dbRDA_plot[["fungi"]] + theme(plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "cm")))) / ((dbRDA_plot[["algae"]] + theme_no_space) + (dbRDA_plot[["diatoms"]] + theme(plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "cm")))) #& theme_no_space
dbRDA_plot_composed <- (dbRDA_plot[["bacteria"]] + dbRDA_plot[["fungi"]]) / (dbRDA_plot[["algae"]] + dbRDA_plot[["diatoms"]]) & theme_no_space

# tag plot
dbRDA_plot_composed <- dbRDA_plot_composed +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# add legend as image on plot
dbRDA_plot_composed <- ggdraw() +
  draw_plot(dbRDA_plot_composed) +
  draw_image(legend_dbrda, x = 0.82, y = 0.6, height = 0.6, width = 0.2, scale = 0.65)

#### (6) EFFECT SIZE ####

dt_effsize <- rbind(effsize, metab_effsize) %>%
  mutate(variable = factor(variable)) %>%
  mutate(group = case_when(
    variable %in% str_subset(levels(variable), "CR|GPP|NCP") ~ "metabolism",
    variable %in% str_subset(levels(variable), "Richness|Shannon") ~ "alpha-diversity",
    variable %in% str_subset(levels(variable), "Abundance") ~ "abundance",
    TRUE ~ NA_character_
  )) %>%
  mutate(group = factor(group)) %>%
  mutate(., contrast = fct_recode(contrast, "transport" = "Migrating - Stationary"))

# select data to plot
plotdat <- dt_effsize %>%
  filter(!grepl("Shannon", variable)) %>%
  filter(!grepl("Initial", contrast)) %>%
  filter(!grepl("GPP", variable))

# only transport
y_reordered <- plotdat %>%
  select(variable, group) %>%
  pull(variable) %>%
  unique() %>%
  as.character()

abundance_labs <- y_reordered[grepl("Abun", y_reordered)]
alpha_labs <- y_reordered[grepl("Richn", y_reordered)]
metab_labs <- y_reordered[grepl("CR|NCP", y_reordered)]

y_reordered <- c(alpha_labs, abundance_labs, metab_labs)
# y_reordered <- c(y_reordered[1:11], "OM", y_reordered[12:13], "BA", y_reordered[14:length(y_reordered)])

plotdat$variable <- factor(plotdat$variable, levels = y_reordered)
plotdat$group <- factor(plotdat$group, levels = levels(plotdat$group)[c(2, 1, 3)])
gr <- plotdat %>%
  group_by(group) %>%
  summarise(n = length(unique(variable))) %>%
  mutate(nn = cumsum(n)) %>%
  mutate(group = str_to_title(group))
gr$group <- recode_factor(gr$group, "Alpha-Diversity" = paste("alpha", "Diversity", sep = "-"))
# gr$group <- factor(gr$group, )
plotdat$dummy <- as.factor(paste(plotdat$contrast, plotdat$signif, sep = "-"))
plotdat$dummy <- factor(plotdat$dummy, levels = rev(levels(plotdat$dummy)))

eff_plot <- ggplot(plotdat) +
  geom_errorbarh(height = 1 / 4, aes(y = variable, xmin = CI_lower, xmax = CI_upper, colour = contrast, group = contrast), position = position_dodge(width = 0.8), show.legend = FALSE) +
  geom_point(aes(y = variable, x = hedges_g, fill = dummy, colour = contrast, group = contrast), size = 4, shape = 21, position = position_dodge(width = 0.8)) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.5) +
  geom_hline(yintercept = gr$nn[gr$group != "Metabolism"] + 0.5, linetype = "longdash", size = 0.5) +
  geom_text(data = gr, aes(x = 3, y = nn - 0.2, label = group), parse = TRUE, size = 5, hjust = 1) +
  scale_y_discrete(labels = c(
    "Algae", "Diatoms", "Bacteria", "Fungi",
    "Algae", "Diatoms", "Bacteria", "Fungi",
    "Initial CR", "Final CR", "Dynamics CR",
    "Initial NCP", "Final NCP", "Dynamics NCP"
  )) +
  # scale_fill_manual(values = c("white", "black"), labels = c("No", "Yes"), guide = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(
    values = c(
      "transport-TRUE" = myfill_no_init[1], "transport-FALSE" = "white",
      "sediment-TRUE" = myfill_no_init[3], "sediment-FALSE" = "white"
    ),
    breaks = c("transport-TRUE", "transport-FALSE"), labels = c("Yes", "No"), guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = c("black", "white")))
  ) +
  scale_colour_manual(
    values = c(myfill_no_init[1], myfill_no_init[3]),
    labels = c("Transport", "Sediment"),
    guide = guide_legend(order = 1, override.aes = list(shape = 21, fill = c(myfill_no_init[1], myfill_no_init[3])))
  ) +
  labs(x = "Hedges' g effect size", y = NULL, fill = "Significant", colour = "Factor") +
  mytheme(no_space = TRUE) +
  theme(legend.key = element_rect(colour = "white", fill = NA)) +
  # theme(legend.position = "bottom")
  theme(
    legend.position = c(0.82, 0.40),
    legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.spacing.x = unit(0, "pt")
  )
eff_plot

#### SAVE OUTPUT ####

# Fig. 1
ggsave(here(dir_out, paste0(filepattern_metab_plot, format)), metab_plot, width = 30, height = 22, units = "cm", dpi = dpi)

# Fig. 2
ggsave(here(dir_out, paste0(filepattern_abundance_plot, format)), abun_plot_composed, width = 20, height = 50, units = "cm", dpi = dpi)
ggsave(here(dir_out, paste0(filepattern_abundance_relab_plot, format)), abun_relab_plot_composed, width = 45, height = 55, units = "cm", dpi = dpi)

# Fig. 3
# ggsave(here(dir_out, paste0(filepattern_richness_plot, format)), richness_plot_composed, width = 20, height = 50, units = "cm", dpi = dpi)
# ggsave(here(dir_out, paste0(filepattern_venn_plot, format)), venn_plot2_composed, width = 20, height = 50, units = "cm", dpi = dpi)
# ggsave(here(dir_out, paste0(filepattern_richness_venn_plot, format)), richness_venn_plot_composed, width = 40, height = 55, units = "cm", dpi = dpi)
ggsave(here(dir_out, paste0(filepattern_richness_plot, format)), richness_plot_composed, width = 22, height = 37, units = "cm", dpi = dpi)
ggsave(here(dir_out, paste0(filepattern_venn_plot, format)), venn_plot2_composed, width = 20, height = 37, units = "cm", dpi = dpi)
ggsave(here(dir_out, paste0(filepattern_richness_venn_plot, format)), richness_venn_plot_composed, width = 40, height = 41, units = "cm", dpi = dpi)
ggsave(here(dir_out, paste0(filepattern_richness_upset_plot, format)), richness_upset_plot_composed, width = 50, height = 55, units = "cm", dpi = dpi)

# Fig. 4
ggsave(here(dir_out, paste0(filepattern_betadiversity_plot, format)), beta_plot_composed, width = 30, height = 30, units = "cm", dpi = dpi)

# Fig. 5
ggsave(here(dir_out, paste0(filepattern_effectsize_plot, format)), eff_plot, width = 20, height = 13, units = "cm", dpi = dpi)

# Fig. Sxx Metabolism dyanmics
ggsave(here(dir_out, paste0(filepattern_metab_plot, "_Dynamics", format)), metab_plot_dyn, width = 30, height = 15, units = "cm", dpi = dpi)

# Fig. Sxx Cyanobacteria abundance
ggsave(here(dir_out, paste0(filepattern_abundance_plot, "_Cyanbacteria", format)), abun_plot_cyanobacteria, width = 20, height = 15, units = "cm", dpi = dpi)

# Fig. Sxx dbRDA
ggsave(here(dir_out, paste0(filepattern_dbRDA_plot, format)), dbRDA_plot_composed, width = 30, height = 30, units = "cm", dpi = dpi)

# Fig. Sxx Shannon index
ggsave(here(dir_out, paste0(filepattern_shannon_plot, format)), shannon_plot_composed, width = 20, height = 37, units = "cm", dpi = dpi)
