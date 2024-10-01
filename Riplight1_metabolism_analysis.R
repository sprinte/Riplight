#### RIPLIGHT1 Analysis of metabolism data ####
#### Last updated: May 10, 2024

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
  "broom",
  "car",
  "cowplot",
  "data.table",
  "dplyr",
  "effectsize",
  "emmeans",
  "ggplot2",
  "ggpubr",
  "grid",
  "gridExtra",
  "lattice",
  "lme4",
  "lmerTest",
  "lmPerm",
  "patternplot",
  "png",
  "RColorBrewer",
  "tidyr"
)

using(packages)

#### USER SETTINGS ####

# determine in- and output directories
dir_in <- paste("D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Results/Metabolism/Tables") # where is all your data
dir_out <- paste("D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Results/Metabolism") # where the results and plots are saved
dir_stats <- "Stats"
dir_table <- "Tables"

# file name 
file <- "Metabolism_alldays_single_RipLight1_noO2outliers.txt"

# set variables
project <- "RipLight1" # to name output files
dates <- c("25.11.2020", "27.11.2020", "29.11.2020", "01.12.2020", "03.12.2020", "05.12.2020") # measurement dates

levels_sediment <- c("aqua" = "Aquatic", "mix" = "Mixed", "blind" = "Blind") # how factor 1 (sediment type) shall be renamed
levels_transport <- c("Roller" = "Migrating", "Wipper" = "Stationary") # how factor 2 (transport regime) shall be renamed

blind_correction <- TRUE # if TRUE, netrate (rate - blindrate) will be analyzed (-> netrate_dw); if FALSE, raw rates will be analyzed (-> rate_dw)
r2.min <- NA # minimum R2 for linear regression (rates) -> if NA or NULL, no restrictions apply. Rates below r2.min will be removed from analysis
implausibles <- 0 # can be 0 or NA; determines handling with implausibles values (positive CR, negative GPP -> put 0 or NA?)
outlier.thresh <- 1.5 # for calculated CR, NCP, and GPP rates: threshold of outlier identification: t*IQR will be subtracted (=lower boundary) from Q1 or added (=upper boundary) from Q3 to identify data points that lie beyond. They will be removed! The higher this factor, the less data will be removed (e.g. 3 = strong outlier). The smaller the number, the more points will be removed (e.g. 1.5 = weak outlier)
remove_outliers <- TRUE # if TRUE, outliers specified with the interquartile range rule (see outlier.thresh) will be removed from the dataset (only O2 values from raw data, no rates will be removed!)
days_as_initial <- c(1)
days_as_final <- c(12)

#### END OF SETTINGS ###

#### FUNCTIONS ####

se <- function(x) sqrt(var(x, na.rm = T) / length(x))

#### DATA IMPORT ####

setwd(dir_in)
activity <- read.table(file, header = TRUE, sep = "\t", dec = ".", stringsAsFactors = TRUE)

# order columns
activity <- activity[, c(1, 4, 2, 3, 5, 14, 13, 12, 6:11)]

# renaming factors and factor levels
# factor 1: sediment type (aquatic, mixed, blind)
activity <- rename(activity, sediment = treatment)
activity$sediment <- dplyr::recode(activity$sediment, !!!levels_sediment)
activity$sediment <- factor(activity$sediment, levels = levels_sediment[c(2, 1, 3)]) # correct order of factor levels

# factor 2: transport regime (Roller=migrating, Wipper=stationary)
activity <- rename(activity, transport = device)
activity$transport <- dplyr::recode(activity$transport, !!!levels_transport)
activity$transport <- factor(activity$transport, levels = levels_transport[c(1, 2)]) # correct order of factor levels

# correct order of time variable
activity$logclass <- factor(activity$logclass, levels = levels(activity$logclass)[c(1, 5:12, 2:4)]) # correct order of days

# remove data below r2.min
if (!is.na(r2.min)) {
  percent_removed <- activity %>%
    filter(variable %in% c("CR", "NCP") & cycle > 2) %>%
    mutate(flag = ifelse(r2 >= r2.min, 1, 0)) %>%
    group_by(variable) %>%
    summarize(
      rel_keep = round(sum(flag == 1, na.rm = T) / n() * 100, 1),
      rel_drop = round(sum(flag == 0, na.rm = T) / n() * 100, 1),
      rel_na = round(sum(is.na(flag)) / n() * 100, 1)
    ) %>%
    as.data.frame()

  activity <- activity[activity$r2 >= r2.min, ]
}

# create new column (ID) to identify replicates: combine the replicate number with the factor level names
# m = migrating, s = stationary, a = aquatic, m = mixed
# activity$ID<-paste(tolower(substr(activity$transport, 1, 1)), tolower(substr(activity$sediment, 1, 1)), activity$replicate, sep="")

#### DATA PREPARATION ####
# remove blind values
metab <- activity
metab <- metab[metab$sediment != levels_sediment[3], ] # without blind
metab <- droplevels(metab)

# remove single replicates (CAUTION: might be biased!)
# metab<-metab[metab$ID!="sa3", ] # remove stationary aquatic #3
# metab<-metab[!(metab$variable=="CR"&metab$ID=="mm5"), ] # remove migrating mixed #5 (CR only)

# create new column (resp) copying the rate which will be analyzed from here on (most likely, netrate_dw (with blind correction) or rate_dw (without blind correction))
if (blind_correction == TRUE) {
  metab$resp <- metab$netrate_dw
} else {
  metab$resp <- metab$rate_dw
}

# # remove implausible CR/GPP (CR must be <0, GPP must be >0, NCP can be either one)
# metab <- metab %>%
#   mutate(resp = ifelse(variable == "CR" & resp > 0, implausibles, resp)) %>% # sets positive CR to value set at "implausibles" (caused by high blank respiration!)
#   mutate(resp = ifelse(variable == "GPP" & resp < 0, implausibles, resp)) %>% # sets negative GPP to value set at "implausibles"
#   as.data.frame()

# select only relevant columns
metab <- metab %>%
  select(c(sediment, transport, replicate, cycle, logclass, variable, totaltimemean, resp)) %>%
  as.data.frame()

# remove outliers (if desired)
if (remove_outliers == TRUE) {
  # identify outliers (per group):
  metab_out <- metab %>%
    group_by(sediment, transport, cycle, variable) %>%
    mutate(
      Q1 = unname(summary(resp))[2],
      Q3 = unname(summary(resp))[5],
      iqr = IQR(resp),
      lower = Q1 - outlier.thresh * iqr,
      upper = Q3 + outlier.thresh * iqr
    ) %>%
    ungroup() %>%
    mutate(flag = ifelse(resp < lower | resp > upper, 0, 1)) %>%
    filter(flag == 0) %>%
    select(-c(Q1, Q3, iqr, lower, upper, flag)) %>%
    as.data.frame()
  nrow(metab_out)
  nrow(metab)
  nrow(metab_out) / nrow(metab) * 100

  # complete data for all cases (unfortunately needed to correctly position the points on the plot)
  all <- metab %>% tidyr::expand(nesting(sediment, transport, replicate), nesting(cycle, logclass), variable)
  metab_out <- metab_out %>% right_join(all)

  metab <- metab %>%
    anti_join(metab_out[which(!is.na(metab_out$resp)), ], by = c("sediment", "transport", "replicate", "cycle", "logclass", "variable", "totaltimemean")) %>%
    as.data.frame()
}

# remove implausible CR/GPP (CR must be <0, GPP must be >0, NCP can be either one)
metab <- metab %>%
  mutate(resp = ifelse(variable == "CR" & resp > 0, implausibles, resp)) %>% # sets positive CR to value set at "implausibles" (caused by high blank respiration!)
  mutate(resp = ifelse(variable == "GPP" & resp < 0, implausibles, resp)) %>% # sets negative GPP to value set at "implausibles"
  as.data.frame()

# calculate means of replicates
meanmetab <- metab %>%
  group_by(sediment, transport, variable, logclass, totaltimemean) %>%
  summarize(mean_resp = mean(resp, na.rm = T), sd_resp = sd(resp, na.rm = T), se_resp = se(resp), replicates = n()) %>%
  as.data.frame()

## Subset data for relevant data: initial, final, and dynamics
# initial activity
metab.init <- metab %>%
  filter(cycle %in% days_as_initial) %>%
  filter(!is.na(resp)) %>%
  mutate(logclass = factor("initial")) %>%
  select(c(sediment, transport, replicate, variable, logclass, resp)) %>%
  as.data.frame()

# final activity
metab.fin <- metab %>%
  filter(cycle %in% days_as_final) %>%
  filter(!is.na(resp)) %>%
  mutate(logclass = factor("final")) %>%
  select(c(sediment, transport, replicate, variable, logclass, resp)) %>%
  as.data.frame()

# dynamics (dO2/dt between start and finish)
metab.regs <- metab %>%
  group_by(sediment, transport, variable, replicate) %>%
  filter(!is.na(resp)) %>%
  do(respchange = tidy(lm(resp ~ totaltimemean, data = .))) %>%
  unnest(respchange) %>%
  rename(resp = estimate)

metab.regs.r2 <- metab %>%
  group_by(sediment, transport, variable, replicate) %>%
  filter(!is.na(resp)) %>%
  do(respchange = glance(lm(resp ~ totaltimemean, data = .))) %>%
  unnest(respchange) %>%
  select(sediment, transport, variable, replicate, adj.r.squared)

metab.regs <- metab.regs %>%
  left_join(metab.regs.r2, by = c("sediment", "transport", "variable", "replicate"))

# calculate fitted rates (linear model)
metab.fitted <- metab %>%
  group_by(sediment, transport, variable, replicate) %>%
  filter(!is.na(resp)) %>%
  do(respchange = augment(lm(resp ~ totaltimemean, data = .))) %>%
  unnest(respchange) %>%
  select(sediment, transport, replicate, variable, totaltimemean, .fitted)

metab.fitted <- metab %>%
  left_join(metab.fitted, by = c("sediment", "transport", "variable", "replicate", "totaltimemean")) %>%
  rename(fitted = .fitted)

# remove outliers (if desired)
if (remove_outliers == TRUE) {
  # identify outliers (per group):
  metab_regs_out <- metab.regs %>%
    group_by(sediment, transport, term, variable) %>%
    mutate(
      Q1 = unname(summary(resp))[2],
      Q3 = unname(summary(resp))[5],
      iqr = IQR(resp),
      lower = Q1 - outlier.thresh * iqr,
      upper = Q3 + outlier.thresh * iqr
    ) %>%
    ungroup() %>%
    mutate(flag = ifelse(resp < lower | resp > upper, 0, 1)) %>%
    filter(flag == 0) %>%
    select(-c(Q1, Q3, iqr, lower, upper, flag)) %>%
    as.data.frame()

  # complete data for all cases (unfortunately needed to correctly position the points on the plot)
  all <- metab.regs %>% tidyr::expand(nesting(sediment, transport, replicate), term, variable)
  metab_regs_out <- metab_regs_out %>% right_join(all)

  metab.regs <- metab.regs %>%
    anti_join(metab_regs_out[which(!is.na(metab_regs_out$resp)), ], by = c("sediment", "transport", "replicate", "variable", "term")) %>%
    as.data.frame()
}

metab.dyn <- metab.regs %>%
  filter(term == "totaltimemean") %>%
  mutate(logclass = factor("dynamics")) %>%
  select(c(sediment, transport, replicate, variable, logclass, resp)) %>%
  as.data.frame()

# combine to single data.frame
metabstat <- rbind.data.frame(metab.init, metab.fin, metab.dyn, stringsAsFactors = TRUE)
# metabstat$logclass<-factor(metabstat$logclass, levels(metabstat$logclass)[c(3,2,1,4)])

meanmetabstat <- metabstat %>%
  group_by(sediment, transport, variable, logclass) %>%
  filter(!is.na(resp)) %>%
  summarise(mean_resp = mean(resp, na.rm = T), sd_resp = sd(resp, na.rm = T), se_resp = se(resp), replicates = n()) %>%
  arrange(sediment, transport, logclass, variable) %>%
  as.data.frame()

#### STATISTICAL ANALYSES ####
## pANOVA for initial respiration, final respiration, and dynamics

aovdat <- metabstat %>%
  mutate(resp = ifelse(variable == "CR", -resp, resp))
head(aovdat)

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
MigMix <- c(1, 0, 0, 0)
MigAqu <- c(0, 1, 0, 0)
StatMix <- c(0, 0, 1, 0)
StatAqu <- c(0, 0, 0, 1)

custom_contrasts <- list(
  "Mixed Migrating - Mixed Stationary" = MigMix - StatMix,
  "Mixed Migrating - Aquatic Migrating" = MigMix - MigAqu,
  "Mixed Stationary - Aquatic Stationary" = StatMix - StatAqu,
  "Aquatic Migrating - Aquatic Stationary" = MigAqu - StatAqu
)

mylevs <- levels(aovdat$variable)

# calculate permANOVA and post-hoc test (selected pairwise comparisons) # here with aovp (lmPerm), can alternatively be done with aovperm (permuco)
for (m in mylevs) {
  aovdat_m <- aovdat[aovdat$variable == m, ] %>% droplevels()

  for (i in seq_along(levels(aovdat_m$logclass))) {
    j <- levels(aovdat_m$logclass)[i]

    aovdat_sub <- aovdat_m %>%
      filter(logclass == j) %>%
      droplevels() %>%
      mutate(treatment = paste(sediment, transport, sep = "-")) %>%
      mutate(treatment = as.factor(treatment))
    aovdatlist[[m]][[j]] <- aovdat_sub

    aovdat_sub_summary <- aovdat_sub %>%
      group_by(sediment, transport, treatment, variable, logclass) %>%
      filter(!is.na(resp)) %>%
      summarize(replicates = n()) %>%
      relocate(sediment, transport, treatment, variable, logclass, replicates) %>%
      as.data.table()

    levelcomb <- combn(aovdat_sub_summary$treatment, 2)
    levelcombnames <- paste(levelcomb[1, ], levelcomb[2, ], sep = " - ")
    replsum <- combn(aovdat_sub_summary$replicates, 2, FUN = sum)
    combinations <- data.frame(contrast = levelcombnames, replicates = replsum)

    # pANOVA
    set.seed(111)
    model_add <- aovp(resp ~ sediment + transport, aovdat_sub, perm = "Exact", seqs = F, maxIter = 9999)
    summary(model_add)
    set.seed(111)
    model_inter <- aovp(resp ~ sediment * transport, aovdat_sub, perm = "Exact", seqs = F, maxIter = 99999)
    summary(model_inter)
    emmip(model_inter, sediment ~ transport) # plot linear prediction

    AIC_model_inter <- AIC(model_inter)
    AIC_model_add <- AIC(model_add)
    if (abs(AIC_model_inter - AIC_model_add) < 2) {
      print("Best model: undecided")
      bestmodellist[[m]][[j]] <- "undecided"
    } else if (AIC_model_inter < AIC_model_add) {
      print("Best model: Interaction")
      bestmodellist[[m]][[j]] <- "interaction"
    } else if (AIC_model_inter > AIC_model_add) {
      print("Best model: Addition")
      bestmodellist[[m]][[j]] <- "addition"
    }

    df_inter <- as.data.frame(summary(model_inter)[[1]])
    df_inter$factor <- row.names(df_inter)
    df_inter <- df_inter %>%
      mutate(
        variable = m,
        logclass = j
      ) %>%
      relocate(variable, logclass, factor)
    row.names(df_inter) <- 1:nrow(df_inter)
    aovlist_inter[[m]][[j]] <- df_inter

    df_add <- as.data.frame(summary(model_add)[[1]])
    df_add$factor <- row.names(df_add)
    df_add <- df_add %>%
      mutate(
        variable = m,
        logclass = j
      ) %>%
      relocate(variable, logclass, factor)
    row.names(df_add) <- 1:nrow(df_add)
    aovlist_add[[m]][[j]] <- df_add

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

    if (all(aovdat_aquatic$resp == 0)) {
      hedges_aquatic <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_aquatic <- hedges_g(resp ~ transport, data = aovdat_aquatic)
    }

    if (all(aovdat_mixed$resp == 0)) {
      hedges_mixed <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_mixed <- hedges_g(resp ~ transport, data = aovdat_mixed)
    }

    if (all(aovdat_migrating$resp == 0)) {
      hedges_migrating <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_migrating <- hedges_g(resp ~ sediment, data = aovdat_migrating)
    }

    if (all(aovdat_stationary$resp == 0)) {
      hedges_stationary <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_stationary <- hedges_g(resp ~ sediment, data = aovdat_stationary)
    }

    df_hedges_inter <- data.frame(
      "contrast" = treatment_contrasts,
      "hedges_g" = c(
        hedges_mixed$Hedges_g,
        hedges_migrating$Hedges_g,
        hedges_stationary$Hedges_g,
        hedges_aquatic$Hedges_g
      ),
      "CI_lower" = c(
        hedges_mixed$CI_low,
        hedges_migrating$CI_low,
        hedges_stationary$CI_low,
        hedges_aquatic$CI_low
      ),
      "CI_upper" = c(
        hedges_mixed$CI_high,
        hedges_migrating$CI_high,
        hedges_stationary$CI_high,
        hedges_aquatic$CI_high
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

    if (all(aovdat_sub$resp == 0)) {
      hedges_transport <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_sediment <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_transport <- hedges_g(resp ~ transport, data = aovdat_sub)
      hedges_sediment <- hedges_g(resp ~ sediment, data = aovdat_sub)
    }

    df_hedges_transport <- data.frame(
      "contrast" = treatment_contrasts_transport,
      "hedges_g" = hedges_transport$Hedges_g,
      "CI_lower" = hedges_transport$CI_low,
      "CI_upper" = hedges_transport$CI_high
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
    emmlist_inter[[m]][[j]] <- emm.contrasts_inter
    hedgeslist_inter[[m]][[j]] <- df_hedges_inter
    emmlist_add[[m]][[j]] <- emm.contrasts_add
    hedgeslist_add[[m]][[j]] <- df_hedges_add

    # add variable
    emmlist_inter[[m]][[j]]$variable <- m
    hedgeslist_inter[[m]][[j]]$variable <- m
    emmlist_add[[m]][[j]]$variable <- m
    hedgeslist_add[[m]][[j]]$variable <- m

    # add logclass
    emmlist_inter[[m]][[j]]$logclass <- j
    hedgeslist_inter[[m]][[j]]$logclass <- j
    emmlist_add[[m]][[j]]$logclass <- j
    hedgeslist_add[[m]][[j]]$logclass <- j
  }
}

# combine to single dt
aovtable <- rbindlist(lapply(aovlist_add, rbindlist))
bestmodeltable <- unlist(bestmodellist)
dt_bestmodel <- data.table(parameter = names(bestmodeltable), best_model = bestmodeltable)

emm_inter <- rbindlist(lapply(emmlist_inter, rbindlist))
emm_add <- rbindlist(lapply(emmlist_add, rbindlist))

hedges_inter <- rbindlist(lapply(hedgeslist_inter, rbindlist)) %>%
  mutate(variable = factor(variable, levels = mylevs)) %>%
  mutate(logclass = factor(logclass)) %>%
  mutate(contrast = factor(contrast))
hedges_add <- rbindlist(lapply(hedgeslist_add, rbindlist)) %>%
  mutate(variable = factor(variable, levels = mylevs)) %>%
  mutate(logclass = factor(logclass)) %>%
  mutate(contrast = factor(contrast)) %>%
  mutate(contrast = recode_factor(contrast,
    "Migrating - Stationary" = "transport",
    "Mixed - Aquatic" = "sediment"
  ))

#### SAVING OUTPUT ####

metab_out <- metab %>%
  arrange(sediment, transport, logclass, variable) %>%
  complete(nesting(sediment, transport, replicate), variable, nesting(cycle, logclass)) %>%
  arrange(sediment, transport, logclass, variable) %>%
  as.data.frame()

metabstat_out <- metabstat %>%
  arrange(sediment, transport, logclass, variable) %>%
  complete(nesting(sediment, transport, replicate), variable, logclass) %>%
  arrange(sediment, transport, logclass, variable) %>%
  as.data.frame()

# metabolism
write.table(metab_out, file = here(dir_out, dir_table, paste0(paste("Metabolism_stats_single_all", project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
write.table(metabstat_out, file = here(dir_out, dir_table, paste0(paste("Metabolism_stats_single", project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
write.table(meanmetabstat, file = here(dir_out, dir_table, paste0(paste("Metabolism_stats_means", project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)

# pANOVA
write.table(aovtable, file = here(dir_out, dir_stats, paste0(paste("Metabolism_stats_pANOVA", project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
write.table(hedges_inter, file = here(dir_out, dir_stats, paste0(paste("Metabolism_stats_effsize_inter", project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
write.table(hedges_add, file = here(dir_out, dir_stats, paste0(paste("Metabolism_stats_effsize_add", project, sep = "_"), ".txt")), sep = "\t", row.names = FALSE, quote = F)
