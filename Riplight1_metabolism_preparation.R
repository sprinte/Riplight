#### RIPLIGHT1 Preparation of metabolism data ####
#### Last updated: May 10, 2024

#### IMPORTANT: DATA ORGANIZATION OUTSIDE R ####
#### If you operated several devices, each device needs its own folder with all necessary data inside (see list of necessary files below). If the measurements ran for several days, and oxy data were logged continuously (can be interrupted and started again, as long as the microcosms stay on their channels!), one folder for each device suffices. If channels were switched, or microcosms were alternated on different devices, then each measurement "cycle" (= microcosm stays on one channel) needs its own data folder.
#### I recommend that the folder names contain the device used and the date of measurement, e.g. "Roller_18.05.2020", "Wipper_18.05.2020"
#### Following files are needed in data directory ("..." may contain whatever name you gave it)
#
# "Calibration... .txt": [tab-separated, "\t"] Calibration file which contains all background data about the samples
#                       "device":       name of device
#                       "treatment":    name of sample treatment
#                       "replicate":    replicate number
#                       "glass_ID":     number of microcosm used
#                       "glass_no":     if a glass broke during the experiment and the sample had to be transferred to a new microcosm, enter a new line (row) with the same information but including the new microcosm number in the column "glass_ID". The first glass (the one that broke) will get a "1" in "glass_no", the new (replaced) glass will get a "2". If no microcosm broke, just put "1" everywhere (default)
#                       "break_time":   ["dd.mm.YYYY HH:MM"] If a microcosm broke, put the time when you replaced it with a new one. If no microcosm broke, just put "NA" everywhere (default)
#                       "density":      [g/cm?], [g/ml] (dry) density of your sediment samples, default of Spree sediment is 2.562 g/cm?
#                       "dw":           [g] dry weight (105 ?C, overnight) of your sediment samples
#                       "afdw":         [g] (OPTIONAL) ash-free dry weight (550 ?C, 5h) of your sediment samples
#                       "glassvolume":  [ml] volume of microcosm, standard of our microcosms is 48 ml
#                       "porevolume":   [ml] (calculated) volume of water in microcosm, to be calculated as the difference between total microcosm volume and volume of dry sediment: porevolume = glassvolume - (dw/density)
#                       "phi_T0":       (calibration data) (OPTIONAL) phase angle at T0 (T0 = temperature at 0% O2 saturation)
#                       "phi100_T100":  (calibration data) (OPTIONAL) phase angle at T100 (T100 = temperature at 100% O2 saturation)
#                       "T0":           (calibration data) (OPTIONAL) [?C] temperature at 0% O2 saturation
#                       "T100":         (calibration data) (OPTIONAL) [?C] temperature at 100% O2 saturation
#                       "p_cal":        (calibration data) (OPTIONAL) [hPa] atmospheric pressure during calibration
#
# "Lights... .txt": [tab-separated, "\t"] Lights file which contains information about the light/dark periods, and possible periods that should not be taken into account
#                       "date":         ["dd.mm.YYYY"] date of measurement
#                       "start":        ["dd.mm.YYYY HH:MM"] start time of period
#                       "end":          ["dd.mm.YYYY HH:MM"] end time of period
#                       "light":        for dark periods, put "off", for light periods, put "on"
#                       "mode":         if oxy values in this period shall be counted into your analysis, put "running". If not, put anything else.
#
# "... -ch1.txt", "... -ch2.txt" etc. (OLD DEVICE) OR "[filename]... .txt" (NEW DEVICE):
#       Oxygen data files created by the oxy device itself. (CAUTION: If measurement were stopped and started again, the Oxy devices generates new files. R will import all files and combine those that belong to the same channel. So R thinks that all files from Ch-xy contain the SAME SAMPLE over the entire meausurement period, regardless of how often the measurement was interrupted and restarted. This means, you CANNOT just swap microcosms and change their channels during the measurement! If you did so, you must import them manually.)
#       OLD DEVICE: -> creates one txt-file for EACH SINGLE channel [semicolon-separated, ";"] Line 1-37: just contains background stuff (e.g. the calibration used). This will be ignored by R. It is important to check that the actual data (i.e., its header) starts in line 38! Sometimes an additional line named "Recalibrated" or so is added - you must remove this prior to importing the file in R!
#                       "Date/dd:mm:yy":Date of measurement
#                       "Time/hh:mm:ss":Time of measurement
#                       "Logtime/h":    Hours since start of meausurement
#                       "Oxy/mg/l":     Measured O2 concentration
#                       "Phase/?":      Raw signal of measurement (phase angle) - important if you need to recalibrate the data afterwards, e.g. because the calibration data during measuring were wrong!
#                       "Amp":          Amplitude of measuring signal
#                       "Temp/?C":      Temperature during measurement - careful, the old device usually does not log temperature, so you need an external temperature logger (-> "Temperature... .txt" file)
#       NEW DEVICE: -> creates one txt-file for ALL channels [tab-separated, "\t"]
#                       "Date":         ["mm/dd/YYYY"] Date of measurement
#                       "Time":         ["mm/dd/YYYY"] Date of measurement
#                       "Channel":      Number of oxygen channel
#                       "User":         Name of user (usually "default")
#                       "SensorID":     quite long number of sensor in use
#                       "Sensor_Name":  Name which you assigned to the sensor/channel, it is necessary (!!) to include the microcosm number somewhere in the name, e.g. "flask20", "glass20", "fl20", "20" - does not really matter how, as long as the number is in there!
#                       "delta_t":      Hours since start of meausurement
#                       "Time_Unit":    Unit of "delta_t"
#                       "Value":        Measured O2 concentration
#                       "O2_Unit":      Unit of "Value"
#                       "Mode":         Mode of meausurement, in water should be "Humid"
#                       "Phase":        Raw signal of measurement (phase angle) - important if you need to recalibrate the data afterwards, e.g. because the calibration data during measuring were wrong!
#                       "Phase_Unit":   Unit of "Phase"
#                       "Amplitude":    Raw signal of measurement (phase angle) - important if you need to recalibrate the data afterwards, e.g. because the calibration data during measuring were wrong!
#                       "Amplitude_Unit":Unit of "Amplitude"
#                       "Temp":         Temperature during measurement - the new device usually logs the temperature, so you do not need an external temperature logger (-> "Temperature... .txt" file)
#                       "Temp_Unit":    Unit of "Temp"
#                       "Pressure":     Atmospheric during measurement - the new device usually logs the pressure, so you do not need the DWD pressure data (-> "Pressure... .txt" file)
#                       "Pressure_Unit":Unit of "Pressure"
#                       "Cal0":         phase angle at T0 (T0 = temperature at 0% O2 saturation)
#                       "Cal0_Unit":    Unit of "Cal0"
#                       "T0":           temperature at 0% O2 saturation
#                       "T0_Unit":      Unit of "T0"
#                       "O2Cal2nd":     O2 concentration at second calibration point (2-point calibration)
#                       "O2_Unit":      Unit of "O2Cal2nd"
#                       "Cal2nd":       phase angle at second calibration point with T2nd
#                       "Cal2nd_Unit":  Unit of "Cal2nd"
#                       "T2nd":         temperature at second calibration point
#                       "T2nd_Unit":    Unit of "T2nd"
#                       "CalPressure":  atmospheric pressure during calibration
#                       "CalPressure_Unit":Unit of "CalPressure"
#             [the remaining columns are additional meta data that we will not use for our calculations, so I will not explain them in detail here]
#
# optional: "Pressure... .txt": [semicolon-separated, ";"] Atmospheric pressure during measurement (e.g., recorded hourly) retrieved from DWD Climate Data Center (ftp://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/pressure/recent/). Not required when your Oxy device is able to log the pressure.
#                       "STATIONS_ID":  weather station of recording (Cottbus: 880, for Bad Saarow use Lindenberg: 3015)
#                       "MESS_DATUM":   ["YYYYmmddHH"] time of recording
#                       "QN_8":         quality index of recorded data
#                       "P":            [hPa] atmospheric pressure at sea level
#                       "P0":           [hPa] atmospheric pressure at weather station level
#                       "eor":          "end of record"
#
# optional: "Temperature... .txt": [tab-separated, "\t"] Temperature file which contains data logged by an external temperature logger. Careful: comes usually in a weird format with a lot of quotes (" "). It is best to open the file first in Excel, clean the file (remove line 2 that contains "?C" and save again to a new .txt-file)
#                       "datetime":     ["dd.mm.YYYY HH:MM" or "dd/mm/YYYY HH:MM"] Timestamp of measurement
#                       "ch1", "ch2",...:[?C] measured temperature values (if more than one column, they will be averaged
#
# "oxyfun.R": contains all functions required to import oxy data and calculate respiration rates
#
# this script here (whatever you named it)

#### USER SETTINGS ####
source("D:/Documents/Studium/Promotion/Skripte/oxyfun.R") # enter the path with the oxyfun script

# determine in- and output directories
path_in <- paste("D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Data/Metabolism") # where is all your data
path_out <- paste("D:/Documents/Studium/Promotion/EcoMigRip/RipLight1/Results/Metabolism") # where the results and plots are saved

# project settings
project <- "RipLight1"
stream <- "Spree" # name of stream
dates <- c("25.11.2020", "27.11.2020", "29.11.2020", "01.12.2020", "03.12.2020", "05.12.2020") # measurement dates
treatments <- c("aqua", "mix") # all treatments used in the experiment, must be named as they appear in the calibration sheet!
devices <- c("Roller", "Wipper") # names of devices used
newoxy_name <- "Oxygen_Wipper" # what is the name of the oxy data in case you used the new device (does not have to be the complete name, it will look for a file that contains this part)
blank_name <- "blind" # if you used a blind that might have a respiration rate, and you want to subtract this from the rest of the rates, put the name of your blind. If you did not use a blind, put NA
myorder <- c(4, 5, 6, 1, 2, 3) # you need to tell R the correct (=chronological) order of the folders, as R will list them alphabetically!

# set parameters
r2.min <- NA # minimum R2 for linear regression (rates) -> if NA or NULL, no restrictions apply
anox <- 2 # anoxic threshold
i.lag <- 6 # incubation lag
l.lag <- 2 # light lag
d.lag <- 2 # dark lag
e.lag <- 0 # water exchange lag
tmax <- Inf # max. running time for experiment

q10 <- TRUE # if TRUE, Q10 temperature adjustion will be made with a factor of 2 (see Winkler et al. 1997)
blind_correction <- TRUE # if TRUE, blind rates will be subtracted from raw rates (-> netrate), if FALSE, raw rates will be used (-> rate)
outlier.thresh_O2 <- 3 # for raw O2 data: threshold t of outlier identification: t*IQR will be subtracted (=lower boundary) from Q1 or added (=upper boundary) from Q3 to identify data points that lie beyond. They will be removed! The higher this factor, the less data will be removed (e.g. 3 = strong outlier). The smaller the number, the more points will be removed (e.g. 1.5 = weak outlier)
outlier.thresh_rate <- 1.5 # for calculated CR, NCP, and GPP rates
remove_outliers <- TRUE # if TRUE, outliers specified with the interquartile range rule (see outlier.thresh) will be removed from the dataset (only O2 values from raw data, no rates will be removed!)

# graphical settings
levels_sediment <- c("aqua" = "Aquatic", "mix" = "Mixed", "blind" = "Blind") # how factor 1 (sediment type) shall be renamed
levels_transport <- c("Roller" = "Migrating", "Wipper" = "Stationary") # how factor 2 (transport regime) shall be renamed

#### END OF SETTINGS ###

#### FUNCTIONS ####
give.n <- function(y) {
  return(data.frame(y = ifelse(quantile(y, 0.75) >= 0, quantile(y, 0.75) * 1.1, quantile(y, 0.75) * 0.9), label = paste("n =", length(y[which(!is.na(y))]))))
}

#### DATA IMPORT ####

setwd(path_in)

# folders with data
folders_r <- list.files(path = path_in, pattern = devices[1], recursive = F, include.dirs = T)[file.info(list.files(pattern = devices[1], recursive = F, include.dirs = T))$isdir == T] # lists all folders from first device
folders_w <- list.files(path = path_in, pattern = devices[2], recursive = F, include.dirs = T)[file.info(list.files(pattern = devices[2], recursive = F, include.dirs = T))$isdir == T] # lists all folders from second device

folders_r <- folders_r[myorder]
folders_w <- folders_w[myorder]

## Oxydata from Roller ##
oxyfiles_r <- vector("list", length = length(folders_r))
names(oxyfiles_r) <- dates
for (i in seq_along(folders_r)) {
  oxyfiles_r[[i]] <- oxyread.old(dir = paste0("./", folders_r[i]), blind = NA, read.patm = TRUE, read.temp = TRUE, mode = "auto", temp_datetimeformat = "%d.%m.%Y %H:%M", oxy_dateformat = "%d.%m.%Y", oxy_timeformat = "%H:%M:%S")
  oxyfiles_r[[i]]$date <- dates[i]
  oxyfiles_r[[i]]$day <- i
}
all_oxyfiles_r <- as.data.frame(do.call(rbind, oxyfiles_r), stringsAsFactors = FALSE)
rownames(all_oxyfiles_r) <- 1:nrow(all_oxyfiles_r)

## Oxydata from Wipper ##
oxyfiles_w <- vector("list", length = length(folders_w))
names(oxyfiles_w) <- dates
for (i in seq_along(folders_w)) {
  oxyfiles_w[[i]] <- oxyread.new(dir = paste0("./", folders_w[i]), filename = newoxy_name, blind = blank_name, read.patm = FALSE, read.temp = FALSE, mode = "auto", oxy_dateformat = "%d.%m.%Y", oxy_timeformat = "%H:%M:%S")
  oxyfiles_w[[i]]$date <- dates[i]
  oxyfiles_w[[i]]$day <- i
}
all_oxyfiles_w <- as.data.frame(do.call(rbind, oxyfiles_w), stringsAsFactors = FALSE)
rownames(all_oxyfiles_w) <- 1:nrow(all_oxyfiles_w)

## COMBINE to single data.frame
oxydat <- rbind(all_oxyfiles_r, all_oxyfiles_w)
oxydat$oxy_mg_l <- as.numeric(oxydat$oxy)
oxydat$oxy100_mg_l <- oxysat(p_atm = oxydat$p_atm, temperature = oxydat$temperature)

## Lights file ##
lights <- read.table(paste0("./", folders_r[1], "/", list.files(paste0("./", folders_r[1]), pattern = "Lights", ignore.case = T, recursive = T)), sep = "\t", dec = ".", header = T, stringsAsFactors = TRUE)
lights$start <- as.POSIXct(lights$start, format = "%d.%m.%Y %H:%M")
lights$end <- as.POSIXct(lights$end, format = "%d.%m.%Y %H:%M")
lights$hours <- round(difftime(lights$end, lights$start, units = "hours"), 2)
lights$cycle <- NA # create a variable that counts the cycles of light and dark
for (i in seq_along(levels(lights$mode))) {
  j <- levels(lights$mode)[i]
  for (k in seq_along(levels(lights$light))) {
    l <- levels(lights$light)[k]
    for (m in 1:(nrow(lights[lights$mode == j & lights$light == l, ]))) {
      lights$cycle[lights$mode == j & lights$light == l][m] <- m
    }
  }
}

### PLOTTING O2 PROFILE ####
o2_profile <- oxyplot(oxydata = oxydat, x.lims = c(0, 48), lightdata = lights, x.breaks = 24, show.sat = TRUE, show.cycle = TRUE, anoxic = anox, incub.lag = i.lag, light.lag = l.lag, dark.lag = d.lag, exchange.lag = e.lag, runtime = tmax, show.outlier = TRUE, threshold.outlier = outlier.thresh_O2)
o2_profile

# O2 plots for single days
o2_singday <- vector("list", length = length(unique(oxydat$cycle)))
names(o2_singday) <- unique(oxydat$cycle)
for (i in unique(oxydat$cycle)) {
  lights_single <- lights[lights$cycle == i, ]
  oxydat_single <- oxydat[oxydat$cycle == i & oxydat$mode == "running", ]
  il <- ifelse(i == 1, i.lag, 0)
  o2_singday[[i]] <- oxyplot(oxydata = oxydat_single, lightdata = lights_single, x.lims = c(min(oxydat_single$totaltime)[1], max(oxydat_single$totaltime)[1]), x.breaks = 3, show.sat = TRUE, show.cycle = TRUE, anoxic = anox, incub.lag = il, light.lag = l.lag, dark.lag = d.lag, exchange.lag = e.lag, runtime = tmax)
}

# O2 plots for single replicates
o2_singrep<-vector("list", length = length(levels(oxydat$replicate))*length(levels(oxydat$treatment)))
names(o2_singrep)<-paste(rep(levels(oxydat$treatment), each=length(levels(oxydat$replicate))), levels(oxydat$replicate), sep="_")
o2_singrep<-o2_singrep[-c(13:15)] # there were only two blinds, therefore removing #3 and #4

for (i in seq_along(o2_singrep)) {
  repname<-names(o2_singrep[i])
  re<-substr(repname, nchar(repname), nchar(repname))
  tr<-substr(repname, 1, nchar(repname)-2)
  oxydat_single<-oxydat[oxydat$treatment==tr&oxydat$replicate==re,]
  #o2_singrep[[i]]<-oxyplot(oxydata=oxydat_single, lightdata=lights, x.lims=c(min(oxydat_single$totaltime)[1], max(oxydat_single$totaltime)[1]), x.breaks=12, show.sat = TRUE, show.cycle = TRUE, anoxic=anox, incub.lag=i.lag, light.lag=l.lag, dark.lag=d.lag, exchange.lag=e.lag, runtime=tmax)
  o2_singrep[[i]]<-oxyplot(oxydata=oxydat_single, lightdata=lights, x.lims=c(min(oxydat_single$totaltime)[1], 48), x.breaks=12, show.sat = TRUE, show.cycle = TRUE, anoxic=anox, incub.lag=i.lag, light.lag=l.lag, dark.lag=d.lag, exchange.lag=e.lag, runtime=tmax)
}

#### ELIMINATE IMPLAUSIBLE VALUES ####
# based on O2 profiles over the measurement time, decide if there is any implausible values, or samples/channels that failed for whatever reason ... and remove those by hand
# manual filtering based on visual inspection (device not working properly):
rn_out <- oxydat %>% 
  rownames_to_column() %>% 
  filter(treatment == "aqua",
         replicate == 4,
         device == "Wipper", 
         light == "off", 
         cycle == 1, 
         oxy_mg_l > oxy100_mg_l) %>% 
  pull(rowname)

oxydat <- oxydat %>%
  filter(!(row.names(.) %in% rn_out))

# find (strong) outliers that are below Q1-[outlier.thresh_O2]*IQR or above Q3+[outlier.thresh_O2]*IQR (rule of interquartile range)
Q1 <- unname(summary(oxydat$oxy_mg_l)[2])
Q3 <- unname(summary(oxydat$oxy_mg_l)[5])
iqr <- IQR(oxydat$oxy_mg_l)
lower_boundary <- Q1 - outlier.thresh_O2 * iqr
upper_boundary <- Q3 + outlier.thresh_O2 * iqr
oxydat_out <- oxydat[oxydat$oxy_mg_l < lower_boundary | oxydat$oxy_mg_l > upper_boundary, ]

if (remove_outliers == TRUE) {
  oxydat <- oxydat[!(row.names(oxydat) %in% row.names(oxydat_out)), ]
}

o2_profile_init <- oxyplot(oxydata = oxydat %>% filter(treatment != "blind") %>% droplevels(), x.lims = c(0, 48), lightdata = lights, x.breaks = 24, show.sat = TRUE, show.cycle = TRUE, anoxic = anox, incub.lag = i.lag, light.lag = l.lag, dark.lag = d.lag, exchange.lag = e.lag, runtime = tmax, show.outlier = TRUE, threshold.outlier = outlier.thresh_O2)
o2_profile_init
o2_profile_fin <- oxyplot(oxydata = oxydat %>% filter(treatment != "blind") %>% droplevels(), x.lims = c(240, 288), lightdata = lights, x.breaks = 24, show.sat = TRUE, show.cycle = TRUE, anoxic = anox, incub.lag = i.lag, light.lag = l.lag, dark.lag = d.lag, exchange.lag = e.lag, runtime = tmax, show.outlier = TRUE, threshold.outlier = outlier.thresh_O2)
o2_profile_fin

#### CALCULATE O2 RATES ####
## Community respiration (CR) = O2 measurements from dark periods (note: CR is ALWAYS negative, i.e. <0)
CR <- oxyrate.batch(oxydata = oxydat, lightdata = lights, blind = blank_name, light = "off", Q10 = q10, reference = "sediment", incub.lag = i.lag, light.lag = l.lag, dark.lag = d.lag, exchange.lag = e.lag, runtime = tmax, anoxic = anox, R2.min = r2.min)
CR <- droplevels(CR)
CR$treatment <- factor(CR$treatment, levels = c(treatments, blank_name))
CR <- CR %>%
  mutate(logclass = factor(paste("Day", cycle), levels = paste("Day", sort(unique(cycle)))))

## Net ecosystem production (NCP) = O2 measurements from light periods (note: NCP can be either positive or negative, but it must be bigger than CR!)
NCP <- oxyrate.batch(oxydata = oxydat, lightdata = lights, blind = blank_name, light = "on", Q10 = q10, reference = "sediment", incub.lag = i.lag, light.lag = l.lag, dark.lag = d.lag, exchange.lag = e.lag, runtime = tmax, anoxic = anox, R2.min = r2.min)
NCP <- droplevels(NCP)
NCP$treatment <- factor(NCP$treatment, levels = c(treatments, blank_name))
NCP <- NCP %>%
  mutate(logclass = factor(paste("Day", cycle), levels = paste("Day", sort(unique(cycle)))))

# Calculate Gross Primary Production (GPP) from NCP+CR (CR <0) (note: GPP is ALWAYS positive, i.e. >0)
GPP <- rbind(CR, NCP) %>%
  group_by(treatment, replicate, device, cycle) %>%
  summarize(
    rate = rate[variable == "NCP"] - rate[variable == "CR"],
    blankrate = blankrate[variable == "NCP"] - blankrate[variable == "CR"],
    netrate = netrate[variable == "NCP"] - netrate[variable == "CR"],
    rate_dw = rate_dw[variable == "NCP"] - rate_dw[variable == "CR"],
    netrate_dw = netrate_dw[variable == "NCP"] - netrate_dw[variable == "CR"]
  ) %>%
  mutate(across(c("rate", "blankrate", "netrate", "rate_dw", "netrate_dw"), ~ ifelse(.x < 0, 0, .x))) %>%
  mutate(variable = "GPP") %>%
  as.data.frame()
GPP <- droplevels(GPP)
GPP$treatment <- factor(GPP$treatment, levels = c(treatments, blank_name))
GPP$totaltimemean <- left_join(GPP, NCP, by = c("treatment", "replicate", "device", "cycle"))$totaltimemean
GPP$glass_ID <- left_join(GPP, NCP, by = c("treatment", "replicate", "device", "cycle"))$glass_ID
GPP$r2 <- NA
GPP <- GPP %>%
  mutate(logclass = factor(paste("Day", cycle), levels = paste("Day", sort(unique(cycle)))))

# # missing observations
# CR_miss<-CR %>%
#   filter(is.na(rate)) %>%
#   group_by(treatment, replicate, device, variable) %>%
#   summarise(missing_cycles=paste(cycle, collapse = ", ")) %>%
#   as.data.frame()
#
# NCP_miss<-NCP %>%
#   filter(is.na(rate)) %>%
#   group_by(treatment, replicate, device, variable) %>%
#   summarise(missing_cycles=paste(cycle, collapse = ", ")) %>%
#   as.data.frame()
#
# GPP_miss<-GPP %>%
#   filter(is.na(rate)) %>%
#   group_by(treatment, replicate, device, variable) %>%
#   summarise(missing_cycles=paste(cycle, collapse = ", ")) %>%
#   as.data.frame()

# calculate means of all replicates
CR_means <- CR %>%
  group_by(treatment, cycle, variable, device, logclass) %>%
  filter(treatment != blank_name) %>%
  summarise(
    rate_dw_mean = mean(rate_dw, na.rm = T),
    rate_dw_sd = sd(rate_dw, na.rm = T),
    netrate_dw_mean = mean(netrate_dw, na.rm = T),
    netrate_dw_sd = sd(netrate_dw, na.rm = T),
    # blankrate=mean(blankrate, na.rm=T),
    replicates = n(),
    meantotaltime = mean(totaltimemean, na.rm = T)
  ) %>%
  distinct() %>%
  ungroup() %>%
  mutate(na_if(., "NaN")) %>%
  as.data.frame()

NCP_means <- NCP %>%
  group_by(treatment, cycle, variable, device, logclass) %>%
  filter(treatment != blank_name) %>%
  summarise(
    rate_dw_mean = mean(rate_dw, na.rm = T),
    rate_dw_sd = sd(rate_dw, na.rm = T),
    netrate_dw_mean = mean(netrate_dw, na.rm = T),
    netrate_dw_sd = sd(netrate_dw, na.rm = T),
    # blankrate=mean(blankrate, na.rm=T),
    replicates = n(),
    meantotaltime = mean(totaltimemean, na.rm = T)
  ) %>%
  distinct() %>%
  ungroup() %>%
  mutate(na_if(., "NaN")) %>%
  as.data.frame()

GPP_means <- GPP %>%
  group_by(treatment, cycle, variable, device, logclass) %>%
  filter(treatment != blank_name) %>%
  summarise(
    rate_dw_mean = mean(rate_dw, na.rm = T),
    rate_dw_sd = sd(rate_dw, na.rm = T),
    netrate_dw_mean = mean(netrate_dw, na.rm = T),
    netrate_dw_sd = sd(netrate_dw, na.rm = T),
    # blankrate=mean(blankrate, na.rm=T),
    replicates = n(),
    meantotaltime = mean(totaltimemean, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(na_if(., "NaN")) %>%
  as.data.frame()

# combine all CR, NCP, GPP values to single data.frame
common.cols1 <- intersect(names(rbind(CR, NCP)), names(GPP))
rates <- rbind(CR[, common.cols1], NCP[, common.cols1], GPP)
common.cols2 <- intersect(names(rbind(CR_means, NCP_means)), names(GPP_means))
meanrates <- rbind(CR_means[, common.cols2], NCP_means[, common.cols2], GPP_means)

#### PLOTTING RATES ####
if (blind_correction == TRUE) {
  plotrate <- "netrate_dw" # which variable shall be plotted on the y-axis of scatter/box plots (rate_dw: without blank correction; netrate_dw: with blank correction; both related to dry weight (dw) of sediment)
} else {
  plotrate <- "rate_dw"
}

tit <- paste(paste(project, stream, paste0("(", ifelse(length(dates) == 1, dates, paste(c(dates[1], dates[length(dates)]), collapse = "-")), ")")), paste0(round(max(oxydat$totaltime), 0), " h incubation"), sep = ", ")
subtit <- paste(paste("Minimum R2:", r2.min), paste("anoxic at:", anox, "mg/l"), paste("incubation lag:", i.lag, "h"), paste("light lag:", l.lag, "h"), paste("dark lag:", d.lag, "h"), paste("exchange lag:", e.lag, "h"), paste("Blind:", ifelse(blind_correction == TRUE, "yes", "no")), paste("Q10:", ifelse(q10 == TRUE, "yes", "no")), paste("Outliers:", paste0(outlier.thresh_rate, "*IQR")), sep = ", ")
subtit_blind <- paste(paste("Minimum R2:", r2.min), paste("anoxic at:", anox, "mg/l"), paste("incubation lag:", i.lag, "h"), paste("light lag:", l.lag, "h"), paste("dark lag:", d.lag, "h"), paste("exchange lag:", e.lag, "h"), sep = ", ")

## BLIND RATES
plotdat <- rates[(rates$treatment %in% blank_name), ]
plotdat <- droplevels(plotdat)
plotdat$treatment <- dplyr::recode(plotdat$treatment, !!!levels_sediment)
plotdat$device <- dplyr::recode(plotdat$device, !!!levels_transport)

scatter_blind <- ggplot(data = plotdat, aes(x = logclass, y = rate, colour = replicate)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, colour = "grey40") +
  labs(title = tit, subtitle = subtit_blind, x = NULL, colour = "Replicate", y = expression(paste(Delta * O[2] ~ "[", mgO[2] ~ h^{
    -1
  }, L^{
    -1
  }, "]"), sep = "")) +
  facet_grid(variable ~ .) +
  theme_bw()
scatter_blind

## CR
plotdat <- CR[!(CR$treatment %in% blank_name), ]
plotdat <- droplevels(plotdat)
plotdat$treatment <- dplyr::recode(plotdat$treatment, !!!levels_sediment)
plotdat$device <- dplyr::recode(plotdat$device, !!!levels_transport)
nrow(plotdat)

# identify outliers (per group):
plotdat_out <- plotdat %>%
  group_by(treatment, device, cycle) %>%
  mutate(
    Q1 = unname(summary(netrate_dw))[2],
    Q3 = unname(summary(netrate_dw))[5],
    iqr = IQR(netrate_dw),
    lower = Q1 - outlier.thresh_rate * iqr,
    upper = Q3 + outlier.thresh_rate * iqr
  ) %>%
  ungroup() %>%
  mutate(flag = ifelse(netrate_dw < lower | netrate_dw > upper, 0, 1)) %>%
  filter(flag == 0) %>%
  select(-c(Q1, Q3, iqr, lower, upper, flag)) %>%
  mutate(na_if(., "NaN")) %>%
  as.data.frame()
nrow(plotdat_out)

# complete data for all cases (unfortunately needed to correctly position the points on the plot)
all <- plotdat %>% expand(nesting(treatment, replicate, glass_ID, device), nesting(cycle, logclass))
plotdat_out <- plotdat_out %>% right_join(all)

scatter_cr <- ggplot(data = plotdat, aes(x = treatment, y = -(!!as.name(plotrate)), colour = device)) +
  geom_point(size = 2.5, aes(group = interaction(device, treatment)), position = position_dodge(0.75, preserve = "total")) +
  geom_point(data = plotdat_out, colour = "orange", inherit.aes = TRUE, size = 2.5, aes(group = interaction(device, treatment)), position = position_dodge(0.75, preserve = "total")) +
  geom_hline(yintercept = 0, colour = "grey40") +
  labs(title = paste0(tit, ", CR"), subtitle = subtit, x = "Sediment", colour = "Transport", y = expression(paste("-" * Delta * O[2] ~ "[", mu, gO[2] ~ h^{
    -1
  }, gDW^{
    -1
  }, "]"), sep = "")) +
  facet_grid(. ~ logclass) +
  theme_bw()
scatter_cr

boxplot_cr <- ggplot(data = plotdat, aes(x = treatment, y = -(!!as.name(plotrate)), fill = device)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  geom_point(data = plotdat_out, colour = "black", fill = "orange", shape = 21, inherit.aes = TRUE, aes(group = interaction(device, treatment)), size = 2.5, position = position_dodge(0.75, preserve = "total")) +
  geom_hline(yintercept = 0, colour = "grey40") +
  stat_summary(fun.data = give.n, size = 3, colour = "grey30", geom = "text", fun = quantile, position = position_dodge(width = 0.75)) +
  labs(title = paste0(tit, ", CR"), subtitle = subtit, x = "Sediment", fill = "Transport", y = expression(paste("-" * Delta * O[2] ~ "[", mu, gO[2] ~ h^{
    -1
  }, gDW^{
    -1
  }, "]"), sep = "")) +
  facet_grid(. ~ logclass) +
  theme_bw()
boxplot_cr

## NCP
plotdat <- NCP[!(NCP$treatment %in% blank_name), ]
plotdat <- droplevels(plotdat)
plotdat$treatment <- dplyr::recode(plotdat$treatment, !!!levels_sediment)
plotdat$device <- dplyr::recode(plotdat$device, !!!levels_transport)
nrow(plotdat)

# identify outliers (per group):
plotdat_out <- plotdat %>%
  group_by(treatment, device, cycle) %>%
  mutate(
    Q1 = unname(summary(netrate_dw))[2],
    Q3 = unname(summary(netrate_dw))[5],
    iqr = IQR(netrate_dw),
    lower = Q1 - outlier.thresh_rate * iqr,
    upper = Q3 + outlier.thresh_rate * iqr
  ) %>%
  ungroup() %>%
  mutate(flag = ifelse(netrate_dw < lower | netrate_dw > upper, 0, 1)) %>%
  filter(flag == 0) %>%
  select(-c(Q1, Q3, iqr, lower, upper, flag)) %>%
  mutate(na_if(., "NaN")) %>%
  as.data.frame()
nrow(plotdat_out)

# complete data for all cases (unfortunately needed to correctly position the points on the plot)
all <- plotdat %>% expand(nesting(treatment, replicate, glass_ID, device), nesting(cycle, logclass))
plotdat_out <- plotdat_out %>% right_join(all)

scatter_ncp <- ggplot(data = plotdat, aes(x = treatment, y = (!!as.name(plotrate)), colour = device)) +
  geom_point(size = 2.5, aes(group = interaction(device, treatment)), position = position_dodge(0.75, preserve = "total")) +
  geom_point(data = plotdat_out, colour = "orange", inherit.aes = TRUE, size = 2.5, aes(group = interaction(device, treatment)), position = position_dodge(0.75, preserve = "total")) +
  geom_hline(yintercept = 0, colour = "grey40") +
  labs(title = paste0(tit, ", NCP"), subtitle = subtit, x = "Sediment", colour = "Transport", y = expression(paste(Delta * O[2] ~ "[", mu, gO[2] ~ h^{
    -1
  }, gDW^{
    -1
  }, "]"), sep = "")) +
  facet_grid(. ~ logclass) +
  theme_bw()
scatter_ncp

boxplot_ncp <- ggplot(data = plotdat, aes(x = treatment, y = (!!as.name(plotrate)), fill = device)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  geom_point(data = plotdat_out, colour = "black", fill = "orange", shape = 21, inherit.aes = TRUE, aes(group = interaction(device, treatment)), size = 2.5, position = position_dodge(0.75, preserve = "total")) +
  geom_hline(yintercept = 0, colour = "grey40") +
  stat_summary(fun.data = give.n, size = 3, colour = "grey30", geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  labs(title = paste0(tit, ", NCP"), subtitle = subtit, x = "Sediment", fill = "Transport", y = expression(paste(Delta * O[2] ~ "[", mu, gO[2] ~ h^{
    -1
  }, gDW^{
    -1
  }, "]"), sep = "")) +
  facet_grid(. ~ logclass) +
  theme_bw()
boxplot_ncp

## GPP
plotdat <- GPP[!(GPP$treatment %in% blank_name), ]
plotdat <- droplevels(plotdat)
plotdat$treatment <- dplyr::recode(plotdat$treatment, !!!levels_sediment)
plotdat$device <- dplyr::recode(plotdat$device, !!!levels_transport)
nrow(plotdat)

# identify outliers (per group):
plotdat_out <- plotdat %>%
  group_by(treatment, device, cycle) %>%
  mutate(
    Q1 = unname(summary(netrate_dw))[2],
    Q3 = unname(summary(netrate_dw))[5],
    iqr = IQR(netrate_dw),
    lower = Q1 - outlier.thresh_rate * iqr,
    upper = Q3 + outlier.thresh_rate * iqr
  ) %>%
  ungroup() %>%
  mutate(flag = ifelse(netrate_dw < lower | netrate_dw > upper, 0, 1)) %>%
  filter(flag == 0) %>%
  select(-c(Q1, Q3, iqr, lower, upper, flag)) %>%
  mutate(na_if(., "NaN")) %>%
  as.data.frame()
nrow(plotdat_out)

# complete data for all cases (unfortunately needed to correctly position the points on the plot)
all <- plotdat %>% expand(nesting(treatment, replicate, glass_ID, device), nesting(cycle, logclass))
plotdat_out <- plotdat_out %>% right_join(all)

scatter_gpp <- ggplot(data = plotdat, aes(x = treatment, y = (!!as.name(plotrate)), colour = device)) +
  geom_point(size = 2.5, aes(group = interaction(device, treatment)), position = position_dodge(0.75, preserve = "total")) +
  geom_point(data = plotdat_out, colour = "orange", inherit.aes = TRUE, size = 2.5, aes(group = interaction(device, treatment)), position = position_dodge(0.75, preserve = "total")) +
  geom_hline(yintercept = 0, colour = "grey40") +
  labs(title = paste0(tit, ", GPP"), subtitle = subtit, x = "Sediment", colour = "Transport", y = expression(paste(Delta * O[2] ~ "[", mu, gO[2] ~ h^{
    -1
  }, gDW^{
    -1
  }, "]"), sep = "")) +
  facet_grid(. ~ logclass) +
  theme_bw()
scatter_gpp

boxplot_gpp <- ggplot(data = plotdat, aes(x = treatment, y = (!!as.name(plotrate)), fill = device)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  geom_point(data = plotdat_out, colour = "black", fill = "orange", shape = 21, inherit.aes = TRUE, aes(group = interaction(device, treatment)), size = 2.5, position = position_dodge(0.75, preserve = "total")) +
  geom_hline(yintercept = 0, colour = "grey40") +
  stat_summary(fun.data = give.n, size = 3, colour = "grey30", geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  labs(title = paste0(tit, ", GPP"), subtitle = subtit, x = "Sediment", fill = "Transport", y = expression(paste(Delta * O[2] ~ "[", mu, gO[2] ~ h^{
    -1
  }, gDW^{
    -1
  }, "]"), sep = "")) +
  facet_grid(. ~ logclass) +
  theme_bw()
boxplot_gpp

#### WRITE OUTPUT FILES ####

## COMMENTS
if (blind_correction == TRUE) {
  comment1 <- "withblind"
} else {
  comment1 <- "noblind"
}

if (remove_outliers == TRUE) {
  comment2 <- "noO2outliers"
} else {
  comment2 <- "withO2outliers"
}

## DATA
setwd(paste(path_out, "Tables", sep = "/"))

write.table(oxydat, file = paste0(paste("Oxygen", project, comment2, sep = "_"), ".txt"), sep = "\t", row.names = FALSE, quote = F)
write.table(rates, file = paste0(paste("Metabolism", "alldays", "single", project, comment2, sep = "_"), ".txt"), sep = "\t", row.names = FALSE, quote = F)
write.table(meanrates, file = paste0(paste("Metabolism", "alldays", "means", project, comment2, sep = "_"), ".txt"), sep = "\t", row.names = FALSE, quote = F)

## PLOTS
setwd(path_out)

ggsave(plot = o2_profile, filename = paste0(paste("O2", project, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")

ggsave(plot = scatter_blind, filename = paste0(paste("Blind", "alldays", "scatterplot", project, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")

ggsave(plot = scatter_cr, filename = paste0(paste("CR", "alldays", "scatterplot", project, comment1, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")
ggsave(plot = boxplot_cr, filename = paste0(paste("CR", "alldays", "boxplot", project, comment1, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")
ggsave(plot = scatter_ncp, filename = paste0(paste("NCP", "alldays", "scatterplot", project, comment1, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")
ggsave(plot = boxplot_ncp, filename = paste0(paste("NCP", "alldays", "boxplot", project, comment1, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")
ggsave(plot = scatter_gpp, filename = paste0(paste("GPP", "alldays", "scatterplot", project, comment1, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")
ggsave(plot = boxplot_gpp, filename = paste0(paste("GPP", "alldays", "boxplot", project, comment1, sep = "_"), ".png"), dpi = 300, width = 10, height = 5, units = "in")

setwd(path_in)
