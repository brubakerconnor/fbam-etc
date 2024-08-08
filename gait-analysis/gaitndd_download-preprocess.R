rm(list = ls())
set.seed(451)

# download files from physionet directly
if (!dir.exists('gait-analysis/data-raw/')) {
  dir.create('gait-analysis/data-raw/')
  dir.create('gait-analysis/data/')
  system('wget -r -N -c -np https://physionet.org/files/gaitndd/1.0.0/ -P gait-analysis/data-raw/')
}
DATA_DIR <- 'gait-analysis/data-raw/physionet.org/files/gaitndd/1.0.0'

# note the following patients were removed due to too many artefacts in the
# stride interval series: als12, als4, park11

# get the path to each of the ts files in the collection
ts_files <- list.files(DATA_DIR, pattern = "\\.ts$", full.names = T)
patient_id <- tools::file_path_sans_ext(basename(ts_files))
labels <- gsub("[0-9]", "", patient_id)

# skip the patient files described above
skip <- c("als12", "als4", "park11") # skip these due to issues in the signal
ts_files <- ts_files[!patient_id %in% skip]
labels <- labels[!patient_id %in% skip]
patient_id <- patient_id[!patient_id %in% skip]

# signal and estimation parameters
duration <- 210 # seconds
sample_rate <- 2 # samples per second
ts_length <- duration * sample_rate
mtfreq <- seq(from = 1 / ts_length, by = 1 / ts_length,
              length = floor(ts_length / 2) - 1)
ntapers <- floor(sqrt(ts_length))

# process each file
xraw <- list()
x <- matrix(nrow = ts_length, ncol = length(ts_files))
mtspec <- matrix(nrow = length(mtfreq), ncol = length(ts_files))
for (i in 1:length(ts_files)) {
  print(paste0("Processing data contained in ", ts_files[i]))

  # read data and define variables
  dat <- readr::read_delim(ts_files[i], col_names = F, show_col_types = F)
  elapsed_time <- dat[[1]]; left_stride <- diff(elapsed_time); elapsed_time <- elapsed_time[-1]
  xraw[[i]] <- list()
  xraw[[i]]$start <- dat[[1]][1]
  xraw[[i]]$elapsed <- dat[[1]]
  xraw[[i]]$left <- dat[[2]]

  # outlier filtering
  # time points that lie outside the percentiles defined below are considered
  # outliers and are replaced by the median value of the series
  median_stride <- median(left_stride)
  quantiles <- quantile(left_stride, c(0.01, 0.96))
  replace <- left_stride < quantiles[1] | left_stride > quantiles[2]
  left_stride[replace] <- median_stride

  # linear interpolation to make regularly sampled series
  start_time <- min(elapsed_time); end_time <- start_time + duration
  time_grid <- seq(from = start_time, to = end_time, by = 1 / sample_rate)
  time_grid <- time_grid[1:ts_length]
  left_stride <- approx(x = elapsed_time, y = left_stride, xout = time_grid,
                        method = "linear")$y

  # remove linear trend
  left_stride <- as.vector(pracma::detrend(left_stride, tt = "linear"))

  # standardize by standard deviation of series
  left_stride <- left_stride / sd(left_stride)
  x[, i] <- left_stride # add to ts data

  # multitaper estimate
  mtspec[, i] <- as.vector(fbam::sine_multitaper(left_stride)$mtspec)
}

# read and clean up the subject-level descriptors file
desc <- read.delim("gait-analysis/data-raw/physionet.org/files/gaitndd/1.0.0/subject-description.txt",
                   row.names = NULL, na.strings = "MISSING")
names(desc) <- c("patient_id", "group", "age", "height", "weight", "gender",
                 "gait_speed", "duration_severity")
desc <- desc[desc$patient_id != "", ] # remove that weird mostly empty row
desc$group[desc$group == "subjects"] <- "als"

# impute those missing als gait speed values with group mean gait speed
als_mean_speed <- mean(desc[desc$group == "als", 'gait_speed'], na.rm = TRUE)
desc$gait_speed[is.na(desc$gait_speed) & desc$group == 'als'] <- als_mean_speed

# The duration/severity column contains condition-specific measures of disease
# duration or severity. For control patients, this value is 0. For Parkinson's
# patients, this is the Hohn and Yahr score which ranges from 1 to 5 with higher
# values indicating more advanced disease. For Huntington's patients, this is
# the total functional capacity measure (lower scores mean more advanced
# functional impairment) which ranges from 0 to 13. Finally, for ALS patients,
# this is simply the number of months since diagnosis.

# in order to provide a standardized scale, the duration or severity within
# each group is re-calculated in the following way:
# PD: 0.25 * score - 0.25
# HD: (1/13) * (13 - score) - 1/13
# ALS: (1/max_duration) * score - (1/max_duration)
# Standardization is done in this way to put all values between 0 and 1 where
# a value of 0 is "best" and 1 is "worst".
desc$duration_severity_std <- desc$duration_severity
desc$duration_severity_std[desc$group == "park"] <- 0.25 *
  desc$duration_severity[desc$group == "park"] - 0.25
desc$duration_severity_std[desc$group == "hunt"] <- (1/13) *
  (13 - desc$duration_severity[desc$group == "hunt"] - (1/13))

max_duration <- max(desc$duration_severity_std[desc$group == "als"])
desc$duration_severity_std[desc$group == "als"] <- (1 / max_duration) *
  desc$duration_severity[desc$group == "als"] - (1 / max_duration)

# plots
## original time series against elapsed time
pdf("gait-analysis/data/all_stride_ts_raw.pdf", width = 16, height = 64)
par(mfrow = c(16, 4))
for (i in 1:length(ts_files)) {
  plot(xraw[[i]]$elapsed, xraw[[i]]$left, type = 'l', main = patient_id[i])
  abline(v = 1:5 * 77 / desc[desc$patient_id == patient_id[i], 'gait_speed'],
         col = 2, lty = 2)
}
dev.off()

## processed time series against elapsed time
pdf("gait-analysis/data/all_stride_ts.pdf", width = 16, height = 64)
par(mfrow = c(16, 4))
for (i in 1:ncol(x)) {
  plot(seq(1, 210, length = nrow(x)), x[,i],
                            type = 'l',
                            main = patient_id[i])
  abline(v = 1:5 * 77 / desc[desc$patient_id == patient_id[i], 'gait_speed'],
         col = 2, lty = 2)
}
dev.off()

# self similarity parameter
ssp <- rep(0, ncol(x))
for (i in 1:ncol(x)) {
  print(i)
  ssp[i] <- DFA::SSP(x[,i])
}

# no clustering and clustering solutions
noclust <- fbam::fbam(gait$x, 1, 2:6, parallel = TRUE)$selected
noclust$rep_collapsed <- data.frame(
  fbam:::rep_collapsed_by_index(gait$mtspec,
                                noclust$endpoints_index,
                                noclust$labels)
)
names(noclust$rep_collapsed) <- c("LF", "HF")

clust <- fbam::fbam(gait$x, 2:6, 2:6, parallel = TRUE)$selected
clust$rep_collapsed <- data.frame(
  fbam:::rep_collapsed_by_index(gait$mtspec,
                                clust$endpoints_index,
                                clust$labels)
)
names(clust$rep_collapsed) <- c("LF", "HF")

# create data object and save
x <- data.frame(x); names(x) <- patient_id
mtspec <- data.frame(mtspec); names(mtspec) <- patient_id
gait <- list(
  x = x,
  mtfreq = mtfreq,
  mtspec = mtspec,
  labels = labels,
  patient_id = patient_id,
  ssp = ssp,
  noclust = noclust,
  clust = clust
)
save(gait, file = "gait-analysis/data/gait.rda")

# create simplified csv file
# add the LF, HF, and LF:HF measures from the no clustering and clustering
# solutions along with cluster labels from the clustering solution
noclust_bap <- fbam:::rep_collapsed_by_index(
  gait$noclust$spec,
  gait$noclust$endpoints_index
)
noclust_df <- data.frame(patient_id = gait$patient_id,
                         noclust_LF = noclust_bap[,1],
                         noclust_HF = noclust_bap[,2],
                         noclust_ratio = noclust_bap[,1] / noclust_bap[,2])

clust_bap <- fbam:::rep_collapsed_by_index(
  gait$clust$spec,
  gait$clust$endpoints_index,
  gait$clust$labels
)
clust_df <- data.frame(patient_id = gait$patient_id,
                       clust_LF = clust_bap[,1],
                       clust_HF = clust_bap[,2],
                       clust_ratio = clust_bap[,1] / clust_bap[,2],
                       clust_label = gait$clust$labels)

desc <- merge(desc, noclust_df, by = "patient_id")
desc <- merge(desc, clust_df, by = "patient_id")

# add endpoint columns
desc$noclust_endpoint <- noclust$endpoints[1]
desc$clust_endpoint <- clust$endpoints[desc$clust_label,]

# add SSP column
desc$ssp <- gait$ssp

# save results to csv
write.csv(desc, file = "gait-analysis/data/gait.csv", row.names = FALSE)

