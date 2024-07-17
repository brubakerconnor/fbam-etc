# utility functions
# convert run time to seconds if not already
convert_to_seconds <- function(s) {
  unit <- sub("^[0-9.]+ ", "", s)
  if (unit == "secs") return(s)
  if (unit == "mins") {
    time <- sub(" mins$", "", s)
    return(as.double(time) * 60)
  }
}

# function to get mode of a vector
get_mode <- function(v) {
  uniq_vals <- unique(v)
  uniq_vals[which.max(tabulate(match(v, uniq_vals)))]
}

# collect all file names to include in the report
fnames <- list.files("~/Data/fbam/sensitivity-analysis", pattern = "\\.rda$",
                     full.names = TRUE)
report <- data.frame()

for (fname in fnames) {
  print(fname)
  load(fname)
  
  # get parameters of the run
  model <- strsplit(tail(strsplit(fname, "/")[[1]], 1), "_")[[1]][1]
  nrep <- dim(sim_data[[1]]$data$x)[2] / length(unique(sim_data[[1]]$data$labels))
  len <- dim(sim_data[[1]]$data$x)[1]
  nsim <- length(sim_data)
  pcrossover <- sim_data[[1]]$fit$selected$params$pcrossover
  pmutate <- sim_data[[1]]$fit$selected$params$pmutate
  endpoint_range <- sim_data[[1]]$fit$selected$params$endpoint_range[1]
  collapsed_range <- sim_data[[1]]$fit$selected$params$collapsed_range[1]
  
  # collect data
  ari <- rep(0, nsim)
  nclust <- nbands <- rep(0, nsim)
  loss <- time <- rep(0, nsim)
  for (n in 1:nsim) {
    ari[n] <- mclust::adjustedRandIndex(
      sim_data[[n]]$fit$selected$labels, 
      sim_data[[n]]$data$labels
    )
    nclust[n] <- sim_data[[n]]$fit$selected$params$nclust
    nbands[n] <- sim_data[[n]]$fit$selected$params$nbands
    loss[n] <- sim_data[[n]]$fit$selected$loss
    time[n] <- convert_to_seconds(sim_data[[n]]$time)
  }
  x <- data.frame(
    model = model,
    nrep = nrep,
    len = len,
    pcrossover = pcrossover,
    pmutate = pmutate,
    endpoint_range = endpoint_range,
    collapsed_range = collapsed_range,
    mode_nclust = get_mode(nclust),
    percent_nclust = sum(nclust == get_mode(nclust)),
    mode_nbands = get_mode(nbands),
    percent_nbands = sum(nclust == get_mode(nbands)),
    mean_loss = mean(loss),
    median_loss = median(loss),
    se_loss = sd(loss),
    mad_loss = mad(loss),
    mean_time_secs = mean(time),
    se_time_secs = sd(time)
  )
  report <- rbind(report, x)
  rm(sim_data)
}
write.csv(report, "report.csv", row.names = FALSE)
