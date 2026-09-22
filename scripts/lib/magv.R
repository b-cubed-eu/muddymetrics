# MAGV saturation reliability, as used in D6.1 Part 2. Reconstructed from the
# D6.1 outputs (magv_by_year.csv, magv_summary.csv), which it reproduces exactly.

#' MAGV saturation reliability (D6.1 Part 2): on the annual total-occurrence
#' series (years with records, in order), relative change r_t = n_t / n_{t-1} - 1;
#' over a right-aligned window of `window` years, mean |r| (lma_abs) and sd(r)
#' (volatility), ignoring the undefined first change; reliability =
#' 1 / (1 + volatility + lma_abs). Years with reliability > threshold are the
#' "Plateau of Stability", all others "Discovery Phase".
magv_series <- function(year, total_occ, window = 5, threshold = 0.8) {
  o <- order(year); year <- year[o]; n <- total_occ[o]
  k <- length(n)
  rc <- c(NA_real_, n[-1] / n[-k] - 1)
  lma <- vol <- rep(NA_real_, k)
  if (k >= window) for (i in window:k) {
    w <- rc[(i - window + 1):i]
    lma[i] <- mean(abs(w), na.rm = TRUE)
    vol[i] <- sd(w, na.rm = TRUE)
  }
  rel <- 1 / (1 + vol + lma)
  data.frame(year = year, total_occ = n, relative_change = rc, lma_abs = lma,
             volatility = vol, reliability = rel,
             phase = ifelse(!is.na(rel) & rel > threshold, "Plateau of Stability", "Discovery Phase"))
}

#' Longest run of consecutive stable rows; usable if at least `min_run` long.
magv_summary <- function(s, min_run = 5) {
  st <- s$phase == "Plateau of Stability"
  r <- rle(st); len <- 0L; start <- end <- NA_real_
  if (any(r$values)) {
    j <- which(r$values)[which.max(r$lengths[r$values])]
    endi <- cumsum(r$lengths)[j]; starti <- endi - r$lengths[j] + 1
    len <- r$lengths[j]; start <- as.numeric(s$year[starti]); end <- as.numeric(s$year[endi])
  }
  data.frame(n_years = nrow(s), first_year = min(s$year), last_year = max(s$year),
             mean_reliability = mean(s$reliability, na.rm = TRUE),
             stable_run_start = start, stable_run_end = end, stable_run_length = len,
             has_usable_window = len >= min_run)
}
