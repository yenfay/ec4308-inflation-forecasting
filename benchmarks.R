
WINDOW <- 240
HORIZONS <- c(1L, 3L, 6L, 12L)
P_MAX <- 12

SAMPLES <- c(
  official = "main_data_official.rds",
  self = "main_data_self.rds"
)

dir.create("results/forecasts", recursive = TRUE, showWarnings = FALSE)

load_sample <- function(filename) {
  df <- readRDS(file.path("data/processed", filename))
  df <- df[order(df$date), ]
  df <- df[!is.na(df$inflation_mom), ]
  
  dates <- as.Date(df$date)
  inflation <- as.numeric(df$inflation_mom)
  
  month_id <- as.integer(format(dates, "%Y")) * 12 +
    as.integer(format(dates, "%m"))
  
  stopifnot(
    length(inflation) > WINDOW,
    !anyNA(dates),
    all(diff(month_id) == 1)
  )
  
  list(dates = dates, inflation = inflation)
}

train_rows <- function(t, h) {
  max(1, t - WINDOW):(t - h)
}

fit_ar_forecast <- function(z, h) {
  n <- length(z)
  s <- (h + P_MAX):n
  
  y <- z[s]
  
  X <- sapply(0:(P_MAX - 1), function(k) {
    z[s - h - k]
  })
  
  colnames(X) <- paste0("lag", 1:P_MAX)
  
  fits <- vector("list", P_MAX)
  bics <- numeric(P_MAX)
  
  for (p in 1:P_MAX) {
    dat <- data.frame(y = y, X[, 1:p, drop = FALSE])
    fits[[p]] <- lm(y ~ ., data = dat)
    bics[p] <- BIC(fits[[p]])
  }
  
  p_best <- which.min(bics)
  
  x_now <- z[n - 0:(p_best - 1)]
  names(x_now) <- paste0("lag", 1:p_best)
  
  forecast <- predict(
    fits[[p_best]],
    newdata = as.data.frame(as.list(x_now))
  )
  
  list(
    forecast = as.numeric(forecast),
    selected_p = p_best
  )
}

run_benchmarks <- function(data, sample_name) {
  pi <- data$inflation
  dates <- data$dates
  n <- length(pi)
  
  results <- list()
  counter <- 1
  
  for (h in HORIZONS) {
    for (t in (WINDOW + 1):n) {
      
      z <- pi[train_rows(t, h)]
      
      rw_pred <- tail(z, 1)
      ar <- fit_ar_forecast(z, h)
      
      results[[counter]] <- data.frame(
        sample = sample_name,
        model = c("RW", "AR_BIC"),
        h = h,
        origin_date = as.character(dates[t - h]),
        target_date = as.character(dates[t]),
        forecast = c(rw_pred, ar$forecast),
        actual = pi[t],
        selected_p = c(NA_integer_, ar$selected_p)
      )
      
      counter <- counter + 1
    }
    
    cat("Completed:", sample_name, "h =", h, "\n")
  }
  
  do.call(rbind, results)
}

all_results <- list()

for (sample_name in names(SAMPLES)) {
  
  data <- load_sample(SAMPLES[[sample_name]])
  res <- run_benchmarks(data, sample_name)
  
  write.csv(
    res[, c("sample", "model", "h", "origin_date",
            "target_date", "forecast", "actual")],
    paste0("results/forecasts/forecasts_benchmarks_",
           sample_name, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    subset(res, model == "AR_BIC",
           select = c(sample, h, origin_date,
                      target_date, selected_p)),
    paste0("results/forecasts/ar_selected_lags_",
           sample_name, ".csv"),
    row.names = FALSE
  )
  
  all_results[[sample_name]] <- res
}

official <- all_results$official

stopifnot(
  nrow(official) == 174 * length(HORIZONS) * 2,
  !anyNA(official$forecast),
  !anyNA(official$actual),
  all(official$selected_p[official$model == "AR_BIC"] %in% 1:P_MAX)
)

table(official$model, official$h)
range(official$target_date)

library(dplyr)

mse_results <- official %>%
  group_by(model, h) %>%
  summarise(
    MSE = mean((actual - forecast)^2),
    .groups = "drop"
  )

print(mse_results, n = Inf)
