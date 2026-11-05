# Adapted from `func-lasso.R` in Lecture 4

#### CHOOSING LAMBDA BY ROLLING-ORIGIN CV ####
# helper function to replace ic.glmnet

cv_indices <- function(n, h, V = 50) {
  # validation rows are the last V rows of the design matrix
  val.rows <- (n - V + 1):n
  
  # fixed training length = number of training rows in the first fold
  L <- n - V + 1 - h
  stopifnot(L >= 1)
  
  # last training row for validation row j = j - h
  # first training row shifts forward with j (preserve shape)
  last.train <- val.rows - h
  first.train <- last.train - L + 1
  
  # return indices for each fold
  return(list("val" = val.rows, "first" = first.train, "last" = last.train, "L" = L))
}

cv_lambda <- function(X, y, h, V = 50, lambda.grid, alpha = 1, na.handling = "fill") {
  # build fold indices
  idx <- cv_indices(nrow(X), h, V)
  
  # blank matrix of squared errors
  # one row per validation point, one column per lambda
  save.error <- matrix(NA, V, length(lambda.grid))
  
  for (i in 1:V) {
    # training block for validation row j (fixed length L)
    train <- idx$first[i]:idx$last[i]
    j <- idx$val[i]
    
    # fit LASSO on the training block using shared lambda grid
    fit <- glmnet(X[train, , drop = FALSE], y[train], alpha = alpha, lambda = lambda.grid)
    
    # predict row j at every lambda (one column per lambda)
    pred <- as.numeric(predict(fit, newx = X[j, , drop = FALSE]))
    
    # fill tail with NA if glmnet terminates early
    pred.full <- rep(NA, length(lambda.grid))
    pred.full[seq_along(pred)] <- pred
    
    # fill missing columns with fold's last fitted prediction
    if(na.handling == "fill") {
      pred.full[is.na(pred.full)] <- pred[length(pred)]
    }
    
    # squared error against realised target
    save.error[i, ] <- (y[j] - pred.full)^2
  }
  # mean squared error for each lambda over the V validation points
  if(na.handling == "rm") {
    mse <- colMeans(save.error, na.rm = TRUE)
  } else {
    mse <- colMeans(save.error)
  }
    
  # return position of minimum
  return(which.min(mse))
}

#### FORMING FORECASTS USING LASSO ####

runlasso <- function(Y, indice, lag, alpha = 1, V = 50, na.handling = "fill") {
  # dropping last h-1 rows so window ends at forecast origin
  Y <- Y[1:(nrow(Y) - lag + 1), ]
  
  # predictor names come from the column names, so they must exist
  stopifnot(!is.null(colnames(Y)))
  
  # compute principal components to add as predictors
  comp <- prcomp(scale(Y, scale = FALSE))
  
  # augment predictors by the first 4 PCs
  Y2 <- cbind(Y, comp$x[,1:4])
  
  # create 4 lags + forecast horizon shift (=lag option)
  aux <- embed (Y2, 4 + lag)
  
  # align Y for missing data due to lags
  y <- aux[,indice]
  
  # lags of Y (predictors) corresponding to forecast horizon
  X <- aux[,-c(1:(ncol(Y2)*lag))]
  
  # name predictors by series and lag relative to the forecast origin (L0 = origin month)
  xnames <- paste0(rep(colnames(Y2), 4), "_L", rep(0:3, each = ncol(Y2)))
  stopifnot(length(xnames) == ncol(X), anyDuplicated(xnames) == 0)
  colnames(X) <- xnames
  
  # retrieve last 4 observations of the window as a one-row matrix
  X.out <- matrix(tail(aux,1)[1:ncol(X)], nrow = 1, dimnames = list(NULL, xnames))
  
  # fit LASSO selected on CV
  # lambda vector is shared grid for every CV fold
  full <- glmnet(X, y, alpha = alpha)
  
  # rolling-origin CV picks the position of lambda.min in the grid
  pos <- cv_lambda(X, y, h = lag, V = V, lambda.grid = full$lambda, 
                   alpha = alpha, na.handling = na.handling)
  lambda.min <- full$lambda[pos]
  
  # named coefficient vector at lambda.min (intercept first)
  beta <- as.matrix(coef(full, s = lambda.min))[, 1]
  
  # LASSO forecast from same fit
  pred.lasso <- as.numeric(predict(full, newx = X.out, s = lambda.min))
  
  # LASSO support (excluding intercept)
  support <- which(beta[-1] != 0)
  
  # check support:
  # forecast training mean if empty support, use LASSO if support above n/2
  empty <- length(support) == 0
  large <- length(support) > nrow(X) / 2
  
  if (empty) {
    # empty support, forecast training mean
    pred.post <- mean(y)
  } else if (large) {
    # large support, use LASSO forecast
    pred.post <- pred.lasso
  } else {
    # run OLS on the LASSO support at the same origin
    ols <- lm.fit(cbind(1, X[, support, drop = FALSE]), y)
    b <- ols$coefficients
    
    b[is.na(b)] <- 0
    pred.post <- sum(b * c(1, X.out[support]))
  }
  
  # return coefficients, both forecasts, chosen lambda, and support diagnostics
  return(list("coef" = beta,
              "pred" = c("lasso" = pred.lasso, "post" = pred.post),
              "lambda" = lambda.min,
              "support" = c("size" = length(support), "empty" = empty, "large" = large)))
}

#### PRODUCING SERIES OF H-STEP LASSO FORECASTS USING ROLLING WINDOW ####

lasso.rolling.window <- function( Y, nprev, dates, indice = 1, lag = 1, alpha = 1, V = 50,
                                na.handling = "fill") {
  
  # dates must be passed explicitly and match rows of Y
  stopifnot(!is.null(dates), length(dates) == nrow(Y))
  
  # blank matrix for coefficients at each iteration
  # intercept + 4 lags of (series + 4 PCs) = 1 + 4*(406 + 4) = 1641
  save.coef <- matrix(NA, nprev, 1 + 4*(ncol(Y) + 4))
  
  # blank forecasts (one column per method)
  save.pred <- matrix(NA, nprev, 2, dimnames = list(NULL, c("lasso", "post")))
  
  # blank support diagnostics
  save.support <- matrix(NA, nprev, 3, dimnames = list(NULL, c("size", "empty", "large")))
  
  # blank chose lambda
  save.lambda <- rep(NA, nprev)
  
  for(i in nprev:1) {
    # define estimation window
    Y.window <- Y[(1 + nprev - i):(nrow(Y) - i),]
    
    # call function to fit LASSO and Post-LASSO
    lasso <- runlasso(Y.window, indice, lag, alpha, V, na.handling)
    
    # save estimated coefficients
    save.coef[(1 + nprev - i),] <- lasso$coef
    
    #save both forecast
    save.pred[(1 + nprev - i),] <- lasso$pred
    
    # save support diagnostics and lambda
    save.support[(1 + nprev - i),] <- lasso$support
    save.lambda[(1 + nprev - i)] <- lasso$lambda
    
    # # display iteration number
    # cat("iteration", (1 + nprev - i), "\n")
  }
  # name coefficient columns (intercept + predictors) for the selected-variable analysis
  colnames(save.coef) <- names(lasso$coef)
  
  # Helpful diagnostics
  # get actual values
  real <- Y[,indice]
  
  # target rows are the last nprev rows at every horizon
  target.rows <- (nrow(Y) - nprev + 1):nrow(Y)
  
  # each origin is h rows before its target
  origin.rows <- target.rows - lag
  
  # count how often the support exceeded n/2 or was empty
  n.large <- sum(save.support[, "large"])
  n.empty <- sum(save.support[, "empty"])
  cat("support above n/2:", n.large, "of", nprev, "| empty:", n.empty, "of", nprev, "\n")
  
  # # plot preds vs actual
  # plot(real, type = "l")
  # lines(c(rep(NA, length(real) - nprev), save.pred), col = "red")
  
  # dates for each forecast
  target_date <- dates[target.rows]
  origin_date <- dates[origin.rows]
  
  # realised values at the target dates
  actual <- real[target.rows]
  
  # calculate RMSE, MAE for each method
  err <- actual - save.pred
  rmse <- sqrt(colMeans(err^2))
  mae <- colMeans(abs(err))
  
  # stack errors in a vector
  errors <- c("RMSE" = rmse, "MAE" = mae)
  
  # return forecasts, actual values, dates, coefficients, lambdas, support info, RMSE, MAE for the period
  return(list("pred" = save.pred, "actual" = actual,
              "target_date" = target_date, "origin_date" = origin_date,
              "coef" = save.coef, "lambda" = save.lambda, "support" = save.support,
              "errors" = errors))
}
