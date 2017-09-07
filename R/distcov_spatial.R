distcov_spatial <- function(X, bias_corr = TRUE, metr.X = "euclidean", isotropic = FALSE) {
    coords <- coordinates(X)
    names_X <- names(X)
    if (length(names_X) != 1)
        stop ("Function works currently only for univariate fields")
    if (isotropic) {
        distances <- gDistance(coords)
    } else {
        lags <- rbind(calculate_lags(coords), c(0, 0))
        n <- nrow(lags)
        distcov_spdf <- SpatialPointsDataFrame(coords = lags, data = data.frame(distcov = rep(NA, n)))
        distcov_spdf_mirror <- SpatialPointsDataFrame(coords = -lags[-n, ], data = data.frame(distcov = rep(NA, n - 1)))
        while (n > 0) {
            coords_copy <- data.table::copy(coords)
            idx_common <- match_coords(coords_copy, lags[n, ])
            sampleX <- X[[names_X]][idx_common[, 2]]
            sampleY <- X[[names_X]][idx_common[, 1]]
            bias_correction <- ifelse(length(sampleX) < 4, FALSE, bias_corr)
            distcov_spdf$distcov[n] <- ifelse(is.null(sampleX), NA, distcov(sampleX, sampleY, bias_corr = bias_correction))
            n <- n - 1
        }
        distcov_spdf_mirror$distcov <- distcov_spdf$distcov[1:nrow(distcov_spdf_mirror)]
    }
    return (rbind(distcov_spdf, distcov_spdf_mirror))
}

distcorr_spatial <- function(X, bias_corr = TRUE, metr.X = "euclidean", isotropic = FALSE) {
  coords <- coordinates(X)
  names_X <- names(X)
  if (length(names_X) != 1)
    stop ("Function works currently only for univariate fields")
  if (isotropic) {
    distances <- coords %>%
      SpatialPoints() %>%
      rgeos::gDistance(X, byid = TRUE)
    
  } else {
    lags <- rbind(calculate_lags(coords), c(0, 0))
    n <- nrow(lags)
    distcorr_spdf <- SpatialPointsDataFrame(coords = lags, data = data.frame(distcorr = rep(NA, n)))
    distcorr_spdf_mirror <- SpatialPointsDataFrame(coords = -lags[-n, ], data = data.frame(distcorr = rep(NA, n - 1)))
    while (n > 0) {
      coords_copy <- data.table::copy(coords)
      idx_common <- match_coords(coords_copy, lags[n, ])
      sampleX <- X[[names_X]][idx_common[, 2]]
      sampleY <- X[[names_X]][idx_common[, 1]]
      bias_correction <- ifelse(length(sampleX) < 4, FALSE, bias_corr)
      distcorr_spdf$distcorr[n] <- ifelse(is.null(sampleX), NA, distcorr(sampleX, sampleY, bias_corr = bias_correction))
      n <- n - 1
    }
    distcorr_spdf_mirror$distcorr <- distcorr_spdf$distcorr[1:nrow(distcorr_spdf_mirror)]
  }
  return (rbind(distcorr_spdf, distcorr_spdf_mirror))
}

calculate_lags <- function(coord_matrix) {
    lag <- 1
    n <- nrow(coord_matrix)
    lag_matrix <- matrix(ncol = 2, nrow = n * (n - 1) / 2)
    start_idx <- 1
    end_idx <- n - 1
    while (n > 1) {
        lag_matrix[start_idx:end_idx, ] <- diff(coord_matrix, lag = lag)
        lag <- lag + 1
        start_idx <- start_idx + n - 1
        end_idx <- end_idx + n - 2
        n <- n - 1
    }
    return (unique(lag_matrix))
}


