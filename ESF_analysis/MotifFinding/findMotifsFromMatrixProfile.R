min_MP_idx <- function (mp, n_dim = NULL, valid = TRUE)  {
  n_dim <- 1
  mp_size <- length(mp$matrix_profile)
  min <- which(mp$matrix_profile == min(mp$matrix_profile) )
  if (any(min == 1) && any(is.infinite(.mp$mp[1, (min == 1)]))) {
    return(NA)
  }
  nn_min <- NULL
  for (i in seq_len(n_dim)) {
    nn_min <- c(nn_min, mp$profile_index[min[i]])
  }
  if (valid) {
    if (all(nn_min > 0 & nn_min <= mp_size) && all(!is.infinite(diag(mp$matrix_profile[nn_min], names = FALSE)))) {
      return(cbind(min, nn_min, deparse.level = 0))
    }
    for (i in seq_len(n_dim)) {
      mp$matrix_profile[min[i]] <- Inf
    }
    stop <- FALSE
    while (!stop) {
      min <- which(mp$matrix_profile == min(mp$matrix_profile))
      if (any(min == 1) && any(is.infinite(mp$matrix_profile[1, (min == 1)]))) {
        stop <- TRUE
      }
      else {
        nn_min <- NULL
        for (i in seq_len(n_dim)) {
          nn_min <- c(nn_min, mp$profile_index[min[i]])
        }
        if (all(nn_min > 0 & nn_min <= mp_size) && all(!is.infinite(diag(mp$matrix_profile[nn_min], 
                                                                         names = FALSE)))) {
          return(cbind(min, nn_min, deparse.level = 0))
        }
        else {
          for (i in seq_len(n_dim)) {
            mp$matrix_profile[min[i]] <- Inf
          }
        }
      }
    }
    return(NA)
  }
  else {
    return(cbind(min, nn_min, deparse.level = 0))
  }
}





findMotifs <- function(MP, data, w_size, n_motifs = 3, n_neighbors = 10, radius = 3, exclusion_zone = NULL) {
  
  # transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data_len <- nrow(data)
    data_dim <- ncol(data)
  } else if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)
    
    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_len <- length(data)
    data_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("`data` must be `matrix`, `data.frame`, `vector` or `list`.")
  }
  
  matrix_profile <- MP # keep mp intact
  matrix_profile_size <- length(matrix_profile$matrix_profile)
  motif_idxs <- list(motifs = list(NULL), neighbors = list(NULL), windows = list(NULL))

  exclusion_zone <- round(w_size * MP$ez + 1e-8)
  window <- w_size
  e_zone <- exclusion_zone
  nn <- NULL
  
  for (i in seq_len(n_motifs)) {
    idxs <- min_MP_idx(matrix_profile)
    
    if (is.na(idxs[1])) {
      break
    }
    
    min_idx <- idxs[1]
    motif_distance <- matrix_profile$matrix_profile[min_idx]
    motif_idxs[[1]][[i]] <- sort(idxs)
    motif_idx <- motif_idxs[[1L]][[i]][1]
    
    # query using the motif to find its neighbors
    nn <- dist_profile(data, data, nn, window_size = window, index = min_idx)
    
    distance_profile <- nn$distance_profile
    
    distance_profile[distance_profile > (motif_distance * radius)^2] <- Inf
    motif_zone_start <- pmax(1, motif_idx - e_zone)
    motif_zone_end <- pmin(matrix_profile_size, motif_idx + e_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    motif_idx <- motif_idxs[[1]][[i]][2]
    motif_zone_start <- pmax(1, motif_idx - e_zone)
    motif_zone_end <- pmin(matrix_profile_size, motif_idx + e_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    st <- sort(distance_profile, index.return = TRUE)
    distance_order <- st$x
    distance_idx_order <- st$ix
    
    motif_neighbor <- vector(mode = "numeric")
    
    for (j in seq_len(n_neighbors)) {
      if (is.infinite(distance_order[1]) || length(distance_order) < j) {
        break
      }
      motif_neighbor[j] <- distance_idx_order[1]
      distance_order <- distance_order[2:length(distance_order)]
      distance_idx_order <- distance_idx_order[2:length(distance_idx_order)]
      distance_order <- distance_order[!(abs(distance_idx_order - motif_neighbor[j]) < e_zone)]
      distance_idx_order <- distance_idx_order[!(abs(distance_idx_order - motif_neighbor[j]) < e_zone)]
    }
    
    motif_neighbor <- motif_neighbor[motif_neighbor != 0]
    motif_idxs[[2]][[i]] <- motif_neighbor
    motif_idxs[[3]][[i]] <- window
    
    remove_idx <- c(motif_idxs[[1]][[i]], motif_idxs[[2]][[i]])
    
    for (j in seq_len(length(remove_idx))) {
      remove_zone_start <- max(1, remove_idx[j] - e_zone)
      remove_zone_end <- min(matrix_profile_size, remove_idx[j] + e_zone)
      matrix_profile$matrix_profile[remove_zone_start:remove_zone_end] <- Inf
    }
  }
  
  if (is.null(motif_idxs[[1]][[1]])) {
    message("No valid motif found.")
    MP <- remove_class(MP, "Motif")
    return(MP)
  }
  
  motifs <- list(motif_idx = motif_idxs[[1]], motif_neighbor = motif_idxs[[2]], motif_window = motif_idxs[[3]])
  return(motifs)
}




for(m in seq_along(motifs$motif_idx)){
  cat(m, "\n")
  for(ml in seq_along(motifs$motif_idx[[m]]) ){
    cat(m, " - ", ml, ": ", motifs$motif_idx[[m]][ml], "\n")
    abline(v = motifs$motif_idx[[m]][ml], col=m)
  }
  for(nl in seq_along(motifs$motif_neighbor[[m]])){
    cat(m, " - ", ml, ": ", motifs$motif_neighbor[[m]][nl], "\n")
    abline(v = motifs$motif_neighbor[[m]][nl], col=m, lty=2)
  }
}

