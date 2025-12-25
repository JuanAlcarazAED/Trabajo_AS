library(av)
library(wavelets)
library(tuneR)
library(RColorBrewer)
# Página que voy a usar de referencia
#https://stackoverflow.com/questions/36756500/visualization-of-wavelets-coefficients-for-different-deconstruction-levels 
# Apuntes del tema 3 de compresión de una onda
# Referencia de cómo se umbraliza y redundancia: https://www.di.ens.fr/~mallat/papiers/WaveletTourChap1-6.pdf

# --------------- FUNCIÓN SEPARAR AUDIO-FRAME DE VIDEO-----------------------


sep_audio_video<-function(main="data/video.mp4",audio="data/audio.wav",video="data/frames",fps=30){
  # Función para separar imagen de audio
  av_audio_convert(main,audio)
  frames<-av_video_images(main,destdir=video,fps=fps)
}
# --------------- FUNCIONES GRÁFICAS-----------------------
expand_to_length <- function(coefs, target_len) {
  coefs <- as.vector(coefs)
  rep_each <- ceiling(target_len / length(coefs))
  expanded <- rep(coefs, each = rep_each)
  expanded[1:target_len]
}

plot_heatmap_wavelet_coef <- function(audio, n.levels, lambda,
                                      plot = "left", threshold = FALSE) {
  vals <- dwt_values(audio, n.levels, lambda)
  
  wt <- if (!threshold) {
    main = "Heatmap coeficientes wavelet (detalle) original"
    if (plot == "left") vals$wt_left else vals$wt_right
    
  } else {
    main = "Heatmap coeficientes wavelet (detalle) umbralado"
    if (plot == "left") vals$thr_left else vals$thr_right
    
  }
  
  L <- length(audio@left)
  
  coef_list <- lapply(1:n.levels, function(i) {
    expanded <- expand_to_length(wt@W[[i]], L)
    log1p(expanded^2)
  })
  
  coef_matrix <- do.call(cbind, coef_list)
  
  image(
    x = 1:L, y = 1:n.levels, z = coef_matrix,
    main = main,
    xlab = "Coeficientes", ylab = "Nivel wavelet",
    yaxt = "n", col = heat.colors(100)
  )
  axis(2, at = 1:n.levels, labels = paste0("d", 1:n.levels))
}



plot_umbralD_energy <- function(audio, n.levels, lambda, i,
                                plot = "left", threshold = TRUE) {
  
  vals <- dwt_values(audio, n.levels, lambda)
  
  wt_orig <- if (plot == "left") vals$wt_left else vals$wt_right
  wt_thr  <- if (plot == "left") vals$thr_left else vals$thr_right
  
  coefs <- coef_to_vector(wt_orig@W[[i]])
  coefs_thr <- coef_to_vector(wt_thr@W[[i]])
  
  idx_keep <- which(coefs_thr != 0)
  
  y_lim <- 1.1 * max(abs(coefs))
  
  plot(
    coefs, type = "l",
    main = paste("Nivel de detalle d", i),
    xlab = "Índice del coeficiente",
    ylab = bquote(w[.(i)]),
    ylim = c(-y_lim, y_lim),
    col = "darkgray"
  )
  
  points(
    idx_keep,
    coefs[idx_keep],
    col = "blue",
    pch = 19,
    cex = 0.6
  )
}



plot_energy_vs_coeff <- function(audio, n.levels, lambda, plot = "left") {
  vals <- dwt_values(audio, n.levels, lambda)
  
  wt  <- if (plot == "left") vals$wt_left  else vals$wt_right
  wtT <- if (plot == "left") vals$thr_left else vals$thr_right
  
  # Energía original
  all_coefs <- unlist(lapply(wt@W, as.vector))
  energies  <- abs(all_coefs)^2
  ord <- order(energies, decreasing = TRUE)
  energies_sorted <- energies[ord]
  
  E_total <- sum(energies_sorted)
  E_cum   <- cumsum(energies_sorted)
  percent_energy <- 100 * E_cum / E_total
  k_vals <- seq_along(percent_energy)
  
  # Energía tras threshold
  all_coefs_thr <- unlist(lapply(wtT@W, as.vector))
  E_after_thr <- sum(abs(all_coefs_thr)^2)
  percent_after_thr <- 100 * E_after_thr / E_total
  
  k_preserved <- sum(all_coefs_thr != 0)
  
  plot(
    k_vals, percent_energy, type = "l",
    col = "darkorange", lwd = 2,
    xlab = "Coeficientes preservados (ordenados por energía)",
    ylab = "Energía capturada (%)",
    main = "Energía capturada vs coeficientes"
  )
  
  
  # Línea horizontal de referencia
  abline(h = percent_after_thr, lty = 5, col = "red")
  
  grid()
}









# -------------- FUNCIONES IMPORTANTES --------------------
hard_threshold_energy <- function(coefs, lambda) {
  # Función para aplicar el umbral de energía
  # lambda=0 => No se comprime la señal
  # lambda=1 => se comprime el 100% señal
  power <- abs(coefs)^2
  total_energy <- sum(power)
  ord <- order(power)
  cum_energy <- cumsum(power[ord])
  idx <- ord[cum_energy <= lambda * total_energy]
  coefs[idx] <- 0
  coefs
}

dwt_values <- function(audio, n.levels, lambda) {
  # Hacer la dwt sobre los dos canales de la señal y aplicar el threshold
  signal_left  <- as.numeric(audio@left)
  signal_right <- as.numeric(audio@right)
  
  wt_left  <- dwt(signal_left,  filter = "la8", n.levels = n.levels)
  wt_right <- dwt(signal_right, filter = "la8", n.levels = n.levels)
  
  thr_left  <- wt_left
  thr_right <- wt_right
  
  thr_left@W  <- lapply(wt_left@W,  hard_threshold_energy, lambda = lambda)
  thr_right@W <- lapply(wt_right@W, hard_threshold_energy, lambda = lambda)
  
  list(
    wt_left  = wt_left,
    wt_right = wt_right,
    thr_left = thr_left,
    thr_right = thr_right
  )
}


coef_to_vector <- function(x) {
  as.vector(x)
}


rle_encode <- function(x) {
  x <- as.vector(x)
  r <- rle(x)
  list(values = r$values, lengths = r$lengths)
}

encode_wt <- function(wt_q) {
  list(
    W = lapply(wt_q@W, rle_encode),
    V = rle_encode(wt_q@V[[length(wt_q@V)]])
  )
}


save_compressed <- function(file, left, right, filter, n.levels) {
  saveRDS(
    list(
      left = left,
      right = right,
      filter = filter,
      n.levels = n.levels
    ),
    file = file
  )
  file.info(file)$size
}


rle_decode <- function(enc) {
  inverse.rle(list(values = enc$values, lengths = enc$lengths))
}


decode_wt <- function(enc, template_wt) {
  wt <- template_wt
  
  wt@W <- Map(
    function(e, ref) {
      v <- rle_decode(e)
      matrix(v, nrow = length(v), ncol = ncol(ref))
    },
    enc$W,
    template_wt@W
  )
  
  v <- rle_decode(enc$V)
  wt@V[[length(wt@V)]] <- matrix(
    v,
    nrow = length(v),
    ncol = ncol(template_wt@V[[length(template_wt@V)]])
  )
  
  wt
}


size_vs_lambda <- function(audio, n.levels, lambda_vals,
                           filter = "la8", tmp_file = "tmp_comp.rds") {
  
  sizes <- numeric(length(lambda_vals))
  
  for (i in seq_along(lambda_vals)) {
    lambda <- lambda_vals[i]
    
    # Wavelet + threshold
    vals <- dwt_values(audio, n.levels, lambda)
    
    wt_left  <- vals[[3]]  # thresholded left delta
    wt_right <- vals[[4]]  # thresholded right
    

    
    # Guardar archivo comprimido
    saveRDS(
      list(
        left   = wt_left,
        right  = wt_right,
        lambda = lambda,
        filter = filter,
        levels = n.levels
      ),
      file = tmp_file
    )
    
    # Tamaño real en bytes
    sizes[i] <- file.info(tmp_file)$size
    
    # Limpieza
    unlink(tmp_file)
  }
  
  data.frame(
    lambda = lambda_vals,
    size_bytes = sizes
  )
}

entropy_redundancy <-function(x){
  nbins<-1000
  energia <- x^2
  
  min_e <- min(energia, na.rm = TRUE)
  max_e <- max(energia, na.rm = TRUE)
  h <- hist(energia, breaks = seq(min_e, max_e, length.out = nbins + 1), plot = FALSE)

  p <- h$counts / sum(h$counts)
  p <- p[p > 0] # Solo usamos los bins que tienen datos
  H<--sum(p*log2(p))
  Hmax<-log2(nbins)

  R<-1-H/Hmax
  
  list(entropy=H,Hmax=Hmax,redundancy=R)
}


entropy_redundancy_wt <- function(wt_q) {
  coefs <- c(
    unlist(wt_q@W),
    wt_q@V[[length(wt_q@V)]]
  )
  
  entropy_redundancy(coefs)
}

redundancy_vs_lambda <- function(audio, n.levels, lambda_vals,
                                 plot = "left") {
  
  H <- numeric(length(lambda_vals))
  R <- numeric(length(lambda_vals))
  
  for (i in seq_along(lambda_vals)) {
    lambda <- lambda_vals[i]
    
    vals <- dwt_values(audio, n.levels, lambda)
    wt_thr <- if (plot == "left") vals[[3]] else vals[[4]]
    
    wt_q <- wt_thr
  
    info <- entropy_redundancy_wt(wt_q)
    
    H[i] <- info$entropy
    R[i] <- info$redundancy
  }
  
  data.frame(
    lambda = lambda_vals,
    entropy = H,
    redundancy = R
  )
}

