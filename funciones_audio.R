library(av)
library(wavelets)
library(tuneR)
library(RColorBrewer)
# Página que voy a usar de referencia
#https://stackoverflow.com/questions/36756500/visualization-of-wavelets-coefficients-for-different-deconstruction-levels 
# Apuntes del tema 3 de compresión de una onda

# --------------- FUNCIONES GRÁFICAS-----------------------


sep_audio_video<-function(main="data/video.mp4",audio="data/audio.wav",video="data/frames",fps=30){
  # Función para separar imagen de audio
  av_audio_convert(main,audio)
  frames<-av_video_images(main,destdir=video,fps=fps)
}

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
    if (plot == "left") vals$wt_left else vals$wt_right
  } else {
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
    main = "Heatmap coeficientes wavelet (detalle)",
    xlab = "Tiempo", ylab = "Nivel wavelet",
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
  
  # Número real de coeficientes preservados
  k_preserved <- sum(all_coefs_thr != 0)
  
  # --- PLOT ---
  plot(
    k_vals, percent_energy, type = "l",
    col = "darkorange", lwd = 2,
    xlab = "Coeficientes preservados (ordenados por energía)",
    ylab = "Energía capturada (%)",
    main = "Energía capturada vs coeficientes"
  )
  
  # Área bajo la curva hasta k_preserved
  polygon(
    x = c(0, k_vals[1:k_preserved], k_preserved),
    y = c(0, percent_energy[1:k_preserved], 0),
    col = rgb(0.2, 0.4, 0.8, 0.35),
    border = NA
  )
  
  # Línea vertical de referencia
  abline(v = k_preserved, lty = 2, col = "blue")
  
  
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


quantize <- function(x, delta) {
  # Cuantización de la señal para la compresión (función usada en quantize_wt)
  round(x / delta)
}


quantize_wt <- function(wt, delta) {
  wt_q <- wt
  wt_q@W <- lapply(wt@W, quantize, delta = delta)
  wt_q@V[[length(wt@V)]] <- quantize(wt@V[[length(wt@V)]], delta)
  wt_q
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


save_compressed <- function(file, left, right, delta, filter, n.levels) {
  saveRDS(
    list(
      left = left,
      right = right,
      delta = delta,
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


dequantize_wt <- function(wt, delta) {
  wt@W <- lapply(wt@W, function(x) x * delta)
  wt@V[[length(wt@V)]] <- wt@V[[length(wt@V)]] * delta
  wt
}