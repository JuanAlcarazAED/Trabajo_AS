library(av)
library(wavelets)
library(tuneR)
library(RColorBrewer)
# Página que voy a usar de referencia
#https://stackoverflow.com/questions/36756500/visualization-of-wavelets-coefficients-for-different-deconstruction-levels 
# Apuntes del tema 3 de compresión de una onda

# Funciones importantes


sep_audio_video<-function(main="data/video.mp4",audio="data/audio.wav",video="data/frames",fps=30){
  # Función para separar imagen de audio
  av_audio_convert(main,audio)
  frames<-av_video_images(main,destdir=video,fps=fps)
}

expand_to_length <- function(coefs, target_len) {
  rep_each <- ceiling(target_len / length(coefs))
  expanded <- rep(coefs, each = rep_each)
  expanded[1:target_len]
}


hard_threshold_energy <- function(coefs, lambda) {
  power <- abs(coefs)^2
  total_energy <- sum(power)
  
  # Ordenar coeficientes por energía creciente
  ord <- order(power)
  cum_energy <- cumsum(power[ord])
  idx <- ord[cum_energy < lambda * total_energy]
  coefs[idx] <- 0
  
  return(coefs)
}


dwt_values<-function(audio,n.levels,lambda){
  
  signal_left<-as.numeric(audio@left)
  signal_rigth<-as.numeric(audio@right)
  
  wt_left <- dwt(signal_left, filter = "la8", n.levels = n.levels)
  wt_right <- dwt(signal_rigth, filter = "la8", n.levels = n.levels)
  
  wt_thresholded_left <- wt_left 
  wt_thresholded_right <- wt_right 
  
  wt_thresholded_left@W  <- lapply(wt_left@W,  hard_threshold_energy, lambda=lambda)
  wt_thresholded_right@W <- lapply(wt_right@W, hard_threshold_energy, lambda=lambda)
  return(list(
    wt_left  = wt_left,
    wt_right = wt_right,
    thresholded_left  = wt_thresholded_left,
    thresholded_right = wt_thresholded_right
  ))
}


plot_heatmap_wavelet_coef<-function(audio,n.levels,lambda,plot="left",threshold=F){
  vals<-dwt_values(audio,n.levels,lambda)
  if (threshold==F){
    if (plot=="left"){
      wt<-vals[[1]]
    }
    else{
      wt<-vals[[2]]
    }
  }
  else{
    if (plot=="left"){
      wt<-vals[[3]]
    }
    else{
      wt<-vals[[4]]
    }
    
  }
  L<- length(audio@left)
  
  coef_list <- lapply(1:n.levels, function(i) {
    expanded <- expand_to_length(wt@W[[i]], L)
    log(1+abs(expanded)^2)
  })
  
  coef_matrix_raw <- do.call(cbind, coef_list)
  coef_matrix <- (coef_matrix_raw[, 1:n.levels])
  
  colores <- brewer.pal(9, "YlOrBr")
  
  image(x=1:L,y=1:n.levels,z=(coef_matrix),main="Wavelet coeficientes de detalle heatmap",xlab="Tiempo",ylab="Nivel wavelet",yaxt="n")
  axis(side=2,at=1:n.levels,labels = paste0("d", 1:n.levels))
}


plot_umbralD<-function(audio,n.levels,lambda,i,plot="left",threshold=F){
  # Coeficientes del nivel i
  vals<-dwt_values(audio,n.levels,lambda)
  if (threshold==F){
    if (plot=="left"){
      wt<-vals[[1]]
    }
    else{
      wt<-vals[[2]]
    }
  }
  else{
    if (plot=="left"){
      wt<-vals[[3]]
    }
    else{
      wt<-vals[[4]]
    }
    
  }
  coefs <- wt@W[[i]]
  
  # Energía logarítmica (solo la calculamos una vez)
  power <- abs(coefs)**2
  # Umbral relativo
  umbral <- lambda * max(power)
  
  # Índices que superan el umbral
  idx_keep <- which(power >= umbral)
  
  # Límite vertical para el plot
  y_lim <- 1.1 * max(abs(coefs))
  
  # Grafica
  plot(
    coefs, type = 'l',
    main = paste("Nivel de Detalle d", i),
    xlab = "Índice del Coeficiente",
    ylab = bquote(w[.(i)]),
    ylim = c(-y_lim, y_lim),
    col = "darkgray"
  )
  
  # Añadimos puntos thresholded
  points(
    x = idx_keep,
    y = coefs[idx_keep],
    col = "blue",
    pch = 19,
    cex = 0.6
  )
}


plot_umbralD_energy <- function(audio, n.levels, perc, i, plot="left", threshold=FALSE) {
  vals <- dwt_values(audio, n.levels, perc)
  
  if (!threshold) {
    wt <- if (plot=="left") vals[[1]] else vals[[2]]
  } else {
    wt <- if (plot=="left") vals[[3]] else vals[[4]]
  }
  
  # Coefs originales (para graficar)
  wt_original <- if (plot=="left") vals[[1]] else vals[[2]]
  coefs <- wt_original@W[[i]]
  
  # Coefs thresholded
  coefs_thr <- wt@W[[i]]
  
  # Determinar índices mantenidos
  idx_keep <- which(coefs_thr != 0)
  
  y_lim <- 1.1 * max(abs(coefs))
  
  plot(
    coefs, type="l",
    main=paste("Nivel de Detalle d", i),
    xlab="Índice del Coeficiente",
    ylab=bquote(w[.(i)]),
    ylim=c(-y_lim, y_lim),
    col="darkgray"
  )
  
  points(
    x=idx_keep,
    y=coefs[idx_keep],
    col="blue",
    pch=19,
    cex=0.6
  )
}
plot_energy_vs_coeff <- function(audio, n.levels, lambda, plot="left") {
  vals <- dwt_values(audio, n.levels, lambda)
  
  wt  <- if (plot=="left") vals$wt_left  else vals$wt_right
  wtT <- if (plot=="left") vals$thresholded_left else vals$thresholded_right
  
  # Energía original
  all_coefs <- unlist(wt@W)
  energies  <- abs(all_coefs)^2
  ord <- order(energies, decreasing = TRUE)
  energies_sorted <- energies[ord]
  E_total <- sum(energies_sorted)
  E_cum   <- cumsum(energies_sorted)
  percent_energy <- 100 * E_cum / E_total
  k_vals <- seq_along(all_coefs)
  
  # Energía tras threshold real
  all_coefs_thr <- unlist(wtT@W)
  E_after_thr   <- sum(abs(all_coefs_thr)^2)
  percent_after_thr <- 100 * E_after_thr / E_total
  
  # Número real de coeficientes preservados:
  coefs_keep <- sum(all_coefs_thr != 0)
  
  # Dibujar curva
  plot(k_vals, percent_energy, type="l",
       col="darkorange", lwd=2,
       xlab="Coeficientes preservados (ordenados por energía)",
       ylab="Energía capturada (%)",
       main="Curva energía capturada vs coeficientes preservados")
  
  # Para situar el punto sobre la curva: encontrar k tal que percent_energy[k] >= percent_after_thr
    k_on_curve <- which(percent_energy >= percent_after_thr)[1]
    if (is.na(k_on_curve)) { # por si ocurre redondeo y percent_after_thr == 100
      k_on_curve <- length(k_vals)
    }
    # Punto sobre la curva
    points(k_on_curve, percent_energy[k_on_curve], pch=19, col="blue", cex=1.3)
    text(k_on_curve, percent_energy[k_on_curve],
         labels = sprintf("  %.2f%% energía", percent_after_thr),
         pos = 4, col="blue")

  
  grid()
}


compress_audio_wavelet_to_aac <- function(audio,
                                          n.levels,
                                          lambda,
                                          ruta_aac = "data_comp/audio_comp.aac",
                                          ruta_wav_tmp = "data_comp/tmp_audio.wav",
                                          borrar_tmp = TRUE) {
  
  # 1. Obtener coeficientes wavelet y umbralizarlos
  vals <- dwt_values(audio, n.levels, lambda)
  
  wt_left  <- vals[[1]]
  wt_right <- vals[[2]]
  thr_left <- vals[[3]]
  thr_right <- vals[[4]]
  
  
  # 2. Reconstrucción
  signal_left  <- idwt(thr_left)
  signal_right <- idwt(thr_right)
  sd(signal_left)
  eps <- 1e-5
  signal_left[abs(signal_left) < eps] <- 0
  signal_right[abs(signal_right) < eps] <- 0
  
  # 4. Crear objeto Wave
  audio_wav <- Wave(
    left = signal_left,
    right = signal_right,
    samp.rate = audio@samp.rate,
    bit = audio@bit
  )
  
  # 5. Guardar WAV temporal
  writeWave(audio_wav, ruta_wav_tmp)
  
  # 6. Convertir a AAC
  av::av_audio_convert(ruta_wav_tmp, ruta_aac)
  
  # 7. Borrar WAV temporal si se pide
  if (borrar_tmp && file.exists(ruta_wav_tmp)) {
    file.remove(ruta_wav_tmp)
  }
  
  cat("Archivo AAC creado en:", ruta_aac, "\n")
  invisible(ruta_aac)
  return(vals)
}

