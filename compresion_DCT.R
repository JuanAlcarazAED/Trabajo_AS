if (!require("imager", quietly = TRUE)) {
  install.packages("imager")
  library(imager)
}

crear_matriz_dct <- function(N) {
  k <- 0:(N-1)
  n <- 0:(N-1)
  
  matriz <- cos(outer(k, n, function(k, n) (pi/N) * (n + 0.5) * k))
  
  matriz[1, ] <- matriz[1, ] * (1 / sqrt(N))
  matriz[2:N, ] <- matriz[2:N, ] * sqrt(2 / N)
  return(matriz)
}

aplicar_dct2d <- function(img_matriz) {
  img_matriz <- as.matrix(img_matriz)
  filas <- nrow(img_matriz)
  cols <- ncol(img_matriz) 
  
  D_filas <- crear_matriz_dct(filas)
  D_cols <- crear_matriz_dct(cols)
  
  return(D_filas %*% img_matriz %*% t(D_cols))
}

aplicar_idct2d <- function(coeficientes) {
  filas <- nrow(coeficientes)
  cols <- ncol(coeficientes)
  
  D_filas <- crear_matriz_dct(filas)
  D_cols <- crear_matriz_dct(cols)
  
  return(t(D_filas) %*% coeficientes %*% D_cols)
}

# se calcula lo de dct "manualmente" porque supuestamente la última versión de pracma no tenía las
# funciones dct e idct, y usando dtt me daba problemas


compresion_DCT <- function(ruta_imagen, nivel_compresion = 0.95) {
  
  img <- load.image(ruta_imagen)
  
  canales <- imsplit(img, "c")
  
  procesar_canal <- function(canal) {
    matriz <- as.matrix(canal)
    
    coeficientes <- aplicar_dct2d(matriz)
    
    umbral <- quantile(abs(coeficientes), nivel_compresion)
    coeficientes[abs(coeficientes) < umbral] <- 0
    
    recuperado <- aplicar_idct2d(coeficientes)
    
    img_recuperada <- as.cimg(recuperado)
    
    img_recuperada <- (img_recuperada + abs(img_recuperada)) / 2 
    img_recuperada <- (img_recuperada - (img_recuperada - 1) * (img_recuperada > 1))
    
    return(img_recuperada)
  }

  canales_comprimidos <- lapply(canales, procesar_canal)
  img_final <- imappend(canales_comprimidos, axis = "c")
  
  nombre_archivo <- tools::file_path_sans_ext(basename(ruta_imagen))
  
  ruta_salida <- file.path("data_comp", paste0(nombre_archivo, "_DCT.jpg"))
  
  save.image(img_final, ruta_salida, quality = 0.8) 

  invisible(img_final)
}