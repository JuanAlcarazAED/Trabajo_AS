if (!require("imager", quietly = TRUE)) {
  install.packages("imager")
  library(imager)
}

if (!require("gsignal", quietly = TRUE)) {
  install.packages("gsignal")
  library(gsignal)
}

compresion_DCT_opt <- function(ruta_imagen, ruta_imagen_comp = "data_comp", nivel_compresion = 0.95) {
  
  img <- load.image(ruta_imagen)
  
  canales <- imsplit(img, "c")
  
  procesar_canal <- function(canal) {
    matriz <- as.matrix(canal)
    
    coeficientes <- gsignal::dct2(matriz)
    
    umbral <- quantile(abs(coeficientes), nivel_compresion)
    coeficientes[abs(coeficientes) < umbral] <- 0
    
    recuperado <- gsignal::idct2(coeficientes)
    
    img_recuperada <- as.cimg(recuperado)
    
    img_recuperada <- (img_recuperada + abs(img_recuperada)) / 2 
    img_recuperada <- (img_recuperada - (img_recuperada - 1) * (img_recuperada > 1))
    
    return(img_recuperada)
  }
  
  canales_comprimidos <- lapply(canales, procesar_canal)
  img_final <- imappend(canales_comprimidos, axis = "c")
  
  nombre_archivo <- tools::file_path_sans_ext(basename(ruta_imagen))
  
  ruta_salida <- file.path(ruta_imagen_comp, paste0(nombre_archivo, "_DCT.jpg"))
  
  save.image(img_final, ruta_salida, quality = 0.8) 
  
  invisible(img_final)
  
  return(img_final)
}