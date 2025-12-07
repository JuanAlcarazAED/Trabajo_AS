if (!require("imager", quietly = TRUE)) {
  message("Instalando paquete 'imager'...")
  install.packages("imager")
  library(imager)
}

if (!require("wavethresh", quietly = TRUE)) {
  message("Instalando paquete 'wavethresh'...")
  install.packages("wavethresh")
  library(wavethresh)
}

obtener_potencia_2 <- function(x) {
  return(2^round(log2(x)))
}

# Esta función existe porque se necesita que la imagen que se ingresa sea potencia de 2 (256, 512, 1024...)
# porque idwt funciona con imágenes cuadradas

compresion_IDWT <- function(ruta_imagen, politica_umbral = "universal") {
  
  img <- load.image(ruta_imagen)
  
  nuevo_tamaño <- obtener_potencia_2(max(width(img), height(img))) #se le asigna el nuevo tamaño a la imagen
  
  img <- resize(img, size_x = nuevo_tamaño, size_y = nuevo_tamaño)

  canales <- imsplit(img, "c")

  procesar_canal <- function(canal) {

    matriz <- as.matrix(canal)

    lwd <- imwd(matriz)
    
    lwdT <- threshold(lwd, policy = politica_umbral, type = "hard")

# political_umbral se define como "universal" que lo que hace es que dependiendo de la imagen
# si tiene mucho ruido o no elimina/comprime más o menos la imagen
# en el método de DCt está establecido en 95%
    
    recuperado <- imwr(lwdT)
    
    img_rec <- as.cimg(recuperado)
    img_rec <- (img_rec + abs(img_rec)) / 2 
    img_rec <- (img_rec - (img_rec - 1) * (img_rec > 1))
    
    return(img_rec)
  }
  
  canales_comprimidos <- lapply(canales, procesar_canal)
  img_final <- imappend(canales_comprimidos, axis = "c")
  
  nombre_archivo <- tools::file_path_sans_ext(basename(ruta_imagen))
  ruta_salida <- file.path("data_comp", paste0(nombre_archivo, "_IDWT.jpg"))
  
  save.image(img_final, ruta_salida, quality = 0.8)
  
  invisible(img_final)
}