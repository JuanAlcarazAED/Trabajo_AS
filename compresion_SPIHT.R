# No hay ninguna librería que implemente de forma sencilla el algoritmo SPIHT para la compresión
# de imágenes. No solo no existe en R si no que en Python tampoco, las implementaciones existentes
# son muy complejas y para su utilización es necesario hacer uso de herramientas muy complicadas
# tanto, que es más fácil hacer nuestra propia implementación

# ------------------------------------------------------------------------------
# LIBRERÍAS

# Para importar imágenes
library(imager)

# Para las transformadas wavelet
library(waveslim)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# El algoritmo necesita que la imagen sea cuadrada y con dimensiones potencias 
# de 4, entonces, vamos a hacer una función que compruebe si una cierta imagen
# cumple las condiciones, y en caso de no hacerlo la recortará para que sí lo haga

recorte <- function(imagen_grises){
  
  # Tomamos el mínimo de las dos 
  n <- min(dim(imagen_grises)[1], dim(imagen_grises)[2])
  
  imagen_recortada <- imagen_grises[1:(n-n%%4),1:(n-n%%4),1,1]
  return(as.cimg(imagen_recortada))
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Ahora la función dwt.2d devuelve un objeto con las subbandas descompuestas hay
# que ensamblarlas en una sola matriz para poder aplicar el algoritmo

Ensamblar <- function(dwt){
  
  # Debemos saber cuántos niveles o subbandas tiene la matriz transformada, para
  # ello tomamos ventaja de la estructura del objeto que devuelve la función 
  # dwt.2d, en la cual los nombres son de la forma LLi o LHi, así que hacemos una
  # lista con los números i que aparecen en los nombres
  niveles <- as.numeric(unique(substr(names(dwt), start=3, stop=3)))
  
  # Los ordenamos de forma descendente
  niveles <- sort(niveles, decreasing = TRUE)
  
  # Tomamos el nivel más bajo para comenzar a ensamblar
  LL_actual <- dwt[[paste0("LL",niveles[1])]]
  
  for(i in niveles){
    HH <- dwt[[paste0("HH", i)]]
    LH <- dwt[[paste0("LH", i)]]
    HL <- dwt[[paste0("HL", i)]]
    
    parte_superior <- rbind(LL_actual, LH)
    parte_inferior <- rbind(HL, HH)
    
    LL_actual <- cbind(parte_superior, parte_inferior)
  }
  return(LL_actual)
}

Desensamblar <- function(W){
  n <- dim(W)[1]
  i <- 1
  Bandas <- list()

  LL_actual <- W
  
  while (!is.null(n) && n > 1){
    
    Bandas[[paste0("LH",i)]] <- LL_actual[((n/2)+1):n, 1:(n/2)]
    Bandas[[paste0("HL",i)]] <- LL_actual[1:(n/2), ((n/2)+1):n]
    Bandas[[paste0("HH",i)]] <- LL_actual[((n/2)+1):n, ((n/2)+1):n]
    LL_actual <- LL_actual[1:(n/2), 1:(n/2)]
    
    n <- dim(LL_actual)[1]
    i <- i +1
  }
  Bandas[[paste0("LL",i-1)]] <- LL_actual
  return(Bandas)
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
################################################################################

# El siguiente es el algoritmo principal, al cual le pasaremos la ruta de la imagen
# y la cantidad máxima de información que se permite guardar

CodificacionSPITH <- function(img, max_bits){

  # A escala de grises
  img_gray <- grayscale(img)
  # En caso de no der cuadrada hay que recortarla
  img_gray <- recorte(img_gray)
  # Se convierte en matriz
  img_mat <- as.matrix(img_gray)
  
  wavelet <- "haar"
  levels <- log2(nrow(img_mat))
  wt <- dwt.2d(img_mat, wf=wavelet, J=levels)
  
  matriz <- Ensamblar(wt)
  
  cod <- Codificar_Transformada(matriz, max_bits = max_bits, max_iter=1000000)
}
################################################################################
