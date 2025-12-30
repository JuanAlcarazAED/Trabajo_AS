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
  
  # Debemos saber cuantos niveles o subbandas tiene las matriz transformada, para
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

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Vamos a necesitar una función que dado un coeficiente (i,j) nos devuelva una
# lista con sus hijos, en caso de tenerlos

Hijos <- function(i, j, size) {
  # Hijos en la siguiente escala
  # Cada padre tiene 4 hijos: (2i,2j), (2i,2j+1), (2i+1,2j), (2i+1,2j+1)
  if ((2*i+1 < size) & (2*j+1 < size) & (2*i<size) & (2*j<size)){
    c1 <- c(2*i,   2*j)
    c2 <- c(2*i,   2*j+1)
    c3 <- c(2*i+1, 2*j)
    c4 <- c(2*i+1, 2*j+1)
    return(rbind(c1, c2, c3, c4))
  }else{
    return(rbind())
  }
}

# Ahora una función para obtener todos los descendientes
Descendientes <- function(i, j, size) {
  queue <- Hijos(i, j, size)
  
  # Si no tiene hijos, no hay descendientes
  if (is.null(queue) || nrow(queue) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }
  
  all_desc <- queue
  
  while (!is.null(queue) && nrow(queue) > 0) {
    new_queue <- NULL
    
    for (k in 1:nrow(queue)) {
      ch <- Hijos(queue[k,1], queue[k,2], size)
      
      if (!is.null(ch) && nrow(ch) > 0) {
        new_queue <- rbind(new_queue, ch)
      }
    }
    
    queue <- new_queue
    
    if (!is.null(new_queue) && nrow(new_queue) > 0) {
      all_desc <- rbind(all_desc, new_queue)
    }
  }
  
  all_desc
}

# ------------------------------------------------------------------------------
# Para comprobar la significancia
Es_significativo <- function(valor, cota) {
  return(abs(valor) >= cota)
}

Es_conjunto_significativo <- function(W, coords, cota) {
  if (is.null(coords) || nrow(coords) == 0) return(FALSE)
  
  # Se recorre el conjunto de las posiciones dadas, si se encuentra una posición
  # con un elemento significativo se para y se devuelve TRUE, si no FALSE
  for (k in 1:nrow(coords)) {
    if (abs(W[coords[k,1] + 1, coords[k,2] + 1]) >= cota) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# ------------------------------------------------------------------------------
# Para obtener la cota inicial
Cota_inicial <- function(W){
  return(floor(log2(max(abs(W)))))
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
}
################################################################################
