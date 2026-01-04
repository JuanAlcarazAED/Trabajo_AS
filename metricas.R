comprobacion_dim <- function(imagen_orig, imagen_comprimida){
  ancho_orig <- width(imagen_orig)
  alto_orig <- height(imagen_orig)
  
  ancho_comp <- width(imagen_comp)
  alto_comp <- height(imagen_comp)
  
  if ((ancho_orig != ancho_comp) || (alto_orig != alto_comp)){
    imagen_comprimida <- resize(imagen_comprimida, size_x = ancho_orig,
                                size_y = alto_orig)
  }
  
  return(imagen_comprimida)
}

MSE <- function(imagen_orig, imagen_comprimida){
  
  comprobacion_dim(imagen_orig, imagen_comprimida)
  
  diferencia <- imagen_orig - imagen_comprimida
  
  canales <- imsplit(diferencia, "c")
  
  mse_canales <- c()
  
  for (canal in canales){
    matriz <- as.matrix(canal)
    mse_canales <- c(mse_canales, c(1/(ancho*alto)*sum(matriz^2)))
  }
  
  mse <- sum(mse_canales)/3
  print(paste('MSE:', mse))
  
  return(mse)
}

# Al elevar al cuadrado, penalizamos los errores grandes (píxeles muy diferentes
# en la imagen comprimida respecto de la original) y perdonamos los pequeños
# errores de redondeo.

MAE <- function(imagen_orig, imagen_comprimida){
  
  ancho_orig <- width(imagen_orig)
  alto_orig <- height(imagen_orig)
  
  ancho_comp <- width(imagen_comp)
  alto_comp <- height(imagen_comp)
  
  if ((ancho_orig != ancho_comp) || (alto_orig != alto_comp)){
    imagen_comprimida <- resize(imagen_comprimida, size_x = ancho_orig,
                                size_y = alto_orig)
  }
  
  diferencia <- imagen_orig - imagen_comprimida
  
  canales <- imsplit(diferencia, "c")
  
  mse_canales <- c()
  
  for (canal in canales){
    matriz <- as.matrix(canal)
    mse_canales <- c(mse_canales, c(1/(ancho*alto)*sum(matriz^2)))
  }
  
  mse <- sum(mse_canales)/3
  print(paste('MSE:', mse))
  
  return(mse)
}


MAPE <- function(imagen, imagen_comprimida){
  
  ancho_orig <- width(imagen_orig)
  alto_orig <- height(imagen_orig)
  
  ancho_comp <- width(imagen_comp)
  alto_comp <- height(imagen_comp)
  
  
  return(mape)
}