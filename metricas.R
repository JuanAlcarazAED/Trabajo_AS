

MSE <- function(imagen_orig, imagen_comprimida){
  
  ancho_orig <- width(imagen_orig)
  alto_orig <- height(imagen_orig)
  
  ancho_comp <- width(imagen_comp)
  alto_comp <- height(imagen_comp)
  
  if ((ancho_orig != ancho_comp) || (alto_orig != alto_comp)){
    imagen_orig <- resize(imagen_orig, size_x = ancho_comp, size_y = alto_comp)
  }
  
  mse <- 1/(ancho*alto)
  return(mse)
}

# Al elevar al cuadrado, penalizamos los errores grandes (píxeles muy diferentes
# en la imagen comprimida respecto de la original) y perdonamos los pequeños
# errores de redondeo.

MAPE <- function(imagen, imagen_comprimida){
  
  ancho_orig <- width(imagen_orig)
  alto_orig <- height(imagen_orig)
  
  ancho_comp <- width(imagen_comp)
  alto_comp <- height(imagen_comp)
  
  
  return(mape)
}