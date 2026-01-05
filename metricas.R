comprobacion_dim <- function(imagen_orig, imagen_comprimida){
  
  ancho_orig <- width(imagen_orig)
  alto_orig <- height(imagen_orig)
  
  ancho_comp <- width(imagen_comprimida)
  alto_comp <- height(imagen_comprimida)
  
  if ((ancho_orig != ancho_comp) || (alto_orig != alto_comp)){
    imagen_comprimida <- resize(imagen_comprimida, size_x = ancho_orig,
                                size_y = alto_orig)
  }
  
  return(imagen_comprimida)
}

## Métricas Basadas en el Error de Píxeles

# Mean Squared Error
MSE <- function(imagen_orig, imagen_comprimida){
  
  imagen_comprimida <- comprobacion_dim(imagen_orig, imagen_comprimida)
  
  ancho <- width(imagen_orig)
  alto <- height(imagen_orig)
  
  diferencia <- imagen_orig - imagen_comprimida
  
  canales <- imsplit(diferencia, "c")
  
  mse_canales <- c()
  
  for (canal in canales){
    matriz <- as.matrix(canal)
    mse_canales <- c(mse_canales, c(1/(ancho*alto)*sum(matriz^2)))
  }
  
  mse <- sum(mse_canales)/3
  # print(paste('MSE:', mse))
  
  return(mse)
}

# Al elevar al cuadrado, penalizamos los errores grandes (píxeles muy diferentes
# en la imagen comprimida respecto de la original) y perdonamos los pequeños
# errores de redondeo. Esta métrica no refleja bien la calidad percibida por el
# ojo humano.

# Mean Absolute Error
MAE <- function(imagen_orig, imagen_comprimida){
  
  imagen_comprimida <- comprobacion_dim(imagen_orig, imagen_comprimida)
  
  ancho <- width(imagen_orig)
  alto <- height(imagen_orig)
  
  diferencia <- imagen_orig - imagen_comprimida
  
  canales <- imsplit(diferencia, "c")
  
  mae_canales <- c()
  
  for (canal in canales){
    matriz <- as.matrix(canal)
    mae_canales <- c(mae_canales, c(1/(ancho*alto)*sum(abs(matriz))))
  }
  
  mae <- sum(mae_canales)/3
  # print(paste('MAE:', mae))
  
  return(mae)
}

# Root Mean Squared Error
RMSE <- function(imagen_orig, imagen_comprimida){
  
  mse <- MSE(imagen_orig, imagen_comprimida)

  rsme <- sqrt(mse)
  
  return(rsme)
}

# Esta métrica es útil para interpretar el error en términos de magnitud de píxeles.

# Peak Signal-to-Noise Ratio
PSNR <- function(imagen_orig, imagen_comprimida){
  
  imagen_comprimida <- comprobacion_dim(imagen_orig, imagen_comprimida)
  
  mse <- MSE(imagen_orig, imagen_comprimida) 
  
  psnr <- 10*log10(1^2/mse)
  
  return(psnr)
}

# Se mide en decibelios. A mayor PSNR mayor calidad. Tiene la misma limitación que
# el MSE: no refleja bien la calidad percibida por el ojo humano.

## Métricas Perceptuales

# Structural Similarity Index
SSMI <- function(imagen_orig, imagen_comprimida){
  return(ssmi)
}


## Métricas Específicas de Compresión

# Razón de Compresión (Compression Ratio)
CR <- function(ruta_imagen_orig, ruta_imagen_comprimida){
  
  cr <- file.size(ruta_imagen_orig)/file.size(ruta_imagen_comprimida)
  
  return(cr)
}


