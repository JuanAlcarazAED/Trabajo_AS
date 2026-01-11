############################################################
# SPIHT COMPLETO OPTIMIZADO (CODIFICACIÓN + DECODIFICACIÓN)
# CON BÚSQUEDA ITERATIVA (SIN RECURSIÓN) DE DESCENDIENTES
############################################################

############################################################
# HIJOS Y PRECOMPUTACIÓN
############################################################

Hijos <- function(i, j, size) {
  coords <- rbind(
    c(2*i,   2*j),
    c(2*i,   2*j+1),
    c(2*i+1, 2*j),
    c(2*i+1, 2*j+1)
  )
  valid <- coords[,1] < size & coords[,2] < size
  
  if (!any(valid)) {
    return(matrix(numeric(0), ncol = 2))
  }
  
  coords[valid, , drop = FALSE]
}

PrecomputeChildren <- function(size) {
  hijos <- vector("list", size*size)
  
  for (i in 0:(size-1)) {
    for (j in 0:(size-1)) {
      idx <- i*size + j + 1
      hijos[[idx]] <- Hijos(i, j, size)
    }
  }
  
  hijos
}

############################################################
# VERSIÓN ITERATIVA: ¿TIENE DESCENDIENTES SIGNIFICATIVOS?
############################################################

TieneDescendientesSignificativos <- function(i, j, size, W, T, hijos) {
  queue <- hijos[[i*size + j + 1]]
  
  if (is.null(queue) || length(queue) == 0) return(FALSE)
  
  while (!is.null(queue) && nrow(queue) > 0) {
    if (any(abs(W[queue[,1]+1, queue[,2]+1]) >= T, na.rm = TRUE)) {
      return(TRUE)
    }
    
    new_queue <- NULL
    for (k in 1:nrow(queue)) {
      idx <- queue[k,1] * size + queue[k,2] + 1
      H <- hijos[[idx]]
      if (!is.null(H) && length(H) > 0) {
        new_queue <- rbind(new_queue, H)
      }
    }
    
    queue <- new_queue
  }
  
  FALSE
}

############################################################
# INICIALIZACIÓN CODIFICADOR
############################################################

InicializarSPIHT <- function(W) {
  size <- nrow(W)
  half <- size / 2
  
  raices <- expand.grid(0:(half-1), 0:(half-1))
  raices <- as.matrix(raices)
  
  LIP <- raices
  LSP <- matrix(numeric(0), ncol = 2)
  
  LIS <- lapply(1:nrow(raices), function(k) {
    list(type = "A", coord = raices[k, ])
  })
  
  n <- floor(log2(max(abs(W))))
  
  # Intervalos [L,H) para coeficientes significativos
  L <- matrix(0, nrow = size, ncol = size)
  H <- matrix(0, nrow = size, ncol = size)
  
  list(LIP = LIP, LIS = LIS, LSP = LSP, n = n, L = L, H = H)
}


############################################################
# PROCESAR LIP (CODIFICACIÓN)
############################################################

ProcesarLIP <- function(W, estado, T) {
  LIP <- estado$LIP
  LSP <- estado$LSP
  L   <- estado$L
  H   <- estado$H
  
  bits <- c()
  eliminar <- c()
  
  if (nrow(LIP) > 0) {
    for (k in 1:nrow(LIP)) {
      i <- LIP[k,1]
      j <- LIP[k,2]
      
      if (abs(W[i+1,j+1]) >= T) {
        bits <- c(bits, 1)
        signo <- ifelse(W[i+1,j+1] < 0, 1, 0)
        bits <- c(bits, signo)
        
        # El coeficiente se vuelve significativo
        LSP <- rbind(LSP, c(i,j))
        
        # Intervalo inicial [T, 2T)
        L[i+1, j+1] <- T
        H[i+1, j+1] <- 2*T
        
        eliminar <- c(eliminar, k)
      } else {
        bits <- c(bits, 0)
      }
    }
  }
  
  if (length(eliminar) > 0)
    LIP <- LIP[-eliminar, , drop = FALSE]
  
  estado$LIP <- LIP
  estado$LSP <- LSP
  estado$L   <- L
  estado$H   <- H
  
  list(estado = estado, bits = bits)
}


############################################################
# PROCESAR LIS (CODIFICACIÓN, ITERATIVO)
############################################################

ProcesarLIS <- function(W, estado, T, hijos) {
  LIS  <- estado$LIS
  LIP  <- estado$LIP
  LSP  <- estado$LSP
  L    <- estado$L
  H    <- estado$H
  size <- nrow(W)
  
  bits <- c()
  nuevos_LIS <- list()
  
  if (length(LIS) > 0) {
    for (elem in LIS) {
      i <- elem$coord[1]
      j <- elem$coord[2]
      idx <- i*size + j + 1
      
      if (elem$type == "A") {
        significativo <- TieneDescendientesSignificativos(i, j, size, W, T, hijos)
        
        if (significativo) {
          bits <- c(bits, 1)
          
          Hhijos <- hijos[[idx]]
          if (!is.null(Hhijos) && length(Hhijos) > 0) {
            for (h in 1:nrow(Hhijos)) {
              hi <- Hhijos[h,1]
              hj <- Hhijos[h,2]
              
              if (abs(W[hi+1,hj+1]) >= T) {
                bits <- c(bits, 1)
                signo <- ifelse(W[hi+1,hj+1] < 0, 1, 0)
                bits <- c(bits, signo)
                LSP <- rbind(LSP, c(hi,hj))
                
                # Intervalo inicial [T, 2T) para este coeficiente
                L[hi+1, hj+1] <- T
                H[hi+1, hj+1] <- 2*T
              } else {
                bits <- c(bits, 0)
                LIP <- rbind(LIP, c(hi,hj))
              }
            }
          }
          
          # Sigue teniendo descendientes por debajo de los hijos
          if (TieneDescendientesSignificativos(i, j, size, W, T, hijos)) {
            nuevos_LIS[[length(nuevos_LIS)+1]] <- list(type="B", coord=c(i,j))
          }
          
        } else {
          bits <- c(bits, 0)
          nuevos_LIS[[length(nuevos_LIS)+1]] <- elem
        }
        
      } else if (elem$type == "B") {
        significativo <- TieneDescendientesSignificativos(i, j, size, W, T, hijos)
        
        if (significativo) {
          bits <- c(bits, 1)
          Hhijos <- hijos[[idx]]
          if (!is.null(Hhijos) && length(Hhijos) > 0) {
            for (h in 1:nrow(Hhijos)) {
              nuevos_LIS[[length(nuevos_LIS)+1]] <- list(type="A", coord=Hhijos[h,])
            }
          }
        } else {
          bits <- c(bits, 0)
          nuevos_LIS[[length(nuevos_LIS)+1]] <- elem
        }
      }
    }
  }
  
  estado$LIP <- LIP
  estado$LSP <- LSP
  estado$L   <- L
  estado$H   <- H
  estado$LIS <- nuevos_LIS
  
  list(estado = estado, bits = bits)
}


############################################################
# PASO DOMINANTE (CODIFICACIÓN)
############################################################

PasoDominante <- function(W, estado, hijos) {
  T <- 2^(estado$n)
  bits <- c()
  
  out1 <- ProcesarLIP(W, estado, T)
  estado <- out1$estado
  bits <- c(bits, out1$bits)
  
  out2 <- ProcesarLIS(W, estado, T, hijos)
  estado <- out2$estado
  bits <- c(bits, out2$bits)
  
  list(estado = estado, bits = bits)
}

############################################################
# PASO SUBORDINADO (CODIFICACIÓN)
############################################################

PasoSubordinado <- function(W, estado) {
  bits <- c()
  
  LSP <- estado$LSP
  L   <- estado$L
  H   <- estado$H
  
  if (nrow(LSP) > 0) {
    for (k in 1:nrow(LSP)) {
      i <- LSP[k,1]
      j <- LSP[k,2]
      
      val_real <- abs(W[i+1,j+1])
      L_ij <- L[i+1, j+1]
      H_ij <- H[i+1, j+1]
      
      if (L_ij == 0 && H_ij == 0) {
        next
      }
      
      mid <- (L_ij + H_ij) / 2
      
      if (val_real >= mid) {
        bits <- c(bits, 1)
        # Intervalo superior [mid, H)
        L[i+1, j+1] <- mid
      } else {
        bits <- c(bits, 0)
        # Intervalo inferior [L, mid)
        H[i+1, j+1] <- mid
      }
    }
  }
  
  estado$L <- L
  estado$H <- H
  
  list(estado = estado, bits = bits)
}


############################################################
# CODIFICADOR SPIHT COMPLETO
############################################################

SPIHT_Codificar <- function(W, max_bits = NULL, max_iter = NULL) {
  size <- nrow(W)
  hijos <- PrecomputeChildren(size)
  
  estado <- InicializarSPIHT(W)
  bitstream <- c()
  iter <- 0
  
  cat("n_inicial:", estado$n, "\n")
  
  repeat {
    iter <- iter + 1
    
    out_dom <- PasoDominante(W, estado, hijos)
    estado <- out_dom$estado
    bitstream <- c(bitstream, out_dom$bits)
    
    out_sub <- PasoSubordinado(W, estado)
    estado <- out_sub$estado
    bitstream <- c(bitstream, out_sub$bits)
    
    estado$n <- estado$n - 1
    
    cat("Iter:", iter, " n:", estado$n, " bits:", length(bitstream), "\n")
    
    if (!is.null(max_bits) && length(bitstream) >= max_bits) { 
      cat("Parada por max_bits\n"); break 
    }
    if (!is.null(max_iter) && iter >= max_iter) { 
      cat("Parada por max_iter\n"); break 
    }
    if (estado$n < 0) { 
      cat("Parada por n<0\n"); break 
    }
  }
  
  bitstream
}


############################################################
# INICIALIZACIÓN DECODIFICADOR
############################################################

InicializarSPIHT_Decode <- function(size, n_inicial) {
  half <- size / 2
  
  raices <- expand.grid(0:(half-1), 0:(half-1))
  raices <- as.matrix(raices)
  
  LIP <- raices
  LSP <- matrix(numeric(0), ncol = 2)
  
  LIS <- lapply(1:nrow(raices), function(k) {
    list(type="A", coord=raices[k,])
  })
  
  L <- matrix(0, nrow = size, ncol = size)
  H <- matrix(0, nrow = size, ncol = size)
  
  list(LIP=LIP, LIS=LIS, LSP=LSP, n=n_inicial, L=L, H=H)
}


############################################################
# DECODIFICAR LIS (USA HIJOS PRECOMPUTADOS)
############################################################
DecodificarLIP <- function(W, estado, bitstream, pos, T) {
  LIP <- estado$LIP
  LSP <- estado$LSP
  L   <- estado$L
  H   <- estado$H
  
  if (nrow(LIP) > 0) {
    for (k in 1:nrow(LIP)) {
      b <- bitstream[pos]; pos <- pos + 1
      
      if (b == 1) {
        signo <- bitstream[pos]; pos <- pos + 1
        sgn <- ifelse(signo == 1, -1, 1)
        
        # Intervalo inicial [T, 2T)
        L[LIP[k,1]+1, LIP[k,2]+1] <- T
        H[LIP[k,1]+1, LIP[k,2]+1] <- 2*T
        
        # Aproximación: punto medio del intervalo
        W[LIP[k,1]+1, LIP[k,2]+1] <- sgn * ( (T + 2*T) / 2 )
        
        LSP <- rbind(LSP, LIP[k,])
        LIP[k,] <- NA
      }
    }
  }
  
  LIP <- LIP[complete.cases(LIP), , drop = FALSE]
  
  estado$LIP <- LIP
  estado$LSP <- LSP
  estado$L   <- L
  estado$H   <- H
  
  list(W=W, estado=estado, pos=pos)
}

DecodificarLIS <- function(W, estado, bitstream, pos, T, hijos) {
  LIS <- estado$LIS
  LIP <- estado$LIP
  LSP <- estado$LSP
  L   <- estado$L
  H   <- estado$H
  size <- nrow(W)
  
  nuevos_LIS <- list()
  
  if (length(LIS) > 0) {
    for (elem in LIS) {
      i <- elem$coord[1]
      j <- elem$coord[2]
      idx <- i*size + j + 1
      
      b <- bitstream[pos]; pos <- pos + 1
      
      if (elem$type == "A") {
        if (b == 1) {
          Hhijos <- hijos[[idx]]
          
          if (!is.null(Hhijos) && length(Hhijos) > 0) {
            for (h in 1:nrow(Hhijos)) {
              hi <- Hhijos[h,1]
              hj <- Hhijos[h,2]
              
              b2 <- bitstream[pos]; pos <- pos + 1
              
              if (b2 == 1) {
                signo <- bitstream[pos]; pos <- pos + 1
                sgn <- ifelse(signo == 1, -1, 1)
                
                # Intervalo inicial [T, 2T)
                L[hi+1, hj+1] <- T
                H[hi+1, hj+1] <- 2*T
                
                # Aproximación: punto medio
                W[hi+1, hj+1] <- sgn * ( (T + 2*T) / 2 )
                
                LSP <- rbind(LSP, c(hi,hj))
              } else {
                LIP <- rbind(LIP, c(hi,hj))
              }
            }
          }
          
          nuevos_LIS[[length(nuevos_LIS)+1]] <- list(type="B", coord=c(i,j))
          
        } else {
          nuevos_LIS[[length(nuevos_LIS)+1]] <- elem
        }
        
      } else if (elem$type == "B") {
        if (b == 1) {
          Hhijos <- hijos[[idx]]
          if (!is.null(Hhijos) && length(Hhijos) > 0) {
            for (h in 1:nrow(Hhijos)) {
              nuevos_LIS[[length(nuevos_LIS)+1]] <- list(type="A", coord=Hhijos[h,])
            }
          }
        } else {
          nuevos_LIS[[length(nuevos_LIS)+1]] <- elem
        }
      }
    }
  }
  
  estado$LIP <- LIP
  estado$LSP <- LSP
  estado$LIS <- nuevos_LIS
  estado$L   <- L
  estado$H   <- H
  
  list(W=W, estado=estado, pos=pos)
}

############################################################
# PASO DOMINANTE (DECODIFICACIÓN)
############################################################

PasoDominante_Decode <- function(W, estado, bitstream, pos, hijos) {
  T <- 2^(estado$n)
  
  out1 <- DecodificarLIP(W, estado, bitstream, pos, T)
  W <- out1$W
  estado <- out1$estado
  pos <- out1$pos
  
  out2 <- DecodificarLIS(W, estado, bitstream, pos, T, hijos)
  W <- out2$W
  estado <- out2$estado
  pos <- out2$pos
  
  list(W=W, estado=estado, pos=pos)
}

############################################################
# PASO SUBORDINADO (DECODIFICACIÓN)
############################################################

PasoSubordinado_Decode <- function(W, estado, bitstream, pos) {
  LSP <- estado$LSP
  L   <- estado$L
  H   <- estado$H
  
  if (nrow(LSP) > 0) {
    for (k in 1:nrow(LSP)) {
      i <- LSP[k,1]
      j <- LSP[k,2]
      
      b <- bitstream[pos]; pos <- pos + 1
      
      L_ij <- L[i+1, j+1]
      H_ij <- H[i+1, j+1]
      
      if (L_ij == 0 && H_ij == 0) {
        next
      }
      
      mid <- (L_ij + H_ij) / 2
      
      if (b == 1) {
        # Intervalo superior [mid, H)
        L_ij <- mid
      } else {
        # Intervalo inferior [L, mid)
        H_ij <- mid
      }
      
      L[i+1, j+1] <- L_ij
      H[i+1, j+1] <- H_ij
      
      # Nueva aproximación: punto medio del nuevo intervalo
      sgn <- ifelse(W[i+1,j+1] < 0, -1, 1)
      W[i+1,j+1] <- sgn * ( (L_ij + H_ij) / 2 )
    }
  }
  
  estado$L <- L
  estado$H <- H
  
  list(W=W, estado=estado, pos=pos)
}



############################################################
# DECODIFICADOR SPIHT COMPLETO
############################################################

SPIHT_Decodificar <- function(bitstream, size, n_inicial) {
  
  hijos <- PrecomputeChildren(size)
  
  W <- matrix(0, size, size)
  estado <- InicializarSPIHT_Decode(size, n_inicial)
  
  pos <- 1
  
  repeat {
    out_dom <- PasoDominante_Decode(W, estado, bitstream, pos, hijos)
    W <- out_dom$W
    estado <- out_dom$estado
    pos <- out_dom$pos
    
    out_sub <- PasoSubordinado_Decode(W, estado, bitstream, pos)
    W <- out_sub$W
    estado <- out_sub$estado
    pos <- out_sub$pos
    
    estado$n <- estado$n - 1
    
    if (estado$n < 0 || pos > length(bitstream)) break
  }
  
  W
}

