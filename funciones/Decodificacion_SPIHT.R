# Ahora las funciones para decodificar
##############################################
# FUNCIONES AUXILIARES (MISMAS QUE CODIFICADOR)
##############################################

Hijos <- function(i, j, size) {
  coords <- rbind(
    c(2*i,   2*j),
    c(2*i,   2*j+1),
    c(2*i+1, 2*j),
    c(2*i+1, 2*j+1)
  )
  valid <- coords[,1] < size & coords[,2] < size
  coords[valid, , drop = FALSE]
}

Descendientes <- function(i, j, size) {
  queue <- Hijos(i, j, size)
  if (is.null(queue) || nrow(queue) == 0)
    return(matrix(numeric(0), ncol = 2))
  
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

##############################################
# INICIALIZACIÓN DEL DECODIFICADOR
##############################################

InicializarSPIHT_Decode <- function(size, n_inicial) {
  half <- size / 2
  
  raices <- expand.grid(0:(half-1), 0:(half-1))
  raices <- as.matrix(raices)
  
  LIP <- raices
  LSP <- matrix(numeric(0), ncol = 2)
  
  LIS <- lapply(1:nrow(raices), function(k) {
    list(type = "A", coord = raices[k, ])
  })
  
  list(LIP = LIP, LIS = LIS, LSP = LSP, n = n_inicial)
}

##############################################
# PASO DOMINANTE (DECODIFICACIÓN)
##############################################

DecodificarLIP <- function(W, estado, bitstream, pos, T) {
  LIP <- estado$LIP
  LSP <- estado$LSP
  
  if (nrow(LIP) > 0) {
    for (k in 1:nrow(LIP)) {
      b <- bitstream[pos]; pos <- pos + 1
      
      if (b == 1) {
        signo <- bitstream[pos]; pos <- pos + 1
        valor <- ifelse(signo == 1, -T, T)
        W[LIP[k,1]+1, LIP[k,2]+1] <- valor
        
        LSP <- rbind(LSP, LIP[k,])
        LIP[k,] <- NA
      }
    }
  }
  
  LIP <- LIP[complete.cases(LIP), , drop = FALSE]
  
  estado$LIP <- LIP
  estado$LSP <- LSP
  
  list(W = W, estado = estado, pos = pos)
}

DecodificarLIS <- function(W, estado, bitstream, pos, T) {
  LIS  <- estado$LIS
  LIP  <- estado$LIP
  LSP  <- estado$LSP
  size <- nrow(W)
  
  nuevos_LIS <- list()
  
  for (elem in LIS) {
    i <- elem$coord[1]
    j <- elem$coord[2]
    
    b <- bitstream[pos]; pos <- pos + 1
    
    if (elem$type == "A") {
      D <- Descendientes(i, j, size)
      
      if (b == 1) {
        H <- Hijos(i, j, size)
        
        for (h in 1:nrow(H)) {
          b2 <- bitstream[pos]; pos <- pos + 1
          
          if (b2 == 1) {
            signo <- bitstream[pos]; pos <- pos + 1
            valor <- ifelse(signo == 1, -T, T)
            W[H[h,1]+1, H[h,2]+1] <- valor
            LSP <- rbind(LSP, H[h,])
          } else {
            LIP <- rbind(LIP, H[h,])
          }
        }
        
        if (nrow(D) > nrow(H)) {
          nuevos_LIS[[length(nuevos_LIS)+1]] <- list(type="B", coord=c(i,j))
        }
        
      } else {
        nuevos_LIS[[length(nuevos_LIS)+1]] <- elem
      }
      
    } else if (elem$type == "B") {
      D <- Descendientes(i, j, size)
      
      if (b == 1) {
        H <- Hijos(i, j, size)
        for (h in 1:nrow(H)) {
          nuevos_LIS[[length(nuevos_LIS)+1]] <- list(type="A", coord=H[h,])
        }
      } else {
        nuevos_LIS[[length(nuevos_LIS)+1]] <- elem
      }
    }
  }
  
  estado$LIP <- LIP
  estado$LSP <- LSP
  estado$LIS <- nuevos_LIS
  
  list(W = W, estado = estado, pos = pos)
}

PasoDominante_Decode <- function(W, estado, bitstream, pos) {
  T <- 2^(estado$n)
  
  out1 <- DecodificarLIP(W, estado, bitstream, pos, T)
  W     <- out1$W
  estado <- out1$estado
  pos    <- out1$pos
  
  out2 <- DecodificarLIS(W, estado, bitstream, pos, T)
  W     <- out2$W
  estado <- out2$estado
  pos    <- out2$pos
  
  list(W = W, estado = estado, pos = pos)
}

##############################################
# PASO SUBORDINADO (DECODIFICACIÓN)
##############################################

PasoSubordinado_Decode <- function(W, estado, bitstream, pos) {
  T  <- 2^(estado$n)
  M  <- T + T/2
  
  if (nrow(estado$LSP) > 0) {
    for (k in 1:nrow(estado$LSP)) {
      i <- estado$LSP[k,1]
      j <- estado$LSP[k,2]
      
      b <- bitstream[pos]; pos <- pos + 1
      
      if (b == 1) {
        if (W[i+1, j+1] >= 0) W[i+1, j+1] <- M
        else                 W[i+1, j+1] <- -M
      } else {
        if (W[i+1, j+1] >= 0) W[i+1, j+1] <- T
        else                 W[i+1, j+1] <- -T
      }
    }
  }
  
  list(W = W, estado = estado, pos = pos)
}

##############################################
# FUNCIÓN PRINCIPAL SPIHT (DECODIFICACIÓN)
##############################################

SPIHT_Decodificar <- function(bitstream, size, n_inicial) {
  W <- matrix(0, size, size)
  estado <- InicializarSPIHT_Decode(size, n_inicial)
  
  pos <- 1
  
  repeat {
    out_dom <- PasoDominante_Decode(W, estado, bitstream, pos)
    W       <- out_dom$W
    estado  <- out_dom$estado
    pos     <- out_dom$pos
    
    out_sub <- PasoSubordinado_Decode(W, estado, bitstream, pos)
    W       <- out_sub$W
    estado  <- out_sub$estado
    pos     <- out_sub$pos
    
    estado$n <- estado$n - 1
    
    if (estado$n < 0 || pos > length(bitstream)) break
  }
  
  W
}
