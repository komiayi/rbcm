############################################
varcovar.estimes <- function(cor_coefs, cor_cste=3){
  
  if(cor_cste==3){
    Delta1 <- NULL; Delt <- NULL; Delta3 <- NULL
    
    Delta1 <- -(prod(cor_coefs[c(1,3,5)])^2+prod(cor_coefs[c(2,4,6)])^2)+
      ((cor_coefs[1])^2-(cor_coefs[2])^2)*((cor_coefs[3])^2-(cor_coefs[4])^2)
    Delt <- (Delta1)^2 - 4*prod(cor_coefs)^2
    
    Delta3 <- 0.5*(-Delta1 + c(sqrt(ifelse(Delt>0,Delt,0)), -sqrt(ifelse(Delt>0,Delt,0))))/(prod(cor_coefs[c(2,3)])^2)
    Delta3 <- Delta3[which(Delta3>=0)]
    r01 <- c(sqrt(Delta3),-sqrt(Delta3))
    indices <- which((r01^2 - cor_coefs[6]^2)*(r01^2 - cor_coefs[5]^2) <= 0)
    if(length(indices) > 0) {
      r01 <- r01[indices]
    }
    r10 <- prod(cor_coefs[c(5,6)])/r01
  }else{
    r01 <- cor_coefs[3]
    r10 <- r01
    indices <- 1
  }
  return(list(r01 = r01, r10 = r10, indices=indices))
}

###################################################################################
varcovarEstimes <- function(cor_coefs, cor_cste = 3){
  
  if(cor_cste == 3){
    Delta1 <- NULL
    Delt <- NULL
    Delta3 <- NULL
    
    # Calcul de Delta1
    Delta1 <- -(prod(cor_coefs[c(1, 3, 5)])^2 + prod(cor_coefs[c(2, 4, 6)])^2) +
      ((cor_coefs[1])^2 - (cor_coefs[2])^2) * ((cor_coefs[3])^2 - (cor_coefs[4])^2)
    
    # Calcul de Delt
    Delt <- (Delta1)^2 - 4 * prod(cor_coefs)^2
    
    # Vérification si Delt est négatif
    if(Delt < 0 || is.na(Delt)){
      return(list(r01 = NA, r10 = NA, indices = NA))  # Retourner directement si Delt < 0
    } else {
      Delta3 <- 0.5 * (-Delta1 + c(sqrt(Delt), -sqrt(Delt))) / (prod(cor_coefs[c(2, 3)])^2)
      # Remplacer les éléments de Delta3 < 0 par NA
      Delta3[Delta3 < 0] <- NA
      
      # Si Delta3 contient des valeurs valides, calculer r01
      
      r01 <- c(sqrt(Delta3), -sqrt(Delta3))  # Calcul de r01 avec le premier élément de Delta3 valide
      
      
      # Vérification de la validité de r01
      indices <- which((r01^2 - cor_coefs[6]^2) * (r01^2 - cor_coefs[5]^2) <= 0)
      # if(length(indices) > 0) {
      #   r01 <- r01[indices]
      # }
      
      # Calcul de r10
      r10 <- prod(cor_coefs[c(5, 6)]) / r01
    }
    
  } else {
    # Cas où cor_cste n'est pas égal à 3
    r01 <- cor_coefs[3]
    r10 <- r01
    indices <- 1
  }
  
  return(list(r01 = r01, r10 = r10, indices = indices))
}
