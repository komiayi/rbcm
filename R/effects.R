effectdirectindirct <- function(alpha, beta, treat, mediators, intmed_cov, sed_cov,
                                inter, vcovar, corC, names_vec, cor_cste, ro, data = NULL){
  
  alp01 = alpha[alpha$name%in%mediators[1], "inter"]
  alp02 = alpha[alpha$name%in%mediators[-1], "inter"]
  alp11 = alpha[alpha$name%in%mediators[1], treat]
  alp12 = alpha[alpha$name%in%mediators[-1], treat]
  alp21 = unlist(alpha[alpha$name%in%mediators[1], intmed_cov])
  alp22 = unlist(alpha[alpha$name%in%mediators[-1], sed_cov])
  
  beta0 <- beta[1,"Beta"]
  beta1 <- beta[beta$name%in%treat,"Beta"]
  beta21 <- beta[beta$name%in%mediators[1],"Beta"]
  beta22 <- beta[beta$name%in%mediators[-1],"Beta"]
  
  if(isTRUE(inter)){
    Bt = beta[beta$name%in%names_vec, "Beta"]
    for (j in seq_along(names_vec)) {
      assign(paste0("beta3", j), Bt[j])
    }
  }else{
    for (j in 1:4) {
      assign(paste0("beta3", j), 0)
    }
  }
  
  DEk <- IEk <- NULL
  
  if(cor_cste==3){
    
    for(k in seq_along(ro)){
      
      termde1 <- beta1
      termde2 <- alp12*beta22
      if(is.null(intmed_cov)){
        termde31 <- alp01
      }else{
        termde31 <- alp01 + as.matrix(data[,intmed_cov])%*%alp21
      }
      termde3 <- termde31*beta31
      if(is.null(sed_cov)){
        termde41 <- alp02 + alp12
      }else{
        termde41 <- alp02 + alp12 + as.matrix(data[,sed_cov])%*%alp22
      }
      termde4 <- termde41*beta32
      
      termde5 <- (alp12*termde31 + prod(corC[c(2,3)])*ro[k] - prod(corC[c(2,4,6)]))*beta33
      termde6 <-(termde31*termde41 + prod(corC[c(2,3)])*ro[k])*beta34
      
      de <- mean(termde1 + termde2 + termde3 + termde4 + termde5 + termde6)
      
      termie1 <- alp11*(beta21 + beta31)
      termie2 <- (alp11*termde41 + prod(corC[c(1,3,5)]) - prod(corC[c(2,3)])*ro[k])*(beta33 + beta34)
      ie <- mean(termie1 + termie2)
      
      DEk <- c(DEk, de)
      IEk <- c(IEk, ie)
    }
    
  }else{
    
    for(k in seq_along(ro)){
      
      termde1 <- beta1
      termde2 <- alp12*beta22
      if(is.null(intmed_cov)){
        termde31 <- alp01
      }else{
        termde31 <- alp01 + as.matrix(data[,intmed_cov])%*%alp21
      }
      termde3 <- termde31*beta31
      if(is.null(sed_cov)){
        termde41 <- alp02 + alp12
      }else{
        termde41 <- alp02 + alp12 + as.matrix(data[,sed_cov])%*%alp22
      }
      termde4 <- termde41*beta32
      
      termde5 <- (alp12*termde31)*beta33
      termde6 <-(termde31*termde41 + prod(corC[c(1,2)])*ro[k])*beta34
      
      de <- mean(termde1 + termde2+ termde3+termde4 +termde5+termde6)
      
      termie1 <- alp11*(beta21 + beta31)
      termie2 <- alp11*termde41*(beta33 + beta34)
      ie <- mean(termie1 + termie2)
      
      DEk <- c(DEk, de)
      IEk <- c(IEk, ie)
    }
  }
  
  return(list(DE1sk = DEk, IE1sk = IEk))
}


############################################
initialParams <- function(treat, mediators,intmed,outcome,intmed_cov, sed_cov, out_cov, 
                          inter_treat_cov = TRUE, cor_cste=1, data=NULL,
                          formula_one3=NULL, names_vec = NULL){
  # initialisation
  coeffs.model.m1 <- NULL
  fit.model.m1 <- NULL; formula_two1 <- NULL; formula_one1 <- NULL
  coeffs.model.m2 <- NULL
  fit.model.m2 <- NULL; formula_two <- NULL; formula_one <- NULL
  formula_two2 <- NULL; formula_one2 <- NULL
  formula_two3 <- NULL
  ################################################################
  
  if(is.null(intmed_cov)){
    formula_one1 <- treat
  }else if(isTRUE(inter_treat_cov)){
    formula_one1 <- paste(treat,intmed_cov, sep = "*")
  }else{
    formula_one1 <- paste(c(treat,intmed_cov), collapse = "+")
  }
  # regressions
  formula_two1 <- as.formula(paste(intmed, formula_one1, sep="~"))
  fit.model.m1 <- lm(formula_two1, data = data)
  
  ### M2
  
  if(is.null(sed_cov)){
    formula_one2 <- treat
  }else if(isTRUE(inter_treat_cov)){
    formula_one2 <- paste(treat,sed_cov, sep = "*")
  }else{
    formula_one2 <- paste(c(treat,sed_cov), collapse = "+")
  }
  # regression
  formula_two2 <- as.formula(paste(mediators[-1], formula_one2, sep="~"))
  fit.model.m2 <- lm(formula_two2, data = data)
  ############ Y
  # regressions
  formula_two3 <- as.formula(paste(outcome, formula_one3, sep="~"))
  fit.outcome <- lm(formula_two3, data = data)
  #########################
  fitcoef <- fit.outcome$coefficients
  # Initialiser un vecteur avec des zéros de la même longueur que names_vec
  result <- setNames(rep(0, length(names_vec)), names_vec)
  
  # Remplacer les éléments existants par leurs valeurs dans fitcoef
  result[names_vec %in% names(fitcoef)] <- fitcoef[names_vec[names_vec %in% names(fitcoef)]]
  
  fitcoef <- as.vector(c(fitcoef[1:(2+length(mediators)+length(out_cov))], result))
  
  #fitcoef[is.na(fitcoef)]
  
  if(cor_cste==3){
    n1 <- nrow(data[data[,treat]==1,])
    n2 <- nrow(data[data[,treat]==0,])
    a11 <- sd(lm(formula_two1, data = data[data[,treat]==1,])$residuals)*sqrt((n1-1)/(n1-1-length(intmed_cov)))
    a10 <- sd(lm(formula_two1, data = data[data[,treat]==0,])$residuals)*sqrt((n2-1)/(n2-1-length(intmed_cov)))
    a21 <- sd(lm(formula_two2, data = data[data[,treat]==1,])$residuals)*sqrt((n1-1)/(n1-1-length(sed_cov)))
    a20 <- sd(lm(formula_two2, data = data[data[,treat]==0,])$residuals)*sqrt((n2-1)/(n2-1-length(sed_cov)))
    
    # corrélations
    
    r11 <- cor(lm(formula_two1, data = data[data[,treat]==1,])$residuals,
               lm(formula_two2, data = data[data[,treat]==1,])$residuals)
    r00 <- cor(lm(formula_two1, data = data[data[,treat]==0,])$residuals,
               lm(formula_two2, data = data[data[,treat]==0,])$residuals)
    return(list(med = c(as.vector(fit.model.m1$coefficients),as.vector(fit.model.m2$coefficients), a11,a10,a21,a20,r11,r00),
                outc = c(fitcoef, sd(fit.outcome$residuals))))
  }else{
    a11 <- sd(lm(formula_two1, data = data)$residuals)*sqrt((nrow(data)-1)/(nrow(data)-2-length(sed_cov)))
    a10 <- a11
    a21 <- sd(lm(formula_two2, data = data)$residuals)*sqrt((nrow(data)-1)/(nrow(data)-2-length(sed_cov)))
    a20 <- a21
    
    # corrélations
    
    r11 <- cor(lm(formula_two1, data = data)$residuals,
               lm(formula_two2, data = data)$residuals)
    r00 <- r11
    return(list(med = c(as.vector(fit.model.m1$coefficients),as.vector(fit.model.m2$coefficients), a11,a21,r11),
                outc = c(fitcoef,sd(fit.outcome$residuals))))
  }
  
}