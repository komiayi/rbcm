
TableCoefficients <- function(J, parmed, pary){
  Mmodel.m <- rmvnorm(J, mean = c(parmed$Cint,parmed$Csed,parmed$Ccor), sigma = solve(parmed$Chess))
  Ymodel.y <- rmvnorm(J, mean = pary$C.y, sigma = solve(pary$Chess))
  return(list(alph = Mmodel.m, beta = Ymodel.y))
}


corCoe <- function(imp_m, treatment,mediators, covariates, outcome, data = NULL){
  residualsData <- residualsData1 <- residualsData0 <-coeffsData <- NULL
  formulesList <- list()  # Stocker les formules utilisées
  
  # Régression de Y sur les médiateurs, le traitement et les covariables
  
  fitY <- lm(as.formula(paste(outcome,paste(c(treatment,mediators,covariates), collapse = "+"),sep="~")), data = data)
  coeffDataY <- summary(fitY)$coefficients[c(treatment,imp_m),c('Estimate','Pr(>|t|)')]
  colnames(coeffDataY) <- c('B', "P-value")
  
  for(i in seq_along(imp_m)){
    formules <- as.formula(paste(mediators[i],paste(c(treatment,covariates), collapse = "+"),sep="~"))
    formulestext <-  deparse(formules)  # Stocker la formule dans la liste
    
    if(all(data[,treatment] %in% c(0,1), na.rm = TRUE)){
      # Ajustement pour x == 1 et x == 0
      
      res1 <- lm(formules, data = filter(data, x == 1))$residuals
      res0 <- lm(formules, data = filter(data, x == 0))$residuals
      
      residualsData1 <- cbind(residualsData1, res1)
      residualsData0 <- cbind(residualsData0, res0)
    }
    # Ajustement général sur l’ensemble des données
    
    res.fit <- lm(formules, data = data)
    res.coef <- summary(res.fit)$coefficients[treatment,c('Estimate','Pr(>|t|)')]
    
    
    res.coef <- data.frame(
      Formula = formulestext, B_treat = res.coef["Estimate"],
      P_value_treat = res.coef["Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
    coeffsData <- rbind(coeffsData, res.coef)
    residualsData <- cbind(residualsData, res.fit$residuals)
  }
  
  #colnames(coeffsData) <- c('B', "P-value")
  rownames(coeffsData) <- NULL
  
  # Calcul des corrélations
  if(all(data[,treatment] %in% c(0,1), na.rm = TRUE)){
    colnames(residualsData1) <- colnames(residualsData0) <- imp_m
    CorelationData1 <- cor(residualsData1, method = 'pearson')
    CorelationData0 <- cor(residualsData0, method = 'pearson')
  }
  
  colnames(residualsData) <- imp_m
  CorelationData <- cor(residualsData, method = 'pearson')
  
  
  # Retourner les résultats avec les formules utilisées
  
  if(all(data[,treatment] %in% c(0,1), na.rm = TRUE)){
    return(list(coeffs = coeffDataY, 
                coeffsMed = coeffsData, 
                CorD = round(CorelationData,3),
                CorD1 = round(CorelationData1,3),
                CorD0 = round(CorelationData0,3)))
  }else{
    return(list(coeffs = coeffDataY, 
                coeffsMed = coeffsData, 
                CorD = round(CorelationData,3)))
  }
}
