library(mvtnorm)
library(Matrix)
##################################################################################
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

###################################################################
Applications_simple_effect <- function(par, treat, outcome, mediators, intmed, out_cov=NULL,
                                       intmed_cov=NULL, sed_cov=NULL, inter=FALSE,
                                       inter_treat_cov = TRUE, cor_cste=1,data=NULL,
                                       B=200,rh=0.5, methode = c("delta", "bootstrap"), formula_one3 = NULL){
  
  DE1s <- IE1s <- NULL
  coefs.intmed.all <- coefs.sedmed.all <- coefs.y.all <-   cor_coefs.all <-  NULL
  
  rho01.all <- rho10.all <- indi.all <- VarbootDE <- VarbootIE <- NULL
  
  DEic <- IEic <- NULL
  
  ################################################
  
  
  if(isTRUE(inter)){
    names_vec <- c( 
      paste(treat, mediators, sep=":"), 
      paste(mediators, collapse = ":"), 
      paste(treat, paste(mediators, collapse = ":"), sep=":"))
  }else{
    names_vec <- NULL
  }
  
  if(is.array(data)){
    J <- dim(data)[3]
  }else if(is.list(data)){
    J <- length(data)
  }else{
    J <- 1
  }
  
  q1 <- 2 + length(intmed_cov)
  q2 <- 2 + length(sed_cov)
  
  
  if(isTRUE(inter)){
    q3 <- 2 + 3*length(mediators)+ length(out_cov)
  }else{
    q3 <- 2 + length(mediators) + length(out_cov)
  }
  
  for(i in 1:J){
    
    donne.i <- NULL
    
    if(is.array(data)){
      donne.i <- as.data.frame(data[,,i])
    }else if(is.list(data)){
      donne.i <- as.data.frame(data[[i]])
    }else{
      donne.i <- as.data.frame(data)
    }
    
    
    Efdi <- calculer_effets(treat = treat, mediators= mediators,intmed = intmed, outcome = outcome,
                            intmed_cov = intmed_cov, sed_cov=sed_cov, out_cov=out_cov,
                            inter_treat_cov = inter_treat_cov, cor_cste=cor_cste, inter = inter,
                            q1 = q1, q2 = q2, names_vec = names_vec, rh =rh, data= donne.i,
                            formula_one3 = formula_one3)
    
    if(cor_cste==3){
      if(!is.null(DE1s) && ncol(DE1s) != 6){
        next
      }
    }
    
    
    ## méthode bootstrap 
    
    if(methode == "bootstrap"){
      
      Efdiboot <- estim.boots(B=B,treat = treat, mediators= mediators,intmed = intmed, outcome = outcome,
                              intmed_cov = intmed_cov, sed_cov=sed_cov, out_cov=out_cov,
                              inter_treat_cov = inter_treat_cov, cor_cste=cor_cste, inter = inter,
                              q1 = q1, q2 = q2, names_vec = names_vec, rh =rh, donne= donne.i, tripl=tripl)
      
      # deIC <- c(as.vector(Efdi$sol.i$DE1sk)-1.96*sqrt(Efdiboot$Vde), as.vector(Efdi$sol.i$DE1sk)+1.96*sqrt(Efdiboot$Vde))
      # ieIC <- c(as.vector(Efdi$sol.i$IE1sk)-1.96*sqrt(Efdiboot$Vie), as.vector(Efdi$sol.i$IE1sk)+1.96*sqrt(Efdiboot$Vie))
      
      VarbootDE <- rbind(VarbootDE, Efdiboot$Vde)
      VarbootIE <- rbind(VarbootIE, Efdiboot$Vie)
      DEic <- rbind(DEic, Efdiboot$deIC)
      IEic <- rbind(IEic, Efdiboot$ieIC)
    }else{
      VarbootDE <- VarbootIE <- DEic <- IEic <- NA
    }
    ################### Effets ##########################
    
    DE1s <- rbind(DE1s, as.vector(Efdi$sol.i$DE1sk))
    IE1s <- rbind(IE1s, as.vector(Efdi$sol.i$IE1sk))
    
    ##################################
    
    coefs.intmed.all <- rbind(coefs.intmed.all, Efdi$Cint)
    coefs.sedmed.all <- rbind(coefs.sedmed.all, Efdi$Csed)
    coefs.y.all <-rbind(coefs.y.all, Efdi$C.y)
    cor_coefs.all <- rbind(cor_coefs.all, Efdi$Ccor)
    
    rho01.all <- rbind(rho01.all,Efdi$var_covar$r01)
    rho10.all <- rbind(rho10.all,Efdi$var_covar$r10)
    indi.all[[i]] <-   Efdi$var_covar$indices
    
    ####################################################
    
    #print(i)
  }
  
  ###########################
  if(!is.null(DE1s) && ncol(DE1s)>2){
    snbsol <- ncol(DE1s)-2
  }else{
    snbsol <- 1
  }
  
  
  if(cor_cste == 3){
    VnamesMed <- c("interM1", paste(treat,1, sep=""), paste(intmed_cov,1, sep=""), 
                   "interM2", paste(treat,2, sep=""), paste(sed_cov,2, sep=""),"sig11", "sig10",
                   "sig21", "sig20", "rh11","rh00" )
    low <- paste0("Low_", 1:(snbsol+2))
    upo <- paste0("Up_", 1:(snbsol+2))
    r01names<- paste0("rh01_",1:snbsol)
    r10names<- paste0("rh10_",1:snbsol)
    
    denames <- paste0("De_", c(paste0("racine",1:snbsol),'moyen','fixe'))
    ienames <- paste0("Ie_", c(paste0("racine",1:snbsol),'moyen','fixe'))
    
    vardenames <- paste0("VarDe_", c(paste0("racine",1:snbsol),'moyen','fixe'))
    varienames <- paste0("VarIe_", c(paste0("racine",1:snbsol),'moyen','fixe'))
    
  }else{
    VnamesMed <- c("interM1", paste(treat,1, sep=""), paste(intmed_cov,1, sep=""),
                   "interM2", paste(treat,2, sep=""), paste(sed_cov,2, sep=""),"sig1", "sig2", "rh1")
    low <- paste0("Low_", 1)
    upo <- paste0("Up_", 1)
    r01names<- "rh01"
    r10names<- "rh10"
    denames <- "De" 
    ienames <- "Ie"
    vardenames <- "VarDe"
    varienames <- "VarIe"
  }
  
  VnamesY <- c("interY", paste(treat,3, sep=""), mediators, paste(out_cov,3, sep=""), names_vec, "sd") 
  
  ##############################
  
  
  colnames(coefs.intmed.all) <- VnamesMed[1:q1]
  colnames(coefs.sedmed.all) <- VnamesMed[(q1+1):(q1+q2)]
  colnames(cor_coefs.all) <- VnamesMed[-(1:(q1+q2))]
  colnames(coefs.y.all) <- VnamesY
  colnames(rho01.all) <- r01names
  colnames(rho10.all) <- r10names
  
  colnames(DE1s) <- denames
  colnames(IE1s) <- ienames
  
  if(methode == "bootstrap"){
    colnames(DEic) <- colnames(IEic) <- c(low, upo)
    colnames(VarbootDE) <- vardenames
    colnames(VarbootIE) <- varienames
  }
  
  
  return(list(IndEff = IE1s, DirEff = DE1s, coefm1 = coefs.intmed.all, 
              coefm2 = coefs.sedmed.all, coefy = coefs.y.all, corcoef = cor_coefs.all,
              rho01 = rho01.all, rho10 = rho10.all, Indir01 = indi.all, 
              VarbootDE = VarbootDE, VarbootIE= VarbootIE, DEIC = DEic, IEIC = IEic))
}

#############################################


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




##################### outcome loglik ############################################

Matr_Design <- function(inter_vec, data = NULL){
  # Calcul des interactions basées sur names_vec
  for (term in inter_vec) {
    # Séparer les variables de l'interaction
    variables <- strsplit(term, ":")[[1]]
    
    # Calculer l'interaction et ajouter à data
    interaction_result <- 1
    for (var in variables) {
      interaction_result <- interaction_result * data[[var]]
    }
    
    # Ajouter la nouvelle colonne d'interaction
    data[[term]] <- interaction_result
  }
  
  der <- data[inter_vec]
  colnames(der) <- inter_vec
  # Afficher les données avec les nouvelles colonnes d'interaction
  return(der)
}



##################### estiamtions des paramètres médiateurs ####################




##################################################################################


calculer_effets <- function(treat, mediators, intmed, outcome, intmed_cov, sed_cov, 
                            out_cov, inter_treat_cov, cor_cste, data= NULL, q1, 
                            q2, inter, names_vec, rh, formula_one3){
  
  # Initialisation des paramètres
  parinit <- initialParams(treat = treat, mediators = mediators, intmed = intmed, outcome = outcome,
                           intmed_cov = intmed_cov, sed_cov = sed_cov, out_cov = out_cov,
                           inter_treat_cov = inter_treat_cov, cor_cste = cor_cste, data = data,
                           formula_one3 = formula_one3, names_vec = names_vec)
  
  # # # Estimation des paramètres des médiateurs
  
  coefs.intmed <- as.vector(parinit$med[1:q1])
  coefs.sedmed <- as.vector(parinit$med[(q1+1):(q1+q2)])
  cor_coefs <- as.vector(parinit$med[-c(1:(q1+q2))])
  
  # Création de la matrice des coefficients
  coef <- data.frame(matrix(c(coefs.intmed, coefs.sedmed), 2, q1, byrow = TRUE), mediators)
  colnames(coef) <- c("inter", treat, intmed_cov, "name")
  
  # Estimation des paramètres pour Y
   
  coefs.y <- parinit$outc
  
  # Création du tableau des coefficients de Y
  if (isTRUE(inter)) {
    Bet <- data.frame(Beta = coefs.y, name = c("inter", treat, mediators, out_cov, names_vec, "sd"))
  } else {
    Bet <- data.frame(Beta = coefs.y, name = c("inter", treat, mediators, out_cov, "sd"))
  }
  
  # Calcul des covariances et corrélations
  var_covar <- varcovarEstimes(cor_coefs, cor_cste = cor_cste)
  
  if (cor_cste == 1) {
    r011 <- as.vector(var_covar$r01)
    r100 <- as.vector(var_covar$r10)
  } else {
    r011 <- c(as.vector(var_covar$r01), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
    r100 <- c(as.vector(var_covar$r10), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
  }
  
  # Calcul des effets naturels
  sol.i <- effectdirectindirct(alpha = coef, beta = Bet, treat = treat, mediators = mediators,
                               intmed_cov = intmed_cov, sed_cov = sed_cov, inter = inter, ro = r011,
                               corC = cor_coefs, names_vec = names_vec, cor_cste = cor_cste, data = data)
  
  return(list(Cint = coefs.intmed, Csed = coefs.sedmed, Ccor = cor_coefs, C.y = coefs.y,
              var_covar = var_covar, 
              sol.i = sol.i))
}


################################### calules varainces bootstrap #################

calculer_effets_bootstrap <- function(treat, mediators, intmed, outcome, intmed_cov, sed_cov, 
                                      out_cov, inter_treat_cov, cor_cste, data= NULL,
                                      q1, q2, inter, names_vec, rh, tripl=TRUE){
  
  # Initialisation des paramètres
  parinit <- initialParams(treat = treat, mediators = mediators, intmed = intmed, outcome = outcome,
                           intmed_cov = intmed_cov, sed_cov = sed_cov, out_cov = out_cov,
                           inter_treat_cov = inter_treat_cov, cor_cste = cor_cste,
                           data = data, tripl = tripl)
  
  coefs.intmed <- as.vector(parinit$med[1:q1])
  coefs.sedmed <- as.vector(parinit$med[(q1+1):(q1+q2)])
  cor_coefs <- as.vector(parinit$med[-c(1:(q1+q2))])
  
  # Création de la matrice des coefficients
  coef <- data.frame(matrix(c(coefs.intmed, coefs.sedmed), 2, q1, byrow = TRUE), mediators)
  colnames(coef) <- c("inter", treat, intmed_cov, "name")
  
  # Estimation des paramètres pour Y
  
  coefs.y <- parinit$outc
  
  # Création du tableau des coefficients de Y
  if (isTRUE(inter)) {
    Bet <- data.frame(Beta = coefs.y, name = c("inter", treat, mediators, out_cov, names_vec, "sd"))
  } else {
    Bet <- data.frame(Beta = coefs.y, name = c("inter", treat, mediators, out_cov, "sd"))
  }
  
  # Calcul des covariances et corrélations
  var_covar <- varcovarEstimes(cor_coefs, cor_cste = cor_cste)
  
  if (cor_cste == 1) {
    r011 <- as.vector(var_covar$r01)
    r100 <- as.vector(var_covar$r10)
  } else {
    r011 <- c(as.vector(var_covar$r01), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
    r100 <- c(as.vector(var_covar$r10), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
  }
  
  # Calcul des effets naturels
  sol.i <- effectdirectindirct(alpha = coef, beta = Bet, treat = treat, mediators = mediators,
                               intmed_cov = intmed_cov, sed_cov = sed_cov, inter = inter, ro = r011,
                               corC = cor_coefs, names_vec = names_vec, cor_cste = cor_cste, data = data)
  
  return(list(sol = sol.i, Vcorr = var_covar$r01))
}


################################################################################
estim.boots <- function(B, treat, mediators, intmed, outcome, intmed_cov, sed_cov, 
                        out_cov, inter_treat_cov, cor_cste, donne = NULL, q1, q2, 
                        inter, names_vec, rh, tripl=TRUE){
  
  DE1s.b <- NULL
  IE1s.b <- NULL
  B_effectif <- 0
  
  
  while (B_effectif < B){
    # Tirage bootstrap
    donne.b <- donne[sample(1:nrow(donne), replace = TRUE),]
    
    # Vérifier que treat == 0 et treat == 1 sont présents avec assez d'observations
    n1 <- sum(donne.b[, treat] == 1)
    n0 <- sum(donne.b[, treat] == 0)
    
    if (n1 > length(intmed_cov) + 1 && n0 > length(sed_cov) + 1 && length(unique(donne.b[, treat])) > 1) {
      # Calcul des effets naturels
      sole <- calculer_effets_bootstrap(treat = treat, mediators = mediators, intmed = intmed, outcome = outcome,
                                        intmed_cov = intmed_cov, sed_cov = sed_cov, out_cov = out_cov,
                                        inter_treat_cov = inter_treat_cov, cor_cste = cor_cste, inter = inter,
                                        q1 = q1, q2 = q2, names_vec = names_vec, rh = rh,
                                        data = donne.b, tripl=tripl)
      
      solb <- sole$sol
      varc <- sole$Vcorr
      
      # Vérifier que solb$DE1sk et solb$IE1sk ne contiennent pas de NA et que le nombre de colonnes correspond
      if (all(!is.na(solb$DE1sk)) && all(!is.na(solb$IE1sk)) && 
          (is.null(DE1s.b) || ncol(DE1s.b) == length(solb$IE1sk))) { 
        DE1s.b <- rbind(DE1s.b, matrix(solb$DE1sk, nrow = 1))
        IE1s.b <- rbind(IE1s.b, matrix(solb$IE1sk, nrow = 1))
        
        # Vérifier si nous avons atteint B lignes valides
        if (nrow(DE1s.b) >= B) {
          break
        }
        
        # Incrémenter B_effectif en fonction du nombre de lignes valides
        B_effectif <- nrow(DE1s.b)
      }
      print(paste("B", B_effectif))
    }
  }
  
  # Calcul des variances en prenant en compte les lignes valides uniquement
  Vde <- ((B-1)/B)*apply(DE1s.b, 2, var, na.rm = TRUE)
  Vie <- ((B-1)/B)*apply(IE1s.b, 2, var, na.rm = TRUE)
  
  ci_lower_DE <- apply(DE1s.b, 2, quantile, probs = 0.025)
  ci_upper_DE <- apply(DE1s.b, 2, quantile, probs = 0.975)
  ci_lower_IE <- apply(IE1s.b, 2, quantile, probs = 0.025)
  ci_upper_IE <- apply(IE1s.b, 2, quantile, probs = 0.975)
  
  return(list(Vde = Vde, Vie = Vie, deIC= c( ci_lower_DE, ci_upper_DE),
              ieIC = c(ci_lower_IE,ci_upper_IE)))
}



############################### calcul des effets ##############################

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


