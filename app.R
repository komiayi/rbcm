# Import Libraries
# Install and load the shinydashboard package, which simultaneously loads shiny.
#install.packages("shinydashboard")
library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tseries)  # Pour le test Jarque-Bera
library(readr)
library(mvtnorm)
library(Matrix)
library(MVN)
#########################


tableau_html <- function(data, data2) {
  effets <- c("\\(\\zeta\\)", "\\(\\delta\\)")
  
  headers <- c(
    "Effets", "Méthode", "Corrélation estimée",
    "Valeur estimée", "Écart-type<sup>1</sup>", "95% IC<sup>1</sup>"
  )
  
  headers_html <- sapply(headers, function(h) {
    words <- unlist(strsplit(h, " "))
    if (length(words) > 1) {
      paste(words[1], "<br>", paste(words[-1], collapse = " "), sep = "")
    } else {
      h
    }
  })
  
  html <- '
  <style>
    table.custom-table {
      margin: auto;
      border-collapse: collapse;
      font-family: Arial, sans-serif;
      font-size: 14px;
      width: 80%;
      border: none;
    }
    .custom-table th,
    .custom-table td {
      text-align: center;
      vertical-align: middle;
      padding: 8px 12px;
      border: none;
    }
    .custom-table thead th {
      border-bottom: 2px solid #999999;
      background-color: #f2f2f2;
    }
    .custom-table tbody tr:nth-child(even) td {
      background-color: #f9f9f9;
    }
    /* Enlever les lignes horizontales sauf celle avant les résultats de data2 */
    .custom-table tbody tr:not(:first-child):not(:nth-child(3)) td {
      border-bottom: none;
    }
    .custom-table tbody tr:nth-child(3) td {
      border-bottom: 1px solid #dddddd; /* Ligne horizontale avant les résultats de data2 */
    }
    .custom-table tbody tr:first-child td {
      border-bottom: none;
    }
    .custom-table tbody tr:last-child td {
      border-bottom: none;
    }
  </style>

  <table class="custom-table">
    <thead>
      <tr>'
  
  for (h in headers_html) {
    html <- paste0(html, '<th>', h, '</th>')
  }
  
  html <- paste0(html, '</tr>
    </thead>
    <tbody>')
  
  # Première méthode : CC
  for (i in seq_along(effets)) {
    html <- paste0(html, '<tr><td>', effets[i], '</td>')
    if (i == 1) {
      html <- paste0(html,
                     '<td rowspan="2" style="vertical-align:middle;">\\(\\texttt{CC}\\)</td>')
    }
    html <- paste0(html, '<td>', "", '</td>')
    html <- paste0(html, '<td>', data[i, "Valeur"], '</td>')
    html <- paste0(html, '<td>', data[i, "ET"], '</td>')
    html <- paste0(html, '<td>', data[i, "IC"], '</td></tr>')
  }
  
  # Méthodes CNCm et CNCr
  effets2 <- c("\\(\\zeta\\)", "\\(\\delta\\)")
  
  for (effet in effets2) {
    sous_data <- data2[data2$Effet == effet, ]
    n <- nrow(sous_data)
    
    for (i in seq_len(n)) {
      ligne <- sous_data[i, ]
      
      if (i == 1) {  # Première ligne => CNCm
        html <- paste0(html, "<tr><td>", effet, "</td>")
        html <- paste0(html, '<td >\\(\\texttt{CNCm}\\)</td>')
      } else {  # Lignes suivantes => CNCr
        html <- paste0(html, "<tr><td></td>")
        html <- paste0(html, '<td rowspan="1" style="vertical-align:middle;">\\(\\texttt{CNCr}\\)</td>')
      }
      html <- paste0(html,
                     "<td>", ligne$rho, "</td>",
                     "<td>", ligne$Valeur, "</td>",
                     "<td>", ligne$ET, "</td>",
                     "<td>", ligne$IC, "</td></tr>")
    }
  }
  
  html <- paste0(html, '</tbody></table>')
  # Ajouter une ligne fine avant "Bootstrap"
  html <- paste0(html, '<hr style="border: 0; border-top: 1px solid #ddd;">')
  
  # Le texte "Bootstrap" à gauche et une ligne après
  html <- paste0(html, '<p style="text-align:left; margin-left: 20px;"><sup>1</sup> Bootstrap</p>')
  html <- paste0(html, '<hr style="border: 0; border-top: 1px solid #ddd;">')
  
  # Ajout de MathJax pour rendre le LaTeX
  html <- paste0(html, '
  <script type="text/javascript" async
    src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
  </script>
  <script type="text/javascript">
    MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
  </script>')
  
  return(HTML(html))
}




##############################################################
# Définir la fonction de test de normalité
test_normalite <- function(vos_donnees) {
  # Sélection du test en fonction de la taille de l'échantillon
  n <- length(vos_donnees)
  if (n <= 50) {
    # Test de Shapiro-Wilk pour les petits échantillons
    test_shapiro <- shapiro.test(vos_donnees)
    return(paste("Test de Shapiro-Wilk:\n", capture.output(print(test_shapiro), type = "output"), collapse = "\n"))
  } else if (n <= 2000) {
    # Test de Kolmogorov-Smirnov pour les échantillons moyens
    test_ks <- ks.test(vos_donnees, "pnorm", mean = mean(vos_donnees), sd = sd(vos_donnees))
    return(print(test_ks))
  } else {
    # Test de Jarque-Bera pour les grands échantillons
    test_jarque <- jarque.bera.test(vos_donnees)
    return(paste("Test de Jarque-Bera:\n", capture.output(print(test_jarque), type = "output"), collapse = "\n"))
  }
}

######################################################################

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
######################################################################




######################

# Set Up UI Components
header <- dashboardHeader(title = "Rbcm", # Customizable header with title 
                          titleWidth = 300)  # Adjusted with to make the title fit
sidebar <- dashboardSidebar(# Sidebar for navigation
  width = 300,  # Adjust sidebar width
  sidebarMenu(
    menuItem("Main Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Méthodes", tabName = "methods", icon = icon("chart-line")),
    
    # Analyse Descriptive avec sous-menus
    menuItem("Analyse Descriptive", icon = icon("th"),
             menuSubItem("Vue d'ensemble", tabName = "overview"),
             menuSubItem("Résumé statistique", tabName = "stat_summary"),
             menuSubItem("Graphiques", tabName = "graphs")
    ),
    
    menuItem("Mediation", tabName = "mediate", icon = icon("chart-bar")),
    menuItem("External Link", href = "https://tilburgsciencehub.com", icon = icon("external-link")),
    
    # Search Form: For enhanced user interaction
    sidebarSearchForm(textId = "search", buttonId = "searchButton", label = "Search...")
  )
)

body <- dashboardBody(

  tabItems(
    tabItem(tabName = "methods",
            h2("Méthodes"),
            withMathJax(),
            # Inclure le fichier Markdown contenant le texte RMarkdown et les formules
            uiOutput("markdown_ui")
    ),
    tabItem(tabName = "overview",
            h2("Descriptive analysis tab content"),
            
            fluidRow(
              column(6, fileInput("data_file", "Télécharger un fichier de données (CSV)", 
                                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".RData"))),
              column(6, textInput("treat", "Exposition")),
              column(6, textInput("outcome", "Réponse")),
              column(6, textInput("mediators", "Médiateurs (séparés par des virgules)")),
              column(6, textInput("intmed", "Médiateur secondaire"))
            ),
            
            fluidRow(
              column(6, textInput("out_cov", "Covariables de la réponse (séparées par des virgules)")),
              column(6, textInput("intmed_cov", "Covariables pour le médiateur primaire (séparées par des virgules)")),
              column(6, textInput("sed_cov", "Covariables pour le médiateur secondaire (séparées par des virgules)")),
              column(6, numericInput("B", "Nombre de rééchantillonnages (Bootstrap)", value = 20)),
              column(6, radioButtons("interaction", "Inclure une interaction entre les médiateurs ?", 
                                     choices = list("TRUE" = TRUE, "FALSE" = FALSE), selected = FALSE)),  # Boolean TRUE/FALSE
              # Section conditionnelle
              conditionalPanel(
                condition = "input.interaction == 'TRUE'",
                column(12,
                       checkboxGroupInput("interaction_types", "Choisir les interactions à inclure :",
                                          choices = c(
                                            "Exposition*Médiateur primaire" = "treat_m1",
                                            "Exposition*Médiateur secondaire" = "treat_m2",
                                            "Médiateur primaire*Médiateur secondaire" = "m1_m2",
                                            "Exposition*Médiateur primaire*Médiateur secondaire" = "treat_m1_m2"
                                          ))
                )
              )
            ),
            
            fluidRow(
              column(12, actionButton("run", "Exécuter", class = "btn-primary")),
              verbatimTextOutput("debug_output")  # Pour afficher les messages de débogage
            )
    ),
    tabItem(tabName = "stat_summary",
            h2("Résumé statistique"),
            fluidRow(
              box(title = "Aperçu des données", width = 12, dataTableOutput("data_table")),
              box(title = "Aperçu des residus", width = 12, dataTableOutput("residual_table")),
              box(title = "Statistiques descriptives des variables quantitatives", width = 12, 
                  dataTableOutput("summary_stats_quant")),
              box(title = "Statistiques descriptives des variables qualitatives", width = 12,
                  dataTableOutput("summary_stats_qual"))
            )
    ),
    tabItem(tabName = "graphs",
            h2("Graphiques et tests statistiques"),
            
            fluidRow(
              tabBox(
                title = "Boxplot",
                side = "right", height = "450px",
                selected = "Réponse",
                tabPanel("Médiateur secondaire", plotOutput("boxplotmedsed")),
                tabPanel("Médiateur primaire", plotOutput("boxplotmedint")),
                tabPanel("Réponse", plotOutput("boxplotoutcome"))
                
              ),
              tabBox(
                title = "Histogramme",
                side = "right", height = "450px",
                selected = "Réponse",
                tabPanel("Médiateur secondaire", plotOutput("histogrammedsed")),
                tabPanel("Médiateur primaire", plotOutput("histogrammedint")),
                tabPanel("Réponse", plotOutput("histogramoutcome"))
                
              )
            ),
            
            fluidRow(
              tabBox(
                title = "QQplot",
                side = "right", height = "450px",
                selected = "Réponse",
                tabPanel("Médiateur secondaire", plotOutput("qqplotmedsed")),
                tabPanel("Médiateur primaire", plotOutput("qqplotmedint")),
                tabPanel("Réponse", plotOutput("qqplotoutcome"))
              ),
              tabBox(
                title = "Test de normalité",
                side = "right", height = "350px",
                selected = "Réponse",
                tabPanel("Multivariée", verbatimTextOutput("testmulti")),
                tabPanel("Univarié (Médiateur secondaire)", verbatimTextOutput("testmedsed")),
                tabPanel("Univarié (Médiateur primaire)", verbatimTextOutput("testmedint")),
                tabPanel("Réponse", verbatimTextOutput("testoutcome"))
              )
            )
    ),
    
    tabItem(tabName = "mediate",
            h2("Résultats Effets de médiation"),
            
            fluidRow(
              box(title = "", width = 12, 
                  uiOutput("resultsCC")),
               box(title = "Corrélation non constante entre les médiateurs", width = 12,
                   dataTableOutput("resultsCNC"))
            )
    )
  )
)  # Main body for content

# Assemble UI
ui <- dashboardPage(header, sidebar, body)  # Combine header, sidebar, and body



# Serveur
server <- function(input, output, session) {
  
  output$markdown_ui <- renderUI({
    includeHTML("Methodes.html")
    #tags$iframe(seamless="seamless", src = "Methodes.html", width = "100%", height = "600px")
  })
  
  # Variable pour stocker les données (persiste pendant la session)
  data <- reactiveVal(NULL)
  
  # Charger les données lorsque l'utilisateur clique sur "Exécuter"
  observeEvent(input$run, {
    req(input$data_file)  # Vérifier qu'un fichier a été téléchargé
    ext <- tools::file_ext(input$data_file$name)
    
    # Message de debug pour vérifier le format du fichier et la réponse
    output$debug_output <- renderPrint({
      paste("Fichier téléchargé :", input$data_file$name, "Extension :", ext)
    })
    
    # Lecture des données selon l'extension
    loaded_data <- NULL  # Initialiser la variable de données
    
    if (ext == "csv") {
      loaded_data <- tryCatch({
        read.csv(input$data_file$datapath)
      }, error = function(e) {
        output$debug_output <- renderPrint({ paste("Erreur lors du chargement du fichier CSV :", e$message) })
        NULL
      })
      
    } else if (ext == "RData") {
      tryCatch({
        load(input$data_file$datapath)  # Charger le fichier RData
        obj_names <- ls()[1]  # Récupérer les objets dans l'environnement
        if (length(obj_names) == 0) {
          stop("Aucun objet trouvé dans le fichier RData.")
        } else if (length(obj_names) == 1) {
          # Si un seul objet, le charger directement
          loaded_data <- get(obj_names[1])
        } else {
          # Si plusieurs objets, informer l'utilisateur
          output$debug_output <- renderPrint({
            paste("Le fichier RData contient plusieurs objets :", paste(obj_names, collapse = ", "), 
                  ". Veuillez spécifier lequel utiliser.")
          })
          return()  # Sortir pour éviter d'essayer de stocker des données non définies
        }
      }, error = function(e) {
        output$debug_output <- renderPrint({ paste("Erreur lors du chargement du fichier RData :", e$message) })
      })
    }
    
    # Conversion en data frame si les données sont une liste ou un tableau
    if (!is.null(loaded_data)) {
      CL <- colnames(loaded_data)
      if (is.array(loaded_data)) {
        loaded_data <- as.data.frame(loaded_data)  # Convertir un tableau en data frame
        colnames(loaded_data) <- CL
      } else if (is.list(loaded_data)) {
        loaded_data <- as.data.frame(loaded_data)  # Convertir une liste en data frame
        colnames(loaded_data) <- CL
      }else {
        loaded_data <- as.data.frame(loaded_data)  # data frame
        colnames(loaded_data) <- CL
      }
      
      data(loaded_data)  # Stocker les données
    }
    
    # Vérification des données chargées
    if (!is.null(data())) {
      output$debug_output <- renderPrint({
        paste("Données chargées avec succès, nombre de lignes :", nrow(data()))
      })
    } else {
      output$debug_output <- renderPrint({ "Aucune donnée chargée." })
    }
  })
  
  # Afficher un aperçu des données dans l'onglet "Résumé statistique"
  output$data_table <- renderDataTable({
    req(data())  # S'assurer que les données existent
    datatable(data(), options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))  # Afficher les premières lignes des données
  })
  
  # Calculer et stocker les résidus dans un data.frame
  get_residuals_df <- reactive({
    req(data())
    
    # Variables à ajuster
    outcome <- input$outcome
    exposure <- input$treat
    mediators <- if (is.character(input$mediators)) {
      strsplit(input$mediators, ",")[[1]] %>% trimws()
    } else {
      input$mediators
    }
    
    # les covariables de la réponse et des médiateurs
    out_cov_vars <- if (input$out_cov == "") {NULL}else{strsplit(input$out_cov, ",")[[1]] %>% trimws()} 
    intmed_cov_vars <- if (input$intmed_cov == "") {NULL}else{strsplit(input$intmed_cov, ",")[[1]] %>% trimws()}
    sed_cov_vars <- if (input$sed_cov == "") {NULL}else{strsplit(input$sed_cov, ",")[[1]] %>% trimws()}
    
    covariate <- if (is.null(intmed_cov_vars) && is.null(sed_cov_vars)) {
      NULL  # Si les deux variables sont NULL, retourne NULL
    } else {
      unique(c(intmed_cov_vars, sed_cov_vars))  # Sinon, combine les deux ensembles de variables et retire les doublons
    }
    
    inter <- as.logical(input$interaction)
    selected_interactions <- input$interaction_types  # vecteur de valeurs ex: c("treat_m1", "m1_m2")
    interaction_terms <- c()
    
    if (isTRUE(inter) && !is.null(selected_interactions)) {
      if ("treat_m1" %in% selected_interactions) {
        interaction_terms <- c(interaction_terms, paste(exposure, mediators[1], sep = ":"))
      }
      if ("treat_m2" %in% selected_interactions && length(mediators) > 1) {
        interaction_terms <- c(interaction_terms, paste(exposure, mediators[2], sep = ":"))
      }
      if ("m1_m2" %in% selected_interactions && length(mediators) > 1) {
        interaction_terms <- c(interaction_terms, paste(mediators[1], mediators[2], sep = ":"))
      }
      if ("treat_m1_m2" %in% selected_interactions && length(mediators) > 1) {
        triple <- paste(exposure, mediators[1], mediators[2], sep = ":")
        interaction_terms <- c(interaction_terms, triple)
      }
      lmfor <- paste(c(exposure, mediators, covariate, interaction_terms), collapse = " + ")
    }else{
      lmfor <- paste(c(exposure, mediators, covariate), collapse = " + ")
    }
    
    # Vérifications
    req(nzchar(outcome), nzchar(exposure), length(mediators) > 0)
    all_vars <- c(outcome, mediators)
    req(all(all_vars %in% names(data())), exposure %in% names(data()))
    resid_mat <- sapply(all_vars, function(var) {
      if (var == outcome) {
        formula_str <- paste(var, "~", lmfor)
      } else {
        predictors <- c(exposure, covariate)
        formula_str <- paste(var, "~", paste(predictors, collapse = " + "))
      }
      
      model <- tryCatch({
        lm(as.formula(formula_str), data = data())
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(model)){
        return(rep(NA, nrow(data())))
      }
      residuals(model)
    })
    #Récupérer les noms des modèles
    # model_names <- sapply(all_vars, function(var) {
    #   predictors <- if (var == outcome) c(exposure, mediators) else exposure
    #   paste(var, "~", paste(predictors, collapse = " + "))  # Formule du modèle
    # })
    # 
    # Mettre les résidus en data.frame
    resid_df <- as.data.frame(resid_mat)
    colnames(resid_df) <- all_vars  # nommer les colonnes avec les noms des variables d'origine
    resid_df
  })
  
  output$residual_table <- renderDataTable({
    df <- get_residuals_df()
    req(nrow(df) > 0)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
  })
  
  # Afficher les statistiques descriptives pour les variables quantitatives
  output$summary_stats_quant <- renderDataTable({
    req(data())

    # Extraire les variables d'intérêt
    treat_var <- input$treat
    outcome_var <- input$outcome
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws() # Méditeurs

    # les covariables de la réponse et des médiateurs
    out_cov_vars <- strsplit(input$out_cov, ",")[[1]]%>% trimws()   # Covariables de la réponse
    intmed_cov_vars <- strsplit(input$intmed_cov, ",")[[1]] %>% trimws()
    sed_cov_vars <- strsplit(input$sed_cov, ",")[[1]]%>% trimws()

    # Créer une liste complète des variables d'intérêt
    selected_vars <- unique(c(treat_var, outcome_var, mediator_vars, out_cov_vars,intmed_cov_vars,sed_cov_vars))
    selected_vars <- selected_vars[selected_vars != ""]  # Retirer les entrées vides

    # Calcul des statistiques descriptives avec quantiles et valeurs manquantes
    stats_dplyr <- data() %>%
      mutate(across(where(is.numeric), ~ if (length(unique(.)) < 4){as.factor(.)}else{.})) %>%
      summarise_if(is.numeric, list(
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        min = ~min(., na.rm = TRUE),
        Q1 = ~quantile(., 0.25, na.rm = TRUE),
        median = ~median(., na.rm = TRUE),
        Q3 = ~quantile(., 0.75, na.rm = TRUE),
        max = ~max(., na.rm = TRUE),
        missing = ~sum(is.na(.))
      ))

    # Vérifiez que stats_dplyr contient des données avant de tenter de transformer
    if (ncol(stats_dplyr) == 0) {
      return(NULL)  # Pas de données à afficher
    }

    # Réorganiser le tableau pour avoir les variables en lignes et les statistiques en colonnes
    stats_tidy <- stats_dplyr %>%
      pivot_longer(cols = everything(), names_to = "variable_statistic", values_to = "value") %>%
      separate(variable_statistic, into = c("variable", "statistic"), sep = "_(?=[^_]+$)", extra = "merge")%>%
      pivot_wider(names_from = statistic, values_from = value)%>%
      mutate(across(where(is.numeric), ~ round(.x, 5)))  # Arrondir les valeurs numériques à 5 décimales


    # # Filtrer pour ne garder que les variables numériques saisis
    stats_filtered <- stats_tidy
    if(length(selected_vars)!= 0){
      stats_filtered <- stats_tidy %>%
        mutate(variable = factor(variable, levels = c(selected_vars, unique(variable[!variable %in% selected_vars]))))%>%
        arrange(variable)
    }
    datatable(stats_filtered, options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
  })

  # Afficher les statistiques descriptives pour les variables qualitatives
  output$summary_stats_qual <- renderDataTable({
    req(data())

    # Extraire les variables d'intérêt
    treat_var <- input$treat
    outcome_var <- input$outcome
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws() # Méditeurs

    # les covariables de la réponse et des médiateurs
    out_cov_vars <- strsplit(input$out_cov, ",")[[1]]%>% trimws()   # Covariables de la réponse
    intmed_cov_vars <- strsplit(input$intmed_cov, ",")[[1]] %>% trimws()
    sed_cov_vars <- strsplit(input$sed_cov, ",")[[1]]%>% trimws()

    # Créer une liste complète des variables d'intérêt
    selected_vars <- unique(c(treat_var, outcome_var, mediator_vars, out_cov_vars,intmed_cov_vars,sed_cov_vars))
    selected_vars <- selected_vars[selected_vars != ""]  # Retirer les entrées vides

    # Calcul des statistiques descriptives avec quantiles et valeurs manquantes
    
    # Sélectionner uniquement les variables qualitatives (caractères, facteurs et facteurs ordonnés)
    
    qualitative_summary <- data() %>%
      mutate(across(where(is.numeric), ~ if (length(unique(.)) < 4){as.factor(.)}else{.})) %>%
      select(where(is.character), where(is.factor), where(is.ordered))
    
    
    
    # Vérifier s'il y a au moins une variable qualitative ou ordinale
    if (ncol(qualitative_summary) == 0) {
      return(NULL)  # Retourner NULL si aucune variable qualitative ou ordinale n'est présente
    }
    
    # Procéder à l'analyse des fréquences
    stats_filtered_quali <- qualitative_summary %>%
      drop_na() %>%  # Exclure les lignes avec des valeurs manquantes
      pivot_longer(cols = everything(), names_to = "Variable", values_to = "Valeur") %>%
      group_by(Variable, Valeur) %>%
      summarise(
        Frequence = n(),
        .groups = 'drop'
      )%>%
      group_by(Variable) %>%  # Re-grouper par variable
      mutate(Pourcentage = round((Frequence / sum(Frequence)) * 100,2)) # Calcul du pourcentage
    
    datatable(stats_filtered_quali, options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
  })

  ###########################################################
  # Boxplot
  output$boxplotoutcome <- renderPlot({
    # Extraire les résidus de l'outcome à partir du tableau de résidus
    residuals_df <- get_residuals_df()
    
    # Extraire les résidus spécifiques à l'outcome
    outcome_resid <- residuals_df[[input$outcome]]
    
    # Vérifier que les résidus existent pour l'outcome
    req(!is.null(outcome_resid), length(outcome_resid) > 0)
    
    # Créer le boxplot des résidus de l'outcome
    boxplot(outcome_resid, main = paste("Boxplot des résidus de", input$outcome),
            col = "lightblue", border = "black", xlab = input$outcome,
            width = 5, height = 4)
  })

  output$boxplotmedint <- renderPlot({
    residuals_df <- get_residuals_df()
    
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[1]]]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)

    boxplot(mediator_resid, main = paste("Boxplot des résidus du médiateur primaire"),
            col = "lightgreen", border = "black", xlab = mediator_vars[1],
            width = 5, height = 4)
  })

  output$boxplotmedsed <- renderPlot({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[2]]]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)
    
    # Créer le boxplot des résidus du médiateur secondaire
    boxplot(mediator_resid, main = paste("Boxplot des résidus du médiateur secondaire"),
            col = "lightcoral", border = "black", xlab = mediator_vars[2],
            width = 5, height = 4)
  })

  # Histogram
  output$histogramoutcome <- renderPlot({
    residuals_df <- get_residuals_df()
    
    # Extraire les résidus spécifiques à l'outcome
    outcome_resid <- residuals_df[[input$outcome]]
    
    # Vérifier que les résidus existent pour l'outcome
    req(!is.null(outcome_resid), length(outcome_resid) > 0)

    hist(outcome_resid , main = paste("Histogramme de", input$outcome),
         xlab = input$outcome, breaks = 20, col = "lightgray", freq=FALSE)
    lines(density(outcome_resid), col = "blue", lwd = 2)

  })

  output$histogrammedint <- renderPlot({
    residuals_df <- get_residuals_df()

    # Extraire les médiateurs
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[1]]]
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)

    hist(mediator_resid, main = paste("Histogramme des résidus du médiateur primaire"),
         xlab = mediator_vars[1], breaks = 20, col = "lightgray", freq=FALSE)
    lines(density(mediator_resid), col = "blue", lwd = 2)
  })

  output$histogrammedsed <- renderPlot({
    residuals_df <- get_residuals_df()

    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[2]]]
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)
    

    hist(mediator_resid, main = paste("Histogramme des résidus du médiateur secondiare"),
         xlab = mediator_vars[2], breaks = 20, col = "lightgray", freq=FALSE)
    lines(density(mediator_resid), col = "blue", lwd = 2)
  })

  # Quantile
  output$qqplotoutcome <- renderPlot({
    residuals_df <- get_residuals_df()
    
    # Extraire les résidus spécifiques à l'outcome
    outcome_resid <- residuals_df[[input$outcome]]
    
    # Vérifier que les résidus existent pour l'outcome
    req(!is.null(outcome_resid), length(outcome_resid) > 0)

    qqnorm(outcome_resid, main = "QQ plot")
    qqline(outcome_resid, col = "red", lwd = 2)

  })

  output$qqplotmedint <- renderPlot({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[1]]]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)
  
    qqnorm(mediator_resid, main = "QQ plot")
    qqline(mediator_resid, col = "red", lwd = 2)
  })

  output$qqplotmedsed <- renderPlot({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[2]]]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)
  
    qqnorm(mediator_resid, main = "QQ plot")
    qqline(mediator_resid, col = "red", lwd = 2)

  })

  ##### Test de normalité
  
  output$testoutcome <- renderPrint({
    # Extraire les résidus de l'outcome à partir du tableau de résidus
    residuals_df <- get_residuals_df()
    
    # Extraire les résidus spécifiques à l'outcome
    outcome_resid <- residuals_df[[input$outcome]]
    
    # Vérifier que les résidus existent pour l'outcome
    req(!is.null(outcome_resid), length(outcome_resid) > 0)
    
    test_normalite(outcome_resid)  # Appeler la fonction de test de normalité
  })

  output$testmedint <- renderPrint({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[1]]]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)
    
    test_normalite(mediator_resid)  # Appeler la fonction de test de normalité
  })

  output$testmedsed <- renderPrint({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[[mediator_vars[2]]]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    req(!is.null(mediator_resid), length(mediator_resid) > 0)
    
    test_normalite(mediator_resid)  # Appeler la fonction de test de normalité
  })
  
  output$testmulti <- renderPrint({
    # Extraire les résidus du médiateur secondaire
    residuals_df <- get_residuals_df()
    
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws()
    mediator_resid <- residuals_df[,mediator_vars]  # Récupérer les résidus du premier médiateur
    # Vérifier que les résidus existent pour le médiateur
    #req(!is.null(mediator_resid), length(mediator_resid) > 0)
   
    test_result <- mvn(mediator_resid, mvnTest = "hz", univariateTest = "SW")  # Appeler la fonction de test de normalité
    print(test_result)
  })

  ####### Corrélation constante
  
  resultsCC <- reactive({
    req(data())
    
    SolVrai <-tableCCVrai <- var_covarVrai <- NULL
    # Extraire les variables d'intérêt
    treat_var <- input$treat
    outcome_var <- input$outcome
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws() # Méditeurs
    ###############################
    inter <- as.logical(input$interaction)
    selected_interactions <- input$interaction_types  # vecteur de valeurs ex: c("treat_m1", "m1_m2")
    interaction_terms <- c()
    
    if (isTRUE(inter) && !is.null(selected_interactions)) {
      if ("treat_m1" %in% selected_interactions) {
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[1], sep = ":"))
      }
      if ("treat_m2" %in% selected_interactions && length(mediator_vars) > 1) {
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[2], sep = ":"))
      }
      if ("m1_m2" %in% selected_interactions && length(mediator_vars) > 1) {
        interaction_terms <- c(interaction_terms, paste(mediator_vars[1], mediator_vars[2], sep = ":"))
      }
      if ("treat_m1_m2" %in% selected_interactions && length(mediator_vars) > 1) {
        triple <- paste(treat_var, mediator_vars[1], mediator_vars[2], sep = ":")
        interaction_terms <- c(interaction_terms, triple)
      }
    }
    ########################################################
    
    
    
    Br <- input$B
    # les covariables de la réponse et des médiateurs
    out_cov_vars <- if (input$out_cov == "") {NULL}else{strsplit(input$out_cov, ",")[[1]] %>% trimws()} 
    intmed_cov_vars <- if (input$intmed_cov == "") {NULL}else{strsplit(input$intmed_cov, ",")[[1]] %>% trimws()}
    sed_cov_vars <- if (input$sed_cov == "") {NULL}else{strsplit(input$sed_cov, ",")[[1]] %>% trimws()}

    covariate <- if (is.null(intmed_cov_vars) && is.null(sed_cov_vars)) {
      NULL  # Si les deux variables sont NULL, retourne NULL
    } else {
      unique(c(intmed_cov_vars, sed_cov_vars))  # Sinon, combine les deux ensembles de variables et retire les doublons
    }
    inter_treat_cov <- FALSE
    if(isTRUE(inter)) {
      lmfor <- paste(c(treat_var, mediator_vars, covariate, interaction_terms), collapse = " + ")
    }else{
     lmfor <- paste(c(treat_var, mediator_vars, covariate), collapse ="+")
    }

    #########################
    if(isTRUE(inter)){
      names_vec <- c(
        paste(treat_var, mediator_vars, sep=":"),
        paste(mediator_vars, collapse = ":"),
        paste(treat_var, paste(mediator_vars, collapse = ":"), sep=":"))
    }else{
      names_vec <- NULL
    }

    cor_cste <- 1
    rh <- 0.5

    q1 <- 2 + length(intmed_cov_vars)
    q2 <- 2 + length(sed_cov_vars)


    if(isTRUE(inter)){
      q3 <- 2 + 3*length(mediator_vars)+ length(out_cov_vars)
    }else{
      q3 <- 2 + length(mediator_vars) + length(out_cov_vars)
    }
    #######################################
  
    tableCC <- initialParams(treat =treat_var, mediators=mediator_vars,intmed = mediator_vars[1],outcome =outcome_var,
                  intmed_cov=intmed_cov_vars, sed_cov=sed_cov_vars, out_cov=out_cov_vars,
                  inter_treat_cov = inter_treat_cov, cor_cste=cor_cste, data=data(),
                  formula_one3=lmfor , names_vec =  names_vec)


    # # # # Estimation des paramètres des médiateurs

    coefs.intmed <- as.vector(tableCC$med[1:q1])
    coefs.sedmed <- as.vector(tableCC$med[(q1+1):(q1+q2)])
    cor_coefs <- as.vector(tableCC$med[-c(1:(q1+q2))])

    # Création de la matrice des coefficients
    coef <- data.frame(matrix(c(coefs.intmed, coefs.sedmed), 2, q1, byrow = TRUE), mediator_vars)
    colnames(coef) <- c("inter", treat_var, intmed_cov_vars, "name")

    # Estimation des paramètres pour Y

    coefs.y <- tableCC$outc

    # Création du tableau des coefficients de Y
    if (isTRUE(inter)) {
      Bet <- data.frame(Beta = coefs.y, name = c("inter", treat_var, mediator_vars, out_cov_vars, names_vec, "sd"))
    } else {
      Bet <- data.frame(Beta = coefs.y, name = c("inter", treat_var, mediator_vars, out_cov_vars, "sd"))
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

    # # Calcul des effets naturels
    sol.i <- effectdirectindirct(alpha = coef, beta = Bet, treat = treat_var, mediators = mediator_vars,
                                 intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, inter = inter, ro = r011,
                                 corC = cor_coefs, names_vec = names_vec, cor_cste = cor_cste, data = data())

    SolVrai <- c(sol.i$DE1sk,sol.i$IE1sk)
    tableCCVrai <- tableCC 
    var_covarVrai <-  var_covar
    bootstrap_list <- vector("list", Br)
    
    set.seed(124)
    b <- 0
    while (b < Br) {  # Tant qu'on n'a pas B échantillons valides

      sample_data <- data()[sample(1:nrow(data()),
                                          size = nrow(data()),
                                          replace = TRUE), ]

      # n1 <- sum(sample_data[, treat_var] == 1)
      # n0 <- sum(sample_data[, treat_var] == 0)
      # 
      # if (n1 > length(intmed_cov_vars) + 1 && n0 > length(sed_cov_vars) + 1 && length(unique(sample_data$x)) > 1) {  # Vérifier la diversité de x
        b <- b + 1  # Incrémenter seulement si l'échantillon est valide
        bootstrap_list[[b]] <- sample_data
      #}
    }
    # 
    SolB <- NULL

    for(j in 1:Br){
      #######################################

      tableCC <- initialParams(treat =treat_var, mediators=mediator_vars,intmed = mediator_vars[1],outcome =outcome_var,
                               intmed_cov=intmed_cov_vars, sed_cov=sed_cov_vars, out_cov=out_cov_vars,
                               inter_treat_cov = inter_treat_cov, cor_cste=cor_cste, data=bootstrap_list[[j]],
                               formula_one3=lmfor , names_vec =  names_vec)


      # # # # Estimation des paramètres des médiateurs

      coefs.intmed <- as.vector(tableCC$med[1:q1])
      coefs.sedmed <- as.vector(tableCC$med[(q1+1):(q1+q2)])
      cor_coefs <- as.vector(tableCC$med[-c(1:(q1+q2))])

      # Création de la matrice des coefficients
      coef <- data.frame(matrix(c(coefs.intmed, coefs.sedmed), 2, q1, byrow = TRUE), mediator_vars)
      colnames(coef) <- c("inter", treat_var, intmed_cov_vars, "name")

      # Estimation des paramètres pour Y

      coefs.y <- tableCC$outc

      # Création du tableau des coefficients de Y
      if (isTRUE(inter)) {
        Bet <- data.frame(Beta = coefs.y, name = c("inter", treat_var, mediator_vars, out_cov_vars, names_vec, "sd"))
      } else {
        Bet <- data.frame(Beta = coefs.y, name = c("inter", treat_var, mediator_vars, out_cov_vars, "sd"))
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

      # # Calcul des effets naturels
      sol.i <- effectdirectindirct(alpha = coef, beta = Bet, treat = treat_var, mediators = mediator_vars,
                                   intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, inter = inter, ro = r011,
                                   corC = cor_coefs, names_vec = names_vec, cor_cste = cor_cste, data = bootstrap_list[[j]])
      SolB <- rbind(SolB, c(sol.i$DE1sk,sol.i$IE1sk))
    }
    
    #SolB <- tt
    # MedCCinit <- Applications_simple_effect(treat= treat_var, outcome= outcome_var,
    #                                         mediators= mediator_vars, intmed = mediator_vars[1],
    #                                         out_cov=out_cov_vars, intmed_cov=intmed_cov_vars,
    #                                         sed_cov= sed_cov_vars, inter=inter, inter_treat_cov = inter_treat_cov,
    #                                         cor_cste = 1, B = Br,rh=0.5,
    #                                         methode = "bootst", data = data(), formula_one3 = lmfor)
    # 
    # # MedCCinit$DirEff
    # # 
    
    SolBi <- ((Br-1)/Br)*apply(SolB, 2, var)
    SolBcant <- apply(SolB, 2, quantile, probs = c(0.025,0.975))
    tibble::tibble(
      #Effets = c("Direct", "Indirect"),
      Valeur = round(SolVrai,2),
      ET = round(SolBi,2),
      IC = paste0("[", round(SolBcant[1, ], 2), ", ", round(SolBcant[2, ], 2), "]")
    )
    # data <- cbind(round(SolVrai,2), round(SolBi,2), paste0("[", round(SolBcant[1, ], 2), ", ", round(SolBcant[2, ], 2), "]"))
    # # Convertir en data frame
    # data <- as.data.frame(data)
    # 
    # # Ajouter les noms des colonnes
    # colnames(data) <- c("Valeur", "ET", "IC")
    # 
    # Convertir les valeurs en caractères (au cas où elles seraient numériques)
    #data$Valeur <- as.character(data$Valeur)
    #data$ET <- as.character(data$ET)
    #$IC <- as.character(data$IC)
    
    # print(class(data))
    # prin
    # str(data)
    # #colnames(data) <- c("Valeurs", "ET","IC")
    # #rownames(data) <- c("Direct", "Indirect")
    # tableau_html(data)  
    # #datatable(data) 
    # #datatable(SolBcant, options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
    #datatable(as.data.frame(t(SolBi)), options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px"))
  })
  
  
  ############## Corrélation non constante ####################################
  
  resultsCNC <- reactive({
    req(data())
    
    SolVrai <-tableCCVrai <- var_covarVrai <- NULL
    # Extraire les variables d'intérêt
    treat_var <- input$treat
    outcome_var <- input$outcome
    mediator_vars <- strsplit(input$mediators, ",")[[1]] %>% trimws() # Méditeurs
    ###############################
    inter <- as.logical(input$interaction)
    selected_interactions <- input$interaction_types  # vecteur de valeurs ex: c("treat_m1", "m1_m2")
    interaction_terms <- c()
    
    if (isTRUE(inter) && !is.null(selected_interactions)) {
      if ("treat_m1" %in% selected_interactions) {
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[1], sep = ":"))
      }
      if ("treat_m2" %in% selected_interactions && length(mediator_vars) > 1) {
        interaction_terms <- c(interaction_terms, paste(treat_var, mediator_vars[2], sep = ":"))
      }
      if ("m1_m2" %in% selected_interactions && length(mediator_vars) > 1) {
        interaction_terms <- c(interaction_terms, paste(mediator_vars[1], mediator_vars[2], sep = ":"))
      }
      if ("treat_m1_m2" %in% selected_interactions && length(mediator_vars) > 1) {
        triple <- paste(treat_var, mediator_vars[1], mediator_vars[2], sep = ":")
        interaction_terms <- c(interaction_terms, triple)
      }
    }
    ########################################################
    
    
    
    Br <- input$B
    # les covariables de la réponse et des médiateurs
    out_cov_vars <- if (input$out_cov == "") {NULL}else{strsplit(input$out_cov, ",")[[1]] %>% trimws()} 
    intmed_cov_vars <- if (input$intmed_cov == "") {NULL}else{strsplit(input$intmed_cov, ",")[[1]] %>% trimws()}
    sed_cov_vars <- if (input$sed_cov == "") {NULL}else{strsplit(input$sed_cov, ",")[[1]] %>% trimws()}
    
    covariate <- if (is.null(intmed_cov_vars) && is.null(sed_cov_vars)) {
      NULL  # Si les deux variables sont NULL, retourne NULL
    } else {
      unique(c(intmed_cov_vars, sed_cov_vars))  # Sinon, combine les deux ensembles de variables et retire les doublons
    }
    inter_treat_cov <- FALSE
    if(isTRUE(inter)) {
      lmfor <- paste(c(treat_var, mediator_vars, covariate, interaction_terms), collapse = " + ")
    }else{
      lmfor <- paste(c(treat_var, mediator_vars, covariate), collapse ="+")
    }
    
    #########################
    if(isTRUE(inter)){
      names_vec <- c(
        paste(treat_var, mediator_vars, sep=":"),
        paste(mediator_vars, collapse = ":"),
        paste(treat_var, paste(mediator_vars, collapse = ":"), sep=":"))
    }else{
      names_vec <- NULL
    }
    
    cor_cste <- 3
    rh <- 0.5
    
    q1 <- 2 + length(intmed_cov_vars)
    q2 <- 2 + length(sed_cov_vars)
    
    
    if(isTRUE(inter)){
      q3 <- 2 + 3*length(mediator_vars)+ length(out_cov_vars)
    }else{
      q3 <- 2 + length(mediator_vars) + length(out_cov_vars)
    }
    #######################################
    
    tableCNC <- initialParams(treat =treat_var, mediators=mediator_vars,intmed = mediator_vars[1],outcome =outcome_var,
                             intmed_cov=intmed_cov_vars, sed_cov=sed_cov_vars, out_cov=out_cov_vars,
                             inter_treat_cov = inter_treat_cov, cor_cste=cor_cste, data=data(),
                             formula_one3=lmfor , names_vec =  names_vec)
    
    
    # # # # Estimation des paramètres des médiateurs
    
    coefs.intmed <- as.vector(tableCNC$med[1:q1])
    coefs.sedmed <- as.vector(tableCNC$med[(q1+1):(q1+q2)])
    cor_coefs <- as.vector(tableCNC$med[-c(1:(q1+q2))])
    
    # Création de la matrice des coefficients
    coef <- data.frame(matrix(c(coefs.intmed, coefs.sedmed), 2, q1, byrow = TRUE), mediator_vars)
    colnames(coef) <- c("inter", treat_var, intmed_cov_vars, "name")
    
    # Estimation des paramètres pour Y
    
    coefs.y <- tableCNC$outc
    
    # Création du tableau des coefficients de Y
    if (isTRUE(inter)) {
      Bet <- data.frame(Beta = coefs.y, name = c("inter", treat_var, mediator_vars, out_cov_vars, names_vec, "sd"))
    } else {
      Bet <- data.frame(Beta = coefs.y, name = c("inter", treat_var, mediator_vars, out_cov_vars, "sd"))
    }
    
    # Calcul des covariances et corrélations
    var_covar <- varcovarEstimes(cor_coefs, cor_cste = cor_cste)
    
    
    r011 <- c(as.vector(var_covar$r01), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
    r100 <- c(as.vector(var_covar$r10), (cor_coefs[5] + cor_coefs[6]) / 2, rh)
    
    # # Calcul des effets naturels
    sol.i <- effectdirectindirct(alpha = coef, beta = Bet, treat = treat_var, mediators = mediator_vars,
                                 intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, inter = inter, ro = r011,
                                 corC = cor_coefs, names_vec = names_vec, cor_cste = 3, data = data())
    
    SolVrai <- c(sol.i$DE1sk,sol.i$IE1sk)
    tableCNCVrai <- tableCNC 
    var_covarVrai <-  var_covar
    bootstrap_list <- vector("list", Br)
    
    set.seed(145)
    b <- 0
    while (b < Br) {  # Tant qu'on n'a pas B échantillons valides
      
      sample_data <- data()[sample(1:nrow(data()),
                                   size = nrow(data()),
                                   replace = TRUE), ]
      
      # n1 <- sum(sample_data[, treat_var] == 1)
      # n0 <- sum(sample_data[, treat_var] == 0)
      # 
      # if (n1 > length(intmed_cov_vars) + 1 && n0 > length(sed_cov_vars) + 1 && length(unique(sample_data$x)) > 1) {  # Vérifier la diversité de x
      b <- b + 1  # Incrémenter seulement si l'échantillon est valide
      bootstrap_list[[b]] <- sample_data
      #}
    }
    # 
    SolB <- NULL
    
    for(j in 1:Br){
      #######################################
      
      tableCNCB <- initialParams(treat =treat_var, mediators=mediator_vars,intmed = mediator_vars[1],outcome =outcome_var,
                               intmed_cov=intmed_cov_vars, sed_cov=sed_cov_vars, out_cov=out_cov_vars,
                               inter_treat_cov = inter_treat_cov, cor_cste=3, data=bootstrap_list[[j]],
                               formula_one3=lmfor , names_vec =  names_vec)
      
      
      # # # # Estimation des paramètres des médiateurs
      
      coefs.intmedB <- as.vector(tableCNCB$med[1:q1])
      coefs.sedmedB <- as.vector(tableCNCB$med[(q1+1):(q1+q2)])
      cor_coefsB <- as.vector(tableCNCB$med[-c(1:(q1+q2))])
      
      # Création de la matrice des coefficients
      coefB <- data.frame(matrix(c(coefs.intmedB, coefs.sedmedB), 2, q1, byrow = TRUE), mediator_vars)
      colnames(coefB) <- c("inter", treat_var, intmed_cov_vars, "name")
      
      # Estimation des paramètres pour Y
      
      coefs.yB <- tableCNCB$outc
      
      # Création du tableau des coefficients de Y
      if (isTRUE(inter)) {
        BetB <- data.frame(Beta = coefs.yB, name = c("inter", treat_var, mediator_vars, out_cov_vars, names_vec, "sd"))
      } else {
        BetB <- data.frame(Beta = coefs.yB, name = c("inter", treat_var, mediator_vars, out_cov_vars, "sd"))
      }
      
      # Calcul des covariances et corrélations
      var_covarB <- varcovarEstimes(cor_coefsB, cor_cste = cor_cste)
     
      r011B <- c(as.vector(var_covarB$r01), (cor_coefsB[5] + cor_coefsB[6]) / 2, rh)
      r100B <- c(as.vector(var_covarB$r10), (cor_coefsB[5] + cor_coefsB[6]) / 2, rh)

      
      # # Calcul des effets naturels
      sol.iB <- effectdirectindirct(alpha = coefB, beta = BetB, treat = treat_var, mediators = mediator_vars,
                                   intmed_cov = intmed_cov_vars, sed_cov = sed_cov_vars, inter = inter, ro = r011B,
                                   corC = cor_coefsB, names_vec = names_vec, cor_cste = 3, data = bootstrap_list[[j]])
      SolB <- rbind(SolB, c(sol.iB$DE1sk,sol.iB$IE1sk))
    }
    
    SolBi <- ((Br-1)/Br)*apply(SolB, 2, var, na.rm=TRUE)
    SolBcant <- apply(SolB, 2, quantile, probs = c(0.025,0.975), na.rm=TRUE)
    dtf <- tibble::tibble(
      Effet = rep(c("\\(\\zeta\\)", "\\(\\delta\\)"), each = 6),
      Méthode = rep(c(rep("CNCr", 4), "CNCm", "CNC"), 2),
      rho = rep(round(r011,2),2),
      Valeur = round(SolVrai,2),
      ET = round(SolBi,2),
      IC = paste0("[", round(SolBcant[1, ], 2), ", ", round(SolBcant[2, ], 2), "]")
    )
    dtf_modifié <- dtf %>%
      slice(c(5,1,2,3,4,11,7,8,9,10))
    #datatable(dtf_modifié)
  })
  
  output$resultsCC <- renderUI({
    df <- resultsCC()
    df2 <- resultsCNC()
    tableau_html(df, df2)
  })
  
}


# Launch Dashboard
shinyApp(ui, server)  # Initialize the shinydashboard




