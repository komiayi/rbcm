

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
