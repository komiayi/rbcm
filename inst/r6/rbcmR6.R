library(Rcpp)
library(RcppArmadillo)
sourceCpp('logmed_rccp.cpp')

RbCorelatedMediators <- R6::R6Class(
  "Rbcm",
  public = list(
    params = list(),
    estimates_params = list(),
    data = NULL,
    mediators_model_fit = NULL,
    outcome_model_fit = NULL,
    method = NULL,
    mediators_interactions =NULL,
    outcome_interactions =NULL,
    direct = NULL,
    indirect = NULL,
    initialize = function(data, mediators, response, exposition, 
                          covariates = NULL) {
      
      if (!is.null(dim(data)) && length(dim(data)) > 2) {
        stop("L'argument 'data' ne peut pas être un Array de dimension 3 ou plus. 
          Veuillez fournir un data.frame ou une matrice 2D.", call. = FALSE)
      }
      
      if (is.list(data) && !is.data.frame(data)) {
        # On vérifie si les éléments de la liste ont tous la même longueur
        lengths <- sapply(data, length)
        if (length(unique(lengths)) > 1) {
          stop("L'argument 'data' est une liste dont les éléments ont des longueurs différentes.
            Elle ne peut pas être convertie en tableau de données.", call. = FALSE)
        }
      }
      
      df_temp <- tryCatch({
        as.data.frame(data)
      }, error = function(e) {
        stop("Format de données invalide : impossible de convertir en data.frame.", call. = FALSE)
      })
      
      self$data <- as.matrix(df_temp)
      self$params$mediators <- mediators
      self$params$response <- response
      self$params$exposition <- exposition
      self$params$covariates <- covariates
      
      private$validate_columns()
      
      private$T_val <- self$data[, self$params$exposition[['exposure']]]
      private$M11   <- self$data[, self$params$mediators[['primary']],   drop = FALSE]
      private$M22   <- self$data[, self$params$mediators[['secondary']], drop = FALSE]
      private$Y     <- self$data[, self$params$response[['outcome']],    drop = FALSE]
      
      private$df_cov_primary <- if(!is.null(self$params$covariates[['cov_primary']])){
        self$data[, self$params$covariates[['cov_primary']], drop = FALSE]
      } else {
        NULL
      }
      
      private$df_cov_secondary <- if(!is.null(self$params$covariates[['cov_secondary']])){
        self$data[, self$params$covariates[['cov_secondary']], drop = FALSE]
      } else { NULL }
      
      private$df_cov_outc <- if(!is.null(self$params$covariates[['cov_outcome']])){
        self$data[, self$params$covariates[['cov_outcome']], drop = FALSE]
      } else { NULL }
      
      private$idx0 <- which(private$T_val == 0)
      private$idx1 <- which(private$T_val == 1)
    },
    parmetors_estimators = function(method = 'CC',
                                    interaction = list(primary = TRUE, secondary = TRUE),
                                    outc_interactions = NULL){
      T_val <- private$T_val
      M1    <- private$M11
      M2    <- private$M22
      idx0  <- private$idx0
      idx1  <- private$idx1
      Y <-  private$Y
      Z1 <- private$build_med_design_matrix(T_val=T_val, cov_data=private$df_cov_primary, interaction = interaction[['primary']])
      Z2 <- private$build_med_design_matrix(T_val=T_val, cov_data=private$df_cov_secondary, interaction = interaction[['secondary']])
      Y1 <- private$build_outc_design_matrix(T_val = T_val, M1= M1, M2=M2,
                                             cov_data = private$df_cov_outc,
                                             interactions=outc_interactions)
      ################### optimisation mediateurs  #############################
      init_results <- private$initial_parms(method = method, Z1=Z1, Z2=Z2, M1=M1,
                                            M2=M2, idx0=idx0,
                                            idx1= idx1)
      theta_init <- init_results$theta_init
      p1 <- ncol(Z1)
      p2 <- ncol(Z2)
      optim_bounds <- private$get_bounds(method = method, p1 = p1, p2 = p2)
      lower_bounds <- optim_bounds$lower
      upper_bounds <- optim_bounds$upper
      
      self$mediators_model_fit <- (function() {
        res <- optim(
          par     = theta_init,
          fn      = log_lik_med_arma,
          M1      = M1,
          M2      = M2,
          Z1      = Z1, Z2 = Z2, T_val = T_val,
          idx0    = idx0,
          idx1    = idx1,
          methods = method,
          return_sum = TRUE,
          lower   = lower_bounds,
          upper   = upper_bounds,
          method  = "L-BFGS-B"#,
          #hessian = TRUE
        )
        list(
          par = res$par,
          #hessian = res$hessian,
          loglik = -res$value,
          converged = res$convergence == 0
        )
      })()
      
      ###################### optimisationoutcome ###############################
      beta_init_results <- get_initials_arma(Y1, Y)
      beta_init <- c(beta_init_results$alphas, log(sd(beta_init_results$res)))
      
      q <- ncol(Y1)
      lower_bounds_outc <- c(rep(-Inf,q),-10)
      upper_bounds_outc <- c(rep(Inf,q),10)
      
      self$outcome_model_fit <- (function() {
        res <- optim(
          par     = beta_init,
          fn      = log_lik_outc_arma,
          Y      = Y, Y1 = Y1,
          return_sum = TRUE,
          lower   = lower_bounds_outc,
          upper   = upper_bounds_outc,
          method  = "L-BFGS-B"#,
          #hessian = TRUE
        )
        list(
          par = res$par,
          #hessian = res$hessian,
          loglik = -res$value,
          converged = res$convergence == 0
        )
      })()
      
      self$method <- method
      self$mediators_interactions <- interaction
      self$outcome_interactions <- outc_interactions
      
      ########################  affichage des coefficients estimés #############
      private$estimates_results()
      
      # NETTOYAGE FINAL
      # rm(med_fit_mle, outc_fit_mle)
      # gc() # On force la libération de la RAM avant de sortir de la fonction
      invisible(self)
    },
    
    ## calcul des effets direct et indirect
    effests_Estimate = function(){
      if(!is.null(private$df_cov_primary)){
        a <-  tcrossprod(private$df_cov_primary, t(self$estimates_params$primary$alpha2))  
        b <-  tcrossprod(private$df_cov_primary, t(self$estimates_params$primary$alpha3)) 
      }else{
        a <- b <- 0 
      }
      
      if(!is.null(private$df_cov_secondary)){
        c <- tcrossprod(private$df_cov_secondary , t(self$estimates_params$secondary$alpha3)) 
        d <- tcrossprod(private$df_cov_secondary , t(self$estimates_params$secondary$alpha2))
      }else{
        c <- d <- 0 
      }
      
      V01 = self$estimates_params$primary$alpha0 + a
      V11 = self$estimates_params$primary$alpha1 + b
      Vbase = self$estimates_params$secondary$alpha1 + c
      V02 = Vbase + self$estimates_params$secondary$alpha0  + d
      
      reslt_cor <- self$estimates_params$Cor
      rho01 <- self$estimates_params$Rho01
      if(self$method == 'CNC'){
        rho01 <- c(rho01, -rho01)
        #self$estimates_params$primary$sigma*self$estimates_params$secondary$sigma
        rho0 <- self$estimates_params$primary$sigma[[1]]*self$estimates_params$secondary$sigma[[1]]*reslt_cor[2]
        rho1 <- self$estimates_params$primary$sigma[[2]]*self$estimates_params$secondary$sigma[[2]]*reslt_cor[1]
        rho22 <- self$estimates_params$primary$sigma[[2]]*self$estimates_params$secondary$sigma[[1]]*rho01
        Rho0 <- rho22- rho0
        Rho1 <- rho1 - rho22
      }else{
        rho22 <- self$estimates_params$primary$sigma*self$estimates_params$secondary$sigma*self$estimates_params$Cor
        Rho0 <- Rho1 <- 0
      }
      
      #############################
      direct  <- sapply(1:length(rho22), function (i){
        self$estimates_params$outcome$beta1 + self$estimates_params$outcome$beta2[[2]]*Vbase + self$estimates_params$outcome$beta3[[1]]*V01+ self$estimates_params$outcome$beta3[[2]]*V02+
          self$estimates_params$outcome$beta3[[3]]*(V01*Vbase + Rho0[i])+ self$estimates_params$outcome$beta3[[4]]*(V01*V02 + rho22[i])
      })
      indirect  <- sapply(1:length(rho22), function (i){
        (self$estimates_params$outcome$beta2[[1]] + self$estimates_params$outcome$beta3[[1]])*V11 + (self$estimates_params$outcome$beta3[[3]]+ self$estimates_params$outcome$beta3[[4]])*(V11*V02 + Rho1[i])
      })
      
      Direct <- if(is.matrix(direct)) colMeans(direct) else direct
      InDirect <- if(is.matrix(indirect)) colMeans(indirect) else indirect
      
      self$direct <- Direct
      self$indirect <- InDirect
      
      # --- AFFICHAGE ---
      cat("========================================================\n")
      cat("       RÉSULTAT EFFETS DE MÉDIATION (Méthode : ", self$method, ") \n")
      cat("========================================================\n\n")
      
      cat("Effet(s) direct  :\n"); print(self$direct)
      cat("\n")
      cat("Effet(s) indirect  :\n"); print(self$indirect)
      
      cat("\n========================================================\n")
      
      invisible(self)
    },
    varcovarH = function(){
      
      T_val <- private$T_val
      M1    <- private$M11
      M2    <- private$M22
      idx0  <- private$idx0
      idx1  <- private$idx1
      Y <-  private$Y
      
      Z1 <- private$build_med_design_matrix(T_val=T_val, cov_data=private$df_cov_primary, interaction = interaction[['primary']])
      Z2 <- private$build_med_design_matrix(T_val=T_val, cov_data=private$df_cov_secondary, interaction = interaction[['secondary']])
      Y1 <- private$build_outc_design_matrix(T_val = T_val, M1= M1, M2=M2,
                                             cov_data = private$df_cov_outc,
                                             interactions=outc_interactions)
      
      skeleton <- ss$estimates_params
      
      if("Rho01" %in% names(skeleton)) skeleton$Rho01 <- NULL
      
      skeleton_med <- skeleton[setdiff(names(skeleton), "outcome")]
      skeleton_outc <- skeleton$outcome
      
      theta0 <- unlist(skeleton)
      theta_med <- unlist(skeleton_med)
      theta_outc <- unlist(skeleton_outc)
      
      psi_mat_med <- t(sapply(1:length(T_val), function(i){
        numDeriv::grad(function(th) log_lik_med_arma(th, M1 = M1, M2 =M2, Z1=Z1,
                                                     Z2=Z2, T_val=T_val, 
                                                     idx0=idx0, idx1=idx1,
                                                     methods=self$method, return_sum=FALSE)[i],
                       theta_med)
      }))
      A_med <- numDeriv::grad(function(th) colMeans(psi_mat_med), theta_med)
      
      
      psi_mat_outc <- t(sapply(1:length(T_val), function(i){
        numDeriv::grad(function(th) log_lik_outc_arma(th, Y =Y, Y1=Y1, return_sum=FALSE)[i],
                       theta_outc)
      }))
      A_outc <- numDeriv::grad(function(th) colMeans(psi_mat_outc), theta_outc)
      
      psi_total <- c(psi_mat_med, psi_mat_outc)
      
      #B <- t(psi_mat) %*% psi_mat / n
      B <- crossprod(psi_total)/ n
      
      A <- c(A_med, A_outc)
      V_theta <- solve(A) %*% B %*% t(solve(A))
      
      ############  direct ###############
      
      grad_val_direct <- if(length(self$direct)>1){
        numDeriv::jacobian(
          func = private$zeta_wrapper_direct,
          x = theta0,
          skeleton = skeleton,
          df_cov_primary = private$df_cov_primary,
          df_cov_secondary = private$df_cov_secondary,
          methode = self$method
        )
      }else{
        numDeriv::grad(
          func = private$zeta_wrapper_direct,
          x = theta0,
          skeleton = skeleton,
          df_cov_primary = private$df_cov_primary,
          df_cov_secondary = private$df_cov_secondary,
          methode = self$method
        )
      }
      
      Var_direct <- grad_val_direct %*% V_theta %*% t(grad_val_direct)
      
      ############  direct ###############
      
      grad_val_indirect <- if(length(self$indirect)>1){
        numDeriv::jacobian(
          func = private$zeta_wrapper_indirect,
          x = theta0,
          skeleton = skeleton,
          df_cov_primary = private$df_cov_primary,
          df_cov_secondary = private$df_cov_secondary,
          methode = self$method
        )
      }else{
        numDeriv::grad(
          func = private$zeta_wrapper_indirect,
          x = theta0,
          skeleton = skeleton,
          df_cov_primary = private$df_cov_primary,
          df_cov_secondary = private$df_cov_secondary,
          methode = self$method
        )
      }
      
      Var_indirect <- grad_val_indirect %*% V_theta %*% t(grad_val_indirect)
      invisible(self)
    }
  ),
  active = list(
    ## mettre ue fonction qui modifie "method"
    ## parce que sur la même base on peut vouloir utiliser les deux méthodes
  ),
  private = list(
    T_val = NULL,
    M11 = NULL,
    M22 = NULL,
    Y = NULL,
    df_cov_primary = NULL,
    df_cov_secondary = NULL,
    df_cov_outc = NULL,
    idx0 = NULL,
    idx1 = NULL,
    validate_columns = function() {
      required_cols <- unlist(list(
        self$params$mediators,
        self$params$response,
        self$params$exposition,
        self$params$covariates
      ))
      
      # Supprimer les doublons éventuels et les valeurs NULL
      required_cols <- unique(required_cols[!is.na(required_cols)])
      
      # Vérifier la présence dans les données
      existing_cols <- colnames(self$data)
      missing_cols <- setdiff(required_cols, existing_cols)
      
      # Si des colonnes manquent, on arrête tout avec un message explicite
      if (length(missing_cols) > 0) {
        stop(paste0(
          "\n[Erreur de validation] Les colonnes suivantes sont introuvables dans la base de données :\n",
          paste("- ", missing_cols, collapse = "\n")
        ), call. = FALSE)
      }
      
      message("✓ Toutes les colonnes ont été vérifiées avec succès.")
    },
    build_med_design_matrix = function(T_val, cov_data = NULL, interaction = FALSE) {
      p <- if(is.null(cov_data)) 0 else ncol(cov_data)
      n <- length(T_val)
      cols <- 1 + 1 + p + if(interaction) p else 0  # intercept + T_val + covs + interactions
      Z <- matrix(0, n, cols)
      
      Z[,1] <- 1                 # intercept
      Z[,2] <- T_val
      if(p > 0) Z[,3:(2+p)] <- cov_data
      if(interaction) Z[,(3+p):(2+2*p)] <- T_val * cov_data
      return(Z)
    },
    build_outc_design_matrix = function(T_val, M1, M2, cov_data = NULL,
                                        interactions=NULL) {
      n <- length(T_val)
      p_cov <- if(is.null(cov_data)) 0 else ncol(cov_data)
      p_inter <- if(is.null(interactions)) 0 else length(interactions) 
      cols <- 4 + p_inter + p_cov 
      Z <- matrix(0, n, cols)
      # Colonnes fixes
      col_names <- character(cols)
      Z[, 1] <- 1
      Z[, 2] <- T_val
      Z[, 3] <- M1
      Z[, 4] <- M2
      
      col_names[1] <- "Intercept"
      col_names[2] <- if(!is.null(names(T_val))) names(T_val) else deparse(substitute(T_val))
      col_names[3] <- if(!is.null(names(M1))) names(M1) else deparse(substitute(M1))
      col_names[4] <- if(!is.null(names(M2))) names(M2) else deparse(substitute(M2))
      col_idx <- 5
      if (!is.null(interactions)) {
        inter_names_map <- c(
          "1" = paste(col_names[2], col_names[3], sep = "_"),
          "2" = paste(col_names[2], col_names[4], sep = "_"),
          "3" = paste(col_names[3], col_names[4], sep = "_"),
          "4" = paste(col_names[2], col_names[3], col_names[4], sep = "_")
        )
        for (choice in interactions) {
          Z[, col_idx] <- switch(as.character(choice),
                                 "1" = T_val * M1,
                                 "2" = T_val * M2,
                                 "3" = M1 * M2,
                                 "4" = T_val * M1 * M2)
          col_names[col_idx] <- inter_names_map[as.character(choice)]
          col_idx <- col_idx + 1
        }
      }
      if (!is.null(cov_data)) {
        Z[, col_idx:(col_idx + p_cov - 1)] <- cov_data
        col_names[col_idx:(col_idx + p_cov - 1)] <- colnames(cov_data)
      }
      
      colnames(Z) <- c(col_names)
      return(Z)
    },
    initial_parms = function(method, Z1,Z2,M1,M2, idx1 = NULL, idx0  = NULL) {
      # res1 <- private$get_initials(Z1, M1)
      # res2 <- private$get_initials(Z2, M2)
      res1 <- get_initials_arma(Z1, M1)
      res2 <- get_initials_arma(Z2, M2)
      
      
      if (method == "CC") {
        # --- CAS CONSTANT (CC) ---
        s1_val <- sd(res1$res)
        s2_val <- sd(res2$res)
        rho_val <- cor(res1$res, res2$res)
        
        # On répète les mêmes valeurs pour T=0 et T=1
        sigmas_init <- log(c(s1_val, s2_val))
        rhos_init   <- atanh(c(rho_val))
        rho01_val <- NULL
      } else {
        # --- CAS NON-CONSTANT (CNC) ---
        # Groupe T=0
        s1_0 <- sd(res1$res[idx0])
        s2_0 <- sd(res2$res[idx0])
        rho_00 <- cor(res1$res[idx0], res2$res[idx0])
        
        # Groupe T=1
        s1_1 <- sd(res1$res[idx1])
        s2_1 <- sd(res2$res[idx1])
        rho_11 <- cor(res1$res[idx1], res2$res[idx1])
        
        sigmas_init <- log(c(s1_1, s1_0, s2_1, s2_0))
        rhos_init   <- atanh(c(rho_00, rho_11))
        rho01_val <- private$calculate_rho01(s1_1, s1_0, s2_1, s2_0, rho_00, rho_11)
        
      }
      
      # Assemblage final
      theta_init <- c(res1$alphas, res2$alphas, sigmas_init, rhos_init)
      
      # Sécurité anti-NA (si un groupe est trop petit pour cor() ou sd())
      #theta_init[is.na(theta_init)] <- 0
      
      return(list(theta_init = theta_init, rho01 = rho01_val))
    },
    calculate_rho01 = function(s1_1, s1_0, s2_1, s2_0, rho_00, rho_11) {
      
      # Calcul des termes intermédiaires (variances)
      v1_0 <- s1_0^2; v1_1 <- s1_1^2
      v2_0 <- s2_0^2; v2_1 <- s2_1^2
      
      # Calcul de D12 (Formule du document)
      D12 <- -(v1_1 * v2_1 * rho_11^2) - (v1_0 * v2_0 * rho_00^2) + (v1_1 - v1_0) * (v2_1 - v2_0)
      
      # Calcul de Delta12
      Delta12 <- D12^2 - 4 * (v1_1 * v1_0 * v2_1 * v2_0 * rho_11^2 * rho_00^2)
      
      # --- Gestion des solutions ---
      
      # Cas 1 : Pas de solution réelle (Delta < 0)
      if (is.na(Delta12) || Delta12 < 0) {
        warning("Delta12 < 0. Initialisation à 0.")
        #return(0)
        Delta12 <- 0
      }
      
      # Cas 2 : Calcul des deux racines potentielles pour rho^2
      sol_plus  <- (-D12 + sqrt(Delta12)) / (2 * v1_0 * v2_1)
      sol_moins <- (-D12 - sqrt(Delta12)) / (2 * v1_0 * v2_1)
      
      # On ne garde que les solutions où rho^2 >= 0 (car on va prendre la racine carrée)
      candidates_sq <- c()
      if (!is.na(sol_plus) && sol_plus >= 0) candidates_sq <- c(candidates_sq, sol_plus)
      if (!is.na(sol_moins) && sol_moins >= 0) candidates_sq <- c(candidates_sq, sol_moins)
      
      if (length(candidates_sq) == 0) return(0)
      
      # On prend la racine carrée pour obtenir rho(0,1)
      # On choisit généralement la solution qui reste < 1
      all_rhos <- sqrt(candidates_sq)
      valid_rhos <- all_rhos[all_rhos <= 1]
      
      if (length(valid_rhos) == 0) {
        return(min(0.99, max(all_rhos))) # Sécurité si tout est > 1
      }
      
      # Par défaut, on prend la plus grande valeur valide (souvent la plus stable)
      return(unique(valid_rhos))
    },
    get_bounds = function(method = "CC", p1, p2) {
      # nl : nombre de sigmas (2 en CC, 4 en CNC)
      nl <- if (method == "CC") 2 else 4
      
      # nr : nombre de rhos (1 en CC, 2 en CNC)
      nr <- if (method == "CC") 1 else 2
      
      # 1. Alphas (p1 pour M1, p2 pour M2)
      low_alphas <- rep(-Inf, p1 + p2)
      upp_alphas <- rep(Inf, p1 + p2)
      
      
      low_sigmas <- rep(-10, nl)
      upp_sigmas <- rep(10, nl)
      
      low_rhos <- rep(-3, nr) 
      upp_rhos <- rep(3, nr)
      
      # Assemblage
      lower <- c(low_alphas, low_sigmas, low_rhos)
      upper <- c(upp_alphas, upp_sigmas, upp_rhos)
      
      return(list(lower = lower, upper = upper))
    },
    split_and_extract_alphas = function(big_alpha_vec, a1,a2, b1,b2,c) {
      # --- MÉDIATEUR 1 (Primary) ---
      p1 <- length(a1)
      # Taille L1 : Intercept(1) + Expo(1) + Cov(p1) + Interaction(p1 si TRUE)
      len_m1 <- 2 + p1 + (if(a2) p1 else 0)
      alpha_m1_raw <- big_alpha_vec[1:len_m1]
      
      # --- MÉDIATEUR 2 (Secondary) ---
      p2 <- length(b1)
      # Taille L2 : Intercept(1) + Expo(1) + Cov(p2) + Interaction(p2 si TRUE)
      len_m2 <- 2 + p2 + (if(b2) p2 else 0)
      # On commence après M1
      alpha_m2_raw <- big_alpha_vec[(len_m1 + 1):(len_m1 + len_m2)]
      
      # --- EXTRACTION DÉTAILLÉE ---
      # On réutilise la logique de découpage par blocs
      return(list(
        primary = private$extract_blocks(alpha_m1_raw, p1, a2, a1, c),
        secondary  = private$extract_blocks(alpha_m2_raw, p2, b2, b1, c),
        corre = big_alpha_vec[(len_m1 +len_m2+ 1):length(big_alpha_vec)]
      ))
    },
    extract_blocks = function(vec, p, has_int, cov_list, a) {
      res <- list(
        alpha0 = vec[1],
        alpha1 = vec[2]
      )
      if (p > 0) {
        res$alpha2 <- setNames(vec[3:(2 + p)], cov_list)
      } else {
        res$alpha2 <- NULL
      }
      if (has_int && p > 0) {
        res$alpha3 <- setNames(vec[(3 + p):(2 + 2 * p)], 
                               paste0(a, ":", cov_list))
      } else {
        res$alpha3 <- NULL
      }
      return(res)
    },
    extract_betas = function(beta_vec, a) {
      res <- list(
        beta0  = beta_vec[1],
        beta1  = beta_vec[2],
        beta21 = beta_vec[3],
        beta22 = beta_vec[4]
      )
      res$beta31 <- 0 
      res$beta32 <- 0 
      res$beta33 <- 0 
      res$beta34 <- 0 
      if (!is.null(a)) {
        for (i in seq_along(a)) {
          choice <- a[i]
          val    <- beta_vec[4 + i]
          if (choice == "1") res$beta31 <- val
          if (choice == "2") res$beta32 <- val
          if (choice == "3") res$beta33 <- val
          if (choice == "4") res$beta34 <- val
        }
      }
      
      return(res)
    },
    estimates_results = function() {
      #if (compute) {
      results_med <- private$split_and_extract_alphas(self$mediators_model_fit$par,
                                                      a1=self$params$covariates$cov_primary,
                                                      a2=self$mediators_interactions$primary,
                                                      b1=self$params$covariates$cov_secondary,
                                                      b2=self$mediators_interactions$secondary,
                                                      c=self$params$exposition$exposure)
      reslt_cor <- results_med$corre 
      
      # --- PRÉPARATION MÉDIATEUR PRIMAIRE (M1) ---
      self$estimates_params$primary$alpha0 <- as.vector(results_med$primary$alpha0)
      self$estimates_params$primary$alpha1 <- as.vector(results_med$primary$alpha1)
      
      alpha2_m1 <- as.vector(results_med$primary$alpha2)
      names(alpha2_m1) <- self$params$covariates$cov_primary
      self$estimates_params$primary$alpha2 <- alpha2_m1
      
      alpha3_m1 <- as.vector(results_med$primary$alpha3)
      names(alpha3_m1) <- if(!is.null(results_med$primary$alpha3)) paste0(self$params$exposition$exposure, ":", self$params$covariates$cov_primary) else NULL
      self$estimates_params$primary$alpha3 <- alpha3_m1
      
      # --- PRÉPARATION MÉDIATEUR SECONDAIRE (M2) ---
      self$estimates_params$secondary$alpha0 <- as.vector(results_med$secondary$alpha0)
      self$estimates_params$secondary$alpha1 <- as.vector(results_med$secondary$alpha1)
      
      alpha2_m2 <- as.vector(results_med$secondary$alpha2)
      names(alpha2_m2) <- self$params$covariates$cov_secondary
      self$estimates_params$secondary$alpha2 <- alpha2_m2
      
      alpha3_m2 <- as.vector(results_med$secondary$alpha3)
      names(alpha3_m2) <- if(!is.null(results_med$secondary$alpha3)) paste0(self$params$exposition$exposure, ":", self$params$covariates$cov_secondary) else NULL
      self$estimates_params$secondary$alpha3 <- alpha3_m2
      
      # --- VENTILATION DES SIGMAS (VARIANCE) ---
      if (self$method == 'CNC') {
        # Ordre supposé : sig11, sig10, sig21, sig20, rho00, rho11
        self$estimates_params$primary$sigma  <- exp(reslt_cor[1:2]) 
        names(self$estimates_params$primary$sigma) <- c("Sigma (T=1)", "Sigma (T=0)")
        
        self$estimates_params$secondary$sigma <- exp(reslt_cor[3:4])
        names(self$estimates_params$secondary$sigma) <- c("Sigma (T=1)", "Sigma (T=0)")
        
        self$estimates_params$Cor <- tanh(reslt_cor[5:6]) # rho00 et rho11
        self$estimates_params$Rho01 <- private$calculate_rho01(
          s1_1=exp(reslt_cor[1]), s1_0=exp(reslt_cor[2]),
          s2_1=exp(reslt_cor[3]), s2_0=exp(reslt_cor[4]), 
          rho_00=tanh(reslt_cor[5]), rho_11=tanh(reslt_cor[6])
        )
      } else {
        # Method CC : sig1, sig2, rho
        self$estimates_params$primary$sigma   <- exp(reslt_cor[1])
        self$estimates_params$secondary$sigma <- exp(reslt_cor[2])
        self$estimates_params$Cor             <- tanh(reslt_cor[3]) # rho unique
        self$estimates_params$Rho01           <- NULL
      }
      
      results_outc <- private$extract_betas(self$outcome_model_fit$par, self$outcome_interactions)
      # --- PRÉPARATION OUTCOME (Betas) ---
      self$estimates_params$outcome$beta0 <- results_outc$beta0
      self$estimates_params$outcome$beta1 <- results_outc$beta1
      
      # 1. Noms pour les Covariables de l'Outcome (beta2)
      beta2_vals <- c(results_outc$beta21, results_outc$beta22)
      names(beta2_vals) <- self$params$mediators
      self$estimates_params$outcome$beta2 <- beta2_vals
      
      # 2. Noms pour les Médiateurs et leurs Interactions (beta3)
      # On récupère les noms des variables pour construire les labels
      
      M1_name <- self$params$mediators$primary
      M2_name <- self$params$mediators$secondary
      T_name  <- self$params$exposition$exposure
      
      # Création du vecteur de noms selon votre logique 1, 2, 3, 4
      # Ordre standard supposé : [M1, M2, Inter_1, Inter_2, Inter_3, Inter_4]
      #beta3_labels <- c(M1_name, M2_name)
      map_inter <- list(
        "1" = paste0(T_name, ":", M1_name),                # Exp-M1
        "2" = paste0(T_name, ":", M2_name),                # Exp-M2
        "3" = paste0(M1_name, ":", M2_name),               # M1-M2
        "4" = paste0(T_name, ":", M1_name, ":", M2_name)      # Exp-M1-M2
      )
      
      active_inter_labels <- unlist(map_inter[self$outcome_interactions])
      beta3_labels <- active_inter_labels
      # On regroupe les résultats de beta3 (assurez-vous que results_outc les contient dans cet ordre)
      beta3_vals <- c(results_outc$beta31, results_outc$beta32, results_outc$beta33, results_outc$beta34)
      
      # Sécurité : on vérifie que la taille correspond avant de nommer
      if(length(beta3_vals) == length(beta3_labels)) {
        names(beta3_vals) <- beta3_labels
      }
      self$estimates_params$outcome$beta3 <- beta3_vals
      
      # --- AFFICHAGE ---
      cat("========================================================\n")
      cat("       RÉSUMÉ DU MODÈLE (Méthode : ", self$method, ") \n")
      cat("========================================================\n\n")
      
      # Section M1
      cat(paste0("--- MÉDIATEUR PRIMAIRE (", self$params$mediators$primary, ") ---\n"))
      cat(sprintf("Intercept  : %f\n", self$estimates_params$primary$alpha0))
      cat(sprintf("Exposition : %f\n", self$estimates_params$primary$alpha1))
      cat("Covariables :\n"); print(self$estimates_params$primary$alpha2)
      cat("Interactions :\n"); print(self$estimates_params$primary$alpha3)
      cat("Variances (Sigma) :\n"); print(self$estimates_params$primary$sigma)
      cat("\n")
      
      # Section M2
      cat(paste0("--- MÉDIATEUR SECONDAIRE (", self$params$mediators$secondary, ") ---\n"))
      cat(sprintf("Intercept  : %f\n", self$estimates_params$secondary$alpha0))
      cat(sprintf("Exposition : %f\n", self$estimates_params$secondary$alpha1))
      cat("Covariables :\n"); print(self$estimates_params$secondary$alpha2)
      cat("Interactions :\n"); print(self$estimates_params$secondary$alpha3)
      cat("Variances (Sigma) :\n"); print(self$estimates_params$secondary$sigma)
      cat("\n")
      
      # Section Corrélation
      cat("--- CORRÉLATION ET PARAMÈTRES DE STRUCTURE ---\n")
      if(self$method == 'CNC') {
        cat(sprintf("Rho_00 (T=0) : %f\n", (self$estimates_params$Cor[1])))
        cat(sprintf("Rho_11 (T=1) : %f\n", (self$estimates_params$Cor[2])))
        cat(sprintf("Rho_01 (cross) : %f\n", c(self$estimates_params$Rho01,-self$estimates_params$Rho01)))
      } else {
        cat(sprintf("Corrélation résiduelle (Rho) : %f\n", (self$estimates_params$Cor)))
      }
      cat("\n")
      
      # Section Outcome
      cat(paste0("--- RÉPONSE / OUTCOME (", self$params$response$outcome, ") ---\n"))
      cat(sprintf("Intercept  : %f\n", self$estimates_params$outcome$beta0))
      cat(sprintf("Exposition : %f\n", self$estimates_params$outcome$beta1))
      cat("Effets des mediateurs  :\n"); print(self$estimates_params$outcome$beta2)
      cat("Effets Médiateurs  :\n"); print(self$estimates_params$outcome$beta3)
      cat("\n")
      
      cat("\n========================================================\n")
    },
    compute_rho01 = function(params) {
      
      s1_0 <- params$primary$sigma[1]
      s1_1 <- params$primary$sigma[2]
      
      s2_0 <- params$secondary$sigma[1]
      s2_1 <- params$secondary$sigma[2]
      
      rho00 <- params$Cor[1]
      rho11 <- params$Cor[2]
      
      # D12
      D12 <- - (s1_1^2 * s2_1^2 * rho11^2) -
        (s1_0^2 * s2_0^2 * rho00^2) +
        (s1_1^2 - s1_0^2) * (s2_1^2 - s2_0^2)
      
      # Delta
      Delta12 <- D12^2 - 4 * s1_1^2 * s1_0^2 * s2_1^2 * s2_0^2 * rho11^2 * rho00^2
      
      sqrt_Delta <- safe_sqrt(Delta12)
      
      num1 <- -D12 + sqrt_Delta
      num2 <- -D12 - sqrt_Delta
      
      denom <- 2 * s1_0^2 * s2_1^2
      
      rho_num1  <-  private$safe_sqrt(num1 / denom)
      rho_num2  <- private$safe_sqrt(num2 / denom)
      
      return(c(rho_num1, rho_num2))
    },
    safe_sqrt = function(x) {
      sqrt(pmax(x, 1e-12))
    },
    
    # ===============================
    # 4. Fonction principale (wrapper)
    # ===============================
    zeta_wrapper_direct = function(theta, skeleton, df_cov_primary, df_cov_secondary, methode='CNC'){
      
      params <- relist(theta, skeleton)
      
      # ---------- PRIMARY ----------
      if(!is.null(df_cov_primary)){
        a <- tcrossprod(df_cov_primary, t(params$primary$alpha2))  
        b <- tcrossprod(df_cov_primary, t(params$primary$alpha3)) 
      } else {
        a <- b <- 0 
      }
      
      # ---------- SECONDARY ----------
      if(!is.null(df_cov_secondary)){
        c <- tcrossprod(df_cov_secondary , t(params$secondary$alpha3)) 
        d <- tcrossprod(df_cov_secondary , t(params$secondary$alpha2))
      } else {
        c <- d <- 0 
      }
      
      V01 <- params$primary$alpha0 + a
      V11 <- params$primary$alpha1 + b
      Vbase <- params$secondary$alpha1 + c
      V02 <- Vbase + params$secondary$alpha0 + d
      
      # ---------- RHO ----------
      if(methode == "CNC"){
        
        rho01 <- private$compute_rho01(params)
        rho01 <- c(rho01, -rho01)
        rho0 <- params$primary$sigma[1] * params$secondary$sigma[1] * params$Cor[2]
        rho1 <- params$primary$sigma[2] * params$secondary$sigma[2] * params$Cor[1]
        
        rho22 <- params$primary$sigma[2] * params$secondary$sigma[1] * rho01
        
        Rho0 <- rho22 - rho0
        Rho1 <- rho1 - rho22
        
      } else {
        
        rho22 <- params$primary$sigma * params$secondary$sigma * params$Cor
        Rho0 <- Rho1 <- 0
      }
      
      # ---------- DIRECT ----------
      direct <- sapply(1:length(rho22), function(i){
        params$outcome$beta1 +
          params$outcome$beta2[2] * Vbase +
          params$outcome$beta3[1] * V01 +
          params$outcome$beta3[2] * V02 +
          params$outcome$beta3[3] * (V01 * Vbase + Rho0[i]) +
          params$outcome$beta3[4] * (V01 * V02 + rho22[i])
      })
      
      Direct <- if(is.matrix(direct)) colMeans(direct) else direct
      
      return(Direct)
    },
    zeta_wrapper_indirect = function(theta, skeleton, df_cov_primary, df_cov_secondary, methode='CNC'){
      
      params <- relist(theta, skeleton)
      
      # ---------- PRIMARY ----------
      if(!is.null(df_cov_primary)){
        a <- tcrossprod(df_cov_primary, t(params$primary$alpha2))  
        b <- tcrossprod(df_cov_primary, t(params$primary$alpha3)) 
      } else {
        a <- b <- 0 
      }
      
      # ---------- SECONDARY ----------
      if(!is.null(df_cov_secondary)){
        c <- tcrossprod(df_cov_secondary , t(params$secondary$alpha3)) 
        d <- tcrossprod(df_cov_secondary , t(params$secondary$alpha2))
      } else {
        c <- d <- 0 
      }
      
      V01 <- params$primary$alpha0 + a
      V11 <- params$primary$alpha1 + b
      Vbase <- params$secondary$alpha1 + c
      V02 <- Vbase + params$secondary$alpha0 + d
      
      # ---------- RHO ----------
      if(methode == "CNC"){
        
        rho01 <- private$compute_rho01(params)
        rho01 <- c(rho01, -rho01)
        rho0 <- params$primary$sigma[1] * params$secondary$sigma[1] * params$Cor[2]
        rho1 <- params$primary$sigma[2] * params$secondary$sigma[2] * params$Cor[1]
        
        rho22 <- params$primary$sigma[2] * params$secondary$sigma[1] * rho01
        
        Rho0 <- rho22 - rho0
        Rho1 <- rho1 - rho22
        
      } else {
        
        rho22 <- params$primary$sigma * params$secondary$sigma * params$Cor
        Rho0 <- Rho1 <- 0
      }
      
      # ---------- INDIRECT ----------
      indirect <- sapply(1:length(rho22), function(i){
        (params$outcome$beta2[1] + params$outcome$beta3[1]) * V11 +
          (params$outcome$beta3[3] + params$outcome$beta3[4]) * (V11 * V02 + Rho1[i])
      })
      
      Indirect <- if(is.matrix(indirect)) colMeans(indirect) else indirect
      
      return(Indirect)
    }
    
  )
)
