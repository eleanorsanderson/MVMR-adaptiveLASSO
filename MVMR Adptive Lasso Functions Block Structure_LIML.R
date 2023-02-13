#Depends: lars, AER
library(lars)
library(AER)

### MVadap.dtblock
### Function: Adaptive Lasso method for instrumental variable (IV) selection and causal effects IV estimation using individual-level data.
###           The penalty parameter for adaptive Lasso is chosen by a downward testing procedure based on the two-step Hansen J-test.
###           The median-of-medians estimators are estimated by the block structure.

### Input: Y:      A numeric vector of continous outcomes (n by 1 vector).
###        D:      A numeric matrix of two endogenous exposures/treatments (n by 2 matrix).
###        Z:      A numeric matrix of instrumental variables, with each column referring to one instrument (n by p_z matrix).
###        index1: A numeric vector of the column identities of instruments for the first endogenous exposures. Default = c(1:ncol(Z)).
###        index2: A numeric vector of the column identities of instruments for the second endogenous exposures. Default = c(1:ncol(Z)).
###        X:      An optional numeric matrix of exogenous explanatory variables including the intercept (n by p_x matrix).
###        alpha:  A numeric scalar between 0 and 1 specifying the significance level for the confidence interval for the causal effect estimate (default = 0.05).
###        tuning: A numeric scalar specifiying the threshold p-value for the Sargan test ((default = 0.1/log(n)).)
###        LIML: Logical. If LIML = TRUE, the downward testing procedure is based on LIML estimation (default = TRUE).

### Output: (1) Valid instruments:Identities of the valid instrumental variables selected by the algorithm.
###         (2) Number of Valid Instruments: The number of the selected valid instrumental variables.
###         (3) Coefficients:The matrix for the post-selection 2SLS estimation results for the coefficients 
###             of the endogenous exposure/treatment variables and exogenous explanatory variables, including the point estimates and standard errors.
###         (4) Confidence interval: The confidence interval for the post-selection 2SLS estimates 
###             for the coefficients of the endogenous exposure/treatment variables with significance level specified by alpha (default = 0.05).
###             If LIML = TRUE, the LIML estimates and their stand errors are also reported.
###         (5) p-value of Sargan: The p-value for the Sargan overidentifying test for the selected valid instruments.


MVadap.dtblock <- function(Y, D, Z, index1 = c(1:ncol(Z)), index2 = c(1:ncol(Z)), X, 
                           alpha = 0.05, tuning = 0.1/log(length(Y)), LIML = TRUE){
  
  # Define Constants
  n <- length(Y); 
  pz <- ncol(Z);
  
  if (!missing(X)){
    if(is.null(colnames(X))){
      colnames(X) = c(paste0('X', 1:(ncol(X))))
    }
    X <- cbind(1, X);
    colnames(X) = c('intercept', colnames(X)[-1]);
  }else{
    X <- matrix(1, n, 1);
    colnames(X) <- "intercept";
  };
  
  # Variable names
  if(!is.null(colnames(D))){colnames(D) <- Dname <- colnames(D)}else{colnames(D) <- Dname <- c("D1", "D2")};
  if(!is.null(colnames(Z))){colnames(Z) <- InstrumentNames <- colnames(Z)}else{colnames(Z) <- InstrumentNames <- paste0("Z", 1:ncol(Z))};
  
  Covariates = cbind(D,X);
  Exogenous = cbind(Z,X);
  
  # Save data
  y <- Y; d <- D; z <- Z;
  
  # Centralize
  Y <- qr.resid(qr(X), Y);
  D <- qr.resid(qr(X), D);
  Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
  U <- solve(crossprod(Z));
  
  Z1 <- Z[, index1]; pz1 <- ncol(Z1);
  Z2 <- Z[, index2]; pz2 <- ncol(Z2);
  
  ## Calculate the median-of-medians estimates
  betam1 = betam2 <- matrix(NA, pz1, pz2);
  for (i in index1){
    for (j in index2){
      if (i != j){
        ivm <- ivreg(Y ~ D + Z[, -c(i,j)] - 1 | Z);
        betam1[which(index1 == i), which(index2 == j)] <- coef(ivm)[1];
        betam2[which(index1 == i), which(index2 == j)] <- coef(ivm)[2];
      }
    }
  }
  
  med1 <- median(apply(t(betam1),2,function(x){median(x, na.rm = TRUE)})); 
  med2 <- median(apply(betam2,   2,function(x){median(x, na.rm = TRUE)}));
  
  # Calculate the penalty weights for the adaptive Lasso
  lm_Reduced <- lm(cbind(Y, D) ~ Z - 1);
  gamma_Y <- coef(lm_Reduced)[, 1];
  gamma_D <- coef(lm_Reduced)[, -1];
  pi_median <- gamma_Y - gamma_D %*% c(med1,med2);
  
  # convert the adaptive Lasso
  QR <- qr(Z); Dhat <- qr.fitted(QR,D); DhatZ <- t(Dhat)%*%Z;  Zt <- Z - Dhat%*%solve(crossprod(Dhat))%*%DhatZ;
  
  # downward testing procedure
  penaltya = 1/abs(pi_median); # weights for the adaptive Lasso
  Ztt <- Zt %*% diag(c(1/penaltya));
  fitall = lars(Ztt,Y,intercept = FALSE,normalize=FALSE);
  pilars <- predict.lars(fitall, type = "coefficients", mode = "step")$coefficients;
  
  # Function for Sargan Test
  Sargan.test <- function(res, Z, U, n) {
    (t(res) %*% Z %*% U %*% t(Z) %*% res) /
      (t(res) %*% res / n)
  };
  
  # downward testing
  Psar <- 0;
  maxs <- pz;
  ind <- min(which (rowSums(pilars != 0) != 0));
  
  while (Psar < tuning && maxs>2){
    
    validiv <- which(pilars[ind, ] == 0);
    maxs <- length(validiv);
    Zv <- Z[, -validiv];
    
    if (LIML){
      Zv <- as.matrix(Z[, -validiv], nrow = n);
      colnames(Zv) <- colnames(Z)[-validiv];
      LIML_fit <- kclass.fit(x = cbind(D, Zv), y = Y, z = Z, model.type="LIML");
      res <- LIML_fit$res;
    }else{
      #res <- resid(AER::ivreg(Y ~ cbind(D, Zv) - 1 | Z));
      res <- resid(AER::ivreg(Y ~ D + Zv - 1 | Z));
    };
    
    sar <- Sargan.test(res, Z, U, n);
    Psar <- pchisq(sar, maxs - 2, lower.tail = FALSE);
    
    ind <- ind + 1;
  }
  
  alphaSuppSize = apply(pilars,1,function(x){sum(x == 0)}); # for each Lasso step, how many valid instruments
  
  if (maxs == 2){ # none valid
    print("Less than two of the instruments are selected as valid, do OLS.")
    validiv <- NULL;
    Psar <- NA;
    maxs <- 0;
  }else if (length(which(alphaSuppSize == maxs)) > 1){ #tie
    rowind <- which(alphaSuppSize == maxs);
    rowind <- rowind[-which(rowind == ind - 1)];
    Psarmax <- rep(NA, length(rowind));
    for (i in 1:length(rowind)){
      validmax <- which(pilars[rowind[i], ] == 0);
      Zv <- Z[, -validmax];
      res <- resid(AER::ivreg(Y ~ cbind(D, Zv) - 1 | Z));
      sar <- Sargan.test(res, Z, U, n);
      Psarmax[i] <- pchisq(sar, maxs - 2, lower.tail = FALSE);
    }
    Psarmax <- c(Psarmax, Psar);
    if (Psar != max(Psarmax)){
      validind <- rowind[which(Psarmax == max(Psarmax))];
      validiv <- which(pilars[validind, ] == 0);
      Psar <- max(Psarmax);
    }
  }
  
  # Selection results
  if(is.null(validiv)){invalid.dt <- seq(pz)}else{invalid.dt = as.numeric(seq(pz)[-validiv])};
  ninvalid.dt <- length(invalid.dt);
  
  #Estimation -- post-selection 2SLS estimator
  if(ninvalid.dt > pz-2){ #underidentified, do OLS
    print("Less than two of the instruments are selected as valid, do OLS.")
    regressor_dt <- cbind(Covariates, z);
    betaa.dt <- qr.coef(qr(regressor_dt), y)[1:ncol(Covariates)];
    res_dt <- c(qr.resid(qr(regressor_dt), y));
    Psar_dt <- NA;
    
    sd_dt <- sqrt(diag(
      mean(res_dt^2) * solve(crossprod(regressor_dt))
    ))[1:ncol(Covariates)];
    ci_dt <- matrix(c(
      betaa.dt[1:2] - qnorm(1-alpha/2) * sd_dt[1:2],
      betaa.dt[1:2] + qnorm(1-alpha/2) * sd_dt[1:2]
    ),2 ,2);
    
    if (LIML){
      betaa_LIML = NA;
      sd_LIML = NA;
      ci_LIML = NA;
      Psar_LIML = NA;
    }
    
  } else {
    z_invalid <- matrix(z[,invalid.dt], ncol = ninvalid.dt, nrow = n);
    colnames(z_invalid) <- colnames(z)[invalid.dt];
    regressor_dt_temp <- cbind(Covariates, z_invalid);
    
    iv_dt <- AER::ivreg(y ~ regressor_dt_temp - 1 | Exogenous);
    betaa.dt <- coef(iv_dt)[1:ncol(Covariates)];
    
    sd_dt <- (summary(iv_dt)$coefficients[, "Std. Error"])[1:ncol(Covariates)];
    ci_dt <- matrix(c(
      betaa.dt[1:2] - qnorm(1-alpha/2) * sd_dt[1:2],
      betaa.dt[1:2] + qnorm(1-alpha/2) * sd_dt[1:2]
    ),2 ,2);
    
    res_dt <- resid(iv_dt);
    sar_dt <- Sargan.test(res_dt, Z, U, n);
    Psar_dt <- pchisq(sar_dt, (pz - ninvalid.dt - 2), lower.tail = FALSE);
    
    #Psar_dt <- (summary(iv_dt, diagnostics = TRUE)$diagnostics)["Sargan", "p-value"];
    
    if(LIML){
      iv_LIML <- kclass.fit(x = regressor_dt_temp, y = y, z = Exogenous, model.type="LIML");
      betaa_LIML <- (iv_LIML$coefficients)[1:ncol(Covariates)];
      
      sd_LIML <- iv_LIML$se[1:ncol(Covariates)];
      ci_LIML <- matrix(c(
        betaa_LIML[1:2] - qnorm(1-alpha/2) * sd_LIML[1:2],
        betaa_LIML[1:2] + qnorm(1-alpha/2) * sd_LIML[1:2]
      ),2 ,2);
      
      res_LIML <- iv_LIML$res;
      sar_LIML <- Sargan.test(res_LIML, Z, U, n);
      Psar_LIML <- pchisq(sar_LIML, (pz - ninvalid.dt - 2), lower.tail = FALSE);
      
      #Psar_LIML <- Psar;
    }
    
  }
  
  # Results
  if (LIML){
    object <- list(
      if(ninvalid.dt > pz-2){
        invalid.dt = paste0("No instruments selected as valid.");
      }else{
        invalid.dt = InstrumentNames[invalid.dt];
      },
      invalid.dt = invalid.dt,
      ninvalid.dt = ninvalid.dt,
      median_estimator = c(med1,med2),
      betaa.dt = betaa.dt,
      sd_dt = sd_dt,
      ci_dt = ci_dt,
      Psar_dt = Psar_dt,
      betaa_LIML = betaa_LIML,
      sd_LIML = sd_LIML,
      ci_LIML = ci_LIML,
      Psar_LIML = Psar_LIML,
      Covariate_Names = colnames(Covariates),
      LIML = TRUE,
      nIV = pz
    )
  }else{
    object <- list(
      if(ninvalid.dt > pz-2){
        invalid.dt = paste0("No instruments selected as valid.");
      }else{
        invalid.dt = InstrumentNames[invalid.dt];
      },
      invalid.dt = invalid.dt,
      ninvalid.dt = ninvalid.dt,
      median_estimator = c(med1,med2),
      betaa.dt = betaa.dt,
      sd_dt = sd_dt,
      ci_dt = ci_dt,
      Psar_dt = Psar_dt,
      Covariate_Names = colnames(Covariates),
      LIML = FALSE,
      nIV = pz
    )
  }
  
  class(object) <- "MVadap.dtblock"
  
  object
  
}

# Internal Function for Printing the Results.
# This is not to be called by the users.
print.MVadap.dtblock <- function(object){
  
  cat("\nInvalid Instruments:\n", object$invalid.dt, 
      "\n","\nNumber of Invalid Instruments:\n", object$ninvalid.dt,"\n");
  
  cat("\nCoefficients:\n");
  names(object$betaa.dt) = names(object$sd_dt) = names(object$median_estimator) = names(object$ci_dt) = NULL;
  
  if(object$LIML){names(object$betaa_LIML) = names(object$sd_LIML) = names(object$ci_LIML) = NULL;}
  
  length(object$median_estimator) = length(object$betaa.dt);
  
  if (object$ninvalid.dt > object$nIV-2){
    coef.MVadp.dt<- cbind(object$betaa.dt, object$sd_dt, object$median_estimator);
    colnames(coef.MVadp.dt) <- c("Post OLS", "OLS-SE", "Median-of-medians");
    rownames(coef.MVadp.dt) <- object$Covariate_Names;
    print(coef.MVadp.dt, quote = FALSE);
    cat(sprintf('\nConfidence Interval-OLS: [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_dt[1,1], object$ci_dt[1,2], object$ci_dt[2,1], object$ci_dt[2,2]));
  }else{
    if (object$LIML){
      coef.MVadp.dt<- cbind(object$betaa.dt, object$sd_dt, object$betaa_LIML, object$sd_LIML, object$median_estimator);
      colnames(coef.MVadp.dt) <- c("Post 2SLS", "2SLS-SE", "Post LIML", "LIML-SE", "Median-of-medians");
    }else{
      coef.MVadp.dt<- cbind(object$betaa.dt, object$sd_dt, object$median_estimator);
      colnames(coef.MVadp.dt) <- c("Post 2SLS", "2SLS-SE", "Median-of-medians");
    };
    rownames(coef.MVadp.dt) <- object$Covariate_Names;
    print(coef.MVadp.dt, quote = FALSE);
    cat(sprintf('\nConfidence Interval-2SLS: [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_dt[1,1], object$ci_dt[1,2], object$ci_dt[2,1], object$ci_dt[2,2]));
    if (object$LIML){
      cat(sprintf('\nConfidence Interval-LIML: [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_LIML[1,1], object$ci_LIML[1,2], object$ci_LIML[2,1], object$ci_LIML[2,2]));
    };
    
    cat("\np-value of Sargan-2SLS:", object$Psar_dt);
    if (object$LIML) {cat("\np-value of Sargan-LIML:", object$Psar_LIML)}
  }
  
}

### NOT to be called by users
### Internal functions to calculate the LIML estimators and variances based on the rkclass package from
### https://github.com/potterzot/rkclass

### Input:
### x: a matrix of ALL the regressors, including the intercept, endogenous and exogenous variables.
### y: a vector of the outcome variable
### z: a matrix of ALL the exogenous variables (exogenous regressors + instruments), including the intercept
### k: default k = NULL, do LIML. If specify k = 1, do 2SLS

### Output:
### (1) Point estimates
### (2) SE estimates
### (3) k value

kclass.fit <- function(x, y, z=NULL, k=NULL,
                       model.type=c("TSLS", "LIML", "FULLER", "KCLASS", "OLS"),
                       eig=c("eigen","geigen"), na.action, ...) {
  
  # ensure as matrices
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(y)) y <- as.matrix(y)
  if(!is.matrix(z) & !is.null(z)) z <- as.matrix(z)
  
  # sizes of each type of regressor
  N = nrow(y)
  
  #Set weights stand in / sqrt because we multiply before the QR decomp.
  w = rep(1,N)
  
  #Offset
  offset = rep(0,N)
  
  #Get names of endogenous, exogenous, and instrument variables
  regressors = colnames(x)
  instruments = colnames(z)
  endogenous <- regressors[!(regressors %in% instruments)] #instrumented variables
  excluded_instruments = instruments[!(instruments %in% regressors)]
  exogenous = regressors[regressors %in% instruments] #non-instrumented, also called included instruments.
  #if no z, then all x are endogenous and there are no exogenous regressors
  if(is.null(z)) {
    n.ins = 0
    n.exc = 0
  }else {
    n.ins = ncol(z) # included and excluded instruments
    n.exc = length(excluded_instruments) #excluded instruments
  }
  n.reg = ncol(x) #endogenous and exogenous regressors
  n.end = length(endogenous) #endogenous/instrumented regressors
  n.exo = n.reg-n.end #exogenous X variables / included instruments
  n.inc = n.exo #included instruments
  
  #Main matrix
  A = cbind(x[,exogenous, drop=FALSE],
            z[,excluded_instruments, drop=FALSE],
            x[,endogenous, drop=FALSE],
            y)
  
  QA = qr(A * w) #multiply by the weights before doing QR
  
  if(!all(QA$pivot==1:ncol(A))) {
    warning("Rank issue.")
    max_piv <- which.max(QA$pivot)
    out_piv <- QA$pivot[(max_piv+1):length(QA$pivot)]
    n_inc <- 1:n.inc
    n_ins <- n.inc + 1:n.ins
    n_end  <- (n.exo + n.ins) + (1:n.end)
    n_y  <- n.exo + n.ins + n.end + 1
    if(any(out_piv %in% n_end)) n_end <- n_end -1
  }
  
  #Calculate model matrix from which all others are derived
  RA = qr.R(QA) [, order(QA$pivot)]
  
  #Define 'R' matrices
  #Following Belsley 1974.
  n1 = NULL
  n2 = NULL
  n3 = n.ins + 1:n.end #endogenous variables (X1)
  n4 = n.ins + n.end + 1 #regressand/response (Y)
  if(n.exo > 0) {
    n1 = 1:n.inc #included instruments / exogenous variables (X2)
    R11 = RA[n1, n1, drop=FALSE] # X2
    R13 = RA[n1, n3, drop=FALSE] # X2X1
    R14 = RA[n1, n4, drop=FALSE] # X2Y
  }
  if(n.exc > 0) {
    n2 = n.inc + 1:n.exc #excluded instruments (Z)
    R23 = RA[n2, n3, drop=FALSE] #ZX1
    R24 = RA[n2, n4, drop=FALSE] #ZY
  }
  R33 = RA[n3, n3, drop=FALSE] #X1
  R34 = RA[n3, n4, drop=FALSE] #X1Y
  
  # Determin k parameter
  # res = cbind(X1,Y) - as.matrix(X2) %*% coef(Mx)
  # coef(Mx) = solve(X2X2) %*% t(X2) %*% X1Y
  # crossprod(X1Y, X1Y) = R34, R34
  # crossprod(R11)
  # crossprod(res) = crossprod(R34) - crossprod(R11) * coef(Mxc)
  if(is.null(k)) k = .set_k(model.type, x[,endogenous], matrix(x[, exogenous], nrow = N), ### !!!add matrix for x[, exogenous]
                            z, y, alpha, N, n.ins, eig="eigen") ### !!! Here might have some problem. Original code: z[, excluded_instruments]
  else model.type = "KCLASS"
  
  #R-based matrices to solve for coefficients
  if(n.exc > 0 & n.inc > 0) { #we have both X2 and Z
    M = rbind(
      cbind(crossprod(R13) + crossprod(R23) + (1-k)*crossprod(R33), crossprod(R13,R11)),
      cbind(crossprod(R11,R13), crossprod(R11)))
    d = rbind(crossprod(R13,R14) + crossprod(R23,R24) + (1-k)*crossprod(R33,R34), crossprod(R11,R14))
    
  }else if(n.exc>0) { #we have Z and no X2 ### !!!! Here might have some problem. Original code: n.inc>0
    M = crossprod(R23) + (1-k)*crossprod(R33)
    d = crossprod(R23,R24) + (1-k)*crossprod(R33,R34)
  }else { #We have no X2 and no Z
    M = crossprod(R33)
    d = crossprod(R33,R34)
  }
  
  #Calculate the coefficients
  Minv = qr.solve(qr(M))
  coef =  Minv %*% d
  rownames(coef) = c(endogenous, exogenous)
  colnames(Minv) = rownames(Minv)
  
  #Calculate the variances: using the formula from the ivmodel package
  
  W = x ## x includes ALL the regressors, endogenous + exogenous variables
  ZXQR = qr(z) ## z includes ALL the exogenous variables (exogenous regressors + instruments)
  degF = N - ncol(x) ## degrees of freedom, N is the sample size
  kVarPointEst = matrix(0,length(k), ncol(x))
  inverseMat = solve(t(W) %*% W - k * t(W) %*% qr.resid(ZXQR,W))
  ## point estimates
  kPointEst = as.numeric(inverseMat %*% ( t(W) %*% Y - k *t(W) %*% qr.resid(ZXQR,Y)))
  ## variance estimates
  kVarPointEst = 1/degF * sum( (Y - W %*% kPointEst)^2) * diag(inverseMat)
  ## se
  sePointEst = sqrt(kVarPointEst)
  
  #Create the returned object
  est = list()
  est$coefficients = coef
  est$se = sePointEst
  est$k = k
  est$res = Y - W %*% kPointEst ## add residuals to output
  
  class(est) = "kclass"
  est
}

### internal function to calculate the k value for LIML
.set_k <- function(model.type, endogenous = NULL, exogenous = NULL, instruments = NULL, y = NULL,
                   alpha=NULL, N=NULL, L=NULL,
                   eig="eigen") {
  if(model.type=="OLS") {
    #print("model.type is set to OLS. Setting k to 0.")
    k = 0
  }
  else if(model.type=="TSLS") {
    #print("model.type is set to TSLS. Setting k to 1.")
    k = 1
  }
  else if(model.type=="LIML" | model.type=="FULLER") {
    if(missing(endogenous) | missing(exogenous) | missing(instruments) | missing(y))
      stop("To determine k for a LIML or FULLER model we must have endogenous and exogenous
           regressors, instruments, and the dependent variable.")
    xy = cbind(endogenous, y)
    ymxy = crossprod(lm.fit(exogenous, xy)$residuals)
    ymzy = crossprod(lm.fit(instruments, xy)$residuals)
    
    if(eig=="eigen") k = min(eigen(ymxy %*% solve(ymzy))$values)
    else k = min(geigen(ymxy, ymzy, symmetric=TRUE, only.values=TRUE)$values)
    
    if(model.type=="FULLER") {
      stopifnot(!is.null(alpha), !is.null(N), !is.null(L))
      k = k-alpha/(N-L)
    }
  }
  #Otherwise, just use the k provided
  k
}

