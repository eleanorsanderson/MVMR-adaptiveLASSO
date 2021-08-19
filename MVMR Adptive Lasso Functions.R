#Depends: lars, AER
library(lars)
library(AER)

### adaptive Lasso method for IV selection with two endogenous variables using (1) 10-fold cross validation (MVadap.cv) and (2) downward testing (MVadap.dt.R).

### (1) MVadap.cv
### Function: Adaptive Lasso method for instrumental variable (IV) selection and causal effects IV estimation using individual level data.
###           The penalty parameter for adaptive Lasso is chosen by 10-fold cross validation.

### Input: Y: A numeric vector of continous outcomes (n by 1 vector).
###        D: A numeric matrix of two endogenous exposures/treatments (n by 2 matrix).
###        Z: A numeric matrix of instrumental variables with each column referring to one instrument (n by p_z matrix).
###        X: An optional numeric matrix of exogenous explanatory variables including the intercept (n by p_x matrix).
###        alpha: A numeric scalar between 0 and 1 specifying the significance level for the confidence intervals for the causal effect estimate (default = 0.05).

### Output: (1) Valid instruments:Identities of the valid instrumental variables selected by the algorithm.
###         (2) Number of Valid Instruments: The number of the selected valid instrumental variables.
###         (3) Coefficients:The matrix for the post-selection 2SLS estimation results for the coefficients 
###             of the endogenous exposure/treatment variables and exogenous explanatory variables, including the point estimates and standard errors.
###         (4) Confidence interval: The confidence interval for the post-selection 2SLS estimates 
###             for the coefficients of the endogenous exposure/treatment variables with significance level specified by alpha (default = 0.05).
###         (5) p-value of Sargan: The p-value for the Sargan overidentifying test for the selected valid instruments.


MVadap.cv <- function(Y, D, Z, X, alpha = 0.05){
  
  # Define Constants
  n <- length(Y); 
  pz <- ncol(Z);
  
  # Save data
  y <- Y; d <- D; z <- Z;
  
  if (!missing(X)){
    X <- cbind(1, X);
    if(is.null(colnames(X))){colnames(X) = c("intercept", paste0('X', 1:(ncol(X)-1)))}else{colnames(X) = c("intercept", colnames(X[,-1]))};
  }else{
    X <- matrix(1, n, 1);
    colnames(X) <- "intercept";
  };
  
  Covariates = cbind(D,X);
  Exogenous = cbind(Z,X);
  
  # Variable names
  if(!is.null(colnames(D))){Dname = colnames(D)}else{Dname = c("D1", "D2")};
  colnames(Covariates)[1:2] <- Dname;
  if(!is.null(colnames(Z))){InstrumentNames = colnames(Z)}else{InstrumentNames = paste0('Z', 1:ncol(Z))};
  
  # Centralize
  Y <- qr.resid(qr(X), Y);
  D <- qr.resid(qr(X), D);
  Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
  
  ## Calculate the median-of-medians estimates
  betam1 = betam2 <- matrix(NA,pz,pz);
  for (i in 1:(pz-1)){
    comb <- c((i+1):pz);
    ivm <- list();
    for(j in 1:length(comb)){
      ivm[[j]] <- ivreg(Y ~ D + Z[,-c(i,comb[j])] - 1 | Z);
      betam1[comb[j],i] <- coef(ivm[[j]])[1];
      betam2[comb[j],i] <- coef(ivm[[j]])[2];
    }
  }
  
  betam1 <- Matrix::forceSymmetric(betam1,uplo="L"); 
  betam2 <- Matrix::forceSymmetric(betam2,uplo="L");
  med1 <- median(apply(betam1,2,function(x){median(x, na.rm = TRUE)})); 
  med2 <- median(apply(betam2,2,function(x){median(x, na.rm = TRUE)}));
  
  # Calculate the penalty weights for the adaptive Lasso
  lm_Reduced <- lm(cbind(Y, D) ~ Z - 1);
  gamma_Y <- coef(lm_Reduced)[, 1];
  gamma_D <- coef(lm_Reduced)[, -1];
  pi_median <- gamma_Y - gamma_D %*% c(med1,med2);
  
  # convert the adaptive Lasso
  QR <- qr(Z); Dhat <- qr.fitted(QR,D); DhatZ <- t(Dhat)%*%Z;  Zt <- Z - Dhat%*%solve(crossprod(Dhat))%*%DhatZ;
  
  # Function of the adaptive Lasso estimators for beta
  adpbeta <- function(Y,D,Z,pi){
    QR <- qr(Z);
    D1hat <- qr.fitted(QR,D[,1]);
    D2hat <- qr.fitted(QR,D[,2]);
    beta <- solve(rbind(cbind(t(D1hat)%*%D1hat, t(D1hat)%*%D2hat), cbind(t(D1hat)%*%D2hat, t(D2hat)%*%D2hat)))%*%rbind(t(D1hat)%*%(Y - Z%*%pi), t(D2hat)%*%(Y - Z%*%pi))
    return(beta)
  }
  
  # Function for Sargan Test
  Sargan.test <- function(res, Z, n) {
    U <- solve(crossprod(Z));
    (t(res) %*% Z %*% U %*% t(Z) %*% res) /
      (t(res) %*% res / n)
  }
  
  # 10-fold cross validation adaptive Lasso for IV selection using lars package
  penaltya = 1/abs(pi_median); # weights for the adaptive Lasso
  Ztt <- Zt %*% diag(c(1/penaltya));
  fitall = lars(Ztt,Y,intercept = FALSE,normalize=FALSE);
  
  lambdaSeq = c(fitall$lambda,seq(from=min(fitall$lambda,na.rm=TRUE),to=2*max(fitall$lambda,na.rm=TRUE),length.out = 100));
  lambdaSeq = sort(unique(lambdaSeq));
  
  K <- 10;
  Kfolds <- split(sample(1:n),rep(1:K, length = n));
  errormat <- matrix(0,K,length(lambdaSeq));
  for(i in seq(K)) {
    testSet <- Kfolds[[i]];
    fit <- lars(Ztt[-testSet, , drop=FALSE], Y[-testSet], intercept = FALSE, normalize = FALSE);
    coefs.pi <- predict.lars(fit, s = lambdaSeq, type = "coefficients", mode = "lambda")$coefficients;
    coefs.beta <- matrix(0, length(lambdaSeq), 2);
    for (ii in 1:length(lambdaSeq)){
      coefs.beta[ii,] <- adpbeta(Y[-testSet], D[-testSet,], Z[-testSet,], coefs.pi[ii,]);
    }
    
    Y.test <- Y[testSet]; D.test <- D[testSet,]; Z.test <- Z[testSet,];
    
    QR <- qr(Z.test);
    residTest <- matrix(0, nrow(Z.test), length(lambdaSeq));
    for (ii in 1:length(lambdaSeq)){residTest[,ii] <- as.numeric(Y.test) - Z.test %*% coefs.pi[ii,] - D.test %*% coefs.beta[ii,]};
    badLambda <- which(is.na(residTest[1,]) == TRUE); goodLambda <- which(is.na(residTest[1,]) == FALSE);
    errormat[i,badLambda] <- NA;
    if(length(goodLambda) == 1) {
      errormat[i,goodLambda] <- sum(qr.fitted(QR,residTest[,goodLambda])^2);
    } else {
      errormat[i,goodLambda] <- colSums(qr.fitted(QR, residTest[,goodLambda])^2); #Sargan statistics as object function for cv?
    }
  }
  cv <- colMeans(errormat);
  
  stderror <- apply(errormat,2,function(x){sd(x)/sqrt(K)});
  mincv.index <- which.min(cv);
  onestderr.rule.index <- max(which(cv <= (cv[mincv.index] + stderror[mincv.index]) & cv[mincv.index] <= cv));
  pialars.cv <- predict.lars(fitall,s = lambdaSeq[mincv.index],type="coefficients", mode = "lambda")$coefficients;
  pialars.cvse <- predict.lars(fitall,s = lambdaSeq[onestderr.rule.index],type="coefficients", mode = "lambda")$coefficients;
  
  ## selection results
  invalid.cv <- as.numeric(which(pialars.cv != 0,arr.ind = T));
  ninvalid.cv <- length(invalid.cv);
  
  invalid.cvse <- as.numeric(which(pialars.cvse != 0,arr.ind = T));
  ninvalid.cvse <- length(invalid.cvse);
  
  #Estimation -- adaptive Lasso estimator
  betalasso.cv <- adpbeta(Y, D, Z, pialars.cv);
  betalasso.cvse <- adpbeta(Y, D, Z, pialars.cvse);
  
  #Estimation -- post-selection 2SLS estimator (cv)
  if(ninvalid.cv > pz-2){ #underidentified, do OLS
    print("Less than two of the instruments are selected as valid, do OLS.")
    regressor_cv <- cbind(Covariates, z);
    betaa.cv <- qr.coef(qr(regressor_cv), y)[1:ncol(Covariates)];
    res_cv <- qr.resid(qr(regressor_cv), y);
    Psar_cv <- NA;
  } else {
    z_invalid <- matrix(z[,invalid.cv], ncol = ninvalid.cv, nrow = n);
    regressor_cv_temp <- cbind(Covariates, z_invalid);
    regressor_cv <- cbind(fitted(lm(Covariates[,1:2] ~ Exogenous)), Covariates[,-c(1,2)], z_invalid);
    iv_cv <- AER::ivreg(y ~ regressor_cv_temp - 1 | Exogenous);
    betaa.cv <- coef(iv_cv)[1:ncol(Covariates)];
    res_cv <- resid(iv_cv);
    Sargan_value_cv <- Sargan.test(res_cv,Exogenous, n);
    Psar_cv <- pchisq(Sargan_value_cv, pz - ninvalid.cv - 2, lower.tail = FALSE);
  }
  
  #Estimation -- post-selection 2SLS estimator (cvse)
  if(ninvalid.cvse > pz-2){ #underidentified, do OLS
    print("Less than two of the instruments are selected as valid, do OLS.")
    regressor_cvse <- cbind(Covariates, z);
    betaa.cvse <- qr.coef(qr(regressor_cvse), y)[1:cv];
    res_cvse <- qr.resid(qr(regressor_cvse), y);
    Psar_cvse <- NA;
  } else {
    z_invalid <- matrix(z[,invalid.cvse], ncol = ninvalid.cvse, nrow = n);
    regressor_cvse_temp <- cbind(Covariates, z_invalid);
    regressor_cvse <- cbind(fitted(lm(Covariates[,1:2] ~ Exogenous)), Covariates[,-c(1,2)], z_invalid);
    iv_cvse <- AER::ivreg(y ~ regressor_cvse_temp - 1 | Exogenous);
    betaa.cvse <- coef(iv_cvse)[1:ncol(Covariates)];
    res_cvse <- resid(iv_cvse);
    Sargan_value_cvse <- Sargan.test(res_cvse, Exogenous, n);
    Psar_cvse <- pchisq(Sargan_value_cvse, pz - ninvalid.cvse - 2, lower.tail = FALSE);
  }
  
  # Standard Errors and Confidence Intervals
  sd_cv <- sqrt(diag(
    solve(crossprod(regressor_cv)) %*%
      crossprod(res_cv * regressor_cv) %*%
      solve(crossprod(regressor_cv))
  ))[1:ncol(Covariates)];
  ci_cv <- matrix(c(
    betaa.cv[1:2] - qnorm(1-alpha/2) * sd_cv[1:2],
    betaa.cv[1:2] + qnorm(1-alpha/2) * sd_cv[1:2]
  ),2 ,2);
  
  sd_cvse <- sqrt(diag(
    solve(crossprod(regressor_cvse)) %*%
      crossprod(res_cvse * regressor_cvse) %*%
      solve(crossprod(regressor_cvse))
  ))[1:ncol(Covariates)];
  ci_cvse <- matrix(c(
    betaa.cvse[1:2] - qnorm(1-alpha/2) * sd_cvse[1:2],
    betaa.cvse[1:2] + qnorm(1-alpha/2) * sd_cvse[1:2]
  ), 2 ,2);
  
  # Results
  object <- list(
    if(ninvalid.cv > pz-2){
      invalid.cv = paste0("No instruments selected as valid.");
    }else{
      invalid.cv = InstrumentNames[invalid.cv];
    },
    invalid.cv = invalid.cv,
    ninvalid.cv = ninvalid.cv,
    if(ninvalid.cvse > pz-2){
      invalid.cvse = paste0("No instruments selected as valid.");
    }else{
      invalid.cvse = InstrumentNames[invalid.cvse];
    },
    invalid.cvse = invalid.cvse,
    ninvalid.cvse = ninvalid.cvse,
    median_estimator = c(med1,med2),
    betalasso.cv = betalasso.cv,
    betalasso.cvse = betalasso.cvse,
    betaa.cv = betaa.cv,
    betaa.cvse = betaa.cvse,
    sd_cv = sd_cv,
    sd_cvse = sd_cvse,
    ci_cv = ci_cv,
    ci_cvse = ci_cvse,
    Psar_cv = Psar_cv,
    Psar_cvse = Psar_cv,
    Covariate_Names = colnames(Covariates)
  )
  
  class(object) <- "MVadap.cv"
  
  object
  
}

# Internal Function for Printing the Results.
# This is not to be called by the users.

print.MVadap.cv <- function(object){
  cat("\nInvalid Instruments (minimum cross-validation Sargan statistics, cv):\n", object$invalid.cv, 
      "\n","\nNumber of Invalid Instruments (cv):\n", object$ninvalid.cv,"\n");
  cat("\nInvalid Instruments (one-standard-error rule,cvse):\n", object$invalid.cvse, 
      "\n","\nNumber of Invalid Instruments (cvse):\n", object$ninvalid.cvse,"\n");
  
  cat("\nCoefficients:\n");
  names(object$median_estimator) = names(object$betalasso.cv) = names(object$betalasso.cvse) = NULL;
  names(object$betaa.cv) = names(object$betaa.cvse) = names(object$sd_cv) = names(object$sd_cvse) = NULL;
  names(object$ci_cv) = names(object$ci_cvse) = NULL;
  
  length(object$median_estimator) = length(object$betalasso.cv) = length(object$betalasso.cvse) = length(object$betaa.cv);
  coef.MVadp.cv<- cbind(object$betaa.cv, object$sd_cv, object$betaa.cvse,  object$sd_cvse, 
                         object$betalasso.cv, object$betalasso.cvse, object$median_estimator);
  rownames(coef.MVadp.cv) <- object$Covariate_Names;
  colnames(coef.MVadp.cv) <- c("Post 2SLS(cv)", "SE(cv)", "Post 2SLS(cvse)", "SE(cvse)",
                               "ALasso(cv)", "ALasso(cvse)", "Median-of-medians");
  print(coef.MVadp.cv, quote = FALSE);
  
  cat(sprintf('\nConfidence Interval (cv): [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_cv[1,1], object$ci_cv[1,2], object$ci_cv[2,1], object$ci_cv[2,2]));
  cat(sprintf('Confidence Interval (cvse): [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_cvse[1,1], object$ci_cvse[1,2], object$ci_cvse[2,1], object$ci_cvse[2,2]));
  
  cat("\np-value of Sargan (cv):", object$Psar_cv, "\n");
  cat("p-value of Sargan (cvse):", object$Psar_cvse);
  
}

### (2) MVadap.dt
### Function: Adaptive Lasso method for instrumental variable (IV) selection and causal effects IV estimation using individual level data.
###           The penalty parameter for adaptive Lasso is chosen by a downward testing procedure based on the Sargan test.

### Input: Y: A numeric vector of continous outcomes (n by 1 vector).
###        D: A numeric matrix of two endogenous exposures/treatments (n by 2 matrix).
###        Z: A numeric matrix of instrumental variables with each column referring to one instrument (n by p_z matrix).
###        X: An optional numeric matrix of exogenous explanatory variables including the intercept (n by p_x matrix).
###        alpha: A numeric scalar between 0 and 1 specifying the significance level for the confidence interval for the causal effect estimate (default = 0.05).
###        tuning: A numeric scalar specifiying the threshold p-value for the Sargan test ((default = 0.1/log(n)).)

### Output: (1) Valid instruments:Identities of the valid instrumental variables selected by the algorithm.
###         (2) Number of Valid Instruments: The number of the selected valid instrumental variables.
###         (3) Coefficients:The matrix for the post-selection 2SLS estimation results for the coefficients 
###             of the endogenous exposure/treatment variables and exogenous explanatory variables, including the point estimates and standard errors.
###         (4) Confidence interval: The confidence interval for the post-selection 2SLS estimates 
###             for the coefficients of the endogenous exposure/treatment variables with significance level specified by alpha (default = 0.05).
###         (5) p-value of Sargan: The p-value for the Sargan overidentifying test for the selected valid instruments.


MVadap.dt <- function(Y, D, Z, X, alpha = 0.05, tuning = 0.1/log(length(Y))){
  
  # Define Constants
  n <- length(Y); 
  pz <- ncol(Z);
  
  # Save data
  y <- Y; d <- D; z <- Z;
  
  if (!missing(X)){
    X <- cbind(1, X);
    if(is.null(colnames(X))){colnames(X) = c("intercept", paste0('X', 1:(ncol(X)-1)))}else{colnames(X) = c("intercept", colnames(X[,-1]))};
  }else{
    X <- matrix(1, n, 1);
    colnames(X) <- "intercept";
  };
  
  Covariates = cbind(D,X);
  Exogenous = cbind(Z,X);
  
  # Variable names
  if(!is.null(colnames(D))){Dname = colnames(D)}else{Dname = c("D1", "D2")};
  colnames(Covariates)[1:2] <- Dname;
  if(!is.null(colnames(Z))){InstrumentNames = colnames(Z)}else{InstrumentNames = paste0('Z', 1:ncol(Z))};
  
  # Centralize
  Y <- qr.resid(qr(X), Y);
  D <- qr.resid(qr(X), D);
  Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
  
  ## Calculate the median-of-medians estimates
  betam1 = betam2 <- matrix(NA,pz,pz);
  for (i in 1:(pz-1)){
    comb <- c((i+1):pz);
    ivm <- list();
    for(j in 1:length(comb)){
      ivm[[j]] <- ivreg(Y ~ D + Z[,-c(i,comb[j])] - 1 | Z);
      betam1[comb[j],i] <- coef(ivm[[j]])[1];
      betam2[comb[j],i] <- coef(ivm[[j]])[2];
    }
  }
  
  betam1 <- Matrix::forceSymmetric(betam1,uplo="L"); 
  betam2 <- Matrix::forceSymmetric(betam2,uplo="L");
  med1 <- median(apply(betam1,2,function(x){median(x, na.rm = TRUE)})); 
  med2 <- median(apply(betam2,2,function(x){median(x, na.rm = TRUE)}));
  
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
  Sargan.test <- function(res, Z, n) {
    U <- solve(crossprod(Z));
    (t(res) %*% Z %*% U %*% t(Z) %*% res) /
      (t(res) %*% res / n)
  }
  
  # downward testing
  Psar <- 0;
  maxs <- pz;
  ind <- min(which (rowSums(pilars != 0) != 0));
  
  while (Psar < tuning && maxs>2){
    
    validiv <- which(pilars[ind, ] == 0);
    maxs <- length(validiv);
    
    Zv <- Z[, -validiv];
    res <- resid(AER::ivreg(Y ~ cbind(D, Zv) - 1 | Z));
    sar <- Sargan.test(res, Z, n);
    
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
      sar <- Sargan.test(res, Z, n);
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
    res_dt <- qr.resid(qr(regressor_dt), y);
    Psar_dt <- NA;
  } else {
    z_invalid <- matrix(z[,invalid.dt], ncol = ninvalid.dt, nrow = n);
    regressor_dt_temp <- cbind(Covariates, z_invalid);
    regressor_dt <- cbind(fitted(lm(Covariates[,1:2] ~ Exogenous)), Covariates[,-c(1,2)], z_invalid);
    iv_dt <- AER::ivreg(y ~ regressor_dt_temp - 1 | Exogenous);
    betaa.dt <- coef(iv_dt)[1:ncol(Covariates)];
    res_dt <- resid(iv_dt);
    Psar_dt <- Psar;
  }
  
  
  # Standard Errors and Confidence Intervals
  sd_dt <- sqrt(diag(
    solve(crossprod(regressor_dt)) %*%
      crossprod(res_dt * regressor_dt) %*%
      solve(crossprod(regressor_dt))
  ))[1:ncol(Covariates)];
  ci_dt <- matrix(c(
    betaa.dt[1:2] - qnorm(1-alpha/2) * sd_dt[1:2],
    betaa.dt[1:2] + qnorm(1-alpha/2) * sd_dt[1:2]
  ),2 ,2);
  
  # Results
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
    Covariate_Names = colnames(Covariates)
  )
  
  class(object) <- "MVadap.dt"
  
  object
  
}

# Internal Function for Printing the Results.
# This is not to be called by the users.

print.MVadap.dt <- function(object){
  cat("\nInvalid Instruments:\n", object$invalid.dt, 
      "\n","\nNumber of Invalid Instruments:\n", object$ninvalid.dt,"\n");
  
  cat("\nCoefficients:\n");
  names(object$betaa.dt) = names(object$sd_dt) = names(object$median_estimator) = names(object$ci_dt) = NULL;
 
  length(object$median_estimator) = length(object$betaa.dt);
  coef.MVadp.dt<- cbind(object$betaa.dt, object$sd_dt, object$median_estimator);
  rownames(coef.MVadp.dt) <- object$Covariate_Names;
  colnames(coef.MVadp.dt) <- c("Post 2SLS", "SE", "Median-of-medians");
  print(coef.MVadp.dt, quote = FALSE);
  
  cat(sprintf('\nConfidence Interval: [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_dt[1,1], object$ci_dt[1,2], object$ci_dt[2,1], object$ci_dt[2,2]));
  
  cat("\np-value of Sargan:", object$Psar_dt);
  
}

