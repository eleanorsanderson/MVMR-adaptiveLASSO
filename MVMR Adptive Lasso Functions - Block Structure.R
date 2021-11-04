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

### Output: (1) Invalid instruments:           Identities of the invalid instrumental variables selected by the algorithm.
###         (2) Number of Invalid Instruments: The number of the selected invalid instrumental variables.
###         (3) Coefficients:                  The matrix for the post-selection 2SLS and two-step GMM estimation results for the coefficients 
###                                            of the endogenous exposure/treatment variables and exogenous explanatory variables, 
###                                            including the point estimates and standard errors.
###         (4) Confidence interval:           The confidence interval for the post-selection 2SLS estimates 
###                                            for the coefficients of the endogenous exposure/treatment variables with significance level 
###                                            specified by alpha (default = 0.05).
###         (5) p-value of the Hansen J-test:  The p-value for the two-step Hansen J-test for the selected valid instruments.


MVadap.dtblock <- function(Y, D, Z, index1 = c(1:ncol(Z)), index2 = c(1:ncol(Z)), X, alpha = 0.05, tuning = 0.1/log(length(Y))){
  
  # Define Constants
  n <- length(Y); 
  pz <- ncol(Z);
  
  # Save data
  y <- Y; d <- D; z <- Z;
  
  if (!missing(X)){
    if(is.null(colnames(X))){
      colnames(X) = c(paste0('X', 1:(ncol(X))))
    }
    X <- cbind(1, X);
    colnames(X) = c('intercept', colnames(X)[-1])
  }else{
    X <- matrix(1, n, 1);
    colnames(X) <- "intercept";
  };
  
  Covariates = cbind(D,X);
  Exogenous = cbind(Z,X);
  
  # Variable names
  if(!is.null(colnames(D))){Dname = colnames(D)}else{Dname = c("D1", "D2")};
  colnames(Covariates)[1:2] <- Dname;
  if(!is.null(colnames(Z))){InstrumentNames = colnames(Z)}else{InstrumentNames = c(1:ncol(Z))};
  
  # Centralize
  Y <- qr.resid(qr(X), Y);
  D <- qr.resid(qr(X), D);
  Z <- scale(qr.resid(qr(X), Z), center = FALSE, scale = TRUE);
  
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
  
  # Function for Hansen J-test
  HansenJTest <- function(Y,X,Z){
    res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z));
    
    Weight_SecondStep <- crossprod(res_FirstStep * Z);
    
    Coef_SecondStep <- solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y;
    
    res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep);
    
    sd_SecondStep <- sqrt(diag(solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
      solve(
        t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
      ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
    )));
    
    HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
      t(Z) %*% res_SecondStep;
    
    list(HansenJ_Stat, Coef_SecondStep,sd_SecondStep)
  }
  
  # downward testing
  PHJ <- 0;
  maxs <- pz;
  ind <- 1;
  
  while (PHJ < tuning && maxs>2){
    
    validiv <- which(pilars[ind, ] == 0);
    maxs <- length(validiv);
    
    Zv <- Z[, -validiv];
    HJtest <- HansenJTest(Y, cbind(D, Zv), Z)[[1]];
    
    PHJ <- pchisq(HJtest, maxs - 2, lower.tail = FALSE);
    
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
  
  #Estimation -- post-selection 2SLS and two-step GMM estimators
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
    
    Coefficients_dt_GMM <- (HansenJTest(y, regressor_dt_temp, Exogenous)[[2]])[1:ncol(Covariates)];
    
    PHJ_dt <- PHJ;
  }
  
  
  # Standard Errors and Confidence Intervals
  sd_dt <- sqrt(diag(
    mean(res_dt^2) * solve(crossprod(regressor_dt))
  ))[1:ncol(Covariates)];
  ci_dt <- matrix(c(
    betaa.dt[1:2] - qnorm(1-alpha/2) * sd_dt[1:2],
    betaa.dt[1:2] + qnorm(1-alpha/2) * sd_dt[1:2]
  ),2 ,2);
  # GMM
  sd_GMM_dt <- (HansenJTest(y, regressor_dt_temp, Exogenous)[[3]])[1:ncol(Covariates)];
  ci_GMM_dt <- matrix(c(
    Coefficients_dt_GMM[1:2] - qnorm(1-alpha/2) * sd_GMM_dt[1:2],
    Coefficients_dt_GMM[1:2] + qnorm(1-alpha/2) * sd_GMM_dt[1:2]
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
    Coefficients_dt_GMM = Coefficients_dt_GMM,
    sd_GMM_dt = sd_GMM_dt,
    ci_GMM_dt = ci_GMM_dt,
    PHJ_dt = PHJ_dt,
    Covariate_Names = colnames(Covariates)
  )
  
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
  names(object$Coefficients_dt_GMM) = names(object$sd_GMM_dt) = names(object$ci_GMM_dt) = NULL;
 
  length(object$median_estimator) = length(object$betaa.dt);
  coef.MVadp.dtblock<- cbind(object$betaa.dt, object$sd_dt, object$Coefficients_dt_GMM, object$sd_GMM_dt, object$median_estimator);
  rownames(coef.MVadp.dtblock) <- object$Covariate_Names;
  colnames(coef.MVadp.dtblock) <- c("2SLS", "2SLS.SE", "GMM", "GMM.SE", "Median-of-medians");
  print(coef.MVadp.dtblock, quote = FALSE);
  
  cat(sprintf('\n2SLS Confidence Interval: [%.4f,%.4f], [%.4f,%.4f]', object$ci_dt[1,1], object$ci_dt[1,2], object$ci_dt[2,1], object$ci_dt[2,2]));
  cat(sprintf('\nGMM Confidence Interval: [%.4f,%.4f], [%.4f,%.4f]\n', object$ci_GMM_dt[1,1], object$ci_GMM_dt[1,2], object$ci_GMM_dt[2,1], object$ci_GMM_dt[2,2]));
  
  cat("\np-value of the Hansen J-test:", object$PHJ_dt);
  
}

