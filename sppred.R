
sppred <- function(object, newdata = NULL, listw = NULL, yobs= object$y,
                   condset= "DEF", blup = NULL, loo = FALSE, power = NULL,
                   zero.policy = NULL, legacy = TRUE, order = 250,
                   tol= .Machine$double.eps^(3/5), ...) {
    require(spdep)
    ## USUAL VERIFICATIONS
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = spdep:::.spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (is.null(power)) power <- object$method != "eigen"
    stopifnot(is.logical(legacy)) ; stopifnot(is.logical(power))
    ## DETERMINING THE MODEL
    if (object$type== "error"){
        mod <- ifelse(object$etype== "error", "sem", "sxm")
    } else {
        mod <- switch(object$type, "lag"= "sar", "mixed"= "sdm",
                                   "sac"= "sac", "sacmixed"= "smc")
    }
    ## DATA SHAPING
    if (mod %in% c("sem", "sxm")) {lab= object$lambda ; rho= 0         }
    if (mod %in% c("sar", "sdm")) {lab= 0             ; rho= object$rho}
    if (mod %in% c("sac", "smc")) {lab= object$lambda ; rho= object$rho}
    Wlg <- substr(names(object$coefficients), 1, 4)== "lag."
    B <- object$coefficients[ !Wlg] ; Bl <- object$coefficients[ Wlg]
    if (is.null(newdata)){
        X   <- object$X[, !Wlg]
    } else {
        frm <- formula(object$call)
        mt  <- delete.response(terms(frm, data = newdata))
        mf  <- model.frame(mt, newdata)
        X   <- model.matrix(mt, mf)
        if (any(object$aliased)) X <- X[, -which(object$aliased)]
    }
    ## WEIGHT MATRIX, add an error message
    if (is.null(listw)) lsw <- eval(object$call$listw) else lsw <- listw
    ## PREDICTORS
    if (is.null(blup)){
        pt <- switch(condset, "X"= 1, "XW"= 2, "DEF"= 3, "XWy"= 4)
    } else {
        pt <- switch(blup, "LSP"= 5, "KP2"= 6, "KP3"= 7, "KPG"= 8)
    }
    prdX <- as.vector(X %*% B)
    if (pt> 1) prdWX   <- prdWX(prdX, X, Bl, mod, lsw)
    if (pt> 2 && pt!= 4) prdKP1  <- prdKP1(prdWX, rho, lsw, power, order, tol)
    if (pt> 3){
        prdWXy <- prdWX+
            rho* lag.listw(lsw, yobs)+ lab* lag.listw(lsw, yobs- prdWX)}
    if (pt==5) prdLSP <- prdLSP(prdKP1, rho, lab, lsw, yobs, loo)
    if (pt> 5 && !loo) stop("Set loo= TRUE for this blup predictor")
    if (pt==6){
        prdKP2 <- prdKP2(prdKP1, prdWXy,
                         rho, lab, lsw, yobs, power, order, tol)}
    if (pt==7){
        prdKP3 <- prdKP3(prdKP1, prdWXy,
                         rho, lab, lsw, yobs, power, order, tol)}
    if (pt==8) stop("not implemented")
    prd <- switch(pt, "1"= prdX  , "2"= prdWX , "3"= prdKP1, "4"= prdWXy,
                      "5"= prdLSP, "6"= prdKP2, "7"= prdKP3, "8"= prdKPG)
    class(prd) <- "sppred" ; as.vector(prd)
}

prdWX <- function(prdX, X= X, Bl= Bl, mod= mod, lsw= lsw){
    if (!mod %in% c("sxm", "sdm", "smc")){
        prdWX <- prdX } else {
            K <- ifelse(colnames(X)[ 1] == "(Intercept)", 2, 1)
            m <- ncol(X) ; WX <- matrix(nrow= length(prdX), ncol= m+ 1- K)
            for (k in K: m){
                WX[, k+ 1- K] <- lag.listw(lsw, X[, k])
            }
            prdWX <- prdX+ (WX %*% Bl)
        } 
    prdWX
}

prdKP1 <- function(prdWX, rho= rho, lsw= lsw,
                   power= power, order= order, tol= tol){
    if (power){
        W <- as(as_dgRMatrix_listw(lsw), "CsparseMatrix")
        prdKP1 <- c(as(powerWeights(W, rho= rho, X= as.matrix(prdWX),
                                    order= order, tol= tol), "matrix"))
    } else {
        prdKP1 <- c(invIrW(lsw, rho) %*% prdWX)
    }
    prdKP1
}

prdLSP <- function(prdKP1, rho= rho, lab= lab,
                   lsw= lsw, yobs= yobs, loo= loo){
    ZL <- diag(length(prdKP1))- (lab* listw2mat(lsw))
    ZR <- diag(length(prdKP1))- (rho* listw2mat(lsw))
    Z  <- ZL %*% ZR ; P22 <- t(Z) %*% Z
    if (loo){
        prdLSP <- matrix(NA, ncol= 1, nrow= length(prdKP1))
        for (i in 1: length(prdKP1)){
            prdLSP[ i] <- prdKP1[ i]-
                (P22[i, -i] %*% (yobs[ -i]- prdKP1[ -i])/ P22[i, i])
        }
    } else {
        P11 <- P22
        prdLSP <- prdKP1+ ((solve(P22) %*% P11 %*% (yobs- prdKP1)))
    }
    prdLSP
}

prdKP2 <- function(prdKP1, prdWXy= prdWXy, rho= rho, lab= lab, lsw= lsw,
                   yobs= yobs, power= power, order= order, tol= tol){
    if (power){
        W <- as(as_dgRMatrix_listw(lsw), "CsparseMatrix")
        GL <- as(powerWeights(W, rho= lab, order= order, tol= tol,
                              X= diag(length(prdWXy))), "matrix")
        GR <- as(powerWeights(W, rho= rho, order= order, tol= tol,
                              X= diag(length(prdWXy))), "matrix")
    } else {
        GL <- invIrW(lsw, rho) ; GR <- invIrW(lsw, lab)
    }
    sum.u <- GL %*% t(GL) ; sum.y <- GR %*% sum.u %*% t(GR)
    prdKP2 <- matrix(NA, ncol= 1, nrow= length(prdWXy))
    for (i in 1: length(prdKP2)){
        WM <- listw2mat(lsw)[i, ]
        rg <- (sum.u[i, ] %*% GR %*% WM)/ (WM %*% sum.y %*% WM)
        prdKP2[ i] <- prdWXy[ i]+ (rg %*% WM %*% (yobs- prdKP1))
    }
    prdKP2
}

prdKP3 <- function(prdKP1, prdWXy= prdWXy, rho= rho, lab= lab, lsw= lsw,
                   yobs= yobs, power= power, order= order, tol= tol){
    if (power){
        W <- as(as_dgRMatrix_listw(lsw), "CsparseMatrix")
        GL <- as(powerWeights(W, rho= lab, order= order, tol= tol,
                              X= diag(length(prdWXy))), "matrix")
        GR <- as(powerWeights(W, rho= rho, order= order, tol= tol,
                              X= diag(length(prdWXy))), "matrix")
    } else {
        GL <- invIrW(lsw, lab) ; GR <- invIrW(lsw, rho)
    }
    sum.u <- GL %*% t(GL) ; sum.y <- GR %*% sum.u %*% t(GR)
    prdKP3 <- matrix(NA, ncol= 1, nrow= length(prdWXy))    
    for (i in 1: length(prdKP3)){
        rg <- sum.u[i, ] %*% GR[, -i] %*% solve(sum.y[-i, -i])
        prdKP3[ i] <- prdWXy[ i]+ (rg %*% (yobs[i ]- prdKP1[-i ]))
    }
    prdKP3
}
