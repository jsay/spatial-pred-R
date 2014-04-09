
sppred <- function(object, newdata = NULL, listw = NULL,
                   zero.policy = NULL, condset = "X", avg = "DEF",
                   legacy= TRUE, power= NULL, order= 250,
                   tol= .Machine$double.eps^(3/5), ...) {
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
    Wlg <- substr(names(object$coefficients), 1, 4)== "lag."
    B <- object$coefficients[ !Wlg]
    if (is.null(newdata)){
        nd  <- FALSE
        X   <- object$X[, !Wlg]
    } else {
        nd  <- TRUE
        frm <- formula(object$call)
        mt  <- delete.response(terms(frm, data = newdata))
        mf  <- model.frame(mt, newdata)
        X   <- model.matrix(mt, mf)
        if (any(object$aliased)) X <- X[, -which(object$aliased)]
    }
    ## WEIGHT MATRIX
    if (!nd) lsw <- eval(object$call$listw) else lsw <- listw
    ## THE PREDICTORS
    if (condset== "X") prd <- as.vector(X %*% B)
    if (condset== "XW")
        prd <- prd1(object, mod, nd, B, X, lsw)
    if (condset== "XW" && !mod %in% c("sem", "sxm") && avg == "INV")
        prd <- prd2(object, prd, mod, lsw, power, order, tol)
    if (condset== "XWe")
        prd <- prd3(object, B, X, listw, power, legacy, order, tol)
    if (condset== "XWy")
        prd <- prd4(object, B, X, listw, power, legacy, order, tol)
    if (condset== "XWc")
        prd <- prd5(object, B, X, listw, power, legacy, order, tol)
    class(prd) <- "sppred"
    prd
}

prd1 <- function(object, mod= mod, nd= nd, B= B, X= X, lsw= lsw){
    if (mod!= "sem" && nd){
        if (is.null(lsw) || !inherits(lsw, "listw"))
            stop("spatial weights list required")
    }
    if (mod %in% c("sxm", "sdm", "smc")){
        m <- ncol(X)
        K <- ifelse(colnames(object$X)[ 1] == "(Intercept)", 2, 1)
        WX <- matrix(nrow= nrow(X), ncol= m+ 1- K)
        for (k in K: m){
            wx <- lag.listw(lsw, X[, k])
            if (any(is.na(wx)))
                stop("NAs in lagged independent variable")
            WX[, k+ 1- K] <- wx
        }
        prdWX <- cbind(X, WX) %*% object$coefficients
    } else {
        prdWX <- X %*% B
    }
    as.vector(prdWX)
}

prd2 <- function(object, prd= prd, mod= mod, lsw= lsw,
                 power= power, order= order, tol= tol){
    if (power){
        W <- as(as_dgRMatrix_listw(lsw), "CsparseMatrix")
        prdWXi <- c(as(powerWeights(W, rho= object$rho, X= prd,
                                    order= order, tol= tol), "matrix"))
    } else {
        prdWXi <- c(invIrW(lsw, object$rho) %*% prd)
    }
    as.vector(prdWXi)
}
