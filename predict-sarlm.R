function (object, newdata = NULL, listw = NULL, zero.policy = NULL, 
    legacy = TRUE, power = NULL, order = 250, tol = .Machine$double.eps^(3/5), 
    ...) 
{
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (object$type == "sac") 
        stop("no predict method for sac")
    if (is.null(power)) 
        power <- object$method != "eigen"
    stopifnot(is.logical(legacy))
    stopifnot(is.logical(power))
    if (is.null(newdata)) {
        res <- fitted.values(object)
        X <- object$X
        B <- object$coefficients
        y <- object$y
        tarX <- object$tarX
        tary <- object$tary
        if (object$type == "error") {
            attr(res, "trend") <- as.vector(X %*% B)
            attr(res, "signal") <- as.vector(-1 * (tary - y) - 
                -1 * (tarX - X) %*% B)
        }
        else {
            attr(res, "trend") <- as.vector(X %*% B)
            attr(res, "signal") <- as.vector(-1 * (tary - y))
        }
    }
    else {
        if (object$type == "error") {
            if (object$etype == "error") {
                B <- object$coefficients
                frm <- formula(object$call)
                mt <- delete.response(terms(frm, data = newdata))
                mf <- model.frame(mt, newdata)
                X <- model.matrix(mt, mf)
                if (any(object$aliased)) 
                  X <- X[, -which(object$aliased)]
                trend <- X %*% B
                signal <- rep(0, length(trend))
                res <- trend + signal
                attr(res, "trend") <- trend
                attr(res, "signal") <- signal
            }
            else if (object$etype == "emixed") {
                if (is.null(listw) || !inherits(listw, "listw")) 
                  stop("spatial weights list required")
                if (nrow(newdata) != length(listw$neighbours)) 
                  stop("mismatch between newdata and spatial weights")
                B <- object$coefficients
                frm <- formula(object$call)
                mt <- delete.response(terms(frm, data = newdata))
                mf <- model.frame(mt, newdata)
                X <- model.matrix(mt, mf)
                K <- ifelse(colnames(X)[1] == "(Intercept)", 
                  2, 1)
                m <- ncol(X)
                if (m > 1) {
                  WX <- matrix(nrow = nrow(X), ncol = (m - (K - 
                    1)))
                  for (k in K:m) {
                    wx <- lag.listw(listw, X[, k], zero.policy = zero.policy)
                    if (any(is.na(wx))) 
                      stop("NAs in lagged independent variable")
                    WX[, (k - (K - 1))] <- wx
                  }
                }
                if (K == 2) {
                  if (!(listw$style == "W")) {
                    intercept <- as.double(rep(1, nrow(X)))
                    wx <- lag.listw(listw, intercept, zero.policy = zero.policy)
                    if (m > 1) {
                      WX <- cbind(wx, WX)
                    }
                    else {
                      WX <- matrix(wx, nrow = nrow(X), ncol = 1)
                    }
                  }
                }
                X <- cbind(X, WX)
                if (any(object$aliased)) 
                  X <- X[, -which(object$aliased)]
                trend <- X %*% B
                signal <- rep(0, length(trend))
                res <- trend + signal
                attr(res, "trend") <- trend
                attr(res, "signal") <- signal
            }
            else stop("unkown error model etype")
        }
        else if (object$type == "mixed") {
            if (is.null(listw) || !inherits(listw, "listw")) 
                stop("spatial weights list required")
            if (nrow(newdata) != length(listw$neighbours)) 
                stop("mismatch between newdata and spatial weights")
            B <- object$coefficients
            frm <- formula(object$call)
            mt <- delete.response(terms(frm, data = newdata))
            mf <- model.frame(mt, newdata)
            X <- model.matrix(mt, mf)
            K <- ifelse(colnames(X)[1] == "(Intercept)", 2, 1)
            m <- ncol(X)
            if (m > 1) {
                WX <- matrix(nrow = nrow(X), ncol = (m - (K - 
                  1)))
                for (k in K:m) {
                  wx <- lag.listw(listw, X[, k], zero.policy = zero.policy)
                  if (any(is.na(wx))) 
                    stop("NAs in lagged independent variable")
                  WX[, (k - (K - 1))] <- wx
                }
            }
            if (K == 2) {
                if (!(listw$style == "W")) {
                  intercept <- as.double(rep(1, nrow(X)))
                  wx <- lag.listw(listw, intercept, zero.policy = zero.policy)
                  if (m > 1) {
                    WX <- cbind(wx, WX)
                  }
                  else {
                    WX <- matrix(wx, nrow = nrow(X), ncol = 1)
                  }
                }
            }
            X <- cbind(X, WX)
            if (any(object$aliased)) 
                X <- X[, -which(object$aliased)]
            trend <- X %*% B
            if (power) {
                W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
                res <- c(as(powerWeights(W, rho = object$rho, 
                  X = trend, order = order, tol = tol), "matrix"))
            }
            else {
                res <- c(invIrW(listw, object$rho) %*% trend)
            }
            if (legacy) {
                signal <- object$rho * lag.listw(listw, res, 
                  zero.policy = zero.policy)
                res <- c(trend + signal)
            }
            else {
                signal <- res - trend
            }
            attr(res, "trend") <- c(trend)
            attr(res, "signal") <- c(signal)
        }
        else {
            if (is.null(listw) || !inherits(listw, "listw")) 
                stop("spatial weights list required")
            if (nrow(newdata) != length(listw$neighbours)) 
                stop("mismatch between newdata and spatial weights")
            B <- object$coefficients
            frm <- formula(object$call)
            mt <- delete.response(terms(frm, data = newdata))
            mf <- model.frame(mt, newdata)
            if (dim(mf)[1] != length(listw$neighbours)) 
                stop("missing values in newdata")
            X <- model.matrix(mt, mf)
            if (any(object$aliased)) 
                X <- X[, -which(object$aliased)]
            trend <- X %*% B
            if (power) {
                W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
                res <- c(as(powerWeights(W, rho = object$rho, 
                  X = trend, order = order, tol = tol), "matrix"))
            }
            else {
                res <- c(invIrW(listw, object$rho) %*% trend)
            }
            if (legacy) {
                signal <- object$rho * lag.listw(listw, res, 
                  zero.policy = zero.policy)
                res <- c(trend + signal)
            }
            else {
                signal <- res - trend
            }
            attr(res, "trend") <- c(trend)
            attr(res, "signal") <- c(signal)
        }
    }
    class(res) <- "sarlm.pred"
    res
}
<environment: namespace:spdep>
