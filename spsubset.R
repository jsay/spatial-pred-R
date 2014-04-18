
nbsubset <- function (x, focal= rep(TRUE, length(x)),
                      target= rep(TRUE, length(x)), ...){
    if (!inherits(x, "nb")) 
        stop("not a neighbours list")
    if (!is.logical(focal)) 
        stop("subset not a logical vector")
    if (!is.logical(target)) 
        stop("target not a logical vector")
    n <- length(x)
    if (n != length(focal)) 
        stop("neighours list and subset focal different lengths")
    if (n != length(target)) 
        stop("neighours list and subset target different lengths")
    tgt.ids <- seq(1, n)[ target]
    x <- sym.attr.nb(x)
    xattrs <- names(attributes(x))
    z <- ifelse(focal, x, 0L)
    for (i in seq(along= 1: n)[ focal]){
        zi <- z[[ i]]
        zn <- zi[zi %in% tgt.ids]
        if (length(zn)== 0) z[[ i]] <- 0L
        else z[[ i]] <- sort(unique(zn))
    }
    attr(z, "region.id") <-attr(x, "region.id")
    for (i in 1:length(xattrs)) {
        if (xattrs[i] != "region.id") 
            attr(z, xattrs[i]) <- attr(x, xattrs[i])
    }
    z <- sym.attr.nb(z)
    z
}

lwsubset <- function (x, focal= rep(TRUE, length(x)),
                      target= rep(TRUE, length(x)), zero.policy = NULL, ...){
    if (!inherits(x, "listw")) 
        stop("not a weights list")
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (!is.logical(focal)) 
        stop("focal not a logical vector")
    if (!is.logical(target)) 
        stop("target not a logical vector")
    nb <- x$neighbours
    vlist <- x$weights
    if (attr(vlist, "mode") != "binary") 
        stop("Not yet able to subset general weights lists")
    style <- x$style
    n <- length(nb)
    if (n != length(focal)) 
        stop("neighbours list and focal vector different lengths")
    if (n != length(target)) 
        stop("neighbours list and target vector different lengths")
    subnb <- nbsubset(nb, focal, target)
    sublistw <- nb2listw(neighbours= subnb, glist= NULL,
                         style= style, zero.policy= zero.policy)
    sublistw
}
