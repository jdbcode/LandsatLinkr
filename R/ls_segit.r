#' Piece-wise regression analysis
#'
#' Piece-wise regression analysis
#' @param obj ?look up segmented
#' @param seg.z ?look up segmented
#' @param psi ?look up segmented
#' @param control ?look up segmented
#' @param model ?look up segmented
#' @import zoo
#' @import segmented
#' @export

ls_segit = function (obj, seg.Z, psi = stop("provide psi"), control = seg.control(), model = TRUE, ...){
  n.Seg <- 1
  if (length(all.vars(seg.Z)) > 1 & !is.list(psi)) 
    stop("`psi' should be a list with more than one covariate in `seg.Z'")
  if (is.list(psi)) {
    if (length(all.vars(seg.Z)) != length(psi)) 
      stop("A wrong number of terms in `seg.Z' or `psi'")
    if (any(is.na(match(all.vars(seg.Z), names(psi), nomatch = NA)))) 
      stop("Variables in `seg.Z' and `psi' do not match")
    n.Seg <- length(psi)
  }
  if (length(all.vars(seg.Z)) != n.Seg) 
    stop("A wrong number of terms in `seg.Z' or `psi'")
  it.max <- old.it.max <- control$it.max
  toll <- control$toll
  visual <- control$visual
  stop.if.error <- control$stop.if.error
  n.boot <- control$n.boot
  size.boot <- control$size.boot
  gap <- control$gap
  random <- control$random
  pow <- control$pow
  visualBoot <- FALSE
  if (n.boot > 0) {
    if (!is.null(control$seed)) {
      set.seed(control$seed)
      employed.Random.seed <- control$seed
    }
    else {
      employed.Random.seed <- eval(parse(text = paste(sample(0:9, 
                                                             size = 6), collapse = "")))
      set.seed(employed.Random.seed)
    }
    if (visual) {
      visual <- FALSE
      visualBoot <- TRUE
    }
    if (!stop.if.error) 
      stop("Bootstrap restart only with a fixed number of breakpoints")
  }
  last <- control$last
  K <- control$K
  h <- min(abs(control$h), 1)
  if (h < 1) 
    it.max <- it.max + round(it.max/2)
  orig.call <- Call <- mf <- obj$call
  orig.call$formula <- mf$formula <- formula(obj)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  if (class(mf$formula) == "name" && !"~" %in% paste(mf$formula)) 
    mf$formula <- eval(mf$formula)
  mfExt <- mf
  mf$formula <- update.formula(mf$formula, paste(seg.Z, collapse = ".+"))
  if (!is.null(obj$call$offset) || !is.null(obj$call$weights) || 
        !is.null(obj$call$subset)) {
    mfExt$formula <- update.formula(mf$formula, paste(".~.+", 
                                                      paste(paste(all.vars(obj$call$offset), collapse = "+"), 
                                                            paste(all.vars(obj$call$weights), collapse = "+"), 
                                                            paste(all.vars(obj$call$subset), collapse = "+"), 
                                                            sep = "+"), sep = ""))
  }
  mf <- eval(mf, parent.frame())
  n <- nrow(mf)
  nomiOff <- setdiff(all.vars(formula(obj)), names(mf))
  if (length(nomiOff) >= 1) 
    mfExt$formula <- update.formula(mfExt$formula, paste(".~.+", 
                                                         paste(nomiOff, collapse = "+"), sep = ""))
  nomiTUTTI <- all.vars(mfExt$formula)
  nomiNO <- NULL
  for (i in nomiTUTTI) {
    r <- try(eval(parse(text = i), parent.frame()), silent = TRUE)
    if (class(r) != "try-error" && length(r) == 1) 
      nomiNO[[length(nomiNO) + 1]] <- i
  }
  if (!is.null(nomiNO)) 
    mfExt$formula <- update.formula(mfExt$formula, paste(".~.-", 
                                                         paste(nomiNO, collapse = "-"), sep = ""))
  mfExt <- eval(mfExt, parent.frame())
  weights <- as.vector(model.weights(mf))
  offs <- as.vector(model.offset(mf))
  mt <- attr(mf, "terms")
  interc <- attr(mt, "intercept")
  y <- model.response(mf, "any")
  XREG <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  namesXREG0 <- colnames(XREG)
  nameLeftSlopeZero <- setdiff(all.vars(seg.Z), names(coef(obj)))
  namesXREG0 <- setdiff(namesXREG0, nameLeftSlopeZero)
  id.duplic <- match(all.vars(formula(obj)), all.vars(seg.Z), 
                     nomatch = 0) > 0
  if (any(id.duplic)) {
    new.mf <- mf[, all.vars(formula(obj))[id.duplic], drop = FALSE]
    new.XREGseg <- data.matrix(new.mf)
    XREG <- cbind(XREG, new.XREGseg)
  }
  n.psi <- length(unlist(psi))
  id.n.Seg <- (ncol(XREG) - n.Seg + 1):ncol(XREG)
  XREGseg <- XREG[, id.n.Seg, drop = FALSE]
  XREG <- XREG[, match(c("(Intercept)", namesXREG0), colnames(XREG), 
                       nomatch = 0), drop = FALSE]
  XREG <- XREG[, unique(colnames(XREG)), drop = FALSE]
  n <- nrow(XREG)
  Z <- lapply(apply(XREGseg, 2, list), unlist)
  name.Z <- names(Z) <- colnames(XREGseg)
  if (length(Z) == 1 && is.vector(psi) && (is.numeric(psi) || 
                                             is.na(psi))) {
    psi <- list(as.numeric(psi))
    names(psi) <- name.Z
  }
  if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) || 
        is.null(names(psi))) 
    stop("Z and psi have to be *named* list")
  id.nomiZpsi <- match(names(Z), names(psi))
  if ((length(Z) != length(psi)) || any(is.na(id.nomiZpsi))) 
    stop("Length or names of Z and psi do not match")
  nome <- names(psi)[id.nomiZpsi]
  psi <- psi[nome]
  initial.psi <- psi
  for (i in 1:length(psi)) {
    if (any(is.na(psi[[i]]))) 
      psi[[i]] <- if (control$quant) {
        quantile(Z[[i]], prob = seq(0, 1, l = K + 2)[-c(1, 
                                                        K + 2)], names = FALSE)
      }
    else {
      (min(Z[[i]]) + diff(range(Z[[i]])) * (1:K)/(K + 
                                                    1))
    }
  }
  a <- sapply(psi, length)
  id.psi.group <- rep(1:length(a), times = a)
  Z <- matrix(unlist(mapply(function(x, y) rep(x, y), Z, a, 
                            SIMPLIFY = TRUE)), nrow = n)
  psi <- unlist(psi)
  psi <- unlist(tapply(psi, id.psi.group, sort))
  k <- ncol(Z)
  PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
  c1 <- apply((Z <= PSI), 2, all)
  c2 <- apply((Z >= PSI), 2, all)
  if (sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) 
    stop("starting psi out of the admissible range")
  colnames(Z) <- nomiZ <- rep(nome, times = a)
  ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))], 
                                          function(xxx) {
                                            1:xxx
                                          })))
  nomiU <- paste("U", ripetizioni, sep = "")
  nomiU <- paste(nomiU, nomiZ, sep = ".")
  nomiV <- paste("V", ripetizioni, sep = "")
  nomiV <- paste(nomiV, nomiZ, sep = ".")
  if (it.max == 0) {
    U <- pmax((Z - PSI), 0)
    colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
    nomiU <- paste("U", colnames(U), sep = "")
    for (i in 1:ncol(U)) mfExt[nomiU[i]] <- mf[nomiU[i]] <- U[, 
                                                              i]
    Fo <- update.formula(formula(obj), as.formula(paste(".~.+", 
                                                        paste(nomiU, collapse = "+"))))
    obj <- update(obj, formula = Fo, evaluate = FALSE)
    obj <- eval(obj, envir = mfExt)
    if (model) 
      obj$model <- mf
    names(psi) <- paste(paste("psi", ripetizioni, sep = ""), 
                        nomiZ, sep = ".")
    obj$psi <- psi
    return(obj)
  }
  if (is.null(weights)) 
    weights <- rep(1, n)
  if (is.null(offs)) 
    offs <- rep(0, n)
  initial <- psi
  obj0 <- obj
  dev0 <- sum(obj$residuals^2)
  list.obj <- list(obj)
  nomiOK <- nomiU
  opz <- list(toll = toll, h = h, stop.if.error = stop.if.error, 
              dev0 = dev0, visual = visual, it.max = it.max, nomiOK = nomiOK, 
              id.psi.group = id.psi.group, gap = gap, visualBoot = visualBoot, 
              pow = pow)
  if (n.boot <= 0) {
    obj <- seg.lm.fit(y, XREG, Z, PSI, weights, offs, opz)
  } else {
    obj <- seg.lm.fit.boot(y, XREG, Z, PSI, weights, offs, 
                           opz, n.boot = n.boot, size.boot = size.boot, random = random)
  }
  if (!is.list(obj)) {
    warning("No breakpoint estimated", call. = FALSE)
    return(obj0)
  }
  if (obj$obj$df.residual == 0) 
    warning("no residual degrees of freedom (other warnings expected)", 
            call. = FALSE)
  id.psi.group <- obj$id.psi.group
  nomiOK <- obj$nomiOK
  it <- obj$it
  psi <- obj$psi
  psi.values <- if (n.boot <= 0) obj$psi.values else obj$boot.restart
  U <- obj$U
  V <- obj$V
  for (jj in colnames(V)) {
    VV <- V[, which(colnames(V) == jj), drop = FALSE]
    sumV <- abs(rowSums(VV))
    if ((any(diff(sumV) >= 2) || any(table(sumV) <= 1)) && 
          stop.if.error) 
      stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
  }
  rangeZ <- obj$rangeZ
  obj <- obj$obj
  k <- length(psi)
  beta.c <- if (k == 1) coef(obj)["U"] else coef(obj)[paste("U", 1:ncol(U), sep = "")]
  psi.values[[length(psi.values) + 1]] <- psi
  id.warn <- FALSE
  if (n.boot <= 0 && it > it.max) {
    warning("max number of iterations attained", call. = FALSE)
    id.warn <- TRUE
  }
  Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
  length.psi <- tapply(as.numeric(as.character(names(psi))), 
                       as.numeric(as.character(names(psi))), length)
  forma.nomiU <- function(xx, yy) paste("U", 1:xx, ".", yy, 
                                        sep = "")
  forma.nomiVxb <- function(xx, yy) paste("psi", 1:xx, ".", 
                                          yy, sep = "")
  nomiU <- unlist(mapply(forma.nomiU, length.psi, name.Z))
  nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, name.Z))
  for (i in 1:ncol(U)) {
    mfExt[nomiU[i]] <- mf[nomiU[i]] <- U[, i]
    mfExt[nomiVxb[i]] <- mf[nomiVxb[i]] <- Vxb[, i]
  }
  nnomi <- c(nomiU, nomiVxb)
  Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", 
                                                       paste(nnomi, collapse = "+"))))
  objF <- update(obj0, formula = Fo, evaluate = FALSE, data = mfExt)
  objF <- eval(objF, envir = mfExt)
  if (any(is.na(objF$coefficients)) && stop.if.error) {
    stop("at least one coef estimate is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", 
         call. = FALSE)
  }
  objF$offset <- obj0$offset
  if (!gap) {
    names.coef <- names(objF$coefficients)
    if (k == 1) {
      names(obj$coefficients)[match(c("U", "V"), names(coef(obj)))] <- nnomi
    }
    else {
      names(obj$coefficients)[match(c(paste("U", 1:k, sep = ""), 
                                      paste("V", 1:k, sep = "")), names(coef(obj)))] <- nnomi
    }
    objF$coefficients[names.coef] <- obj$coefficients[names.coef]
    objF$fitted.values <- obj$fitted.values
    objF$residuals <- obj$residuals
  }
  if (any(is.na(objF$coefficients))) {
    stop("some estimate is NA: premature stopping with a large number of breakpoints?", 
         call. = FALSE)
  }
  Cov <- vcov(objF)
  id <- match(nomiVxb, names(coef(objF)))
  if (any(is.na(id))){return(control$K-1)}
  if (ncol(Cov) < max(id)){return(control$K-1)}
  vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
  a <- tapply(id.psi.group, id.psi.group, length)
  initial <- unlist(mapply(function(x, y) {
    if (is.na(x)[1]) 
      rep(x, y)
    else x
  }, initial.psi, a))
  psi <- cbind(initial, psi, sqrt(vv))
  rownames(psi) <- colnames(Cov)[id]
  colnames(psi) <- c("Initial", "Est.", "St.Err")
  objF$rangeZ <- rangeZ
  objF$psi.history <- psi.values
  objF$psi <- psi
  objF$it <- (it - 1)
  objF$epsilon <- obj$epsilon
  #objF$call <- match.call()
  objF$nameUV <- list(U = nomiU, V = rownames(psi), Z = name.Z)
  objF$id.group <- if (length(name.Z) <= 1) -rowSums(as.matrix(V))
  objF$id.psi.group <- id.psi.group
  objF$id.warn <- id.warn
  objF$orig.call <- orig.call
  if (model) 
    objF$model <- mf
  if (n.boot > 0) 
    objF$seed <- employed.Random.seed
  class(objF) <- c("segmented", class(obj0))
  list.obj[[length(list.obj) + 1]] <- objF
  class(list.obj) <- "segmented"
  if (last) 
    list.obj <- list.obj[[length(list.obj)]]
  return(list.obj)
}