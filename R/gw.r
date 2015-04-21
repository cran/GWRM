gw <- function(formula, data, weights, k = NULL, subset, na.action, 
               kstart = 1, rostart = 2, betastart = NULL,offset, control = list(...), method = NULL, hessian = TRUE, model = TRUE, x = FALSE, y = TRUE, ...){
  
  if (is.null(control$trace)){
    control$trace = 0 
  }
  call <- match.call()
  if (missing(data))  
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  
  fitted<-FALSE
  
  Terms <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  nobs<-nrow(as.matrix(Y))
  X <- if(!is.empty.model(Terms)){ 
    model.matrix(Terms, mf, contrasts)
  }
  else{
    matrix(1, nobs, 0)
  }
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  
  if (!is.null(offset)) {
    if (length(offset) != nobs) 
      stop(gettextf("Number of offsets is %d should equal %d (number of observations)", length(offset), nobs), domain = NA)
  }
  if (is.null(method)){
    if (control$trace > 0) {
      cat("Trying 'nlm' initial fit", "\n")
    }
    nlm.fit <-try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = kstart, rostart = rostart, betastart = betastart, offset = offset, control = control, method = "nlm", hessian = hessian)), silent = TRUE)
    if ('try-error' %in% class(nlm.fit)){
      if (control$trace > 0) {
        cat("Crashed 'nlm' initial fit", "\n")
      }
      nlm.aic <- Inf
    }
    else {
      nlm.aic <- AIC(nlm.fit)
      if (control$trace > 0) {
        cat(paste(round(nlm.aic, 2), "AIC value in 'nlm' initial fit"), "\n")
      }
    }
    if (control$trace > 0) {
      cat("Trying 'Nelder-Mead' initial fit", "\n")
    }
    neldermead.fit <- try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = kstart, rostart = rostart, betastart = betastart, offset = offset, control = control, method = "Nelder-Mead", hessian = hessian)), silent = TRUE)
    
    if ('try-error' %in% class(neldermead.fit)){
      if (control$trace > 0) {
        cat("Crashed 'Nelder-Mead' initial fit", "\n")
      }
      neldermead.aic <- Inf
    }
    else {
      neldermead.aic <- AIC(neldermead.fit)
      if (control$trace > 0) {
        cat(paste(round(neldermead.aic,2), "AIC value in 'Nelder-Mead' initial fit"), "\n")
      }
    }
    if (neldermead.aic == Inf & nlm.aic == Inf){
      warning("No 'nlm' neither 'Nelder-Mead' provide fits. Try to change initial values")
    } 
    else { 
      if (neldermead.aic < nlm.aic | nlm.aic <= 0){
        fit2 <- neldermead.fit
      }
      else{
        fit2<-nlm.fit          
      }
      if (control$trace > 0) {
        cat("L-BFGS-B fitting","\n")
      }
      fit <- try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = fit2$betaIIpars[1], rostart = fit2$betaIIpars[2], betastart = fit2$betascoefs, offset=offset, control = control, method = "L-BFGS-B", hessian = hessian)),silent=TRUE)
      if ('try-error' %in% class(fit)){
        if (control$trace > 0) {
          cat("Crashed 'L-BFGS-B' final fit", "\n")
        }
        fit <- fit2
      }
      
      fitted <- TRUE
      
    }
  }
  else{
    if (!any(method == c("nlm", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
      stop("method must be in c('nlm', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN')")
    } 
    else{
      fit <- try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = kstart, rostart = rostart, betastart = betastart, control = control, method = method, hessian = hessian)), silent = TRUE)
      if ('try-error' %in% class(fit)){
        warning("Crashed fit", "\n")
      }
      else{
        fitted<-TRUE
      }
    }
  }
  if (fitted){
    if (x) 
      fit$X <- X
    if (!y) 
      fit$Y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = Terms, data = data, offset = offset, control = control, method = method, contrasts = attr(X, "contrasts"), xlevels = .getXlevels(Terms, mf)))
    class(fit) <- "gw"
  } 
  else{
    fit$aic<- -Inf
  }  
  fit
}


gw.fit <-function (x, y, weights = NULL, k = NULL, kstart = 1, rostart = 2, betastart = NULL, offset = NULL, control = list(), method = "L-BFGS-B", hessian=TRUE, intercept = TRUE){
    
  control <- do.call("gw.control", control)
  X <- x
  xnames<-dimnames(X)[[2L]]
  Y <- y    
  conv <- FALSE
  nobs <- nrow(as.matrix(X))
  ncovars <- ncol(as.matrix(X))
  if (is.null(weights)) 
    w <- rep(1, nobs)
  else{
    w <- weights
  }
  if (is.null(offset)){ 
    offset <- rep(0, nobs)
    covoffset <- FALSE
  }
  else{
    covoffset <- TRUE
  }
  if (is.null(kstart))
    kstart <- 1
  if (is.null(rostart))
    rostart <- 2
  if(!is.null(k)){ 
    kBool <- TRUE   #Boolean about k parameter 
    gCorrect <- 1   #To correct df, since there is a parameter less
  }
  else{  
    kBool <- FALSE
    gCorrect <- 0
  }
  
  if (kstart <= 0){
    stop("kstart must be positive")
  }
  
  if (rostart <= 1){
    stop("rostart must be greater than 1")
  }
  
  if (!is.logical(hessian)){
    stop("Hessian must be logical")
  }
  
  if (!any(method == c("nlm", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
    stop("method must be in c('nlm', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN')")
  }
  
  if (is.null(betastart)){
    betastart <- rep(0, ncovars)
  }
  
  #Log-likelihood and initial value
  if(!kBool){
    logL<-function(p){
      beta <- p[1:(ncovars)]
      betak <- p[ncovars + 1]
      betaro <- p[ncovars + 2]
      if (method != "L-BFGS-B"){ 
        k <- exp(betak)
        ro <- 1 + exp(betaro)
      }
      else{
        k <- betak
        ro <- betaro
      }
      mu <- exp(offset + X %*% beta)
      a <- mu * (ro - 1) / k
      gama <- a + k + ro
      -sum(w * (lgamma(a + Y) - lgamma(a) + lgamma(k + Y) - lgamma(k) - lgamma(gama + Y) + lgamma(a + ro) + lgamma(k + ro) - lgamma(ro)))
    }
    if (method != "L-BFGS-B"){
      p0 <- c(betastart, log(kstart), log(rostart - 1)) 
    }
    else{
      p0 <- c(betastart, kstart, rostart) 
    }
  }
  else{
    logL <- function(p){
      beta <- p[1:(ncovars)]
      betaro<-p[ncovars + 1] 
      if (method != "L-BFGS-B"){
        ro <- 1 + exp(betaro)
      }
      else{
        ro <- betaro
      }
      mu <- exp(offset + X %*% beta)
      a <- mu * (ro - 1) / k
      gama <- a + k + ro
      #If k=1 there is a simpler expresion
      if (k == 1){
        -sum(w * (lbeta(a + Y, ro + 1) - lbeta(a, ro)))
      }
      else{
        -sum(w * (lgamma(a + Y) - lgamma(a) + lgamma(k + Y) - lgamma(k) - lgamma(gama + Y) + lgamma(a + ro) + lgamma(k + ro) - lgamma(ro)))
      }
    }
    if (method != "L-BFGS-B"){ 
      p0 <- c(betastart, log(rostart - 1)) 
    }
    else{
      p0 <- c(betastart, rostart) 
    }
  }
  
  #Optimizing log-likelihood
  if (method == "nlm"){
    fit <- nlm(logL, p = p0, hessian = hessian, iterlim = control$maxit, print.level = control$trace)
    fit$value <- fit$minimum
    fit$par <- fit$estimate
    fit$convergence <- fit$code
    methodText = "nlm"
  }
  else if (any(method == c("Nelder-Mead", "BFGS", "CG","SANN"))){
    fit <- optim(p0, logL, method = method, hessian = hessian, control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
  }
  else if (any(method == c("L-BFGS-B"))){    
    #The constraints depends on whether k is known so
    if(!kBool){
      lower <- c(rep(-Inf, ncovars), 0.0000001, 1.0000001)
    }
    else{
      lower <- c(rep(-Inf, ncovars), 1.0000001)
    }
    fit <- optim(p0, logL, method = method, hessian = hessian, lower = lower, control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
  }
  
  if(!kBool){
    if (any(method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(xnames, "log(k)", "log(ro-1)"))
      betaIIpars <- c(exp(fit$par[ncovars + 1]), 1 + exp(fit$par[ncovars + 2]))
    }
    else{
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(xnames, "k", "ro"))
      betaIIpars <- c(fit$par[ncovars + 1], fit$par[ncovars + 2])
    }
  }
  else{
    if (any(method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(rep("beta", ncovars), "log(ro-1)"))
      dimnames(coef.table) <- list("", c(xnames, "log(ro-1)"))
      betaIIpars <- c(k, 1 + exp(fit$par[ncovars + 1]))
    }
    else{
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(xnames, "ro"))
      betaIIpars <- c(k, fit$par[ncovars + 1])
    }
  }
  
  
  results <- list(
    Y = Y,
    W = w,
    covars = dimnames(X)[[2]],
    nobs = sum(w),
    covoffset = covoffset,
    loglik = -(fit$value + sum(w * lfactorial(Y))),
    aic = 2 * (fit$value + sum(w * lfactorial(Y))) + (length(betastart) + 2 - gCorrect) * 2,
    bic = 2 * (fit$value + sum(w * lfactorial(Y))) + (length(betastart) + 2 - gCorrect) * log(sum(w)),
    df.residual = sum(w) - (length(betastart) + 2 - gCorrect),
    residuals = Y - exp(offset + X %*% fit$par[1:ncovars]),
    coefficients = coef.table,
    betaIIpars = betaIIpars,
    betascoefs = fit$par[1:ncovars],
    fitted.values = exp(offset + X %*% fit$par[1:ncovars]),
    hessian = fit$hessian,
    cov = solve(fit$hessian),
    se = sqrt(diag(solve(fit$hessian))),
    corr = solve(fit$hessian) / (sqrt(diag(solve(fit$hessian))) %o% sqrt(diag(solve(fit$hessian)))),
    code = fit$convergence,    
    method = methodText,
    k = k,
    kBool = kBool
  )
  
  class(results) <- "gw"
  return(results)
} 



print.gw<-function (x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts)) {
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    }
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  
  cat("\nDegrees of Freedom:", x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }
  cat("AIC:", format(signif(x$aic, digits)))
  cat("\n")
  invisible(x)
}

summary.gw <- function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...){
  
  df.r <- object$df.residual
  if (is.null(object$covars)) object$covars <- "(Intercept)"
  coef.p <- object$betascoefs
  s.err <- object$se[1:length(object$betascoefs)]
  tvalue <- (object$betascoefs) / (object$se[1:length(object$betascoefs)])
  pvalue <- 2 * pnorm(abs((object$betascoefs) / (object$se[1:length(object$betascoefs)])), lower.tail = FALSE)              
  dn <- c("Estimate", "Std. Error")
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(object$covars, c(dn, "z value", "Pr(>|z|)"))                     
  fitted <- cbind(object$loglik, object$aic, object$bic, object$df)
  dimnames(fitted) <- list("", c("log-likelihood", "AIC", "BIC", "df"))
  if (!object$kBool){
    if (any(object$method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coefk <- object$coefficients[length(object$betascoefs) + 1]
      coefro <- object$coefficients[length(object$betascoefs) + 2]
      k <- object$betaIIpars[1]
      ro <- object$betaIIpars[2]
      #Std. error aproximated by delta method
      SE.k <- object$se[length(object$betascoefs) + 1] * exp(coefk)
      SE.ro <- object$se[length(object$betascoefs) + 2] * exp(coefro)
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("",""), c("par", "Estimate", "Std. Error (Delta method)"))
    }
    else{
      k <- object$betaIIpars[1]
      ro <- object$betaIIpars[2]
      SE.k <- object$se[length(object$betascoefs) + 1]
      SE.ro <- object$se[length(object$betascoefs) + 2]
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("", ""), c("par", "Estimate", "Std. Error"))
    }
  }
  else{
    if (any(object$method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coefro <- object$coefficients[length(object$betascoefs) + 1]
      k <- 1
      ro <- object$betaIIpars[2]
      #Std. Error aproximated by delta method
      SE.k <- NA
      SE.ro <- object$se[length(object$betascoefs) + 1] * exp(coefro)
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("",""), c("par", "Estimate", "Std. Error"))
    }
    else{
      k <- 1
      ro <- object$betaIIpars[2]
      SE.k <- NA
      SE.ro <- object$se[length(object$betascoefs) + 1]
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("",""), c("par", "Estimate", "Std. Error"))
    }
  }
  
  keep <- match(c("call", "terms", "deviance", "aic", "contrasts", "df.residual", "na.action"), names(object), 0L)
  ans <- c(object[keep], list(coefficients = coef.table, fitted = fitted, betaII = betaII, method = object$method, convergence = object$code))
  class(ans) <- "summary.gw"
  return(ans)
}

print.summary.gw <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))){
    cat("Coefficients")
    if (is.character(co <- x$contrasts)){
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    }
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  if (length(x$fitted)){
    cat("\n")
    cat("Fit:\n")
    print.default(format(x$fitted, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No fits\n\n")
  if (length(x$betaII)){
    cat("\n")
    cat("betaII:\n")
    print.default(format(x$betaII, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No betaII\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
  cat("\nCode of convergence:", x$convergence, "\n")
  cat("\nMethod:", x$method, "\n")  
  if (nzchar(mess <- naprint(x$na.action))){
    cat("  (", mess, ")\n", sep = "")
  }
  invisible(x)
}

gw.control <- function (maxit = 10000, epsilon = 1e-08, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0) 
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit, trace = trace)
}

model.matrix.gw<-function (object, ...) {
  if (n_match <- match("x", names(object), 0L)) object[[n_match]]
  else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
  }
}

predict.gw <- function(object = NULL, newdata = NULL, ...){
  tt <- terms(object)
  if (!inherits(object, "gw")) 
    warning("calling predict.gw(<fake-gw-object>) ...")
  
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    nobs<-nrow(as.matrix(m))
    X <- model.matrix(Terms, m, offset <- rep(0, nobs))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  
  ncovars <- ncol(X)
  beta <- object$coefficients[1:(ncovars)]
  if(!object$kBool){
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
  }
  else{
    k <- object$k
    ro <- object$betaIIpars[2]
  }
  
  if (is.null(offset))
    fits <- exp(X %*% beta)
  else
    fits <- exp(offset + X %*% beta)
  predictor <- cbind(fits)
  colnames(predictor) <- c("fit")
  ans <- data.frame(predictor)
  return(ans)
}


partvar <- function(object, newdata = NULL, ...){
  UseMethod("partvar")
}

partvar.default <- function(object, newdata = NULL, ...){
  tt <- terms(object)
  if (!inherits(object, "gw")) 
    warning("calling predict.gw(<fake-gw-object>) ...")
  
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    nobs<-nrow(as.matrix(m))
    X <- model.matrix(Terms, m, offset <- rep(0, nobs))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  
  ncovars <- ncol(X)
  beta <- object$coefficients[1:(ncovars)]
  if(!object$kBool){
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
  }
  else{
    k <- object$k
    ro <- object$betaIIpars[2]
  }
  
  
  if (is.null(offset)){
    mus <- exp( X %*% beta)
  }
  else{
    mus <- exp(offset +X %*% beta)
  }
  
  a <- mus * (ro - 1) / k
  var <- mus * ((a + ro - 1) * (k + ro - 1)) / ((ro - 1) * (ro - 2))
  prand <- mus / var
  pliabi <- ((ro - 1) * (k + 1)) / ((a + ro - 1) * (k + ro - 1))
  pprone <- a / (a + ro - 1)
  rand <- mus
  liabi <- pliabi * var
  prone <- pprone * var
  partvar_rate = data.frame(Levels = X, Randomness = prand, Liability = pliabi, Proneness = pprone)
  partvar = data.frame(Levels = X, Randomness = rand, Liability = liabi, Proneness = prone)
  
  ans <- list(Prop.Variance.Components = partvar_rate[, c("Randomness", "Liability", "Proneness")], Variance.Components = partvar[, c("Randomness", "Liability", "Proneness")])
  return(ans)
}


logLik.gw <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$nobs
  attr(val, "df") <- length(object$coefficients)
  class(val) <- "logLik"
  val
}

residuals.gw <- function(object, type = "pearson", rep = 19, envelope = TRUE, title = "Simulated Envelope of Residuals", trace = FALSE, ...){ 
  
  if (envelope == TRUE){
    elapsedtime <- function (st, et){
      time <- et[3] - st[3]
      h <- trunc(time/3600)
      if (h < 10) 
        hs <- paste("0", h, sep = "", collapse = " ")
      else hs <- h
      time <- time - h * 3600
      min <- trunc(time/60)
      if (min < 10) 
        mins <- paste("0", min, sep = "", collapse = " ")
      else mins <- min
      time <- time - min * 60
      sec <- trunc(time)
      if (sec < 10) 
        secs <- paste("0", sec, sep = "", collapse = " ")
      else secs <- sec
      retval <- paste(hs, ":", mins, ":", secs, sep = "", collapse = " ")
      return(retval)
    }  
  }
  
  #if (is.null(object$W)) object$W <- rep(1, length(object$Y))
  
  st <- proc.time()
  
  # Model parameters
  k <- object$betaIIpars[1]
  ro <- object$betaIIpars[2]
  mu <- object$fitted.values
  a <- mu * (ro - 1) / k
  if (rep < 19){ 
    stop("Number of replications must equal or be greater than 19")    
  }
  
  if (type == "pearson"){
    if (k > 0 & ro > 2){
      varianza <- ((k + ro - 1) / (ro - 2)) * (mu + (mu ^ 2) / k)
      residuos <- sort(object$residuals / sqrt(varianza))
      if (envelope == TRUE){
        residuos.sim <- matrix(rep(NA, rep * sum(object$W)), nrow = sum(object$W), ncol = rep)
        residuos.sim[, 1] <- rep(residuos, object$W)
        i <- 1
        while (i <= (rep - 1)){
          datos <- object$data[rep(1:nrow(object$data), object$W),]
          datosSimulados <- as.matrix(rghyper(object$W, -a, -k, ro - 1))
          datos <- object$data[rep(1:nrow(object$data), object$W),]
          datos[getResponse(object$formula)] <- datosSimulados
          fit <- try(gw(object$formula, data = datos, k = object$k, method = object$methodCode), silent = TRUE)
          if ('try-error' %in% class(fit)){
            cat("Crashed fit", "\n")
          }
          else {
            if (trace > 0) {
              cat(paste(c("Fit", i, "out of", rep)), "\n")
            }
            media <- object$fitted.values
            k <- fit$betaIIpars[1]
            ro <- fit$betaIIpars[2]
            varianza <- ((k + ro - 1) / (ro - 2)) * (media + (media ^ 2) / k)
            residuos.i <- try(sort(fit$residuals / sqrt(varianza)), silent=TRUE)
            if ('try-error' %in% class(residuos.i)){
              if (trace > 0) {
                cat("Crashed fit", "\n")
              }
            }
            else{
              residuos.i <- sort(fit$residuals / sqrt(varianza))
              residuos.i <- rep(residuos.i, object$W)
              if (length(residuos.i) == sum(object$W)){
                residuos.sim[, i + 1] <- residuos.i
                i <- i + 1
              }
            }
          }
        }
      }
      no.res <- 0
    }
    else{
      warning("Residuals cannot be calculated")
      no.res <- 1
    }
  }
  
  if (type == "response"){
    no.res <- 0
    media <- object$fitted.values
    residuos <- sort(object$residuals)
    if (envelope == TRUE){
      residuos.sim <- matrix(rep(NA, rep * sum(object$W)), nrow = sum(object$W), ncol = rep)
      residuos.sim[, 1] <- rep(residuos, object$W)
      i <- 1
      while (i <= (rep - 1)){
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datosSimulados <- as.matrix(rghyper(object$W, -a, -k, ro - 1))
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datos[getResponse(object$formula)] <- datosSimulados
        fit <- try(gw(object$formula, data = datos, k = object$k, method = object$methodCode), silent = TRUE)#Le he quitado control = object$control
        if ('try-error' %in% class(fit)){
          cat("Crashed fit", "\n")
        }
        else {
          if (trace > 0) {
            cat(paste(c("Fit", i, "out of", rep)), "\n")
          }
          residuos.i <- sort(fit$residuals)
          residuos.i <- rep(residuos.i, object$W)
          if (length(residuos.i) == sum(object$W)){
            residuos.sim[, i + 1] <- residuos.i
            i <- i + 1
          }
        }
      }
    }
  }
  
  if (type == "deviance"){
    no.res <- 0
    media <- object$fitted.values
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
    a <- media * (ro - 1) / k
    gama <- a + k + ro
    y <- object$Y
    residuos <- sort(2 * (lgamma(a + y) + lgamma(k + y) - lgamma(gama + y) - (lgamma(a + media) + lgamma(k + media) - lgamma(gama + media))))
    if (envelope == TRUE){
      residuos.sim <- matrix(rep(NA, rep * sum(object$W)), nrow = sum(object$W), ncol = rep)
      residuos.sim[, 1] <- rep(residuos, object$W)
      i <- 1
      while (i <= (rep - 1)){
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datosSimulados <- as.matrix(rghyper(object$W, -a, -k, ro - 1))
        datos <- object$data[rep(1:nrow(object$data), object$W),]
        datos[getResponse(object$formula)] <- datosSimulados
        fit <- try(gw(object$formula, data = datos, k = object$k, method = object$methodCode), silent = TRUE)
        if ('try-error' %in% class(fit)){
          cat("Crashed fit", "\n")
        }
        else {
          if (trace > 0) {
            cat(paste(c("Fit", i, "out of", rep)), "\n")
          }
          media <- fit$fitted.values
          k <- fit$betaIIpars[1]
          ro <- fit$betaIIpars[2]
          a <- media * (ro - 1) / k
          gama <- a + k + ro
          y <- fit$Y
          residuos.i <- sort(2 * (lgamma(a + y) + lgamma(k + y) - lgamma(gama + y) - (lgamma(a + media) + lgamma(k + media) - lgamma(gama + media))))
          residuos.i <- rep(residuos.i, object$W)
          if (length(residuos.i) == sum(object$W)){
            residuos.sim[, i + 1] <- residuos.i
            i <- i + 1
          }
        }
      }
    }
  }
  
  if (envelope==TRUE & no.res != 1){
    minimos <- apply(residuos.sim[, 2:rep], 1, min)
    maximos <- apply(residuos.sim[, 2:rep], 1, max)
    n <- sum(object$W)
    t <- 1:n
    normal.score <- qnorm(t / (n + 1))
    xx <- c(normal.score, rev(normal.score))
    yy <- c(minimos, rev(maximos))
    plot(normal.score, residuos, type = "l", xlab = "Standard normal quantiles", ylab = paste("Residuals ","(",type[1],")", sep = ""), main = title)
    polygon(xx, yy, col = "gray", border = NA)
    lines(normal.score, residuos)
  }
  
  if (envelope == TRUE){
    et <- proc.time()
    if (trace > 0) 
      cat(paste("\nOverall envelope simulation process took", elapsedtime(st, et)), "(hh:mm:ss)\n")
    else cat("\n")
  }
  
  if (envelope == FALSE) residuos.sim <- NULL
  
  if (no.res != 1){
    ans <-  list(type = type, residuals = residuos, sim.residuals = residuos.sim)
    class(ans) <- "residuals.gw"
    return(ans)
  }
  else
    return("No residuals calculated")
}

getResponse <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] 
  response <- attr(tt, "response")
  vars[response] 
}

extractAIC.gw <- function (fit, scale, k = 2, ...){
  if (fit$aic < 0 || (fit$method != "nlm" && fit$code != 0)){
    c(Inf, Inf)
  }
  else{
    n <- length(fit$residuals)
    edf <- n - fit$df.residual
    aic <- fit$aic
    c(edf, aic + (k - 2) * edf)
  }
}

formula.gw<-function (x, ...) {
  form <- x$formula
  if (!is.null(form)) {
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
  }
  else formula(x$terms)
}

add1.gw <- function (object, scope, test = c("none", "Chisq"), k = 2, trace = FALSE, ...) {
  
  safe_pchisq <- function(q, df, ...){
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
  
  if (missing(scope) || is.null(scope)) 
    stop("no terms in scope")
  if (!is.character(scope)) 
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope)) 
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", scope), c("df", "AIC")))
  ans[1L, ] <- extractAIC(object, k = k, ...)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq(ns)) {
    tt <- scope[i]
    if (trace > 0) {
      cat("trying +", tt, "\n", sep = "")
      utils::flush.console()
    }
    nfit <- update(object, as.formula(paste("~ . +", tt)), evaluate = FALSE)
    nfit <- eval(nfit, envir = env)
    ans[i + 1L, ] <- extractAIC(nfit, k = k, ...)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if (all(is.finite(c(n0, nnew))) && nnew != n0) 
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[, 1L] - ans[1L, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2L])
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- ans[, 2L] - k * ans[, 1L]
    dev <- dev[1L] - dev
    dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term additions", "\nModel:", deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

drop1.gw<-function (object, scope, test = c("none", "Chisq"), k = 2, trace = FALSE, ...) {
  
  safe_pchisq <- function(q, df, ...){
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
  }
  
  tl <- attr(terms(object), "term.labels")
  if (missing(scope)) 
    scope <- drop.scope(object)
  else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), "term.labels")
    if (!all(match(scope, tl, 0L) > 0L)) 
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", scope), c("df", "AIC")))
  ans[1, ] <- extractAIC(object, k = k, ...)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq(ns)) {
    tt <- scope[i]
    if (trace > 0) {
      cat("trying -", tt, "\n", sep = "")
      utils::flush.console()
    }
    nfit <- update(object, as.formula(paste("~ . -", tt)), evaluate = FALSE)
    nfit <- eval(nfit, envir = env)
    ans[i + 1, ] <- extractAIC(nfit, k = k, ...)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if (all(is.finite(c(n0, nnew))) && nnew != n0) 
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[1L, 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2])
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- ans[, 2L] - k * ans[, 1L]
    dev <- dev - dev[1L]
    dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

model.frame.gw <- function (formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) 
      env <- parent.frame()
    eval(fcall, env, parent.frame())
  }
  else formula$model
}
