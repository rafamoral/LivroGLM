## Método Delta
confint.glm.dose <- function(object, parm, level = 0.95, ...) {
  if(missing(parm)) parm <- seq_along(object)
  if(!is.numeric(parm)) stop("only numeric 'parm' indicators are supported")
  nam <- names(object)[parm]
  se <- attr(object, "SE")[parm]
  p <- attr(object, "p")[parm]
  object <- as.vector(object[parm])
  z <- sqrt(qchisq(level, 1))
  res <- cbind(lower = object - z*se,
               estimate = object,
               upper = object + z*se)
  row.names(res) <- paste("p = ", format(p), ":", sep = "")
  res
}

## Método de Fieller
Fieller <- local({
  ofun <- function(xt, b, xi, v, chi2)
    (b[1] + b[2]*xt - xi)^2/(v[1,1] + 2*xt*v[1,2] + xt^2*v[2,2]) - chi2
  
  ci <- function(b, xi, v, chi2) {
    xhat <- (xi - b[1])/b[2]
    xl <- xhat-1
    while(ofun(xl, b, xi, v, chi2) < 0) xl <- xl - 1
    low <- uniroot(ofun, interval = c(xl, min(xhat, xl+1)), b=b, xi=xi, v=v, chi2=chi2)$root
    xu <- xhat+1
    while(ofun(xu, b, xi, v, chi2) < 0) xu <- xu + 1
    upp <- uniroot(ofun, interval = c(max(xhat, xu-1), xu), b=b, xi=xi, v=v, chi2=chi2)$root
    c(lower = as.vector(low), estimate = as.vector(xhat), upper = as.vector(upp))
  }
  
  function(object, cf = 1:2, p = 0.5, level = 0.95) {
    b <- coef(object)[cf]
    V <- vcov(object)[cf, cf]
    xiv <- family(object)$linkfun(p)
    chi2 <- qchisq(level, df = 1)
    
    R <- NULL
    for(xi in xiv)
      R <- rbind(R, ci(b, xi, V, chi2))
    row.names(R) <- paste("p = ", format(p), ":", sep = "")
    structure(R, p = p, class = "Fieller")
  }
})

print.Fieller <- function(x, ...) {
  attr(x, "p") <- class(x) <- NULL
  NextMethod("print", x, ...)
}

## Método da razão de verossimilhanças
LR_confint <- local({
  ofun <- Vectorize(function(theta_p) {
    glm.fit(x = cbind(X0, x-theta_p), y = Y, weights = wts,
            etastart = etastart, offset = off, family = fam,
            control = control, intercept = TRUE)$deviance - D0 - chi2
  }, "theta_p")
  
  ci <- function(b, xi) {
    theta_p <- (xi - b[1])/b[2]
    theta_p_l <- theta_p-1
    while(ofun(theta_p_l) < 0) theta_p_l <- theta_p_l - 1
    low <- uniroot(ofun, interval = c(theta_p_l, min(theta_p, theta_p_l+1)))$root
    theta_p_u <- theta_p+1
    while(ofun(theta_p_u) < 0) theta_p_u <- theta_p_u + 1
    upp <- uniroot(ofun, interval = c(max(theta_p, theta_p_u-1), theta_p_u))$root
    c(lower = as.vector(low), estimate = as.vector(theta_p), upper = as.vector(upp))
  }
  
  prof <- function(b, xi, R, ...) {
    R <- R[1,]
    theta_p <- (xi - b[1])/b[2]
    inc <- diff(range(R))
    dc <- D0 + chi2
    dev.x <- seq(theta_p - inc, theta_p + inc, length=50)
    dev.y <- ofun(dev.x) + dc
    plot(dev.x, dev.y, type="l", ylab="Deviance", las=1, ...)
    abline(h=dc, lty=2)
    arrows(R, c(0,0,0), R, c(dc,min(dev.y),dc), 
           code = 3, length = 0, lty = 2)
    points(R[c(1,3)], rep(par("usr")[3], 2), xpd = TRUE, pch = 16, cex = .8)
    points(R[2], par("usr")[3], xpd = TRUE, pch = "*", cex = 1.8)
    #axis(1, at = R[c(1,3)], labels = FALSE)
    text(x = R[c(1,3)], par("usr")[3], cex = .7,
         labels = c("limite inferior","limite superior"), 
         xpd = TRUE, srt = 45, pos = 1)
  }
  
  function(object, cf = 1:2, p = 0.5, level = 0.95, profile = FALSE, ...) {
    fam <<- family(object)
    Y <<- object$y
    X <- model.matrix(object)
    X0 <<- X[, -cf]
    x <<- X[, cf[2]]
    b <- as.vector(coef(object)[cf])
    etastart <<- object$linear.predictors
    wts <<- weights(object)
    originalOffset <<- if(is.null(o <- object$offset)) 0 else o
    control <<- object$control
    xiv <- fam$linkfun(p)
    chi2 <<- qchisq(level, 1)
    D0 <<- deviance(object)
    R <- NULL
    for(xi in xiv) {
      off <<- originalOffset + xi
      R <- rbind(R, ci(b, xi))
    }
    if(profile) {
      if(length(p) > 1) cat("Profile produced only for p = ", p[1], sep="", "\n")
      xi <- xiv[1]
      off <<- originalOffset + xi
      prof(b, xi, R, ...)
    }
    row.names(R) <- paste("p = ", format(p), ":", sep = "")
    structure(R, p = p, class = "LR_glm_dose")
  }
})

print.LR_glm_dose <- function(x, ...) {
  attr(x, "p") <- class(x) <- NULL
  NextMethod("print", x, ...)
}