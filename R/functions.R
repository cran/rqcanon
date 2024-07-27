#' Canonical Quantile Regression.
#'
#' Given multivariate data matrices X (explanatory variables) and Y (response variables), the function fits coefficients of the Y-variables that are best fit by a quantile regression on X. These are analogous to the coefficients given by a classical canonical correlations analysis, but replace the implicit L2 norm by an L1 norm. See: "Method" and "Reference" below.
#'
#' @param X input design matrix of explanatory variables (without intercept)
#' @param Y input matrix of response variables (nrow(Y) = nrow(X))
#' @param tau desired quantile, default = .5
#' @param a.pos for first component: non-empty vector of indices of Y-variables (Y-columns) whose alpha coefficients are constrained to be positive (to provide the direction of increasing responses); default = 1
#' @param ap for subsequent components j = 2, 3, . . . , na: vector whose (j-1)-th element is the Y-variable (column) index whose alpha coefficients are constrained to be positive; default = rep(1,na-1)
#' @param na number of components desired (1 <= na <= ncol(Y))
#' @param wts used only for use with the bootstrap methods. If weighting is desired for the sample observations, rqcan will multiply the columns of X and Y by the vector wts; but the bootstrap methods will apply to the unweighted data, and so will be incorrect; default = rep(1,nrow(X)) (unweighted analysis)
#'
#' @details
#' Finds orthogonal alpha coefficients and corresponding best-fitting beta coefficients to minimize sum|x_i' beta - y_i' alpha| subject to sum|alpha| = 1 (where x_i and y_i are the i-th rows of X and Y). The intercept is included (X should not include intercept). Need ncol(Y) > 1. For first component: if length(a.pos) < ncol(Y), sum(|alpha|) = 1  is constrained by going through all sign choices (s_j = sign(alpha_j)) and setting Y1_j =  s_j Y_j  (j not in a.pos). A constrained regression quantile fit is applied from quantreg: rq.fit.fnc(cbind(1,X,Y1),y0=0,R,r,tau). where (R,r) constrains all alpha_j >= 0 and sum(alpha_j) >= 1 (sum = 1 at min). Note: rq.fit.fnc solves by generating a sequence of quadratic approximations. The matrix defining one quadratic problem may be singular (and stop the computation) even if he input design matrices are of full rank. If a singularity stop occurs, jittering the data (see jitter()) sometimes helps.For the subsequent j-th component, only the index given by ap(j-1) is constrained to be positive. Alpha coefficients for subsequent components are constrained to be orthogonal to previous alpha coefficients. 
#'
#' @return object of class "rqcan"; a list of matrices of the alpha and beta coefficients: the j-th row of each matrix is the coefficients for the j-th component; input data and the constraint matrices R an r are also returned in the list
#' @references S. Portnoy, 2022. Canonical quantile regression, J. Multivar. Anal., 192, 105071.
#' @seealso See \link{summary.rqcan} for a description of the summary function.
#' @examples
#' X <- as.matrix(example_data[,1:3])
#' Y <- as.matrix(example_data[,4:7])
#'
#' a <- rqcan(X,Y,tau=.75,a.pos=2)
#' summary(a)
#' @export
rqcan <- function(X,Y,tau=.5,a.pos=1,ap=rep(1,na),na=ncol(Y),wts=rep(1,nrow(X)) ) {
  # finds orthog alphas and corresponding best-fitting betas s.t. sum|alpha| = 1
  # intercept included (X should not include 1st col = 1) ; need ncol(Y) > 1
  # a.pos = indices of Y-variables constrained to be >= 0 (ie, cols of Y matrix)
  # after 1st can component sign determined by constraining alpha[ap] > 0
  # na = number of orthog alphas to find; must be <= ncol(Y) & > 1
  # length(a.pos) should be at least 1 to specify coef signs if tau=.5
  #   if tau = .5 and a.pos = NULL,  coef and -coef give the same solution
  # for length(a.pos) < ncol(Y), constrain  sum(|alpha|) = 1  by setting
  # Y1_j =  s_j Y_j  (j !in a.pos) where s_j = sgn(alpha_j) ; use all sign choices
  # apply constrained rq.fit.fnc( cbind(1,X,Y1),y0=0,R,r,tau)
  # (R,r) constrains all alpha_j >= 0 and sum(alpha_j) >= 1 (makes sum = 1)
  #
  # initialization
  py <- ncol(Y) ; px <- ncol(X)
  if(nrow(X) != nrow(Y)) { warning("Row dimensions of X and Y must be equal.") }
  if(na > py) { warning("can not request more components than nrow(Y)") }
  if(na < 2 & na != 1) { na <- 1 ; warning("na must be a pos integer") }
  alps <- matrix(0,na,py) ; bets <- matrix(0,na,px+1)
  rownames(alps) <- rownames(bets) <- paste(1:na)
  colnames(alps) <- colnames(Y) ; colnames(bets) <- c("int",colnames(X))
  #
  #  use a output from rqcan1
  a <- rqcan1(X,Y,tau,a.pos,wts)
  if( !is.list(a) ) {
    warning("singular matrix input or generated, a='sing' returned; perhaps jitter X and Y")
    return(a) }
  R <- a$R ; r <- a$r
  alps[1,] <- a$alpha ; bets[1,] <- a$beta
  if(max(a$alpha) < .0000001) {
    warning("need at least 1 alpha coef > 0, Y set to -Y")
    a <- rqcan1(X,-Y,tau,a.pos,wts)
    R <- a$R ; r <- a$r
    alps[1,] <- a$alpha ; bets[1,] <- a$beta }
  if(na == 1)
  {rqcanresults<-list(alpha=alps,beta=bets,X=X,Y=Y,tau=tau,ap=ap,na=na,a.pos=a.pos,R=R,r=r)
  class(rqcanresults)<-"rqcan"
  return(rqcanresults)}
  # for subsequent component, constrain alpha[ap] > 0
  # old: if(length(a.pos) == 0) ap <- which.max(a$alpha)[1]
  #     else ap <- which.max(a$alpha[a.pos])[1]
  #
  # compute orthog canonical rq's ; need to orthogonalize and solve by rqcan1
  for(j in 2:na) {
    g <- qr.Q(qr(t(alps)),complete=TRUE)
    if(j == py) { # use orthog alpha and best-fitting beta
      b <- c(g[,py])
      b <- b/sum(abs(b))
      if(b[ap[j-1]] < 0) b <- -b
      alps[j,] <- b
      Y1 <- Y %*% b
      bets[j,] <- quantreg::rq.fit.br(cbind(1,X),Y1,tau)$coef
    } else { # orthogonalize and use rqcan1 with a.pos = ap
      g <- g[,j:py]
      Y1 <- Y %*% g
      a <- rqcan1(X,Y1,tau,a.pos = ap[j-1],wts)
      if(!is.list(a)) {
       warning("singular matrix input or generated, a='sing' returned; perhaps jitter X and Y")
       return(a) }
      b <- c(g %*% a$alpha)
      if(b[ap[j-1]] < 0) b <- -b
      b <- b/sum(abs(b))
      alps[j,] <- b
      Y1 <- Y %*% b
      bets[j,] <- quantreg::rq.fit.br(cbind(1,X),Y1,tau)$coef }
  }
  rqcanresults<-list(alpha=alps,beta=bets,X=X,Y=Y,tau=tau,ap=ap,na=na,a.pos=a.pos,R=R,r=r)
  class(rqcanresults)<-"rqcan"
  return(rqcanresults) }

NULL

# Document rqcan1
#' @export rqcan1
#' @title First Component
#' @description Internal function to find the first component \cr
#'  It is not intended for general use, but the documentation may be helpful if errors occur or if one whishes to modify the algorithms 
#' \cr
#'
#' @param X input X-matrix
#' @param Y input Y-matrix
#' @param tau probability for qualtile (default = .5)
#' @param a.pos indices of Y-variable whose coefficient is constrained to be positive (default = 1)
#' @param wts case weights (default = rep(1,nrow(X)) ) \cr
#' @details
#' The function finds the leading pair of indices. Notes: an intercept is added  (X should not include 1st col = 1) ;  ncol(Y) should be > 1 ; length(a.pos) should be at least 1 to specify coef signs (if tau = .5 and a.pos = NULL,  coef and -coef give the same solution) ; for length(a.pos) < ncol(Y), the constraint  sum(|alpha|) = 1  is set by setting Y1_j =  s_j Y_j  (j !in a.pos) where s_j = sgn(alpha_j) ;  all sign choices are used and then constrained rq.fit.fnc( cbind(1,X,Y1),y0=0,R,r,tau) is applied \cr (R,r) contrains all alpha_j >= 0 and sum(alpha_j) >= 1 (makes sum = 1)
#'
#' @return Returns list(a,X,Y,a.pos,R,r,rho1): a = output from rq.fit.fnc(XY,y0,R,r,tau) ; X,Y,a.pos = input data ; R,r = constraint matrices for rq.fit.fnc ; rho1 = rq objective fct. ; if rq.fit.fnc generates a singular matrix, returns "sing" \cr
#' \cr

rqcan1 <- function(X,Y,tau=.5,a.pos = 1,wts=rep(1,nrow(X)) ) {
  # finds leading pair of indices
  # intercept added (X should not include 1st col = 1) ; ncol(Y) > 1
  # a.pos = indices of Y-variables constrained to be >= 0 (ie, cols of Y matrix)
  # length(a.pos) should be at least 1 to specify coef signs
  # if tau = .5 and a.pos = NULL,  coef and -coef give the same solution
  # for length(a.pos) < ncol(Y), constrain  sum(|alpha|) = 1  by setting
  # Y1_j =  s_j Y_j  (j !in a.pos) where s_j = sgn(alpha_j) ; use all sign choices
  # apply constrained rq.fit.fnc( cbind(1,X,Y1),y0=0,R,r,tau)
  # (R,r) contrains all alpha_j >= 0 and sum(alpha_j) >= 1 (makes sum = 1)
  # if rq.fit.fnc generates a singular matrix, return "sing"
  px <- ncol(X) ; py <- ncol(Y) ; m <- length(a.pos) ; m1 <- py - m
  y0 <- rep(0,nrow(X))
  for(i in 1:(2^m1)) { # check each set of signs, note: if py = m, need s1=rep(1,py)
    s1 <- rep(1,py)
    if(m1 > 0) { s <- 1 - 2*as.integer(intToBits(i-1))[1:m1]
    s1 <- replace(s1,setdiff(1:py,a.pos),s) }
    R <- cbind(matrix(0,py,px+1),diag(s1))
    R <- rbind(R,c(rep(0,px+1),s1))
    r <- c(rep(0,py),1)
    XY <- wts * cbind(1,X,-Y)
    a <- tryCatch(quantreg::rq.fit.fnc(XY,y0,R,r,tau), error = function(e) TRUE)
    if(!is.list(a)) return("sing")
    I <- a$res < 0 ; rho <- sum((1-tau)*a$res[!I]) - sum(tau*a$res[I])
    if(i == 1) { as <- a ; rho1 <- rho }
    else if(rho < rho1) { as <- a ; rho1 <- rho }  }
  a <- as
  a$alpha <- a$coef[(px+2):(px+py+1)]
  a$beta <- a$coef[1:(px+1)]
  a$coef <- c(a$beta,a$alpha)
  names(a$beta) <- c("int",colnames(X))
  names(a$alpha) <- colnames(Y)
  #class(a) <- "rqcan"
  # change a$object :  a <- class("can") ?
  m <- length(a)
  a[[m+1]] <- X ; a[[m+2]] <- Y ; a[[m+3]] <- a.pos
  a[[m+4]] <- R ; a[[m+5]] <- r ; a[[m+6]] <- rho1
  names(a)[(m+1):(m+6)] <- c("X","Y","a.pos","R","r","rho")
  return(a)
  }

NULL

# Document boot.can
#' @export boot.can
#' @title Bootstrap
#' @description Internal function to carry out the bootstrap ; It is a sub-module of summary.rqcan; not intended for general use. \cr
#' The parameters may be passed by summary.rqcan (see below)
#' \cr
#' @details See help(summary.rqcan) ; If errors occur or modification is wanted, see the routine boot.can \cr
#' @param a output from rqcan \cr
#' @param Rep number of bootstrap replications (default=200) \cr
#' @param method "Andrews" (default) or ""xy" \cr
#' @param msub parameter defining the size of the bootstrap subsample for developmental work only: see the boot.can function \cr
#' @param seed a starting seed (default: missing: no new seed set) \cr
#' @param nsing number of consecutive singular replicatios to ignore (default = 5) \cr
#' @param prb if TRUE (default = FALSE), print every time 10 percent of the bootstrap samples are done \cr
#' @return Returns list(As, Bs, sdc): As (Bs) are N by dim(alpha) (dim(beta)) arrays of all bootstrap alphs and beta values; sdc = sqrt(m/n): SD adjustment for m-choose-n bootstrap, or 1 for "xy" bootstrap
#' \cr

boot.can <- function(a,Rep=200,method="Andrews",msub=.9,seed,nsing=5,prb=FALSE) {
  # returns As = alpha coef for Y, Bs = bet coef for X, sdc (SD = bootSD/sdc)
  if(!missing(seed)) set.seed(seed)
  X <- a$X ; Y <- a$Y ; tau <- a$tau ; a.pos <- a$a.pos ; ap = a$ap ; na <- a$na
  n <- nrow(X) ; px <- ncol(X) ; py <- ncol(Y)
  m <- ceiling( min(n, max(log(n)*(px+py+1), n^msub)) )
  As <- array(0,c(Rep,na,py)) ; Bs <- array(0,c(Rep,na,px+1))
  for(j in 1:Rep) {
    J <- 0
    b <- TRUE
    while( !is.list(b)) {
    if(method=="xy") { wts <- rep(1,n) ; I <- sample(1:n,replace=TRUE) ; sdc <- 1 }
    else if(method=="Andrews")  {
      I <- sample(n, m) ; wts <- -log(runif(m)) ; sdc <- sqrt(m/n) }
    else stop("Bootstrap method must be either Andrews or xy")
     b <- rqcan(X[I,],Y[I,],tau,a.pos,ap,na,wts)
     J <- J + 1  
     if(J > nsing) break }
    if(!is.list(b)) { 
      warning(paste("Bootstrap stopped at ",j," replications (",
        nsing," consecutive singularities occurred)",sep=""))
      return(list(As= As[1:j,,], Bs = Bs[1:j,,], sdc=sdc)) }
    As[j,,] <- b$alpha
    Bs[j,,] <- b$beta
    if(prb == TRUE & j %% (Rep/10) == 0) print(paste("Progress: #rep =",j)) }
  dimnames(Bs) <- list(NULL,paste(1:na),c("Intercept",colnames(X)))
  dimnames(As) <- list(NULL,paste(1:na),colnames(Y))
  return(list(As=As, Bs=Bs, sdc=sdc)) }


# Create a new S3 class
# new_rqcan <- function(x = list()) {
#   structure(x, class = "rqcan")
# }

# Document the rqcan class
#' @export
#' @name rqcan
#' @title My S3 Class
#' @description This is a simple S3 class for display formatting purposes.
#' @field list A list.

NULL

# Document the summary method
#' @export summary.rqcan
#' @title Summary of rqcan function results.
#' @description
#' Uses one of two bootstrap methods to provide Standard Error and confidence intervals for the alpha and beta coefficients for all components. \cr
#' \cr
#'
#' @param object rqcan object returned by rqcan
#' @param pr print tables if TRUE (default)
#' @param ci type of 95% confidence interval \cr
#' ci=1: use (adjusted) .025 and .975 percentiles of bootstrap distribution \cr
#' ci=2: use normal approx with adjusted SE based on interquartile range \cr
#' ci may be a vector to provide more than one type of interval; default=1 \cr
#' @param fact a factor to adjust conf ints for components 2:na; used only for development \cr
#' @param ... parameters that are sent to the bootstrap function: \cr
#' Rep: number of bootstrap replications (default=200) \cr
#' method: "Andrews" (default) or ""xy" \cr
#' msub: parameter defining the size of the bootstrap subsample for developmental work only: see the boot.can function \cr
#' seed: a starting seed (default: missing: no new seed set) \cr
#' nsing: number of consecutive singular replicatios to ignore (default = 5) \cr
#' prb: if TRUE (default = FALSE), print every time 10 percent of the bootstrap samples are done \cr
#' @details
#' The Portnoy reference showed that a subsample bootstrap (as described by Andrews) gives consistent estimates of SE's and confidence intervals. The subsample size is m = ceiling( min(n, max(log(n)*(px+py+1), n^msub)) ) (where n = nrow(X), px = ncol(X), py = ncol(Y)), msub is as above). Some simulations and examples suggest that this is OK. The usual "xy" bootstrap (sampling rows independently with replacement) can be specified. It seems to give similar confidence intervals to "Andrews", but the SE estimates may be wrong; and no form of consistency has been proven. Note: as noted in help(rqcan), the quantreg function rq.fit.fnc may generate singular matrices even if the input design matrix is of full rank. In simulation examples, this can happen for some bootstrap replications (perhaps less than 1/1000 times). When this occurs, a new bootstrap replication is drawn. If more than nsing consecutive singularities are produced, the bootstrap function returns with those replications that it has already found (a number less than Rep), with a warning. If a singularity warning occurs, using "xy", or changing the seed or "jittering" the data (see jitter()) sometimes helps. \cr
#'
#' @return Returns list(As,Bs,sdc): As and Bs are matrices with Rep rows giving alpha beta coefficients for each bootstrap replication; and sdc is a standard error adjustment based on the subsample bootstrap:   sdc = sqrt(1 - m/n). \cr
#' \cr
summary.rqcan <- function(object, pr = TRUE, ci = 1, fact = 1, ...)  {
  # object = result from rqcan(X,Y,tau,a.pos,na) ; wts should be the default: rep(1,n)
  # print if pr == TRUE
  # ci: vector of one or more ci types: 1 = %-tile, 2 = adj se from IQR
  #     3 = %-tile using norm-adj quartiles; deleted
  # ... = input params for boot.can if not default: Rep,method,seed,np,na,msub,pr
  a <- object
  cis <- ci
  y <- a$Y ; x <- a$X ; q <- ncol(y)
  tau <- a$tau ; na <- a$na
  vnames <- c(paste("alpha:",colnames(a$alpha),sep=""),
              paste("beta:",colnames(a$beta),sep="") )
  p <- ncol(x) + ncol(y) + 1
  rdf <- nrow(x) - p
  coef <- vector("list", length=na)
  names(coef) <- paste("RQcan",1:na)
  B <- boot.can(a, ...)
# p1 <- function(v) c(v[2] + (v[3]-v[2])*2.90585,v[2] - (v[2]-v[1])*2.90585)
  sdc <- B$sdc
  fac <- c(1,rep(fact,na-1))
  plim <- c(.025,rep(.025/fact,na-1))
  # form coef array for each can vector
  for(j in 1:na) {
    serr <- sqrt(c(diag(cov(B$As[,j,])),fac[j] * diag(cov(B$Bs[,j,])) ))
    est <- c(a$alpha[j,],a$beta[j,])
    nci <- length(ci)
    cs <- matrix(0, p, 2+2*nci)
    dimnames(cs) <- list(vnames, c("Value", "Std. Error",
                      paste("ci",c(rbind(ci,ci)),rep(c(":low", ":up"),nci),sep="")))
    cs[,1] <- est
    cs[,2] <- serr
    # compute ci's for all ci-types
    for(k in 1:nci) 
      if(ci[k] != 1 & ci[k] != 2) { mess <- paste("ci = conf int type must be 1 or 2; if not, ourput reflects ci[",k,"] = 1",sep="")
        warning(mess)
        cis[k] <- 1 }
      if(cis[k] == 1) {
        bd <- apply(B$As[,j,],2,quantile,prob=c(.025,.975),type=4)
        bd <- cbind(bd, apply(B$Bs[,j,],2,quantile,prob=c(plim[j],1-plim[j]),type=4))
        low <- est*(1 + sdc) - sdc*bd[2,]   #  old: 2*est - bd[2,]
        up <-  est*(1 + sdc) - sdc*bd[1,]   #  old: 2*est - bd[1,]
        cs[,1+2*k] <- low
        cs[,2+2*k] <- up }
      else if(cis[k] == 2) {
        bd <- apply(B$As[,j,],2,IQR,type=4)
        bd <- c(bd, fac[j] * apply(B$Bs[,j,],2,IQR,type=4) ) / 1.349
        low <- est - 2*bd
        up <-  est + 2*bd
        cs[,1+2*k] <- low
        cs[,2+2*k] <- up }
    coef[[j]] <- cs }
  if(pr) { options(digits = 4) ; print(coef) }
  return(invisible( list(coef=coef,parnames=vnames) ))
}
