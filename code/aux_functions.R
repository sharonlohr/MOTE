# Auxiliary functions for analyzing benchmark poll data

group_mean <- function(varname,groupname,data,remove_NA=TRUE) {
  vname <- substitute(varname)
  gname <- substitute(groupname)
  mean_all <- mean(data[[vname]],na.rm=remove_NA)
  mean_group <- tapply(data[[vname]],data[[gname]],mean,na.rm=remove_NA)
  mean_group <- c(mean_all,mean_group)
  names(mean_group)[1] <- "All"
  return(mean_group)
}

emp_coverage <- function(y,trueval = 0,estvar,alpha=.05) {
  # calculates vector of indicators, = 1 if trueval is in y +/- z_critical*sqrt(estvar), 0 otherwise
  # estvar is the estimated total variance

  m <- length(y)
  if (length(estvar) != m) {return("length of y unequal to length of estvar")}
  MOTE <- qnorm(1-alpha/2)*sqrt(estvar)
  coverage <- rep(0,m)
  coverage[y-MOTE < trueval & y+MOTE > trueval] <- 1
  return(coverage)
}

logdet <- function(vec,multiplier = 1) {
   logdeterm <- sum(log(vec)) + log(1 + multiplier*sum(1/vec))
   return(logdeterm)
}


multiplicative_model <- function(fitdata,y_var,sampling_var,benchmark_var,tol=1e-25,ub_interval=400) {

 # Uses method of moments to estimate multiplicative factor (multmodel$root) for inflating sampling variance
 # Finds the root f giving \sum( (\hat{p} - \hat{p0})^2 )/ (sampling_var * f + benchmark_var) ) = m
 # All variable names should be in quotes
 # y_var = name of variable containing \hat{p} - \hat{p0}
 # sampling_var = variance of \hat{p}
 # benchmark_var = variance of \hat{p0}, should be set to 0 if benchmark has no sampling error
 # tol = convergence tolerance 
 # ub_interval = upper bound of interval in which to search for root 
   multmodel <- uniroot(f = function(x,y,sv,bv) {mean(y^2/(sv*x+bv))-1}, 
          interval = c(0,ub_interval),extendInt="yes",tol=tol,y=fitdata[[y_var]],sv=fitdata[[sampling_var]],
          bv=fitdata[[benchmark_var]])
   return(multmodel$root)
}


Gfcn_1 <- function(delta) {
   l <- length(delta)
   G <- matrix(0,ncol=l,nrow=l)
   diag(G) <- delta
   return(G)
}


partialG_1 <- function(delta) {
   #returns array of partial derivatives of G with respect to components of delta
   # Partial derivatives for delta_l are in G[l,,]
   # This function is for categorical variable with l categories
   l <- length(delta)
   dG <- array(0,dim=c(l,l,l))
   for (j in 1:l) {dG[j,j,j] <- 1}
   return(dG)
}


