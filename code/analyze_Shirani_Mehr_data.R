# Compute estimates for Shirani-Mehr data

####### PRELIMINARIES TO LOAD DATA AND FUNCTIONS ########

# SET PATH # PATH is directory containing auxiliary functions

PATH <- "C:"

# Input auxiliary functions
source(paste(PATH,"aux_functions.R",sep="//"))

# Input functions to estimate parameters

source(paste(PATH,"MLE.R",sep="//"))
source(paste(PATH,"MLE_poll_effect.R",sep="//"))

# Input and pre-process data

source(paste(PATH,"input_Shirani_Mehr_data.R",sep="//"))

names(SMdata)

# Create subset data for elections before 2007, after 2007

SMdata1 <- SMdata[SMdata$before2007,]
SMdata2 <- SMdata[!SMdata$before2007,]


######## ESTIMATE EMSE, COVERAGE WITH SAMPLING VARIABILITY #######

### Table 3.1

# Average EMSE

group_mean(emse,election,SMdata)

# Average variance v-bar

group_mean(vi,election,SMdata)

# Average squared bias b-bar

group_mean(emp_biassq_vi,election,SMdata)

# Empirical coverage with sampling variability

group_mean(sampling_error_coverage_ind,election,SMdata)

# Empirical coverage for sampling variability with 2008-2014 data, Table 3.4

group_mean(sampling_error_coverage_ind,election,SMdata2)

######## ESTIMATE EMSE, COVERAGE WHEN MULTIPLY SE BY 2, ADD BBAR #######

## Coverage if multiply standard error by 2

SMdata$mult2_coverage_ind <- emp_coverage(SMdata$poll_minus_vote,estvar=SMdata$vi*4)
group_mean(mult2_coverage_ind,election,SMdata)

# Data after 2007

group_mean(mult2_coverage_ind,election,SMdata[!SMdata$before2007,])

## Coverage if add bbar

SMdata$add_bbar_coverage_ind <- emp_coverage(SMdata$poll_minus_vote,estvar=SMdata$vi+mean(SMdata$emp_biassq_vi))
group_mean(add_bbar_coverage_ind,election,SMdata)

# Data after 2007

SMdata2$add_bbar_coverage_ind <- emp_coverage(SMdata2$poll_minus_vote,estvar=SMdata2$vi+mean(SMdata1$emp_biassq_vi))
group_mean(add_bbar_coverage_ind,election,SMdata2)


########################### MULTIPLICATIVE ADJUSTMENTS #################################

mult_factor_variance <- group_mean(U_i,election,SMdata)
mult_factor_variance
sqrt(mult_factor_variance)

# Evaluate coverage if one factor applied to all estimates

SMdata$cov_multmodel_one_factor <- emp_coverage(SMdata$poll_minus_vote, trueval=0,estvar=SMdata$vi*mult_factor_variance[1])
group_mean(cov_multmodel_one_factor,election,SMdata)

# Estimate coverage if separate factors by election type

SMdata$multfactor_by_type <- mult_factor_variance[match(SMdata$election, names(mult_factor_variance))]

SMdata$cov_multmodel_mult_factor <- emp_coverage(SMdata$poll_minus_vote, trueval=0,estvar=SMdata$vi*SMdata$multfactor_by_type)
group_mean(cov_multmodel_mult_factor,election,SMdata)

###### Use factors obtained from 1998 to 2006 data on data from 2007 to 2014

mult_factor_variance1 <- group_mean(U_i,election,SMdata1)
mult_factor_variance1

# Coverage, one factor 

SMdata2$cov_multmodel_one_factor <- emp_coverage(SMdata2$poll_minus_vote, trueval=0,estvar=SMdata2$vi*mult_factor_variance1[1])
group_mean(cov_multmodel_one_factor,election,SMdata2)

# Coverage, separate multiplicative factors by election type

SMdata2$multfactor_by_type <- mult_factor_variance1[match(SMdata2$election, names(mult_factor_variance1))]

SMdata2$cov_multmodel_mult_factor <- emp_coverage(SMdata2$poll_minus_vote, trueval=0,estvar=SMdata2$vi*SMdata2$multfactor_by_type)
group_mean(cov_multmodel_mult_factor,election,SMdata2)

##################################################################################################

####### Fit regression models of empirical squared bias as function of covariates.

linmodel1 <- lm(emp_biassq_vi ~ election + p0_1minusp0 + days_to_election-1,data=SMdata)
summary(linmodel1)
round(linmodel1$coef,6)
plot(linmodel1$fitted.values,linmodel1$residuals)

# Calculate coverage

mean(linmodel1$fitted.values < 0 )
linmodel1$gu <- linmodel1$fitted.values
linmodel1$gu[linmodel1$gu < 0] <- 0
linmodel1$emp_coverage <- emp_coverage(SMdata$poll_minus_vote,estvar = SMdata$vi + linmodel1$gu)

mean(linmodel1$emp_coverage)
tapply(linmodel1$emp_coverage,SMdata$election,mean)


#### Try with log tranformation, does not help

model1_log <- lm(log(emp_biassq_vi + .005)~ election + p0_1minusp0 + days_to_election-1,data=SMdata)
summary(model1_log)

plot(model1_log$fitted.values,model1_log$residuals)

# Calculate predicted squared bias in original scale

model1_log$fitted.values.orig.scale <- exp(model1_log$fitted.values) - .005
mean(model1_log$fitted.values.orig.scale < 0 ) 

##### Fit to before2007, evaluate coverage for rest of data

model1s <- lm(emp_biassq_vi ~ election + p0_1minusp0 + days_to_election-1,data=SMdata1)
summary(model1s)

# Calculate coverage
model1s$gu <- predict(model1s,SMdata2)
model1s$gu[model1s$gu < 0] <- 0
model1s$emp_coverage <- emp_coverage(SMdata2$poll_minus_vote,estvar = SMdata2$vi + model1s$gu)
mean(model1s$emp_coverage)
tapply(model1s$emp_coverage,SMdata2$election,mean)


##################### LINEAR MIXED MODELS #########################


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
   # Or for any setup where gammas are independent
   l <- length(delta)
   dG <- array(0,dim=c(l,l,l))
   for (j in 1:l) {dG[j,j,j] <- 1}
   return(dG)
}

#################################FIT MODELS################################################

## Model (3.1) with one G for all data, no fixed effects

ML_fit_no_intercept =G_est(formula=poll_minus_vote~1,sampling_var="vi",is.nullX = TRUE,REformula = poll_minus_vote~1,,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=SMdata)

ML_fit_no_intercept$fit

# Evaluate coverage

mean(ML_fit_no_intercept$coverage_ind)
tapply(ML_fit_no_intercept$coverage_ind,SMdata[["election"]],mean)

# Fit with SMdata1, eval with SMdata2

ML_fit_no_intercepts =G_est(formula=poll_minus_vote~1,sampling_var="vi",is.nullX = TRUE,REformula = poll_minus_vote~1,,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=SMdata1)

ML_fit_no_intercepts$fit

cov2 <- emp_coverage(SMdata2$poll_minus_vote,estvar = SMdata2$vi + as.vector(ML_fit_no_intercepts$fit$estG))
mean(cov2)
tapply(cov2,SMdata2[["election"]],mean)



##### Model (3.2) ###############################

ML_fit_days_type =G_est(poll_minus_vote~days_to_election,sampling_var="vi",is.nullX = FALSE,REformula=poll_minus_vote~election - 1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,method="ML",data=SMdata)

ML_fit_days_type$fit


mean(ML_fit_days_type$coverage_ind)
tapply(ML_fit_days_type$coverage_ind,SMdata[["election"]],mean)


# Fit with SMdata1, eval with SMdata2

ML_fit_days_type_s =G_est(poll_minus_vote~days_to_election,sampling_var="vi",is.nullX = FALSE,REformula=poll_minus_vote~election - 1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=SMdata1)

ML_fit_days_type_s$fit

estsqbias_fixed <- (ML_fit_days_type_s$fit$estcoef[1,1] + ML_fit_days_type_s$fit$estcoef[2,1]*SMdata2[["days_to_election"]])^2
estsqbias_random <- rep(ML_fit_days_type_s$fit$estdelta[1],nrow(SMdata2))
estsqbias_random[SMdata2$election=="Pres"] <- ML_fit_days_type_s$fit$estdelta[2]
estsqbias_random[SMdata2$election=="Sen"] <- ML_fit_days_type_s$fit$estdelta[3]

cov2 <- emp_coverage(SMdata2$poll_minus_vote,estvar = SMdata2$vi +  estsqbias_fixed + estsqbias_random)
mean(cov2)
tapply(cov2,SMdata2[["election"]],mean)

###############################################################

# Model (3.3) scaling each G by  z = (pi0*(1-pi0))


ML_fit_no_intercept_zscale =G_est(poll_minus_vote~1,sampling_var="vi",is.nullX = TRUE,REformula=poll_minus_vote~p0_1minusp0-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,method="ML",data=SMdata)

ML_fit_no_intercept_zscale$fit

Model33_cov <- emp_coverage(SMdata$poll_minus_vote,estvar = SMdata$vi +  (SMdata$poll_prop_rep*(1-SMdata$poll_prop_rep))^2*ML_fit_no_intercept_zscale$fit$estdelta[1])

mean(Model33_cov)
tapply(Model33_cov,SMdata[["election"]],mean)


# Fit with SMdata1, eval with SMdata2

ML_fit_no_intercept_zscale_s =G_est(poll_minus_vote~1,sampling_var="vi",is.nullX = TRUE,REformula=poll_minus_vote~p0_1minusp0-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=SMdata1)

ML_fit_no_intercept_zscale_s$fit

estsqbias_random <- (SMdata2$poll_prop_rep*(1-SMdata2$poll_prop_rep))^2*ML_fit_no_intercept_zscale_s$fit$estdelta[1]
cov2 <- emp_coverage(SMdata2$poll_minus_vote,estvar = SMdata2$vi +   estsqbias_random)
mean(cov2)
tapply(cov2,SMdata2[["election"]],mean)










