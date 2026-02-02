#### Set up directories containing auxiliary functions (PATH) and data (datadir)

PATH <- "C:"
datadir <- "C:"

# Input smoothing and empirical coverage functions
source(paste(PATH,"aux_functions.R",sep="//"))

# Input functions to estimate parameters
source(paste(PATH,"MLE.R",sep="//"))
source(paste(PATH,"MLE_poll_effect.R",sep="//"))

# Read and process Pew data

mercer <- data.frame(read.csv(paste(datadir,"pew_benchmarking_data.csv",sep="//"),header=T))

### If phat \approx 0 or 1, then use Agresti-Coull to modify point estimate and calculate vi

modify_p <- mercer$phat > 0.95 | mercer$phat < 0.05

# Check min(n*p*(1-p)) for estimates not subject to AC modification
min(mercer$phat[!modify_p]*(1-mercer$phat[!modify_p])*mercer$n[!modify_p])

mercer$phat[modify_p] <- (mercer$phat[modify_p]*mercer$n[modify_p] + 2)/(mercer$n[modify_p] + 4)
mercer$vi <- mercer$phat*(1-mercer$phat)/mercer$n
mercer$vi[modify_p] <- mercer$phat[modify_p]*(1-mercer$phat[modify_p])/(mercer$n[modify_p] + 4)

mercer$phat_minus_p0 <- mercer$phat - mercer$p0
mercer$emse <- mercer$phat_minus_p0^2
mercer$emp_biassq_vi <- mercer$emse - mercer$vi - mercer$v0
mercer$vi_plus_v0 <- mercer$vi + mercer$v0

mercer$topic <- rep("Other",nrow(mercer))
mercer$topic[mercer$benchmark_category %in% c("Received food stamps","Received Social Security","Received unemployment compensation","Received worker's compensation")] <- "Received government benefits"
mercer$topic[mercer$benchmark_category %in% c("Never smoked 100 cigarettes","Never vaped")] <- "Never smoked, never vaped"

mercer$type_topic <- paste(mercer$type,mercer$topic,sep=" - ")

mercer$type <- factor(mercer$type,levels=c("prob","nonprob"))
mercer$topic <- factor(mercer$topic,levels=c("Other","Never smoked, never vaped","Received government benefits"))

mercer$type_topic <- factor(mercer$type_topic,levels=c("prob - Other","prob - Never smoked, never vaped","prob - Received government benefits","nonprob - Other","nonprob - Never smoked, never vaped","nonprob - Received government benefits")) 

mercer$subgroup_category <- factor(mercer$subgroup_category,levels=c("All adults","Hispanic","Black non-Hispanic","White non-Hispanic","18-29","30-64","65+","HS or less","Some college","College grad")  )

mercer$sampling_error_coverage_ind <- emp_coverage(mercer$phat_minus_p0,estvar=mercer$vi+mercer$v0)

# Write datafile to csv for use with other software
write.csv(mercer,paste(datadir,"mercer.csv",sep="//"),row.names=F)
                 

############### Calculate coverage probability, overall and by type and topic #########################

group_mean(sampling_error_coverage_ind,type,data=mercer)
group_mean(sampling_error_coverage_ind,response_rate,data=mercer)
group_mean(sampling_error_coverage_ind,topic,data=mercer)
group_mean(sampling_error_coverage_ind,type_topic,data=mercer)


############### MULTIPLICATIVE MOTES #############################################

# One factor

multfactor_all <- multiplicative_model(fitdata=mercer,y_var="phat_minus_p0",sampling_var="vi",benchmark_var="v0")  

multfactor_all
sqrt(multfactor_all)


# By type

multfactor_type <- lapply(split(mercer,mercer$type),multiplicative_model,
                   y_var="phat_minus_p0",sampling_var="vi",benchmark_var="v0")
multfactor_type
lapply(multfactor_type,sqrt)

# By topic

multfactor_topic <- lapply(split(mercer,mercer$topic),multiplicative_model,
                   y_var="phat_minus_p0",sampling_var="vi",benchmark_var="v0")
multfactor_topic
lapply(multfactor_topic,sqrt)

# By type and topic


multfactor_type_topic <- lapply(split(mercer,mercer$type_topic),multiplicative_model,
                   y_var="phat_minus_p0",sampling_var="vi",benchmark_var="v0")
multfactor_type_topic
lapply(multfactor_type_topic,sqrt)


# Look at coverage by sample size

mercer$multmote_coverage <- emp_coverage(mercer$phat_minus_p0,estvar=mercer$vi*multfactor_all+mercer$v0)
mean(mercer$multmote_coverage)
mean(mercer$multmote_coverage[mercer$n <= 2000])
mean(mercer$multmote_coverage[mercer$n >= 4000])


############# CALCULATE COVERAGE IF MULTIPLY SE BY 2 ##################################

mercer$factor2_coverage_ind <- emp_coverage(mercer$phat_minus_p0,estvar=mercer$vi*4+mercer$v0)


group_mean(factor2_coverage_ind,type,data=mercer)
group_mean(factor2_coverage_ind,topic,data=mercer)
group_mean(factor2_coverage_ind,type_topic,data=mercer)


############# LINEAR REGRESSION MODEL (SECTION 2.3) #########################################

# See SAS code for PROC MIXED in Pew_linear_model.sas
# Can also fit using lme4 package

################# LINEAR MIXED MODEL ###########################################


ML_fit <-G_est(phat_minus_p0~type*topic,sampling_var="vi_plus_v0",is.nullX = FALSE,REformula = phat_minus_p0~type-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=mercer)

ML_fit$fit


ML_fit_dep <-G_est_dep(phat_minus_p0~type*topic,sampling_var="vi_plus_v0",poll_ID = "sampid",is.nullX = FALSE,REformula = phat_minus_p0~type-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=mercer)

ML_fit_dep$fit

# Look at expected squared bias for each type/topic combination

tapply(ML_fit_dep$estsqbias,mercer$type_topic1,mean)


ML_fit_type =G_est_dep(phat_minus_p0~type-1,sampling_var="vi_plus_v0",poll_ID = "sampid",is.nullX = FALSE,REformula = phat_minus_p0~type-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=mercer)


ML_fit_type$fit


ML_fit_topic =G_est_dep(phat_minus_p0~topic,sampling_var="vi_plus_v0",poll_ID = "sampid",is.nullX = FALSE,REformula = phat_minus_p0~type-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=mercer)


ML_fit_topic$fit

### What happens if we fit separate variance inflations for the types and topics, with no fixed effects?

ML_fit_G_only <-G_est_dep(phat_minus_p0~1,sampling_var="vi_plus_v0",poll_ID = "sampid",is.nullX = TRUE,REformula = phat_minus_p0~type_topic-1,Gfcn = Gfcn_1,partialG = partialG_1, PRECISION = 1e-10, MAXITER = 40,data=mercer)

ML_fit_G_only$fit

tapply(ML_fit_G_only$estsqbias,mercer$type_topic,mean)







