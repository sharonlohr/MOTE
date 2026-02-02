# Read and preprocess Shirani-Mehr data

# Data source: Shirani-Mehr, H., Rothschild, D., Goel, S., and Gelman, A. (2018) 
# “Disentangling Bias and Variance in Election Polls,” Journal of the American Statistical Association, 113, 607–614.
# Dataset polls_main_dataset.tsv can be downloaded from https://github.com/5harad/polling-errors

# SET PATH AND DATADIR
# PATH is directory containing auxiliary functions; datadir is directory containing the data
PATH <- "C:"
datadir <- "C:" 

source(paste(PATH,"aux_functions.R",sep="//"))

SMdata <- read.delim(paste(datadir,"polls_main_dataset.tsv",sep="//"))

SMdata$poll_prop_rep <- SMdata$republican/(SMdata$republican+SMdata$democratic)  # p-hat
SMdata$ni <- round((SMdata$republican+SMdata$democratic)*SMdata$numberOfRespondents/100) # sample size
SMdata$final_prop_rep <- SMdata$finalTwoPartyVSRepublican/100 # p_0
SMdata$poll_minus_vote <- SMdata$poll_prop_rep - SMdata$final_prop_rep
SMdata$emse <- SMdata$poll_minus_vote^2
SMdata$vi <- SMdata$poll_prop_rep*(1-SMdata$poll_prop_rep)/SMdata$ni  # vi-hat based on p-hat
SMdata$vi_p0 <- SMdata$final_prop_rep*(1-SMdata$final_prop_rep)/SMdata$ni # vi-hat based on p0
SMdata$vi0 <- 0   # Variance of benchmark proportions is 0 because known exactly
SMdata$emp_biassq_vi <- SMdata$emse - SMdata$vi
SMdata$emp_biassq_vi_p0 <- SMdata$emse - SMdata$vi_p0
SMdata$p0_1minusp0 <- SMdata$final_prop_rep*(1-SMdata$final_prop_rep)
SMdata$election_cat <- paste(SMdata$election,SMdata$year,SMdata$state)
SMdata$days_to_election <- as.integer(as.Date(SMdata$electionDate)-as.Date(SMdata$endDate))
SMdata$sampling_error_coverage_ind <- emp_coverage(SMdata$poll_minus_vote,estvar=SMdata$vi)
SMdata$before2007 <- (SMdata$year <= 2006)
SMdata$U_i <- SMdata$poll_minus_vote^2/(SMdata$vi)

write.csv(SMdata,paste(datadir,"SMdata.csv",sep="//"))




