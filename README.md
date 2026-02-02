# MOTE
R functions for estimating a margin of total error using benchmark poll data

This repository has programs (in R and SAS&reg; software) and data used to calculate estimates in

Lohr, SL, Mercer, A, Kennedy, C, and Brick, JM. (2026). Toward a Margin of Total Error. Journal of Survey Statistics and Methodology.

The MIT license applies to these programs (see LICENSE.txt). You are free to adapt them provided you give attribution to the original source. Please cite the above paper when you use the software.

SAS&reg; and all other SAS Institute Inc. product or service names are registered trademarks or trademarks of SAS Institute Inc., Cary, NC, USA. &reg; indicates USA registration.


**R functions called in analysis programs**

|   Name    | Description  |
|--------------------------|-----------------------------------------------------------------------------------------------------|
|   MLE.R		| estimate parameters for the linear mixed model |
|   MLE_poll_effect.R	| estimate parameters for the linear mixed model with extra random effect for poll (model in Appendix A) |
|   aux_functions.R	| auxiliary functions used to calculate estimates |

**Analysis of Shirani-Mehr (2018) data**

The data are described in

Shirani-Mehr, H, Rothschild, D, Goel, S, and Gelman, A (2018). Disentangling Bias and Variance in Election Polls. Journal of the American Statistical Association, 113, 607–614. 
Data are available at https://github.com/5harad/polling-errors.

|   Name    | Description  |
|--------------------------|------------------------------------------------------------------|
| input_Shirani_Mehr_data.R   | code to read and process data |
| analyze_Shirani_Mehr_data.R | code to analyze data |

**Analysis of Pew data**

The Pew data are described in

Mercer, A and Lau, A (2023). Comparing Two Types of Online Survey Samples. Pew Research, https://www.pewresearch.org/methods/2023/09/07/comparing-two-typesof-online-survey-samples/. 
Data are available in file "Pew_benchmarking_data.csv" in the Data folder of this directory. The questionnaire for the survey is at https://www.pewresearch.org/methods/wp-content/uploads/sites/10/2023/09/pm_09.07.23_benchmarking_questionnaire.pdf


|   Name    | Description  |
|--------------------------|------------------------------------------------------------------|
| analyze_Pew_data.R  | R code to analyze the Pew data |
| Pew_linear_model.sas | SAS&reg; software code to calculate the coefficients in Table 3.7 |
