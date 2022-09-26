# contact_and_NegVE

This repository is for the paper, "Higher contact among vaccinated can be a mechanism for negative vaccine effectiveness", currently available here: https://www.medrxiv.org/content/10.1101/2022.04.25.22274266v2 with a plain-language summary of the paper found here: https://mishra-lab.ca/2022/05/12/higher-contact-among-vaccinated-can-be-a-mechanism-for-negative-vaccine-effectiveness/

The R scripts, "contact_and_NegVE_SEIR.R" and "contact_and_NegVE_SIR.R" are the two main coding scripts. These scripts perform the following two functions:

1) run the main simulations with an SEIR and SIR model, respectively, and;
2) generate all figures in the manusript (including those in the supplementary material). 

The main R script, "contact_and_NegVE_SEIR.R", uses "VES_and_ContactIncrease_SEIR.csv" and "VES_and_VEI_SEIR.csv" to create the main contour plots (Figure 2a and b) and "VES_and_ContactIncrease_SIR.csv" and "VES_and_VEI_SIR.csv" as well as "VEI_and_ContactIncrease_Under_VES0_*.csv" to produce contour plots for the supplementary material. 

Data for the contour plots (the .csv files listed above) are generated by the R scripts,"contour_data_generation_SEIR.R" and "contour_data_generation_SIR.R".
