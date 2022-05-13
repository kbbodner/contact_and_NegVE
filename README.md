# contact_and_NegVE

This repository is for the paper, "Higher contact among vaccinated can be a mechanism for negative vaccine effectiveness", currently available here: https://www.medrxiv.org/content/10.1101/2022.04.25.22274266v2 with a plain-language summary of the paper found here: https://mishra-lab.ca/2022/05/12/higher-contact-among-vaccinated-can-be-a-mechanism-for-negative-vaccine-effectiveness/

The two R scripts, "contact_and_NegVE.R" and "contour_data_generation.R" contain code to recreate all simulations and to produce Figure1 (a,b,c and d) in the main document and Figure S2 in the supplementary material.

The R script, "contact_and_NegVE.R", contains code to run the simulations for Figure 1a and b and to generate all figures in the manusript. It calls in "VES_and_ContactIncrease.csv" and "VES_and_VEI.csv", which contain the scenario data used to produce the contour plots (Figure 1c and d). "VES_and_ContactIncrease.csv" and "VES_and_VEI.csv" were generated with code from "contour_data_generation.R".


