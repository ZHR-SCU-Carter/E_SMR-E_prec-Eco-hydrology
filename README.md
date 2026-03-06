# E_SMR-E_prec-Eco-hydrology
This is the reproduction code of the results of the article "freeze thaw driven soil moisture return significant contributions to spring physiology on the warming Qinghai Tibet Plateau", including the processing, drawing, analysis and other steps of the core results.

This research code is mainly divided into MATLAB code. We used MATLAB for data processing, analysis of θcri, implementation of sliding segmentation, modeling of non-linear mixed effects model and quantify the contribution rate of meteorological and hydrological element variables to soil moisture content at different depths in the main research.


                      The overall workflow for the MATLAB R2022b scripts used in this study.

                  Step 1: Identify the site-specific θ_cri (critical soil moisture threshold)
                     using the procedures provided in "1 - theta_cri_determination.m".
 
                Step 2: Visualize the soil moisture data and perform manual outlier checks using
                          "2 - extractAndVisualizeSMData.m". This script also highlights
                                       periods when soil moisture exceeds θ_cri.
 
                Step 3: Based on the outcomes from Step 2, apply "3 - analyzeESmrprecEffect.m" 
                    to perform a segmented analysis of spring soil moisture enhancement,
                                    and quantify both E_SMR and E_prec.
 
                Step 4: Evaluate the impacts of E_SMR, E_prec, and near-surface air temperature (airT) 
                    on the start of season (SOS) using the nonlinear mixed-effects model, 
                                and quantify their relative contributions.
          
                       We uploaded a Demon.mat database to test the usability of the code

                                  
