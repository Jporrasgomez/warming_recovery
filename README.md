# warming_recovery


R-code for the work "Warming alters post-disturbance temporal dynamics of grassland plant communities"


All data required to run this code can be found in: public repository will be indicated upon acceptance of the manuscript


Empty folders (data, data>processed_data and code) are included to enhance code running. 


The code is structured with sections from 0 to 13. They are intended to be run subsequently, 
minding some exceptions. 



- 0.OTC_effect.R (0.R) - Contains the code to open, process and visualize the effect of open top chambers (OTC) on air temperature. 
  Data was recorded using TMS-4 sensors (Tomst (R)). 
  
- 1.0.processing_data.R (1.R) - Contains the code to open and process data compiled in the field experiment. 
  Several subsections are found in this code: 
    -  1.1.0.biomass_imputation_level_MICE.R (1.1.R) - The code used to generate the data opened in line 238 of 1.R. This code requires computation time. 
    -  1.1.1.loop_reliability. R - Code to check imputation reliability of mice in 1.1.0.
      
- 1.2.biomass_imputation_level_LM.R (1.2.R) - Contains the imputation level of samplings 0, 1, 2 and 12 using linear regression models. 

- 2.Evenness_RAD_models.R (2.R) - Code for Rank Abundance Distribution models (i.e. Evenness analysis). 

- 3.0.species_composition_NMDSbray.R (3.R) - Code for species composition analysis and visualization (Sorensen, NMDS, PERMANOVA).

- 3.1.species_composition_new_species.R (3.1.R) - Code focused on sub-analysis of colonizers in perturbation and combined treatments.

- 4.functional_traits_LES.R (4.R) - Code for functional traits data analysis, CWM estimation and visualization. 

- 5.merging_database.R (5.R) - Processed data generated in 0.R, 1.R, 1.2.R, 2.R and 4.R is merged into a single database for effect size analysis

- 6.choosing_GLMMS (6.R) - Building and optimizing Generalized Linear Mixed Models for each variable using the merged database created in 5.R

- 7.running_GLMMs (7.R) - Running the chosen GLMMs in 6.R and storing the results.

- 8.0.Log_Response_Ratio_analysis - Code for the effect size using log response ratio (LRR)
  - 8.1.effect_size_metafor_checking - Checking the reliability of our LRR code using a specialized package for effect size analysis

- 9.LRR_vs_GLMM_results_visualization - Some plots to visualize the consistency between GLMM and LRR results

- 10.Main_Results_Visualization - Joininf GLMM and LRR results to be plot for main Figures (Figures 2 and 4). Here we use the functions stored in the folder 
  "functions". 

- 11.Sampling_0 - Analysis for sampling 0 (The sampling before perturbation were applied)

- 12.Biomass~Richness - Code for the biomass ~ richness linear fitted relationships 
      
- 13.0.sensitivity_analysis_Z_biomass.R - Code to perform sensitivity analysis of Z coefficient (2/3) used in allometric equation of 
  biomass and its interaction with imputation levels in 1.1.R and 1.2.R.
  - 13.1.GLMM_biomass_Z_sensitivity - Choosing the best GLMMs for biomass iterations. 

- palettes_labels.R - Code where objects for visualization colors and labels are created to ensure visual consistency.















