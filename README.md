# warming_recovery
R-code for the work "Experimental warming delays ecosystem recovery in grasslands"




All data required to run this code can be found in: public repository will be indicated upon acceptance of the manuscript


Empty folders (data, data>processed_data and code) are included to enhance code running. 


The code is structured with sections from 0 to 9. They are intended to be run subsequently, 
minding some exceptions. 



- 0.OTC_effect.R (0.R) - Contains the code to open, process and visualize the effect of open top chambers (OTC) on air temperature. 
  Data was recorded using TMS-4 sensors (Tomst (R)). 
  
  
- 1.0.processing_data.R (1.R) - Contains the code to open and process data compiled in the field experiment. 
  Several subsections are found in this code: 
      · 1.1.0.biomass_imputation_level_MICE.R (1.1.R) - The code used to generate the data opened in line 238 of 1.R. This code requires computation time. 
      · 1.1.1.loop_reliability. R - Code to check imputation reliability of mice in 1.1.0.
    
      
- 1.2.biomass_imputation_level_LM.R (1.2.R) - Contains the imputation level of samplings 0, 1, 2 and 12 using linear regression models. 


- 2.RADmodel.R (2.R) - Code for Rank Abundance Distribution models (i.e. Evenness analysis). 


- 3.0.species_composition_NMDSbray.R (3.R) - Code for species composition analysis and visualization (Sorensen, NMDS, PERMANOVA).


- 3.1.species_composition_new_species.R (3.1.R) - Code focused on sub-analysis of colonizers in perturbation and combined treatments.


- 4.functional_traits_LES.R (4.R) - Code for functional traits data analysis, CWM estimation and visualization. 


- 5.merging_database.R (5.R) - Processed data generated in 0.R, 1.R, 1.2.R, 2.R and 4.R is merged into a single database for effect size analysis


- 6.EFFECT_SIZE.R (6.R) - Here we perform Geary test and the Log Response Ratio (LRR) analysis and visualization of main figures of the manuscript
  (Figures 2 and 3) as well as Supplementary Figures 10 and 11. To run this code we need the source the functions contained in the folder
  "functions". Here we find: 
      · eff_size_LRR_function.R - function needed to perform LRR analysis at aggregated level
      · eff_size_dynamics_LRR_function.R - function needed to perform LRR analysis of temporal series
      · gg_aggregated_function.R - function needed to generate plots at aggregated level (left panel of Figures 2 and 3)
      · gg_dynamics_function. R - function needed to generate plots of temporal series (right panel of Figures 2 and 3) 
      
      
- 7.sensitivity_analysis_Z_biomass.R (7.R) - Code to perform sensitivity analysis of Z coefficient (2/3) used in allometric equation of 
  biomass and its interaction with imputation levels in 1.1.R and 1.2.R. Here we generate Supplementary Figures 8 and 9. 
  
  
  
- 8.effect_size_metafor_checking.R (8.R) - Code where we use package metafor to check our own LRR analysis. 



- palettes_labels.R - Code where objects for visualization colors and labels are created to ensure visual consistency.















