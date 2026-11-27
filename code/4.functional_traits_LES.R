



rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #se quitan todos los paquetes (limpiamos R)
pacman::p_load(dplyr, tidyr, tidyverse, ggplot2, BIEN, ape, maps, sf, rtry, ggrepel)




flora_abrich <- read.csv("data/processed_data/flora_abrich.csv")
traits <- read.csv('data/traits_data.csv')
source("code/palettes_labels.R")


flora_abrich <- flora_abrich %>% 
  filter(!is.na(code))

traits <- traits %>%
  mutate(across(where(is.character), as.factor)) %>% 
  mutate(species = recode(species, "Anagallis arvensis" = "Lysimachia arvensis"))%>%
  filter(species != "Medicago polymorpha") %>%
  mutate(species = droplevels(species))

species_code <- read.csv("data/species_code.csv") %>% 
  
  mutate(species = recode(species, "CAPSELLA BURSA-PASTORIS" = "Capsella bursa-pastoris")) %>% 
  mutate(across(where(is.character), as.factor))


## Checking which species are absent
checking <- anti_join(traits, species_code, by = "species") %>% 
  distinct(species, .keep_all = T) %>% 
  print()#Anagallis arvensis = Lysimachia arvensis

traits <- traits  %>% 
  left_join(species_code, by = "species") %>% 
  select(c("code", "species", "trait_name", "trait_ID", "database", "trait_value")) %>% 
  group_by(species, trait_name) %>% 
  mutate(n_observations = n()) %>% 
  ungroup() %>% 
  mutate(
    code_n_obs = paste(code, n_observations, sep = ", ")
  )




# Remove species that are "sp". This is: working just with identified species up to species level. 
flora_abrich <- flora_abrich %>% 
  filter(!code %in% c("chsp", "amsp", "cisp", "casp"))

traits <- traits %>% 
  filter(!code %in% c("chsp", "amsp")) %>% 
  filter(!species %in% c("Erophila verna", "Stellaria media",
                      "Lotus corniculatus")) # I remove also Erophila, Stellaria and Lotus because we did not measured them in the field


######## WE WILL FOCUS ON LEAF ECONOMIC SPECTRUM #########

traits <- traits %>% 
  filter(
    trait_name %in% c("LDMC", "leafN", "SLA.ex", "SLA.inc")
  )


# Remove outliers based on Z-scores
traits_cleaned <- traits %>%
  group_by(trait_name, species) %>% 
  mutate(
    mean_trait = mean(trait_value, na.rm = TRUE),
    sd_trait = sd(trait_value, na.rm = TRUE)
  ) %>%
  mutate(
    z_score = (trait_value - mean_trait) / sd_trait
  ) %>% 
  filter(
    n_observations <= 5 | abs(z_score) <= 4  # Do not filter outliers for n_observations <= 5
  ) %>%
  ungroup() %>%
  select(-mean_trait, -sd_trait, -z_score)  # Remove temporary columns if not needed

# View the cleaned data
outliers <- anti_join(traits, traits_cleaned)

########################################
## CALCULATING MEAN VALUES FOR TRAITS
#######################################


#Calculating the average traits per taxonomic groups: torilis sp, poaceae, asteraceae and orchidaceae
sp_poaceae <- c("Avena sterilis", "Bromus hordeaceus", "Bromus sterilis", "Cynosurus echinatus",
                "Elymus repens", "Hordeum murinum", "Poa annua", "Poa bulbosa",
                "Lolium perenne", "Gaudinia fragilis")

sp_asteraceae <- c("Crepis capillaris", "Hypochaeris radicata", "Leontodon hispidus")

sp_orchidaceae <- c("Anacamptis pyramidalis", "Ophrys apifera")

sp_torilis <- c("Torilis nodosa", "Torilis arvensis")


labels_sp <- c("poaceae", "asteraceae", "orchidaceae", "tosp")

list_sp <- list(sp_poaceae, sp_asteraceae, sp_orchidaceae, sp_torilis)
# 1          # 2            # 3             # 4
trait_values_species <- list()
trait_values_aggtax <- list()


# In the case of species that we have aggregated at certain taxon levels, we have 
# to avoid bias based on the number of trait datapoints available per each species in the online repositories (TRY and others).
# Thus, we will FIRST calculate an average value of each trait per species and then calculate
# the average value per trait and taxonomical aggregation. 

for (i in c(1:4)){
  
  # First: calculating average value of trait per species
  data1 <- traits_cleaned %>% 
    filter(species %in% list_sp[[i]]) %>% 
    group_by(trait_name, species, code) %>% 
    summarize(species_mean_trait = mean(trait_value, na.rm = T),
              species_sd_trait   = sd(trait_value, na.rm = T), 
              n_obs_sp           = n())
  
  
  trait_values_species[[i]] <- data1 
  
  # Second: calculating average value of trait at aggregated level using previous calculation (means of means)
  data2 <- data1 %>% 
    group_by(trait_name) %>% 
    summarize(trait_mean       = mean(species_mean_trait, na.rm = T),
              trait_sd = sd(species_mean_trait, na.rm = T),
    ) %>% 
    mutate(code = labels_sp[i])
  
  trait_values_aggtax[[i]] <- data2
  
  
}


trait_means_poaceae_splevel <- trait_values_species[[1]] %>%  print()
trait_means_asteraceae_splevel <- trait_values_species[[2]] %>%  print()
trait_means_orchidaceae_splevel <- trait_values_species[[3]] %>%  print()
trait_means_tosp_splevel <- trait_values_species[[4]] %>%  print()


trait_means_poaceae <- trait_values_aggtax[[1]] %>%  print()
trait_means_asteraceae <- trait_values_aggtax[[2]] %>%  print()
trait_means_orchidaceae <- trait_values_aggtax[[3]] %>%  print()
trait_means_tosp <- trait_values_aggtax[[4]] %>%  print()




sp_all <- c(sp_poaceae, sp_asteraceae, sp_orchidaceae, sp_torilis)

sp_others <- setdiff(unique(traits_cleaned$species), sp_all)

# Here we do not need to perform 2 steps, since we directly calculate one trait value per species using
# all available data in the repositories

trait_means_others0 <- traits_cleaned %>% 
  filter(species %in% sp_others) %>% 
  group_by(trait_name, species, code) %>% 
  summarize(trait_mean = mean(trait_value, na.rm = T),
            trait_sd   = sd(trait_value, na.rm = T), 
            n_obs_sp           = n(), 
  )


trait_means_others <- trait_means_others0 %>% 
  ungroup() %>%
  select(trait_name, trait_mean, trait_sd, code)




trait_means <- bind_rows(trait_means_others,
                         trait_means_poaceae,
                         trait_means_asteraceae,
                         trait_means_orchidaceae,
                         trait_means_tosp)


#The function checks which species are present in traits_mean$species but not present in flora$species.
setdiff(unique(trait_means$code), unique(flora_abrich$code))
setdiff(unique(flora_abrich$code), unique(trait_means$code))



setdiff(unique(species_code$code), unique(trait_means$code))
setdiff(unique(species_code$code), unique(traits$code))

# I have to check on Myosotis, Sonchus and Cirsium. 

#Changes that have been made: 
# - Sonchus sp has been substituted by "asteraceae" in the flora_abrich database
# - For the species Cirsium sp, Myosotis discolor and Cardamine sp we do not have functional traits. 
# What can I do about Cirsium? it is one of the most abudance species at the end of the samplings. 



#######################################
########## CWM analysis RAW ###########
#######################################

library(FD)
library(psych)


####### TRAIT MATRIX: 

traits_mean_wide <- trait_means %>%
  ungroup() %>%                           # Remove grouping structure
  select(-trait_sd) %>%  # Remove unnecessary columns
  pivot_wider(
    names_from = trait_name,              # Columns will be based on trait_name levels
    values_from = trait_mean              # Values will come from trait_mean
  )


##|  Treating SLA and LA traits. 
##| We have 2 traits for SLA: SLA.ex and SLA.inc
##| We have 3 traits for LA: LA.ex, LA.un and LA.inc. 
##| They all differ if they include or no the pedunculum of the leaf. 
##| When there is one trait available for a species, the others are not. So we will keep just one value per SLA and LA and species
##| 

traits_mean_wide <- traits_mean_wide %>% 
  mutate(SLA.inc = ifelse(is.na(SLA.inc), SLA.ex, SLA.inc)) %>% 
  select(!"SLA.ex") %>% 
  rename(SLA = SLA.inc ) 


traits_data_final <- traits_mean_wide %>% 
  left_join(species_code) %>% 
  select(species, family, sampling_level, code, SLA, LDMC, leafN)

#traits_data_final %>%  write.csv("results/traits_data_final.csv")

##  ---REMOVING SPECIES WITH LOW NUMBER OF TRAITS AVAILABLE---

traits_mean_wide %>%
  mutate(
    missing_data = rowSums(is.na(.))
  ) %>%
  select(code, missing_data) %>%
  mutate(
    code = factor(code, levels = code[order(missing_data, decreasing = TRUE)])
  ) %>%
  ggplot(aes(x = code, y = missing_data)) +
  geom_col(fill = "steelblue") +
  labs(x = "Species code", y = "% Missing traits data", 
       title = "Number of traits without data per species (n traits = 10)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

## According to Carmona et al. 2021, we have to remove the species with less than the 50% of the traits. 
# These are: libi and rapa. 
# Before deleting them, we can check the relevance of these species (in checking.results we can do it)
# rapa is a species with a lot of abundance and n_observations, but not one of the highest

traits_mean_wide <- traits_mean_wide %>% 
  filter(code != "libi") %>% 
  droplevels()

## REMOVING TRAITS WITH HIGH NUMBER OF NA

missing_data <- colSums(is.na(traits_mean_wide)) / nrow(traits_mean_wide) * 100
print(missing_data)

data.frame(
  variable = names(missing_data),
  missing_perc = as.numeric(missing_data)
) %>%
  mutate(
    variable = factor(variable, levels = variable[order(missing_perc, decreasing = TRUE)])
  ) %>% 
  ggplot(aes(x = variable, y = missing_perc)) +
  geom_col(fill = "seagreen") +
  geom_hline(yintercept = 50, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Trait", y = "% Missing data", 
       title = "Percentage of Missing Trait Data per Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


################ CHOOSING TRAITS ##############################################################################

##
traits_final <- traits_mean_wide %>% 
  select(c("code",
           "SLA", 
           "LDMC",
           "leafN",
           #"LA", 
           #"seed.mass", 
           #"vegetation.height"
  ))  ### Choose here the FUNCTIONAL TRAITS WE WANT

#################################################################################################################


#Now we have 2 traits above 50 % NA, 1 = 50% and 1 around 45%. I would delete them all
# given that the following one is araound 15%

# We now have only 5 NA in trait leafN 

# Prepare the database for the FD function: 

traits_matrix <- as.data.frame(traits_final)

rownames(traits_matrix) <- traits_matrix$code

traits_matrix <- traits_matrix %>% 
  arrange(code)

traits_matrix <- traits_matrix %>% 
  select(!code)






### ABUNDANCE MATRIX #####

flora_abrich <- flora_abrich %>% 
  filter(!code %in% c("cisp", "casp", "mydi",  #Deleting the species for which we do not have traits
                      "libi", "amsp", "rapa"))  #Deleting the species for which we had more than 50% of NAs




########### COMMUNITY WEIGHTED MEAN (CWM) ############
# Abundance matrix 2: plot level (dynamics)

abundance_matrix_dynamics <- flora_abrich %>% 
  ungroup() %>% 
  select(sampling, date, treatment, plot, abundance_s, code) %>% 
  mutate(
    com = paste(sampling, treatment, plot, sep = "/")
  ) %>% 
  select(com, abundance_s, code) %>% 
  mutate(abundance_s = abundance_s/100) %>% # abundance = (0,1)
  pivot_wider(
    names_from = code, 
    values_from = abundance_s) %>%
  select(com, sort(setdiff(names(.), "com"))) %>%
  column_to_rownames("com") %>% 
  mutate(across(everything(), ~ replace_na(., 0))) %>%  # NA = 0 
  as.data.frame() 


# Calculation of cwm with function functcomp

#cwm at plot level (dynamics)

cwm_plot <- FD::functcomp(as.matrix(traits_matrix), as.matrix(abundance_matrix_dynamics))
cwm_plot$com <- rownames(abundance_matrix_dynamics); rownames(cwm_plot) <- NULL


cwm_plot <- cwm_plot %>% 
  mutate(
    sampling = sapply(strsplit(com, "/"), function(x) x[1]),
    treatment = sapply(strsplit(com, "/"), function(x) x[2]),
    plot = sapply(strsplit(com, "/"), function(x) x[3])
  ) %>% 
  select(-com)
  
  

sampling_dates <- read.csv("data/sampling_dates.csv") %>% 
  mutate(sampling = as.factor(sampling),
         date = ymd(date), 
         month = month(date, label = TRUE),
         day = day(date), 
         year = year(date)) %>% 
  select(sampling, date, year 
         #, day, month, one_month_window, omw_date
  ) %>% 
  mutate(across(where(is.character), as.factor))




cwm_plot_db <- cwm_plot %>% 
  merge(sampling_dates)

# Storing results for GLMM and LRR results

cwm_plot_db %>%  write.csv("data/processed_data/cwm_plot_db.csv")





