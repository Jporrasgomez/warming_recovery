





rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, tidyverse, DT, viridis, ggrepel, codyn, vegan, eulerr, ggplot2, ggthemes, ggpubr, ggforce )#


flora_abrich <- read.csv("data/processed_data/flora_abrich.csv")



source("code/palettes_labels.R")
palette <- palette_CB
labels <- labels1



####### DIFFERENCES AT SAMPLING x TREATMENT LEVEL #############

{
species_ab_sampling <- flora_abrich |> 
  group_by(code, date, sampling, treatment) |> 
  summarise(abundance = mean(abundance_s, na.rm = T),
            .groups = "drop_last") |> 
  ungroup() |> 
  mutate(id = paste0(as.character(treatment), "/" , as.character(sampling))) |> 
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))



treats <- c("c", "w", "p", "wp")
palette_treats <- c("#12D08C", "#E4345A", "#1CC4FD", "#842F8A")
labels_treats <- c("Control", "Warmed-only", "Perturbed-only", "Combined")
  
list1 <- list()
gglist1 <- list()
gglist2 <- list()
count = 0

for(i in 1:length(treats)){
  
  count = count + 1
  
  list1[[count]] <- subset(species_ab_sampling, treatment == treats[i])
  
  sp_wide <- list1[[count]] |>
    pivot_wider(id_cols = sampling,
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) |> 
    arrange(as.numeric(sampling)) |> 
    column_to_rownames(var = "sampling")
  
  
  #Relative abundance
  
  # Perform NMDS using Bray-Curtis distance
  nmds_res <- metaMDS(sp_wide, distance = "bray", k = 2, trymax = 250, maxit = 999) 
  
  # Extract NMDS sample scores
  nmds_samples <- as.data.frame(scores(nmds_res, display = "sites"))
  
  # Extract NMDS species scores (optional)
  nmds_species <- as.data.frame(scores(nmds_res, display = "species"))
  
  gglist1[[count]] <- ggplot() +
    geom_text_repel(data = nmds_species |> 
                      rownames_to_column(var = "sp"),
                    aes(x = NMDS1, y = NMDS2, label = sp),
                    color = "grey",
                    max.overlaps = 30) +
    geom_point(data = nmds_samples |> 
                 rownames_to_column(var = "sampling"),
               aes(x = NMDS1, y = NMDS2),
               size = 3,
               color = palette_treats[i]) +
    geom_text_repel(data = nmds_samples |> 
                      rownames_to_column(var = "sampling"),
                    aes(x = NMDS1, y = NMDS2, label = sampling),
                    max.overlaps = 9,
                    color = palette_treats[i]) +
    geom_path(data = nmds_samples |> 
                rownames_to_column(var = "sampling"),
              aes(x = NMDS1, y = NMDS2),
              color = palette_treats[i]) +
    geom_hline(aes(yintercept = 0), color = "gray52", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), color = "gray52", linetype = "dashed") +
    labs(title = paste(labels_treats[i]),
         subtitle = paste0("Stress = ", round(nmds_res$stress, 3)),
         x = "NMDS1",
         y = "NMDS2") + 
    theme1
  
  
}


gg_sp_comp_dynamics <- 
ggarrange(
  gglist1[[1]],
  gglist1[[2]],
  gglist1[[3]],
  gglist1[[4]],
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2)


#Stress-level of NMDS:
#  
#  Higher than 0.2 is poor (risks for false interpretation).
#0.1 - 0.2 is fair (some distances can be misleading for interpretation).
#0.05 - 0.1 is good (can be confident in inferences from plot).
#Less than 0.05 is excellent (this can be rare).





##############################################################################
#### 1. SPECIES COMPOSITION DIFFERENCES AT TREATMENT LEVEL #######################
##############################################################################

# Abundance matrix at treatment level (4 dots : 4 replicates per treatment)

species_ab_treatment <- flora_abrich |> 
  filter(sampling != "0") |>  # Post-disturbance samplings
  group_by(code, plot, treatment) |> 
  summarise(abundance = mean(abundance_s, na.rm = T),
            .groups = "drop_last") |>
  ungroup() |> 
  mutate(id_treat_plot = paste0(as.character(treatment), "/" , as.character(plot)))


sp_wide_treatment <- species_ab_treatment |>
  pivot_wider(id_cols = c(treatment, id_treat_plot),
              names_from = code,
              values_from = abundance,
              values_fill = list(abundance = 0)) |> 
  ungroup() |> 
  select(-`NA`) |> 
  as_tibble()

abudance_data_treatment <- sp_wide_treatment |> 
  select(-treatment, -id_treat_plot)

abundance_data_treatment_hellinger <- vegan::decostand(abudance_data_treatment, method = "hellinger")


abundance_treatment_list <- list(abudance_data_treatment, abundance_data_treatment_hellinger)
# Compute Bray-Curtis distance matrix: i = 1
# Compute with  Hellinger transformation: i = 2
 i = 1
 
 
distance_matrix_treatment_bc <- vegan::vegdist(as.data.frame(abundance_treatment_list[[i]]), method = "bray")

# Run NMDS (2 dimensions, 100 tries)
nmds_bc_treatment <- metaMDS(distance_matrix_treatment_bc, k = 2, trymax = 250, maxit = 999)
nmds_df_treatment <- data.frame(
  NMDS1 = nmds_bc_treatment$points[, 1],
  NMDS2 = nmds_bc_treatment$points[, 2],
  treatment = sp_wide_treatment$treatment
)

gg_sp_comp_treatment <- 
  
  ggplot(nmds_df_treatment, aes(x = NMDS1, y = NMDS2)) +  # <- IMPORTANTE: Quita el color general de aquí
  
  # Stat ellipse
  stat_ellipse(geom = "polygon", aes(fill = treatment, color = treatment),
               alpha = 0.20, show.legend = FALSE, level = 0.95) + 
  

  geom_point(aes(fill = treatment, shape = treatment, color = treatment),
             size = 8, 
             color = "black",
             stroke = 1, 
             show.legend = TRUE) +
  
  # Axis lines
  geom_hline(yintercept = 0, color = "gray52", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray52", linetype = "dashed") +
  
  
  # Color and shape codes
  
  scale_fill_manual(values = palette_CB, 
                    labels = c("c" = "Control", "w" = "Warmed-only", 
                               "p" = "Perturbed-only", "wp" = "Combined")) +
  
  scale_color_manual(values = palette_CB,
                     labels = c("c" = "Control", "w" = "Warmed-only", 
                                "p" = "Perturbed-only", "wp" = "Combined")) +
  
  scale_shape_manual(values = c("c"  = 21,  # Circle
                                "w"  = 22,  # Square
                                "p"  = 24,  # Triangle
                                "wp" = 23), # Rhombus
                     labels = c("c" = "Control", "w" = "Warmed-only", 
                                "p" = "Perturbed-only", "wp" = "Combined")) +
  
  labs(
    subtitle = paste0("Stress = ", round(nmds_bc_treatment$stress, 3)),
    x = "NMDS1", y = "NMDS2", 
    fill = NULL, color = NULL, shape = NULL) +
  theme1


}


print(list_sorensen)
#View(sorensen_df)
gg_sp_comp_dynamics
#ggsave("results/NMDS_sp_composition_dyn.png", plot = gg_sp_comp_dynamics, dpi = 600)

gg_sp_comp_treatment
ggsave("results/NMDS_sp_composition_agg.png", plot = gg_sp_comp_treatment, dpi = 600)













########### 1.2. STATISTICAL ANALYSIS: PERMANOVA ##########################

adonis_treatment <- adonis2(
  distance_matrix_treatment_bc ~ treatment,  
  data = sp_wide_treatment,                 
  permutations = 999,                      
  method = "bray"
)

print(adonis_treatment)
adonis_treatment$`Pr(>F)`[1]
adonis_treatment$R2[1] # R2 reflects how much variance the treatment explain

# Dispersion test
bd <- betadisper(distance_matrix_treatment_bc, sp_wide_treatment$treatment)

test_dispersion <- permutest(bd) 
test_dispersion
test_dispersion[[1]]$`Pr(>F)`[1] #p-value is bigger than 0.05. 
plot(bd)
boxplot(bd)

permutest(bd, pairwise = TRUE)
TukeyHSD(bd)  # Differences of dispersion across treatment comparisons 

library(pairwiseAdonis)

# Pairwise comparisons
pw_adonis <- pairwise.adonis(
  x           = distance_matrix_treatment_bc,           
  factors     = sp_wide_treatment$treatment,  
  perm        = 999,                         
  p.adjust.m  = "BH")                        

print(pw_adonis)

pw_adonis |>  write.csv("results/adonis_pairwise.csv")


# Data for colonizers

species_ab_plot <- flora_abrich |> 
  select(date, code, sampling, plot, treatment, family, genus_level, species_level, abundance_s)

species_ab_plot <- species_ab_plot |> 
  mutate(id = paste0(treatment, "/", sampling, "/", plot))

species_ab_plot <- species_ab_plot |> 
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))

sp_wide_plot <- species_ab_plot |>
  pivot_wider(id_cols = c(sampling, date, treatment, plot),
              names_from = code,
              values_from = abundance_s,
              values_fill = list(abundance = 0)) |>
  mutate(across(everything(), ~ replace_na(., 0)))

sp_wide_plot |> write.csv("data/processed_data/relative_abudance_species_time.csv")





############### SORENSEN  ################################

list_sorensen <- list()
list_sorensen_df <- list()
samps <- sort(unique(flora_abrich$sampling))

for(i in 1:length(samps)){
  
  sp_wide <- species_ab_sampling |> 
    filter(sampling == samps[i]) |> 
    pivot_wider(id_cols = c(sampling, date, treatment, id),
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) |> 
    ungroup() |> 
    as_tibble()
  
  
  abundance_data <- sp_wide |>
    select(-treatment, -sampling, -date, -id)
  
  sorensen <- vegdist(abundance_data, method = "bray", binary = TRUE)
  sorensen <- as.matrix(sorensen)  
  rownames(sorensen) <- sp_wide$id
  colnames(sorensen) <- sp_wide$id
  sorensen[upper.tri(sorensen)] <- NA
  
  sorensen_df <- sorensen |> 
    as.data.frame() |>
    rownames_to_column(var = "row_name") |>
    pivot_longer(-row_name, names_to = "col_name", values_to = "value") |> 
    filter(!is.na(value))
  
  list_sorensen[[i]] <- sorensen
  list_sorensen_df[[i]] <- sorensen_df
  
  
}



sorensen_df <- bind_rows(list_sorensen_df) |> 
  filter(!value == "0") |>
  separate(col_name, into = c("treatment_x", "sampling_x"), sep = "/") |>
  separate(row_name, into = c("treatment_y", "sampling_y"), sep = "/") |> 
  select(-sampling_x) |> 
  rename(sampling = sampling_y) |> 
  mutate(sampling = factor(as.numeric(sampling), levels = sort(unique(as.numeric(sampling))))) |> 
  arrange(sampling) |> 
  mutate(comparison = paste0(treatment_x, "-", treatment_y)) |> 
  select(-treatment_x, -treatment_y) |> 
  mutate(comparison = ifelse(comparison == "p-c", "c-p", comparison),
         comparison = ifelse(comparison == "w-c", "c-w", comparison),
         comparison = ifelse(comparison == "wp-c", "c-wp", comparison))





##############################################################################
#### SPECIES COMPOSITION DIFFERENCES AT PLOT X SAMPLING LEVEL ###########################
##############################################################################

species_ab_plot <- flora_abrich |> 
  select(date, code, sampling, plot, treatment, family, genus_level, species_level, abundance_s) |> 
  mutate(id = paste0(treatment, "/", sampling, "/", plot)) |> 
  filter(!(sampling == "1" & treatment %in% c("p", "wp"))) # Removing wp and p from sampling 1 because there is no information


# Initialize a list to store NMDS results
samps <- sort(unique(species_ab_plot$sampling))
nmds_list <- list()
stress_values <- data.frame()

for (i in seq_along(samps)) {
  sp_wide_subset <- species_ab_plot |>
    filter(sampling == samps[i]) |>    # Filter for each sampling event
    pivot_wider(id_cols = c(sampling, date, treatment, plot),
                names_from = code,
                values_from = abundance_s,
                values_fill = list(abundance = 0))|>
    mutate(across(everything(), ~ replace_na(., 0)))
  
  # Prepare abundance data
  abundance_data_subset <- sp_wide_subset |> select(-treatment, -sampling, -date, -plot)
  
  # Compute Bray-Curtis distance
  if (nrow(abundance_data_subset) > 1) {  # Avoid NMDS failure for single-row cases
    distance_matrix_bc <- vegan::vegdist(abundance_data_subset, method = "bray")
    
    # Run NMDS
    nmds_result <- metaMDS(distance_matrix_bc, k = 2, trymax = 250, maxit = 999)
    
    # Store NMDS scores and stress values
    nmds_df_temp <- data.frame(
      NMDS1 = nmds_result$points[, 1],
      NMDS2 = nmds_result$points[, 2],
      treatment = sp_wide_subset$treatment,
      sampling = sp_wide_subset$sampling,
      plot = sp_wide_subset$plot,
      date = sp_wide_subset$date
    )
    
    # Store stress values
    stress_values <- rbind(stress_values, data.frame(sampling = samps[i], stress = round(nmds_result$stress, 3)))
    
    # Store NMDS results
    nmds_list[[i]] <- nmds_df_temp
  }
}

# Combine NMDS results into one dataframe
nmds_df_1x1sampling <- bind_rows(nmds_list)

# Merge stress values into NMDS dataframe
nmds_df_1x1sampling <- left_join(nmds_df_1x1sampling, stress_values, by = "sampling")

nmds_df_1x1sampling <- nmds_df_1x1sampling |> 
  mutate(date_sampling = as.factor(date))

# Create NMDS plot with stress values in facet labels


gg_samplings <- 
    ggplot(nmds_df_1x1sampling, aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment, shape = treatment)) +
    #facet_wrap(~ paste0(sampling, " (Stress = ", stress, ")"), ncol = 5, nrow = 5) +
    facet_wrap(~ date_sampling, ncol = 5, nrow = 5, scales = "free") + 
    geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +

  stat_ellipse(geom = "polygon", aes(fill = treatment),
               alpha = 0.12, show.legend = FALSE, level = 0.95) +
  geom_point(size = 2, alpha = 0.8) +
    labs(
      x = "NMDS1",
      y = "NMDS2"
    ) +
    
  scale_fill_manual(values = palette_CB, 
                    labels = c("c" = "Control", "w" = "Warmed-only", 
                               "p" = "Perturbed-only", "wp" = "Combined")) +
  
  scale_color_manual(values = palette_CB,
                     labels = c("c" = "Control", "w" = "Warmed-only", 
                                "p" = "Perturbed-only", "wp" = "Combined")) +
  
  scale_shape_manual(values = c("c"  = 21,  # Circle
                                "w"  = 22,  # Square
                                "p"  = 24,  # Triangle
                                "wp" = 23), # Rhombus
                     labels = c("c" = "Control", "w" = "Warmed-only", 
                                "p" = "Perturbed-only", "wp" = "Combined")) +
  
  labs(
    x = "NMDS1", y = "NMDS2", 
    fill = NULL, color = NULL, shape = NULL) +
  theme1 + 
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"), # Adjust facet label size
      axis.text.y = element_text(size = 10), 
      axis.text.x = element_text(size = 10), 
    
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      #legend.axis.line = element_text(size = 12)
    ) 
    
   
 


 
  print(gg_samplings)
  print(stress_values)
  mean(stress_values$stress)
  sd(stress_values$stress)

ggsave("results/species_composition_allsamplings.png", plot = gg_samplings, dpi = 600)



