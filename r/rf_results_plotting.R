# script for plotting random forest model outputs

library(tidyverse)

rf_performance = read.csv('E:/SDSU_GEOG/Thesis/Data/RandomForest/outputs/rf_performance_caravan_output.csv')
rf_performance_plus = read.csv('E:/SDSU_GEOG/Thesis/Data/RandomForest/outputs/rf_performance_caravan_plus_output.csv')
rf_performance_camels = read.csv('E:/SDSU_GEOG/Thesis/Data/RandomForest/outputs/rf_performance_camels_plus_output.csv')

rf_performance_wet = read.csv('E:/SDSU_GEOG/Thesis/Data/RandomForest/outputs/rf_performance_caravan_plus_wet_output.csv')
rf_performance_eco = read.csv('E:/SDSU_GEOG/Thesis/Data/RandomForest/outputs/rf_performance_caravan_plus_eco_eastern_forests_output.csv')

#### Performance Plotting ####

rf_performance_r2 = rf_performance %>% 
  group_by(sig) %>% 
  summarize(cv_r2 = median(cv_r2)) %>% 
  rename(cv_r2_caravan = cv_r2)

# sorting  high to low predictability
rf_performance_r2$sig <- factor(rf_performance_r2$sig, levels = rf_performance_r2$sig[order(rf_performance_r2$cv_r2_caravan, decreasing = TRUE)])


rf_performance_plus_r2 = rf_performance_plus %>% 
  group_by(sig) %>% 
  summarize(cv_r2 = median(cv_r2)) %>% 
  rename(cv_r2_plus = cv_r2)


rf_performance_camels_r2 = rf_formance_camels %>% 
  group_by(sig) %>% 
  summarize(cv_r2 = median(cv_r2)) %>% 
  rename(cv_r2_camels = cv_r2)

rf_performance_final =  rf_performance_r2 %>% 
  left_join(rf_performance_plus_r2, by = "sig") %>%
  # left_join(rf_performance_camels, by = "sig") %>% 
  pivot_longer(-sig, names_to = "r2_variable", values_to = "value") %>% 
  arrange(r2_variable) %>% 
  filter(sig != "Spearmans_rho") %>% 
  filter(sig != "MRC_num_segments") %>% 
  filter(sig != "RecessionParameters_a")
  
  

# eventually, want a bar plot with two different R2 for each signature
# R2 when inlcude or remove the new catchment attributes
# to do so, have the data long, with another column labeling R2with and R2without
# then, geom_col(data = df, aes(x = sigs, y = value, fill = r2_variable)

# boxplot, R2 on y axis, signature on x axis
ggplot(rf_performance_final, aes(x = sig, y = value, fill = r2_variable)) +
  geom_col(position = "dodge", width = 0.6) +
  labs(
    y = "Mean R2",
    x = NULL,
    title = NULL
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 24),
    text = element_text(size = 14),   # Increase text (axis labels, title) size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis labels diagonally
    axis.text.y = element_text(size = 14),  # Adjust y-axis text size
    legend.position = "left",  # Move legend to the left side
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),  # Adjust legend text size
    legend.key.size = unit(2, "lines")  # Adjust the size of the legend color key
  )+
  scale_fill_manual(name = "RF Model Performance",  # Changing legend name
                      labels = c("Caravan attributes", "Plus new attributes"),
                      values = c("blue", "green"))


# ggsave("E:/SDSU_GEOG/Thesis/Data/RandomForest/figures_final/rf_performance_r2.png", width = 12, height = 6, dpi = 300,bg = "white")


# trying slightly different aesthetics
ggplot(rf_performance_r2, aes(x = sig, y = cv_r2)) +
  geom_col(width = 0.6) +  # Adjust the width of the bars
  labs(
    y = "R2",
    x = NULL,
    title = "Model Performance, Cross Validation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),  # Adjust axis title size
    legend.position = "none"  # Remove legend
  )



#### VARIABLE IMPORTANCE PLOTTING (incMSE) ####
rf_performance_test = rf_performance %>% 
  filter(cv_r2 > 0.4) %>% 
  filter(IncMSE > 0)

rf_performance_2 = rf_performance_plus %>% 
  mutate(IncMSE_Category = case_when(
    IncMSE <= 5 ~ "<5",
    IncMSE > 5 & IncMSE <= 15 ~ "5-15",
    IncMSE > 15 & IncMSE <= 25 ~ "15-25",
    IncMSE > 25 & IncMSE <= 35 ~ "25-35",
    IncMSE > 35 ~ "35+",
  ))

# Define the color and size scales
color_scale <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", 'pink')
size_scale <- c(3, 6, 9, 12, 15)

# Define the categories and their corresponding colors and sizes
categories <- c("<5", "5-15", "15-25", "25-35","35+" )

# Create the heatmap plot with categorized colors and sizes
ggplot(rf_performance_2, aes(x = sig, y = variable, color = IncMSE_Category, size = IncMSE_Category)) +
  geom_point() +
  scale_color_manual(values = color_scale, name = "IncMSE") +
  scale_size_manual(values = size_scale) +
  labs(
    y = NULL,  # No y-axis label
    x = NULL,
    title = "Variable Importance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24),
    text = element_text(size = 20),   # Increase text (axis labels, title) size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels diagonally
    axis.text.y = element_text(size = 14),  # Adjust y-axis text size
    legend.position = "left",  # Move legend to the left side
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.key.size = unit(2, "lines")  # Adjust the size of the legend color key
  ) +
  guides(size = FALSE)

# trying the continuous scale plotting
# I still don't love the color scales since a lot of the data is middle range, but ok for now


rf_variables = rf_performance_eco %>% 
  filter(sig != "Spearmans_rho") %>% 
  filter(sig != "MRC_num_segments") %>% 
  filter(sig != "RecessionParameters_a") %>% 
  mutate(IncMSE = ifelse(IncMSE <0, 0, IncMSE))

rf_variables_2 = rf_performance_plus %>% 
  filter(sig != "Spearmans_rho") %>% 
  filter(sig != "MRC_num_segments") %>% 
  filter(sig != "RecessionParameters_a") %>% 
  mutate(IncMSE = ifelse(IncMSE <0, 0, IncMSE))

desired_order <- rev(c("p_mean", "pet_mean", "aridity", "frac_snow","moisture_index", "seasonality", "high_prec_freq", "high_prec_dur",
                   "low_prec_freq", "low_prec_dur",
                   "ele_mt_smn", "slp_dg_sav", 'area',
                   'for_pc_sse',
                   "swc_pc_syr", "cly_pc_sav", "slt_pc_sav", "snd_pc_sav", "soc_th_sav", "kar_pc_sse",
                   'fresh_no_giw', 'geol_major_age_ma' ))

rf_variables$variable <- factor(rf_variables$variable, levels = desired_order)

ggplot(rf_variables, aes(x = sig, y = variable, fill = IncMSE, size = IncMSE)) +
  geom_point(shape = 21, color = 'black') +
  # scale_color_viridis_c()+
  scale_fill_gradient2(low = "white",high = "black", name = "IncMSE",
                        midpoint = median(rf_variables_2$IncMSE),  limits = c(0, 35)) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    y = NULL,  # No y-axis label
    x = NULL,
  ) + # Add plot title
  theme_minimal() +                    # Use a minimal theme
  theme(
    plot.title = element_text(size=24),
    text = element_text(size = 20),   # Increase text (axis labels, title) size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels diagonally
    axis.text.y = element_text(size = 14),  # Adjust y-axis text size
    legend.position = "left",  # Move legend to the left side
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),  # Adjust legend text size
    legend.key.size = unit(2, "lines")  # Adjust the size of the legend color key
  ) +
  guides(size = FALSE)

ggsave("E:/SDSU_GEOG/Thesis/Data/RandomForest/figures_final/rf_plus_eco_eastern_incmse.png", width = 10, height = 8, dpi = 300,bg = "white")
