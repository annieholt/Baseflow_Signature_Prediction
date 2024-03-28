# script for plotting random forest model outputs

library(tidyverse)

rf_performance = read.csv('E:/SDSU_GEOG/Thesis/Data/RandomForest/rf_performance_output.csv')

rf_performance_test = rf_performance %>% 
  filter(cv_r2 > 0.4) %>% 
  filter(IncMSE > 0)

rf_performance_2 = rf_performance %>% 
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
ggplot(rf_performance_test, aes(x = sig, y = variable, fill = IncMSE, size = IncMSE)) +
  geom_point(shape = 21, color = 'black') +
  # scale_color_viridis_c()+
  scale_fill_gradient2(low = "grey",high = "darkgreen", name = "IncMSE",
                        midpoint = median(rf_performance$IncMSE)) +
  scale_size_continuous(range = c(4, 12)) +
  labs(
    y = NULL,  # No y-axis label
    x = NULL,
  ) +
  ggtitle("Variable Importance") +  # Add plot title
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
