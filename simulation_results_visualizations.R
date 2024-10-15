library(ggplot2)

load("results/Simulation_Design_1_results_20_reps.rda")


list_to_analze <- c("value_ml_fr_indv_set", "RU_ml_fr_indv_set", "AU_ml_fr_indv_set",
                    "value_ml_fr_indv_cluster_set", "RU_ml_fr_indv_cluster_set", "AU_ml_fr_indv_cluster_set",
                    "FURG_indv_set", "FUTR_indv_set", "FURG_indv_cluster_set", "FUTR_indv_cluster_set",
                    "RU_ml_fr_indv_CATE_set", "AU_ml_fr_indv_CATE_set",
                    "RU_ml_fr_indv_cluster_CATE_set", "AU_ml_fr_indv_cluster_CATE_set")

average_values_delta0 <- list()
average_values_delta20 <- list()

for (element in list_to_analze) {
  average_values_delta0[[element]] <- mean(sapply(results, function(x) x[[element]][1]))
  average_values_delta20[[element]] <- mean(sapply(results, function(x) x[[element]][27]))
}

# Convert the list to a named vector for better readability
average_values_delta0 <- unlist(average_values_delta0)
average_values_delta20 <- unlist(average_values_delta20)

# Print the average values
print(round(average_values_delta0,3))
print(average_values_delta20)


# --------------------------------------------
# Convert the list-elements into one dataframe
# --------------------------------------------

# Define the names of the elements you want to convert into a dataframe
elements_to_convert <- list_to_analze

# Initialize a list to store the extracted elements
extracted_elements <- list()

# Extract each element from all list elements
for (element in elements_to_convert) {
  # Use sapply to extract the element from each list and store as a column in a dataframe
  extracted_elements[[element]] <- t(sapply(results, function(x) x[[element]]))
}


# extract each element from the list and convert it to a dataframe, and calculate the column mean
value_ml_fr_indv_set_mean <- apply(t(sapply(results, function(x) x$value_ml_fr_indv_set)),2, mean)
RU_ml_fr_indv_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_set)),2, mean)
AU_ml_fr_indv_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_set)),2, mean)

value_ml_fr_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$value_ml_fr_indv_cluster_set)),2, mean)
RU_ml_fr_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_cluster_set)),2, mean)
AU_ml_fr_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_cluster_set)),2, mean)

FURG_indv_set_mean <- apply(t(sapply(results, function(x) x$FURG_indv_set)),2, mean)
FUTR_indv_set_mean <- apply(t(sapply(results, function(x) x$FUTR_indv_set)),2, mean)
FURG_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$FURG_indv_cluster_set)),2, mean)
FUTR_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$FUTR_indv_cluster_set)),2, mean)

RU_ml_fr_indv_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_CATE_set)),2, mean)
AU_ml_fr_indv_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_CATE_set)),2, mean)

RU_ml_fr_indv_cluster_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_cluster_CATE_set)),2, mean)
AU_ml_fr_indv_cluster_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_cluster_CATE_set)),2, mean)

# Combine all the means into one dataframe
means_df <- data.frame(value_ml_fr_indv_set_mean, RU_ml_fr_indv_set_mean, AU_ml_fr_indv_set_mean,
                       value_ml_fr_indv_cluster_set_mean, RU_ml_fr_indv_cluster_set_mean, AU_ml_fr_indv_cluster_set_mean,
                       FURG_indv_set_mean, FUTR_indv_set_mean, FURG_indv_cluster_set_mean, FUTR_indv_cluster_set_mean,
                       RU_ml_fr_indv_CATE_set_mean, AU_ml_fr_indv_CATE_set_mean,
                       RU_ml_fr_indv_cluster_CATE_set_mean, AU_ml_fr_indv_cluster_CATE_set_mean)
head(means_df)

# ---------------------------------------------------------------
#                    Design 1: 
#               Plot 1: RU vs AU
#
# ---------------------------------------------------------------

# Define plot parameters
line_thick <- 1.5
x_max <- max(AU_ml_fr_indv_set_mean, AU_ml_fr_indv_cluster_set_mean)
y_max <- round(max(RU_ml_fr_indv_set_mean, RU_ml_fr_indv_cluster_set_mean), 3)
y_min <- round(min(RU_ml_fr_indv_set_mean, RU_ml_fr_indv_cluster_set_mean), 3)
y_interval <- round((y_max - y_min) / 5, 3)
x_interval <- round(x_max / 8, 2)

# Create data frames for plotting
data1 <- data.frame(AU = AU_ml_fr_indv_set_mean, 
                    RU = RU_ml_fr_indv_set_mean, 
                    Group = "Individual")
data2 <- data.frame(AU = AU_ml_fr_indv_cluster_set_mean, 
                    RU = RU_ml_fr_indv_cluster_set_mean,
                    Group = "Individual + Cluster")
data <- rbind(data1, data2)

# Create the plot
p <- ggplot(data, aes(x = AU, y = RU, color = Group, shape = Group)) +
  # Add the lines
  geom_line(size = line_thick) +
  # Add the points
  geom_point(size = 3) +
  # Custom axis labels and limits
  labs(x = "Mean Unfairness", y = "Relative Utility") +
  scale_x_continuous(limits = c(0.00, x_max), 
                     breaks = seq(0, x_max, by = x_interval)) +
  scale_y_continuous(limits = c(y_min - 2*y_interval, y_max + 2*y_interval), 
                     breaks = seq(y_min - 2*y_interval, y_max + 2*y_interval, y_interval)) +
  # Add text annotations
  annotate("text", 
           x = AU_ml_fr_indv_cluster_set_mean[1] + 0.035, 
           y = RU_ml_fr_indv_cluster_set_mean[1], 
           label = bquote(delta == .("0.0001")), 
           color = "black", hjust = 0,  alpha = 0.28, size = 6) +
  annotate("text", 
           x = AU_ml_fr_indv_set_mean[length(AU_ml_fr_indv_set_mean) - 1], 
           y = RU_ml_fr_indv_set_mean[length(RU_ml_fr_indv_set_mean) - 1] + 0.001, 
           label = expression(delta == infinity), 
           color = "black", hjust = 1, alpha = 0.6, size = 6) +
  # Add custom grid lines
  theme_minimal(base_size = 18) +  # Set the base font size for the entire plot
  theme(
    panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted"),
    legend.position = c(0.9985, 0.0795),
    legend.justification = "right",
    legend.background = element_rect(color = "black", size=0.6, fill = alpha('white', 0.7)),
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border
  ) +
  # Customize legend labels
  scale_color_manual(values = c("Individual" = "blue", "Individual + Cluster" = "red"), 
                     labels = c(expression(uf[1]~with~S[1*ij]^{(1)}), 
                                expression(uf[1]~with~S[1*ij]^{(1)}~and~uf[2]~with~S[1*j]^{(2)}))) +
  scale_shape_manual(values = c("Individual" = 16, "Individual + Cluster" = 17),
                     guide = "none") +  # Hide the shape legend
  guides(color = guide_legend(override.aes = list(shape = c(16, 17))))



ggsave(
  filename = "results/Design_1_plot_1.jpeg",  # File name
  plot = p,  # The plot object to save
  width = 593 * 3.5 / 300,  # Convert width to inches (pixels / DPI)
  height = 639 * 3 / 300,  # Convert height to inches (pixels / DPI)
  dpi = 300  # Resolution in dots per inch (DPI)
)

# ---------------------------------------------------------------
#                    Design 1: 
#               Plot 2: MSE vs Difference in CATE
#
# ---------------------------------------------------------------

x_max <- max(AU_ml_fr_indv_CATE_set_mean, AU_ml_fr_indv_cluster_CATE_set_mean)
y_max <- round(max(RU_ml_fr_indv_CATE_set_mean, RU_ml_fr_indv_cluster_CATE_set_mean),3)
y_min <- round(min(RU_ml_fr_indv_CATE_set_mean, RU_ml_fr_indv_cluster_CATE_set_mean),3)
y_interval <- round((y_max - y_min) / 5,3)
x_interval <- round(x_max / 5,1)

# Create data frames for plotting
data1 <- data.frame(AU = AU_ml_fr_indv_CATE_set_mean, 
                    RU = RU_ml_fr_indv_CATE_set_mean, 
                    Group = "Individual")
data2 <- data.frame(AU = AU_ml_fr_indv_cluster_CATE_set_mean,
                    RU = RU_ml_fr_indv_cluster_CATE_set_mean,
                    Group = "Individual + Cluster")

data <- rbind(data1, data2)

# Create the plot
p2 <- ggplot(data, aes(x = AU, y = RU, color = Group, shape = Group)) +
  # Add the lines
  geom_line(size = line_thick) +
  # Add the points
  geom_point(size = 3) +
  # Custom axis labels and limits
  labs(x = "Mean Unfairness", y = "Mean Squared Error") +
  scale_x_continuous(limits = c(0.00, x_max), 
                     breaks = seq(0, x_max, by = x_interval)) +
  scale_y_continuous(limits = c(y_min - 2*y_interval, y_max + 2*y_interval), 
                     breaks = seq(y_min - 2*y_interval, y_max + 2*y_interval, y_interval)) +
  # Add text annotations
  annotate("text", 
           x = AU_ml_fr_indv_cluster_CATE_set_mean[1] + 0.1, 
           y = RU_ml_fr_indv_cluster_CATE_set_mean[1], 
           label = bquote(delta == .("0.0001")), 
           color = "black", hjust = 0,  alpha = 0.28, size = 6) +
  annotate("text", 
           x = AU_ml_fr_indv_CATE_set_mean[length(AU_ml_fr_indv_CATE_set_mean) - 1], 
           y = RU_ml_fr_indv_CATE_set_mean[length(RU_ml_fr_indv_CATE_set_mean) - 1] + 0.03, 
           label = expression(delta == infinity), 
           color = "black", hjust = 1, alpha = 0.6, size = 6) +
  # Add custom grid lines
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted"),
    legend.position = c(0.9985, 0.079),
    legend.justification = "right",
    legend.background = element_rect(color = "black", size=0.6, fill = alpha('white', 0.7)),
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border
  ) +
  # Customize legend labels
  scale_color_manual(values = c("Individual" = "blue", "Individual + Cluster" = "red"), 
                     labels = c(expression(uf[1]~with~S[1*ij]^{ (1)}), 
                                expression(uf[1]~with~S[1*ij]^{ (1) }~and~uf[2]~with~S[1*j]^{ (2) }))) +
  scale_shape_manual(values = c("Individual" = 16, "Individual + Cluster" = 17),
                     guide = "none") +  # Hide the shape legend
  guides(color = guide_legend(override.aes = list(shape = c(16, 17))))

print(p2)

ggsave(
  filename = "results/Design_1_plot_2.jpeg",  # File name
  plot = p2,  # The plot object to save
  width = 593 * 3.5 / 300,  # Convert width to inches (pixels / DPI)
  height = 639 * 3 / 300,  # Convert height to inches (pixels / DPI)
  dpi = 300  # Resolution in dots per inch (DPI)
)


# ---------------------------------------------------------------
#                    Design 2: 
#               Plot 1: RU vs AU
#
# ---------------------------------------------------------------

# remove all the variables from the environment
rm(list = ls())
gc()


load("results/Simulation_Design_2_results_20_reps.rda")

list_to_analze <- c("value_ml_fr_indv_set", "RU_ml_fr_indv_set", "AU_ml_fr_indv_set",
                    "value_ml_fr_indv_cluster_set", "RU_ml_fr_indv_cluster_set", "AU_ml_fr_indv_cluster_set",
                    "FURG_indv_set", "FUTR_indv_set", "FURG_indv_cluster_set", "FUTR_indv_cluster_set",
                    "value_ml_fr_is_set","RU_ml_fr_is_set","AU_ml_fr_is_set","FURG_is_set","FUTR_is_set",
                    "RU_ml_fr_indv_CATE_set", "AU_ml_fr_indv_CATE_set",
                    "RU_ml_fr_indv_cluster_CATE_set", "AU_ml_fr_indv_cluster_CATE_set",
                    "RU_ml_fr_is_CATE_set", "AU_ml_fr_is_CATE_set")

# Define the names of the elements you want to convert into a dataframe
elements_to_convert <- list_to_analze

# Initialize a list to store the extracted elements
extracted_elements <- list()

# Extract each element from all list elements
for (element in elements_to_convert) {
  # Use sapply to extract the element from each list and store as a column in a dataframe
  extracted_elements[[element]] <- t(sapply(results, function(x) x[[element]]))
}

# extract each element from the list and convert it to a dataframe, and calculate the column mean
value_ml_fr_indv_set_mean <- apply(t(sapply(results, function(x) x$value_ml_fr_indv_set)),2, mean)
RU_ml_fr_indv_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_set)),2, mean)
AU_ml_fr_indv_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_set)),2, mean)

value_ml_fr_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$value_ml_fr_indv_cluster_set)),2, mean)
RU_ml_fr_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_cluster_set)),2, mean)
AU_ml_fr_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_cluster_set)),2, mean)

FURG_indv_set_mean <- apply(t(sapply(results, function(x) x$FURG_indv_set)),2, mean)
FUTR_indv_set_mean <- apply(t(sapply(results, function(x) x$FUTR_indv_set)),2, mean)
FURG_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$FURG_indv_cluster_set)),2, mean)
FUTR_indv_cluster_set_mean <- apply(t(sapply(results, function(x) x$FUTR_indv_cluster_set)),2, mean)

value_ml_fr_is_set_mean <- apply(t(sapply(results, function(x) x$value_ml_fr_is_set)),2, mean)
RU_ml_fr_is_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_is_set)),2, mean)
AU_ml_fr_is_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_is_set)),2, mean)
FURG_is_set_mean <- apply(t(sapply(results, function(x) x$FURG_is_set)),2, mean)
FUTR_is_set_mean <- apply(t(sapply(results, function(x) x$FUTR_is_set)),2, mean)

RU_ml_fr_indv_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_CATE_set)),2, mean)
AU_ml_fr_indv_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_CATE_set)),2, mean)

RU_ml_fr_indv_cluster_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_cluster_CATE_set)),2, mean)
AU_ml_fr_indv_cluster_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_cluster_CATE_set)),2, mean)

RU_ml_fr_is_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_is_CATE_set)),2, mean)
AU_ml_fr_is_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_is_CATE_set)),2, mean)

# create the dataframe
data1 <- data.frame(AU = AU_ml_fr_indv_set_mean, 
                    RU = RU_ml_fr_indv_set_mean, 
                    Group = "Individual")

data2 <- data.frame(AU = AU_ml_fr_indv_cluster_set_mean,
                    RU = RU_ml_fr_indv_cluster_set_mean,
                    Group = "Individual + Cluster")

data3 <- data.frame(AU = AU_ml_fr_is_set_mean,
                    RU = RU_ml_fr_is_set_mean,
                    Group = "Intersectional")

data <- rbind(data1, data2, data3)

# define the range of x-y axis
x_max <- max(AU_ml_fr_indv_set_mean, AU_ml_fr_indv_cluster_set_mean, AU_ml_fr_is_set_mean)
y_max <- round(max(RU_ml_fr_indv_set_mean, RU_ml_fr_indv_cluster_set_mean, RU_ml_fr_is_set_mean),3)
y_min <- round(min(RU_ml_fr_indv_set_mean, RU_ml_fr_indv_cluster_set_mean, RU_ml_fr_is_set_mean),3)
y_interval <- round((y_max - y_min) / 5,3)
x_interval <- round(x_max / 5,2)

line_thick <- 1.5
# Create the plot
p3 <- ggplot(data, aes(x = AU, y = RU, color = Group, shape = Group)) +
  # Add the lines
  geom_line(size = line_thick) +
  # Add the points
  geom_point(size = 3) +
  # Custom axis labels and limits
  labs(x = "Mean Unfairness", y = "Relative Utility") +
  scale_x_continuous(limits = c(0.00, x_max), 
                     breaks = seq(0, x_max, by = x_interval)) +
  scale_y_continuous(limits = c(y_min - 2*y_interval, y_max + 2*y_interval), 
                     breaks = seq(y_min - 2*y_interval, y_max + 2*y_interval, y_interval)) +
  # Add text annotations
  annotate("text", 
           x = AU_ml_fr_indv_cluster_set_mean[1]-0.025, 
           y = RU_ml_fr_indv_cluster_set_mean[1]-0.0015, 
           label = bquote(delta == .("0.0001")), 
           color = "black", hjust = 0,  alpha = 0.28, size = 6)  +
  annotate("text", 
           x = AU_ml_fr_indv_set_mean[length(AU_ml_fr_indv_set_mean) - 1]-0.025, 
           y = RU_ml_fr_indv_set_mean[length(RU_ml_fr_indv_set_mean) - 1] + 0.0015, 
           label = expression(delta == infinity), 
           color = "black", hjust = 0,  alpha = 0.6, size = 6)  +
  # Add custom grid lines
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted"),
    legend.position = c(0.998, 0.1105),
    legend.justification = "right",
    legend.background = element_rect(color = "black", size=0.6, fill = alpha('white', 0.7)),
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border
  ) +
  # Customize legend labels
  scale_color_manual(values = c("Individual" = "blue", "Individual + Cluster" = "red", "Intersectional" = "black"), 
                     labels = c(expression(uf[1]~with~S[1*ij]^{ (1)}), 
                                expression(uf[1]~with~S[1*ij]^{ (1) }~and~uf[2]~with~S[1*j]^{ (2) }),
                                expression(uf[k]~with~(S[1*ij]^{ (1) }*","~S[1*j]^{ (2)})~"for"~k==1*","*2*","*3))) +
  scale_shape_manual(values = c("Individual" = 16, "Individual + Cluster" = 17, "Intersectional" = 18),
                     guide = "none") +  # Hide the shape legend
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 18))))



ggsave(
  filename = "results/Design_2_plot_1.jpeg",  # File name
  plot = p3,  # The plot object to save
  width = 593 * 3.5 / 300,  # Convert width to inches (pixels / DPI)
  height = 639 * 3 / 300,  # Convert height to inches (pixels / DPI)
  dpi = 300  # Resolution in dots per inch (DPI)
)

# ---------------------------------------------------------------
#                    Design 2: 
#               Plot 2: MSE vs. Difference in CATE
#
# ---------------------------------------------------------------

x_max <- max(AU_ml_fr_indv_CATE_set_mean, 
             AU_ml_fr_indv_cluster_CATE_set_mean,
             AU_ml_fr_is_CATE_set_mean)
y_max <- round(max(RU_ml_fr_indv_CATE_set_mean, 
                   RU_ml_fr_indv_cluster_CATE_set_mean,
                   RU_ml_fr_is_CATE_set_mean),3)
y_min <- round(min(RU_ml_fr_indv_CATE_set_mean, 
                   RU_ml_fr_indv_cluster_CATE_set_mean,
                   RU_ml_fr_is_CATE_set_mean),3)
y_interval <- round((y_max - y_min) / 5,3)
x_interval <- round(x_max / 5,2)

# Create data frames for plotting
data1 <- data.frame(AU = AU_ml_fr_indv_CATE_set_mean, 
                    RU = RU_ml_fr_indv_CATE_set_mean, 
                    Group = "Individual")
data2 <- data.frame(AU = AU_ml_fr_indv_cluster_CATE_set_mean,
                    RU = RU_ml_fr_indv_cluster_CATE_set_mean,
                    Group = "Individual + Cluster")
data3 <- data.frame(AU = AU_ml_fr_is_CATE_set_mean,
                    RU = RU_ml_fr_is_CATE_set_mean,
                    Group = "Intersectional")

data <- rbind(data1, data2, data3)

# Create the plot
p4 <- ggplot(data, aes(x = AU, y = RU, color = Group, shape = Group)) +
  # Add the lines
  geom_line(size = line_thick) +
  # Add the points
  geom_point(size = 3) +
  # Custom axis labels and limits
  labs(x = "Mean Unfairness", y = "Mean Squared Error") +
  scale_x_continuous(limits = c(0.00, x_max), 
                     breaks = seq(0, x_max, by = x_interval)) +
  scale_y_continuous(limits = c(y_min - 2*y_interval, y_max + 2*y_interval), 
                     breaks = seq(y_min - 2*y_interval, y_max + 2*y_interval, y_interval)) +
  # Add text annotations
  annotate("text", 
           x = AU_ml_fr_indv_cluster_CATE_set_mean[1] - 0.06, 
           y = RU_ml_fr_indv_cluster_CATE_set_mean[1]  + 0.03, 
           label = bquote(delta == .("0.0001")), 
           color = "black", hjust = 0,  alpha = 0.28, size = 6)  +
  annotate("text", 
           x = AU_ml_fr_indv_CATE_set_mean[length(AU_ml_fr_indv_CATE_set_mean) - 1] - 0.06, 
           y = RU_ml_fr_indv_CATE_set_mean[length(RU_ml_fr_indv_CATE_set_mean) - 1] + 0.03, 
           label = expression(delta == infinity), 
           color = "black", hjust = 0,  alpha = 0.6, size = 6)  +
  # Add custom grid lines
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted"),
    legend.position = c(0.998, 0.1105),
    legend.justification = "right",
    legend.background = element_rect(color = "black", size=0.6,
                                     fill = alpha('white', 0.7)),
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border
  ) +
  # Customize legend labels
  scale_color_manual(values = c("Individual" = "blue", 
                                "Individual + Cluster" = "red", 
                                "Intersectional" = "black"), 
                     labels = c(expression(uf[1]~with~S[1*ij]^{ (1)}), 
                                expression(uf[1]~with~S[1*ij]^{ (1) }~and~uf[2]~with~S[1*j]^{ (2) }),
                                expression(uf[k]~with~(S[1*ij]^{ (1) }*","~S[1*j]^{ (2)})~"for"~k==1*","*2*","*3))) +
  scale_shape_manual(values = c("Individual" = 16,
                                "Individual + Cluster" = 17,
                                "Intersectional" = 18),
                     guide = "none") +  # Hide the shape legend
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 18))))

print(p4)

ggsave(
  filename = "results/Design_2_plot_2.jpeg",  # File name
  plot = p4,  # The plot object to save
  width = 593 * 3.5 / 300,  # Convert width to inches (pixels / DPI)
  height = 639 * 3 / 300,  # Convert height to inches (pixels / DPI)
  dpi = 300  # Resolution in dots per inch (DPI)
)