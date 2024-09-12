# load the results

design_1 <- FALSE

load("results/Simulation_Design_2_results.rda")

length(results)

results[[2]]

# Assuming your list is named 'results'
# Define the names of the elements you want to average
# PAN: 2024.08.03
# add the CATE unfairness and MSE cate
elements_to_average <- c("value_true", "value_random", "value_BART", "value_cf", 
                         "RU_true", "RU_BART", "RU_cf","RU_cf_ps", "average_unfairness_true",
                         "average_unfairness_bart","average_unfairness_cf","average_unfairness_cf_ps",
                         "RU_BART_CATE", "RU_cf_CATE", "RU_cf_CATE_ps",
                         "average_unfairness_bart_CATE","average_unfairness_cf_CATE","average_unfairness_cf_CATE_ps")

# Initialize a list to store the average values
average_values <- list()

# Loop through each element name and calculate its average value
for (element in elements_to_average) {
  average_values[[element]] <- mean(sapply(results, function(x) x[[element]]))
}

# Convert the list to a named vector for better readability
average_values <- unlist(average_values)

# Print the average values
print(round(average_values,3))

# ----------------------------------
# extract the rest lists in the results
# ----------------------------------
if (design_1) {
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
  print(round(average_values_delta20,3))
  
  
  # ----------------------------------
  # Convert the list-elements into one dataframe
  # ----------------------------------
  
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
  
}else {
  
  list_to_analze <- c("value_ml_fr_indv_set", "RU_ml_fr_indv_set", "AU_ml_fr_indv_set",
                      "value_ml_fr_indv_cluster_set", "RU_ml_fr_indv_cluster_set", "AU_ml_fr_indv_cluster_set",
                      "FURG_indv_set", "FUTR_indv_set", "FURG_indv_cluster_set", "FUTR_indv_cluster_set",
                      "value_ml_fr_is_set","RU_ml_fr_is_set","AU_ml_fr_is_set","FURG_is_set","FUTR_is_set",
                      "RU_ml_fr_indv_CATE_set", "AU_ml_fr_indv_CATE_set",
                      "RU_ml_fr_indv_cluster_CATE_set", "AU_ml_fr_indv_cluster_CATE_set",
                      "RU_ml_fr_is_CATE_set", "AU_ml_fr_is_CATE_set")
  
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
  print(round(average_values_delta20,3))
  
  
  # ----------------------------------
  # Convert the list-elements into one dataframe
  # ----------------------------------
  
  # Define the names of the elements you want to convert into a dataframe
  elements_to_convert <- list_to_analze
  
  # Initialize a list to store the extracted elements
  extracted_elements <- list()
  
  # Extract each element from all list elements
  for (element in elements_to_convert) {
    # Use sapply to extract the element from each list and store as a column in a dataframe
    extracted_elements[[element]] <- t(sapply(results, function(x) x[[element]]))
  }
  
  names(extracted_elements)
  class(extracted_elements[[1]])
  
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
  
  # PAN:2024.08.04
  RU_ml_fr_indv_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_CATE_set)),2, mean)
  AU_ml_fr_indv_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_CATE_set)),2, mean)
  
  RU_ml_fr_indv_cluster_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_indv_cluster_CATE_set)),2, mean)
  AU_ml_fr_indv_cluster_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_indv_cluster_CATE_set)),2, mean)
  
  RU_ml_fr_is_CATE_set_mean <- apply(t(sapply(results, function(x) x$RU_ml_fr_is_CATE_set)),2, mean)
  AU_ml_fr_is_CATE_set_mean <- apply(t(sapply(results, function(x) x$AU_ml_fr_is_CATE_set)),2, mean)
  
  
  # Combine all the means into one dataframe
  means_df <- data.frame(value_ml_fr_indv_set_mean, RU_ml_fr_indv_set_mean, AU_ml_fr_indv_set_mean,
                         value_ml_fr_indv_cluster_set_mean, RU_ml_fr_indv_cluster_set_mean, AU_ml_fr_indv_cluster_set_mean,
                         FURG_indv_set_mean, FUTR_indv_set_mean, FURG_indv_cluster_set_mean, FUTR_indv_cluster_set_mean,
                         value_ml_fr_is_set_mean, RU_ml_fr_is_set_mean, AU_ml_fr_is_set_mean, FURG_is_set_mean, FUTR_is_set_mean,
                         RU_ml_fr_indv_CATE_set_mean, AU_ml_fr_indv_CATE_set_mean,
                         RU_ml_fr_indv_cluster_CATE_set_mean, AU_ml_fr_indv_cluster_CATE_set_mean,
                         RU_ml_fr_is_CATE_set_mean, AU_ml_fr_is_CATE_set_mean)
  head(means_df)
  
}