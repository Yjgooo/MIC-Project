# Load functions
source("simulation_functions/MIC_utils.R")
source("simulation_functions/np_simulation.R")
source("simulation_functions/ph_simulation.R")
source("simulation_functions/po_simulation.R")
source("simulation_functions/real_simulation.R")
source("simulation_functions/generic_simulation.R")

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
setting <- as.numeric(args[1])
algorithm <- as.numeric(args[2])

# Define model names
sim_models <- c("simulate_data_po", "simulate_data_ph", "simulate_data_np", "simulate_data_real")
models <- c("np", "ph", "po", "roc", "pred", "pred.adj")
n_values <- c(100, 500, 1000)

results_list <- list()

# Main simulation loop
for (n_val in n_values) {
  sim_model <- sim_models[setting]
  model <- models[algorithm]
  
  sim_result <- generic_sim(n = n_val, m = 100, setup = 0, cov = NULL,
                            sim_model = sim_model, model = model, boot = TRUE)
  
  key <- paste("n", n_val, sim_model, "fit_model", model, sep = "_")
  results_list[[key]] <- sim_result 
}

# Prepare output directory
output_dir <- "/mnt/home/gongxi/my-project/MIC_simulation/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare output filename
output_file <- file.path(output_dir, paste0(sim_models[setting], "_Algorithm_", models[algorithm], ".txt"))

# Write output
sink(output_file)
print(results_list)
sink()   # just sink(), not sink(file = NULL)

