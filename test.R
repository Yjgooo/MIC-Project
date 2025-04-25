


source("simulation_functions/MIC_utils.R")
source("simulation_functions/np_simulation.R")
source("simulation_functions/ph_simulation.R")
source("simulation_functions/po_simulation.R")
source("simulation_functions/real_simulation.R")
source("simulation_functions/generic_simulation.R")


args = commandArgs(T)


setting <- as.numeric(args[1])
algorithm <- as.numeric(args[2])

sim_models <- c("simulate_data_po", "simulate_data_ph", "simulate_data_np", "simulate_data_real")
models <- c("np", "ph", "po", "roc", "pred", "pred.adj")
n_values <- c(100, 500, 1000)

sink(file = paste0("/mnt/home/gongxi/my-project/MIC_simulation", setting, "_Algorithm", algorithm, ".text"))


for (n_val in n_values) {
  for (j in seq_along(models)) {
    sim_model <- sim_models[setting]
    model <- models[j]
    
    # Call the generic simulation function
    sim_result <- generic_sim(n = n_val, m = m, setup = setup, cov = cov,
                              sim_model = sim_model, model = model, boot = boot)
    
    # Create a key for the result list, e.g. "n_100_po"
    key <- paste("n", n_val, sim_model, "fit_model", model, sep = "_")
    results_list[[key]] <- sim_result 
  }
}

sink(file = paste0("mnt/home/gongxi/my-project/MIC_simulation/HPCC_result", sim_models[setting], "_Algorithm", models[algorithm], ".text"))
print(sim_result)
sink(file = NULL)

