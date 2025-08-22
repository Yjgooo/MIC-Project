library(here) # here() always resolves paths relative to the project root (where your .Rproj file lives).

source(here("R", "simulation_functions", "MIC_utils.R"))
source(here("R", "simulation_functions", "np_simulation.R"))
source(here("R", "simulation_functions", "ph_simulation.R"))
source(here("R", "simulation_functions", "po_simulation.R"))
source(here("R", "simulation_functions", "real_simulation.R"))
source(here("R", "simulation_functions", "generic_simulation.R"))



#system.time({
  print("entered")
  p_result <- run_simulations(setup = 0, m = 100, boot = FALSE)
  print(p_result)
  print("what?")
#})



print("done")


