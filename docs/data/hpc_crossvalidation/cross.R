library(magrittr)

args <- unlist(strsplit(commandArgs(trailingOnly = TRUE), " "))
run <- args[1]
ds_for_this_run <- as.numeric(args[2])
dt_for_this_run <- as.numeric(args[3])

samples_projected <- readr::read_csv(
  "docs/data/samples_projected.csv"
)

ind <- mobest::create_spatpos(
  id = samples_projected$Sample_ID,
  x  = samples_projected$x,
  y  = samples_projected$y,
  z  = samples_projected$Date_BC_AD_Median
)
dep <- mobest::create_obs(
  C1 = samples_projected$MDS_C1,
  C2 = samples_projected$MDS_C2
)

kernel_for_this_run <- mobest::create_kernset_multi(
  mobest::create_kernset(
    C1 = mobest::create_kernel(
      dsx = ds_for_this_run*1000,
      dsy = ds_for_this_run*1000,
      dt  = dt_for_this_run,
      g   = 0.071
    ),
    C2 = mobest::create_kernel(
      dsx = ds_for_this_run*1000,
      dsy = ds_for_this_run*1000,
      dt  = dt_for_this_run,
      g   = 0.059
    )
  ),
  .names = paste0("kernel_", run)
)

set.seed(123)

interpol_comparison <- mobest::crossvalidate(
  independent = ind,
  dependent   = dep,
  kernel      = kernel_for_this_run,
  iterations  = 10,
  groups      = 10,
  quiet       = F
)

kernel_grid <- interpol_comparison %>%
  dplyr::group_by(
    dependent_var_id, ds = dsx, dt) %>%
  dplyr::summarise(
    mean_squared_difference = mean(difference^2),
    .groups = "drop"
  )

readr::write_csv(
  kernel_grid,
  file = paste0(
    "docs/data/hpc_crossvalidation/kernel_grid_",
    sprintf("%06d", as.integer(run)),
    ".csv"
  )
)
