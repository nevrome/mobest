vars <- setdiff(names(iris), "Species")

pageWithSidebar(
  NULL,
  sidebarPanel(
    numericInput("kernel_spatial_size", "Kernel: spatial lengthscale parameter in [m]", 500000),
    numericInput("kernel_temporal_size", "Kernel: temporal lengthsclae parameter in [a]", 1000),
    numericInput("kernel_nugget", "Kernel: nugget", 0.1),
    numericInput("pred_grid_spatial_cell_size", "Prediction grid: spatial cell size in [m]", 200000),
    numericInput("pred_grid_temporal_distance", "Prediction grid: temporal distance in [a]", 500)
  ),
  mainPanel(
    plotOutput('plot1', height = "1000px")
  )
)
