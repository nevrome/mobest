vars <- setdiff(names(iris), "Species")

pageWithSidebar(
  NULL,
  sidebarPanel(
    numericInput("kernel_spatial_size", "Kernel: spatial lengthscale parameter in [m]", 500000),
    numericInput("kernel_temporal_size", "Kernel: temporal lengthscale parameter in [a]", 1000),
    numericInput("kernel_nugget", "Kernel: nugget", 0.1),
    hr(style="border-color: black;"),
    numericInput("pred_grid_spatial_cell_size", "Prediction grid: spatial cell size in [m]", 200000),
    numericInput("pred_grid_temporal_distance", "Prediction grid: temporal distance in [a]", 500),
    hr(style="border-color: black;"),
    selectInput("mobility_algorithm", "Mobility algorithm", choices = c("clemens", "stephan")),
    sliderInput("clemens_mobility_steps", "Clemens' algorithm: steps in the past", 1, 10, 1),
    numericInput("stephan_mobility_delta_x", "Stephan's algorithm: delta_x in [m]", 10000),
    numericInput("stephan_mobility_delta_y", "Stephan's algorithm: delta_y in [m]", 10000),
    numericInput("stephan_mobility_delta_z", "Stephan's algorithm: delta_z in [a]", 10),
    hr(style="border-color: black;"),
    selectInput("plot_type", "Plot type", choices = c(
      "C1_comic", "C2_comic", "C1", "C2",
      "mobility_clemens_comic", "mobility_clemens",
      "mobility_stephan_comic", "mobility_stephan"
    )),
    uiOutput("time_slider_input")
  ),
  mainPanel(
    plotOutput('plot1', height = "1000px")
  )
)
