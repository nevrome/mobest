library(magrittr)
library(ggplot2)

print(getwd())

load("../../coest.interpol.2020/data/poseidon_data/janno_final.RData")
load("../../coest.interpol.2020/data/spatial/area.RData")
load("../../coest.interpol.2020/data/spatial/mobility_regions.RData")

function(input, output, session) {

  # Combine the selected variables into a new data frame
  interpol_grid <- reactive({

    model_grid <- mobest::create_model_grid(
      independent = list(
        tibble::tibble(
          x = janno_final$x,
          y = janno_final$y,
          z = janno_final$Date_BC_AD_Median_Derived
        )
      ) %>% stats::setNames("age_median"),
      dependent = list(
        C1 = janno_final$C1,
        C2 = janno_final$C2
      ),
      kernel = list(
        ds600_dt300_g01 = list(d = c(500000, 500000, 1000), g = 0.1, on_residuals = T, auto = F)
      ),
      prediction_grid = list(
        scs100_tl100 = mobest::create_prediction_grid(
          area,
          mobility_regions,
          spatial_cell_size = 200000,
          time_layers = seq(-7500, 1500, 500)
        )
      )
    )

    interpol_grid <- mobest::run_model_grid(model_grid)

    interpol_grid
  })

  # clusters <- reactive({
  #   kmeans(selectedData(), input$clusters)
  # })

  output$plot1 <- renderPlot({

    interpol_grid() %>%
      dplyr::filter(
        #kernel_setting_id == "ds400_dt700_g001",
        dependent_var_id == "C1",
        z %in% seq(-7000, -2000, 500)
        #z %% 500 == 0
      ) %>%
      ggplot() +
      geom_raster(aes(x, y, fill = mean))+#, alpha = sd)) +
      facet_wrap(~z) +
      scale_fill_viridis_c() +
      scale_alpha_continuous(range = c(1, 0), na.value = 0)

  })

}
