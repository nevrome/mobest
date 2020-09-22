library(magrittr)
library(ggplot2)

print(getwd())

load("../../coest.interpol.2020/data/poseidon_data/janno_final.RData")
load("../../coest.interpol.2020/data/spatial/area.RData")
load("../../coest.interpol.2020/data/spatial/mobility_regions.RData")

function(input, output, session) {

  # GPR
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
        shiny_kernel = list(d = c(
            input$kernel_spatial_size,
            input$kernel_spatial_size,
            input$kernel_temporal_size
          ), g = input$kernel_nugget,
          on_residuals = T, auto = F
        )
      ),
      prediction_grid = list(
        scs100_tl100 = mobest::create_prediction_grid(
          area,
          mobility_regions,
          spatial_cell_size = input$pred_grid_spatial_cell_size,
          time_layers = seq(-7500, 1500, input$pred_grid_temporal_distance)
        )
      )
    )

    withProgress(message = "GPR", {
      interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)
    })

    interpol_grid
  })

  # Spatial origin
  mobility_clemens <- reactive({
    withProgress(message = "Spatial origin search", {
      origin_grid <- mobest::search_spatial_origin(interpol_grid(), steps = 1)
    })
    mobest::estimate_mobility(origin_grid)
  })


  # dynamic inputs
  output$time_slider_input <- renderUI({
    sliderInput(
      "plot_z", "Plot: time to show in [a calBC/AD]", -7500, 1500, -3000,
      step = input$pred_grid_temporal_distance
    )
  })

  # plots
  output$plot1 <- renderPlot({

    if (input$plot_type == "C1_comic") {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C1") %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c() +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "C2_comic") {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C2") %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "plasma") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "C1") {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C1", z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c() +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "C2") {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C2", z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "plasma") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "mobility_clemens_comic") {

      mobility_clemens() %>%
        ggplot() +
        geom_raster(aes(x, y, fill = speed_km_per_decade)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "cividis") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "mobility_clemens") {

      mobility_clemens() %>% dplyr::filter(z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = speed_km_per_decade)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "cividis") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    }

  })

}
