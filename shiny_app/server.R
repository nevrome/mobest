library(magrittr)
library(ggplot2)

print(getwd())

load("../../coest.interpol.2020/data/poseidon_data/janno_final.RData")
load("../../coest.interpol.2020/data/spatial/area.RData")
load("../../coest.interpol.2020/data/spatial/mobility_regions.RData")

function(input, output, session) {

  # dynamic inputs
  output$time_slider_input <- renderUI({
    if (!input$comic) {
      sliderInput(
        "plot_z", "Plot: time to show in [a calBC/AD]", -7500, 1500, -3000,
        step = input$pred_grid_temporal_distance
      )
    } else {
      HTML("Comic mode active")
    }
  })

  # prediction grid preparation
  pred_grids <- reactive({

    main_pred_grid <- mobest::create_prediction_grid(
      area,
      mobility_regions,
      spatial_cell_size = input$pred_grid_spatial_cell_size,
      time_layers = seq(-7500, 1500, input$pred_grid_temporal_distance)
    )

    if (input$mobility_algorithm == "clemens") {
      list(
        main = main_pred_grid
      )
    } else if (input$mobility_algorithm == "stephan") {
      list(
        main = main_pred_grid,
        offset_x = main_pred_grid %>% dplyr::mutate(x = x + input$stephan_mobility_delta_x),
        offset_y = main_pred_grid %>% dplyr::mutate(y = y + input$stephan_mobility_delta_y),
        offset_z = main_pred_grid %>% dplyr::mutate(z = z + input$stephan_mobility_delta_z)
      )
    }

  })

  # GPR
  interpol_grid <- reactive({

    req(pred_grids())

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
      prediction_grid = pred_grids()
    )

    withProgress(message = "GPR", {
      interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)
    })

    interpol_grid
  })

  # mobility calculation
  mobility_clemens <- reactive({
    # hacky trick to reduce calculcation time
    if (input$mobility_algorithm == "stephan") {
      interpol_grid <- interpol_grid() %>% dplyr::filter(pred_grid_id == "main")
    } else {
      interpol_grid <- interpol_grid()
    }
    # updateSelectInput(session, "mobility_algorithm", selected = "clemens")
    withProgress(message = "Mobility estimation (Clemens)", {
      origin_grid <- mobest::search_spatial_origin(interpol_grid, steps = input$clemens_mobility_steps)
      mobility <- mobest::estimate_mobility(origin_grid)
    })
    mobility
  })

  mobility_stephan <- reactive({
    updateSelectInput(session, "mobility_algorithm", selected = "stephan")
    withProgress(message = "Mobility estimation (Stephan)", {
      mobility <- mobest::estimate_mobility(
        interpol_grid(),
        input$stephan_mobility_delta_x,
        input$stephan_mobility_delta_y,
        input$stephan_mobility_delta_z
      )
    })
    mobility
  })

  # plots
  output$plot1 <- renderPlot({

    if (input$plot_type == "C1" && input$comic) {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C1", pred_grid_id == "main") %>%
        ggplot() +
        {
          if (input$sd_as_alpha) {
            geom_raster(aes(x, y, fill = mean, alpha = sd))
          } else {
            geom_raster(aes(x, y, fill = mean))
          }
        } +
        facet_wrap(~z) +
        scale_fill_viridis_c() +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "C2" && input$comic) {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C2", pred_grid_id == "main") %>%
        ggplot() +
        {
          if (input$sd_as_alpha) {
            geom_raster(aes(x, y, fill = mean, alpha = sd))
          } else {
            geom_raster(aes(x, y, fill = mean))
          }
        } +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "plasma") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "C1" && !input$comic) {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C1", pred_grid_id == "main", z == input$plot_z) %>%
        ggplot() +
        {
          if (input$sd_as_alpha) {
            geom_raster(aes(x, y, fill = mean, alpha = sd))
          } else {
            geom_raster(aes(x, y, fill = mean))
          }
        } +
        scale_fill_viridis_c() +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "C2" && !input$comic) {

      interpol_grid() %>% dplyr::filter(dependent_var_id == "C2", pred_grid_id == "main", z == input$plot_z) %>%
        ggplot() +
        {
          if (input$sd_as_alpha) {
            geom_raster(aes(x, y, fill = mean, alpha = sd))
          } else {
            geom_raster(aes(x, y, fill = mean))
          }
        } +
        scale_fill_viridis_c(option = "plasma") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "mobility_clemens" && input$comic) {

      mobility_clemens() %>%
        ggplot() +
        geom_raster(aes(x, y, fill = speed_km_per_decade)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "cividis") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "mobility_clemens" && !input$comic) {

      mobility_clemens() %>% dplyr::filter(z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = speed_km_per_decade)) +#, alpha = sd)) +
        scale_fill_viridis_c(option = "cividis") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "mobility_stephan" && input$comic) {

      mobility_stephan() %>%
        ggplot() +
        geom_raster(aes(x, y, fill = J_final_outlier_removed)) +#, alpha = sd)) +
        facet_wrap(~z) +
        scale_fill_viridis_c(option = "cividis") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "mobility_stephan" && !input$comic) {

      mobility_stephan() %>% dplyr::filter(z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = J_final_outlier_removed)) +#, alpha = sd)) +
        scale_fill_viridis_c(option = "cividis") +
        scale_alpha_continuous(range = c(1, 0), na.value = 0)

    } else if (input$plot_type == "clemens_origin_segments" && input$comic) {

      mobility_clemens() %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean_C1)) +
        geom_segment(aes(x, y, xend = x_origin, yend = y_origin), arrow = arrow(length = unit(0.1,"cm"))) +
        facet_wrap(~z)

    } else if (input$plot_type == "clemens_origin_segments" && !input$comic) {

      mobility_clemens() %>% dplyr::filter(z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean_C1)) +
        geom_segment(aes(x, y, xend = x_origin, yend = y_origin), arrow = arrow(length = unit(0.1,"cm")))

    } else if (input$plot_type == "clemens_regional_curves") {

      mobility_clemens() %>%
        # main
        dplyr::group_by(region_id, z, independent_table_id, kernel_setting_id) %>%
        dplyr::summarise(
          mean_speed_km_per_decade = mean(speed_km_per_decade)#,
          #mean_angle_deg = mobest::mean_deg(angle_deg)
        ) %>%
        ggplot() +
        geom_line(
          aes(
            x = z, y = mean_speed_km_per_decade,
            group = interaction(independent_table_id, kernel_setting_id)#,
            #color = mean_angle_deg
          ),
          alpha = 0.5
        ) +
        facet_wrap(dplyr::vars(region_id))
    }

  })

}
