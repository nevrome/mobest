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

  # GPR
  interpol_grid <- reactive({

    #req(pred_grids())

    cache_file_path <- paste0(
      "cache/interpol_",
      "cs", input$pred_grid_spatial_cell_size, "_",
      "td", input$pred_grid_temporal_distance, "_",
      "ds", input$kernel_spatial_size, "_",
      "dt", input$kernel_temporal_size, "_",
      "g", input$kernel_nugget, "_",
      "smdx", input$stephan_mobility_delta_x, "_",
      "smdy", input$stephan_mobility_delta_y, "_",
      "smdz", input$stephan_mobility_delta_z,
      ".RData"
    )

    if (file.exists(cache_file_path)) {
      load(cache_file_path)
      interpol_grid
    } else {

      main_pred_grid <- mobest::create_prediction_grid(
        area,
        mobility_regions,
        spatial_cell_size = input$pred_grid_spatial_cell_size,
        time_layers = seq(-7500, 1500, input$pred_grid_temporal_distance)
      )

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
          main = main_pred_grid,
          offset_x = main_pred_grid %>% dplyr::mutate(x = x + input$stephan_mobility_delta_x),
          offset_y = main_pred_grid %>% dplyr::mutate(y = y + input$stephan_mobility_delta_y),
          offset_z = main_pred_grid %>% dplyr::mutate(z = z + input$stephan_mobility_delta_z)
        )
      )

      withProgress(message = "GPR", {
        interpol_grid <- mobest::run_model_grid(model_grid, quiet = F)
      })

      save(interpol_grid, file = cache_file_path)

      interpol_grid

    }

  })

  # mobility calculation
  mobility_clemens <- reactive({
    interpol_grid <- interpol_grid() %>% dplyr::filter(pred_grid_id == "main")
    # updateSelectInput(session, "mobility_algorithm", selected = "clemens")
    withProgress(message = "Mobility estimation (Clemens)", {
      origin_grid <- mobest::search_spatial_origin(interpol_grid, steps = input$clemens_mobility_steps)
      mobility <- origin_grid#mobest::estimate_mobility(origin_grid)
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

    req(input$plot_z)

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
        geom_segment(aes(x, y, xend = x_origin, yend = y_origin), arrow = arrow(length = unit(0.3,"cm"))) +
        facet_wrap(~z)

    } else if (input$plot_type == "clemens_origin_segments" && !input$comic) {

      mobility_clemens() %>% dplyr::filter(z == input$plot_z) %>%
        ggplot() +
        geom_raster(aes(x, y, fill = mean_C1)) +
        geom_segment(aes(x, y, xend = x_origin, yend = y_origin), arrow = arrow(length = unit(0.3,"cm")))

    } else if (input$plot_type == "clemens_directed_mobility_regional_curves") {

      mobility_clemens() %>%
        # main
        dplyr::group_by(region_id, z, independent_table_id, kernel_setting_id) %>%
        dplyr::summarise(
          mean_directed_distance = sqrt(mean(x - x_origin)^2 + mean(y - y_origin)^2),
          mean_angle_deg = mobest::vec2deg(c(mean(x_origin - x), mean(y_origin - y)))
        ) %>%
        ggplot() +
        geom_line(
          aes(
            x = z, y = mean_directed_distance,
            group = interaction(independent_table_id, kernel_setting_id),
            color = mean_angle_deg
          )
        ) +
        facet_wrap(dplyr::vars(region_id)) +
        scale_color_gradientn(
          colours = c("#F5793A", "#85C0F9", "#85C0F9", "#A95AA1", "#A95AA1", "#33a02c", "#33a02c", "#F5793A")
        )
    } else if (input$plot_type == "clemens_absolute_mobility_regional_curves") {

      mobility_clemens() %>%
        # main
        dplyr::group_by(region_id, z, independent_table_id, kernel_setting_id) %>%
        dplyr::summarise(
          mean_absolute_distance = mean(sqrt((x - x_origin)^2 + (y - y_origin)^2))
        ) %>%
        ggplot() +
        geom_line(
          aes(
            x = z, y = mean_absolute_distance,
            group = interaction(independent_table_id, kernel_setting_id)#,
            #color = mean_angle_deg
          ),
          color = "red"
        ) +
        facet_wrap(dplyr::vars(region_id))
    }

  })

}
