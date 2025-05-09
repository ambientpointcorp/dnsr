#' @title Dynamic Nelson-Siegel Model Visualization Functions
#' @description Functions for visualizing DNS model results and yield curves.

#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_bar facet_wrap labs theme_minimal
#' @importFrom ggplot2 theme scale_color_manual scale_linetype_manual geom_tile
#' @importFrom ggplot2 scale_fill_gradientn scale_x_continuous scale_y_continuous scale_y_date geom_text
#' @importFrom ggplot2 annotate geom_ribbon sec_axis element_text element_blank
#' @importFrom reshape2 melt
#' @importFrom dplyr filter select mutate arrange
#' @importFrom graphics par boxplot matplot
#' @importFrom grDevices rainbow
#' @importFrom stats quantile smooth.spline predict

#' @title Plot Yield Curves and Fitted DNS Models
#'
#' @description
#' Creates a plot comparing actual yield curves with DNS model fits.
#'
#' @param yield_data Data frame. The original yield curve data.
#' @param dns_result List. The DNS fitting result from fit_dns_model().
#' @param dates_to_plot Date vector. Specific dates to plot. If NULL, uses the first few dates.
#' @param max_plots Integer. Maximum number of plots to create if dates_to_plot is NULL.
#' @param show_rmse Logical. If TRUE, display RMSE values on the plots.
#'
#' @return A ggplot object with the yield curve plots.
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' dns_result <- fit_dns_model(yield_data, n_days = 30)
#' p <- plot_yield_curves(yield_data, dns_result)
#' print(p)
#' }
#'
#' @export
plot_yield_curves <- function(yield_data, dns_result, dates_to_plot = NULL,
                              max_plots = 5, show_rmse = TRUE) {
  # Parameter validation
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_list(dns_result)
  checkmate::assert_names(
    names(dns_result),
    must.include = c("parameters", "fitted_yields", "maturities")
  )
  if (!is.null(dates_to_plot)) {
    checkmate::assert_date(dates_to_plot)
  }
  checkmate::assert_count(max_plots, positive = TRUE)
  checkmate::assert_logical(show_rmse, len = 1)

  # If no dates specified, use the first few dates
  if (is.null(dates_to_plot)) {
    all_dates <- intersect(yield_data$Date, dns_result$parameters$Date)
    dates_to_plot <- all_dates[1:min(max_plots, length(all_dates))]
  }

  # Handle case with no dates to plot
  if (length(dates_to_plot) == 0) {
    warning("No matching dates found between yield data and DNS results")
    return(NULL)
  }

  # Prepare data for plotting
  plot_data <- data.frame()

  # Get maturities
  maturities <- as.numeric(colnames(yield_data)[-1])

  # For each date, add actual and fitted yields
  for (date in dates_to_plot) {
    # Find row index for this date
    date_idx <- which(yield_data$Date == date)

    if (length(date_idx) > 0) {
      # Actual yields
      actual_yields <- as.numeric(yield_data[date_idx, -1])

      for (i in 1:length(maturities)) {
        if (!is.na(actual_yields[i])) {
          plot_data <- rbind(
            plot_data,
            data.frame(
              date = as.character(date),
              maturity = maturities[i],
              yield = actual_yields[i] * 100,  # Convert to percentage
              type = "Actual"
            )
          )
        }
      }

      # Fitted yields (if available)
      dns_date_idx <- which(dns_result$parameters$Date == date)
      if (length(dns_date_idx) > 0) {
        # Get RMSE for this date
        rmse_bp <- dns_result$parameters$rmse[dns_date_idx] * 10000  # Convert to basis points

        # Get fitted yields
        fitted_yields <- dns_result$fitted_yields[dns_date_idx, ]

        for (i in 1:length(maturities)) {
          if (!is.na(fitted_yields[i])) {
            plot_data <- rbind(
              plot_data,
              data.frame(
                date = as.character(date),
                maturity = maturities[i],
                yield = fitted_yields[i] * 100,  # Convert to percentage
                type = "Fitted"
              )
            )
          }
        }

        # Add RMSE annotation if requested
        if (show_rmse) {
          plot_data <- rbind(
            plot_data,
            data.frame(
              date = as.character(date),
              maturity = max(maturities, na.rm = TRUE),
              yield = max(c(actual_yields, fitted_yields), na.rm = TRUE) * 100,
              type = sprintf("RMSE: %.1f bp", rmse_bp)
            )
          )
        }
      }
    }
  }

  # Check if we have data to plot
  if (nrow(plot_data) == 0) {
    warning("No valid data to plot yield curves")
    return(NULL)
  }

  # Create plot
  if (show_rmse) {
    # Create helper column for plotting
    plot_data$plot_group <- plot_data$type
    plot_data$plot_group[grepl("RMSE", plot_data$type)] <- "Annotation"

    # Base plot
    p <- ggplot(
      plot_data,
      aes(
        x = maturity,
        y = yield,
        color = type,
        linetype = plot_group,
        group = interaction(date, type)
      )
    ) +
      geom_line() +
      geom_point(data = plot_data[plot_data$type %in% c("Actual", "Fitted"), ]) +
      # Add RMSE as text annotation
      geom_text(
        data = plot_data[grepl("RMSE", plot_data$type), ],
        aes(label = type),
        hjust = 1,
        vjust = 1,
        size = 3,
        show.legend = FALSE
      ) +
      facet_wrap(~ date) +
      labs(
        title = "Treasury Yield Curves and DNS Model Fits",
        x = "Maturity (years)",
        y = "Yield (%)"
      ) +
      scale_color_manual(
        values = c(
          "Actual" = "blue",
          "Fitted" = "red",
          "RMSE: 0.1 bp" = "transparent"
        ),
        name = "Type"
      ) +
      scale_linetype_manual(
        values = c(
          "Actual" = "solid",
          "Fitted" = "dashed",
          "Annotation" = "blank"
        ),
        guide = "none"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  } else {
    # Simpler plot without RMSE
    p <- ggplot(
      plot_data,
      aes(
        x = maturity,
        y = yield,
        color = type,
        linetype = type,
        group = interaction(date, type)
      )
    ) +
      geom_line() +
      geom_point() +
      facet_wrap(~ date) +
      labs(
        title = "Treasury Yield Curves and DNS Model Fits",
        x = "Maturity (years)",
        y = "Yield (%)",
        color = "Type",
        linetype = "Type"
      ) +
      scale_color_manual(values = c("Actual" = "blue", "Fitted" = "red")) +
      scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "dashed")) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  return(p)
}

#' @title Plot DNS Factor Time Series
#'
#' @description
#' Creates plots of DNS factor time series and fit quality.
#'
#' @param dns_result List. The DNS fitting result from fit_dns_model().
#' @param show_changes Logical. If TRUE, include plots of factor changes.
#'
#' @return A list of ggplot objects:
#'   - factors: Plot of factor time series
#'   - rmse: Plot of RMSE time series
#'   - changes: Plot of factor changes (if show_changes = TRUE)
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' dns_result <- fit_dns_model(yield_data, n_days = 100)
#' plots <- plot_factor_time_series(dns_result)
#' print(plots$factors)
#' print(plots$rmse)
#' }
#'
#' @export
plot_factor_time_series <- function(dns_result, show_changes = FALSE) {
  # Parameter validation
  checkmate::assert_list(dns_result)
  checkmate::assert_names(names(dns_result), must.include = "parameters")
  checkmate::assert_logical(show_changes, len = 1)

  # Extract parameters data frame
  params_df <- dns_result$parameters

  # Create empty list to store plots
  plots <- list()

  # Convert to long format for plotting factors
  params_long <- reshape2::melt(
    params_df,
    id.vars = "Date",
    measure.vars = c("beta1", "beta2", "beta3"),
    variable.name = "Factor",
    value.name = "Value"
  )

  # Create factor labels
  factor_labels <- c(
    beta1 = "Level (β₁)",
    beta2 = "Slope (β₂)",
    beta3 = "Curvature (β₃)"
  )

  # Plot factors
  plots$factors <- ggplot(params_long, aes(x = Date, y = Value, color = Factor)) +
    geom_line() +
    labs(
      title = "DNS Model Factor Time Series",
      x = "Date",
      y = "Factor Value",
      color = "Factor"
    ) +
    scale_color_manual(
      values = c("blue", "red", "green"),
      labels = factor_labels
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Plot RMSE
  plots$rmse <- ggplot(params_df, aes(x = Date, y = rmse * 10000)) +  # Convert to basis points
    geom_line(color = "purple") +
    labs(
      title = "DNS Model Fit Error (RMSE)",
      x = "Date",
      y = "RMSE (basis points)"
    ) +
    theme_minimal()

  # Plot factor changes if requested
  if (show_changes && all(c("beta1_change", "beta2_change", "beta3_change") %in% names(params_df))) {
    # Convert to long format for plotting changes
    changes_long <- reshape2::melt(
      params_df,
      id.vars = "Date",
      measure.vars = c("beta1_change", "beta2_change", "beta3_change"),
      variable.name = "Factor",
      value.name = "Change"
    )

    # Create change labels
    change_labels <- c(
      beta1_change = "Level (β₁) Change",
      beta2_change = "Slope (β₂) Change",
      beta3_change = "Curvature (β₃) Change"
    )

    # Plot changes
    plots$changes <- ggplot(changes_long, aes(x = Date, y = Change, color = Factor)) +
      geom_line() +
      labs(
        title = "DNS Model Factor Changes",
        x = "Date",
        y = "Daily Factor Change",
        color = "Factor"
      ) +
      scale_color_manual(
        values = c("blue", "red", "green"),
        labels = change_labels
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  return(plots)
}

#' @title Plot Simulated DNS Factors
#'
#' @description
#' Creates plots of simulated DNS factors from a DNS simulation.
#'
#' @param sim_result List. The simulation result from simulate_dns_factors().
#' @param path_subset Integer vector. Indices of paths to plot. If NULL, plots summary statistics.
#' @param include_quantiles Logical. If TRUE, include quantile paths.
#' @param plot_trajectories Logical. If TRUE, plot individual trajectories instead of summary statistics.
#'
#' @return A list of ggplot objects, one for each factor.
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' sim_result <- simulate_dns_factors(params, initial_factors)
#' plots <- plot_simulated_factors(sim_result)
#' print(plots$level)
#' }
#'
#' @export
plot_simulated_factors <- function(sim_result, path_subset = NULL,
                                   include_quantiles = TRUE,
                                   plot_trajectories = FALSE) {
  # Parameter validation
  checkmate::assert_list(sim_result)
  checkmate::assert_names(
    names(sim_result),
    must.include = c("factors", "time", "params")
  )

  # Extract relevant information
  factors <- sim_result$factors
  time <- sim_result$time
  n_times <- dim(factors)[1]
  n_paths <- dim(factors)[2]

  # Factor names
  factor_names <- if ("names" %in% names(sim_result$params)) {
    sim_result$params$names
  } else {
    c("Level", "Slope", "Curvature")
  }

  # Handle path subset
  if (is.null(path_subset) && plot_trajectories) {
    # If no subset specified but trajectories requested, use first 10 paths
    path_subset <- 1:min(10, n_paths)
  }

  # Create a list to store plots
  plots <- list()
  factor_colors <- c("blue", "red", "green")

  # Generate plots for each factor
  for (f in 1:3) {
    factor_name <- factor_names[f]

    if (plot_trajectories) {
      # Plot individual trajectories
      plot_data <- data.frame()

      for (p in path_subset) {
        path_data <- data.frame(
          Time = time,
          Value = factors[, p, f],
          Path = paste("Path", p)
        )
        plot_data <- rbind(plot_data, path_data)
      }

      plots[[tolower(factor_name)]] <- ggplot(
        plot_data,
        aes(x = Time, y = Value, color = Path)
      ) +
        geom_line(alpha = 0.7) +
        labs(
          title = paste("Simulated", factor_name, "Factor Paths"),
          x = "Time (years)",
          y = factor_name
        ) +
        theme_minimal()

    } else {
      # Calculate summary statistics for each time point
      means <- apply(factors[, , f], 1, mean)
      medians <- apply(factors[, , f], 1, median)
      q10 <- apply(factors[, , f], 1, function(x) quantile(x, 0.1))
      q90 <- apply(factors[, , f], 1, function(x) quantile(x, 0.9))

      # Create data for the mean path
      mean_data <- data.frame(
        Time = time,
        Value = means,
        Type = "Mean"
      )

      # Create plotting data
      plot_data <- mean_data

      if (include_quantiles) {
        # Add quantile data
        quantile_data <- data.frame(
          Time = rep(time, 3),
          Value = c(medians, q10, q90),
          Type = rep(c("Median", "10% Quantile", "90% Quantile"), each = length(time))
        )
        plot_data <- rbind(plot_data, quantile_data)
      }

      # Create the plot
      plots[[tolower(factor_name)]] <- ggplot(
        plot_data,
        aes(x = Time, y = Value, color = Type, linetype = Type)
      ) +
        geom_line() +
        labs(
          title = paste("Simulated", factor_name, "Factor"),
          x = "Time (years)",
          y = factor_name,
          color = "Statistic",
          linetype = "Statistic"
        ) +
        theme_minimal()
    }
  }

  return(plots)
}

#' @title Plot Simulated Yield Curves
#'
#' @description
#' Creates plots of simulated yield curves from a DNS simulation.
#'
#' @param sim_result List. The simulation result from simulate_dns_yield_curves().
#' @param time_points Numeric vector. Time points to plot (in years).
#' @param path_subset Integer vector. Indices of paths to plot. If NULL, plots summary statistics.
#' @param show_forwards Logical. If TRUE, also plot forward curves.
#'
#' @return A list of ggplot objects:
#'   - yield_curves: Plot of yield curves at specified time points
#'   - term_structure: Plot of the term structure evolution (if time_points contains a single value)
#'   - forwards: Plot of forward curves (if show_forwards = TRUE)
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' sim_result <- simulate_dns_yield_curves(params, initial_factors, maturities)
#' plots <- plot_simulated_yields(sim_result, time_points = c(0, 1, 5, 10))
#' print(plots$yield_curves)
#' }
#'
#' @export
plot_simulated_yields <- function(sim_result, time_points = c(0, 1, 5, 10),
                                  path_subset = NULL, show_forwards = FALSE) {
  # Parameter validation
  checkmate::assert_list(sim_result)
  checkmate::assert_names(
    names(sim_result),
    must.include = c("yields", "time", "maturities")
  )
  checkmate::assert_numeric(time_points)

  # Extract relevant information
  yields <- sim_result$yields
  time <- sim_result$time
  maturities <- sim_result$maturities

  # Find indices of requested time points
  time_indices <- numeric(length(time_points))
  for (i in seq_along(time_points)) {
    # Find closest time point
    time_indices[i] <- which.min(abs(time - time_points[i]))
  }

  # Create a list to store plots
  plots <- list()

  # Get dimensions
  n_times <- dim(yields)[1]
  n_maturities <- dim(yields)[2]
  n_paths <- dim(yields)[3]

  # Handle path subset
  if (is.null(path_subset)) {
    # If no subset specified, use first 5 paths
    path_subset <- 1:min(5, n_paths)
  }

  # Plot yield curves at specified time points
  plot_data <- data.frame()

  for (t_idx in seq_along(time_indices)) {
    t <- time_indices[t_idx]
    t_value <- time[t]

    for (p in path_subset) {
      for (m in seq_along(maturities)) {
        plot_data <- rbind(
          plot_data,
          data.frame(
            Time = paste("Year", round(t_value, 1)),
            Maturity = maturities[m],
            Yield = yields[t, m, p] * 100,  # Convert to percentage
            Path = paste("Path", p)
          )
        )
      }
    }
  }

  # Create yield curve plot
  plots$yield_curves <- ggplot(
    plot_data,
    aes(
      x = Maturity,
      y = Yield,
      color = Path,
      group = interaction(Time, Path)
    )
  ) +
    geom_line() +
    facet_wrap(~ Time) +
    labs(
      title = "Simulated Yield Curves",
      x = "Maturity (years)",
      y = "Yield (%)",
      color = "Path"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # If only one time point is requested, also create a term structure evolution plot
  if (length(time_points) == 1) {
    # Calculate statistics across paths
    t_idx <- time_indices[1]

    # Create data for term structure evolution
    term_data <- data.frame()

    for (m in seq_along(maturities)) {
      yields_at_m <- yields[, m, ]

      means <- apply(yields_at_m, 1, mean)
      q10 <- apply(yields_at_m, 1, function(x) quantile(x, 0.1))
      q90 <- apply(yields_at_m, 1, function(x) quantile(x, 0.9))

      for (t in 1:n_times) {
        term_data <- rbind(
          term_data,
          data.frame(
            Time = time[t],
            Maturity = paste(maturities[m], "yr"),
            Mean = means[t] * 100,
            Q10 = q10[t] * 100,
            Q90 = q90[t] * 100
          )
        )
      }
    }

    # Create term structure evolution plot
    plots$term_structure <- ggplot(
      term_data,
      aes(
        x = Time,
        y = Mean,
        group = Maturity,
        color = Maturity
      )
    ) +
      geom_line() +
      geom_ribbon(
        aes(ymin = Q10, ymax = Q90, fill = Maturity),
        alpha = 0.2,
        color = NA
      ) +
      labs(
        title = "Term Structure Evolution",
        x = "Time (years)",
        y = "Yield (%)",
        color = "Maturity",
        fill = "Maturity"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  # Plot forward curves if requested
  if (show_forwards) {
    # Calculate forward rates
    forwards <- calculate_dns_forwards(
      factors = sim_result$factors,
      maturities = maturities,
      lambda = if ("lambda" %in% names(sim_result)) sim_result$lambda else 0.0609
    )

    # Prepare data
    forward_data <- data.frame()

    for (t_idx in seq_along(time_indices)) {
      t <- time_indices[t_idx]
      t_value <- time[t]

      for (p in path_subset) {
        for (m in seq_along(maturities)) {
          forward_data <- rbind(
            forward_data,
            data.frame(
              Time = paste("Year", round(t_value, 1)),
              Maturity = maturities[m],
              Forward = forwards[t, m, p] * 100,  # Convert to percentage
              Path = paste("Path", p)
            )
          )
        }
      }
    }

    # Create forward curve plot
    plots$forwards <- ggplot(
      forward_data,
      aes(
        x = Maturity,
        y = Forward,
        color = Path,
        group = interaction(Time, Path)
      )
    ) +
      geom_line() +
      facet_wrap(~ Time) +
      labs(
        title = "Simulated Forward Rate Curves",
        x = "Maturity (years)",
        y = "Forward Rate (%)",
        color = "Path"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  return(plots)
}

#' @title Create Interactive Yield Curve Animation
#'
#' @description
#' Creates an interactive animation of yield curve evolution over time.
#'
#' @param sim_result List. The simulation result from simulate_dns_yield_curves().
#' @param path_idx Integer. Index of the path to animate.
#' @param time_subset Integer vector. Indices of time points to include.
#' @param title Character. Title for the animation.
#'
#' @return An animation object that can be displayed in the R viewer.
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' sim_result <- simulate_dns_yield_curves(params, initial_factors, maturities)
#' anim <- animate_yield_curves(sim_result)
#' print(anim)
#' }
#'
#' @importFrom gganimate transition_states shadow_mark ease_aes animate view_static
#' @importFrom gifski save_gif
#'
#' @export
animate_yield_curves <- function(sim_result, path_idx = 1, time_subset = NULL,
                                 title = "Yield Curve Evolution") {
  # Check if required packages are installed
  if (!requireNamespace("gganimate", quietly = TRUE) ||
      !requireNamespace("gifski", quietly = TRUE)) {
    stop(
      "Packages 'gganimate' and 'gifski' are required for animation. ",
      "Please install them with install.packages(c('gganimate', 'gifski'))"
    )
  }

  # Extract relevant information
  yields <- sim_result$yields
  time <- sim_result$time
  maturities <- sim_result$maturities

  # Check path index
  n_paths <- dim(yields)[3]
  if (path_idx < 1 || path_idx > n_paths) {
    stop("path_idx must be between 1 and ", n_paths)
  }

  # Handle time subset
  n_times <- dim(yields)[1]
  if (is.null(time_subset)) {
    # If no subset specified, use all time points or a subset if there are many
    if (n_times > 50) {
      time_subset <- seq(1, n_times, length.out = 50)
    } else {
      time_subset <- 1:n_times
    }
  }

  # Prepare data for animation
  anim_data <- data.frame()

  for (t_idx in time_subset) {
    for (m in seq_along(maturities)) {
      anim_data <- rbind(
        anim_data,
        data.frame(
          Time = time[t_idx],
          TimeIndex = t_idx,
          Maturity = maturities[m],
          Yield = yields[t_idx, m, path_idx] * 100  # Convert to percentage
        )
      )
    }
  }

  # Create base plot
  p <- ggplot(anim_data, aes(x = Maturity, y = Yield, group = TimeIndex)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    labs(
      title = title,
      subtitle = "Time: {closest_state} years",
      x = "Maturity (years)",
      y = "Yield (%)"
    ) +
    theme_minimal() +
    ylim(min(anim_data$Yield) * 0.9, max(anim_data$Yield) * 1.1)

  # Create animation
  anim <- gganimate::transition_states(
    states = Time,
    transition_length = 2,
    state_length = 1
  ) %>%
    gganimate::shadow_mark(alpha = 0.3, size = 0.5) %>%
    gganimate::ease_aes('cubic-in-out')

  # Apply animation
  anim_p <- p + anim

  # Return animation object
  return(anim_p)
}

#' @title Plot DNS Factor Loadings
#'
#' @description
#' Creates a plot showing the factor loadings of the DNS model.
#'
#' @param maturities Numeric vector. Maturities (in years) for which to calculate loadings.
#' @param lambda Numeric. DNS decay parameter lambda.
#' @param show_labels Logical. If TRUE, add labels to the plot.
#'
#' @return A ggplot object with the factor loadings plot.
#'
#' @examples
#' \dontrun{
#' maturities <- seq(0, 30, length.out = 100)
#' p <- plot_dns_loadings(maturities)
#' print(p)
#' }
#'
#' @export
plot_dns_loadings <- function(maturities = seq(0, 30, length.out = 100),
                              lambda = 0.0609, show_labels = TRUE) {
  # Parameter validation
  checkmate::assert_numeric(maturities, lower = 0)
  checkmate::assert_numeric(lambda, lower = 0, len = 1)
  checkmate::assert_logical(show_labels, len = 1)

  # Calculate factor loadings
  loading_data <- data.frame(
    Maturity = maturities,
    Level = rep(1, length(maturities)),
    Slope = numeric(length(maturities)),
    Curvature = numeric(length(maturities))
  )

  for (i in seq_along(maturities)) {
    m <- maturities[i]

    # Handle numerical stability for very small maturities
    if (m * lambda < 1e-10) {
      loading_2 <- 1.0
      loading_3 <- 0.0
    } else {
      loading_2 <- (1 - exp(-lambda * m)) / (lambda * m)
      loading_3 <- loading_2 - exp(-lambda * m)
    }

    loading_data$Slope[i] <- loading_2
    loading_data$Curvature[i] <- loading_3
  }

  # Convert to long format for plotting
  loading_long <- reshape2::melt(
    loading_data,
    id.vars = "Maturity",
    variable.name = "Factor",
    value.name = "Loading"
  )

  # Create plot
  p <- ggplot(loading_long, aes(x = Maturity, y = Loading, color = Factor)) +
    geom_line(size = 1) +
    labs(
      title = "DNS Model Factor Loadings",
      subtitle = paste("λ =", lambda),
      x = "Maturity (years)",
      y = "Loading",
      color = "Factor"
    ) +
    scale_color_manual(
      values = c("Level" = "blue", "Slope" = "red", "Curvature" = "green")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Add labels if requested
  if (show_labels) {
    # Find points for labels
    label_points <- data.frame(
      Maturity = c(
        max(maturities) * 0.9,
        0.1,
        maturities[which.max(loading_data$Curvature)]
      ),
      Loading = c(
        1,
        loading_data$Slope[which(maturities >= 0.1)[1]],
        max(loading_data$Curvature)
      ),
      Factor = c("Level", "Slope", "Curvature"),
      Label = c("Level (β₁)", "Slope (β₂)", "Curvature (β₃)")
    )

    p <- p + geom_text(
      data = label_points,
      aes(label = Label, color = Factor),
      hjust = -0.1,
      vjust = 1.2,
      size = 4
    )
  }

  return(p)
}

#' @title Create Heat Map of Yield Curve Evolution
#'
#' @description
#' Creates a heat map visualization of yield curve evolution over time.
#'
#' @param yield_data Data frame. Historical yield curve data with 'Date' column.
#' @param reverse_time Logical. If TRUE, display most recent dates at the top.
#' @param plot_type Character. One of "yields", "changes", or "curvature".
#' @param smooth_data Logical. If TRUE, apply smoothing to the heat map data.
#' @param color_palette Character vector. Colors for the heat map.
#'
#' @return A ggplot object with the heat map.
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' p <- plot_yield_heatmap(yield_data)
#' print(p)
#' }
#'
#' @importFrom ggplot2 geom_tile scale_fill_gradientn scale_x_continuous scale_y_date
#'
#' @export
plot_yield_heatmap <- function(yield_data, reverse_time = FALSE,
                               plot_type = "yields", smooth_data = FALSE,
                               color_palette = NULL) {
  # Parameter validation
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_logical(reverse_time, len = 1)
  checkmate::assert_choice(plot_type, c("yields", "changes", "curvature"))
  checkmate::assert_logical(smooth_data, len = 1)

  # Get maturities from column names
  maturities <- as.numeric(colnames(yield_data)[-1])

  # Prepare data for heatmap
  if (plot_type == "changes") {
    # Calculate daily changes in yields
    yield_changes <- yield_data
    for (col in names(yield_data)[-1]) {
      yield_changes[[col]] <- c(NA, diff(yield_data[[col]]))
    }
    plot_df <- yield_changes
    title_prefix <- "Daily Yield Changes"
    legend_label <- "Change (%)"
  } else if (plot_type == "curvature") {
    # Calculate curvature (2nd derivative approximation)
    yield_curvature <- yield_data
    for (i in 2:(length(maturities) - 1)) {
      col <- names(yield_data)[i + 1]
      col_prev <- names(yield_data)[i]
      col_next <- names(yield_data)[i + 2]
      m <- maturities[i]
      m_prev <- maturities[i - 1]
      m_next <- maturities[i + 1]

      # Centered finite difference approximation
      # f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h^2
      # Adjust for non-uniform grid
      h1 <- m - m_prev
      h2 <- m_next - m
      h_avg <- (h1 + h2) / 2

      yield_curvature[[col]] <- (
        yield_data[[col_next]] - 2*yield_data[[col]] + yield_data[[col_prev]]
      ) / (h_avg^2)
    }

    # Set edge columns to NA for curvature
    yield_curvature[[names(yield_data)[2]]] <- NA
    yield_curvature[[names(yield_data)[ncol(yield_data)]]] <- NA

    plot_df <- yield_curvature
    title_prefix <- "Yield Curve Curvature"
    legend_label <- "Curvature"
  } else {
    # Just use yields
    plot_df <- yield_data
    title_prefix <- "Yield Curve Evolution"
    legend_label <- "Yield (%)"
  }

  # Make sure all values are in percentage for display
  for (col in names(plot_df)[-1]) {
    if (max(plot_df[[col]], na.rm = TRUE) < 1 && plot_type == "yields") {
      plot_df[[col]] <- plot_df[[col]] * 100
    }
  }

  # Convert to long format for ggplot
  plot_data <- reshape2::melt(
    plot_df,
    id.vars = "Date",
    variable.name = "Maturity",
    value.name = "Value"
  )

  # Convert maturity to numeric
  plot_data$Maturity <- as.numeric(as.character(plot_data$Maturity))

  # Apply smoothing if requested
  if (smooth_data && requireNamespace("mgcv", quietly = TRUE)) {
    # Create a grid for prediction
    date_seq <- seq(min(plot_data$Date), max(plot_data$Date), length.out = 100)
    maturity_seq <- seq(min(plot_data$Maturity), max(plot_data$Maturity), length.out = 50)
    grid <- expand.grid(Date = date_seq, Maturity = maturity_seq)

    # Fit a smooth surface
    smooth_model <- mgcv::gam(Value ~ s(as.numeric(Date), Maturity), data = plot_data)

    # Predict on the grid
    grid$Value <- predict(smooth_model, newdata = grid)

    # Use the smoothed data
    plot_data <- grid
  }

  # Set default color palette if not provided
  if (is.null(color_palette)) {
    if (plot_type == "changes") {
      color_palette <- c("blue", "white", "red")
    } else if (plot_type == "curvature") {
      color_palette <- c("purple", "white", "orange")
    } else {
      color_palette <- c("navy", "blue", "cyan", "yellow", "orange", "red")
    }
  }

  # Create the heatmap
  p <- ggplot(plot_data, aes(x = Maturity, y = Date, fill = Value)) +
    geom_tile() +
    labs(
      title = paste(title_prefix, "Heatmap"),
      x = "Maturity (years)",
      y = "Date",
      fill = legend_label
    ) +
    scale_fill_gradientn(colors = color_palette) +
    scale_x_continuous(breaks = maturities) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  # Reverse time axis if requested
  if (reverse_time) {
    p <- p + scale_y_date()
  }

  return(p)
}

#' @title Plot 3D Yield Surface
#'
#' @description
#' Creates a 3D visualization of the yield curve surface.
#'
#' @param yield_data Data frame. Historical yield curve data with 'Date' column.
#' @param interactive Logical. If TRUE, create an interactive plot with plotly.
#' @param sample_dates Integer. Number of dates to sample for plotting.
#' @param color_scale Character. Color scale for the surface.
#'
#' @return A 3D plot object.
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' plot_yield_surface(yield_data)
#' }
#'
#' @export
plot_yield_surface <- function(yield_data, interactive = TRUE,
                               sample_dates = 50, color_scale = "Viridis") {
  # Check if required packages are installed
  if (interactive && !requireNamespace("plotly", quietly = TRUE)) {
    warning("Package 'plotly' is required for interactive plots. Using static plot instead.")
    interactive <- FALSE
  }

  if (!interactive && !requireNamespace("lattice", quietly = TRUE)) {
    stop("Package 'lattice' is required for static 3D plots. Please install it.")
  }

  # Parameter validation
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_logical(interactive, len = 1)
  checkmate::assert_count(sample_dates, positive = TRUE)

  # Get maturities from column names
  maturities <- as.numeric(colnames(yield_data)[-1])

  # Sample dates if needed
  if (nrow(yield_data) > sample_dates) {
    idx <- round(seq(1, nrow(yield_data), length.out = sample_dates))
    sampled_data <- yield_data[idx, ]
  } else {
    sampled_data <- yield_data
  }

  # Convert data to long format
  plot_data <- reshape2::melt(
    sampled_data,
    id.vars = "Date",
    variable.name = "Maturity",
    value.name = "Yield"
  )

  # Convert maturity to numeric
  plot_data$Maturity <- as.numeric(as.character(plot_data$Maturity))

  # Convert to percentage if needed
  if (max(plot_data$Yield, na.rm = TRUE) < 1) {
    plot_data$Yield <- plot_data$Yield * 100
  }

  # Create plot based on type
  if (interactive) {
    # Check if plotly is available again (just to be extra safe)
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive plots. Please install it.")
    }

    # Create interactive 3D surface with plotly
    p <- plotly::plot_ly(
      x = ~Maturity,
      y = ~Date,
      z = ~Yield,
      data = plot_data,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = 3,
        color = ~Yield,
        colorscale = color_scale,
        opacity = 0.8
      )
    ) %>%
      plotly::add_trace(
        type = "mesh3d",
        alphahull = 0
      ) %>%
      plotly::layout(
        title = "3D Yield Curve Surface",
        scene = list(
          xaxis = list(title = "Maturity (years)"),
          yaxis = list(title = "Date"),
          zaxis = list(title = "Yield (%)")
        )
      )
  } else {
    # Create static 3D surface with lattice
    p <- lattice::wireframe(
      Yield ~ Maturity * as.numeric(Date),
      data = plot_data,
      scales = list(arrows = FALSE),
      drape = TRUE,
      colorkey = TRUE,
      screen = list(z = 30, x = -60),
      xlab = "Maturity (years)",
      ylab = "Time",
      zlab = "Yield (%)",
      main = "3D Yield Curve Surface",
      shade = TRUE
    )
  }

  return(p)
}

#' @title Plot Principal Component Analysis of Yield Curve Movements
#'
#' @description
#' Performs Principal Component Analysis (PCA) on yield curve movements and
#' visualizes the results.
#'
#' @param yield_data Data frame. Historical yield curve data with 'Date' column.
#' @param n_components Integer. Number of principal components to extract.
#' @param plot_weights Logical. If TRUE, plot the loadings (eigenvectors).
#' @param plot_scores Logical. If TRUE, plot the scores (projections).
#' @param plot_variance Logical. If TRUE, plot the explained variance.
#'
#' @return A list containing:
#'   - pca: The PCA result object
#'   - plots: A list of created plots
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' pca_result <- plot_yield_pca(yield_data)
#' print(pca_result$plots$weights)
#' }
#'
#' @importFrom stats prcomp
#'
#' @export
plot_yield_pca <- function(yield_data, n_components = 3,
                           plot_weights = TRUE, plot_scores = TRUE,
                           plot_variance = TRUE) {
  # Parameter validation
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_count(n_components, positive = TRUE)
  checkmate::assert_logical(plot_weights, len = 1)
  checkmate::assert_logical(plot_scores, len = 1)
  checkmate::assert_logical(plot_variance, len = 1)

  # Get yields only (exclude Date)
  yields <- as.matrix(yield_data[, -1])

  # Handle missing values
  if (any(is.na(yields))) {
    # Simple imputation - fill with column means
    col_means <- colMeans(yields, na.rm = TRUE)
    for (j in 1:ncol(yields)) {
      na_idx <- is.na(yields[, j])
      if (any(na_idx)) {
        yields[na_idx, j] <- col_means[j]
      }
    }
  }

  # Get maturities from column names
  maturities <- as.numeric(colnames(yield_data)[-1])

  # Calculate yield changes (first differences)
  yield_changes <- apply(yields, 2, diff)

  # Perform PCA on yield changes
  pca_result <- prcomp(yield_changes, scale. = TRUE)

  # Limit to requested number of components
  n_components <- min(n_components, ncol(pca_result$rotation))

  # Create a list to store plots
  plots <- list()

  # Plot component weights (loadings/eigenvectors)
  if (plot_weights) {
    # Prepare data for plotting
    weights_data <- data.frame(
      Maturity = rep(maturities, n_components),
      Weight = as.vector(pca_result$rotation[, 1:n_components]),
      Component = factor(
        rep(1:n_components, each = length(maturities)),
        labels = paste("PC", 1:n_components)
      )
    )

    # Create component labels with variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    component_labels <- sprintf(
      "PC%d (%.1f%%)",
      1:n_components,
      var_explained[1:n_components] * 100
    )

    # Create plot
    plots$weights <- ggplot(
      weights_data,
      aes(x = Maturity, y = Weight, color = Component)
    ) +
      geom_line(size = 1) +
      geom_point() +
      labs(
        title = "PCA Components of Yield Curve Movements",
        x = "Maturity (years)",
        y = "Component Weight",
        color = "Component"
      ) +
      scale_color_manual(
        values = c("red", "blue", "green", "purple", "orange")[1:n_components],
        labels = component_labels
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  # Plot component scores (projections)
  if (plot_scores) {
    # Calculate scores
    scores <- yield_changes %*% pca_result$rotation[, 1:n_components]

    # Prepare data for plotting
    scores_data <- data.frame(
      Date = yield_data$Date[-1],  # Remove first date since we used diff()
      scores
    )
    names(scores_data)[-1] <- paste0("PC", 1:n_components)

    # Convert to long format for plotting
    scores_long <- reshape2::melt(
      scores_data,
      id.vars = "Date",
      variable.name = "Component",
      value.name = "Score"
    )

    # Create component labels with variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    component_labels <- sprintf(
      "PC%d (%.1f%%)",
      1:n_components,
      var_explained[1:n_components] * 100
    )

    # Create plot
    plots$scores <- ggplot(
      scores_long,
      aes(x = Date, y = Score, color = Component)
    ) +
      geom_line() +
      labs(
        title = "Time Series of PCA Component Scores",
        x = "Date",
        y = "Component Score",
        color = "Component"
      ) +
      scale_color_manual(
        values = c("red", "blue", "green", "purple", "orange")[1:n_components],
        labels = component_labels
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  # Plot explained variance
  if (plot_variance) {
    # Calculate explained variance
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    cum_var_explained <- cumsum(var_explained)

    # Prepare data for plotting
    variance_data <- data.frame(
      Component = factor(
        1:length(var_explained),
        labels = paste("PC", 1:length(var_explained))
      ),
      Variance = var_explained * 100,
      Cumulative = cum_var_explained * 100
    )

    # Create plot
    plots$variance <- ggplot(variance_data, aes(x = Component)) +
      geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue") +
      geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
      geom_point(aes(y = Cumulative), color = "red", size = 3) +
      labs(
        title = "Explained Variance by Principal Components",
        x = "Principal Component",
        y = "Variance Explained (%)"
      ) +
      scale_y_continuous(
        sec.axis = sec_axis(~., name = "Cumulative Variance Explained (%)")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(color = "steelblue"),
        axis.title.y.right = element_text(color = "red")
      )
  }

  return(list(
    pca = pca_result,
    plots = plots,
    var_explained = var_explained,
    n_components = n_components
  ))
}

#' @title Compare Nelson-Siegel Model to Other Term Structure Models
#'
#' @description
#' Compares the Nelson-Siegel model fit against other term structure models.
#'
#' @param yield_data Data frame. Historical yield curve data with 'Date' column.
#' @param date_to_compare Date. The specific date to compare models.
#' @param models Character vector. Models to compare, from "NS", "NSS", "Svensson", "Spline".
#' @param spline_df Integer. Degrees of freedom for spline model.
#'
#' @return A list containing:
#'   - plot: A ggplot object with model comparisons
#'   - metrics: A data frame with fit metrics for each model
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2023-01-01")
#' comparison <- compare_term_structure_models(yield_data, Sys.Date() - 30)
#' print(comparison$plot)
#' }
#'
#' @export
compare_term_structure_models <- function(yield_data, date_to_compare = NULL,
                                          models = c("NS", "NSS", "Spline"),
                                          spline_df = 5) {
  # Parameter validation
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_subset(models, c("NS", "NSS", "Svensson", "Spline"))
  checkmate::assert_count(spline_df, positive = TRUE)

  # Select date to compare (default to most recent)
  if (is.null(date_to_compare)) {
    date_to_compare <- max(yield_data$Date)
  } else {
    date_to_compare <- as.Date(date_to_compare)
    if (!(date_to_compare %in% yield_data$Date)) {
      closest_date <- yield_data$Date[which.min(abs(yield_data$Date - date_to_compare))]
      warning(
        sprintf(
          "Date %s not found in data. Using closest date: %s",
          date_to_compare, closest_date
        )
      )
      date_to_compare <- closest_date
    }
  }

  # Extract yield curve for the selected date
  date_idx <- which(yield_data$Date == date_to_compare)
  if (length(date_idx) == 0) {
    stop("Selected date not found in data")
  }

  # Extract yield data for the selected date
  yields <- as.numeric(yield_data[date_idx, -1])
  maturities <- as.numeric(colnames(yield_data)[-1])

  # Create a data frame for storing model fits
  model_fits <- data.frame(
    Maturity = maturities,
    Actual = yields
  )

  # Create a data frame for storing fit metrics
  fit_metrics <- data.frame(
    Model = character(),
    RMSE = numeric(),
    MAE = numeric(),
    MaxError = numeric(),
    stringsAsFactors = FALSE
  )

  # Fit Nelson-Siegel model if requested
  if ("NS" %in% models) {
    ns_result <- fit_ns_model(maturities, yields)
    model_fits$NS <- ns_result$fitted_yields

    # Calculate fit metrics
    rmse <- ns_result$rmse
    mae <- mean(abs(yields - ns_result$fitted_yields), na.rm = TRUE)
    max_error <- max(abs(yields - ns_result$fitted_yields), na.rm = TRUE)

    fit_metrics <- rbind(
      fit_metrics,
      data.frame(
        Model = "Nelson-Siegel",
        RMSE = rmse * 100,  # Convert to basis points
        MAE = mae * 100,
        MaxError = max_error * 100
      )
    )
  }

  # Fit Nelson-Siegel-Svensson model if requested
  if ("NSS" %in% models || "Svensson" %in% models) {
    # Check if SmithWilson package is available
    if (!requireNamespace("SmithWilson", quietly = TRUE)) {
      warning("Package 'SmithWilson' is required for NSS model. Skipping.")
    } else {
      # Fit NSS model
      tryCatch({
        nss_fit <- SmithWilson::Nelson.Siegel.Svensson(
          rate = yields * 100,  # Convert to percentage
          maturity = maturities,
          lambda1 = 0.5,
          lambda2 = 3
        )

        # Extract fitted values
        nss_fitted <- nss_fit$Spot.Rate / 100  # Convert back to decimal
        model_fits$NSS <- nss_fitted

        # Calculate fit metrics
        rmse <- sqrt(mean((yields - nss_fitted)^2, na.rm = TRUE))
        mae <- mean(abs(yields - nss_fitted), na.rm = TRUE)
        max_error <- max(abs(yields - nss_fitted), na.rm = TRUE)

        fit_metrics <- rbind(
          fit_metrics,
          data.frame(
            Model = "Nelson-Siegel-Svensson",
            RMSE = rmse * 100,  # Convert to basis points
            MAE = mae * 100,
            MaxError = max_error * 100
          )
        )
      }, error = function(e) {
        warning("Error fitting NSS model: ", e$message)
      })
    }
  }

  # Fit cubic spline model if requested
  if ("Spline" %in% models) {
    # Fit smoothing spline
    spline_fit <- smooth.spline(maturities, yields, df = spline_df)
    spline_fitted <- predict(spline_fit, maturities)$y
    model_fits$Spline <- spline_fitted

    # Calculate fit metrics
    rmse <- sqrt(mean((yields - spline_fitted)^2, na.rm = TRUE))
    mae <- mean(abs(yields - spline_fitted), na.rm = TRUE)
    max_error <- max(abs(yields - spline_fitted), na.rm = TRUE)

    fit_metrics <- rbind(
      fit_metrics,
      data.frame(
        Model = "Cubic Spline",
        RMSE = rmse * 100,  # Convert to basis points
        MAE = mae * 100,
        MaxError = max_error * 100
      )
    )
  }

  # Convert to long format for plotting
  model_fits_long <- reshape2::melt(
    model_fits,
    id.vars = "Maturity",
    variable.name = "Model",
    value.name = "Yield"
  )

  # Convert to percentage for display
  if (max(model_fits_long$Yield, na.rm = TRUE) < 1) {
    model_fits_long$Yield <- model_fits_long$Yield * 100
  }

  # Create plot
  p <- ggplot(
    model_fits_long,
    aes(x = Maturity, y = Yield, color = Model, linetype = Model)
  ) +
    geom_line(size = 1) +
    geom_point(data = subset(model_fits_long, Model == "Actual")) +
    labs(
      title = paste("Term Structure Model Comparison -", format(date_to_compare, "%Y-%m-%d")),
      x = "Maturity (years)",
      y = "Yield (%)",
      color = "Model",
      linetype = "Model"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Add annotation with fit metrics
  if (nrow(fit_metrics) > 0) {
    # Create annotation text
    annotation_text <- ""
    for (i in 1:nrow(fit_metrics)) {
      annotation_text <- paste0(
        annotation_text,
        fit_metrics$Model[i], ": RMSE = ", round(fit_metrics$RMSE[i], 2), " bp\n"
      )
    }

    # Add annotation to plot
    p <- p + annotate(
      "text",
      x = max(maturities) * 0.8,
      y = min(model_fits_long$Yield, na.rm = TRUE) * 1.1,
      label = annotation_text,
      hjust = 0,
      vjust = 0,
      size = 3
    )
  }

  return(list(
    plot = p,
    metrics = fit_metrics
  ))
}
