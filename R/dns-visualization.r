#' @title Dynamic Nelson-Siegel Model Visualization Functions
#' @description Functions for visualizing DNS model results and yield curves.

#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_bar facet_wrap labs theme_minimal
#' @importFrom ggplot2 theme scale_color_manual scale_linetype_manual geom_tile
#' @importFrom ggplot2 scale_fill_gradientn scale_x_continuous scale_y_continuous scale_y_date geom_text
#' @importFrom ggplot2 annotate geom_ribbon sec_axis element_text element_blank
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
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
#' data(treasury_yields)
#' dns_result <- fit_dns_model(treasury_yields, n_days = 30, show_progress = FALSE)
#' plot_yield_curves(treasury_yields, dns_result, max_plots = 4)
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

  # Maturity columns, located by name rather than by position
  mat_cols <- setdiff(names(yield_data), "Date")
  maturities <- suppressWarnings(as.numeric(mat_cols))

  # If no dates specified, use the first few dates common to both inputs.
  # intersect() drops the Date class, so restore it explicitly: formatting a
  # bare numeric date gives the day number, not a date.
  if (is.null(dates_to_plot)) {
    all_dates <- as.Date(
      intersect(as.Date(yield_data$Date), as.Date(dns_result$parameters$Date)),
      origin = "1970-01-01"
    )
    dates_to_plot <- all_dates[seq_len(min(max_plots, length(all_dates)))]
  }
  dates_to_plot <- as.Date(dates_to_plot)

  # Handle case with no dates to plot
  if (length(dates_to_plot) == 0) {
    warning("No matching dates found between yield data and DNS results")
    return(NULL)
  }

  # Prepare data for plotting
  curve_rows <- list()
  rmse_rows <- list()

  # For each date, add actual and fitted yields
  for (k in seq_along(dates_to_plot)) {
    this_date <- dates_to_plot[k]
    date_label <- format(this_date, "%Y-%m-%d")

    date_idx <- which(as.Date(yield_data$Date) == this_date)
    if (length(date_idx) == 0) next
    date_idx <- date_idx[1]

    actual_yields <- as.numeric(yield_data[date_idx, mat_cols])
    keep <- !is.na(actual_yields)
    if (any(keep)) {
      curve_rows[[length(curve_rows) + 1L]] <- data.frame(
        date = date_label,
        maturity = maturities[keep],
        yield = actual_yields[keep] * 100,  # Convert to percentage
        type = "Actual",
        stringsAsFactors = FALSE
      )
    }

    # Fitted yields (if available)
    dns_date_idx <- which(as.Date(dns_result$parameters$Date) == this_date)
    if (length(dns_date_idx) == 0) next
    dns_date_idx <- dns_date_idx[1]

    rmse_bp <- dns_result$parameters$rmse[dns_date_idx] * 10000  # basis points
    fitted_yields <- as.numeric(dns_result$fitted_yields[dns_date_idx, ])

    keep_fit <- !is.na(fitted_yields)
    if (any(keep_fit)) {
      curve_rows[[length(curve_rows) + 1L]] <- data.frame(
        date = date_label,
        maturity = maturities[keep_fit],
        yield = fitted_yields[keep_fit] * 100,  # Convert to percentage
        type = "Fitted",
        stringsAsFactors = FALSE
      )
    }

    # RMSE annotations live in their own layer, so they do not become a
    # colour level that the manual scale has to carry.
    if (show_rmse) {
      all_vals <- c(actual_yields, fitted_yields)
      if (any(!is.na(all_vals))) {
        rmse_rows[[length(rmse_rows) + 1L]] <- data.frame(
          date = date_label,
          maturity = max(maturities, na.rm = TRUE),
          yield = max(all_vals, na.rm = TRUE) * 100,
          label = sprintf("RMSE: %.1f bp", rmse_bp),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  plot_data <- if (length(curve_rows)) do.call(rbind, curve_rows) else data.frame()
  rmse_data <- if (length(rmse_rows)) do.call(rbind, rmse_rows) else NULL

  # Check if we have data to plot
  if (nrow(plot_data) == 0) {
    warning("No valid data to plot yield curves")
    return(NULL)
  }

  # Facet in calendar order, not alphabetically by label
  date_levels <- format(sort(unique(dates_to_plot)), "%Y-%m-%d")
  plot_data$date <- factor(plot_data$date, levels = date_levels)
  plot_data$date <- droplevels(plot_data$date)
  if (!is.null(rmse_data)) {
    rmse_data$date <- factor(rmse_data$date, levels = levels(plot_data$date))
  }

  # Create plot
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

  if (show_rmse && !is.null(rmse_data)) {
    p <- p + geom_text(
      data = rmse_data,
      aes(x = maturity, y = yield, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3,
      show.legend = FALSE
    )
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
#' data(treasury_yields)
#' dns_result <- fit_dns_model(treasury_yields, n_days = 250, show_progress = FALSE)
#' plots <- plot_factor_time_series(dns_result)
#' plots$factors
#' plots$rmse
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
    beta1 = "Level (\u03b2\u2081)",
    beta2 = "Slope (\u03b2\u2082)",
    beta3 = "Curvature (\u03b2\u2083)"
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
      beta1_change = "Level (\u03b2\u2081) Change",
      beta2_change = "Slope (\u03b2\u2082) Change",
      beta3_change = "Curvature (\u03b2\u2083) Change"
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
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' sim_result <- simulate_dns_factors(params, c(0.035, -0.015, 0.005),
#'                                    n_paths = 50, engine = "exact", seed = 1)
#' plots <- plot_simulated_factors(sim_result)
#' plots$level
#'
#' plot_simulated_factors(sim_result, plot_trajectories = TRUE)$level
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

  # Generate plots for each factor
  for (f in 1:3) {
    factor_name <- factor_names[f]

    if (plot_trajectories) {
      # Plot individual trajectories. Built in one shot: as.vector() on the
      # path slice is column-major, so time varies fastest within each path,
      # which is the order the per-path loop produced.
      plot_data <- data.frame(
        Time = rep(time, length(path_subset)),
        Value = as.vector(factors[, path_subset, f, drop = FALSE]),
        Path = rep(paste("Path", path_subset), each = length(time)),
        stringsAsFactors = FALSE
      )

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
      # Calculate summary statistics for each time point. The explicit
      # matrix() keeps a single-path simulation two-dimensional: factors[, , f]
      # drops to a vector when n_paths == 1 and apply() rejects it.
      factor_paths <- matrix(factors[, , f], nrow = n_times, ncol = n_paths)
      means <- apply(factor_paths, 1, mean)
      medians <- apply(factor_paths, 1, stats::median)
      q10 <- apply(factor_paths, 1, function(x) quantile(x, 0.1))
      q90 <- apply(factor_paths, 1, function(x) quantile(x, 0.9))

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
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' sim_result <- simulate_dns_yield_curves(params, c(0.035, -0.015, 0.005),
#'                                         maturities, n_paths = 50,
#'                                         engine = "exact", seed = 1)
#'
#' # Facets follow the calendar, so Year 10 comes after Year 5
#' plot_simulated_yields(sim_result, time_points = c(0, 1, 5, 10))$yield_curves
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

  # Facet labels, as a factor in time_points order. Left as character these
  # sort lexically, which puts "Year 10" between "Year 1" and "Year 5".
  time_labels <- paste("Year", round(time[time_indices], 1))
  time_levels <- unique(time_labels)

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

  # Plot yield curves at specified time points. expand.grid() varies its first
  # argument fastest, so this reproduces the maturity-inner, path-middle,
  # time-outer order of the nested loops it replaces.
  grid_y <- expand.grid(
    m = seq_along(maturities), p = seq_along(path_subset),
    t = seq_along(time_indices), KEEP.OUT.ATTRS = FALSE
  )
  plot_data <- data.frame(
    Time = factor(time_labels[grid_y$t], levels = time_levels),
    Maturity = maturities[grid_y$m],
    Yield = yields[cbind(time_indices[grid_y$t], grid_y$m,
                         path_subset[grid_y$p])] * 100,  # Convert to percentage
    Path = paste("Path", path_subset[grid_y$p]),
    stringsAsFactors = FALSE
  )

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

    # Create data for term structure evolution. Statistics are collected into
    # (time x maturity) matrices, then the frame is built once.
    stat_means <- matrix(NA_real_, n_times, length(maturities))
    stat_q10 <- stat_means
    stat_q90 <- stat_means

    for (m in seq_along(maturities)) {
      # matrix() keeps this two-dimensional when n_paths == 1
      yields_at_m <- matrix(yields[, m, ], nrow = n_times, ncol = n_paths)

      stat_means[, m] <- apply(yields_at_m, 1, mean)
      stat_q10[, m] <- apply(yields_at_m, 1, function(x) quantile(x, 0.1))
      stat_q90[, m] <- apply(yields_at_m, 1, function(x) quantile(x, 0.9))
    }

    grid_ts <- expand.grid(t = seq_len(n_times), m = seq_along(maturities),
                           KEEP.OUT.ATTRS = FALSE)
    term_data <- data.frame(
      Time = time[grid_ts$t],
      Maturity = paste(maturities[grid_ts$m], "yr"),
      Mean = as.vector(stat_means) * 100,
      Q10 = as.vector(stat_q10) * 100,
      Q90 = as.vector(stat_q90) * 100,
      stringsAsFactors = FALSE
    )

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

    # Prepare data, on the same grid as the yield panel above
    forward_data <- data.frame(
      Time = factor(time_labels[grid_y$t], levels = time_levels),
      Maturity = maturities[grid_y$m],
      Forward = forwards[cbind(time_indices[grid_y$t], grid_y$m,
                               path_subset[grid_y$p])] * 100,  # to percentage
      Path = paste("Path", path_subset[grid_y$p]),
      stringsAsFactors = FALSE
    )

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
#' plot_dns_loadings()
#'
#' # A faster decay concentrates the curvature loading at shorter maturities
#' plot_dns_loadings(lambda = 0.4)
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

  loadings <- ns_loadings(maturities, lambda)
  loading_data$Slope <- loadings[, 1]
  loading_data$Curvature <- loadings[, 2]

  # Convert to long format for plotting
  loading_long <- reshape2::melt(
    loading_data,
    id.vars = "Maturity",
    variable.name = "Factor",
    value.name = "Loading"
  )

  # Create plot
  p <- ggplot(loading_long, aes(x = Maturity, y = Loading, color = Factor)) +
    geom_line(linewidth = 1) +
    labs(
      title = "DNS Model Factor Loadings",
      subtitle = paste("\u03bb =", lambda),
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
      Label = c("Level (\u03b2\u2081)", "Slope (\u03b2\u2082)", "Curvature (\u03b2\u2083)")
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
#' data(treasury_yields)
#' plot_yield_heatmap(treasury_yields)
#'
#' # Most recent dates at the top
#' plot_yield_heatmap(treasury_yields, reverse_time = TRUE)
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
    # A centred second difference needs an interior point to sit on.
    if (length(maturities) < 3) {
      stop("plot_type = \"curvature\" needs at least 3 maturities; ",
           "yield_data has ", length(maturities), ".")
    }
    # Calculate curvature (2nd derivative approximation)
    yield_curvature <- yield_data
    for (i in seq.int(2, length(maturities) - 1)) {
      col <- names(yield_data)[i + 1]
      col_prev <- names(yield_data)[i]
      col_next <- names(yield_data)[i + 2]
      m <- maturities[i]
      m_prev <- maturities[i - 1]
      m_next <- maturities[i + 1]

      # Centered finite difference approximation
      # f''(x) ~= (f(x+h) - 2f(x) + f(x-h)) / h^2
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

  # Rescale to percent for display. The decision is made once, from the whole
  # frame: taken per column it could scale the short end and not the long end,
  # putting two different units on one heatmap.
  if (plot_type == "yields") {
    value_cols <- setdiff(names(plot_df), "Date")
    overall_max <- suppressWarnings(
      max(unlist(plot_df[value_cols]), na.rm = TRUE)
    )
    if (is.finite(overall_max) && overall_max < 1) {
      for (col in value_cols) plot_df[[col]] <- plot_df[[col]] * 100
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

  # Reverse the date axis if requested. scale_y_reverse() cannot be used here:
  # it negates the variable, and unary minus is undefined for Date. Reversed
  # limits on the coordinate system flip the axis and keep date labels.
  if (reverse_time) {
    date_range <- range(plot_data$Date, na.rm = TRUE)
    p <- p + ggplot2::coord_cartesian(ylim = rev(date_range))
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
#' # Needs plotly (interactive) or lattice (static)
#' data(treasury_yields)
#' plot_yield_surface(treasury_yields)
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
#' @param scale Logical. Whether to standardise each maturity's change series
#'   before the decomposition. \code{TRUE} is the default, preserving earlier
#'   behaviour. \code{FALSE} decomposes the covariance matrix instead, which is
#'   the usual convention for yield curves: it lets the more volatile
#'   maturities carry more weight, and is what gives the familiar level, slope
#'   and curvature shapes.
#'
#' @return A list containing:
#'   - pca: The PCA result object
#'   - plots: A list of created plots
#'
#' @examples
#' data(treasury_yields)
#'
#' # scale = FALSE decomposes the covariance matrix, the usual convention for
#' # yield curves, which is what produces the level/slope/curvature shapes
#' pca_result <- plot_yield_pca(treasury_yields, scale = FALSE)
#' round(pca_result$var_explained[1:3], 4)
#' pca_result$plots$weights
#'
#' @importFrom stats prcomp
#'
#' @export
plot_yield_pca <- function(yield_data, n_components = 3,
                           plot_weights = TRUE, plot_scores = TRUE,
                           plot_variance = TRUE, scale = TRUE) {
  # Parameter validation
  checkmate::assert_logical(scale, len = 1)
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
  pca_result <- prcomp(yield_changes, scale. = scale)

  # Limit to requested number of components
  n_components <- min(n_components, ncol(pca_result$rotation))

  # Computed once here: it is part of the return value, so it must exist even
  # when every plot_* flag is FALSE.
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

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
      geom_line(linewidth = 1) +
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
      geom_line(aes(y = Cumulative, group = 1), color = "red", linewidth = 1) +
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
#' @param models Character vector. Models to compare, from "NS", "NSS",
#'   "Svensson", "Spline". The NSS column needs the \pkg{SmithWilson} package,
#'   which is not on CRAN; without it that model is skipped with a warning.
#' @param spline_df Integer. Degrees of freedom for spline model.
#' @param units Character. Units for the returned fit metrics: \code{"bp"}
#'   (the default, basis points) or \code{"decimal"} (the units of the input
#'   yields).
#'
#' @return A list containing:
#'   - plot: A ggplot object with model comparisons
#'   - metrics: A data frame with RMSE, MAE and MaxError for each model, in
#'     \code{units}
#'
#' @details
#' Decimal yields scale to basis points by 10,000, matching
#' \code{\link{fit_dns_model}} and \code{\link{plot_factor_time_series}}.
#' Maturities with no observed yield on \code{date_to_compare} are excluded
#' from every model's fit and from the metrics.
#'
#' @examples
#' data(treasury_yields)
#'
#' # Errors are in basis points
#' comparison <- compare_term_structure_models(treasury_yields,
#'                                             models = c("NS", "Spline"))
#' comparison$metrics
#' comparison$plot
#'
#' @export
compare_term_structure_models <- function(yield_data, date_to_compare = NULL,
                                          models = c("NS", "NSS", "Spline"),
                                          spline_df = 5,
                                          units = c("bp", "decimal")) {
  # Parameter validation
  units <- match.arg(units)
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_subset(models, c("NS", "NSS", "Svensson", "Spline"))
  checkmate::assert_count(spline_df, positive = TRUE)

  # Decimal yields scale to basis points by 10,000.
  err_scale <- if (identical(units, "bp")) 10000 else 1

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

  # Extract yield data for the selected date. Maturities with no observation
  # are dropped up front: every fitter below would otherwise either error or
  # silently return a vector shorter than the maturity grid.
  mat_cols <- setdiff(names(yield_data), "Date")
  yields <- as.numeric(yield_data[date_idx[1], mat_cols])
  maturities <- suppressWarnings(as.numeric(mat_cols))

  observed <- !is.na(yields) & !is.na(maturities) & maturities > 0
  if (sum(observed) < 3) {
    stop("Fewer than 3 observed maturities on ", format(date_to_compare),
         "; nothing to compare.")
  }
  if (any(!observed)) {
    warning(sprintf("Dropping %d maturit%s with no observation on %s.",
                    sum(!observed), if (sum(!observed) == 1) "y" else "ies",
                    format(date_to_compare)))
  }
  yields <- yields[observed]
  maturities <- maturities[observed]

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

    # fit_ns_model() drops points it could not use, so place the fitted values
    # back onto the full maturity grid by index.
    ns_fitted <- rep(NA_real_, length(maturities))
    ns_fitted[ns_result$valid_idx] <- ns_result$fitted_yields
    model_fits$NS <- ns_fitted

    # Calculate fit metrics
    rmse <- ns_result$rmse
    mae <- mean(abs(yields - ns_fitted), na.rm = TRUE)
    max_error <- max(abs(yields - ns_fitted), na.rm = TRUE)

    fit_metrics <- rbind(
      fit_metrics,
      data.frame(
        Model = "Nelson-Siegel",
        RMSE = rmse * err_scale,
        MAE = mae * err_scale,
        MaxError = max_error * err_scale
      )
    )
  }

  # Fit Nelson-Siegel-Svensson model if requested
  if ("NSS" %in% models || "Svensson" %in% models) {
    # Check if SmithWilson package is available
    if (!requireNamespace("SmithWilson", quietly = TRUE)) {
      # Not listed in Suggests: SmithWilson is not available from CRAN, so
      # declaring it there makes the package's dependencies unsolvable. The
      # NSS comparison is simply unavailable unless it is installed by hand.
      warning("Package 'SmithWilson' is needed for the NSS model and is not ",
              "on CRAN; install it manually to include that column. Skipping.")
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
            RMSE = rmse * err_scale,
            MAE = mae * err_scale,
            MaxError = max_error * err_scale
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
        RMSE = rmse * err_scale,
        MAE = mae * err_scale,
        MaxError = max_error * err_scale
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
    geom_line(linewidth = 1) +
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
        fit_metrics$Model[i], ": RMSE = ",
        signif(fit_metrics$RMSE[i], 3),
        if (identical(units, "bp")) " bp\n" else "\n"
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
