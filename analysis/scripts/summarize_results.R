# ============================================
# Summarize Simulation Results
# Reads and summarizes full_pmsim_analysis_hyb_versus_co.RData
# ============================================

library(tidyverse)

# Load the results
cat("Loading results from ../output/full_pmsim_analysis_hyb_versus_co.RData\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

load("../output/full_pmsim_analysis_hyb_versus_co.RData")

# List all loaded objects
loaded_objects <- ls()
cat("Objects loaded:\n")
for (obj_name in loaded_objects) {
  obj <- get(obj_name)
  cat(sprintf("  - %s: %s\n", obj_name, class(obj)[1]))
}
cat("\n")

# Function to summarize an object
summarize_object <- function(obj_name) {
  obj <- get(obj_name)

  cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
  cat(sprintf("Object: %s\n", obj_name))
  cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

  # Basic info
  cat(sprintf("Class: %s\n", paste(class(obj), collapse = ", ")))

  if (is.data.frame(obj) || is_tibble(obj)) {
    cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(obj), ncol(obj)))
    cat("\nColumn names:\n")
    cat(paste("  -", names(obj), collapse = "\n"), "\n")

    cat("\nColumn types:\n")
    for (col in names(obj)) {
      cat(sprintf("  - %s: %s\n", col, class(obj[[col]])[1]))
    }

    cat("\nFirst few rows:\n")
    print(head(obj))

    # If it looks like simulation results, provide more detail
    if ("power" %in% names(obj) || "Power" %in% names(obj)) {
      cat("\n--- Power Analysis Summary ---\n")

      # Find power column (case-insensitive)
      power_col <- names(obj)[tolower(names(obj)) == "power"][1]

      if (!is.na(power_col)) {
        cat(sprintf("\nPower statistics:\n"))
        cat(sprintf("  Mean:   %.3f\n", mean(obj[[power_col]], na.rm = TRUE)))
        cat(sprintf("  Median: %.3f\n", median(obj[[power_col]], na.rm = TRUE)))
        cat(sprintf("  Min:    %.3f\n", min(obj[[power_col]], na.rm = TRUE)))
        cat(sprintf("  Max:    %.3f\n", max(obj[[power_col]], na.rm = TRUE)))
        cat(sprintf("  SD:     %.3f\n", sd(obj[[power_col]], na.rm = TRUE)))
      }

      # Group by design if present
      if ("design" %in% names(obj) || "design_name" %in% names(obj)) {
        design_col <- if ("design" %in% names(obj)) "design" else "design_name"

        cat("\nPower by design:\n")
        obj %>%
          group_by(!!sym(design_col)) %>%
          summarize(
            n = n(),
            mean_power = mean(!!sym(power_col), na.rm = TRUE),
            sd_power = sd(!!sym(power_col), na.rm = TRUE),
            .groups = "drop"
          ) %>%
          print()
      }

      # Group by carryover if present
      if ("carryover_t1half" %in% names(obj)) {
        cat("\nPower by carryover condition:\n")
        obj %>%
          group_by(carryover_t1half) %>%
          summarize(
            n = n(),
            mean_power = mean(!!sym(power_col), na.rm = TRUE),
            .groups = "drop"
          ) %>%
          print()
      }

      # Group by biomarker correlation if present
      if ("biomarker_correlation" %in% names(obj)) {
        cat("\nPower by biomarker correlation:\n")
        obj %>%
          group_by(biomarker_correlation) %>%
          summarize(
            n = n(),
            mean_power = mean(!!sym(power_col), na.rm = TRUE),
            .groups = "drop"
          ) %>%
          print()
      }

      # Check for model_carryover comparison
      if ("model_carryover" %in% names(obj)) {
        cat("\nPower by modeling approach:\n")
        obj %>%
          group_by(model_carryover) %>%
          summarize(
            n = n(),
            mean_power = mean(!!sym(power_col), na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(approach = if_else(model_carryover, "WITH carryover", "WITHOUT carryover")) %>%
          select(approach, n, mean_power) %>%
          print()
      }
    }

  } else if (is.list(obj)) {
    cat(sprintf("Length: %d elements\n", length(obj)))
    if (!is.null(names(obj))) {
      cat("\nElement names:\n")
      cat(paste("  -", names(obj)[1:min(10, length(names(obj)))], collapse = "\n"), "\n")
      if (length(names(obj)) > 10) {
        cat(sprintf("  ... and %d more\n", length(names(obj)) - 10))
      }
    }

  } else if (is.matrix(obj)) {
    cat(sprintf("Dimensions: %d x %d\n", nrow(obj), ncol(obj)))
    cat("\nFirst few rows/cols:\n")
    print(obj[1:min(5, nrow(obj)), 1:min(5, ncol(obj))])

  } else if (is.vector(obj)) {
    cat(sprintf("Length: %d\n", length(obj)))
    cat(sprintf("Type: %s\n", typeof(obj)))
    if (length(obj) <= 10) {
      cat("Values:", paste(obj, collapse = ", "), "\n")
    } else {
      cat("First 10 values:", paste(head(obj, 10), collapse = ", "), "...\n")
    }

  } else {
    cat("Structure:\n")
    str(obj)
  }

  cat("\n")
}

# Summarize each object
for (obj_name in loaded_objects) {
  summarize_object(obj_name)
}

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("Summary complete.\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
