# High-Performance Computing Setup
# Configure parallel processing for computationally intensive tasks

# Load required packages
library(parallel)
library(doParallel)
library(foreach)
library(future)
library(furrr)

# Detect available cores
n_cores <- parallel::detectCores()
cat("Detected", n_cores, "CPU cores\\n")

# Set up parallel backend (choose one)

# 1. Using doParallel (for foreach loops)
cl <- makeCluster(max(1, n_cores - 1))  # Leave one core free
registerDoParallel(cl)

# 2. Using future (for purrr-style functions)
plan(multisession, workers = max(1, n_cores - 1))

# Configuration for different computing environments

# Local development (conservative)
if (Sys.getenv("COMPUTING_ENV") == "local" || Sys.getenv("COMPUTING_ENV") == "") {
  n_workers <- min(4, max(1, n_cores - 1))
  options(mc.cores = n_workers)
  cat("Local environment: Using", n_workers, "cores\\n")
}

# High-performance cluster
if (Sys.getenv("COMPUTING_ENV") == "cluster") {
  # Read from SLURM environment variables or config
  n_workers <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", n_cores))
  options(mc.cores = n_workers)
  cat("Cluster environment: Using", n_workers, "cores\\n")
}

# Cloud computing (AWS, GCP, Azure)
if (Sys.getenv("COMPUTING_ENV") == "cloud") {
  # Configure based on instance type
  n_workers <- n_cores  # Use all available cores in cloud
  options(mc.cores = n_workers)
  cat("Cloud environment: Using", n_workers, "cores\\n")
}

# Memory management for large datasets
if (Sys.getenv("LARGE_DATA") == "true") {
  # Increase memory limits
  if (.Platform$OS.type == "unix") {
    system("ulimit -v unlimited", ignore.stderr = TRUE)
  }
  
  # Configure garbage collection
  options(expressions = 500000)  # Increase expression limit
  
  # Use memory-efficient data structures
  options(datatable.fwrite.sep = ",")
  options(datatable.optimize = 2)
  
  cat("Large data mode: Optimized for memory efficiency\\n")
}

# Progress reporting setup
options(future.progress = TRUE)

# Cleanup function
cleanup_parallel <- function() {
  if (exists("cl")) {
    stopCluster(cl)
  }
  plan(sequential)  # Reset future plan
}

# Register cleanup on exit
on.exit(cleanup_parallel(), add = TRUE)

cat("Parallel computing setup complete\\n")
