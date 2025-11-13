# ==========================================
# ZZCOLLAB .Rprofile - Three-Part Structure
# ==========================================
# Part 1: User Personal Settings (from ~/.Rprofile)
# Part 2: renv Activation + Reproducibility Options
# Part 3: Auto-Snapshot on Exit
# ==========================================

# ==========================================
# Part 1: User Personal Settings
# ==========================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
q <- function(save="no", ...) quit(save=save, ...)

# Package installation behavior (non-interactive)
# Prevents prompts during install.packages()
options(
install.packages.check.source = "no",
install.packages.compile.from.source = "never",

# Parallel installation (faster package installs)
  Ncpus = parallel::detectCores()
)

# ==========================================
# Part 2: ZZCOLLAB Template - renv + Options
# ==========================================
# Activate renv (set project-local library paths)
# renv is pre-installed in Docker system library
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# renv consent (skips first-time prompts)
options(
  renv.consent = TRUE,
  renv.config.install.prompt = FALSE  # Skip "Do you want to proceed?" prompts
)

# Helper function for initializing renv without prompts
renv_init_quiet <- function() {
  renv::init(
    settings = list(snapshot.type = "explicit"),
    force = TRUE,
    restart = FALSE
  )
}

# ==========================================
# Auto-Initialize renv (New Projects)
# ==========================================
# Automatically initialize renv for new projects on first R session
# Only runs if:
# 1. renv not yet initialized (no renv.lock)
# 2. This is a project directory (has DESCRIPTION or mounted in container)
# 3. ZZCOLLAB_AUTO_INIT not disabled
#
# Disable with: docker run -e ZZCOLLAB_AUTO_INIT=false ...

if (!file.exists("renv.lock")) {
  # Check if auto-init is enabled (default: true)
  auto_init <- Sys.getenv("ZZCOLLAB_AUTO_INIT", "true")

  # Check if this looks like a project directory
  is_project <- file.exists("DESCRIPTION") ||
                getwd() == "/home/analyst/project"  # Container path

  if (tolower(auto_init) %in% c("true", "t", "1") && is_project) {
    message("\n🔧 ZZCOLLAB: Auto-initializing renv for new project...")

    tryCatch({
      renv_init_quiet()
      message("✅ renv initialized successfully")
      message("   Install packages with: install.packages('package')")
    }, error = function(e) {
      warning("⚠️  Auto-init failed: ", conditionMessage(e),
              "\n   Run manually: renv_init_quiet()", call. = FALSE)
    })
  } else if (!is_project) {
    message("\n💡 Tip: Initialize renv with renv_init_quiet()")
  }
}

# ==========================================
# Critical Reproducibility Options
# See: docs/COLLABORATIVE_REPRODUCIBILITY.md Pillar 3
# ==========================================
# These options affect computational results and should not be modified
# without team review. Changes are monitored by check_rprofile_options.R

options(
  # Character vector treatment in data.frames
  # FALSE ensures characters stay as characters (R >= 4.0.0 default)
  stringsAsFactors = FALSE,

  # Statistical contrasts for factor variables in models
  # Treatment contrasts for unordered factors, polynomial for ordered
  contrasts = c("contr.treatment", "contr.poly"),

  # Missing data handling in modeling functions
  # na.omit removes rows with any NA values
  na.action = "na.omit",

  # Numeric precision in printed output
  # 7 significant digits (R default)
  digits = 7,

  # Decimal separator for output
  # Period (US standard) ensures consistency across locales
  OutDec = ".",

  # CRAN mirror for package installation
  repos = c(CRAN = "https://cloud.r-project.org"),

  # Package installation behavior (non-interactive)
  # Prevents prompts during install.packages()
  install.packages.check.source = "no",
  install.packages.compile.from.source = "never",

  # Parallel installation (faster package installs)
  Ncpus = parallel::detectCores()
)

# ==========================================
# Auto-Snapshot on R Exit
# ==========================================
# Automatically updates renv.lock when exiting R session
# This captures any packages installed during the session

.Last <- function() {
  # Check if auto-snapshot is enabled (default: true)
  auto_snapshot <- Sys.getenv("ZZCOLLAB_AUTO_SNAPSHOT", "true")

  if (tolower(auto_snapshot) %in% c("true", "t", "1")) {
    # Check if we're in an renv project
    if (file.exists("renv.lock") && file.exists("renv/activate.R")) {
      message("\n📸 Auto-snapshot: Updating renv.lock...")

      snapshot_result <- tryCatch({
        # Snapshot with prompt disabled (non-interactive)
        # Uses default snapshot type (captures all installed packages)
        renv::snapshot(prompt = FALSE)
        TRUE
      }, error = function(e) {
        warning("Auto-snapshot failed: ", conditionMessage(e), call. = FALSE)
        FALSE
      })

      if (snapshot_result) {
        message("✅ renv.lock updated successfully")
        message("   Commit changes: git add renv.lock && git commit -m 'Update packages'")
      }
    }
  }

  # Call any user-defined .Last function from .Rprofile.local
  if (exists(".Last.user", mode = "function", envir = .GlobalEnv)) {
    tryCatch(
      .Last.user(),
      error = function(e) warning("User .Last failed: ", conditionMessage(e))
    )
  }
}

# ==========================================
# Personal/Team Customizations (Optional)
# ==========================================
# Load personal settings from git-ignored file
# This allows team members to have personal preferences without
# affecting version-controlled reproducibility settings

if (file.exists(".Rprofile.local")) {
  tryCatch(
    source(".Rprofile.local"),
    error = function(e) {
      warning(".Rprofile.local failed to load: ", conditionMessage(e))
    }
  )
}
