# Database Connection Setup
# Template for connecting to various database systems

# Load database packages
suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(RPostgres)  # For PostgreSQL
  library(RMySQL)     # For MySQL
  library(odbc)       # For ODBC connections
  library(here)
})

# 1. SQLite Database (local file-based) ====
setup_sqlite <- function(db_path = here("analysis", "data", "project.db")) {
  con <- dbConnect(RSQLite::SQLite(), db_path)
  cat("Connected to SQLite database at:", db_path, "\\n")
  return(con)
}

# 2. PostgreSQL Database ====
setup_postgresql <- function() {
  con <- dbConnect(RPostgres::Postgres(),
    dbname = Sys.getenv("POSTGRES_DB", "research_db"),
    host = Sys.getenv("POSTGRES_HOST", "localhost"),
    port = as.numeric(Sys.getenv("POSTGRES_PORT", "5432")),
    user = Sys.getenv("POSTGRES_USER"),
    password = Sys.getenv("POSTGRES_PASSWORD")
  )
  cat("Connected to PostgreSQL database\\n")
  return(con)
}

# 3. MySQL Database ====
setup_mysql <- function() {
  con <- dbConnect(RMySQL::MySQL(),
    dbname = Sys.getenv("MYSQL_DB", "research_db"),
    host = Sys.getenv("MYSQL_HOST", "localhost"),
    port = as.numeric(Sys.getenv("MYSQL_PORT", "3306")),
    user = Sys.getenv("MYSQL_USER"),
    password = Sys.getenv("MYSQL_PASSWORD")
  )
  cat("Connected to MySQL database\\n")
  return(con)
}

# 4. ODBC Connection (for various databases) ====
setup_odbc <- function(dsn_name) {
  con <- dbConnect(odbc::odbc(), dsn = dsn_name)
  cat("Connected via ODBC to:", dsn_name, "\\n")
  return(con)
}

# Example usage:
# con <- setup_sqlite()  # For local SQLite database
# con <- setup_postgresql()  # For PostgreSQL (requires environment variables)
# dbDisconnect(con)  # Always disconnect when done

cat("Database connection functions loaded\\n")
cat("Set environment variables for database credentials\\n")
cat("Use setup_sqlite(), setup_postgresql(), setup_mysql(), or setup_odbc()\\n")
