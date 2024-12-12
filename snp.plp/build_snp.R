# Step 1: Unload and clean any existing package installation
clean_environment <- function(package_name) {
  if (paste0("package:", package_name) %in% search()) {
    detach(paste0("package:", package_name), unload = TRUE, force = TRUE)
  }
  # Remove the package from the library
  if (package_name %in% installed.packages()) {
    remove.packages(package_name)
  }
}

# Step 2: Clear old compiled files
clean_compiled_files <- function() {
  if (dir.exists("src")) {
    tryCatch({
      devtools::clean_dll()
    }, error = function(e) {
      message("Error during cleaning compiled files: ", e$message)
    })
  }
}

# Step 3: Regenerate Rcpp attributes and documentation
regenerate_attributes_and_docs <- function() {
  tryCatch({
    Rcpp::compileAttributes()
    roxygen2::roxygenize()
    devtools::document()
  }, error = function(e) {
    message("Error during attribute or documentation generation: ", e$message)
  })
}

# Step 4: Build and install the package
build_and_install_package <- function(package_name) {
  tryCatch({
    devtools::install(force = TRUE, build_vignettes = FALSE, upgrade = "never")
    message("Package '", package_name, "' built and installed successfully!")
  }, error = function(e) {
    message("Error during package build or installation: ", e$message)
  })
}

# Step 5: Test package loading and a function
test_package <- function(package_name, test_function) {
  tryCatch({
    library(package_name, character.only = TRUE)
    message("Package '", package_name, "' loaded successfully!")
    if (!missing(test_function)) {
      do.call(test_function, list())
      message("Test function '", test_function, "' executed successfully!")
    }
  }, error = function(e) {
    message("Error during package test: ", e$message)
  })
}

# Main execution
tryCatch({
  package_name <- "snp.plp"  # Replace with your package name
  test_function <- "htslib_version"     # Replace with a function from your package
  
  if ("package:snp.plp" %in% search()) {
    detach("package:snp.plp", unload = TRUE, force = TRUE)
  }
  if ("snp.plp" %in% installed.packages()[, "Package"]) {
    remove.packages("snp.plp")
  }
  

  # Clean up the environment (do not delete all user variables, only package-specific)
  clean_environment(package_name)
  
  # Clean compiled files
  clean_compiled_files()
  
  # Regenerate attributes and documentation
  regenerate_attributes_and_docs()
  
  # Build and install the package
  build_and_install_package(package_name)
  
  # Test the package
  test_package(package_name, test_function)
}, error = function(e) {
  message("An error occurred during the process: ", e$message)
})

