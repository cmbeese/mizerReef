# Developer tools for building, checking, and documenting the mizerReef package

# Load required packages for development
library(devtools) # Tools for package development (build, install, check, document, etc.)
library(pkgdown) # For building a static website from package documentation

# Generate documentation from roxygen2 comments
# This updates the man/ directory and NAMESPACE file based on your roxygen2 comments in R/ scripts.
devtools::document()

# Build a source tarball of the package (like R CMD build)
# This creates a .tar.gz file in the parent directory, which can be distributed or submitted to CRAN
devtools::build()

# (Re-)generate documentation again in case build modified any files (optional, but sometimes useful)
devtools::document()

# Run package checks (like R CMD check)
# This checks for errors, warnings, and notes, and runs all tests and examples.
devtools::check()

# Run all tests
# This runs all tests in the tests folder
devtools::test()

# Build the pkgdown website
# This generates a static website (in /docs) from your documentation, vignettes, and README.
pkgdown::build_site()

# Installing for local testing
devtools::install() # Install the package in your local R library
remove.packages("mizerReef") # Unload the package if already loaded
