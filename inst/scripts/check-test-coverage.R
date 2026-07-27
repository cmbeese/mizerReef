# Test coverage report for mizerReef
#
# Runs the package's test suite under covr and reports coverage, both
# package-wide and per file. Intended for reviewers/readers checking the
# "code testing" claims in the paper's code-availability/checklist section,
# and for periodic maintainer use -- it's the same tool the
# .github/workflows/test-coverage.yaml CI job uses, run locally with a
# more readable summary.
#
# Usage: from the package root, `Rscript inst/scripts/check-test-coverage.R`
# (or source it interactively). Writes:
#   - a coverage percentage summary to the console (by file and overall)
#   - covr-test-coverage.html, an interactive line-by-line coverage report,
#     in the current working directory (gitignored, not part of the package)

if (!requireNamespace("covr", quietly = TRUE)) {
    stop("Package 'covr' is required. Install it with install.packages('covr').")
}

cov <- covr::package_coverage(
    path = ".",
    quiet = FALSE,
    clean = FALSE
)

cat("\n==================== Coverage by file ====================\n")
print(covr::coverage_to_list(cov)$filecoverage[order(covr::coverage_to_list(cov)$filecoverage)])

overall <- covr::percent_coverage(cov)
cat(sprintf("\n==================== Overall coverage: %.1f%% ====================\n", overall))

report_path <- "covr-test-coverage.html"
covr::report(cov, file = report_path, browse = interactive())
cat("\nLine-by-line HTML report written to:", normalizePath(report_path), "\n")
