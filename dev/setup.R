## AXtools — Developer Workflow
## Run these after editing functions in R/
## Source this file or run sections interactively.

# ── Document & reload ────────────────────────────────────────────────────────
# Run after adding/editing roxygen2 headers or exported functions.
# Updates NAMESPACE and man/ pages from roxygen comments.
devtools::document()

# Reload the package in your current session (no reinstall needed).
devtools::load_all()

# ── Check ────────────────────────────────────────────────────────────────────
# Full R CMD check — run before pushing to GitHub.
devtools::check()

# ── Install locally ──────────────────────────────────────────────────────────
# Install from your local source so library(axtools) picks up changes.
devtools::install()

# ── Quick reference ──────────────────────────────────────────────────────────
# Typical edit cycle:
#   1. Edit R/*.R (add/change functions and roxygen headers)
#   2. devtools::document()
#   3. devtools::load_all()
#   4. Test interactively
#   5. devtools::check() before committing
