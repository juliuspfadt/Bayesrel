# Bayesrel - R Session Startup for Claude Code
# Run this script in your interactive R session (RStudio/Positron/radian)
# to prepare and hand over the session to Claude Code.
#
# Usage: source(".claude/session_startup.R")

# Fix cli::get_spinner() conflict with testthat/tinytest in btw/evaluate context.
# BOTH the option AND the monkey-patch are needed:
#   - options(cli.spinner = "line") sets a sane default
#   - But btw:::local_reproducible_output() overrides it to FALSE on every
#     btw_tool_run_r() call, causing cli::get_spinner() to hit:
#       FALSE$frames -> "$ operator is invalid for atomic vectors"
#   - The monkey-patch intercepts logical values and coerces to "line"
# Upstream fixes pending: posit-dev/btw and r-lib/cli
options(cli.spinner = "line")
local({
  original_get_spinner <- cli::get_spinner
  patched_get_spinner <- function(which = NULL) {
    if (is.null(which)) {
      opt <- getOption("cli.spinner")
      if (identical(opt, FALSE) || identical(opt, TRUE)) {
        options(cli.spinner = "line")
      }
    } else if (identical(which, FALSE) || identical(which, TRUE)) {
      which <- "line"
    }
    original_get_spinner(which = which)
  }
  utils::assignInNamespace("get_spinner", patched_get_spinner, ns = "cli")
})

# MCP server dependencies
for (pkg in c("btw", "mcptools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# Load Bayesrel from source (recompiles src/ when C++ changed).
# Use pkgbuild::clean_dll() first if you hit stale-object link errors.
pkgload::load_all(".", export_all = FALSE)

# Hand the session over to Claude Code.
btw::btw_mcp_session()
