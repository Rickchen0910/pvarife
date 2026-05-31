# Run all tests under a sequential plan so tests are deterministic
# and do not accidentally use parallel workers.
# This mirrors the convention in the pgfe package.
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}
