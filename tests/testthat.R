Sys.setenv("R_TESTS" = "")

# automatic unit tests

library("testthat")
# library("abn")

test_check("abn")
