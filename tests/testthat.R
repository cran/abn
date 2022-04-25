Sys.setenv("R_TESTS" = "")

# automatic unit tests

library("testthat")
# library("abn", lib='../../lib')

test_check("abn")
