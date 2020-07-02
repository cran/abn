Sys.setenv("R_TESTS" = "")

# automatic unit tests

library("testthat")
# library("abn")
# abn cannot be tested on old version
if (as.numeric(sessionInfo()$R.version$major)>=3 & as.numeric(sessionInfo()$R.version$minor)>6.3){
test_check("abn")
}