context("validation")

### distr 
dd <- list(a="binomial",b="gaussian", c="multinomial",d="poisson")

test_that("dists",{
  expect_identical( abn:::validate_dists(c(a="bin",b="gaus",c="mu",d="pois")), dd)
  expect_error( abn:::validate_dists(c(a="bbin",b="gaus",c="mu",d="pois")))
})


### formula
test_that("formula", {
  expect_error( abn:::validate_abnDag(~a))   # data.df missing
  
  expect_equal( abn:::validate_abnDag(~a, c(a=1)), matrix(0,1,1, dimnames=list("a","a")))
  
  expect_error( abn:::validate_abnDag(~a|b+b|a, c(a=1,b=1)))  # 1-diag(2), cyclic     
})

m1 <- matrix( 0, 2,2, dimnames=list(c("a","b"), c("a","b")))
m2 <- matrix( 0, 4,4, dimnames=list(c("a","b","c","d"), c("a","b","c","d")))
m2[c(2,4,12,8)] <- 1

test_that("matrix", {
  expect_error( abn:::validate_abnDag( m1[2:1,]))
  expect_equal( abn:::validate_abnDag(m2), m2)
})


### all other cases

test_that("All other cases",{
  expect_error( abn:::validate_abnDag("~a"))  # string not a formula
  
})




context("Markov Blanket")
m2 <- matrix( 0, 4,4, dimnames=list(c("a","b","c","d"), c("a","b","c","d")))
m2[c(2,4,12,8)] <- 1

test_that("Trivial blankets",{
  
  expect_equal( sort(mb(m2, "c")), sort(c("a","b","d")))
  expect_error( mb(m2))
  
})



context("Compare Graph")
m2 <- matrix( 0, 4,4, dimnames=list(c("a","b","c","d"), c("a","b","c","d")))
m0 <- m2
m2[c(2,4,12,8)] <- 1

m2a <- create_abnDag( ~a+b|a+c+d|a:b:c, c(a=1,b=1,c=1,d=1))



#create_abnDag(~a+b+c+d, c(a=1,b=1,c=1,d=1))

test_that("part dag",{
  expect_equal(m2, create_abnDag( ~a+b|a+c+d|a:b:c, c(a=1,b=1,c=1,d=1))$dag)
  
  expect_equal(m0, create_abnDag( ~a+b+c+d, c(a=1,b=1,c=1,d=1))$dag)
  expect_equal(m0, create_abnDag( ~d+a+b+c, c(a=1,b=1,c=1,d=1))$dag)    # here name perturbation ok
  
  
})


test_that("part compareDag",{
  expect_equal( compareDag(m2, m2), compareDag(m2,   ~a+b|a+c+d|a:b:c, c(a=1,b=1,c=1,d=1)))  

 expect_equal( suppressWarnings(compareDag(m0, m0)), suppressWarnings(compareDag(m0,   ~a+b+c+d, c(a=1,b=1,c=1,d=1))))  
 expect_warning( compareDag(ref = m0,test =    ~a+b+c+d,node.names =  c(a=1,b=1,c=1,d=1)))
  
  
})

