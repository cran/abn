#Rdevel -d "valgrind --tool=memcheck --leak-check=full --log-file=test-valgrind.Rout" --no-save < test-valgrind.R 

require(abn, lib='../lib') # to have latest version ready...

mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
mydists <- list(b1="binomial",  b2="binomial",  b3="binomial", g1="gaussian",
                      b4="binomial",  p2="poisson", p4="poisson")
mydag <- matrix(0, nrow=7, ncol=7)
colnames(mydag) <- rownames(mydag) <- names(mydat)


mydag["b1","b2"] <- 1; # b1 <- b2
mydag["b2","p4"] <- 1; # b2 <- p4
mydag["b2","g1"] <- 1; # b2 <- g1
mydag["g1","p2"] <- 1; # g1 <- p2
mydag["b3","g1"] <- 1; # b3 <- g1
mydag["b4","b1"] <- 1; # b4 <- b1
mydag["p4","g1"] <- 1; # p4 <- g1

myres.inla <- fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                     compute.fixed=TRUE,max.mode.error=100);




