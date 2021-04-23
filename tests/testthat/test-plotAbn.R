
context("plotAbn")
suppressPackageStartupMessages(require(testthat))
suppressPackageStartupMessages(require(abn))
suppressPackageStartupMessages(require(Rgraphviz))


# construct example
dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="binomial")
edge.strength <- matrix(c(0,0.5,0.5,0.7,0.1,0,   #Define a matrix formulation
                          0,0,0.3,0.1,0,0.8,
                          0,0,0,0.35,0.66,0,
                          0,0,0,0,0.9,0,
                          0,0,0,0,0,0.8,
                          0,0,0,0,0,0),nrow = 6L, ncol = 6L, byrow = TRUE)
colnames(edge.strength) <- rownames(edge.strength) <- names(dist)  #Naming of the matrix

#Plot from a formula
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.fillcolor.list= "e")
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.fillcolor.list= "e", markov.blanket.node = "b")
# mb has precedance!

# with diamond
expect_warning( plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.shape=c('diamond','box','circle')))


#Plot form a matrix
expect_warning( plotabn(dag = edge.strength, data.dist = dist))

# No plotting:
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, plot = FALSE)

# Markov blanket
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = "e")
mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = "e")

plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = "c")
mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = "c")
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = c("d"))
mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = "d")
plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = c("d","f"))
mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = c("d","f"))




# edge.strength
plotAbn(edge.strength, edge.strength=edge.strength, data.dist = dist)
vals <- c(t(edge.strength))
expect_error(plotAbn(edge.strength, edge.strength=edge.strength, data.dist = dist,
                     fitted.values=vals[vals>0]), "argument is of length zero")
tmp <- edge.strength
tmp[1,6] <- 1
expect_error(plotAbn(tmp, edge.strength=-edge.strength-4, data.dist = dist),
  "'edge.strength' should be positive")
plotAbn(tmp, edge.strength=edge.strength, data.dist = dist) # zeros allowd
expect_error(plotAbn(edge.strength, edge.strength=tmp, data.dist = dist),"'edge.strength' does not match dag")



# testing fitted.values
mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
mydists <- list(b1="binomial", b2="binomial", b3="binomial", g1="gaussian",
                b4="binomial", p2="poisson", p4="poisson")
mydag.empty <- matrix(0, nrow=7, ncol=7)
colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat)
(myres <- fitAbn(dag = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1+p4|g1, data.df = mydat, data.dists = mydists))
g <- plotAbn(myres$abnDag$dag, fitted.values = myres$modes, data.dists = mydists, edge.direction = 'pc')
myres

(myres <- fitAbn(dag = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1:g1+p4|g1:b3:p2, data.df = mydat, data.dists = mydists))
g <- plotAbn(myres$abnDag$dag, fitted.values = myres$modes, data.dists = mydists, edge.direction = 'pc')
myres


mydat1 <- cbind(m1=as.factor(as.numeric(mydat[,1])*2-as.numeric(mydat[,2])),
                mydat[3:7])
mydists1 <- list(m1="multinomial", b3="binomial", g1="gaussian",
                b4="binomial", p2="poisson", p4="poisson")
# cored dump!!
(myres1 <- fitAbn(dag = ~m1|b3:g1:b4+b3|p2:p4+p2|p4, data.df = mydat1, data.dists = mydists1, method='mle'))
plotAbn(dag = ~m1|b3:g1:b4+b3|p2:p4+p2|p4, data.dist=mydists1, node.shape=rep('box',4))

#  str(AgEdge(g))[[1]]
#
#  getMethod('plot', 'Ragraph')  # entire graph
#  getMethod('lines', 'AgEdge')  # single edge
## labels are done with `drawTxtLabel`
# getMethod("labelLoc",'AgTextLabel')


## PLAIN WRONG
## abn:::plotabn(myres$abnDag$dag, fitted.values.abn.mle = myres$modes, data.dist = mydists)






