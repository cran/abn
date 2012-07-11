##################################################################################################################
## Estimate all posterior estimates from a given DAG
##################################################################################################################
# Step 1. Define the DAG
##################################################################################################################
library(abn);
# define a DAG and call it mydag
mydag<-matrix(c(
# D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 Year Loc.x Loc.y
  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D1
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D2
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D3
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   1,    0,    # D4
  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0,   0,    0,    # D5
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D6
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,   0,    0,    # D7
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0,   0,    0,    # D8
  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D9
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0,   0,    0,    # D10
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # Year
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,    0,    # Loc.x
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,    0     # Loc.y
                                               ),byrow=TRUE,ncol=13); 
colnames(mydag)<-rownames(mydag)<-names(pigs.1par);#set names
##################################################################################################################
# Step 2. Go through all the different parameters and calculate and store the densities
# note that the ranges (post.x) were all set manually with some trial and error
##################################################################################################################
#D1|D2
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D1",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=1,to=1.7,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D1.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D1",
                   whichvar="D2",#this is the intercept
                     post.x=seq(from=0.6,to=1.5,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# add to existing grid x,fx,x,fx
D1.p<-cbind(D1.p,marg[,1],marg[,2]);
##################################################################
#D2|D3
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D2",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=0.6,to=1.1,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D2.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D2",
                   whichvar="D3",#this is the intercept
                     post.x=seq(from=1.2,to=1.9,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D2.p<-cbind(D2.p,marg[,1],marg[,2]);
#################################################################
#D3|D4
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D3",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=0.5,to=0.8,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D3.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D3",
                   whichvar="D4",#this is the intercept
                     post.x=seq(from=0.7,to=1.5,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D3.p<-cbind(D3.p,marg[,1],marg[,2]);
#################################################################
#D4|Loc.x
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D4",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-1.8,to=-1.4,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D4.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D4",
                   whichvar="Loc.x",#this is the intercept
                     post.x=seq(from=0.25,to=0.6,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D4.p<-cbind(D4.p,marg[,1],marg[,2]);
#################################################################
#D5|D6
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D5",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-1.7,to=-1.3,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D5.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D5",
                   whichvar="D6",#this is the intercept
                     post.x=seq(from=0.5,to=1.0,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D5.p<-cbind(D5.p,marg[,1],marg[,2]);
##################################################################
#D6|D4
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D6",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-0.05,to=0.2,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D6.p<-cbind(marg[,1],marg[,2]);
#D6|D4
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D6",
                   whichvar="D4",#this is the intercept
                     post.x=seq(from=0.8,to=1.4,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D6.p<-cbind(D6.p,marg[,1],marg[,2]);
###################################################################
#D7|Year
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D7",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-1.6,to=-1.3,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D7.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D7",
                   whichvar="Year",#this is the intercept
                     post.x=seq(from=-0.4,to=-0.1,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D7.p<-cbind(D7.p,marg[,1],marg[,2]);
######################################################################


#D8|D10
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D8",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-2.8,to=-2.2,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D8.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D8",
                   whichvar="D10",#this is the intercept
                     post.x=seq(from=1,to=2,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D8.p<-cbind(D8.p,marg[,1],marg[,2]);
######################################################################
#D9|D2
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D9",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-3,to=-2,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D9.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D9",
                   whichvar="D2",#this is the intercept
                     post.x=seq(from=0.8,to=1.8,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D9.p<-cbind(D9.p,marg[,1],marg[,2]);
#######################################################################
#D10|D9
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D10",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-3,to=-2.4,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D10.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="D10",
                   whichvar="D9",#this is the intercept
                     post.x=seq(from=0.5,to=1.2,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
D10.p<-cbind(D10.p,marg[,1],marg[,2]);
########################################################################
#Year
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Year",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=-0.05,to=0.05,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
Year.p<-cbind(marg[,1],marg[,2]);

#Year.pec
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Year",
                   whichvar="precision",#this is the intercept
                     post.x=seq(from=0.92,to=1.07,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
prec.Year.p<-cbind(marg[,1],marg[,2]);
##########################################################################
#Loc.x
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Loc.x",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=0,to=0.15,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
Loc.x.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Loc.x",
                   whichvar="D7",#this is the intercept
                     post.x=seq(from=-0.5,to=-0.2,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
Loc.x.p<-cbind(Loc.x.p,marg[,1],marg[,2]);
#Loc.x.prec
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Loc.x",
                   whichvar="precision",#this is the intercept
                     post.x=seq(from=0.95,to=1.1,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
prec.Loc.x.p<-cbind(marg[,1],marg[,2]);
#############################################################################
#Loc.y
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Loc.y",
                   whichvar="constant",#this is the intercept
                     post.x=seq(from=0.02,to=0.15,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
Loc.y.p<-cbind(marg[,1],marg[,2]);
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Loc.y",
                   whichvar="D7",#this is the intercept
                     post.x=seq(from=-0.6,to=-0.3,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
Loc.y.p<-cbind(Loc.y.p,marg[,1],marg[,2]);
#Loc.y.prec
marg<-getmarginal(data.df=pigs.1par,
                   dag.m=mydag,
                   whichnode="Loc.y",
                   whichvar="precision",#this is the intercept
                     post.x=seq(from=0.95,to=1.1,len=1000),
                     verbose=FALSE,std=TRUE);

#plot(marg[,1],marg[,2],type="l");
# get grid of discrete values d and p(d)
prec.Loc.y.p<-cbind(marg[,1],marg[,2]);
###############################################################################################
# Step 3. dump all the densities into a format which can be read be JAGS
# dump this into a file called post_params.R
###############################################################################################
## now dump it all to a format which JAGS then import
Nints<-1000;
dump(c("Nints","D1.p","D2.p","D3.p","D4.p","D5.p","D6.p","D7.p","D8.p","D9.p","D10.p",
       "Year.p","prec.Year.p","Loc.x.p","prec.Loc.x.p","Loc.y.p","prec.Loc.y.p"),file="post_params.R");

