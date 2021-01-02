##################  DMC Lesson 4: Multiple Subjects


### Lesson 4.8 Fitting real data

###################### GETTING REAL DATA INTO DMC ##############################
#
# First, create a model and simulate data for exactly the same design as your
# experiment. Note that the simulated data doesnt have to look like your real
# data you have to fit the data first to figure out the parameters (if any) that
# can do that. At this point all you want to do is make sure you design is 
# identical.

# Now setup your data frame like the output of h.simulate.dmc.
#
# 1) Make sure you have an "s" column identifying the data for each subject that
#    is a factor
# 2) Follow that by columns specifying the factors in your design. Again each 
#    must be a factor
# 3) The second last column give the response (again a factor)
# 4) The last columns is called RT, make sure the unit is seconds.
#
# Column names and factor levels for (2) and (3) depend on your model 
# specification so make sure they are compatible. data.model.dmc attempts to 
# check this.
#
# Note that model.dmc is very picky about overlap in factor levels, particularly
# where one level is a substring of another. This means you should avoid short
# and cryptic level names (this also helps with reading output!).
#
# Avoid superfluous columns, they are usually ignored but can sometimes cause
# problems.
#
# Make sure that your factors not only have the same levels but that the levels
# are in the same order as specified in the model (including the response 
# factor).
#
# See advanced lesson dmc_5_1_Factorial.R for setting up models for the complex
# designs that often occur in real experiments.
# 
#
############ A SCRIPT (THAT TRYS) TO AUTOMATICALLY FIT YOUR DATA ###############

# With the usual caveats that automatic procedures dont always work you can 
# often use the following script first fit fixed effects then use it as a basis
# for random effects. The script automatically constructs hyper mu priors 
# from the data level priors used in the fixed effect fits and uses exponentials
# with a scale of 1 for the the hyper sigma priors. These choices (particularly
# the scale parameters of the mu and sigma priors) might need tweaking.
#
# The following code assumes and lba_B model and that you have a data.model 
# object called "dmX" and have created an object ready to fit called sX with 
# h.samples.dmc (I recommend you set nmc=120 and dont thin, the script then sets 
# thin=10 for the hierarchical fitting, again this might need tweaking) that is 
# saved to a file called modelmodelX.RData, e.g.,
# 
# save(dmX,sX,file="modelmodelX.RData")
#
# You can save the following code off into a script and then run it on the
# command line. For example on a Linux system with a script called modelX.R
#
# nohup R CMD BATCH modelX.R &
#
# The "nohup" means you can log out of the terminal and it will keep running. 
# The "&" means it runs in background so you get the command prompt back. The
# console output that is created will be saved to a file "modelX.Rout" which 
# you can look at to see how your if is progressing. Note that you can change
# all of these names (dmX, sX, modelX) as required. You can fit several
# different models conveniently by just saving the files with new names and 
# searching and replaceing "X" with corresponding model designators in the 
# script.
#
# Note that the script assumes that you have at least as many cores as subjects.
# If you dont then edit the line "cores=length(sX)".
#
# Note also that the results of individual fitting and the data.model object are 
# saved off to "donesmodelX.RData" and then deleted in order to minimize memory 
# usage when the hierarhcial run starts.  
#
# After h.converge.dmc finishes the hierarhcial fits (doing an extra final run
# after convergence is satistifed with nmc=250, a number that you can tweak) 
# the script runs a set of diagnotics that will be saved in the .Rout file with
# diagonstic graphs saved in modelX.pdf. Look at these to figure out if your
# sampling is satisfactory.
#
# If it is not satisfactory you will need to do extra runs of the script,  
# commening out some sections (as indicated below) and adding extra 
# h.run.dmc commands. An example is provided commented out below. If you see
# that there are problems with some participant fits developing stuck chains 
# then you may want to add the argument replace.bad.chains=TRUE to the 
# h.samples.dmc function (this throws out chains identified as bad by a 
# heuristic and replaces them with slightly perturbed versions of good chains). 
# If you do this then you it is best to do a short run with to clean out the 
# stuck chains followed by a longer run without this  argument to ensure that 
# you are getting respresentative posteriors.
#

rm(list=ls()) 
source ("dmc/dmc.R")
load_model("LBA","lba_B.R")
tmp=load("modelX.RData")
cores=length(sX)

# Fixed effects, comment out on subsequent runs
sX <- h.RUN.dmc(sX,cores=cores)
save(dmX,sX,file="donesmodelX.RData")
# Make hierarchical, comment out on subsequent runs
p.prior <- sX[[1]]$p.prior
p1 <- get.p.vector(sX[[1]])[names(p.prior)]
s.prior <- prior.p.dmc(p1=p1,p2=p1,
  dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)
hstart <- make.hstart(sX)
theta1 <- make.theta1(sX)
hX <- h.samples.dmc(nmc=100,p.prior,dmX,pp.prior,
  hstart.prior=hstart,theta1=theta1,thin=10)
rm(dmX,sX); gc()

# Hierarhcial fit, comment out all but the next line on subsequent runs
pp.prior <- attr(hX,"hyper")$pp.prior
hX  <- h.run.unstuck.dmc(hX, p.migrate = .05, cores = cores)
save(hX,file="modelX.RData")
hX <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hX),
  thorough=TRUE,nmc=50,cores=cores,finalrun=TRUE,finalI=250)
save(hX,file="modelX.RData")
# Uncomment the following and make more copies to run more fits.
# hX <- h.run.dmc(h.samples.dmc(samples=hX,nmc=250),cores=cores,report=10)
# save(hX,file="modelX.RData")


# Make diagnostics, tweak layout as desired
pdf("modelX_chains.pdf",height=6,width = 8)
plot.dmc(hX,hyper=TRUE,pll.chain=TRUE) 
plot.dmc(hX,hyper=TRUE,layout=c(3,3))
plot.dmc(hX,hyper=TRUE,layout=c(3,3),p.prior=pp.prior)
dev.off()

ppX <- h.post.predict.dmc(hX,cores=20)
save(hX,ppX,file="modelX.RData")
pdf("X_fit.pdf",height=6, width = 8)
plot.pp.dmc(ppX,model.legend = FALSE,layout=c(2,2))
dev.off()

h.IC.dmc(hX,DIC=TRUE)


gelman.diag.dmc(hX,hyper=TRUE)
gelman.diag.dmc(hX)

parsX <- summary.dmc(hX,hyper=TRUE)
round(parsX$quantiles[,c(1,3,5)],3)
parsX <- summary.dmc(hX)
save(hX,ppX,parsX,file="modelX.RData")
