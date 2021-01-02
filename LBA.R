### Script for fitting the LBA model to choice-RT data ###

# Clear workspace
rm(list = ls())

# Set working directory to the top-level folder containing the dmc and tutorial subfolders 
getwd()  # shows current working directory
setwd("/Users/rboag/Documents/Course material - Methods/dmc")

# Load packages and model functions
source("dmc/dmc.R")
load_model(dir_name = "LBA", model_name = "lba_B.R")

# Load data
dat <- read.csv("dat.csv", header=TRUE)

# Show the first few rows of the dataframe
head(dat)

# Show how the data are structured 
str(dat)
# So we have several factors each with two levels:
# "S" = stimulus type (left/right)
# "D" = stimulus difficulty (easy/hard)
# "I" = speed/accuracy instruction (acc/spd)
# "R" = response (L/R)
# The numeric "RT" column contains the response time on each trial in seconds

# Set categorical variables to "factor" type and RT variable to "numeric" type
dat$S <- as.factor(dat$S)
dat$D <- as.factor(dat$D)
dat$I <- as.factor(dat$I)
dat$R <- as.factor(dat$R)
dat$RT <- as.numeric(dat$RT)
str(dat)

# Add factor to code accuracy (i.e., correct = 1, error = 0)
dat$C <- NA  # create empty column
dat$C[(dat$S == "left" & dat$R == "L")|(dat$S == "right" & dat$R == "R")] <- 1  # code correct responses "1"
dat$C[(dat$S == "right" & dat$R == "L")|(dat$S == "left" & dat$R == "R")] <- 0  # code error responses "0"
dat$C <- as.numeric(dat$C)  # set as "numeric" type
head(dat)
str(dat)

# Now we can look at accuracy by cell
dat %>% 
  group_by(I, D) %>%             # factors to group by
  summarise(accuracy = mean(C))  # compute accuracy

# We can also look at mean RT, standard deviation, and cell counts
dat %>% 
  group_by(I, D, C) %>%          # factors to group by
  summarise(mean.RT = mean(RT),  # compute mean RT
            sd.RT = sd(RT),      # compute std. dev.
            n.obs = n())         # get cell counts

# Plot RT histograms
par(mfcol = c(2,2))
hist(dat$RT[(dat$D == "easy" & dat$C == 1) & dat$I == "spd"], 
     breaks = 20, xlab = "RT (s)", 
     main = "Easy, Spd", xlim = c(.2,1.5))  # easy, speed condition
hist(dat$RT[(dat$D == "easy" & dat$C == 1) & dat$I == "acc"], 
     breaks = 20, xlab = "RT (s)", 
     main = "Easy, Acc", xlim = c(.2,1.5))  # easy, accuracy condition
hist(dat$RT[(dat$D == "hard" & dat$C == 1) & dat$I == "spd"], 
     breaks = 20, xlab = "RT (s)", 
     main = "Hard, Spd", xlim = c(.2,2.2))  # hard, speed condition
hist(dat$RT[(dat$D == "hard" & dat$C == 1) & dat$I == "acc"], 
     breaks = 20, xlab = "RT (s)", 
     main = "Hard, Acc", xlim = c(.2,2.2))  # hard, accuracy condition

# We will now build an LBA and fit it to these data


### Model 1: LBA with drift rate effect

# We will start with a model that allows drift rate to vary by stimulus 
# difficulty/discriminability (i.e., low versus high motion coherence)
# We would expect perceptually difficult stimuli (low coherence) to have 
# lower drift rates than perceptually easy stimuli (high coherence)

## Step 1 - Build model and set prior parameter values
# Build model
model <- model.dmc(p.map     = list(A = "1",
                                    B = "1",
                                    t0 = "1",
                                    mean_v = c("D", "M"),
                                    sd_v = "M",
                                    st0 = "1"), 
                   match.map = list(M = list(left = "L", right = "R")),
                   factors = list(S = c("left", "right"),
                                  D = c("easy", "hard"),
                                  I = c("acc", "spd")),
                   constants = c(sd_v.false = 1, st0 = 0), 
                   responses = c("L", "R"), 
                   type = "norm"
)

# Set parameter mean values
p.vector = c(A = .3, B = .7, t0 = .2,
             mean_v.easy.true = 1.2, mean_v.hard.true = .9,
             mean_v.easy.false = 0, mean_v.hard.false = 0,
             sd_v.true = .2
)
length(p.vector)

# Show mapping from parameters to design cells
print.cell.p(p.vector, model)

# Creat priors object
p.prior <- prior.p.dmc(
  dists = rep("tnorm", 8),
  p1 = p.vector,
  p2 = c(rep(1, 8)),
  lower = c(0, rep(0, 1), .1, rep(NA, 4), 0),
  upper = c(NA, rep(NA, 1), 1, rep(NA, 4), 1)
)

# Plot prior distributions
par(mfcol = c(2, 4)); for (i in names(p.prior)) plot.prior(i, p.prior)

# Combine model and data into a data.model object
data.model <- data.model.dmc(dat, model)

# See the likelihoods of the model given the data 
likelihood.dmc(p.vector, data.model)  # the sampler will use these to find the best parameters


## Step 2 - Run sampling and check sampling diagnostics
# Create samples object to pass to the sampler
samples <- samples.dmc(nmc = 100, p.prior = p.prior, data = data.model, thin = 5, n.chains = 24)

# Run sampling
samples <- run.dmc(samples, cores = 4, report = 5, p.migrate = 0.05)

# Check whether chains have converged
gelman.diag.dmc(samples)  # is Gelman's R < 1.1?

# Run automatic sampling (this will continue sampling until several diagnostic criteria are met)
samples1 <- RUN.dmc(samples.dmc(nmc = 100, samples = samples, thin = 5, n.chains = 24),
                    cores = 4, report = 5)  # this will take a few minutes

# Chains should have converged by now
gelman.diag.dmc(samples1)  # is Gelman's R < 1.1?

# Save samples (so that you don't have to re-run the sampling)
save(samples, samples1, file = "samples_mod1.RData")
load(file = "samples_mod1.RData")  # load saved samples

# Plot sampling chains to check for convergence (they should look like 'big fuzzy caterpillars')
plot.dmc(samples1)  # chains by parameter
par(mfcol = c(1, 1)); plot.dmc(samples1, pll.chain = TRUE, start = 1)  # overall chain log-likelihood


## Step 3 - Assess and plot model fit 
# Sample posterior predictives
posts <- post.predict.dmc(samples1, n.post = 100)

# Plot model fit
plot.pp.dmc(posts, "pdf", layout = c(2, 2), model.legend = FALSE)  # prob. density functions
plot.pp.dmc(posts, "cdf", layout = c(2, 2), model.legend = FALSE)  # cumulative prob. density functions

# Show deviance stats (goodness of fit between model and data)
Dstats.dmc(samples1)  # smaller is better!


## Step 4 - Check parameter recovery and plot parameter effects
# Check parameter recovery
check.recovery.dmc(samples1, digits = 3)  # shows the uncertainty surrounding parameter estimates

# Plot correlations between parameters
pairs.dmc(samples1, start = 50)  # correlated parameters can be a problem

# Get parameter summary statistics (means & SDs) for plotting
parms <- summary.dmc(samples1)
parms$statistics

# Before we can plot, we need to make a dataframe containing only the parameters
# that we want to plot
# First we grab the mean and SDs for drift rates from the summary.dmc object and 
# put it into a data.frame called 'drift'
drift <- data.frame(parms$statistics[grep("mean_v", rownames(parms$statistics)), c("Mean","SD")])

# Then we create a factor for the effect we are interested in
drift$D <- NA; drift$C <- NA
drift$D[grep("easy", rownames(drift))] <- "Easy"
drift$D[grep("hard", rownames(drift))] <- "Hard"

# We will also create a factor so that we can plot correct and error drift rates separately
drift$C[grep("true", rownames(drift))] <- "Correct"
drift$C[grep("false", rownames(drift))] <- "Error"
drift$D <- factor(drift$D); drift$C <- factor(drift$C)
drift

# Now we can plot drift rates by stimulus discriminability
ggplot(drift, aes(x = D, y = Mean, shape = C)) +
  geom_point(stat = "identity", aes(), size = 3) +
  geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD, width = 0.2)) +
  geom_line(aes(group = C, y = Mean), linetype = 2) + 
  ylim(-.8, 2.4) +
  xlab("Stimulus difficulty (discriminability)") + 
  ylab("Drift rate") +
  theme_classic()

# So we see that drift rates toward the correct response are lower for perceptually
# difficult stimuli
# However, this model didn't seem to fit the data very well, so it is likely too simple
# We will now build a model that lets thresholds vary by speed/accuracy 
# condition and see if that captures the data better


### Model 2: LBA with threshold effect

# In this model we will let thresholds vary by speed/accuracy instructions (acc/spd)
# We would expect people to lower their thresholds (i.e., gather less evidence 
# before responding) under speed instructions and raise their thresholds (i.e., 
# gather more evidence before responding) under accuracy instructions

## Step 1 - Build model and set prior parameter values
# Build model
model <- model.dmc(p.map     = list(A = "1", 
                                    B = "I",
                                    t0 = "1",
                                    mean_v = c("M"),
                                    sd_v = "M",
                                    st0 = "1"), 
                   match.map = list(M = list(left = "L", right = "R")),
                   factors = list(S = c("left", "right"),
                                  D = c("easy", "hard"),
                                  I = c("acc", "spd")),
                   constants = c(sd_v.false = 1, st0 = 0), 
                   responses = c("L", "R"),
                   type="norm"
)

# Set parameter mean values
p.vector <- c(A = .3, 
              B.acc = 1.3, B.spd = .9, 
              t0 = .2,
              mean_v.true = 1.2, mean_v.false = 0,
              sd_v.true = .2
)

# Show mapping from parameters to design cells
print.cell.p(p.vector, model)

# Creat priors object
p.prior <- prior.p.dmc(
  dists = rep("tnorm", 7),
  p1 = p.vector,
  p2 = c(rep(1, 7)),
  lower = c(0, rep(0, 2), .1, rep(NA, 2), 0),
  upper = c(NA, rep(NA, 2), 1, rep(NA, 2), 1)
)

# Plot prior distributions
par(mfcol = c(2, 4)); for (i in names(p.prior)) plot.prior(i, p.prior)

# Combine model and data into a data.model object
data.model <- data.model.dmc(dat, model)

# See the likelihoods of the model given the data 
likelihood.dmc(p.vector, data.model)  # the sampler will use these to find the best parameters


## Step 2 - Run sampling and check sampling diagnostics
# Create samples object to pass to the sampler
samples <- samples.dmc(nmc = 100, p.prior = p.prior, data = data.model, thin = 5, n.chains = 21)

# Run sampling
samples <- run.dmc(samples, cores = 4, report = 5, p.migrate = 0.05)

# Check whether chains have converged
gelman.diag.dmc(samples)  # is Gelman's R < 1.1?

# Run automatic sampling (this will continue sampling until several diagnostic criteria are met)
samples1 <- RUN.dmc(samples.dmc(nmc = 100, samples = samples, thin = 5, n.chains = 21),
                    cores = 4, report = 5)  # this will take a few minutes

# Chains should have converged by now
gelman.diag.dmc(samples1)  # is Gelman's R < 1.1?

# Save samples (so that you don't have to re-run the sampling)
save(samples, samples1, file = "samples_mod2.RData")
load(file = "samples_mod2.RData")  # load saved samples

# Plot sampling chains to check for convergence (they should look like 'big fuzzy caterpillars')
plot.dmc(samples1)  # chains by parameter
par(mfcol = c(1, 1)); plot.dmc(samples1, pll.chain = TRUE, start = 1)  # overall chain log-likelihood


## Step 3 - Assess and plot model fit 
# Sample posterior predictives
posts <- post.predict.dmc(samples1, n.post = 100)

# Plot model fit
plot.pp.dmc(posts, "pdf", layout = c(2, 2), model.legend = FALSE)  # prob. density functions
plot.pp.dmc(posts, "cdf", layout = c(2, 2), model.legend = FALSE)  # cumulative prob. density functions

# Show deviance stats (goodness of fit between model and data)
Dstats.dmc(samples1)  # smaller is better!


## Step 4 - Check parameter recovery and plot parameter effects
# Check parameter recovery
check.recovery.dmc(samples1, digits = 3)  # shows the uncertainty surrounding parameter estimates

# Plot correlations between parameters
pairs.dmc(samples1, start = 50)  # correlated parameters can be a problem

# Get parameter summary statistics (means & SDs) for plotting
parms <- summary.dmc(samples1)
parms$statistics

# Before we can plot, we need to make a dataframe containing only the parameters
# that we want to plot
# First we grab the mean and SDs for thresholds from the summary.dmc object and 
# put them into a dataframe called 'thresholds'
threshold <- data.frame(parms$statistics[grep("B", rownames(parms$statistics)), c("Mean","SD")])

# Then we create a factor for the effect we are interested in 
threshold$I <- NA
threshold$I[grep("acc", rownames(threshold))] <- "Accuracy"
threshold$I[grep("spd", rownames(threshold))] <- "Speed"
threshold$I <- factor(threshold$I)
threshold

# Now we can plot thresholds by speed/accuracy instructions
ggplot(threshold, aes(x = I, y = Mean)) +
  geom_point(stat = "identity", aes(), size = 3) +
  geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD, width = 0.2)) +
  geom_line(aes(group = 1, y = Mean), linetype = 2) + 
  ylim(0.3, 2) +
  xlab("Speed-accuracy instruction") + 
  ylab("Threshold") +
  theme_classic()

# So we see that thresholds tend to be lower under speed instructions (as expected)
# However this model also didn't fit very well, so we probably need to also include 
# a drift rate effect 
# We will now build a model with both drift and threshold effects and see if that 
# improves the fit


### Model 3: LBA with drift rate and threshold effects

# Now we will build a model that lets both drift rate and thresholds vary
# We will let drift rate vary by stimulus difficulty (as in the first model) 
# and let thresholds vary by speed/accuracy instructions

## Step 1 - Build model and set prior parameter values
# Build model
model <- model.dmc(p.map     = list(A = "1",
                                    B = "I",
                                    t0 = "1",
                                    mean_v = c("D", "M"),
                                    sd_v = "M",
                                    st0 = "1"), 
                   match.map = list(M = list(left = "L", right = "R")),
                   factors = list(S = c("left", "right"),
                                  D = c("easy", "hard"),
                                  I = c("acc", "spd")),
                   constants = c(sd_v.false = 1, st0 = 0), 
                   responses = c("L", "R"),
                   type = "norm")

# Set parameter mean values
p.vector <- c(A = .3,
              B.acc = 1.3, B.spd = .9,
              t0 = .2,
              mean_v.easy.true = 1.4, mean_v.hard.true = .9,
              mean_v.easy.false = 0, mean_v.hard.false = 0,
              sd_v.true = .2
)

# Show mapping from parameters to design cells
print.cell.p(p.vector, model)

# Creat priors object
p.prior <- prior.p.dmc(
  dists = rep("tnorm", 9),
  p1 = p.vector,
  p2 = c(rep(1, 9)),
  lower = c(0, rep(0, 2), .1, rep(NA, 4), 0),
  upper = c(NA, rep(NA, 2), 1, rep(NA, 4), 1)
)

# Plot prior distributions
par(mfcol = c(3, 3)); for (i in names(p.prior)) plot.prior(i, p.prior)

# Combine model and data into a data.model object
data.model <- data.model.dmc(dat, model)

# See the likelihoods of the model given the data 
likelihood.dmc(p.vector, data.model)  # the sampler will use these to find the best parameters


## Step 2 - Run sampling and check sampling diagnostics
# Create samples object to pass to the sampler
samples <- samples.dmc(nmc = 100, p.prior = p.prior, data = data.model, thin = 5, n.chains = 27)

# Run sampling
samples <- run.dmc(samples, cores = 4, report = 5, p.migrate = 0.05)

# Check whether chains have converged
gelman.diag.dmc(samples)  # is Gelman's R < 1.1?

# Run automatic sampling (this will continue sampling until several diagnostic criteria are met)
samples1 <- RUN.dmc(samples.dmc(nmc = 100, samples = samples, thin = 5, n.chains = 27),
                    cores = 4, report = 5)  # this will take a few minutes

# Chains should have converged by now
gelman.diag.dmc(samples1)  # is Gelman's R < 1.1?

# Save samples (so that you don't have to re-run the sampling)
save(samples, samples1, file = "samples_mod3.RData")
load(file = "samples_mod3.RData")  # load saved samples

# Plot sampling chains to check for convergence (they should look like 'big fuzzy caterpillars')
plot.dmc(samples1)  # chains by parameter
par(mfcol = c(1, 1)); plot.dmc(samples1, pll.chain = TRUE, start = 1)  # overall chain log-likelihood


## Step 3 - Assess and plot model fit 
# Sample posterior predictives
posts <- post.predict.dmc(samples1, n.post = 100)

# Plot model fit
plot.pp.dmc(posts, "pdf", layout = c(2, 2), model.legend = FALSE)  # prob. density functions
plot.pp.dmc(posts,"cdf", layout = c(2, 2), model.legend = FALSE)  # cumulative prob. density functions

# Show deviance stats (goodness of fit between model and data)
Dstats.dmc(samples1)  # smaller is better!


## Step 4 - Check parameter recovery and plot parameter effects
# Check parameter recovery
check.recovery.dmc(samples1, digits = 3)  # shows the uncertainty surrounding parameter estimates

# Plot correlations between parameters
pairs.dmc(samples1, start = 50)  # correlated parameters can be a problem

# Get parameter summary statistics (means & SDs) for plotting
parms <- summary.dmc(samples1)
parms$statistics

# Before we can plot, we need to make a dataframe containing only the parameters
# that we want to plot
# First we grab the mean and SDs for thresholds from the summary.dmc object and 
# put them into a dataframe called 'thresholds'
threshold <- data.frame(parms$statistics[grep("B", rownames(parms$statistics)), c("Mean", "SD")])

# Then we create a factor for the effect we are interested in 
threshold$I <- NA
threshold$I[grep("acc", rownames(threshold))] <- "Accuracy"
threshold$I[grep("spd", rownames(threshold))] <- "Speed"
threshold$I <- factor(threshold$I)
threshold

# Now we can plot thresholds by speed/accuracy instructions
ggplot(threshold, aes(x = I, y = Mean)) +
  geom_point(stat = "identity", aes(), size = 3) +
  geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD, width = 0.2)) +
  geom_line(aes(group = 1, y = Mean), linetype = 2) + 
  ylim(0.3, 1) +
  xlab("Speed-accuracy instruction") + 
  ylab("Threshold") +
  theme_classic()

# Now to plot drift rates, we will grab the mean and SDs for drift rates 
# from the summary.dmc object and put it into a data.frame called 'drift'
drift <- data.frame(parms$statistics[grep("mean_v", rownames(parms$statistics)), c("Mean", "SD")])

# Then we create a factor for the effect we are interested in
drift$D <- NA; drift$C <- NA
drift$D[grep("easy", rownames(drift))] <- "Easy"
drift$D[grep("hard", rownames(drift))] <- "Hard"

# We will also create a factor so that we can plot correct and error drift rates separately
drift$C[grep("true", rownames(drift))] <- "Correct"
drift$C[grep("false", rownames(drift))] <- "Error"
drift$D <- factor(drift$D); drift$C <- factor(drift$C)
drift

# Plot drift rates by stimulus discriminability
ggplot(drift, aes(x = D, y = Mean, shape = C)) +
  geom_point(stat = "identity", aes(), size = 3) +
  geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD, width = 0.2)) +
  geom_line(aes(group = C, y = Mean), linetype = 2) + 
  ylim(-1, 2) +
  xlab("Stimulus difficulty (discriminability)") + 
  ylab("Drift rate") +
  theme_classic()

# So we see both effects present, i.e., lower drift rates for perceptually difficult 
# stimuli and lower thresholds under speed instructions
# This model also fits much better than either the drift-only or threshold-only models 
# (and we should expect it to since this model was used to generate the data!)

