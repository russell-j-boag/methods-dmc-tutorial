### Script for making simulated LBA data ###

# Clear workspace
rm(list = ls())

# Set working directory to the top-level folder containing the dmc and tutorial subfolders 
getwd()

# Load packages and model functions
source("dmc/dmc.R")
load_model(dir_name = "LBA", model_name = "lba_B.R")

# Build model
model <- model.dmc(p.map     = list(A = "1", 
                                    B = "I", 
                                    t0 = "1", 
                                    mean_v = c("D", "M"), 
                                    sd_v = "M", 
                                    st0="1"), 
                   match.map = list(M = list(left = "L", right = "R")),
                   factors = list(S = c("left", "right"), 
                                  D = c("easy", "hard"), 
                                  I = c("acc", "spd")),
                   constants = c(sd_v.false = 1, st0 = 0), 
                   responses = c("L", "R"),
                   type = "norm"
)

# Set parameter mean values
p.vector <- c(A = .4,
              B.acc = .8, B.spd = .6,
              t0 = .3,
              mean_v.easy.true = 1.5, mean_v.hard.true = 1,
              mean_v.easy.false = 0, mean_v.hard.false = 0,
              sd_v.true = .25
)

# Show mapping from parameters to design cells
print.cell.p(p.vector, model)

# Creat priors object
p.prior <- prior.p.dmc(
  dists = rep("tnorm", 9),
  p1 = p.vector,
  p2 = c(rep(.2, 9)),
  lower = c(0, rep(0,2), .1, rep(NA,4), 0),
  upper = c(NA, rep(NA,2), 1, rep(NA,4), 1)
)
length(p.prior)

# Simulate data (single subject)
dat <- simulate.dmc(p.vector = p.vector, model = model, n = 80)
head(dat)

# Simulate data (multiple subjects)
hdat <- h.simulate.dmc(model = model, ps = p.vector, ns = 30, n = 80)
head(hdat)

# Summarise accuracy
dat %>% 
  group_by(I, D) %>%      # factors to group by
  summarise(accuracy = mean(C))  # compute accuracy

# Summarise mean RT, standard deviation, and cell counts
dat %>% 
  group_by(I, D, C) %>%    # factors to group by
  summarise(mean.RT = mean(RT),  # compute mean RT
            sd.RT = sd(RT),      # compute std. dev.
            n.obs = n())         # get cell counts

# Save data
write.csv(dat, file="dat.csv", row.names=FALSE)
write.csv(hdat, file="hdat.csv", row.names=FALSE)


