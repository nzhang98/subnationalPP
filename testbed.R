library(renv)
renv::init()
# renv::snapshot(type = "implicit")
renv::settings$snapshot.type("explicit")
renv::status()
renv::dependencies()

######################

install.packages(c("bayesTFR", "bayesLife", "bayesPop", "MortCast", 
                   "bayesMig", "devtools"), dependencies = TRUE)

library(devtools)
# options(timeout = 600)
# install_github("PPgp/wpp2024")

library(bayesPop)
library(wpp2024)
library(bayesTFR)
library(bayesLife)
library(bayesMig)
library(bayesPop)
# sessionInfo()

getwd()
setwd('C:/Users/nakaz/OneDrive/Desktop/Research/PopulationProjection/SubnationalPP/bayespop_workshop')

my.tfr.file <- "data/tfr.txt"
head(read.delim(my.tfr.file, check.names = FALSE))


nat.tfr.dir <- "data/wpp2024_projections/TFR1unc/sim20241101"

subnat.tfr.dir <- "results/tfr"

start_year = 2014
tfr.preds <- tfr.predict.subnat(840, my.tfr.file = my.tfr.file,
                                sim.dir = nat.tfr.dir, output.dir = subnat.tfr.dir,
                                annual = TRUE, start.year = start_year, end.year = 2100,
                                nr.traj = 1000)


tfr.preds <- get.regtfr.prediction(subnat.tfr.dir)

names(tfr.preds)
tfr.subnat <- tfr.preds[["840"]]
# identical to
# tfr.subnat <- get.tfr.prediction(file.path(subnat.tfr.dir, "subnat/c840"))
# or
# tfr.subnat <- get.regtfr.prediction(subnat.tfr.dir, 840)
get.countries.table(tfr.subnat)

county <- "King"
tfr.trajectories.plot(tfr.subnat, county, nr.traj = 50)
tfr.trajectories.table(tfr.subnat, county)
summary(tfr.subnat, county)

nat.tfr.pred <- get.tfr.prediction(nat.tfr.dir)
tfr.trajectories.plot(nat.tfr.pred, 840, nr.traj = 50, 
                      uncertainty = TRUE)

county <- "Adams"
tfr.trajectories.plot(tfr.subnat, county, pi = 80, 
                      half.child.variant = FALSE, nr.traj = 50)
tfr.trajectories.plot(nat.tfr.pred, 840, pi = 80,
                      half.child.variant = FALSE, nr.traj = 0, 
                      add = TRUE, col = rep("darkblue", 5), 
                      show.legend = FALSE)

trajs <- get.tfr.trajectories(tfr.subnat, county)
dim(trajs)
summary(t(trajs))

my.e0M.file <- "data/e0M.txt"  # male data
my.e0F.file <- "data/e0F.txt"  # female data
head(read.delim(my.e0F.file, check.names = FALSE))

nat.e0.dir <- "data/wpp2024_projections/e01/sim20241101"


subnat.e0.dir <- "results/e0"

e0.preds <- e0.predict.subnat(840, my.e0.file = my.e0F.file, 
                              sim.dir = nat.e0.dir, output.dir = subnat.e0.dir,
                              annual = TRUE, start.year = start_year, end.year = 2100,
                              predict.jmale = TRUE, my.e0M.file = my.e0M.file,
                              nr.traj = 1000)

e0.subnat <- get.rege0.prediction(subnat.e0.dir, 840)
e0.trajectories.plot(e0.subnat, "Spokane", both.sexes = TRUE)

nat.e0.pred <- get.e0.prediction(nat.e0.dir)
e0.trajectories.plot(e0.subnat, "King", pi = 80, nr.traj = 50)
e0.trajectories.plot(nat.e0.pred, 840, pi = 80, add = TRUE, 
                     col = rep("darkblue", 5), nr.traj = 0, 
                     show.legend = FALSE)

e0.trajectories.table(e0.subnat, "King", both.sexes = TRUE, 
                      pi = c(60, 80, 90))

e0.joint.plot(e0.subnat, "Pierce", years = c(2023, 2030, 2040, 2050))

e0M.subnat <- get.e0.jmale.prediction(e0.subnat)
trajF <- get.e0.trajectories(e0.subnat, "Clark")
trajM <- get.e0.trajectories(e0M.subnat, "Clark")

my.mig.file <- "data/mig_rates.txt"
head(read.delim(my.mig.file, check.names = FALSE))
exclude.codes <- c(53001, 53003, 53013, 53019, 53023, 53039, 
                   53043, 53049, 53051, 53055, 53059, 53069, 53)


subnat.mig.dir <- "results/mig"
mc <- run.mig.mcmc(nr.chains = 3, iter = 5000, 
                   output.dir = subnat.mig.dir,
                   present.year = 2013, start.year = 1990, 
                   my.mig.file = my.mig.file,  annual = TRUE,
                   exclude.from.world = exclude.codes,
                   replace.output = TRUE, verbose.iter = 1000)

mig.subnat <- mig.predict(sim.dir = subnat.mig.dir, nr.traj = 100,
                          burnin = 1000, end.year = 2050, save.as.ascii = 100)

mig.subnat <- get.mig.prediction(sim.dir = subnat.mig.dir)
mig.trajectories.plot(mig.subnat, "Walla Walla", nr.traj = 30)

location.file <- "data/wafips.txt"
head(read.delim(location.file))

popM0.file <- "data/popM.txt"
popF0.file <- "data/popF.txt"
head(read.delim(popM0.file, check.names = FALSE))
popF_df = read.delim(popM0.file, check.names = FALSE)

mxM.file <- "data/mxM.txt"
mxF.file <- "data/mxF.txt"
head(read.delim(mxM.file, check.names = FALSE))

mig.file <- "data/mig_counts.txt"
head(read.delim(mig.file, check.names = FALSE))

gqm.file <- "data/gqM.txt"
gqf.file <- "data/gqF.txt"
head(read.delim(gqm.file))

mig.traj.file <- file.path(subnat.mig.dir, 
                           "predictions/ascii_trajectories.csv")

subnat.tfr.results <- file.path(subnat.tfr.dir, "subnat/c840")
subnat.e0.results <- file.path(subnat.e0.dir, "subnat_ar1/c840")


subnat.pop.dir <- "results/pop"

pop.subnat <- pop.predict.subnat(output.dir = subnat.pop.dir,
                                 locations = location.file, default.country = 840, 
                                 verbose = TRUE, annual = TRUE, wpp.year = 2024, 
                                 present.year = 2013, end.year = 2050,
                                 nr.traj = 100, replace.output = TRUE,
                                 inputs = list(
                                   popM = popM0.file, popF = popF0.file,
                                   mxM = mxM.file, mxF = mxF.file,
                                   mig = mig.file, migtraj = mig.traj.file,
                                   tfr.sim.dir = subnat.tfr.results,
                                   e0F.sim.dir = subnat.e0.results,
                                   e0M.sim.dir = "joint_",
                                   GQpopM = gqm.file, GQpopF = gqf.file
                                 ),
                                 mig.age.method = "rc", 
                                 mig.is.rate = c(FALSE, TRUE), 
                                 keep.vital.events = TRUE, 
                                 pasfr.ignore.phase2 = TRUE
)


pop.aggr <- pop.aggregate.subnat(pop.subnat, regions = 53, 
                                 locations = location.file)

pop.subnat <- get.pop.prediction(subnat.pop.dir)
pop.aggr <- get.pop.aggregation(sim.dir = subnat.pop.dir) 

get.countries.table(pop.subnat)
get.countries.table(pop.aggr)

pop.trajectories.plot(pop.subnat, "King", nr.traj = 20)
pop.trajectories.plot(pop.aggr, 53, nr.traj = 20)

pop.trajectories.table(pop.subnat, "King")

pop.byage.plot(pop.subnat, 53061, year = 2050, nr.traj = 20)

pop.pyramid(pop.aggr, 53, year = 2050)

pop.pyramid(pop.aggr, 53, year = c(2050, 2020), proportion = TRUE)

pop.trajectories.plot(pop.subnat, 
                      expression = "P53033[20:64] / P53033[65:130]", 
                      nr.traj = 20,
                      main = "Potential Support Ratio for King")


pop.byage.plot(pop.subnat, expression = "log(M53031_M{})", nr.traj = 20, 
               year = 2050,
               main = "2050 Male mortality rates in Jefferson")

# extract birth trajectories for Snohomish + Pierce
trajSP <- get.pop.ex("B53061 + B53053", pop.subnat) 
# extract birth trajectories for King
trajK <- get.pop.ex("B53033", pop.subnat)

dim(trajSP)
# frequency of the event happening in 2050
sum(trajSP["2050", ] > trajK["2050",])/ncol(trajSP) * 100
# or equivalently, extracting the comparison as an expression
sum(get.pop.ex("(B53061 + B53053) > B53033", pop.subnat)["2050", ]
) / pop.subnat$nr.traj * 100


pop.trajectories.plot(pop.subnat, expression = "B53033", 
                      nr.traj = 20, pi = 80,
                      main = "Number of births in King vs. Snohomish + Pierce")
pop.trajectories.plot(pop.subnat, expression = "B53061 + B53053", 
                      nr.traj = 0, pi = 80, show.legend = FALSE,
                      add = TRUE, col = rep("blue",4))

sum(get.pop.ex("P53033[20:64] / P53033[65:130] < 3", pop.subnat)["2050", ]
) / pop.subnat$nr.traj * 100


get.pop.ex("B53033[30]", pop.subnat)
get.pop.ex("D53033_M[131]", pop.subnat)
get.pop.exba("G53033_F{10:20}", pop.subnat)
get.pop.ex("G53033_F[10]", pop.subnat)
