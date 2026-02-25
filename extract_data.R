library(hts)
library(bayesRecon)
library(bayesPop)
library(wpp2024)
library(bayesTFR)
library(bayesLife)
library(bayesMig)
library(bayesPop)

subnat.pop.dir <- "bayespop_workshop/results_14/pop"
pop.subnat <- get.pop.prediction(subnat.pop.dir)
pop.aggr <- get.pop.aggregation(sim.dir = subnat.pop.dir) 

location.file <- "bayespop_workshop/data/wafips.txt"
county.codes = read.delim(location.file)$reg_code[-1] #Remove code '53' for Washington State
county.names = read.delim(location.file)$name

library(readxl)
components.file = "data/ofm_1960_2024_components.xlsx"

agg.births = data.frame(read_excel(components.file, col_names = TRUE, sheet = 5, skip = 3)[-40,-(2:3)])
row.names(agg.births) = agg.births[[1]]
agg.births = agg.births[,-1]

agg.deaths = data.frame(read_excel(components.file, sheet = 6, skip = 3)[-40,-(2:3)])
row.names(agg.deaths) = agg.deaths[[1]]
agg.deaths = agg.deaths[,-1]

agg.mig = data.frame(read_excel(components.file, sheet = 8, skip = 3)[-40,-(2:3)])
row.names(agg.mig) = agg.mig[[1]]
agg.mig = agg.mig[,-1]


res_list <- vector("list", length = 0)
for(county.code in county.codes){
  res_list[[length(res_list) + 1]] = apply(get.pop.exba(paste0("B", county.code, "_{10:54}"), pop.subnat, nr.traj = 100), c(1,2), mean)[,2:11] |> t()
  
  res_list[[length(res_list) + 1]] = apply(get.pop.exba(paste0("D", county.code, "_{0:130}"), pop.subnat, nr.traj = 100), c(1,2), mean)[,2:11] |> t()
  
  res_list[[length(res_list) + 1]] = apply(get.pop.exba(paste0("G", county.code, "_{0:100}"), pop.subnat, nr.traj = 100), c(1,2), mean)[,2:11] |> t()
}

col_names = c()
for (county.code in county.codes){
  for (age in 10:54){
    col_names <- c(col_names, paste0("B", county.code, "_",age))
  }
  for (age in 0:130){
    col_names <- c(col_names, paste0("D", county.code, "_",age))
  }
  for (age in 0:100){
    col_names <- c(col_names, paste0("G", county.code, "_",age))
  }
}

df = as.data.frame(do.call(cbind, res_list))
colnames(df) = col_names

S <- rbind(
  c(1, 1, 1, 1),
  c(1, 1, 0, 0),
  c(0, 0, 1, 1),
  c(1, 0, 0, 0),
  c(0, 1, 0, 0),
  c(0, 0, 1, 0),
  c(0, 0, 0, 1)
)

S %*% ginv(t(S)%*%S) %*% t(S)

library(MASS)
S_ginv = ginv(S)

S_ginv %*% S %*% S_ginv


