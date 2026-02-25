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





























