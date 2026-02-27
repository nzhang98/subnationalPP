library(hts)
library(bayesRecon)
library(bayesPop)
library(wpp2024)
library(bayesTFR)
library(bayesLife)
library(bayesMig)
library(bayesPop)

### Script containing steps to extract relevant datasets to perform Forecast Reconciliation

# WA county code and names
location.file = "bayespop_workshop/data/wafips.txt"
county.codes = read.delim(location.file)$reg_code[-1] #Remove code '53' for Washington State
county.names = read.delim(location.file)$name[-1]

library(readxl)
# Regional aggregate historical and forecasts from the OFM for WA
components.file = "data/ofm_1960_2024_components.xlsx"

agg.births = data.frame(read_excel(components.file, col_names = TRUE, sheet = 5, skip = 3)[-40,-(2:3)])
row.names(agg.births) = agg.births[[1]]
agg.births = agg.births[,-1]
colnames(agg.births) = seq(1961,2025)

agg.deaths = data.frame(read_excel(components.file, sheet = 6, skip = 3)[-40,-(2:3)])
row.names(agg.deaths) = agg.deaths[[1]]
agg.deaths = agg.deaths[,-1]
colnames(agg.births) = seq(1961,2025)

agg.mig = data.frame(read_excel(components.file, sheet = 8, skip = 3)[-40,-(2:3)])
row.names(agg.mig) = agg.mig[[1]]
agg.mig = agg.mig[,-1]
colnames(agg.mig) = seq(1961,2025)

agg.pop = data.frame(read_excel(components.file, sheet = 3, skip = 3)[-40,-(2:3)])
row.names(agg.pop) = agg.pop[[1]]
agg.pop = agg.pop[,-1]
colnames(agg.pop) = seq(1961,2025)

### OFM county aggregate forecasts from 2007, used to compute in-sample forecast errors
agg.pop.fcast07 = data.frame(read_excel("data/ofm_2007_pop_county_forecasts_medium.xls", col_names = TRUE,
                                        sheet = 1, skip = 3)[-(41:43),]) 
row.names(agg.pop.fcast07) = agg.pop.fcast07[[1]]
agg.pop.fcast07 = agg.pop.fcast07[,-1]
colnames(agg.pop.fcast07) = c(2000, 2005, seq(2010, 2030))

### OFM county aggregate forecasts from 2012, used as aggregate forecasts from 2016 to reconcile with
agg.pop.fcast12 = data.frame(read_excel("data/ofm_2012_pop_county_forecasts_medium.xls", col_names = TRUE,
                                        sheet = 1, skip = 4)[-(41:46),]) 
row.names(agg.pop.fcast12) = agg.pop.fcast12[[1]]
agg.pop.fcast12 = agg.pop.fcast12[,-1]
colnames(agg.pop.fcast12) = c(2010, seq(2015, 2040))

###################

#### Retrieve Observations/Historical data from bayesPop
subnat.pop.dir = "bayespop_workshop/results_24/pop"
pop.subnat = get.pop.prediction(subnat.pop.dir)
pop.aggr = get.pop.aggregation(sim.dir = subnat.pop.dir) 

res_list = vector("list", length = 0)
for(county.code in county.codes){
  res_list[[length(res_list) + 1]] = get.pop.exba(paste0("B", county.code), observed = TRUE, pop.subnat)[, as.character(("2014":"2023"))]
  
  res_list[[length(res_list) + 1]] = get.pop.exba(paste0("D", county.code), observed = TRUE, pop.subnat)[, as.character(("2014":"2023"))]
  
  res_list[[length(res_list) + 1]] = get.pop.exba(paste0("G", county.code), observed = TRUE, pop.subnat)[, as.character(("2014":"2023"))]
}

b_obs = do.call(cbind, res_list)
colnames(b_obs) = col_names

#### Retrieve forecasts from bayesPop, made at year 2013
subnat.pop.dir = "bayespop_workshop/results_14/pop"
pop.subnat = get.pop.prediction(subnat.pop.dir)
pop.aggr = get.pop.aggregation(sim.dir = subnat.pop.dir) 

res_list = vector("list", length = 0)
for(county.code in county.codes){
  res_list[[length(res_list) + 1]] = rowMeans(get.pop.exba(paste0("B", county.code), pop.subnat, nr.traj = 100))[2:11]
  
  res_list[[length(res_list) + 1]] = rowMeans(get.pop.exba(paste0("D", county.code), pop.subnat, nr.traj = 100))[2:11]
  
  res_list[[length(res_list) + 1]] = rowMeans(get.pop.exba(paste0("G", county.code), pop.subnat, nr.traj = 100))[2:11]
}

b_base = do.call(cbind, res_list)

col_names = c()
for (county.name in county.names){
    col_names = c(col_names, paste0("B_", county.name))
    col_names = c(col_names, paste0("D_", county.name))
    col_names = c(col_names, paste0("G_", county.name))
}

colnames(b_base) = col_names



################### Retrieve data disaggregate into age groups

res_list = vector("list", length = 0)
for(county.code in county.codes){
  res_list[[length(res_list) + 1]] = apply(get.pop.exba(paste0("B", county.code, "_{10:54}"), pop.subnat, nr.traj = 100), c(1,2), mean)[,2:11] |> t()
  
  res_list[[length(res_list) + 1]] = apply(get.pop.exba(paste0("D", county.code, "_{0:130}"), pop.subnat, nr.traj = 100), c(1,2), mean)[,2:11] |> t()
  
  res_list[[length(res_list) + 1]] = apply(get.pop.exba(paste0("G", county.code, "_{0:100}"), pop.subnat, nr.traj = 100), c(1,2), mean)[,2:11] |> t()
}

col_names = c()
for (county.code in county.codes){
  for (age in 10:54){
    col_names = c(col_names, paste0("B_", county.code, "_",age))
  }
  for (age in 0:130){
    col_names = c(col_names, paste0("D_", county.code, "_",age))
  }
  for (age in 0:100){
    col_names = c(col_names, paste0("G_", county.code, "_",age))
  }
}

df = as.data.frame(do.call(cbind, res_list))
colnames(df) = col_names


