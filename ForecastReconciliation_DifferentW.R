library(quadprog) # to solve QP
library(MASS) # for generalised inverse

# Script to perform Forecast Reconciliation: Example 1
# Forecast Reconciliation with different reconciliation solutions, simplified example with no age-sex groupings
# Datasets for this script are obtained with "extract_data.R"

# Aggregates
# 
# N_rc = t(agg.pop[, as.character(2014:2023)])
# colnames(N_rc) = paste0("dN_", colnames(N_rc))
# N_rc = cbind(dN_tot = rowSums(N_rc), N_rc)
# 
# Y_base = cbind(N_rc, b_base)[as.character(2016:2023),]
# colnames(Y_base) = colnames(Y_obs)
################################################################
# Fcast error (in-sample)

diff(as.numeric(agg.pop.fcast07))

fcast07 = t((agg.pop.fcast07[, -1] - agg.pop.fcast07[, -ncol(agg.pop.fcast07)])[, c(6:10)])
obs = t(rbind(colSums(agg.pop), agg.pop)[, as.character(2014:2018)])

Y_base_insample = cbind(fcast07, b_base[as.character(2014:2018),])
Y_obs = cbind(obs, b_obs[as.character(2014:2018),])

fcast_error = Y_base_insample - Y_obs
#################################################################

fcast12 = t((agg.pop.fcast12[, -1] - agg.pop.fcast12[, -ncol(agg.pop.fcast12)])[, as.character(2016:2023)])
colnames(fcast12) = paste0("dN_", colnames(fcast12))
Y_base = cbind(fcast12, b_base[as.character(2016:2023),])

obs = t(rbind(colSums(agg.pop), agg.pop)[, as.character(2016:2023)])
Y_obs = cbind(obs, b_obs[as.character(2016:2023),])
colnames(Y_obs) = colnames(Y_base)

# OLS - unconstrained

n_r = length(county.names)
m = 3*n_r
n = 1 + n_r + m

dim(S)      # n x m
dim(Amat)   # m x (2*n_r)
dim(b_base) # T x m
dim(Y_base) # T x n

Ymat = as.matrix(Y_base)

P  = ginv(t(S) %*% S) %*% t(S)
SP = S %*% P

Y_rec_OLS = t(SP %*% t(Ymat))
colnames(Y_rec_OLS) = colnames(Y_base)
Y_rec_OLS = as.data.frame(Y_rec_OLS)

# OLS - non-negative constraints

Ymat = as.matrix(Y_base)

Tn = nrow(Ymat)

Dmat = t(S) %*% S
Dmat = Dmat + 1e-20 * diag(m)   # small ridge for numerical stability

bvec = rep(0, ncol(Amat))

## Compute bottom-level reconciled
b_rec_OLSnn = matrix(NA, nrow = Tn, ncol = m)

for (t in 1:Tn) {
  
  yhat = Ymat[t, ]
  
  dvec = t(S) %*% yhat
  
  sol = solve.QP(
    Dmat = Dmat,
    dvec = dvec,
    Amat = Amat,
    bvec = bvec,
    meq  = 0
  )
  
  b_rec_OLSnn[t, ] = sol$solution
}

## Recover full hierarchy
Y_rec_OLSnn = t(S %*% t(b_rec_OLSnn))

colnames(Y_rec_OLSnn) = colnames(Y_base)
colnames(b_rec_OLSnn) = colnames(b_base)

Y_rec_OLSnn = as.data.frame(Y_rec_OLSnn)
b_rec_OLSnn = as.data.frame(b_rec_OLSnn)

# checks
all(b_rec_OLSnn[, grep("^B_", colnames(b_rec_OLSnn))] >= -0.000001)
all(b_rec_OLSnn[, grep("^D_", colnames(b_rec_OLSnn))] >= -0.000001)
all(Y_rec_OLS[, grep("^B_", colnames(Y_rec_OLS))] >= 1e-8)
all(Y_rec_OLS[, grep("^D_", colnames(Y_rec_OLS))] >= 1e-8)

View(round(rbind(Y_obs,
                 rep(NA, n),
                 Y_base,
                 rep(NA, n),
                 Y_rec_OLS,
                 rep(NA, n),
                 Y_rec_OLSnn), 1))

# write.csv(round(rbind(Y_obs,
#                       rep(NA, n),
#                       Y_base,
#                       rep(NA, n),
#                       Y_rec_OLS,
#                       rep(NA, n),
#                       Y_rec_OLSnn), 1),
#           "OLS.csv")

##########------- WLS v
#######################################
sigma2 = apply(fcast_error, 2, var)
Winv = diag(1 / sigma2)


## Unconstrained
Ymat = as.matrix(Y_base)

Tn = nrow(Ymat)
n  = ncol(Ymat)
m  = ncol(S)

P  = ginv(t(S) %*% Winv %*% S) %*% t(S) %*% Winv
SP = S %*% P

Y_rec_WLS = t(SP %*% t(Ymat))

colnames(Y_rec_WLS) = colnames(Y_base)
Y_rec_WLS = as.data.frame(Y_rec_WLS)

## Non-negative constraint 
Dmat = t(S) %*% Winv %*% S
Dmat = Dmat + 1e-8 * diag(ncol(S))   # numerical stability

b_rec_WLSnn = matrix(NA, nrow = nrow(Y_base), ncol = ncol(S))

for (t in 1:nrow(Y_base)) {
  
  yhat = as.numeric(Y_base[t, ])
  
  dvec = t(S) %*% Winv %*% yhat
  
  sol = solve.QP(
    Dmat = Dmat,
    dvec = dvec,
    Amat = Amat,
    bvec = rep(0, ncol(Amat)),
    meq  = 0
  )
  
  b_rec_WLSnn[t, ] = sol$solution
}

Y_rec_WLSnn = t(S %*% t(b_rec_WLSnn))

colnames(Y_rec_WLSnn) = colnames(Y_base)
Y_rec_WLSnn = as.data.frame(Y_rec_WLSnn)

# checks
all(b_rec_WLSnn[, grep("^B_", colnames(b_rec_WLSnn))] >= 0)
all(b_rec_WLSnn[, grep("^D_", colnames(b_rec_WLSnn))] >= 0)
all(Y_rec_WLS[, grep("^B_", colnames(Y_rec_WLS))] >= 0)
all(Y_rec_WLS[, grep("^D_", colnames(Y_rec_WLS))] >= 0)

View(round(rbind(Y_obs,
                 rep(NA, n),
                 Y_base,
                 rep(NA, n),
                 Y_rec_WLS,
                 rep(NA, n),
                 Y_rec_WLSnn), 1))

# write.csv(round(rbind(Y_obs,
#                       rep(NA, n),
#                       Y_base,
#                       rep(NA, n),
#                       Y_rec_WLS,
#                       rep(NA, n),
#                       Y_rec_WLSnn), 1),
#           "WLS.csv")
########## MinT

W_sample = cov(fcast_error)
W_inv_mint = ginv(W_sample)

## Unconstrained
Ymat = as.matrix(Y_base)

Tn = nrow(Ymat)
n  = ncol(Ymat)
m  = ncol(S)

P  = ginv(t(S) %*% W_inv_mint %*% S) %*% t(S) %*% W_inv_mint
SP = S %*% P

Y_rec_MinT = t(SP %*% t(Ymat))

colnames(Y_rec_MinT) = colnames(Y_base)
Y_rec_MinT = as.data.frame(Y_rec_MinT)

## Non-negative constraint
Dmat = t(S) %*% W_inv_mint %*% S
Dmat = Dmat + 1e-8 * diag(ncol(S))   # numerical stability

b_rec_Y_rec_MinTnn = matrix(NA, nrow = nrow(Y_base), ncol = ncol(S))

for (t in 1:nrow(Y_base)) {
  
  yhat = as.numeric(Y_base[t, ])
  
  dvec = t(S) %*% W_inv_mint %*% yhat
  
  sol = solve.QP(
    Dmat = Dmat,
    dvec = dvec,
    Amat = Amat,
    bvec = rep(0, ncol(Amat)),
    meq  = 0
  )
  
  b_rec_Y_rec_MinTnn[t, ] = sol$solution
}

Y_rec_MinTnn = t(S %*% t(b_rec_Y_rec_MinTnn))

colnames(Y_rec_MinTnn) = colnames(Y_base)
Y_rec_MinTnn = as.data.frame(Y_rec_MinTnn)

# checks
all(b_rec_Y_rec_MinTnn[, grep("^B_", colnames(b_rec_Y_rec_MinTnn))] >= 0)
all(b_rec_Y_rec_MinTnn[, grep("^D_", colnames(b_rec_Y_rec_MinTnn))] >= 0)
all(Y_rec_MinT[, grep("^B_", colnames(Y_rec_MinT))] >= 0)
all(Y_rec_MinT[, grep("^D_", colnames(Y_rec_MinT))] >= 0)

View(round(rbind(Y_obs,
                 rep(NA, n),
                 Y_base,
                 rep(NA, n),
                 Y_rec_MinT,
                 rep(NA, n),
                 Y_rec_MinTnn), 1))

# write.csv(round(rbind(Y_obs,
#                       rep(NA, n),
#                       Y_base,
#                       rep(NA, n),
#                       Y_rec_MinT,
#                       rep(NA, n),
#                       Y_rec_MinTnn), 1),
#           "MinT.csv")































