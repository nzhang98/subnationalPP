install.packages("quadprog")  # uncomment if not installed
library(quadprog)

set.seed(123)

# Bottom-level series: B, D, G for 2 regions
m <- 6

# Random "base forecasts" for bottom-level (all counts >=0)
B_r1 <- 100
D_r1 <- 50
G_r1 <- 10

B_r2 <- 20
D_r2 <- 3
G_r2 <- -10

b_base <- c(B_r1, D_r1, G_r1, B_r2, D_r2, G_r2)

Delta_r1_base <- 100   # independent base forecast for region 1
Delta_r2_base <- 30   # independent base forecast for region 2
Delta_c_base  <- 130  # independent base forecast for country

# Combined base forecast vector (top to bottom)
y_base <- c(Delta_c_base, Delta_r1_base, Delta_r2_base, b_base)
y_base

S <- matrix(c(
  1, -1,  1,  1, -1, 1,   # ΔN_c
  1, -1,  1,  0,  0, 0,   # ΔN_r1
  0,  0,  0,  1, -1, 1,   # ΔN_r2
  1,  0,  0,  0,  0, 0,   # B_r1
  0,  1,  0,  0,  0, 0,   # D_r1
  0,  0,  1,  0,  0, 0,   # G_r1
  0,  0,  0,  1,  0, 0,   # B_r2
  0,  0,  0,  0,  1, 0,   # D_r2
  0,  0,  0,  0,  0, 1    # G_r2
), nrow=9, byrow=TRUE)


# OLS weights = identity
W <- diag(nrow(S))
# -------------------
# WLS_s: structural scaling
# -------------------
k_h <- 1  # assume unit variance
Lambda <- diag(as.numeric(S %*% rep(1,6)))  # number of bottom-level series per aggregate
W <- k_h * Lambda
# -------------------
# WLS_v: variance scaling
# -------------------
# simulate variance for each series (for demo)
var_vec <- runif(9, 1, 5)
W <- diag(var_vec)
# -------------------
# MinT (sample covariance) - demo with small random covariance
# -------------------
# Here we simulate a covariance matrix for bottom-level forecasts
cov_b <- matrix(runif(36,0,10),6,6)
cov_b <- cov_b %*% t(cov_b)/100  # make it symmetric positive-definite
# Expand to all levels via S
W <- S %*% cov_b %*% t(S)

H <- t(S) %*% ginv(W) %*% S    # same as t(S) %*% S for W=I
dvec <- t(S) %*% ginv(W) %*% y_base

# Constraints:
# B >=0 , D >=0, G free
# Bottom-level order: B_r1, D_r1, G_r1, B_r2, D_r2, G_r2
Amat <- matrix(0, nrow=6, ncol=4)
Amat[1,1] <- 1  # B_r1 >=0
Amat[2,2] <- 1  # D_r1 >=0
Amat[4,3] <- 1  # B_r2 >=0
Amat[5,4] <- 1  # D_r2 >=0

b0 <- rep(0, 4)

sol <- solve.QP(Dmat = 2*H, dvec = dvec, Amat = Amat, bvec = b0)
b_hat <- sol$solution
b_hat

y_tilde <- S %*% b_hat
y_tilde

df_OLS = data.frame(
  Series = c("ΔN_c","ΔN_r1","ΔN_r2","B_r1","D_r1","G_r1","B_r2","D_r2","G_r2"),
  BaseForecast = round(y_base,1),
  Reconciled = round(as.numeric(y_tilde),1)
)

df_WLSs = data.frame(
  Series = c("ΔN_c","ΔN_r1","ΔN_r2","B_r1","D_r1","G_r1","B_r2","D_r2","G_r2"),
  BaseForecast = round(y_base,1),
  Reconciled = round(as.numeric(y_tilde),1)
)

df_WLSv = data.frame(
  Series = c("ΔN_c","ΔN_r1","ΔN_r2","B_r1","D_r1","G_r1","B_r2","D_r2","G_r2"),
  BaseForecast = round(y_base,1),
  Reconciled = round(as.numeric(y_tilde),1)
)

df_MinT = data.frame(
  Series = c("ΔN_c","ΔN_r1","ΔN_r2","B_r1","D_r1","G_r1","B_r2","D_r2","G_r2"),
  BaseForecast = round(y_base,1),
  Reconciled = round(as.numeric(y_tilde),1)
)

cbind(df_OLS, df_WLSs[,3], df_WLSv[,3], df_MinT[,3])
