

location.file = "bayespop_workshop/data/wafips.txt"
county.codes = read.delim(location.file)$reg_code[-1] #Remove code '53' for Washington State
county.names = read.delim(location.file)$name[-1]

subnat.pop.dir = "bayespop_workshop/results/pop"
pop.subnat = get.pop.prediction(subnat.pop.dir)
pop.aggr = get.pop.aggregation(sim.dir = subnat.pop.dir) 
n.traj = 1000
n_years = 10

agg_dN = get.pop.exba("B53", pop.aggr, nr.traj = n.traj) + 
  get.pop.exba("G53", pop.aggr, nr.traj = n.traj) -
  get.pop.exba("D53", pop.aggr, nr.traj = n.traj)
agg_dN = agg_dN[5:(4+n_years), ]
dens_list <- apply(agg_dN, 1, density)

df <- as.data.frame(agg_dN)
df$year <- rownames(agg_dN)

long_df <- pivot_longer(df, -year, values_to = "value")
p_agg = ggplot(long_df, aes(x = value, y = factor(year))) +
  geom_density_ridges() +
  labs(x = "dN",
       y = "Year",
       title = "Predictive distribution of yearly population change (Washington state), 1000 samples per year") +
  theme_minimal()

county_dN = array(NA, dim = c(n_years, n.traj, length(county.codes)))
for (c in (1:length(county.codes))){
  county.code = county.codes[c]
  temp_dN = (get.pop.exba(paste0("B", county.code), pop.subnat, nr.traj = n.traj) + 
               get.pop.exba(paste0("G", county.code), pop.subnat, nr.traj = n.traj) -
               get.pop.exba(paste0("D", county.code), pop.subnat, nr.traj = n.traj)
  )[5:(4+n_years),]
  county_dN[, , c] = temp_dN
}

results.bhm.dir = "data/bhm_results_2014/"
county_dN = saveRDS(county_dN, paste0(results.bhm.dir, "county_dN.rds"))

county_dN = readRDS(paste0(results.bhm.dir, "county_dN.rds"))
K = 100000

b.tilde_array = array(NA, dim = c(n_years, K, length(county.codes)))
u.tilde_mat = matrix(NA, nrow = n_years, ncol = K)

for (t in 1:n_years){
  # pi_U_hat = dnorm
  
  Ab_hat = apply(county_dN[t, , ], 1, sum)
  
  w_unnorm = dnorm(Ab_hat, as.numeric(mu[1,t]), sd = as.numeric(sigma[1,t]))
  
  w = w_unnorm / sum(w_unnorm)
  
  sample_idx = sample(1:length(w), size = K, replace = TRUE, prob = w)
  
  b.tilde = county_dN[t, sample_idx, ]
  b.tilde_array[t, ,] = b.tilde
  
  u.tilde_mat[t, ] = rowSums(b.tilde)
}


# for (t in 1:n_years){
#   pi_U_hat = density(agg_dN[t,])
#   
#   Ab_hat = apply(county_dN[t, , ], 1, sum)
#   
#   w_unnorm = approx(pi_U_hat$x, pi_U_hat$y, xout = Ab_hat, rule = 1)$y
#   
#   w = w_unnorm / sum(w_unnorm)
#   
#   sample_idx = sample(1:length(w), size = K, replace = TRUE, prob = w)
#   
#   b.tilde = county_dN[t, sample_idx, ]
#   b.tilde_array[t, ,] = b.tilde
#   
#   u.tilde_mat[t, ] = rowSums(b.tilde)
# }


############ Obs

subnat.pop.dir = "bayespop_workshop/results_24/pop"
pop.subnat = get.pop.prediction(subnat.pop.dir)
pop.aggr = get.pop.aggregation(sim.dir = subnat.pop.dir) 

obs_dN = (get.pop.exba("B53", pop.aggr, observed = TRUE) + 
  get.pop.exba("G53", pop.aggr, observed = TRUE) -
  get.pop.exba("D53", pop.aggr, observed = TRUE))[, as.character(("2014":"2023"))]

obs_county_dN = matrix(NA, nrow = n_years, ncol = length(county.codes))
for (c in (1:length(county.codes))){
  county.code = county.codes[c]
  temp_dN = (get.pop.exba(paste0("B", county.code), observed = TRUE, pop.subnat) + 
               get.pop.exba(paste0("G", county.code), observed = TRUE, pop.subnat) -
               get.pop.exba(paste0("D", county.code), observed = TRUE, pop.subnat)
  )[,as.character(("2017":"2026"))]
  obs_county_dN[, c] = temp_dN
}


############Plot

df <- as.data.frame(agg_dN)
df$year <- rownames(agg_dN)

long_df <- pivot_longer(df, -year, values_to = "value")

vline_df <- data.frame(
  year = factor(unique(long_df$year)),
  xintercept = obs_dN  # one per year
)

p1 = ggplot(long_df, aes(x = value, y = factor(year))) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.5),
                      vline_width = 0.5,
                      alpha = 0.7) +
  # geom_segment(
  #   data = vline_df,
  #   aes(
  #     x = xintercept,
  #     xend = xintercept,
  #     y = as.numeric(year) - 0.1,  # adjust vertical placement
  #     yend = as.numeric(year) + 0.8
  #   ),
  #   color = "red",
  #   linetype = "dashed"
  # ) +
  scale_x_continuous(limits = c(25000, 125000)) +
  labs(x = "dN",
       y = "Year",
       title = "Predictive distribution (yearly dN) (Washington state), 1000 samples") +
  theme_minimal()


df2 <- as.data.frame(u.tilde_mat)
df2$year <- rownames(agg_dN)

long_df2 <- pivot_longer(df2, -year, values_to = "value")
p2 = ggplot(long_df2, aes(x = value, y = factor(year))) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.5),
                      vline_width = 0.5,
                      alpha = 0.7) +
  # geom_segment(
  #   data = vline_df,
  #   aes(
  #     x = xintercept,
  #     xend = xintercept,
  #     y = as.numeric(year) - 0.1,  # adjust vertical placement
  #     yend = as.numeric(year) + 0.8
  #   ),
  #   color = "red",
  #   linetype = "dashed"
  # ) +
  scale_x_continuous(limits = c(25000, 125000)) +
  labs(x = "dN",
       y = "Year",
       title = "Coherent distribution (yearly dN) (Washington state), 1000 samples") +
  theme_minimal()

p1+p2
p_ofm + p2
### for b hat
county_c = 17
df <- as.data.frame(county_dN[, , county_c])
df$year <- rownames(agg_dN)

long_df <- pivot_longer(df, -year, values_to = "value")

vline_df <- data.frame(
  year = factor(unique(long_df$year)),
  xintercept = obs_county_dN[, county_c]  # one per year
)

p1 = ggplot(long_df, aes(x = value, y = factor(year))) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.5),
                      vline_width = 0.5,
                      alpha = 0.7) +
  geom_segment(
    data = vline_df,
    aes(
      x = xintercept,
      xend = xintercept,
      y = as.numeric(year) - 0.1,  # adjust vertical placement
      yend = as.numeric(year) + 0.8
    ),
    color = "red",
    linetype = "dashed"
  ) +
  scale_x_continuous(limits = c(-10000, 60000)) +
  labs(x = "dN",
       y = "Year",
       title = paste0("Predictive distribution (yearly dN) (", county.names[county_c], " county), 1000 samples")) +
  theme_minimal()


df2 <- as.data.frame(b.tilde_array[, , county_c])
df2$year <- rownames(agg_dN)

long_df2 <- pivot_longer(df2, -year, values_to = "value")
p2 = ggplot(long_df2, aes(x = value, y = factor(year))) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.5),
                      vline_width = 0.5,
                      alpha = 0.7) +
  geom_segment(
    data = vline_df,
    aes(
      x = xintercept,
      xend = xintercept,
      y = as.numeric(year) - 0.1,  # adjust vertical placement
      yend = as.numeric(year) + 0.8
    ),
    color = "red",
    linetype = "dashed"
  ) +
  scale_x_continuous(limits = c(-10000, 60000)) +
  labs(x = "dN",
       y = "Year",
       title = paste0("Coherent distribution (yearly dN) (", county.names[county_c], " county), 1000 samples")) +
  theme_minimal()
library(patchwork)
p1+p2


######
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
county_id <- 2

mu_vec <- as.numeric(mu[county_id, ])
sd_vec <- as.numeric(sigma[county_id, ])

lower_bounds <- qnorm(0.005, mean = mu_vec, sd = sd_vec)
upper_bounds <- qnorm(0.995, mean = mu_vec, sd = sd_vec)

x_min <- min(lower_bounds, na.rm = TRUE)
x_max <- max(upper_bounds, na.rm = TRUE)

x_grid <- seq(x_min, x_max, length.out = 300)

years <- colnames(mu)

plot_df <- do.call(rbind, lapply(seq_along(years), function(j) {
  data.frame(
    x = x_grid,
    year = years[j],
    density = dnorm(x_grid, mean = mu_vec[j], sd = sd_vec[j])
  )
}))

vline_df <- data.frame(
  year = years,
  low  = as.numeric(low_diff[county_id, ]),
  mid  = as.numeric(medium_diff[county_id, ]),
  high = as.numeric(high_diff[county_id, ])
) %>%
  pivot_longer(cols = c(low, mid, high),
               names_to = "scenario",
               values_to = "xintercept")

p_ofm = ggplot(plot_df, aes(x = x, y = factor(year), height = density)) +
  geom_density_ridges(
    stat = "identity",
    scale = 1.5,
    alpha = 0.7
  ) +
  geom_segment(
    data = vline_df,
    inherit.aes = FALSE,
    aes(
      x = xintercept,
      xend = xintercept,
      y = as.numeric(factor(year)) - 0.2,
      yend = as.numeric(factor(year)) + 0.6,
      color = scenario
    ),
    linetype = "dashed"
  ) +
  labs(
    x = "Yearly population change",
    y = "Year",
    color = "Scenario"
  ) +
  theme_minimal()


######### Counties_Comp

compcounty_dN = array(NA, dim = c(n_years, n.traj, length(county.codes), 3))
for (c in (1:length(county.codes))){
  county.code = county.codes[c]
  compcounty_dN[, , c, 1] = get.pop.exba(paste0("B", county.code), pop.subnat, nr.traj = n.traj)[5:(4+n_years),]
  compcounty_dN[, , c, 2] = get.pop.exba(paste0("G", county.code), pop.subnat, nr.traj = n.traj)[5:(4+n_years),]
  compcounty_dN[, , c, 3] = get.pop.exba(paste0("D", county.code), pop.subnat, nr.traj = n.traj)[5:(4+n_years),]
  
}

results.bhm.dir = "data/bhm_results_2014/"
# saveRDS(compcounty_dN, paste0(results.bhm.dir, "compcounty_dN.rds"))
compcounty_dN = readRDS(paste0(results.bhm.dir, "compcounty_dN.rds"))

K = 1000

b.tilde_array = array(NA, dim = c(n_years, K, length(county.codes), 3))
u.tilde_mat = matrix(NA, nrow = n_years, ncol = K, length(county.codes))

for (t in 1:n_years){
  # Lower level
  for (c in 1:length(county.codes)){
    county.code = county.codes[c]
    Ab_hat = compcounty_dN[t, , c, 1] +
             compcounty_dN[t, , c, 2] -
             compcounty_dN[t, , c, 3]
    
    w_unnorm = dnorm(Ab_hat, as.numeric(mu[c+1,t]), sd = as.numeric(sigma[c+1,t]))
    
    w = w_unnorm / sum(w_unnorm)
    
    sample_idx = sample(1:length(w), size = K, replace = TRUE, prob = w)
    
    b.tilde = compcounty_dN[t, sample_idx, c, ]
    b.tilde_array[t, , c, ] = compcounty_dN[t, sample_idx, c, ]
  }
  
  Ab_hat = apply(county_dN[t, , ,], 1, sum)
  
  w_unnorm = dnorm(Ab_hat, as.numeric(mu[1,t]), sd = as.numeric(sigma[1,t]))
  
  w = w_unnorm / sum(w_unnorm)
  
  sample_idx = sample(1:length(w), size = K, replace = TRUE, prob = w)
  
  b.tilde = county_dN[t, sample_idx, ]
  b.tilde_array[t, ,] = b.tilde
  
  u.tilde_mat[t, ] = rowSums(b.tilde)
}

######### Counties

subcounty_dN = array(NA, dim = c(n_years, n.traj, length(county.codes)))
for (c in (1:length(county.codes))){
  county.code = county.codes[c]
  temp_dN = (get.pop.exba(paste0("B", county.code), pop.subnat, nr.traj = n.traj) + 
               get.pop.exba(paste0("G", county.code), pop.subnat, nr.traj = n.traj) -
               get.pop.exba(paste0("D", county.code), pop.subnat, nr.traj = n.traj)
  )[5:(4+n_years),]
  county_dN[, , c] = temp_dN
}

get.pop.exba(paste0("B", county.code), pop.subnat, nr.traj = 1)



test = get.pop.exba(paste0("B", county.code,"_{10:55}"), pop.subnat, nr.traj = 1)
colSums(test)


colSums(get.pop.exba(paste0("D", county.code,"_{0:130}"), pop.subnat, nr.traj = 1))

get.pop.exba(paste0("D", county.code,"_{0:130}"), pop.subnat, nr.traj = 2)


colSums(get.pop.exba(paste0("G", county.code,"_{0:100}"), pop.subnat, nr.traj = 1))
get.pop.exba(paste0("G", county.code), pop.subnat, nr.traj = 1)






















