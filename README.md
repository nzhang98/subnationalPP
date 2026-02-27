Repository Forecast Reconciliation of Bayesian sub-national probabilistic population projections.

"synthetic_data_test.R" contains a simple example with synthetic data of Forecast Reconciliation methods
"extract_data.R" contains the steps to extract the relevant data from bayesPop and OFM projections/historical data from the folder data/.
          Note that to extract data from bayesPop, need to have results on subnational data as explained in https://github.com/PPgp/CSDEworkshop
"utils.R" contains useful functions used throughout;
"ForecastReconciliation_DifferentW.R" contains the first attempt at reconciling aggregate B,D,G (posterior means) from the BHM with aggregate county forecasts from the OFM.

After cloning, to match dependencies use

renv::activate()
renv::restore()
