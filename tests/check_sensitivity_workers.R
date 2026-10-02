args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
setwd(root);source('03_hyperparameter-sensitivity/R/run_mixed_sensitivity.R')
config<-yaml::read_yaml('03_hyperparameter-sensitivity/config/mixed_hyperparameter_path.yml')
inputs<-load_mixed_dgp43_inputs(config$dgp_id,config$variables)
row<-as.list(make_grid_from_config(config)[1,])
run<-function(x)run_one_condition(x,inputs,config)$metrics
future::plan(future::sequential)
a<-future.apply::future_lapply(list(row),run,future.seed=TRUE)[[1]]
future::plan(future::multisession,workers=2)
b<-future.apply::future_lapply(list(row),run,future.seed=TRUE)[[1]]
future::plan(future::sequential)
stopifnot(isTRUE(all.equal(a,b,tolerance=0)))
# Frozen first-condition values; reference design uses L'Ecuyer-CMRG in future workers.
stopifnot(abs(a$true_local_mean-1.7637210776935606)<1e-12,
          abs(a$false_local_mean-1.0194619677011947)<1e-12,
          abs(a$top10_true_jaccard-0.5151515151515151)<1e-12)
cat('PASS: sensitivity first condition matches v6 and is invariant to worker count\n')
