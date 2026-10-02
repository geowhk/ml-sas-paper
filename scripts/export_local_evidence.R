args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: export_local_evidence.R REPO_DIR REPEATED_RUN_DIR NEW_OUTPUT_DIR")
root<-normalizePath(args[1]);run<-normalizePath(args[2]);out<-args[3]
if(file.exists(out)) stop("Output must not exist")
dir.create(out,recursive=TRUE)
library(sf)
source(file.path(root,'analysis/repeated/R/generate_professor_dgps.R'))
grid <- make_professor_grid_data(20,20,'hex',TRUE,1)
fixed_files <- c('sim_dgp1.rds','sim_dgp2.rds','sim_dgp3.1.rds','sim_dgp4.3.rds')
functions <- list(simulate_professor_dgp1,simulate_professor_dgp2,simulate_professor_dgp3,simulate_professor_dgp4)
checks <- list()
for (d in 1:4) {
 generated <- functions[[d]](grid,seed=2025L+d)
 saved <- readRDS(file.path(root,'02_revision-simulation/sim_dgp',fixed_files[d]))
 stopifnot(setequal(names(generated),names(saved)))
 for(nm in setdiff(names(generated),'geometry')) {
  a<-generated[[nm]]; b<-saved[[nm]]
  if(is.numeric(a)) stopifnot(isTRUE(all.equal(a,b,tolerance=5e-15))) else stopifnot(identical(a,b))
 }
 stopifnot(all(vapply(seq_len(nrow(saved)),function(i) i %in% st_equals_exact(generated[i,],saved,par=1e-12)[[1]],logical(1))))
 write.csv(st_drop_geometry(saved),file.path(out,paste0('fixed_dgp',d,'_attributes.csv')),row.names=FALSE)
 if(d==1) {
  coords <- do.call(rbind,lapply(seq_len(nrow(saved)),function(i) data.frame(id=saved$id[i],st_coordinates(saved[i,])[,1:2])))
  write.csv(coords,file.path(out,'hex_vertices.csv'),row.names=FALSE)
 }
 checks[[d]]<-data.frame(dgp=d,fixed_file=fixed_files[d],seed=2025+d,n=nrow(saved),targets=sum(saved$true_cluster),matched=TRUE)
}
write.csv(do.call(rbind,checks),file.path(out,'fixed_reproduction_checks.csv'),row.names=FALSE)
cat('PASS: fixed DGP 1–4 values, categories, truth labels and geometry match saved RDS files.\n')
truth <- make_professor_true_cluster(grid$sf)$true_cluster
stopifnot(sum(truth)==60L)
regions <- st_drop_geometry(grid$sf)
regions$target <- truth
regions$degree <- lengths(grid$nb)
regions$target_neighbors <- vapply(grid$nb,function(j)sum(truth[j]),integer(1))
regions$target_neighbor_fraction <- regions$target_neighbors/regions$degree
write.csv(regions,file.path(out,'lattice_and_target.csv'),row.names=FALSE)
all <- list(); k<-0L
for(d in 3:4) for(rep in 1:100) {
 sfdata <- generate_professor_dgp(d,rep)
 attrs <- st_drop_geometry(sfdata)
 methods <- if(d==3)c('local_geary','path_continuous') else c('local_geary','path_continuous','path_mixed','terminal_mixed')
 for(method in methods) {
  local <- read.csv(file.path(run,paste0('dgp',d),sprintf('rep_%04d',rep),paste0(method,'_local.csv')))
  stopifnot(identical(as.integer(local$id),as.integer(attrs$id)))
  local$dgp<-d;local$replication_id<-rep;local$method<-method
  local$truth<-attrs$true_cluster;local$x1<-attrs$x1;local$x2<-attrs$x2
  local$cat1<-if(d==4)as.character(attrs$cat1) else NA_character_
  local$target_neighbor_fraction<-regions$target_neighbor_fraction
  if(!'permutation_mean'%in%names(local)) local$permutation_mean<-NA_real_
  local$p_raw <- if(method=='local_geary') local$p_value else local$p_upper
  k<-k+1L
  all[[k]]<-local[,c('dgp','replication_id','method','id','truth','x1','x2','cat1','target_neighbor_fraction','statistic','permutation_mean','p_raw','p_adjusted','detected')]
 }
}
write.csv(do.call(rbind,all),file.path(out,'repeated_local_evidence.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('PASS: exported 600 method-dataset local results with regenerated attributes, without rerunning forests or permutations.\n')
