args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: prepare_empirical.R REPO_DIR NEW_OUTPUT_DIR")
root<-normalizePath(args[1]);out<-args[2]
if(file.exists(out))stop("Output must not exist")
source(file.path(root,"05_empirical-study/R/prepare_inputs.R"))
x<-prepare_empirical_inputs();fingerprint<-check_empirical_inputs(x,root)
dir.create(out,recursive=TRUE)
saveRDS(x,file.path(out,"analysis_inputs.rds"))
y<-x$sf[,c("code","population","dwelling_value","employment_rate","urbanity")]
sf::st_write(y,file.path(out,"empirical_inputs.geojson"),quiet=TRUE)
jsonlite::write_json(fingerprint,file.path(out,"fingerprint.json"),auto_unbox=TRUE,pretty=TRUE)
cat("PASS: package-sourced data match frozen v6 attributes, geometry and adjacency (337 units)\n")
