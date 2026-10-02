args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else '.')
source(file.path(root,'analysis/repeated/R/production_common.R'))
reference <- read.csv(file.path(root,'reference/v6/all_replication_metrics.csv'))
for(i in which(reference$dgp %in% c(3,4))) {
 z<-reference[i,]; tp<-z$true_positive;fp<-z$false_positive;fn<-z$false_negative
 truth<-c(rep(TRUE,tp+fn),rep(FALSE,400-tp-fn))
 detected<-c(rep(TRUE,tp),rep(FALSE,fn),rep(TRUE,fp),rep(FALSE,400-tp-fn-fp))
 m<-make_method_metrics(data.frame(detected=detected),truth,z$dgp,z$replication_id,z$method,as.numeric(truth))
 stopifnot(abs(m$f1-z$f1)<1e-12, identical(is.na(m$precision),is.na(z$precision)))
}
stopifnot(safe_rank_auc(orient_rank_auc_score(c(1,2,3,4),'lower'),c(TRUE,TRUE,FALSE,FALSE))==1)
cat('PASS: all 600 saved non-null confusion matrices use v6 F1; conditional precision and AUC direction\n')
