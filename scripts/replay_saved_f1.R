# Reaggregate saved results only. Original run outputs remain immutable.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: Rscript scripts/replay_saved_f1.R REPO_DIR NEW_OUTPUT_DIR")
root <- normalizePath(args[1], mustWork=TRUE)
src <- file.path(root, 'reference/legacy_summary')
out <- args[2]
if (file.exists(out)) stop("Output must not exist")
dir.create(out, recursive=TRUE)
m <- read.csv(file.path(src,'all_replication_metrics.csv'))
old <- m$f1
idx <- m$dgp %in% c(3,4)
m$f1[idx] <- with(m[idx,], 2*true_positive/(2*true_positive+false_positive+false_negative))
changed <- which(idx & (is.na(old) | abs(m$f1-old)>1e-12))
stopifnot(length(changed)==1, m$dgp[changed]==3, m$replication_id[changed]==91, m$f1[changed]==0)
write.csv(m,file.path(out,'all_replication_metrics.csv'),row.names=FALSE,na='')
s <- read.csv(file.path(src,'dgp3_dgp4_method_summaries.csv'))
for(i in seq_len(nrow(s))) {v<-m$f1[m$dgp==s$dgp[i]&m$method==s$method[i]];s$mean_f1[i]<-mean(v);s$f1_defined_reps[i]<-sum(is.finite(v))}
write.csv(s,file.path(out,'dgp3_dgp4_method_summaries.csv'),row.names=FALSE)
# Recreate original bootstrap RNG stream, and reuse the DGP3/F1 entry state.
# Restore the original post-entry state to avoid changing unrelated CIs by RNG drift.
p <- read.csv(file.path(src,'paired_differences_by_replication.csv'))
c <- read.csv(file.path(src,'prespecified_paired_contrasts.csv'))
set.seed(5200001L)
for(i in seq_len(nrow(c))) {
 d<-p[p$contrast==c$contrast[i],];d<-d[order(d$replication_id),]
 v<-d[[c$metric[i]]];v<-v[is.finite(v)];before<-.Random.seed
 b<-replicate(9999L,mean(sample(v,length(v),replace=TRUE)));after<-.Random.seed
 ci<-quantile(b,c(.025,.975),names=FALSE,type=7)
 stopifnot(abs(mean(v)-c$mean_difference[i])<1e-12, max(abs(ci-c(c$bootstrap_ci_lower[i],c$bootstrap_ci_upper[i])))<1e-12)
 if(c$dgp[i]==3 && c$metric[i]=='f1') {
  a<-m[m$dgp==3&m$method==c$lhs[i],];z<-m[m$dgp==3&m$method==c$rhs[i],]
  a<-a[order(a$replication_id),];z<-z[order(z$replication_id),];stopifnot(identical(a$replication_id,z$replication_id))
  v<-a$f1-z$f1;stopifnot(length(v)==100,all(is.finite(v)))
  .Random.seed<-before;b<-replicate(9999L,mean(sample(v,length(v),replace=TRUE)));ci<-quantile(b,c(.025,.975),names=FALSE,type=7)
  c$n_paired_datasets[i]<-100;c$mean_difference[i]<-mean(v);c$bootstrap_ci_lower[i]<-ci[1];c$bootstrap_ci_upper[i]<-ci[2]
  p$f1[p$contrast==c$contrast[i]]<-v
  .Random.seed<-after
 }
}
write.csv(c,file.path(out,'prespecified_paired_contrasts.csv'),row.names=FALSE)
write.csv(p,file.path(out,'paired_differences_by_replication.csv'),row.names=FALSE,na='')
print(s[,c('dgp','method','mean_f1','f1_defined_reps')]);print(c[c$dgp==3 & c$metric=='f1',])
# Optional audit interval for the corrected mean; not added to the paper's mean-only table.
set.seed(5200002L);v<-m$f1[m$dgp==3&m$method=='path_continuous'];ci<-quantile(replicate(9999L,mean(sample(v,length(v),replace=TRUE))),c(.025,.975))
write.csv(data.frame(n=100,mean=mean(v),lower=ci[1],upper=ci[2],bootstrap_seed=5200002),file.path(out,'dgp3_f1_mean_audit.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
