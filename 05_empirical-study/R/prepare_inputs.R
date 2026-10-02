prepare_empirical_inputs <- function() {
  e<-new.env();data("NLD_muni",package="tmap",envir=e);x<-e$NLD_muni
  stopifnot(nrow(x)==345L,!anyNA(x),all(sf::st_is_valid(x)))
  nb<-suppressWarnings(spdep::poly2nb(x,queen=TRUE))
  components<-spdep::n.comp.nb(nb);sizes<-table(components$comp.id)
  x<-x[which(components$comp.id==as.integer(names(sizes)[which.max(sizes)])),]
  nb<-spdep::poly2nb(x,queen=TRUE)
  stopifnot(nrow(x)==337L,all(spdep::card(nb)>0))
  frame<-sf::st_drop_geometry(x)
  frame$urbanity<-ordered(as.character(frame$urbanity),levels=levels(e$NLD_muni$urbanity))
  continuous<-as.data.frame(frame[,c("population","dwelling_value","employment_rate"),drop=FALSE])
  mixed<-as.data.frame(frame[,c("population","dwelling_value","employment_rate","urbanity"),drop=FALSE])
  rownames(continuous)<-rownames(mixed)<-frame$code
  list(sf=x,queen=nb,continuous=continuous,mixed=mixed)
}
empirical_fingerprint <- function(x) {
  h<-function(z) digest::digest(z,algo="sha256",serialize=TRUE,serializeVersion=2)
  list(n=nrow(x$sf),attributes=h(x$mixed),geometry=h(sf::st_as_binary(sf::st_geometry(x$sf))),
       adjacency=h(lapply(x$queen,as.integer)),crs=sf::st_crs(x$sf)$epsg)
}
check_empirical_inputs <- function(x,root) {
  expected<-jsonlite::read_json(file.path(root,"reference/empirical_fingerprint.json"),simplifyVector=TRUE)
  actual<-empirical_fingerprint(x)
  if(!isTRUE(all.equal(actual,expected,check.attributes=FALSE))) stop("Empirical input fingerprint differs from v6; do not continue")
  invisible(actual)
}
