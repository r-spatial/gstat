require(tools)
packages_to_check <- function(dep, which = c("Depends", "Imports", "LinkingTo", "Suggests"), recursive = FALSE){

download.file("http://cran.R-project.org/web/packages/packages.rds", "packages.rds", mode="wb")
    x <- readRDS("packages.rds")
    x <- x[!duplicated(x[,1]),]
    packages <- x[,1]
    rdeps <- package_dependencies(packages = dep, x,
                        which = which,
                        recursive = recursive, reverse = TRUE)
    paste(apply(x[x[,1] %in% rdeps[[1]], 1:2], 1, paste, collapse="_"), ".tar.gz", sep="")
}

#RCheck = function(x, URL = "http://ftp5.gwdg.de/pub/misc/cran/src/contrib/") {
RCheck = function(x, URL = "http://cran.r-project.org/src/contrib/") {
	if (!file.exists(x))
		download.file(paste(URL, x, sep=""), x)
	cmd = paste("R CMD check --as-cran ", x, " > ", x, ".log", sep = "")
	print(cmd)
	ret = system(cmd)
	print(ret)
	ret
}

result <-  packages_to_check("gstat")
result
sel = TRUE
library(parallel)
ncores_to_use = 2
cl <- makeCluster(getOption("cl.cores", ncores_to_use))
clusterExport(cl, c("RCheck", "sel", "result"))
out = parLapply(cl, result[sel], function(x) RCheck(x))
succ = unlist(out)
x = which(succ != 0)
result[x]
bla = lapply(result[x], function(y) {
		cat(paste(y, ":\n"))
		system(paste("tail -20 ",y,".log", sep=""))
	}
)

#result <-  packages_to_check("sp", recursive=TRUE) 
