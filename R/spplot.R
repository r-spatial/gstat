spplot.vcov = function(x, ...) {
    basenames = sub(".pred", "", names(x)[grep(".pred",names(x))])
    n = length(basenames)
    names = NULL
    skip = NULL
    for (i in 1:n) {
        skp = rep(TRUE, n)
        pos = 1
        if (i > 1) {
            for (j in 1:(i-1)) {
                names = c(names, paste("cov", basenames[j], basenames[i], sep = "."))
                skp[pos] = FALSE
                pos = pos + 1
            }
        }
        names = c(names, paste(basenames[i], ".var", sep = ""))
        skp[pos] = FALSE
        skip = c(skip, skp)
    }
    spplot(x, names, skip = skip, layout = c(n,n), as.table = TRUE, ...)
}
