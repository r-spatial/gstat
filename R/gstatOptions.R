.gstatOptions <- new.env(FALSE, globalenv())

assign("gstat_progress", TRUE, envir = .gstatOptions)

get_gstat_progress <- function() {
    get("gstat_progress", envir = .gstatOptions)
}

set_gstat_progress <- function(value) {
        stopifnot(is.logical(value))
        stopifnot(length(value) == 1)
        assign("gstat_progress", value, envir = .gstatOptions)
        get_gstat_progress()
}
