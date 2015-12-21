library(sp)

data(meuse)
data(meuse.grid)

set.seed(13131) # reproduce results

# select 10 random rows;
# create two missing values in the coordinates:
m = meuse.grid[sample(nrow(meuse.grid), 10), ]
m[c(2,8), "x"] = NA

library(gstat)
## this is not allowed anymore:
try(krige(log(zinc)~1,~x+y,meuse,m, na.action = na.pass))
try(krige(log(zinc)~1,~x+y,meuse,m, na.action = na.omit))
try(krige(log(zinc)~1,~x+y,meuse,m, na.action = na.exclude))
try(krige(log(zinc)~1,~x+y,meuse,m, na.action = na.fail))

# select 10 random rows;
# create two missing values in the regressor variable:
m = meuse.grid[sample(nrow(meuse.grid), 10), ]
m[c(3,7), "dist"] = NA
krige(log(zinc)~dist,~x+y,meuse,m, na.action = na.pass)
krige(log(zinc)~dist,~x+y,meuse,m, na.action = na.omit)
krige(log(zinc)~dist,~x+y,meuse,m, na.action = na.exclude)
try(krige(log(zinc)~dist,~x+y,meuse,m, na.action = na.fail))
