########################################################
# note that demo(stkrige) and demo(stkrige-prediction) #
# need to be run before this one                       #
########################################################

## cross-validation
crossStat <- function(var1, var2="PM10", STxDF=DE_RB_2005, digits=NA) {
  diff <- STxDF[,,var1,drop=F]@data[[1]] - STxDF[,,var2,drop=F]@data[[1]]
  RMSE <- sqrt(mean(diff^2))
  MAE <- mean(abs(diff))
  ME <- mean(diff)
  COR <- cor(STxDF[,,var1,drop=F]@data[[1]], STxDF[,,var2,drop=F]@data[[1]])
  res <- c(RMSE, MAE, ME, COR)
  names(res) <- c("RMSE", "MAE", "ME", "COR")
  if(is.na(digits))
    return(res)
  else
    return(round(res, digits))
}

# purely spatial:
DE_RB_2005 <- as(DE_RB_2005, "STSDF")

pureSp <- NULL
for(i in 1:365) { # i <- 1
  pureSp <- c(pureSp, krige.cv(PM10~1,DE_RB_2005[,i,"PM10"],
                               model=spVgmMod, nmax=10)$var1.pred)
}

DE_RB_2005@data$pureSp10Nghbr <- pureSp

pureSp <- NULL
for(i in 1:365) { # i <- 1
  pureSp <- c(pureSp, krige.cv(PM10~1,DE_RB_2005[,i,"PM10"],model=spVgmMod,nmax=50)$var1.pred)
}

DE_RB_2005@data$pureSp50Nghbr <- pureSp

## spatio-temporal LOOCV
target <- as(DE_RB_2005[,,"PM10"],"STFDF")

## seprable model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) {
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSepModel, nmax=10,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$sepModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSepModel, nmax=50,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$sepModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

## product-sum model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitProdSumModel, nmax=10,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$psModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitProdSumModel, nmax=50,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$psModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

## metric model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitMetricModel, nmax=10,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$metricModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitMetricModel, nmax=50,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$metricModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

## sum-metric model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSumMetricModel, nmax=10,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$sumMetricModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- array(NA, c(length(DE_RB_2005@sp), 365,2))
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 38
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"],] <- as.matrix(krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                                     newdata=DE_RB_2005[loc,,drop=F], 
                                                     fitSumMetricModel, nmax=50,
                                                     computeVar=T,
                                                     stAni=linStAni*1000/24/3600)@data[,c("var1.pred","var1.var")])
}

DE_RB_2005@data$sumMetricModel50Nghbr <- as.vector(res[,,1])[!is.na(target@data)]
DE_RB_2005@data$sumMetricModel50NghbrVar <- as.vector(res[,,2])[!is.na(target@data)]
DE_RB_2005@data$sumMetricModel50Nghbr95u <- apply(DE_RB_2005@data, 1, 
                                                  function(x) {
                                                    qnorm(0.975, x["sumMetricModel50Nghbr"],
                                                          sqrt(x["sumMetricModel50NghbrVar"]))
                                                  })
DE_RB_2005@data$sumMetricModel50Nghbr95l <- apply(DE_RB_2005@data, 1, 
                                                  function(x) {
                                                    qnorm(0.025, x["sumMetricModel50Nghbr"],
                                                          sqrt(x["sumMetricModel50NghbrVar"]))
                                                  })

## simple sum-metric model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSimpleSumMetricModel, nmax=10,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$simpleSumMetricModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSimpleSumMetricModel, nmax=50,
                                           stAni=fitMetricModel$stAni/24/3600)$var1.pred
}

DE_RB_2005@data$simpleSumMetricModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

###
# cross-stats

rbind(
crossStat("pureSp10Nghbr", digits=2),
crossStat("pureSp50Nghbr", digits=2),

crossStat("sepModel10Nghbr", digits=2),
crossStat("sepModel50Nghbr", digits=2),

crossStat("psModel10Nghbr", digits=2),
crossStat("psModel50Nghbr", digits=2),

crossStat("metricModel10Nghbr", digits=2),
crossStat("metricModel50Nghbr", digits=2),

crossStat("sumMetricModel10Nghbr", digits=2),
crossStat("sumMetricModel50Nghbr", digits=2))

if(paper) {
  texRow <- function(x) {
    paste(paste(x,collapse=" & ")," \\\\ \n")
  }
  
  cat(apply(round(rbind(crossStat("pureSp10Nghbr"),
                        
                        crossStat("sepModel10Nghbr"),
                        crossStat("psModel10Nghbr"),
                        crossStat("metricModel10Nghbr"),
                        crossStat("sumMetricModel10Nghbr"),
                        crossStat("simpleSumMetricModel50Nghbr"),
                        
                        crossStat("pureSp50Nghbr"),
                        
                        crossStat("sepModel50Nghbr"),
                        crossStat("psModel50Nghbr"),
                        crossStat("metricModel50Nghbr"),
                        crossStat("sumMetricModel50Nghbr"),
                        crossStat("simpleSumMetricModel50Nghbr")
                        ),
                  2),1,texRow))

  loc <- 38 # sample(68,1) # 15
  tw <- "2005-01-15/2005-04-15"
  
  png("vignettes/figures/singleStationTimeSeries.png", 9, 4, "in", bg="white", res = 149)
  plot(DE_RB_2005[loc,tw][,"sumMetricModel50Nghbr"], main=paste("Location", DE_RB_2005@sp@data$station_european_code[loc]),
       ylim=c(0,70))
  points(DE_RB_2005[loc,tw][,"PM10"], type="l", col="darkgreen", lty=1)
  points(DE_RB_2005[loc,tw][,"sumMetricModel50Nghbr95u"], type="l", col="darkgrey", lty=2)
  points(DE_RB_2005[loc,tw][,"sumMetricModel50Nghbr95l"], type="l", col="darkgrey", lty=2)
  legend("topright",legend = c("observed","sum-metric","95 % prediction band"), lty=c(1,1,2),
         col=c("darkgreen", "black", "darkgrey") )
  dev.off()
  
  DE_RB_2005@data$diffPM10 <- DE_RB_2005@data$sumMetricModel50Nghbr - DE_RB_2005@data$PM10

  stpl <- stplot(as(DE_RB_2005[,smplDays, "diffPM10"],"STFDF"),
                 col.regions=bpy.colors(5),
                 sp.layout = list("sp.polygons", DE_NUTS1), scales=list(draw=F), 
                 key.space="right", colorkey=T, cuts=c(-25,-15,-5,5,15,25),
                 main=NULL) #expression(paste("daily mean ","PM"[10]," concentration")))
  
  png("vignettes/figures/diffs_daily_means_PM10.png", width=9, height=6, "in", res=150)
  print(stpl)
  dev.off()
}