\name{DE_RB_2005}
\alias{DE_RB_2005}
\docType{data}
\title{
Spatio-temporal data set with rural background PM10 concentrations in Germany 2005
}
\description{
Spatio-temporal data set with rural background PM10 concentrations in Germany 2005 (airbase v6). 
}
\usage{data("DE_RB_2005")}
\format{
  The format is:
Formal class 'STSDF' [package "spacetime"] with 5 slots
  ..@ data   :'data.frame':	23230 obs. of  2 variables:
  .. ..$ PM10   : num [1:23230] 16.7 31.7 5 22.4 26.8 ...
  .. ..$ logPM10: num [1:23230] 2.82 3.46 1.61 3.11 3.29 ...
  ..@ index  : int [1:23230, 1:2] 1 2 3 4 5 6 7 8 9 10 ...
  ..@ sp     :Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
  .. .. ..@ data       :'data.frame':	69 obs. of  9 variables:
  .. .. .. ..$ station_altitude     : int [1:69] 8 3 700 15 35 50 343 339 45 45 ...
  .. .. .. ..$ station_european_code: Factor w/ 7965 levels "AD0942A","AD0944A",..: 1991 1648 1367 2350 1113 1098 1437 2043 1741 1998 ...
  .. .. .. ..$ country_iso_code     : Factor w/ 39 levels "AD","AL","AT",..: 10 10 10 10 10 10 10 10 10 10 ...
  .. .. .. ..$ station_start_date   : Factor w/ 2409 levels "1900-01-01","1951-04-01",..: 152 1184 1577 1132 744 328 1202 1555 1148 407 ...
  .. .. .. ..$ station_end_date     : Factor w/ 864 levels "","1975-02-06",..: 1 1 1 579 1 1 1 1 1 1 ...
  .. .. .. ..$ type_of_station      : Factor w/ 5 levels "","Background",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ station_type_of_area : Factor w/ 4 levels "rural","suburban",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ street_type          : Factor w/ 5 levels "","Canyon street: L/H < 1.5",..: 4 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ annual_mean_PM10     : num [1:69] 20.9 21.8 16.5 20.3 23.3 ...
  .. .. ..@ coords.nrs : num(0) 
  .. .. ..@ coords     : num [1:69, 1:2] 538709 545414 665711 551796 815738 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:69] "DESH001" "DENI063" "DEBY109" "DEUB038" ...
  .. .. .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  .. .. ..@ bbox       : num [1:2, 1:2] 307809 5295752 907375 6086661
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  .. .. .. .. ..$ : chr [1:2] "min" "max"
  .. .. ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot
  .. .. .. .. ..@ projargs: chr "+init=epsg:32632 +proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  ..@ time   :An ?xts? object on 2005-01-01/2005-12-31 containing:
  Data: int [1:365, 1] 5115 5116 5117 5118 5119 5120 5121 5122 5123 5124 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "..1"
  Indexed by objects of class: [POSIXct,POSIXt] TZ: GMT
  xts Attributes:  
 NULL
  ..@ endTime: POSIXct[1:365], format: "2005-01-02" "2005-01-03" "2005-01-04" "2005-01-05" ...
}
\source{
EEA, airbase v6
}
\examples{
data(DE_RB_2005)
str(DE_RB_2005)
}
\keyword{datasets}
