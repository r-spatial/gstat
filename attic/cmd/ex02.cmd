#
# Two variables with (initial estimates of) variograms,
# start the variogram modelling user interface
#
data(zinc): 'zinc.eas', x=1, y=2, v=3;  
data(ln_zinc): 'zinc.eas', x=1, y=2, v=3, log;

variogram(zinc): 10000 Nug() + 140000 Sph(800);
variogram(ln_zinc): 1 Nug() + 1 Sph(800);
