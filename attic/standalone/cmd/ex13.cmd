#
# Local universal kriging, using one continuous variable
#
data(ln_zinc): 'zincmap.eas', x=1, y=2, v=3, log, 
  X=4, # sqrtdist values at data locations
  min=20, max=40, radius=1000;  # apply model locally
# the variogram should be that of the residual:
variogram(ln_zinc): 0.0674 Nug(0) + 0.149 Sph(700);
mask: 'sqrtdist';  # sqrtdist values at prediction locations
predictions(ln_zinc): 'lzn_ukpr'; 
variances(ln_zinc):   'lzn_ukvr';
