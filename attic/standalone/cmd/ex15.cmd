#
# Local linear model, using one continuous variable
#
data(ln_zinc): 'zincmap.eas', x=1, y=2, v=3, X=4, log,
  min=20, max=40, radius=1000;  # apply linear model locally
# no variogram definition: assume residual to be IID.
mask: 'sqrtdist';
predictions(ln_zinc): 'lzn_trpr';
# prediction variance for point locations:
variances(ln_zinc):   'lzn_trvr';
