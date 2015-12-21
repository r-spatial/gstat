#
# Multivariable kriging: ordinary local cokriging of two variables
#
data(ln_zinc): 'zincmap.eas', x=1, y=2, v=3, log,
  min=20, max=40, radius=1000;
data(sq_dist): 'zincmap.eas', x=1, y=2, v=4,
  min=20, max=40, radius=1000;
variogram(ln_zinc): 0.0554 Nug(0) + 0.581 Sph(900);
variogram(sq_dist): 0.0001 Nug(0) + 0.0631 Sph(900);
variogram(ln_zinc, sq_dist): 0 Nug(0)-0.156 Sph(900);
  # NOTE: the 0 Nug(0)'s are added to make gstat recognize 
  # the Linear Model of Coregionalization
mask: 'mask_map.tif';
predictions(ln_zinc): 'lzn_ckpr';
variances(ln_zinc):   'lzn_ckvr';
predictions(sq_dist): 'sqd_ckpr';
variances(sq_dist):   'sqd_ckvr';
# the next map holds the prediction error covariances:
covariances(sq_dist,ln_zinc): 'znsqdcov';
