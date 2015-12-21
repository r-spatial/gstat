#
# Multiple kriging: prediction on more than one variable
# (ordinary kriging of two variables)
# (note that zinc_map.eas wass obtained through ex09.gst)
#
data(ln_zinc): 'zincmap.eas', x=1, y=2, v=3, log,
  min=20, max=40, radius=1000;
data(sq_dist): 'zincmap.eas', x=1, y=2, v=4,
  min=20, max=40, radius=1000;
variogram(ln_zinc): 0.0554 Nug(0) + 0.581 Sph(900);
variogram(sq_dist): 0.0631 Sph(900);
mask: 'mask_map.tif';
predictions(ln_zinc): 'lzn_okpr';
variances(ln_zinc):   'lzn_okvr';
predictions(sq_dist): 'sqd_okpr';
variances(sq_dist):   'sqd_okvr';
