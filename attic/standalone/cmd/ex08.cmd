#
# Change of support: local ordinary block kriging on a mask
#
data(ln_zinc): 'zinc.eas', x=1, y=2, v=3, log,
  min=20, max=40, radius=1000;
variogram(ln_zinc): 0.0554 Nug(0) + 0.581 Sph(900);
mask: 'mask_map.tif';
predictions(ln_zinc): 'lzn_okbp';
variances(ln_zinc):   'lzn_okbv';
blocksize: dx=40, dy=40; # define block dimensions
