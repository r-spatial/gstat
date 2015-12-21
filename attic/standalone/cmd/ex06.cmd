#
# Unconditional Gaussian simulation on a mask
# (local neigbourhoods, simple kriging)
#
# defines empty variable:
data(ln_zn_dummy): dummy, sk_mean=5.9, max=20;
variogram(ln_zn_dummy): 0.0554 Nug(0) + 0.581 Sph(900);
mask:  'mask_map.tif';
method:  gs; # Gaussian simulation instead of kriging
predictions(ln_zn_dummy): 'lzn_uspr';
set nsim = 10;
