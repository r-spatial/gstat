#
# Gaussian simulation, conditional upon data
# (local neighbourhoods, simple and ordinary kriging)
#
data(ln_zinc): 'zinc.eas', x=1, y=2, v=3, log,
  sk_mean=5.9, max=20;
variogram(ln_zinc): 0.0554 Nug(0) + 0.581 Sph(900);
mask: 'mask_map.tif';
method: gs;
predictions(ln_zinc): 'lzn_cspr';
set nsim=10;
