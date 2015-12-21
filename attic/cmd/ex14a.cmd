#
# Stratified universal kriging (within-categorie universal kriging)
#
data(zinc.at.0): 'zincat0.eas', x=1, y=2, v=3, X=4, log,
  min=20, max=40, radius=1000;  # where part_a = 0
data(zinc.at.1): 'zincat1.eas', x=1, y=2, v=3, X=4, log,
  min=20, max=40, radius=1000;  # where part_a = 1

# residual variograms:
variogram(zinc.at.0): 0.096572 Nug(0) + 0.226367 Sph(1069.33);
variogram(zinc.at.1): 0.115766 Sph(237.257);

mask: 'part_a', # 0 for zinc.at.0 locations, 1 for zinc.at.1 locs.
  'sqrtdist', # predictor values corresp. to col. 4 for zinc.at.0
  'sqrtdist'; # predictor values corresp. to col. 4 for zinc.at.1
predictions: 'lzn_stup';
variances:   'lzn_stuv';
