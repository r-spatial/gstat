#
# Stratified ordinary kriging (within-categorie ordinary kriging)
#
data(zinc.at.0): 'zincat0.eas', x=1, y=2, v=3, log,
  min=20, max=40, radius=1000;  # where part_a = 0
data(zinc.at.1): 'zincat1.eas', x=1, y=2, v=3, log,
  min=20, max=40, radius=1000;  # where part_a = 1
variogram(zinc.at.0): 0.0654 Nug(0) + 0.548 Sph(900);
variogram(zinc.at.1): 0.716 Sph(900);
# the mask map is 0 for zinc.at.0 locations, 1 for zinc.at.1
mask: 'part_a';
# stratified mode: one map holds predictions for all vars:
predictions: 'lzn_stpr';
# another the prediction variances for all vars:
variances:   'lzn_stvr';
