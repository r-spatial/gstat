#
# Universal kriging, using one continuous and
# two binary variables.
#
data(ln_zinc): 'zincmap.eas', x=1, y=2, v=3, log,
  X=-1&4&5&6;
# -1: no default intercept (col. 5 and 6 form an intercept)
# use global kriging: local kriging would lead to a singularity
# the variogram of e is:
variogram(ln_zinc): 0.0698 Nug(0) + 0.147 Sph(709);
  # mask maps holding the independent variable values
  # at prediction locations, their order corresponding
  # that of the X-columns:
mask: 'sqrtdist.tif', 'part_a.tif', 'part_b.tif';
predictions(ln_zinc): 'lzn_vkpr'; 
variances(ln_zinc):   'lzn_vkvr';
