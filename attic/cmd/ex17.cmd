#
# global coordinate polynomial trend surfaces
# trend orders 0-3.
#
data(zinc.0): 'zinc.eas', x=1, y=2, v=3, d=0, log;
data(zinc.1): 'zinc.eas', x=1, y=2, v=3, d=1, log;
data(zinc.2): 'zinc.eas', x=1, y=2, v=3, d=2, log;
data(zinc.3): 'zinc.eas', x=1, y=2, v=3, d=3, log;
mask:                'mask_map.tif'; 
# predict block averages for very small blocks:
blocksize: dx=1, dy=1;
# variances apply to mean values,
# not for single observations
predictions(zinc.0): 'lzn_tr0';
variances  (zinc.0): 'lzn_vr0';
predictions(zinc.1): 'lzn_tr1';
variances  (zinc.1): 'lzn_vr1';
predictions(zinc.2): 'lzn_tr2';
variances  (zinc.2): 'lzn_vr2';
predictions(zinc.3): 'lzn_tr3';
variances  (zinc.3): 'lzn_vr3';
