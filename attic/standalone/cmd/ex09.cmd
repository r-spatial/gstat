#
# Obtain map values at data() locations
# (Point-map overlay)
#
data(a): dummy; # define n dummy data variable (n=n_masks/2)
data(b): dummy; 
data(): 'zinc.eas', x=1, y=2, v=3; # prediction locations
method: map;                       # mapvalues as `predictions'
masks: 'sqrtdist', 'part_a', 'part_b'; # the maps
set output = 'zincmap.eas';        # ascii output file.
