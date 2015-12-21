#
# Local ordinary block kriging at non-gridded locations
#
data(zinc): 'zinc.eas', x=1, y=2, v=3, 
  min=20, max=40, radius=1000; # local neighbourhood
variogram(zinc): 2.42e+04 Nug(0) + 1.34e+05 Sph(800);
data(): 'locs.eas', x=1, y=2; # prediction locations
blocksize: dx=40, dy=40;      # 40 $\times$ 40 block averages
set output = 'zincok.out';    # ascii output file
