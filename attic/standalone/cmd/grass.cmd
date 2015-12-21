#
# Inverse distance interpolation on a mask map
#
#Was, for grass5: data(zinc): 'zincpts', v=1, I=200;
#is, now, for grass 6.0beta1:
data(zinc): 'zincpts', v=1;
mask:       'mask_map.tif'; # the prediction locations
predictions(zinc): 'id_pr'; # result map
