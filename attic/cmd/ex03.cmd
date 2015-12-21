#
# Inverse distance interpolation on a mask map
#
data(zinc): 'zinc.eas', x=1, y=2, v=3;
mask:       'mask_map.tif';     # the prediction locations
predictions(zinc): 'id_pr'; # result map
