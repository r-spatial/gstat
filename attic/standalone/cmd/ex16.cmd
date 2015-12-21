#
# Multivariable indicator cosimulation 
#
data(i200): 'zinc.eas', x=1, y=2, v=3, max=20, radius=1000,
    I=200, sk_mean = 0.28;
data(i400): 'zinc.eas', x=1, y=2, v=3, max=20, radius=1000,
    I=400, sk_mean = 0.56;
data(i800): 'zinc.eas', x=1, y=2, v=3, max=20, radius=1000,
    I=800, sk_mean = 0.85;

# define an LMC:
variogram(i200): 0.0490637 Nug(0) + 0.182814 Exp(300);
variogram(i400): 0.0608225 Nug(0) + 0.21216 Exp(300);
variogram(i200, i400): 0 Nug() + 0.14806 Exp(300);
variogram(i800): 0.0550284 Nug(0) + 0.0842966 Exp(300);
variogram(i200, i800): 0 Nug() + 0.0525584 Exp(300);
variogram(i400, i800): 0 Nug() + 0.102852 Exp(300);

method: is;
mask: 'mask_map.tif';
# apply order corrections for cumulative indicators:
set order = 4;

predictions(i200): 'i200pr';
predictions(i400): 'i400pr';
predictions(i800): 'i800pr';
# uncomment next line to get 5 simulations:
# set nsim = 5;
