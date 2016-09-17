This includes 4 files. You probably want to run 'forest_reverb' first

-forest_reverb.m: Top level script to calculate reverberation in a forest or hard cylinder. Adjust N_trees and N_sc_max to adjust number of trees and scatterings. With the default values, it takes about 200 s to run on a Core2Duo 6400 @ 2.13 GHz. Increase the size to fill whatever time you have.

-scatter_impulse.m: Function to calculate the impulse response of one scattering.

-scatter_cyl.m: Function to calculate one scattered wave

-mpf.m: optional function to put impulses into minimum phase. You probably don't want to do this unless you have a reason.

Travis Wiens
travis.mlfx@nutaksas.com