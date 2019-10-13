% wave 7 and uniform 4 for wavelet analysis

load('makemovies_9_19_w18_u15.mat')
wave_plot = wave{7};
uniform_plot = uniform{4};

[sl_cell,keepme,cells2trash] = order_shoreline_bwbound(lake);