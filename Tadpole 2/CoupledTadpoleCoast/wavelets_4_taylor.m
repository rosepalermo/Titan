% wave 7 and uniform 4 for wavelet analysis

load('makemovies_9_19_w18_u15.mat')
for i = 1:length(wave)
    w_l(i) = length(find(wave{i}));
end
for i = 1:17
    u_l(i) = length(find(uniform{i}));
end

for i = 1:length(wave)-1
    wld(i) = w_l(i+1)-w_l(1);
end
for i = 1:length(u_l)-1
    uld(i) = u_l(i+1)-u_l(1);
end

wave_plot = wave{14};
uniform_plot = uniform{4};
wave_plot = wave{1};
uniform_plot = uniform{1};

[sl_cell_wave,~,cells2trash_wave] = order_shoreline_bwbound(wave_plot);
[sl_cell_uniform,~,cells2trash_uniform] = order_shoreline_bwbound(uniform_plot);