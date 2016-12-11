close all
[x, fs] = audioread('../annotations/when_a_fire_starts_to_burn.wav');
sample_frame = 1:1024;
%% Zero crossing
tz = odf_zero();
tz.initialise(x, fs, 'zero crossings')
tz.compute_feature();
tz.plot_sample_intermediate_data(sample_frame)
%% Weighted phase deviation
twpd = odf_wpd();
twpd.initialise(x, fs, 'wpd')
twpd.compute_feature();
twpd.plot_sample_intermediate_data(sample_frame)
%% Energy flux
te = odf_energyflux();
te.initialise(x, fs, 'ef');
te.compute_feature();
te.plot_sample_intermediate_data(sample_frame)
%% Spectral flux
ts = odf_spectralflux();
ts.initialise(x, fs, 'sf');
ts.compute_feature();
ts.plot_sample_intermediate_data(sample_frame)
%% Complex spectral diff
tcsd = odf_complexspecdiff();
tcsd.initialise(x, fs, 'csd');
tcsd.compute_feature();
tcsd.plot_sample_intermediate_data(sample_frame)