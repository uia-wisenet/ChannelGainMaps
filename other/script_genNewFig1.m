clear
close all

gsinfo;
load savedResults/results_2515.mat
vs_titles = ["locF", "locB", ...
"MoE.locF", "MoE.locB", "MoE", "true"];

intercept = mean(trueGains(:));
g0 = 0.5;
v_zm_locBasedH = str_mapEstimates.locBased_fromHybrid' - intercept;
v_corrected_locBasedH = g0*v_zm_locBasedH + intercept;

str_dataset.animator.ch_colormap = 'jet';
str_dataset.animator.c_plotSource_options{1} = 'pk';
str_dataset.animator.c_plotSource_options{3} = 'r';
str_dataset.animator.c_plotSource_options{5} = 5;

figure;
[h_cb] = str_dataset.animator.storyBoard([...
str_mapEstimates.locFree'...
str_mapEstimates.locBased' ...
str_mapEstimates.locFree_fromHybrid' ...
v_corrected_locBasedH ...
str_mapEstimates.hybrid(:) ...
trueGains(:)], ...
vec(vs_titles), [2 5 9]);

h_cb.Position = [0.8333 0.1235 0.0119 0.7840];
h_cb.Limits   = [-80 -40];
