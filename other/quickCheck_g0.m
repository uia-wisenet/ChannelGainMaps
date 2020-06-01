clear
gsinfo;
load savedResults/results_2515.mat

intercept = mean(trueGains(:));

v_zm_locBasedH = str_mapEstimates.locBased_fromHybrid' - intercept;
v_zm_trueGains = trueGains(:) - intercept;
v_g0 = linspace(0.1, 1, 100);
v_rms = zeros(size(v_g0));

for i_g = 1:100
    v_rms(i_g) = rms(v_zm_trueGains - v_g0(i_g)*v_zm_locBasedH);
end
figure;
plot(v_g0, v_rms)