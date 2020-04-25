function channel_gain = calculateCGBetween2(obj, xt, xr)

Rx = calculatePaths(obj, xt, xr);

%% Reflection & Line Of Signt Propagation Map

channel_LOS=Rx.LosCG;
channel_first_refl= obj.reflectExaggerationFac*sum(Rx.reflecjCG);
channel_second_refl= obj.reflectExaggerationFac*sum(Rx.SeconReflWallJCG);

channel_gain=10*log10(abs(channel_LOS+channel_first_refl+channel_second_refl));