function H_D = calculateImpulse_Resp2(obj, v_xy_rx)
    n_sources=length(obj.xt_loc);
    H_D=zeros(n_sources,obj.maxSamplesPerPilot);
    for ind_source=1:n_sources
        v_xy_tx = [obj.xt_loc(ind_source) obj.yt_loc(ind_source)];
        Rx = calculatePaths(obj, v_xy_tx, v_xy_rx);
        
        % This section comes from the street canyon generator
        
        distanceNray=[RxTx.dist  firstOrderRef  secondOrderReflec]; %  Rx.distFirstRefl  TxSecondRef2Rx.dist
        delays=distanceNray/obj.c;
        
        powerDelayProfile=10*log10([abs(Rx.LosRssi) abs(Rx.reflecjRssi) abs(Rx.SeconReflWallJRSSI)]);   % abs(Rx.ReflecRssi) abs(Rx.SecondRefRSSI)
        
        %
        alpha_coeff = sqrt(db2pow(powerDelayProfile));
        alpha_coeff2 = sqrt([abs(Rx.LosRssi) abs(Rx.reflecjRssi) abs(Rx.SeconReflWallJRSSI)]);
        norm(alpha_coeff - alpha_coeff2);
        keyboard
        error 'TODO: add the transmit power of each source!'
        
        [ h_D, ~] = obj.digitalImpulseResponse(alpha_coeff, delays);
        
        %%
        H_D(ind_source, 1:length(h_D))=h_D; % Get the impulse responses from all sources
    end
    
end