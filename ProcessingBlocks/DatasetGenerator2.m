classdef DatasetGenerator2
    properties
        generator        ChannelGainGenerator
        locEstimator     LocationEstimator_tdoa
        featureExtractor
        sampler
        
        b_syntheticLocError = 0
        b_syntheticLocEstimate = 0;
        std_syntheticLocationNoise = 5;
        snr_locErrors = 50;
        
        n_realizations = 1;
        b_inParallel = 0;

    end
    
    % Note: the input to "generate" could be in the form:
    % t_grid_xy is a 3-way array of size s1 x s2 x 2.
    
    methods
        function [str_out] = generate (obj, m_locations, m_pairs)
            % m_locations is an Nl x 2 matrix, where each row contains the x
            % and y coordinates of each location to have into account
            % m_pairs is an Np x 2 matrx, where each row contains
            % the indices (respective to the ordering in m_locations) for 
            % the transmitter and receiver
            
            [n_l, the_dim] = size(m_locations); % n_l = number of locations
            assert(the_dim==2)
            assert(size(m_pairs,2) == 2);
            n_p = size(m_pairs, 1); % n_p = number of pairs
            assert(  isequal(m_pairs, round(m_pairs)) ...
                && all(m_pairs(:) >= 1   ) ...
                && all(m_pairs(:) <= n_l )  );
            
            str_out = struct();
            str_out.m_locations = m_locations;
            str_out.m_pairs     = m_pairs;
            
            n_s = size(obj.generator.xt_loc, 2); % number of location sources
            
            % generate pilot signals
            % 3-way array: t3_nps(sourceIndex, timeIndex, locationIndex)
            % t3_nps = noiseless pilot signals
            t3_nps = zeros(n_s, obj.generator.maxSamplesPerPilot, n_l); 
            ltc = LoopTimeControl(n_l);
            for i_location = 1:n_l
                t3_nps(:,:,i_location) = ...
                    obj.generator.calculateImpulse_Resp( ...
                    m_locations(i_location,:) );
                ltc.go(i_location)
            end            
            % add noise
            t4_nps = repmat(t3_nps, [1 1 1 obj.n_realizations]);
            t4_receivedPilotSignals = t4_nps + ...
                obj.sampler.pilotNoiseSTD/sqrt(2) *(...
                randn(size(t4_nps)) + 1i*randn(size(t4_nps)) );
            
            % extract features
            for i_r = obj.n_realizations:-1:1
                str_out.am_features_LF(:,:,i_r) = obj.featureExtractor.locFreeExtract(...
                    t4_receivedPilotSignals(:,:,:,i_r)); % Center of mass of xcorr
            end
            
            % channel gain for each pair
            v_gains = zeros(n_p,1);
            ltc2 = LoopTimeControl(n_p);
            if obj.b_inParallel
                disp ('Calculating channel gains in parallel')
                parfor i_pair = 1:n_p
                    xy_transmitter = m_locations(m_pairs(i_pair,1), :);
                    xy_receiver    = m_locations(m_pairs(i_pair,2), :);
                    v_gains(i_pair) = obj.generator.calculateCGBetween(...
                        xy_transmitter, xy_receiver);
                end
            else
                for i_pair = 1:n_p
                    xy_transmitter = m_locations(m_pairs(i_pair,1), :);
                    xy_receiver    = m_locations(m_pairs(i_pair,2), :);
                    v_gains(i_pair) = obj.generator.calculateCGBetween(...
                        xy_transmitter, xy_receiver);
                    ltc2.go(i_pair);
                end
            end
            str_out.v_trueGains = v_gains;
            av_channelGains = repmat(v_gains, [1 1 obj.n_realizations]);
            str_out.av_noisy_channelGains = av_channelGains + ...
                obj.sampler.powerNoiseSTD*randn(size(av_channelGains));
            
            if obj.b_syntheticLocEstimate
                warning('just simulating location algorithm by adding random noise')
                am_noiseless_locations = repmat(m_locations, [1 1 obj.n_realizations]);
                am_estimatedLocations = am_noiseless_locations ...
                    + obj.std_syntheticLocationNoise*randn(size(am_noiseless_locations));
                assert (obj.b_syntheticLocError==1)
            else
                disp 'extracting LocBased features (TDoA)...'
                for i_r = 1:obj.n_realizations
                    str_out.am_features_LB(:,:,i_r) = ...
                        obj.featureExtractor.locBasedExtract(...
                        t4_receivedPilotSignals(:,:,:,i_r));  
                       %^ Estim. time diff. of arrival
                end
                disp 'estimating locations at train set...'
                am_estimatedLocations = zeros(n_l, 2, obj.n_realizations);
                av_locUncertainties = zeros(n_l, 1, obj.n_realizations);
                if obj.b_inParallel
                    parfor i_r = 1:obj.n_realizations
                        [am_estimatedLocations(:,:,i_r), av_locUncertainties(:,:,i_r)] = ...
                            obj.locEstimator.estimateLocations(str_out.am_features_LB(:,:,i_r));
                    end   
                else
                    for i_r = 1:obj.n_realizations                        
                        [am_estimatedLocations(:,:,i_r), av_locUncertainties(:,:,i_r)] = ...
                            obj.locEstimator.estimateLocations(str_out.am_features_LB(:,:,i_r));
                    end
                end
            end
            if obj.b_syntheticLocError
                v_locErrors_noiseless = vecnorm(m_locations - am_estimatedLocations, 2, 2);
                av_locUncertainties = v_locErrors_noiseless ...
                    + randn(size(v_locErrors_noiseless))...
                    * std(am_estimatedLocations(:)/obj.snr_locErrors);
            end
            str_out.am_estimatedLocations = am_estimatedLocations;
            str_out.av_locUncertainties   = av_locUncertainties;
        end
    end
end   