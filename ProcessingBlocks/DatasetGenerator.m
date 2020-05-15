classdef DatasetGenerator
    properties
        generator        ChannelGainGenerator
        locEstimator     LocationEstimator_tdoa
        featureExtractor
        sampler
        
        b_syntheticLocError = 0
        b_syntheticLocEstimate = 0;
        std_syntheticLocationNoise = 5;
        snr_locErrors = 50;

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
            t3_receivedPilotSignals = t3_nps + ...
                obj.sampler.pilotNoiseSTD/sqrt(2) *(...
                randn(size(t3_nps)) + 1i*randn(size(t3_nps)) );
            
            % extract features
            str_out.m_features_LF = obj.featureExtractor.locFreeExtract(...
                t3_receivedPilotSignals); % Center of mass of xcorr
            
            % channel gain for each pair
            v_channels = zeros(n_p,1);
            ltc2 = LoopTimeControl(n_p);
            for i_pair = 1:n_p
                xy_transmitter = m_locations(m_pairs(i_pair,1), :);
                xy_receiver    = m_locations(m_pairs(i_pair,2), :);
                v_channels(i_pair) = obj.generator.calculateCGBetween(...
                    xy_transmitter, xy_receiver);
                ltc2.go(i_pair);
            end
            str_out.v_channelGains = v_channels;
            
            if obj.b_syntheticLocEstimate
                warning('just simulating location algorithm by adding random noise')
                m_estimatedLocations = m_locations ...
                    + obj.std_syntheticLocationNoise*randn(size(m_locations));
                assert (obj.b_syntheticLocError==1)
            else
                disp 'extracting LocBased features (TDoA)...'
                str_out.m_features_LB = obj.featureExtractor.locBasedExtract(...
                t3_receivedPilotSignals);  % Estim. time diff. of arrival

                disp 'estimating locations at train set...'
                [m_estimatedLocations, v_locUncertainties] = ...
                    obj.locEstimator.estimateLocations(str_out.m_features_LB);
            end
            if obj.b_syntheticLocError
                locErrors_noiseless = vecnorm(m_locations - m_estimatedLocations, 2, 2);
                v_locUncertainties = locErrors_noiseless ...
                    + randn(size(locErrors_noiseless))...
                    * std(m_estimatedLocations(:)/obj.snr_locErrors);
            end
            str_out.m_estimatedLocations = m_estimatedLocations;
            str_out.v_locUncertainties = v_locUncertainties;
        end
    end
end   