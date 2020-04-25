classdef BiasErrorSimulator
    properties
        generator MapGenerator = PrecomputedMapGenerator();
        locEstimator LocationEstimator_tdoa = WangLocationEstimator();
        featureExtractor
        sampler
    end
    
    methods
        function [m_featDiff, v_locErrors, m_trueTDoAs, m_allLocBasedFeatures,...
                m_estim_loc_tm_rx, m_trueLoc_all_points] = ...
                simulate ( obj, m_grid_x, m_grid_y )
            %inputs:
            % m_grid_x: matrix containing the x-coordinates of all grid points
            % m_grid_y: matrix containing the y-coordinates of all grid points
      
            assert(all(size(m_grid_x)==size(m_grid_y)));
            num_all_points = numel(m_grid_x);
            source_loc=[obj.generator.xt_loc; obj.generator.yt_loc];
            n_sources = size(source_loc,2);
            
            %% generate pilot signals and extract features
            %  for all points in the grid
            
            disp('Precomputing pilot signals...')
            m_trueLoc_all_points  = zeros([2, num_all_points]);
            t4_nps = zeros([n_sources, ... % noiseless pilot signals
                obj.generator.maxSamplesPerPilot,size(m_grid_x)]);
            ltc = LoopTimeControl(numel(m_grid_x));
            for ind_gridPoint = 1:numel(m_grid_x)
                m_trueLoc_all_points(:, ind_gridPoint) =...
                    [m_grid_x(ind_gridPoint), m_grid_y(ind_gridPoint)]';
                [ind_row, ind_col] = ind2sub(size(m_grid_x), ind_gridPoint);
                t4_nps(:,:,ind_row, ind_col) = ...
                    obj.generator.calculateImpulse_Resp(...
                    [m_grid_x(ind_gridPoint), m_grid_y(ind_gridPoint)]);
                ltc.go(ind_gridPoint);
            end
            
            disp('Extracting location-based features...')
            t4_receivedPilotSignals = t4_nps + ... % taken from Sampler
                obj.sampler.pilotNoiseSTD/sqrt(2) *( ...
                randn(size(t4_nps)) + 1i*randn(size(t4_nps)) );
            % Reshape tensor to be used by featureExtractor:
            v_size = size(t4_receivedPilotSignals);
            t3_reshaped_receivedPS = reshape(t4_receivedPilotSignals, ...
                [v_size(1:2), prod(v_size(3:4))]);
            
            m_allLocBasedFeatures = obj.featureExtractor.locBasedExtract(...
                t3_reshaped_receivedPS); %Time differences of arrival
            
            % obtain true TDoAS
            m_trueTDoAs = obj.featureExtractor.obtainRangeDiffer(...
                source_loc,m_trueLoc_all_points);   
            disp('Computing TDoA-based position estimates...')
            m_estim_loc_tm_rx = obj.locEstimator.estimateLocations(...
                m_allLocBasedFeatures);
           
            % compute localization errors and feature differences
            m_featDiff = m_trueTDoAs - m_allLocBasedFeatures;
            v_locErrors = vecnorm(m_trueLoc_all_points - m_estim_loc_tm_rx);
            disp('Done');
            
            
        end
    end
end
