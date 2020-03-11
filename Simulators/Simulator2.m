classdef Simulator2
    properties
        generator MapGenerator = PrecomputedMapGenerator();
        locFreeEstimator LocationFreeEstimator
        locBasedEstimator LocationBasedEstimator
        hybridEstimator HybridEstimator
        locEstimator LocationEstimator
        featureExtractor
        sampler
    end
    
    methods
        function [NMSE_locFree, NMSE_locBased, NMSE_hybrid, ...
                locFreeMapEstimate,locBasedMapEstimate,...
                hybridMapEstimate,channelForPairs_tracemap] = ...
                simulate ( obj, m_train_pairs, m_tracemap_pairs, ...
                m_grid_x, m_grid_y )
            %inputs: 
            %    m_train_pairs:   N_train x 2 matrix containing grid indices
            %    m_tracemap_pairs: N_test x 2 matrix containing grid indices
            
            assert(size(m_train_pairs, 2)    == 2);
            assert(size(m_tracemap_pairs, 2) == 2);
            assert(all(size(m_grid_x)==size(m_grid_y)));
            
            N_pair_train    = size(m_train_pairs,    1);
            N_pair_tracemap = size(m_tracemap_pairs, 1);
            
            source_loc=[obj.generator.xt_loc; obj.generator.yt_loc];
            n_sources = size(source_loc,2);
            
            %% generate pilot signals and extract features 
            %  for all points in the grid
            
            disp('Precomputing pilot signals...')
            t4_nps = zeros([n_sources, ... % noiseless pilot signals
                obj.generator.maxSamplesPerPilot,size(m_grid_x)]);
            ltc = LoopTimeControl(numel(m_grid_x));
            for ind_gridPoint = 1:numel(m_grid_x)
                [ind_row, ind_col] = ind2sub(size(m_grid_x), ind_gridPoint);
                t4_nps(:,:,ind_row, ind_col) = ...
                    obj.generator.calculateImpulse_Resp(...
                    [m_grid_x(ind_gridPoint), m_grid_y(ind_gridPoint)]);
                ltc.go(ind_gridPoint);
            end
            
            disp('Extracting features...')
            t4_receivedPilotSignals = t4_nps + ... % taken from Sampler
                obj.sampler.pilotNoiseSTD/sqrt(2) *( ...
                randn(size(t4_nps)) + 1i*randn(size(t4_nps)) ); 
            
            v_size = size(t4_receivedPilotSignals);
            t3_reshaped_receivedPS = reshape(t4_receivedPilotSignals, ...
                [v_size(1:2), prod(v_size(3:4))]);

            t_allLocFreeFeatures  = obj.featureExtractor.locFreeExtract(...
                t3_reshaped_receivedPS);
            t_allLocBasedFeatures = obj.featureExtractor.locBasedExtract(...
                t3_reshaped_receivedPS); %Time differences of arrival
            
            %% this section can be refactored by defining a static method
            channelForPairsTr=zeros(1,N_pair_train);
            n_sp_CoM  = size(t_allLocFreeFeatures,  1); % no. of source pairs
            n_sp_TDoA = size(t_allLocBasedFeatures, 1);
            locFreeFeatures_train  = zeros([n_sp_CoM, size(m_train_pairs)]);
            locBasedFeatures_train = zeros([n_sp_TDoA, size(m_train_pairs)]);
            
            disp('Computing channel gains (for training)...')
            ltc = LoopTimeControl(N_pair_train);
            for ind_pair = 1:N_pair_train
                ind_tx = m_train_pairs(ind_pair, 1);
                ind_rx = m_train_pairs(ind_pair, 2);
                xy_transmitter = [m_grid_x(ind_tx), m_grid_y(ind_tx)];
                xy_receiver    = [m_grid_x(ind_rx), m_grid_y(ind_rx)];
                channelForPairsTr(ind_pair) = ...
                    obj.generator.calculateCGBetween(...
                    xy_transmitter, xy_receiver);
                locFreeFeatures_train(:, ind_pair, :) = ...
                    t_allLocFreeFeatures(:, [ind_tx ind_rx]);
                locBasedFeatures_train(:, ind_pair, :) = ...
                    t_allLocBasedFeatures(:, [ind_tx ind_rx]);
                ltc.go(ind_pair);
            end
            channelForPairs_tracemap = zeros(1, N_pair_tracemap);    
            locFreeFeat_tracemap  = zeros([n_sp_CoM, size(m_tracemap_pairs)]);
            locBasedFeat_tracemap = zeros([n_sp_TDoA, size(m_tracemap_pairs)]);
            disp('Computing channel gains (for map tracing)...')
            ltc = LoopTimeControl(N_pair_tracemap);
            for ind_pair = 1:N_pair_tracemap
                ind_tx = m_tracemap_pairs(ind_pair, 1);
                ind_rx = m_tracemap_pairs(ind_pair, 2);
                xy_transmitter = [m_grid_x(ind_tx), m_grid_y(ind_tx)];
                xy_receiver    = [m_grid_x(ind_rx), m_grid_y(ind_rx)];
                channelForPairs_tracemap(ind_pair) = ...
                    obj.generator.calculateCGBetween(...
                    xy_transmitter, xy_receiver);
                locFreeFeat_tracemap(:, ind_pair, :) = ...
                    t_allLocFreeFeatures(:, [ind_tx ind_rx]);
                locBasedFeat_tracemap(:, ind_pair, :) = ...
                    t_allLocBasedFeatures(:, [ind_tx ind_rx]);
                ltc.go(ind_pair);
            end
                     
            %%
            % here, the original simulator had some lines defining which
            % features to consider. This is skipped here.        
            
            disp ('Training LocFree...')
            my_locFreeTrainer = LocFreeTrainer;
            my_locFreeTrainer.estimator = obj.locFreeEstimator;
            %TODO: reshape?
            [locFreeFKM, probMissingFeat] = my_locFreeTrainer.train(...
                locFreeFeatures_train, channelForPairsTr);
            disp('Evaluating LocFree at validation set...')
            [meanErrOnEvalFeat,locFreeMapEstimate] = locFreeFKM.evaluate(...
                locFreeFeat_tracemap);
            disp('Done');
            
            disp ('Training LocBased...')
            my_locBasedTrainer = LocBasedTrainer;
            my_locBasedTrainer.estimator = obj.locBasedEstimator;
            estim_loc_tx =obj.locEstimator.estimateLocationIRWSRDLS(...
                locBasedFeatures_train(:,:,1)); %estimated locations
            estim_loc_rx = obj.locEstimator.estimateLocationIRWSRDLS(...
                locBasedFeatures_train(:,:,2));
            estim_loc_train = cat(3, estim_loc_tx, estim_loc_rx);
            locBasedFKM = my_locBasedTrainer.train(...
                estim_loc_train,channelForPairsTr);
            
            disp('Evaluating LocBased at validation set...')
            estim_loc_tm_tx = obj.locEstimator.estimateLocationIRWSRDLS(...
                locBasedFeat_tracemap(:,:,1)); %estimated loc (tracemap)
            estim_loc_tm_rx = obj.locEstimator.estimateLocationIRWSRDLS(...
                locBasedFeat_tracemap(:,:,2));
            estim_loc_tracemap = cat(3, estim_loc_tm_tx, estim_loc_tm_rx);
            locBasedMapEstimate = locBasedFKM.evaluate(estim_loc_tracemap);
            
             disp ('Training Hybrid...')
            my_hybridTrainer = HybridTrainer;
            my_hybridTrainer.estimator = obj.hybridEstimator;
            %TODO: reshape?
            hybridFKM = my_hybridTrainer.train(...
                locFreeFeatures_train, estim_loc_train, channelForPairsTr);
            disp('Evaluating Hybrid at validation set...')
            hybridMapEstimate = hybridFKM.evaluate(...
                locFreeFeat_tracemap, estim_loc_tracemap);
            disp('Done');

            %% compute error measures
            SqErr_locFree = Simulator.compute_SqErr(...
                channelForPairs_tracemap, locFreeMapEstimate);
            NMSE_locFree = SqErr_locFree/sum((channelForPairs_tracemap...
                -mean(channelForPairs_tracemap)).^2);
            SqErr_locBased = Simulator.compute_SqErr(...
                channelForPairs_tracemap, locBasedMapEstimate);
            NMSE_locBased = SqErr_locBased/sum((channelForPairs_tracemap...
                -mean(channelForPairs_tracemap)).^2);
            SqErr_hybrid = Simulator.compute_SqErr(...
                channelForPairs_tracemap, hybridMapEstimate);
            NMSE_hybrid = SqErr_hybrid/sum((channelForPairs_tracemap...
                -mean(channelForPairs_tracemap)).^2);
        end
    end
end
