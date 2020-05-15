classdef Simulator2
    properties
        generator MapGenerator = PrecomputedMapGenerator();
        locFreeEstimator LocationFreeEstimator
        locBasedEstimator LocationBasedEstimator
        hybridEstimator HybridEstimator = HybridEstimator();
        locEstimator LocationEstimator_tdoa = DiverseLocationEstimator('');
        featureExtractor
        sampler % the only property we use from sampler is pilotNoiseSTD
        
        PCAenable = 0;
        
        b_trainLocFree = 1;
        b_trainLocBased = 1;
        b_trainHybrid = 0;
        
        b_cv = 0;
        b_cvLambdas_hybrid = 0;
        
        b_syntheticLocError = 0
        b_syntheticLocEstimate = 0;
        std_syntheticLocationNoise = 9;
        snr_locErrors = 20;
    end
    
    methods
        function [str_NMSE, str_mapEstimates, channelForPairs_tracemap] = ...
                simulate ( obj, m_train_pairs, m_tracemap_pairs, ...
                m_grid_x, m_grid_y )
            %inputs: 
            %    m_train_pairs:   N_train x 2 matrix containing grid indices
            %    m_tracemap_pairs: N_test x 2 matrix containing grid indices
            
            assert(size(m_train_pairs, 2)    == 2);
            assert(size(m_tracemap_pairs, 2) == 2);
            assert(all(size(m_grid_x)==size(m_grid_y)));
            
            assert(obj.PCAenable==0, 'tpca not implemented yet in Simulator2')
            
            N_pair_train    = size(m_train_pairs,    1);
            N_pair_tracemap = size(m_tracemap_pairs, 1);
            
            source_loc=[obj.generator.xt_loc; obj.generator.yt_loc];
            n_sources = size(source_loc,2);
            
            %% generate pilot signals 
            %  for all points in the grid
            
            disp('Precomputing pilot signals...')
            % noiseless pilot signals
            % 4-way tensor: t4_nps(source, pilotSample, ind_x, ind_y)
            t4_nps = zeros([n_sources, ... 
                obj.generator.maxSamplesPerPilot,size(m_grid_x)]);
            ltc = LoopTimeControl(numel(m_grid_x));
            for ind_gridPoint = 1:numel(m_grid_x)
                [ind_row, ind_col] = ind2sub(size(m_grid_x), ind_gridPoint);
                t4_nps(:,:,ind_row, ind_col) = ...
                    obj.generator.calculateImpulse_Resp(...
                    [m_grid_x(ind_gridPoint), m_grid_y(ind_gridPoint)]);               
                ltc.go(ind_gridPoint);
            end
            % Add noise:            
            t4_receivedPilotSignals = t4_nps + ... % taken from Sampler
                obj.sampler.pilotNoiseSTD/sqrt(2) *( ...
                randn(size(t4_nps)) + 1i*randn(size(t4_nps)) ); 
            % Reshape tensor to be used by featureExtractor:
            v_size = size(t4_receivedPilotSignals);
            t3_reshaped_receivedPS = reshape(t4_receivedPilotSignals, ...
                [v_size(1:2), prod(v_size(3:4))]);
            
            %% Extract features for locFree and locBased algs
            
            disp('Extracting features...')
            t_allLocFreeFeatures  = obj.featureExtractor.locFreeExtract(...
                t3_reshaped_receivedPS); %Center of mass of xcorr
            t_allLocBasedFeatures = obj.featureExtractor.locBasedExtract(...
                t3_reshaped_receivedPS); %Time differences of arrival
            
            %% Compute channel gains
            % this section can be refactored by defining a static method
            channelForPairsTr_noiseless=zeros(N_pair_train,1);
            n_sp_CoM  = size(t_allLocFreeFeatures,  1); % no. of source pairs
            n_sp_TDoA = size(t_allLocBasedFeatures, 1); % no. of source pairs
            t3_trueLoc_train  = zeros([2, size(m_train_pairs)]);
            locFreeFeatures_train  = zeros([n_sp_CoM,  size(m_train_pairs)]);
            locBasedFeatures_train = zeros([n_sp_TDoA, size(m_train_pairs)]);
            
            % channel gains at training pairs
            disp('Computing channel gains (for training)...')
            ltc = LoopTimeControl(N_pair_train);
            for ind_pair = 1:N_pair_train
                ind_tx = m_train_pairs(ind_pair, 1);
                ind_rx = m_train_pairs(ind_pair, 2);
                xy_transmitter = [m_grid_x(ind_tx), m_grid_y(ind_tx)];
                xy_receiver    = [m_grid_x(ind_rx), m_grid_y(ind_rx)];
                channelForPairsTr_noiseless(ind_pair) = ...
                    obj.generator.calculateCGBetween(...
                    xy_transmitter, xy_receiver);
                locFreeFeatures_train(:, ind_pair, :) = ...
                    t_allLocFreeFeatures(:, [ind_tx ind_rx]);
                locBasedFeatures_train(:, ind_pair, :) = ...
                    t_allLocBasedFeatures(:, [ind_tx ind_rx]);
                t3_trueLoc_train(:, ind_pair, 1) = xy_transmitter';
                t3_trueLoc_train(:, ind_pair, 2) = xy_receiver';
                ltc.go(ind_pair);
            end
            % add noise:
            channelForPairsTr = channelForPairsTr_noiseless + ...
                obj.sampler.powerNoiseSTD* ...
                randn(size(channelForPairsTr_noiseless));
            
            % compute channel gains at test pairs
            channelForPairs_tracemap = zeros(1, N_pair_tracemap); 
            t3_trueLoc_tracemap   = zeros([2, size(m_tracemap_pairs)]);
            locFreeFeat_tracemap  = zeros([n_sp_CoM, size(m_tracemap_pairs)]);
            locBasedFeat_tracemap = zeros([n_sp_TDoA, size(m_tracemap_pairs)]);
            disp('Computing channel gains (noiseless, for map tracing)...')
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
                t3_trueLoc_tracemap(:, ind_pair, 1) = xy_transmitter';
                t3_trueLoc_tracemap(:, ind_pair, 2) = xy_receiver';
                ltc.go(ind_pair);
            end
            
            str_NMSE = struct();
            str_mapEstimates = struct();
                     
            %%
            % here, the original simulator had some lines defining which
            % features to consider. This is skipped here.        
            
            %% Training and testing: locFree
            
            if obj.b_trainLocFree
                disp ('Training LocFree...')
                my_locFreeTrainer = LocFreeTrainer;
                my_locFreeTrainer.estimator = obj.locFreeEstimator;
                %
                v_lambdas_toTryLF = logspace(-4, -2, 14);    %TODO: move to properties
                v_sigmas_toTryLF  = linspace(85, 200, 15);   %TODO: move to properties
                if obj.b_cv
                    disp ('Cross-validation...')
                    [m_crossValScoresLF, best_lambdaLF, best_sigmaLF] = ...
                        my_locFreeTrainer.crossValidate(...
                        locFreeFeatures_train, channelForPairsTr, ...
                        v_lambdas_toTryLF, v_sigmas_toTryLF);
                    my_locFreeTrainer.estimator.regularizationParameter = best_lambdaLF;
                    my_locFreeTrainer.estimator.kernel =  @(x, y) ...
                        exp(-norms(x-y, 2, 1).^2/(best_sigmaLF^2));
                    figure(999); clf
                    mesh(v_sigmas_toTryLF, v_lambdas_toTryLF, m_crossValScoresLF); 
                    ax = gca;
                    ax.YAxis.Scale = 'log';
                    xlabel \sigma
                    ylabel \lambda
                    title 'LOCATION FREE'
                end
                %
                [locFreeFKM, probMissingFeat] = my_locFreeTrainer.train(...
                    locFreeFeatures_train, channelForPairsTr); %#ok<ASGLU>
                %
                disp('Evaluating LocFree at validation set...')
                [str_mapEstimates.locFree, meanErrOnEvalFeat] = locFreeFKM.evaluate(...
                    locFreeFeat_tracemap); %#ok<ASGLU>
                disp('Done');
                SqErr_locFree = Simulator.compute_SqErr(...
                    channelForPairs_tracemap, str_mapEstimates.locFree);
                str_NMSE.locFree = SqErr_locFree/sum((channelForPairs_tracemap...
                    -mean(channelForPairs_tracemap)).^2);                
            else
                str_NMSE.locFree = 1;
                str_mapEstimates.locFree = zeros(size(channelForPairs_tracemap));
            end
            
            
            %% If locationBased or Hybrid is used, 
            % call the locationEstimator
            if obj.b_trainLocBased || obj.b_trainHybrid
                if obj.b_syntheticLocEstimate
                    warning('just simulating location algorithm by adding random noise')
                    
                    t3_estim_loc_train    = t3_trueLoc_train ...
                        + obj.std_syntheticLocationNoise*randn(size(t3_trueLoc_train));
                    t3_estim_loc_tracemap = t3_trueLoc_tracemap ...
                        + obj.std_syntheticLocationNoise*randn(size(t3_trueLoc_train));                   
                else
                    disp 'estimating locations at train set...'
                    m_estim_loc_tx = obj.locEstimator.estimateLocations(...
                        locBasedFeatures_train(:,:,1));
                    m_estim_loc_rx = obj.locEstimator.estimateLocations(...
                        locBasedFeatures_train(:,:,2));
                    % t3_estimatedLocation is a 3-way tensor:
                    % (space_dim_index, sample_index, rx_or_tx)
                    t3_estim_loc_train = cat(3, m_estim_loc_tx, m_estim_loc_rx);
                    
                    disp('estimating locations at tracemap (test) set...')
                    estim_loc_tm_tx = obj.locEstimator.estimateLocations(...
                        locBasedFeat_tracemap(:,:,1)); %estimated loc (tracemap)
                    estim_loc_tm_rx = obj.locEstimator.estimateLocations(...
                        locBasedFeat_tracemap(:,:,2));
                    % t3_estim_loc_tracemap is a 3-way tensor:
                    % (space_dim_index, sample_index, rx_or_tx)
                    t3_estim_loc_tracemap = cat(3, estim_loc_tm_tx, estim_loc_tm_rx);
                end
                if obj.b_syntheticLocError
                    locErrors_noiseless_train = ...
                        squeeze(vecnorm(t3_trueLoc_train - t3_estim_loc_train));
                    locErrors_train = locErrors_noiseless_train ...
                        + randn(size(locErrors_noiseless_train))...
                        * std(locErrors_noiseless_train(:))/obj.snr_locErrors;

                    locErrors_noiseless_val = squeeze(...
                        vecnorm(t3_trueLoc_tracemap - t3_estim_loc_tracemap));
                    locErrors_val = locErrors_noiseless_val ...
                        + randn(size(locErrors_noiseless_val))...
                        * std(locErrors_noiseless_train(:))/obj.snr_locErrors;
                else
                    warning ('location uncertainty not implemented yet...')
                    locErrors_train = [];
                    locErrors_val   = [];
                end
                    

            end
            
            %% Train and test location Based estimator
            
            if obj.b_trainLocBased
                disp ('Training LocBased...')
                my_locBasedTrainer = LocBasedTrainer;
                my_locBasedTrainer.estimator = obj.locBasedEstimator;
                %
                v_lambdas_toTryLB = logspace(-5, -2,  15);   %TODO: move to properties
                v_sigmas_toTryLB  = linspace(20, 100, 12);   %TODO: move to properties
                if obj.b_cv
                    disp ('Cross-validation...')
                    [m_crossValScoresLB, best_lambdaLB, best_sigmaLB] = ...
                        my_locBasedTrainer.crossValidate(...
                        t3_estim_loc_train, channelForPairsTr, ...
                        v_lambdas_toTryLB, v_sigmas_toTryLB);
                    my_locBasedTrainer.estimator.regularizationParameter = best_lambdaLB;
                    my_locBasedTrainer.estimator.kernel =  @(x, y) ...
                        exp(-norms(x-y, 2, 1).^2/(best_sigmaLB^2));
                    figure(998); clf
                    mesh(v_sigmas_toTryLB, v_lambdas_toTryLB, m_crossValScoresLB); 
                    ax = gca;
                    ax.YAxis.Scale = 'log';
                    xlabel \sigma
                    ylabel \lambda
                    title 'LOCATION BASED'
                end % save por_siacaso
                
                %
                locBasedFKM = my_locBasedTrainer.train(...
                    t3_estim_loc_train,channelForPairsTr);
                
                str_mapEstimates.locBased = locBasedFKM.evaluate(t3_estim_loc_tracemap);
                str_NMSE.locBased = obj.compute_NMSE(...
                    channelForPairs_tracemap, str_mapEstimates.locBased);
            else
                str_NMSE.locBased = 1;
                str_mapEstimates.locBased = zeros(size(channelForPairs_tracemap));
            end
            
            %% Train and test hybrid estimator
            if obj.b_trainHybrid
                
                my_hybridTrainer = HybridTrainer;
                my_hybridTrainer.hybridEstimator = obj.hybridEstimator;
                
                % TODO: using the true location to compute the locErrors is "cheating"
                % we must do it using only the positioning pilots
                
                if obj.b_cv
                    disp(['Hybrid model inherits sigma from '...
                        'the pure loc-Free and the pure locBased,']);
                    disp ('and sets the lambdas with an extra round of crossval');
                    my_hybridTrainer.hybridEstimator.h_kernelLB = @(x, y) ...
                        exp(-norms(x-y, 2, 1).^2/(best_sigmaLB^2));
                    
                    my_hybridTrainer.hybridEstimator.h_kernelLF = @(x, y) ...
                        exp(-norms(x-y, 2, 1).^2/(best_sigmaLF^2));
                    
                    if obj.b_cvLambdas_hybrid
                        v_lambdas_toTryLB_hyb = best_lambdaLB *...
                            logspace(-1, 1,  5);  %TODO: set based on properties
                        v_lambdas_toTryLF_hyb = best_lambdaLF *...
                            logspace(-1, 1,  5);  %TODO: set based on properties
                        [m_crossValScoresHybrid, best_lambdaLF_hyb, ...
                            best_lambdaLB_hyb] = my_hybridTrainer.crossValidateLambdas(...
                            locFreeFeatures_train, t3_estim_loc_train, locErrors_train, channelForPairsTr, ...
                            v_lambdas_toTryLF_hyb, v_lambdas_toTryLB_hyb);
                        my_hybridTrainer.hybridEstimator.regularizationParameterLF...
                            = best_lambdaLF_hyb;
                        my_hybridTrainer.hybridEstimator.regularizationParameterLB...
                            = best_lambdaLB_hyb;
                        %
                        figure(997); clf
                        mesh(v_lambdas_toTryLB_hyb, v_lambdas_toTryLF_hyb, m_crossValScoresHybrid);
                        ax = gca;
                        ax.XAxis.Scale = 'log';
                        ax.YAxis.Scale = 'log';
                        xlabel '\lambda LocBased'
                        ylabel '\lambda LocFree'
                        title 'HYBRID'
                    else
                        disp(['Hybrid model inherits lambda and sigma from '...
                            'the pure loc-Free and the pure locBased']);
                        my_hybridTrainer.hybridEstimator.regularizationParameterLB = best_lambdaLB;
                        my_hybridTrainer.hybridEstimator.regularizationParameterLF = best_lambdaLF;
                    end
                
                end
                
                disp ('Training Hybrid...')
                
                my_hybridTrainer.hybridEstimator.b_debugPlots = 1;
                hybridFKM = my_hybridTrainer.train(...
                locFreeFeatures_train, t3_estim_loc_train, locErrors_train, channelForPairsTr);
                
                disp('Evaluating Hybrid at validation set...')
                
                str_mapEstimates.hybrid = hybridFKM.evaluate(...
                    locFreeFeat_tracemap, t3_estim_loc_tracemap, locErrors_val);
                %
                str_mapEstimates.locFree_fromHybrid = hybridFKM.fkm_lf.evaluate(...
                    locFreeFeat_tracemap);
                str_mapEstimates.locBased_fromHybrid = hybridFKM.fkm_lb.evaluate(...
                    t3_estim_loc_tracemap);
                %
                disp('Done');
                str_NMSE.hybrid =             obj.compute_NMSE(...
                    channelForPairs_tracemap, str_mapEstimates.hybrid);
                str_NMSE.locFree_fromHybrid = obj.compute_NMSE(...
                    channelForPairs_tracemap, str_mapEstimates.locFree_fromHybrid);
                str_NMSE.locBased_fromHybrid = obj.compute_NMSE(...
                    channelForPairs_tracemap, str_mapEstimates.locBased_fromHybrid);
            else
                str_NMSE.hybrid              = 1;
                str_NMSE.locFree_fromHybrid  = 1;
                str_NMSE.locBased_fromHybrid = 1;
                [str_mapEstimates.hybrid, str_mapEstimates.locBased_fromHybrid,  ...
                    str_mapEstimates.locFree_fromHybrid  ] ...
                    = deal(zeros(size(channelForPairs_tracemap)));                
            end
        end
    end
    methods(Static)
        function NMSE = compute_NMSE(trueMap, estimatedMap)
            sqErr = Simulator.compute_SqErr(...
                    trueMap(:), estimatedMap(:));
                NMSE = sqErr/sum((trueMap - mean(trueMap)).^2);
        end
    end
end
