classdef Simulator3
    properties
        locFreeEstimator LocationFreeEstimator
        locBasedEstimator LocationBasedEstimator
        hybridEstimator HybridEstimator = HybridEstimator();
        PCAenable = 0;
        
        v_lambdas_toTryLF = []
        v_sigmas_toTryLF  = []
        
        v_lambdas_toTryLB = []
        v_sigmas_toTryLB  = []
        
        b_trainLocFree = 1;
        b_trainLocBased = 1;
        b_trainHybrid = 0;
        
        b_cvLambdas_hybrid = 0;
        
        b_augmentTraining = 0;
        
        b_meshCVSurfaces = 0;
        
        b_cv_inParallel = 0;
    end
    
    methods
        function [str_NMSE, str_mapEstimates, v_channelGains_test,...
                str_NMSE_num, str_NMSE_den] = ...
                simulate ( obj, str_dataset, v_train_pairs_in, v_test_pairs_in)
            %inputs: 
            %    m_train_pairs:   N_train x 2 matrix containing grid indices
            %    m_test_pairs: N_test x 2 matrix containing grid indices
            
            m_allPairs = str_dataset.m_pairs;
            assert(size(m_allPairs,2) == 2); 
            n_allPairs = size(m_allPairs,1);
            
            assert( iscolumn(v_train_pairs_in) && ...
                isequal(v_train_pairs_in, round(v_train_pairs_in)) && ...
                max(v_train_pairs_in <= n_allPairs) )
            assert( iscolumn(v_test_pairs_in) && ...
                isequal(v_test_pairs_in, round(v_test_pairs_in)) && ...
                max(v_test_pairs_in <= n_allPairs) )
            
            assert(obj.PCAenable==0, 'tpca not implemented yet in Simulator3')
            

            % Augment the training set with the symmetric pairs
            
            n_pair_train    = numel(v_train_pairs_in);
            n_pair_test     = numel(v_test_pairs_in);
            
            source_loc=[str_dataset.datasetGen.generator.xt_loc; ...
                        str_dataset.datasetGen.generator.yt_loc];
            n_sources = size(source_loc,2);
            
            str_NMSE = struct();
            str_mapEstimates = struct();
            
            m_pairs_augmented = ... % adding symmetric pairs:
                [m_allPairs ; fliplr(m_allPairs)];
            v_channelGains_augmented = ...
                [str_dataset.v_noisy_channelGains; str_dataset.v_noisy_channelGains];        
            v_trainAug_pairs = [v_train_pairs_in; v_train_pairs_in + n_allPairs];
            
            if isequal(...
                    sortrows(m_pairs_augmented(v_trainAug_pairs,:)), ...
                    unique(m_pairs_augmented(v_trainAug_pairs,:), 'rows') )
                v_trainAug_noDupes = v_trainAug_pairs;
            else
                warning('Some training instances are duplicated. Removing...')
                b_aux = true(length(v_trainAug_pairs),1);
                for ii = 1:length(v_train_pairs_in)
                    if b_aux(ii)
                        my_pair = m_pairs_augmented(v_trainAug_pairs(ii),:);
                        b_found = m_pairs_augmented(:,1)==my_pair(1) & ...
                            m_pairs_augmented(:,2)==my_pair(2);
                        b_found(v_trainAug_pairs(ii)) = 0;
                        for iii = 1:length(b_found)
                            if b_found(iii)
                                b_aux(v_trainAug_pairs==iii) = 0;
                            end
                        end
                    end
                end
                v_trainAug_noDupes = v_trainAug_pairs(b_aux);
            end
                       
            v_channelGains_train = v_channelGains_augmented(v_trainAug_noDupes);
            % TODO: maybe remove duplicates from the set of map/test
            % pairs, for efficiency (probably not critical)
            v_channelGains_test  = str_dataset.v_trueGains(v_test_pairs_in);
            
            %%
            % here, the original simulator had some lines defining which
            % features to consider. This is skipped here.  
                      
            %% Training and testing: locFree
            
            t3_locFreeFeatures_train = obj.build_set(str_dataset.m_features_LF, ...
                v_trainAug_noDupes, m_pairs_augmented);
            t3_locFreeFeatures_test =  obj.build_set(str_dataset.m_features_LF, ...
                v_test_pairs_in,  m_allPairs);
            
            if obj.b_trainLocFree
                disp ('Training LocFree...')
                my_locFreeTrainer = LocFreeTrainer;
                my_locFreeTrainer.b_cv_inParallel = obj.b_cv_inParallel;
                my_locFreeTrainer.estimator = obj.locFreeEstimator;
                %
                
                n_values_toTry = numel(obj.v_lambdas_toTryLF)*numel(obj.v_sigmas_toTryLF);
                if n_values_toTry > 1
                    disp ('Cross-validation...')
                    [m_crossValScoresLF, best_lambdaLF, best_sigmaLF] = ...
                        my_locFreeTrainer.crossValidate(...
                        t3_locFreeFeatures_train, v_channelGains_train, ...
                        obj.v_lambdas_toTryLF, obj.v_sigmas_toTryLF);
                    my_locFreeTrainer.estimator.regularizationParameter = best_lambdaLF;
                    my_locFreeTrainer.estimator.kernel =  @(x, y) ...
                        exp(-vecnorm(x-y, 2, 1).^2/(best_sigmaLF^2));
                    if obj.b_meshCVSurfaces
                        figure(999); clf
                        mesh(obj.v_sigmas_toTryLF, obj.v_lambdas_toTryLF, m_crossValScoresLF);
                        ax = gca;
                        ax.YAxis.Scale = 'log';
                        xlabel \sigma
                        ylabel \lambda
                        title 'LOCATION FREE'
                    end
                elseif n_values_toTry == 1
                    best_lambdaLF = obj.v_lambdas_toTryLF;
                    best_sigmaLF = obj.v_sigmas_toTryLF;
                else
                    warning 'v_lambdas_toTryLB or v_sigmas_toTryLB properties should not be empty'
                    keyboard
                end
                my_locFreeTrainer.estimator.regularizationParameter = best_lambdaLF;
                my_locFreeTrainer.estimator.kernel =  @(x, y) ...
                    exp(-vecnorm(x-y, 2, 1).^2/(best_sigmaLF^2)); 
                
                [locFreeFKM, probMissingFeat] = my_locFreeTrainer.train(...
                    t3_locFreeFeatures_train, v_channelGains_train); %#ok<ASGLU>
                %
                disp('Evaluating LocFree at validation set...')
                [str_mapEstimates.locFree, meanErrOnEvalFeat] = locFreeFKM.evaluate(...
                    t3_locFreeFeatures_test); %#ok<ASGLU>
                disp('Done');
                [str_NMSE.locFree, str_NMSE_num.locFree, ...
                    str_NMSE_den.locFree] = obj.compute_NMSE(...
                    v_channelGains_test, str_mapEstimates.locFree);
            else
                str_NMSE.locFree = 1;
                str_mapEstimates.locFree = zeros(size(v_channelGains_test));
            end
                        
            %% Train and test location Based estimator
            
            t3_estimatedLocations_train = obj.build_set(str_dataset.m_estimatedLocations', ...
                v_trainAug_noDupes, m_pairs_augmented);
            t3_estimatedLocations_test =  obj.build_set(str_dataset.m_estimatedLocations', ...
                v_test_pairs_in,  m_allPairs);
            
            if obj.b_trainLocBased
                disp ('Training LocBased...')
                my_locBasedTrainer = LocBasedTrainer;
                my_locBasedTrainer.b_cv_inParallel = obj.b_cv_inParallel;
                my_locBasedTrainer.estimator = obj.locBasedEstimator;
                %
                n_values_toTry = numel(obj.v_lambdas_toTryLB)*numel(obj.v_sigmas_toTryLB);
                if n_values_toTry > 1
                    disp ('Cross-validation...')
                    [m_crossValScoresLB, best_lambdaLB, best_sigmaLB] = ...
                        my_locBasedTrainer.crossValidate(...
                        t3_estimatedLocations_train, v_channelGains_train, ...
                        obj.v_lambdas_toTryLB, obj.v_sigmas_toTryLB);
                    if obj.b_meshCVSurfaces
                        figure(998); clf
                        mesh(obj.v_sigmas_toTryLB, obj.v_lambdas_toTryLB, m_crossValScoresLB);
                        ax = gca;
                        ax.YAxis.Scale = 'log';
                        xlabel \sigma
                        ylabel \lambda
                        title 'LOCATION BASED'
                    end
                elseif n_values_toTry == 1
                    best_lambdaLB = obj.v_lambdas_toTryLB;
                    best_sigmaLB  = obj.v_sigmas_toTryLB;
                else
                    warning 'v_lambdas_toTryLB or v_sigmas_toTryLB properties should not be empty'
                    keyboard
                end
                my_locBasedTrainer.estimator.regularizationParameter = best_lambdaLB;
                my_locBasedTrainer.estimator.kernel =  @(x, y) ...
                    exp(-vecnorm(x-y, 2, 1).^2/(best_sigmaLB^2));
                
                locBasedFKM = my_locBasedTrainer.train(...
                    t3_estimatedLocations_train, v_channelGains_train);
                
                disp('Evaluating LocBased at validation set...')
                str_mapEstimates.locBased = locBasedFKM.evaluate(t3_estimatedLocations_test);
                [str_NMSE.locBased, str_NMSE_num.locBased, str_NMSE_den.locBased] = ...
                    obj.compute_NMSE(...
                    v_channelGains_test, str_mapEstimates.locBased);
            else
                str_NMSE.locBased = 1;
                str_mapEstimates.locBased = zeros(size(v_channelGains_test));
            end
            
            %% Train and test hybrid estimator
 
            m_locUncertainties_train = squeeze(obj.build_set(...
                str_dataset.v_locUncertainties', ...
                v_trainAug_noDupes, m_pairs_augmented));
            m_locUncertainties_test =  squeeze(obj.build_set(...
                str_dataset.v_locUncertainties', ...
                v_test_pairs_in,  m_allPairs));
            
            if obj.b_trainHybrid
                
                my_hybridTrainer = HybridTrainer;
                my_hybridTrainer.hybridEstimator = obj.hybridEstimator;
                
                % TODO: using the true location to compute the locErrors is "cheating"
                % we must do it using only the positioning pilots
                
                disp(['Hybrid model inherits sigma from '...
                    'the pure loc-Free and the pure locBased,']);
                disp ('and sets the lambdas with an extra round of crossval');
                my_hybridTrainer.hybridEstimator.h_kernelLB = @(x, y) ...
                    exp(-vecnorm(x-y, 2, 1).^2/(best_sigmaLB^2));
                
                my_hybridTrainer.hybridEstimator.h_kernelLF = @(x, y) ...
                    exp(-vecnorm(x-y, 2, 1).^2/(best_sigmaLF^2));
                                
                if obj.b_cvLambdas_hybrid
                    v_lambdas_toTryLB_hyb = best_lambdaLB *...
                        logspace(-1, 1,  5);  %TODO: set based on properties
                    v_lambdas_toTryLF_hyb = best_lambdaLF *...
                        logspace(-1, 1,  5);  %TODO: set based on properties
                    [m_crossValScoresHybrid, best_lambdaLF_hyb, ...
                        best_lambdaLB_hyb] = my_hybridTrainer.crossValidateLambdas(...
                        t3_locFreeFeatures_train, t3_estimatedLocations_train, ...
                        m_locUncertainties_train, v_channelGains_train, ...
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
                
                
                disp ('Training Hybrid...')
                
%                 my_hybridTrainer.hybridEstimator.b_debugPlots = 1;
                hybridFKM = my_hybridTrainer.train(...
                t3_locFreeFeatures_train, t3_estimatedLocations_train, ...
                m_locUncertainties_train, v_channelGains_train);
                
                disp('Evaluating Hybrid at validation set...')
                
                str_mapEstimates.hybrid = hybridFKM.evaluate(...
                    t3_locFreeFeatures_test, t3_estimatedLocations_test, ...
                    m_locUncertainties_test);
                %
                str_mapEstimates.locFree_fromHybrid = hybridFKM.fkm_lf.evaluate(...
                    t3_locFreeFeatures_test);
                str_mapEstimates.locBased_fromHybrid = hybridFKM.fkm_lb.evaluate(...
                    t3_estimatedLocations_test);
                str_mapEstimates.locBased_fromHybrid_corrected = ...
                    hybridFKM.fkm_lb.intercept +...
                    hybridFKM.wcu.evalWeights([0 0])*( +...
                    hybridFKM.fkm_lb.evaluate(t3_estimatedLocations_test)...
                    -hybridFKM.fkm_lb.intercept);
                %
                disp('Done');
                [str_NMSE.hybrid, str_NMSE_num.hybrid, str_NMSE_den.hybrid]...
                    =             obj.compute_NMSE(...
                    v_channelGains_test, str_mapEstimates.hybrid);
                [str_NMSE.locFree_fromHybrid, str_NMSE_num.locFree_fromHybrid,...
                    str_NMSE_den.locFree_fromHybrid] = obj.compute_NMSE(...
                    v_channelGains_test, str_mapEstimates.locFree_fromHybrid);
                [str_NMSE.locBased_fromHybrid, str_NMSE_num.locBased_fromHybrid, ...
                    str_NMSE_den.locBased_fromHybrid] = obj.compute_NMSE(...
                    v_channelGains_test, str_mapEstimates.locBased_fromHybrid);
                [str_NMSE.locBased_fromHybrid_corrected, ...
                    str_NMSE_num.locBased_fromHybrid_corrected,...
                    str_NMSE_den.locBased_fromHybrid_corrected] = obj.compute_NMSE(...
                    v_channelGains_test, str_mapEstimates.locBased_fromHybrid_corrected);
            else
                str_NMSE.hybrid              = 1;
                str_NMSE.locFree_fromHybrid  = 1;
                str_NMSE.locBased_fromHybrid = 1;
                str_NMSE.locBased_fromHybrid_corrected = 1;
                [str_mapEstimates.hybrid, str_mapEstimates.locBased_fromHybrid,  ...
                    str_mapEstimates.locFree_fromHybrid  ] ...
                    = deal(zeros(size(v_channelGains_test)));                
            end
        end
    end
    methods(Static)
        function [NMSE, num, den] = compute_NMSE(trueMap, estimatedMap)
            sqErr = Simulator.compute_SqErr(...
                    trueMap(:), estimatedMap(:));               
                num = sqErr;
                den = sum((trueMap - mean(trueMap)).^2);
                NMSE = num/den;
        end
        
        function [t_features_out] = build_set(m_features_in, v_pairs_in, m_allPairs)
            n_pair_train = numel(v_pairs_in);            
            n_feat_CoM = size(m_features_in, 1);
            t_features_out = zeros(n_feat_CoM, n_pair_train, 2);
            for i_pair = 1:n_pair_train
                ind_tx = m_allPairs(v_pairs_in(i_pair), 1);
                ind_rx = m_allPairs(v_pairs_in(i_pair), 2);
                t_features_out(:, i_pair, :) = ...
                    m_features_in(:, [ind_tx, ind_rx]);
            end
        end
        
    end
end
