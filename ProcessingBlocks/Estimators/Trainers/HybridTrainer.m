classdef HybridTrainer
    % Trains a kernel machine for hybrid (LocF + LocB) Cartography.
    
    properties
        hybridEstimator HybridEstimator
        n_folds = 5;
        b_cv_inParallel = 0;
    end
    
    methods
        function fkm_hybrid = train(obj, m_features_LF, m_features_LB,...
                m_locErrors, v_c2m_metric)
            
            [v_coefficients_LF, v_coefficients_LB, intercept, wcu] ...
                = obj.hybridEstimator.train(m_features_LF, ...
                m_features_LB, m_locErrors, v_c2m_metric);
            
            
            lf = ForwardLFKernelMachine;
            %lines adapted from LocFreeTrainer
            lf.estimator = LocationFreeEstimator;
            lf.estimator.kernel = obj.hybridEstimator.h_kernelLF;
            lf.estimator.regularizationParameter = ...
                obj.hybridEstimator.regularizationParameterLF;
            % some lines defining properties necessary for missing feature
            % estimation can be inserted here
            lf.coefficients = v_coefficients_LF;
            lf.m_trainFeatures = m_features_LF;
            lf.v_trainMetric = v_c2m_metric;
            lf.intercept = intercept;
            
            lb = ForwardLBKernelMachine;
            lb.estimator = LocationBasedEstimator;
            lb.estimator.kernel = obj.hybridEstimator.h_kernelLB;
            lb.estimator.regularizationParameter = ...
                obj.hybridEstimator.regularizationParameterLB;
            % some lines defining properties necessary for missing feature
            % estimation can be inserted here
            lb.coefficients = v_coefficients_LB;
            lb.m_trainFeatures = m_features_LB;
            lb.intercept = intercept;
            
            fkm_hybrid = ForwardHybridKernelMachine;
            fkm_hybrid.fkm_lf = lf;
            fkm_hybrid.fkm_lb = lb;
            fkm_hybrid.wcu = wcu;
            
        end
        
        function [m_crossValScores, best_lambda_LF, best_lambda_LB] = ...
                crossValidateLambdas(obj, m_features_LF, m_features_LB, ...
                m_locErrors, v_c2m_metric, v_lambdas_LF, v_lambdas_LB)
            cv_obj = cvpartition(v_c2m_metric, 'KFold', obj.n_folds);
            [m_lambdas_LF, m_lambdas_LB] = ndgrid(v_lambdas_LF, v_lambdas_LB);
            m_crossValScores = zeros(size(m_lambdas_LF));
            assert( not(ishandle(obj.hybridEstimator)));
            if obj.b_cv_inParallel, error 'not implemented yet', end
            ltc = LoopTimeControl(numel(m_crossValScores));
            for ii = 1:numel(m_crossValScores)                                    
                my_estimator = obj.hybridEstimator;
                my_estimator.regularizationParameterLF = m_lambdas_LF(ii);
                my_estimator.regularizationParameterLB = m_lambdas_LB(ii);
%                 my_estimator.kernel = @(x, y) ...
%                     exp(-norms(x-y, 2, 1).^2/(m_sigmas(ii)^2));
                v_scoresFolds = zeros(obj.n_folds, 1);
                for i_fold = 1:obj.n_folds
                    obj_now = obj;
                    obj_now.hybridEstimator = my_estimator;
                    fkm = obj_now.train(...
                        m_features_LF( :, cv_obj.training(i_fold), : ), ...
                        m_features_LB( :, cv_obj.training(i_fold), : ), ...
                        m_locErrors  ( cv_obj.training(i_fold),: ),    ...
                        v_c2m_metric(  cv_obj.training(i_fold))       );
                    v_c2m_test = fkm.evaluate(  ...
                        m_features_LF( :, cv_obj.test(i_fold), : ), ...
                        m_features_LB( :, cv_obj.test(i_fold), : ), ...
                        m_locErrors     ( cv_obj.test(i_fold),: )    );
                    v_scoresFolds(i_fold) = mean(  (v_c2m_test(:) - ...
                         v_c2m_metric(cv_obj.test(i_fold))).^2  );
                end
                m_crossValScores(ii) = mean(v_scoresFolds)/...
                    mean((v_c2m_metric - mean(v_c2m_metric)).^2);
                ltc.go(ii)
            end
            [min_value, best_index] = min(m_crossValScores(:));
            best_lambda_LF = m_lambdas_LF (best_index);
            best_lambda_LB = m_lambdas_LB (best_index);
        end
        
    end
end