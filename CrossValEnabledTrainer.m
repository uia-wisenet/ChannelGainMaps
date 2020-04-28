classdef CrossValEnabledTrainer
    properties
        n_folds = 5; % for cross validation
        b_cv_inParallel = 0;
 
    end
    methods
        
        [v_estimate, meanErrOnEvalFeat] = evaluate(obj, m_queryFeatures)
        
        function [m_crossValScores, best_lambda, best_sigma] = ...
                crossValidate(obj, m_features, v_c2m_metric, ...
                v_lambdas, v_sigmas)
            cv_obj = cvpartition(v_c2m_metric, 'KFold', obj.n_folds);
            [m_lambdas, m_sigmas] = ndgrid(v_lambdas, v_sigmas);
            m_crossValScores = zeros(size(m_lambdas));
            assert( not(ishandle(obj.estimator)));
            if obj.b_cv_inParallel, error 'not implemented yet', end
            ltc = LoopTimeControl(numel(m_lambdas));
            for ii = 1:numel(m_lambdas)                                    
                my_estimator = obj.estimator;
                my_estimator.regularizationParameter = m_lambdas(ii);
                my_estimator.kernel = @(x, y) ...
                    exp(-norms(x-y, 2, 1).^2/(m_sigmas(ii)^2));
                v_scoresFolds = zeros(obj.n_folds, 1);
                for i_fold = 1:obj.n_folds
                    obj_now = obj;
                    obj_now.estimator = my_estimator;
                    fkm = obj_now.train(...
                        m_features( :, cv_obj.training(i_fold), : ), ...
                        v_c2m_metric(  cv_obj.training(i_fold)) );
                    v_c2m_test = fkm.evaluate(  ...
                        m_features( :, cv_obj.test(i_fold), : )  );
                    v_scoresFolds(i_fold) = mean(...
                    (v_c2m_test(:) - v_c2m_metric(cv_obj.test(i_fold))).^2);
                end
                m_crossValScores(ii) = mean(v_scoresFolds)/...
                    mean((v_c2m_metric - mean(v_c2m_metric)).^2);
                ltc.go(ii)
            end
            [min_value, best_index] = min(m_crossValScores(:));
            best_lambda = m_lambdas(best_index);
            best_sigma  = m_sigmas (best_index);
        end
        
    end
end

