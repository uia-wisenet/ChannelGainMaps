classdef ForwardHybridKernelMachine_old
% An object of this class just contains the results of learning from a
% single training set, and provides a method to evaluate the estimated
% function for a query point.

properties        
    estimator HybridEstimator
    intercept
    coefficients
    m_trainFeatures_LF
    m_trainFeatures_LB
    funct_weights
    v_trainMetric
end

methods
    function v_estimate_hybrid = evaluate(obj, ...
            m_queryFeatures_LF, m_queryFeatures_LB, v_locationErrors)
        %input check
        assert(iscolumn(v_locationErrors));
        n_ues = length(v_locationErrors);
        assert(size(m_queryFeatures_LB, 2) == n_ues);
        assert(size(m_queryFeatures_LF, 2) == n_ues);
        assert(size(m_queryFeatures_LB, 1) == size(obj.m_trainFeatures_LB, 1));
        assert(size(m_queryFeatures_LF, 1) == size(obj.m_trainFeatures_LF, 1));
        
%         % obtain LF estimate
%         v_estimate_LF = obj.FKM_LF.evaluate(m_queryFeatures_LF);
%         % obtain LB estimate
%         v_estimate_LB = obj.FKM_LB.evaluate(m_queryFeatures_LB);
%         % obtain eval weights from interpolating the input errors
%         v_weights = obj.eval_weights(v_locationErrors);
%         v_estimate_hybrid = v_weights  .* v_estimate_LB + ...
%                          (1-v_weights) .* v_estimate_LF;
       

        v_estimate = obj.estimator.estimateGivenEstimatedDistancesAndLoc(... 
                    obj.coefficients, obj.m_trainFeatures_LF,...
                    obj.m_trainFeatures_LB, obj.v_trainMetric, ...
                    obj.funct_weights,...
                    m_queryFeatures_LF, m_queryFeatures_LB);
                % sugg: change method name to: evaluateEstimatedMap
                % [given  features]
    end
end
end