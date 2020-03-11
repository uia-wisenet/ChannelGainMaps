classdef ForwardHybridKernelMachine
    % An object of this class just contains the results of learning from a
    % single training set, and provides a method to evaluate the estimated
    % function for a query point.
    
    properties        
        estimator HybridEstimator
        coefficients
        m_trainFeatures_LF
        m_trainFeatures_LB
        funct_weights
        v_trainMetric
    end
    
    methods
        function v_estimate = ...
                evaluate(obj, m_queryFeatures_LF, m_queryFeatures_LB)
            
            v_estimate = obj.estimator.estimateGivenEstimatedDistancesAndLoc(... % sugg: change method name to: evaluateEstimatedMap[given features]
                        obj.coefficients, obj.m_trainFeatures_LF,obj.m_trainFeatures_LB, obj.v_trainMetric, obj.funct_weights,...
                        m_queryFeatures_LF,m_queryFeatures_LB);
            
            
        end
    end
end