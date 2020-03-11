classdef ForwardLBKernelMachine
    % An object of this class just contains the results of learning from a
    % single training set, and provides a method to evaluate the estimated
    % function for a query point.
    
    properties        
        estimator LocationBasedEstimator

        coefficients
        intercept % was set as mean of v_trainMetric
        m_trainFeatures 
    end
    
    methods
        function v_estimate = evaluate(obj, m_queryFeatures)
            
            v_estimate = obj.estimator.estimate(obj.coefficients, obj. m_trainFeatures,...
                                   m_queryFeatures, obj.intercept);
            
            
        end
    end
end
