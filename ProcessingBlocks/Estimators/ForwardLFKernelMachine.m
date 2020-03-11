classdef ForwardLFKernelMachine
    % An object of this class just contains the results of learning from a
    % single training set, and provides a method to evaluate the estimated
    % function for a query point.
    
    properties        
        estimator LocationFreeEstimator

        coefficients
        combin_sources
        orthBasis
        meanFeature
        m_featCov 
        completedmeasurements
        
        m_trainFeatures
        v_trainMetric
        intercept % was set as mean of v_trainMetric
    end
    
    methods
        function [meanErrOnEvalFeat,v_estimate] = ...
                evaluate(obj, m_queryFeatures)
            
            [meanErrOnEvalFeat,v_estimate] = obj.estimator.estimateGivenEstimatedDistances(... % sugg: change method name to: evaluateEstimatedMap[given features]
                        obj.coefficients, obj.m_trainFeatures, obj.v_trainMetric,...
                        obj.completedmeasurements, m_queryFeatures,...
                        obj.combin_sources,obj.orthBasis, obj.intercept, ...
                        obj.meanFeature, obj.m_featCov);
            
            
        end
    end
end

