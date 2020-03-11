classdef LocFreeTrainer
    % Trains a kernel machine for Location Free Cartography.
    
    properties
        estimator LocationFreeEstimator
    end
    
    methods
        function [fkm_out, avgNumberOfMissingFeat] = train(obj, m_features, v_c2m_metric)
            
            [locFreeCoefficients,~,ues_many_misses,avgNumberOfMissingFeat,...
                combin_sources,orthBasis,meanFeature,...
                featCovMat,completedmeasurements] ...
                = obj.estimator.train(m_features, v_c2m_metric);
            
            fkm_out = ForwardLFKernelMachine;
            fkm_out.estimator = obj.estimator;
            
            fkm_out.coefficients = locFreeCoefficients;
            fkm_out.combin_sources = combin_sources;
            fkm_out.orthBasis = orthBasis;
            fkm_out.meanFeature = meanFeature;
            fkm_out.m_featCov = featCovMat;
            fkm_out.completedmeasurements = completedmeasurements;
            
            fkm_out.m_trainFeatures = m_features;
            fkm_out.v_trainMetric    = v_c2m_metric;
            fkm_out.intercept       = mean(v_c2m_metric);
        end
    end
end

