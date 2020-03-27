classdef HybridTrainer
    % Trains a kernel machine for hybrid (LocF + LocB) Cartography.
    
    properties
        estimator HybridEstimator
    end
    
    methods
        function fkm_out = train(obj, m_features_LF, m_features_LB,...
                local_errors, v_c2m_metric)
            
            [hybridCoefficients,~, sig_ens] ...
                = obj.estimator.train(m_features_LF, m_features_LB,...
                local_errors, v_c2m_metric);
            
            fkm_out = ForwardHybridKernelMachine;
            fkm_out.estimator = obj.estimator;
            
            fkm_out.coefficients = hybridCoefficients;
            fkm_out.m_trainFeatures_LF = m_features_LF;
            fkm_out.m_trainFeatures_LB = m_features_LB;
            fkm_out.funct_weights    = sig_ens;
            fkm_out.v_trainMetric    = v_c2m_metric;
        end
    end
end