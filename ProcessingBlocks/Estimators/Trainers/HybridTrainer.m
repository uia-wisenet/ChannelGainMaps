classdef HybridTrainer
    % Trains a kernel machine for hybrid (LocF + LocB) Cartography.
    
    properties
        hybridEstimator HybridEstimator
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
        
        % TODO: create here a crossValidate method
        
    end
end