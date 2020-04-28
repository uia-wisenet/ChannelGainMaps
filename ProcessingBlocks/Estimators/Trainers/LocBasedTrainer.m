classdef LocBasedTrainer < CrossValEnabledTrainer
    % Trains a kernel machine for Location Based Cartography.
    
    properties
        estimator LocationBasedEstimator
    end
    
    methods
        function fkmLB_out = train(obj, t_features, v_c2m_metric)
            
            [locBasedCoefficients, ~] = obj.estimator.train(t_features,v_c2m_metric);
            
            fkmLB_out = ForwardLBKernelMachine;
            fkmLB_out.estimator = obj.estimator;
            fkmLB_out.coefficients = locBasedCoefficients;
            fkmLB_out.m_trainFeatures = t_features;
            fkmLB_out.intercept = mean(v_c2m_metric);
        end        
    end
end

