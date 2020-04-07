classdef RobustLocationExperiments < ExperimentFunctionSet
    % Experimenting with Robust optimization to solve the
    % minimax problem for the source localization
    
    properties
        N = 8;        % number of sources
        
        % upper bound on the NLOS range error
        rho = 25; % meters
        
        % std of the noise on range difference
        sigma_range_d_noise = 10; % meters
        
        % std of the Gaussian distributed sources and sensor
        sigma_source_spread = 100; %meters
        sigma_sensor_spread = 10;
                
        lambda = 0
        
    end
    
    methods
        function [m_sourceLocs, v_x_sensor] = generateSources_Sensor(obj)
            m_sourceLocs = obj.sigma_source_spread*randn(2, obj.N);
            v_x_sensor   = obj.sigma_sensor_spread*randn(2, 1);
        end
        
        function v_d = rangeDifferences(obj, v_x,  m_s, rho)
            v_d = zeros(obj.N-1, 1);
            v_r0 = norm(v_x - m_s(:, 1)) + rho*rand;
            for n = 2:obj.N
                v_d(n-1) = norm(v_x - m_s(:, n)) + rho*rand - v_r0 ... 
                    + obj.sigma_range_d_noise*randn;
            end
        end
               
        % first experiment with ADMMLocationEstimator
        function F = experiment_201(obj, niter)
            rng(2)
            m_s      = obj.generateSources;
            v_x_true = randn(2, 1);
            v_d      = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            
%             %sweep
%             [v_e0, v_values] = obj.sweep_e(v_d, m_s, obj.rho);
%             [worst_v, idx] = max(v_values);
%             e_worst_case = v_e0(idx);
%             

            e_0 = obj.rho;
            [v_x_0, ~]  = obj.estimatePosition(v_d, m_s, obj.rho, e_0);
            %v_x_0 =  [-0.2641; 0.5308];
            ale = ADMMLocationEstimator;
            v_x_hat = ale.solve_admm_loc(v_d, m_s, obj.rho, e_0, v_x_0);
            
            %% plot
%             figure(998); clf
%                 plot(v_e0, v_values);hold on
%                 plot(e_worst_case, worst_v, '*r');

            figure(999); clf
                plot(v_x_true(1), v_x_true(2), '+k'); hold on
                plot(m_s(1,:), m_s(2,:), 'xg')
                plot(v_x_0(1), v_x_0(2), '^r');
                plot(v_x_hat(1,:), v_x_hat(2,:))

%                 viscircles(m_s', v_r)
                
            F = GFigure.captureCurrentFigure;
        end
        
        %successive convex approximation
        function F = experiment_301(obj, niter)
            %rng(2)
            m_s      = obj.generateSources;
            v_x_true = 100*randn(2, 1);
            v_d      = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            
%             %sweep
%             [v_e0, v_values] = obj.sweep_e(v_d, m_s, obj.rho);
%             [worst_v, idx] = max(v_values);
%             e_worst_case = v_e0(idx);
%             

            e_0 = obj.rho;
            %[v_x_0, v_r_0]  = obj.estimatePosition(v_d, m_s, obj.rho, e_0);
            v_x_0 = zeros(2,1);
            v_r_0 = 100*ones(obj.N, 1);
            mle = SCALocationEstimator;
            [m_x_hat, m_r] = mle.solve_sca_loc(v_d, m_s, obj.rho, e_0, v_x_0, v_r_0);
            
            %% plot
%             figure(998); clf
%                 plot(v_e0, v_values);hold on
%                 plot(e_worst_case, worst_v, '*r');

            figure(999); clf
                plot(v_x_true(1), v_x_true(2), '+k'); hold on
                plot(m_s(1,:), m_s(2,:), 'xg')
                plot(v_x_0(1), v_x_0(2), '^r');
                plot(m_x_hat(1,:), m_x_hat(2,:), '<b')

                % viscircles(m_s', v_r_0)
                % viscircles(m_s', m_r(:,end), 'Color', 'blue')
                
            F = GFigure.captureCurrentFigure;
        end

    end
end

