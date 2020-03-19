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
                        
        n_corod_grid_default  = 500
        
        ch_sourceSymbol = 'xk'
        ch_sensorSymbol = 'ok'
        ch_wangEstimate = '^m'
        ch_altEstimate  = 'vr'
        ch_gsEstimate   = 'sm'
    end
    
    methods % auxiliary functions: data generation and range differences
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
    end
    
    methods % plotting routines
        function legend_entry = plot_sources(obj, m_s)
            plot(m_s(1,:), m_s(2,:), obj.ch_sourceSymbol);
            legend_entry = 'source location';
        end
        function legend_entry = plot_sensor(obj, v_x)
            plot(v_x(1), v_x(2), obj.ch_sensorSymbol);
            legend_entry = 'sensor location';
        end
        function plot_estimates(obj, m_x_hat, ch_symbol)
            plot(m_x_hat(1,:), m_x_hat(2,:), ch_symbol);
        end
        
    end
    
    methods % experiments
        
        function F = experiment_101(obj, niter)
            rng(3);
            [m_s, v_x_true] = obj.generateSources_Sensor();
            v_d             = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            
            %sweep
            mrle = MyRobustLocationEstimator(obj.N);
            mrle.lambda = 1e-5;
%             [v_e0, v_values, v_x_hat, v_r, v_wc] = ...
%                 mrle.sweep_e(v_d, m_s, obj.rho);
            [v_x_hat, v_r, value] = mrle.solve_robust_loc(...
                v_d, m_s, obj.rho);
            
            wle = WangLocationEstimator;
            [v_x_wang, v_r_wang] = wle.solve_robust_loc(...
                v_d, m_s, obj.rho);

            %% plot
            figure(999); clf
                c_legend{1} = obj.plot_sensor(v_x_true); hold on
                c_legend{2} = obj.plot_sources(m_s);
                obj.plot_estimates(v_x_hat, obj.ch_altEstimate); 
                c_legend{3} = 'convex relax. 1';
                obj.plot_estimates(v_x_wang, obj.ch_wangEstimate);
                c_legend{4} = 'Wang et al.';
                %viscircles(m_s', v_r, 'Color', 'blue')
            legend(c_legend);
            F = GFigure.captureCurrentFigure;
        end
        
        % visualize the estimates when different sensors are taken as the
        % reference source
        function F = experiment_110(obj, niter)
            rng(3);
            [m_s, v_x_true] = obj.generateSources_Sensor();
            v_d             = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            
            mrle = MyRobustLocationEstimator(obj.N);
            wle = WangLocationEstimator;

            ltc = LoopTimeControl(obj.N);
            for k = 1:obj.N
                % take source k as the reference source
                v_d_all = [0; v_d];
                v_d_k = v_d_all - v_d_all(k);
                v_d_k(k) = [];
                m_my_s = [m_s(:,k), m_s];
                m_my_s(:,k+1) = [];
                [c_v_x_hat{k}, v_r_mrle, v_value_mrle(k)] = ...
                    mrle.solve_robust_loc(v_d_k, m_my_s, obj.rho);
            
                [c_v_x_wang{k}, v_r_wang, v_value_wang(k)] = ...
                    wle.solve_robust_loc(v_d_k, m_my_s, obj.rho);
                ltc.go(k);
            end
            
            v_value_wc_normalized   = v_value_mrle  /max(v_value_mrle  );
            v_value_wang_no_inf = v_value_wang(not(isinf(v_value_wang)));
            v_value_wang_normalized = v_value_wang/max(v_value_wang_no_inf);

            
            %% plot           
            figure(998); clf
                plot3(v_x_true(1), v_x_true(2), 0,           '+k'); hold on
                plot3(m_s(1,:),    m_s(2,:), zeros(1,obj.N), 'xg');
                stem3(c_v_x_hat{1}(1),  c_v_x_hat{1}(2),  ...
                        v_value_wc_normalized  (1), 'xr');
                stem3(c_v_x_wang{1}(1), c_v_x_wang{1}(2), ...
                            v_value_wang_normalized(1), 'xb')    
                for k = 1:obj.N
                    stem3(c_v_x_hat{k}(1),  c_v_x_hat{k}(2),  ...
                        v_value_wc_normalized  (k), '^r');
                    if not(isinf(v_value_wang(k)))
                        stem3(c_v_x_wang{k}(1), c_v_x_wang{k}(2), ...
                            v_value_wang_normalized(k), 'vb')
                    end
                end
            F = GFigure.captureCurrentFigure;
            
            figure(999); clf
                c_legend{1} = obj.plot_sensor(v_x_true); hold on
                c_legend{2} = obj.plot_sources(m_s);
                
                obj.plot_estimates(c_v_x_hat{1},  ['x' obj.ch_altEstimate(2) ]); 
                c_legend{3} = 'convex relax. 1, def. ref. source';
                obj.plot_estimates(c_v_x_wang{1}, ['x' obj.ch_wangEstimate(2)]);
                c_legend{4} = 'Wang et al., def. ref. source';
                for k = 1:obj.N
                    plot(c_v_x_hat{k}(1), c_v_x_hat{k}(2), obj.ch_altEstimate);
                    plot(c_v_x_wang{k}(1), c_v_x_wang{k}(2), obj.ch_wangEstimate)
                end
                c_legend{5} = 'convex relax. 1';
                c_legend{6} = 'Wang et al.';
            legend(c_legend);
                                
                %viscircles(m_s', v_r)               
            F(2) = GFigure.captureCurrentFigure;
        end
        
        % trace value function as a map for visualization
        function F = experiment_120(obj, n_coord_grid)
            rng(3);
            [m_s, v_x_true] = obj.generateSources_Sensor();
            v_d             = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            
            mrle = MyRobustLocationEstimator(obj.N);
            mrle.N = obj.N;
            [v_x_hat, v_r_mrle, v_value_mrle] = ...
                mrle.solve_robust_loc(v_d, m_s, obj.rho);
            value_mrle = mrle.calculate_value(v_x_hat, v_d, m_s, obj.rho);
            
            wle = WangLocationEstimator;
            [v_x_wang, v_r_wang] = ...
                wle.solve_robust_loc(v_d, m_s, obj.rho);
            value_wang = mrle.calculate_value(v_x_wang, v_d, m_s, obj.rho);
                
            % trace value map
            if isempty(n_coord_grid)
                n_coord_grid = obj.n_corod_grid_default;
            end
            m_points = [m_s v_x_true];
            range1 = linspace(min(m_points(1,:)), max(m_points(1,:)));
            range2 = linspace(min(m_points(2,:)), max(m_points(2,:)));
            [m_grid_x1, m_grid_x2] = meshgrid(range1, range2);
            disp 'calculating value map...'
            m_value = mrle.value_map(m_grid_x1, m_grid_x2, v_d, m_s, obj.rho); 
            [value_gridsearch, idx] = min(m_value(:));
            v_x_grid_s = [m_grid_x1(idx); m_grid_x2(idx)];
            
            figure(998); clf
                mesh(m_grid_x1, m_grid_x2, m_value);
                hold on
                stem3(v_x_hat(1),    v_x_hat(2),    value_mrle, obj.ch_altEstimate);
                stem3(v_x_wang(1),   v_x_wang(2),   value_wang, obj.ch_wangEstimate)
                stem3(v_x_grid_s(1), v_x_grid_s(2), value_gridsearch, obj.ch_gsEstimate);
                legend ('objective function', 'convex relax. 1', ...
                    'Wang et al.', 'grid search')
            F = GFigure.captureCurrentFigure; 
            
            contour_levels = linspace(value_gridsearch, 10*value_gridsearch, 30);
            figure(999); clf
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);             
                contour(range1, range2, m_value, contour_levels);
                c_legend{3} = 'objective function';
                obj.plot_estimates(v_x_hat, obj.ch_altEstimate)
                c_legend{4} = 'convex relax. 1';
                obj.plot_estimates(v_x_wang, obj.ch_wangEstimate)
                c_legend{5} = 'Wang et al.';
                obj.plot_estimates(v_x_grid_s, obj.ch_gsEstimate)
                c_legend{6} = 'grid search';  
            legend(c_legend)
            F(2) = GFigure.captureCurrentFigure;
        end
        
        % trace value function as a map for visualization
        function F = experiment_130(obj, n_coord_grid)
            rng(3);
            [m_s, v_x_true] = obj.generateSources_Sensor();
            v_d             = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            v_d_LOS         = obj.rangeDifferences(v_x_true, m_s, 0 );
            
             mrle = MyRobustLocationEstimator(obj.N);
%             [v_x_hat, v_r_mrle, v_value_mrle] = ...
%                 mrle.solve_robust_loc(v_d, m_s, obj.rho);
%             value_mrle = mrle.calculate_value(v_x_hat, v_d, m_s, obj.rho);
%             
%             wle = WangLocationEstimator;
%             [v_x_wang, v_r_wang] = ...
%                 wle.solve_robust_loc(v_d, m_s, obj.rho);
%             value_wang = mrle.calculate_value(v_x_wang, v_d, m_s, obj.rho);
                
            % create value map
            if isempty(n_coord_grid)
                n_coord_grid = obj.n_corod_grid_default;
            end
            m_points = [m_s v_x_true];
            range1 = linspace(min(m_points(1,:)), max(m_points(1,:)));
            range2 = linspace(min(m_points(2,:)), max(m_points(2,:)));
            [m_grid_x1, m_grid_x2] = meshgrid(range1, range2);
            disp 'calculating value map...'
            m_value =     mrle.value_map( m_grid_x1, m_grid_x2, ...
                v_d,     m_s, obj.rho); 
            m_value_LOS = mrle.value_map( m_grid_x1, m_grid_x2, ...
                v_d_LOS, m_s, obj.rho); 
            [value_gridsearch, idx] = min(m_value(:));
            value_gs_LOS =  min(m_value_LOS(:));
%             v_x_grid_s = [m_grid_x1(idx); m_grid_x2(idx)];
            
            

            figure(997); clf
            contour_levels = linspace(...
                value_gridsearch, 10*value_gridsearch, 30);
            contour_levels_LOS = linspace(...
                value_gs_LOS, 10*value_gs_LOS, 30);
            subplot(1, 2, 1)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, m_value, contour_levels)
                c_legend{3} = 'objective function'
                title 'worst case, some NLOS range meas.'
                legend(c_legend);
            subplot(1, 2, 2)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, m_value_LOS, contour_levels)
                c_legend{3} = 'objective function'
                title 'worst case, LOS-free range meas.'
                legend(c_legend);
            F = GFigure.captureCurrentFigure;
               
            figure(998); clf
            ax (1) = subplot(1, 2, 1);
                mesh(m_grid_x1, m_grid_x2, m_value);
            ax (2) = subplot(1, 2, 2);
                mesh(m_grid_x1, m_grid_x2, m_value_LOS);
                % stem3(v_x_hat(1),    v_x_hat(2),    value_mrle, '^m');
                % stem3(v_x_wang(1),   v_x_wang(2),   value_wang, 'vr')
                % stem3(v_x_grid_s(1), v_x_grid_s(2), value_gridsearch, 'sm')
            linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'XLim', 'YLim', 'ZLim'})
            F(2) = GFigure.captureCurrentFigure;          
                                      
        end


    end
end

