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
            Ns = size(m_s, 2);
            v_d = zeros(Ns-1, 1);
            v_r0 = norm(v_x - m_s(:, 1)) + rho*rand;
            for n = 2:Ns
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
                c_legend{3} = 'objective function';
                title 'worst case, some NLOS range meas.'
                legend(c_legend);
            subplot(1, 2, 2)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, m_value_LOS, contour_levels)
                c_legend{3} = 'objective function';
                title 'worst case,NLOS-free range meas.'
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

        % interpret value function as a negative log likelihood
        % visualize the associated PDF(s)
        function F = experiment_140(obj, n_coord_grid)
            rng(4);
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
            range1 = obj.exagg_linspace(min(m_points(1,:)), ...
                max(m_points(1,:)),  0, 1, n_coord_grid);
            range2 = obj.exagg_linspace(min(m_points(2,:)), ...
                max(m_points(2,:)), 0, 1, n_coord_grid);
            [m_grid_x1, m_grid_x2] = meshgrid(range1, range2);
            disp 'calculating value map...'
            m_value =     mrle.value_map( m_grid_x1, m_grid_x2, ...
                v_d,     m_s, obj.rho); 
            m_value_LOS = mrle.value_map( m_grid_x1, m_grid_x2, ...
                v_d_LOS, m_s, obj.rho); 
            [value_gridsearch, idx] = min(m_value(:));
            value_gs_LOS =  min(m_value_LOS(:));
%             v_x_grid_s = [m_grid_x1(idx); m_grid_x2(idx)];
            factor = 5000;
            unnormalized_pdf = exp((value_gridsearch-m_value)/factor);
            unnormalized_pdf_LOS = exp((value_gs_LOS-m_value_LOS)/factor);
            normalized_pdf     = ...
                unnormalized_pdf./sum(unnormalized_pdf(:));
            normalized_pdf_LOS = ...
                unnormalized_pdf_LOS./sum(unnormalized_pdf_LOS(:));

            n_cl = 20;
            v_pdfLevels = linspace(0.001, 1, n_cl)./sum(unnormalized_pdf(:));
            v_confidenceLevels =     zeros(n_cl, 1);
            v_confidenceLevels_LOS = zeros(n_cl, 1);
            for k_cl = 1:n_cl
                m_b = normalized_pdf >= v_pdfLevels(k_cl);
                m_b_LOS = normalized_pdf_LOS >= v_pdfLevels(k_cl);
                v_confidenceLevels(k_cl)     = sum(normalized_pdf(m_b(:)));
                v_confidenceLevels_LOS(k_cl) = sum(normalized_pdf_LOS(m_b_LOS(:)));
            end
            
            figure(996); clf
            plot(v_confidenceLevels, v_pdfLevels);
            hold on
            plot(v_confidenceLevels_LOS, v_pdfLevels);
            
            v_standardConfidenceLevels = linspace(0.6, 0.8, 5);
            contour_levels = interp1(v_confidenceLevels, v_pdfLevels, ...
                v_standardConfidenceLevels);
            contour_levels_LOS = interp1(v_confidenceLevels_LOS, v_pdfLevels, ...
                v_standardConfidenceLevels);
                
            figure(997); clf         
            subplot(1, 2, 1)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, normalized_pdf, contour_levels)
                c_legend{3} = 'objective function';
                title 'worst case, some NLOS range meas.'
                legend(c_legend);
            subplot(1, 2, 2)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, normalized_pdf_LOS, contour_levels_LOS)
                c_legend{3} = 'objective function';
                title 'worst case,NLOS-free range meas.'
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
               
            figure(999); clf
            ax (1) = subplot(1, 2, 1);
                mesh(m_grid_x1, m_grid_x2, unnormalized_pdf);
            ax (2) = subplot(1, 2, 2);
                mesh(m_grid_x1, m_grid_x2, unnormalized_pdf_LOS);
            linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'XLim', 'YLim', 'ZLim'})
            F(3) = GFigure.captureCurrentFigure;  
            
            figure(990); clf
                mesh(m_grid_x1, m_grid_x2, normalized_pdf); hold on
                mesh(m_grid_x1, m_grid_x2, normalized_pdf_LOS);            
            F(4) = GFigure.captureCurrentFigure;  
                                      
        end
        
        % interpret value function as a negative log likelihood
        % visualize the associated PDF(s)
        % Using ARSRobustLocationEstimator
        function F = experiment_145(obj, n_coord_grid)
            rng(4);
            obj.rho = 60;
            obj.sigma_range_d_noise = 25;
            [m_s, v_x_true] = obj.generateSources_Sensor();
            v_d             = obj.rangeDifferences(v_x_true, m_s, obj.rho);
            v_d_LOS         = obj.rangeDifferences(v_x_true, m_s, 0 );
            
             mrle = MyRobustLocationEstimator(obj.N);
             arle = ARSRobustLocationEstimator();
             arle.param_rho = obj.rho;
%             [v_x_hat, v_r_mrle, v_value_mrle] = ...
%                 mrle.solve_robust_loc(v_d, m_s, obj.rho);
%             value_mrle = mrle.calculate_value(v_x_hat, v_d, m_s, obj.rho);
%             
%             wle = WangLocationEstimator;
%             [v_x_wang, v_r_wang] = ...
%                 wle.solve_robust_loc(v_d, m_s, obj.rho);
%             value_wang = mrle.calculate_value(v_x_wang, v_d, m_s, obj.rho);
            tic
            [v_x_hat, v_r_mrle, v_value_mrle] = ...
                arle.solveRelaxed(v_d, m_s);     
            v_x_hat_LOS = ...
                arle.solveRelaxed(v_d_LOS, m_s); 
            time_relaxed = toc
            % create value map
            if isempty(n_coord_grid) % if the numer of grid points was not
                % passed as an input argument, take it from the property
                n_coord_grid = obj.n_corod_grid_default;
            end
            m_points = [m_s v_x_true];
            tic
            v_range1 = obj.exagg_linspace(min(m_points(1,:)), ...
                max(m_points(1,:)),  0, 1, n_coord_grid);
            v_range2 = obj.exagg_linspace(min(m_points(2,:)), ...
                max(m_points(2,:)), 0, 1, n_coord_grid);
            [m_grid_x1, m_grid_x2] = meshgrid(v_range1, v_range2);
            disp 'calculating value map...'
%             m_valueMap =     arle.valueMap( m_grid_x1, m_grid_x2, ...
%                 v_d,     m_s); 
%             m_valueMap_LOS = arle.valueMap( m_grid_x1, m_grid_x2, ...
%                 v_d_LOS, m_s); 
%             [value_gs, idx]    = min(m_valueMap(:));           
%             v_x_hat_gs = [m_grid_x1(idx); m_grid_x2(idx)];
%             [value_gs_LOS, idx] =  min(m_valueMap_LOS(:));
%             v_x_hat_gs_LOS = [m_grid_x1(idx); m_grid_x2(idx)];
            [v_x_hat_gs, ~, value_gs, m_valueMap] = ...
                arle.solveGridSearch(v_d, m_s, m_grid_x1, m_grid_x2);
            [v_x_hat_gs_LOS, ~, value_gs_LOS, m_valueMap_LOS] = ...
                arle.solveGridSearch(v_d_LOS, m_s, m_grid_x1, m_grid_x2);
            time_gridSearch = toc
            %
            factor = 5000;
            unnormalized_pdf = exp((value_gs-m_valueMap)/factor);
            unnormalized_pdf_LOS = exp((value_gs_LOS-m_valueMap_LOS)/factor);
            normalized_pdf     = ...
                unnormalized_pdf./sum(unnormalized_pdf(:));
            normalized_pdf_LOS = ...
                unnormalized_pdf_LOS./sum(unnormalized_pdf_LOS(:));

            n_cl = 20;
            v_pdfLevels = linspace(0.001, 1, n_cl)./sum(unnormalized_pdf(:));
            v_confidenceLevels =     zeros(n_cl, 1);
            v_confidenceLevels_LOS = zeros(n_cl, 1);
            for k_cl = 1:n_cl
                m_b = normalized_pdf >= v_pdfLevels(k_cl);
                m_b_LOS = normalized_pdf_LOS >= v_pdfLevels(k_cl);
                v_confidenceLevels(k_cl)     = sum(normalized_pdf(m_b(:)));
                v_confidenceLevels_LOS(k_cl) = sum(normalized_pdf_LOS(m_b_LOS(:)));
            end
            
            figure(996); clf
            plot(v_confidenceLevels, v_pdfLevels);
            hold on
            plot(v_confidenceLevels_LOS, v_pdfLevels);
            
            v_standardConfidenceLevels = linspace(0.6, 0.8, 5);
            contour_levels_pdf = interp1(v_confidenceLevels, v_pdfLevels, ...
                v_standardConfidenceLevels);
            contour_levels_pdf_LOS = interp1(v_confidenceLevels_LOS, v_pdfLevels, ...
                v_standardConfidenceLevels);
                
            figure(995); clf         
            subplot(1, 2, 1)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(v_range1, v_range2, normalized_pdf , contour_levels_pdf)
                c_legend{3} ='interpreted PDF';
                obj.plot_estimates(v_x_hat, obj.ch_altEstimate)
                c_legend{4} = 'Relaxed ARS';
                obj.plot_estimates(v_x_hat_gs, obj.ch_gsEstimate)
                c_legend{5} = 'ARS, grid search';
                title 'worst case, some NLOS range meas.'
                legend(c_legend);
            subplot(1, 2, 2)
                c_legend{1} = obj.plot_sources(m_s); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(v_range1, v_range2, normalized_pdf_LOS, ...
                    contour_levels_pdf_LOS)
                c_legend{3} = 'interpreted PDF';
                obj.plot_estimates(v_x_hat_LOS, obj.ch_altEstimate)
                c_legend{4} = 'Relaxed ARS';
                obj.plot_estimates(v_x_hat_gs_LOS, obj.ch_gsEstimate)
                c_legend{5} = 'ARS, grid search';
                title 'worst case, NLOS-free range meas.'
                legend(c_legend);
            F = GFigure.captureCurrentFigure;
            
            figure(998); clf
            ax (1) = subplot(1, 2, 1);
                mesh(m_grid_x1, m_grid_x2, m_valueMap);
            ax (2) = subplot(1, 2, 2);
                mesh(m_grid_x1, m_grid_x2, m_valueMap_LOS);
                % stem3(v_x_hat(1),    v_x_hat(2),    value_mrle, '^m');
                % stem3(v_x_wang(1),   v_x_wang(2),   value_wang, 'vr')
                % stem3(v_x_grid_s(1), v_x_grid_s(2), value_gridsearch, 'sm')
            linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'XLim', 'YLim', 'ZLim'})
            F(2) = GFigure.captureCurrentFigure;
               
            figure(999); clf
            ax (1) = subplot(1, 2, 1);
                mesh(m_grid_x1, m_grid_x2, unnormalized_pdf);
            ax (2) = subplot(1, 2, 2);
                mesh(m_grid_x1, m_grid_x2, unnormalized_pdf_LOS);
            linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'XLim', 'YLim', 'ZLim'})
            F(3) = GFigure.captureCurrentFigure;  
            
            figure(990); clf
                mesh(m_grid_x1, m_grid_x2, normalized_pdf); hold on
                mesh(m_grid_x1, m_grid_x2, normalized_pdf_LOS);            
            F(4) = GFigure.captureCurrentFigure;  
                                      
        end
        
        function F = experiment_150(obj, n_coord_grid)
            rng(4);
            [m_s{1}, v_x_true] = obj.generateSources_Sensor();
            [m_s{2}] = obj.generateSources_Sensor();
            v_d{1}   = obj.rangeDifferences(v_x_true, m_s{1}, obj.rho);
            v_d{2}   = obj.rangeDifferences(v_x_true, m_s{2}, obj.rho);
            %v_d_LOS = obj.rangeDifferences(v_x_true, m_s{1}, 0 );
            
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
            m_points = [m_s{1} v_x_true];
            range1 = obj.exagg_linspace(min(m_points(1,:)), ...
                max(m_points(1,:)),  0, 1, n_coord_grid);
            range2 = obj.exagg_linspace(min(m_points(2,:)), ...
                max(m_points(2,:)), 0, 1, n_coord_grid);
            [m_grid_x1, m_grid_x2] = meshgrid(range1, range2);
            disp 'calculating value map...'
            factor = 5000; % for the time being, a magic number :(
            for k = 1:2
                cm_value{k} = mrle.value_map( m_grid_x1, m_grid_x2, ...
                    v_d{k}, m_s{k}, obj.rho);
                v_value_gs(k) = min(cm_value{k}(:));
                cm_unPdf{k} = exp((v_value_gs(k)-cm_value{k})/factor);
                cm_nPdf{k}  = cm_unPdf{k}./sum(cm_unPdf{k}(:)); %normalized pdf
                c_sumUP(k)  = sum(cm_unPdf{k}(:));
            end

            n_cl = 20;
            v_pdfLevels = linspace(0, 1/min(c_sumUP), n_cl);
            for k = 1:2
                cv_confidenceLevels{k} = zeros(n_cl, 1);
                for k_cl = 1:n_cl
                    m_b = (cm_nPdf{k} >= v_pdfLevels(k_cl));
                    cv_confidenceLevels{k}(k_cl) = sum(cm_nPdf{k}(m_b(:)));
                end
            end
            
            figure(996); clf
            plot(cv_confidenceLevels{1}, v_pdfLevels);
            hold on
            plot(cv_confidenceLevels{2}, v_pdfLevels);
            
            v_standardConfidenceLevels = linspace(0.6, 0.8, 5);
            for k = 1:2
                [v_cLevels, v_idx] = unique(cv_confidenceLevels{k});                
                contour_levels{k} = interp1(v_cLevels, v_pdfLevels(v_idx), ...
                    v_standardConfidenceLevels);                
            end
            figure(997); clf 
            for k = 1:2
               subplot(1, 2, k)
                c_legend{1} = obj.plot_sources(m_s{k}); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, cm_nPdf{k}, contour_levels{k})
                c_legend{3} = 'objective function';
                title(sprintf('worst case, source set %d', k))
                legend(c_legend);
            end
            F = GFigure.captureCurrentFigure;
            
%             figure(998); clf
%             ax (1) = subplot(1, 2, 1);
%                 mesh(m_grid_x1, m_grid_x2, m_value);
%             ax (2) = subplot(1, 2, 2);
%                 mesh(m_grid_x1, m_grid_x2, m_value_LOS);
%                 % stem3(v_x_hat(1),    v_x_hat(2),    value_mrle, '^m');
%                 % stem3(v_x_wang(1),   v_x_wang(2),   value_wang, 'vr')
%                 % stem3(v_x_grid_s(1), v_x_grid_s(2), value_gridsearch, 'sm')
%             linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
%                 'CameraTarget', 'XLim', 'YLim', 'ZLim'})
%             F(2) = GFigure.captureCurrentFigure;
%                
%             figure(999); clf
%             ax (1) = subplot(1, 2, 1);
%                 mesh(m_grid_x1, m_grid_x2, unnormalized_pdf);
%             ax (2) = subplot(1, 2, 2);
%                 mesh(m_grid_x1, m_grid_x2, unnormalized_pdf_LOS);
%             linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
%                 'CameraTarget', 'XLim', 'YLim', 'ZLim'})
%             F(3) = GFigure.captureCurrentFigure;  
%             
%             figure(990); clf
%                 mesh(m_grid_x1, m_grid_x2, normalized_pdf); hold on
%                 mesh(m_grid_x1, m_grid_x2, normalized_pdf_LOS);            
%             F(4) = GFigure.captureCurrentFigure;  
                                      
        end

        function F = experiment_160(obj, n_coord_grid)
            rng(8);
            obj.N = 10;
            [m_s{1}, v_x_true] = obj.generateSources_Sensor();
            m_s{2} = obj.generateSources_Sensor();
            m_s{2}(:,end-3:end) = [];
            v_d{1}   = obj.rangeDifferences(v_x_true, m_s{1}, obj.rho);
            v_d{2}   = obj.rangeDifferences(v_x_true, m_s{2}, obj.rho);
            %v_d_LOS = obj.rangeDifferences(v_x_true, m_s{1}, 0 );
            
            for k = 1:2
                mrle{k} = MyRobustLocationEstimator(size(m_s{k},2));
            end
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
            m_points = [cat(2, m_s{:}) v_x_true];
            range1 = obj.exagg_linspace(min(m_points(1,:)), ...
                max(m_points(1,:)),  0, 1, n_coord_grid);
            range2 = obj.exagg_linspace(min(m_points(2,:)), ...
                max(m_points(2,:)), 0, 1, n_coord_grid);
            [m_grid_x1, m_grid_x2] = meshgrid(range1, range2);
            disp 'calculating value map...'
            factor = 5000; % for the time being, a magic number :(
            for k = 1:2
                cm_value{k} = mrle{k}.value_map( m_grid_x1, m_grid_x2, ...
                    v_d{k}, m_s{k}, obj.rho);
                v_value_gs(k) = min(cm_value{k}(:));
                cm_unPdf{k} = exp((v_value_gs(k)-cm_value{k})/factor);
                cm_nPdf{k}  = cm_unPdf{k}./sum(cm_unPdf{k}(:)); %normalized pdf
                c_sumUP(k)  = sum(cm_unPdf{k}(:));
            end

            n_cl = 20;
            v_pdfLevels = linspace(0, 1/min(c_sumUP), n_cl);
            for k = 1:2
                cv_confidenceLevels{k} = zeros(n_cl, 1);
                for k_cl = 1:n_cl
                    m_b = (cm_nPdf{k} >= v_pdfLevels(k_cl));
                    cv_confidenceLevels{k}(k_cl) = sum(cm_nPdf{k}(m_b(:)));
                end
            end
            
            figure(996); clf
            plot(cv_confidenceLevels{1}, v_pdfLevels);
            hold on
            plot(cv_confidenceLevels{2}, v_pdfLevels);
            
            v_standardConfidenceLevels = linspace(0.8, 0.95, 5);
            for k = 1:2
                [v_cLevels, v_idx] = unique(cv_confidenceLevels{k});                
                contour_levels{k} = interp1(v_cLevels, v_pdfLevels(v_idx), ...
                    v_standardConfidenceLevels);                
            end
            figure(997); clf 
            for k = 1:2
               subplot(1, 2, k)
                c_legend{1} = obj.plot_sources(m_s{k}); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, cm_nPdf{k}, contour_levels{1})%!
                c_legend{3} = 'objective function';
                title(sprintf('worst case, source set %d', k))
                legend(c_legend);
            end
            F = GFigure.captureCurrentFigure;
            
            figure(998); clf
            
%             contour_levels_LOS = linspace(...
%                 value_gs_LOS, 10*value_gs_LOS, 30);
            for k = 1:2
               subplot(1, 2, k)
                contour_levels_value = linspace(...
                    v_value_gs(k), v_value_gs(k)+ v_value_gs(1), 5);
                c_legend{1} = obj.plot_sources(m_s{k}); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, cm_value{k}, contour_levels_value)
                c_legend{3} = 'objective function';
                title(sprintf('worst case, source set %d', k))
                legend(c_legend);
            end
            F(2) = GFigure.captureCurrentFigure;
            
%             figure(998); clf
%             ax (1) = subplot(1, 2, 1);
%                 mesh(m_grid_x1, m_grid_x2, m_value);
%             ax (2) = subplot(1, 2, 2);
%                 mesh(m_grid_x1, m_grid_x2, m_value_LOS);
%                 % stem3(v_x_hat(1),    v_x_hat(2),    value_mrle, '^m');
%                 % stem3(v_x_wang(1),   v_x_wang(2),   value_wang, 'vr')
%                 % stem3(v_x_grid_s(1), v_x_grid_s(2), value_gridsearch, 'sm')
%             linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
%                 'CameraTarget', 'XLim', 'YLim', 'ZLim'})
%             F(2) = GFigure.captureCurrentFigure;
%                
%             figure(999); clf
%             ax (1) = subplot(1, 2, 1);
%                 mesh(m_grid_x1, m_grid_x2, unnormalized_pdf);
%             ax (2) = subplot(1, 2, 2);
%                 mesh(m_grid_x1, m_grid_x2, unnormalized_pdf_LOS);
%             linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
%                 'CameraTarget', 'XLim', 'YLim', 'ZLim'})
%             F(3) = GFigure.captureCurrentFigure;  
%             
%             figure(990); clf
%                 mesh(m_grid_x1, m_grid_x2, normalized_pdf); hold on
%                 mesh(m_grid_x1, m_grid_x2, normalized_pdf_LOS);            
%             F(4) = GFigure.captureCurrentFigure;  
                                      
        end

        % First experiment with ARSRobustLocationEstimator
        function F = experiment_200(obj, n_coord_grid)
            rng(5);
            obj.N = 15;
            N_fewer = 6;
            [m_s{1}, v_x_true] = obj.generateSources_Sensor();
            m_s{2} = obj.generateSources_Sensor();
            m_s{2} = m_s{2}(:,1:N_fewer);
            v_d{1}   = obj.rangeDifferences(v_x_true, m_s{1}, obj.rho);
            v_d{2}   = obj.rangeDifferences(v_x_true, m_s{2}, obj.rho);
            %v_d_LOS = obj.rangeDifferences(v_x_true, m_s{1}, 0 );
            
            for k = 1:2
                arle{k} = ARSRobustLocationEstimator();
                arle{k}.param_rho = obj.rho;
            end
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
            m_points = [cat(2, m_s{:}) v_x_true];
            range1 = obj.exagg_linspace(min(m_points(1,:)), ...
                max(m_points(1,:)),  0, 1, n_coord_grid);
            range2 = obj.exagg_linspace(min(m_points(2,:)), ...
                max(m_points(2,:)), 0, 1, n_coord_grid);
            [m_grid_x1, m_grid_x2] = meshgrid(range1, range2);
            disp 'calculating value map...'
            factor = 5000; % for the time being, a magic number :(
            for k = 1:2
                cm_value{k} = arle{k}.value_map( m_grid_x1, m_grid_x2, ...
                    v_d{k}, m_s{k});
                v_value_gs(k) = min(cm_value{k}(:));
                cm_unPdf{k} = exp((v_value_gs(k)-cm_value{k})/factor);
                cm_nPdf{k}  = cm_unPdf{k}./sum(cm_unPdf{k}(:)); %normalized pdf
                c_sumUP(k)  = sum(cm_unPdf{k}(:));
            end

            n_cl = 20;
            v_pdfLevels = linspace(0, 1/min(c_sumUP), n_cl);
            for k = 1:2
                cv_confidenceLevels{k} = zeros(n_cl, 1);
                for k_cl = 1:n_cl
                    m_b = (cm_nPdf{k} >= v_pdfLevels(k_cl));
                    cv_confidenceLevels{k}(k_cl) = sum(cm_nPdf{k}(m_b(:)));
                end
            end
            
            figure(996); clf
            plot(cv_confidenceLevels{1}, v_pdfLevels);
            hold on
            plot(cv_confidenceLevels{2}, v_pdfLevels);
            
            v_standardConfidenceLevels = linspace(0.8, 0.95, 5);
            for k = 1:2
                [v_cLevels, v_idx] = unique(cv_confidenceLevels{k});                
                contour_levels{k} = interp1(v_cLevels, v_pdfLevels(v_idx), ...
                    v_standardConfidenceLevels);                
            end
            figure(997); clf 
            for k = 1:2
               subplot(1, 2, k)
                c_legend{1} = obj.plot_sources(m_s{k}); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, cm_nPdf{k}, contour_levels{1})%!
                c_legend{3} = 'objective function';
                title(sprintf('worst case, source set %d', k))
                legend(c_legend);
            end
            F = GFigure.captureCurrentFigure;
            
            figure(998); clf
            
%             contour_levels_LOS = linspace(...
%                 value_gs_LOS, 10*value_gs_LOS, 30);
            for k = 1:2
               subplot(1, 2, k)
                contour_levels_value = linspace(...
                    v_value_gs(k), v_value_gs(k)+ v_value_gs(1), 5);
                c_legend{1} = obj.plot_sources(m_s{k}); hold on
                c_legend{2} = obj.plot_sensor(v_x_true);
                contour(range1, range2, cm_value{k}, contour_levels_value)
                c_legend{3} = 'objective function';
                title(sprintf('worst case, source set %d', k))
                legend(c_legend);
            end
            F(2) = GFigure.captureCurrentFigure;
            
%             figure(998); clf
%             ax (1) = subplot(1, 2, 1);
%                 mesh(m_grid_x1, m_grid_x2, m_value);
%             ax (2) = subplot(1, 2, 2);
%                 mesh(m_grid_x1, m_grid_x2, m_value_LOS);
%                 % stem3(v_x_hat(1),    v_x_hat(2),    value_mrle, '^m');
%                 % stem3(v_x_wang(1),   v_x_wang(2),   value_wang, 'vr')
%                 % stem3(v_x_grid_s(1), v_x_grid_s(2), value_gridsearch, 'sm')
%             linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
%                 'CameraTarget', 'XLim', 'YLim', 'ZLim'})
%             F(2) = GFigure.captureCurrentFigure;
%                
%             figure(999); clf
%             ax (1) = subplot(1, 2, 1);
%                 mesh(m_grid_x1, m_grid_x2, unnormalized_pdf);
%             ax (2) = subplot(1, 2, 2);
%                 mesh(m_grid_x1, m_grid_x2, unnormalized_pdf_LOS);
%             linkprop(ax, {'CameraUpVector', 'CameraPosition', ...
%                 'CameraTarget', 'XLim', 'YLim', 'ZLim'})
%             F(3) = GFigure.captureCurrentFigure;  
%             
%             figure(990); clf
%                 mesh(m_grid_x1, m_grid_x2, normalized_pdf); hold on
%                 mesh(m_grid_x1, m_grid_x2, normalized_pdf_LOS);            
%             F(4) = GFigure.captureCurrentFigure;  
                                      
        end

    end
    
    methods (Static)
        function v_out = exagg_linspace(x_0,  x_1, ...
                min_z, max_z, npoints)
            v_z = linspace(min_z, max_z, npoints);
            v_out = x_0 + (x_1-x_0)*v_z;
        end
    end
end

