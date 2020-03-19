classdef MyRobustLocationEstimator
    properties
        N % number of sources
        ngrid_e = 2 % number of values for e_0 (NLOS bias on the reference 
        % source) that we try: for the simplified relaxed robust estimator,
        % it should be enough with just 0 and rho
        lambda = 0 % regularization parameter: in case the relaxed problem
        % has multiple solutions, if \lambda is small and > 0 we retrieve
        % the one with minimal r
    end
    methods
        function obj = MyRobustLocationEstimator(N) %constructor
            obj.N = N;
        end
        function[x_out, v_r_out, value_out] = solve_robust_loc(obj, ...
                v_d, m_s, rho)
            obj.N = numel(v_d)+1;
            [~, ~, x_out, v_r_wc, value_wc, e0_wc] = ...
                obj.sweep_e( v_d, m_s, rho);
            
            [value_old, v_r_old] = ...
                obj.calculate_value_old(x_out, v_d, m_s, rho, e0_wc);   
            [value_out, v_r_out] = obj.calculate_value(x_out, v_d, m_s, rho);
            % TODO!
            %norm(value_out-value_old) 
            %keyboard
        end
        
        function [v, v_r] = calculate_value(obj, v_x, v_d, m_s, rho)
            v_r = nan(obj.N,1);
            v1 = 0;
            v2 = 0;
            for k = 1:obj.N
                v_r(k) = norm(v_x-m_s(:,k));
            end
            for n = 2:obj.N
                my_expr = v_d(n-1) - v_r(n) + v_r(1) + 0;
                v1 = v1 + my_expr^2       + rho*max(0, rho-2* my_expr     );
                v2 = v2 + (my_expr+rho)^2 + rho*max(0, rho-2*(my_expr+rho));
            end
            v = max(v1, v2);
        end
        
        function [v, v_r] = calculate_value_old(obj, x_out, v_d, m_s, rho, e_0)
            v_r = nan(obj.N,1);
            v = 0;
            for k = 1:obj.N
                v_r(k) = norm(x_out-m_s(:,k));
            end
            for n = 2:obj.N
                my_expr = v_d(n-1) - v_r(n) + v_r(1) + e_0;
                v = v + my_expr^2 + rho*max(0, rho-2*my_expr);
            end
        end
        
        function m_map = value_map(obj, m_grid_x1, m_grid_x2, v_d, m_s, rho)
            assert(isequal(size(m_grid_x1), size(m_grid_x2)), ...
                'size of grids not coherent')
            m_map = nan(size(m_grid_x1));
            %val_aux = zeros(2, 1);
            ltc = LoopTimeControl(numel(m_map));
            for n = 1:numel(m_map)
%                 val_aux(1) = obj.calculate_value_old(...
%                     [m_grid_x1(n);m_grid_x2(n)], v_d, m_s, rho, 0);
%                 val_aux(2) = obj.calculate_value_old(...
%                     [m_grid_x1(n);m_grid_x2(n)], v_d, m_s, rho, rho);
                m_map(n) = obj.calculate_value(...
                    [m_grid_x1(n);m_grid_x2(n)], v_d, m_s, rho);
                %norm(m_map(n)-max(val_aux))
                ltc.go(n);
            end
            
        end
        
        function [v_x_hat, v_r_out, value] = estimatePosition(obj, ...
                v_d, m_s, rho, e_0)
            
            assert(obj.N-1 == length(v_d));
            cvx_begin quiet
                variable x(2)
                variable v_r(obj.N)
                objective = obj.lambda*sum_square(v_r);
                for n = 2:obj.N
                    my_expr = v_d(n-1) - v_r(n) + v_r(1) + e_0;
                    objective = objective + ...
                        my_expr^2 + rho*max(0, rho-2*my_expr);                   
                end
                minimize(objective)
                subject to
                for n = 1:obj.N
                    norm(x-m_s(:,n)) <= v_r(n); %relaxed constraint
                end
            cvx_end
            v_x_hat = x;
            v_r_out = v_r;
            value = cvx_optval;
        end
        
        function [v_e0, v_values, x_wc, r_wc, value_wc, e0_wc] = ...
                sweep_e(obj, v_d, m_s, rho)
            v_e0 = linspace(0, rho, obj.ngrid_e);
            v_values = zeros(size(v_e0));
            c_x= cell(obj.ngrid_e,1);
            c_r = c_x;
            for k = 1:obj.ngrid_e
                [c_x{k}, c_r{k}, v_values(k)] = obj.estimatePosition(...
                    v_d, m_s, rho, v_e0(k));
            end
            [value_wc, idx] = max(v_values);
            e0_wc = v_e0(idx);
            x_wc = c_x{idx};
            r_wc = c_r{idx};
        end
    end
end
