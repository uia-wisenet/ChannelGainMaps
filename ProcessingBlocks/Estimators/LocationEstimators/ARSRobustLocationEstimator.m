classdef ARSRobustLocationEstimator < RobustLocationEstimator
    %Averaging Reference Sources to obtain a Robust Location Estimate
properties
    %param_rho       % upper bound on the NLOS bias - already defined in
    %superclass
    
    ch_maxType = 'relaxed' % options: 'greedy', 'relaxed'
end

methods
    
    function [v_x_hat, v_r_out, value, m_valueMap, uncertainty_measure...
            ] = solveGridSearch(obj, v_d_in, m_s, m_grid_x1, m_grid_x2)
        % Solve the convex relaxation of the problem of minimizing the
        % cost function implemented by the method computeValue
        m_valueMap = obj.valueMap(v_d_in, m_s, m_grid_x1, m_grid_x2);
        [value, idx] = min(m_valueMap(:));           
        v_x_hat = [m_grid_x1(idx); m_grid_x2(idx)];
        v_r_out = norms(v_x_hat - m_s);
        uncertainty_measure = nan;
        if nargout ==5
            warning('Uncertainty measure not implemented yet.')
        end
    end
    
    function varargout = computeValue(obj, varargin)
        varargout = cell(nargout, 1);
        [varargout{:}] = feval(...
            ['computeValue_' obj.ch_maxType 'Max'], obj, varargin{:});
    end
    
    function [value_out, v_e] = computeValue_greedyMax(obj, ...
            v_x, v_d_in, m_s, v_e_0)
        % Given a point in space, range difference measurements, and
        % source positions, computes an approximation (lower bound) to 
        % the maximization of the cost via a greedy algorithm
        N = size(m_s, 2);
        if not(exist('v_e_0', 'var'))
            % initialize at the optimum of the relaxed maximization:
            [~, ~, v_e_0] = obj.computeValue_relaxedMax(...
                v_x, v_d_in, m_s);
        end
        v_d1 = [0; v_d_in(:)];
        m_D = v_d1' - v_d1;
        v_r = norms(v_x - m_s);
                    
        %v_e_prev = [];
        max_greedy_iter = 1000;
        v_valueSeq = -inf(max_greedy_iter, 1);
        v_e = v_e_0;
        keyboard
        for k_iter = 2:max_greedy_iter
            m_neighbors = diag(obj.param_rho-2*v_e(:)) + v_e(:);
            v_valNeigh = nan(N, 1);
            for k = 1:N
                v_e_now = m_neighbors(:,k);
                m_toSum = (v_r(:)+v_e_now)' - (v_r(:)+v_e_now) - m_D;
                %set m_toSum's diag to 0:
                m_toSum = m_toSum -diag(diag(m_toSum));
                v_valNeigh(k) = sum(m_toSum(:).^2);
            end
            [v_value_now, idx_max] = max(v_valNeigh);
            if v_value_now > v_valueSeq(k_iter-1)
                %v_e_prev = v_e;
                v_e = m_neighbors(:,idx_max);
                v_valueSeq(k_iter) = v_value_now;
            else
                % v_e is a local maximum
                value_out = v_valueSeq(k_iter-1);
                return
            end            
        end
        warning('max num. of iter. exceeded')
    end
    
    function [value, v_r, argmax_e] = computeValue_relaxedMax(obj, ...
            v_x, v_d_in, m_s)
        % Given a point in space, range difference measurements, and
        % source positions, computes the objective value of the 
        % relaxed maximization (upper bound of the exact cost)
        N = size(m_s, 2);
        rho = obj.param_rho;
        zero_or_rho = [0 rho];
        assert(size(m_s, 1)==2);
        assert(length(v_d_in)==N-1);
        v_d1 = [0; v_d_in(:)];
        m_D = v_d1' - v_d1;
        v_r = norms(v_x - m_s);
        
        v_valueForRefSource = nan(N, 1);
        v_argmax_e = nan(N, 1);
        for i_s = 1:N % index pointing at the reference source
            val_0 = 0;
            val_rho = 0;
            for k = 1:N % index pointing at the source whose TDOA we add
                if k~=i_s
                    expr_0 = v_r(i_s) - v_r(k) - m_D(k, i_s);
                    expr_rho = expr_0 - rho;
                    val_0   = val_0   +expr_0^2   +2*rho*max(0, rho+expr_0);
                    val_rho = val_rho +expr_rho^2 +2*rho*max(0, rho+expr_rho);
                end
            end
            [v_valueForRefSource(i_s), idx_max] = max([val_0, val_rho]);
            argmax_e(i_s) = zero_or_rho(idx_max);
        end
        value = mean(v_valueForRefSource);
    end
    
    function m_map = valueMap(obj, v_d, m_s,m_grid_x1, m_grid_x2, m_grid_x3)
        % Given a grid of points as produced by meshgrid, returns a matrix 
        % containing the objective value for each point in the grid.
        
        % This function simply iterates over every point in the grid, calls
        % the method calculate_value and stores the output in a matrix.
        if obj.dim == 3
            error 'not implemented yet'
        elseif obj.dim ~= 2
            error 'wrong dimension, must be 2 or 3'
        end
        assert(isequal(size(m_grid_x1), size(m_grid_x2)), ...
            'size of grids not coherent')
        
        m_map = nan(size(m_grid_x1));
        
        ltc = LoopTimeControl(numel(m_map));
        for n = 1:numel(m_map)
            m_map(n) = obj.computeValue(...
                [m_grid_x1(n);m_grid_x2(n)], v_d, m_s);
            ltc.go(n);
        end
        
    end
    
    function [v_x_hat, v_r_out, value] = solveRelaxed(obj, ...
                v_d_in, m_s)
        % Solve the convex relaxation of the problem of minimizing the 
        % cost function implemented by the method computeValue_relaxedMax
        
        assert(strcmp(obj.ch_maxType, 'relaxed'))
            
        N = size(m_s, 2);
        rho = obj.param_rho;
        assert(size(m_s, 1)==2);
        assert(length(v_d_in)==N-1);
        v_d1 = [0; v_d_in(:)];
        m_D = v_d1' - v_d1;
        
        cvx_begin quiet
            variable x(2)
            variable v_r(N)
            for n = 1:N
                norm(x-m_s(:,n)) <= v_r(n);  %#ok<VUNUS>
                %relaxed constraint
            end
            variable v_valueForRefSource(N)
            disp 'Conforming constraints...'
            b_faster = 1;
            if b_faster
                expression m_value_0(N,N)
                expression m_value_rho(N,N)
                m_expr_0 = ones(N, 1)*v_r' - v_r*ones(1, N) - m_D;
                m_value_0   = m_expr_0.^2 + 2*rho*max(0, rho+m_expr_0);
                m_value_rho = (m_expr_0-rho).^2 + 2*rho*max(0, m_expr_0);
                for i = 1:N
                    m_value_0(i,i) = 0;
                    m_value_rho(i,i) = 0;
                end
                value_0 = sum(m_value_0);
                value_rho = sum(m_value_rho);
            else
                expression value_0(N); %#ok<UNRCH>
                expression value_rho(N);                
                for i_s = 1:N % index pointing at the reference source
                    for k = 1:N % index pointing at the source whose TDOA we add
                        if k~=i_s
                            expr_0 = v_r(i_s) - v_r(k) - m_D(k, i_s);
                            expr_rho = expr_0 - rho;
                            value_0(i_s)   = value_0(i_s)   + expr_0^2   + 2*rho*max(0, rho+expr_0);
                            value_rho(i_s)= value_rho(i_s) + expr_rho^2 + 2*rho*max(0, rho+expr_rho);
                        end
                    end
                end
            end
            v_valueForRefSource >= value_0';   %#ok<VUNUS>
            v_valueForRefSource >= value_rho'; %#ok<VUNUS>

            disp 'Solving relaxed problem...'
            v_precValues = cvx_precision;
            lambda = v_precValues(3)./trace(cov(m_s'));
            minimize(sum(v_valueForRefSource)/N + lambda*sum_square(v_r));
        cvx_end
        v_x_hat = x;
        v_r_out = v_r;
        value = cvx_optval;
    end
    
end

end