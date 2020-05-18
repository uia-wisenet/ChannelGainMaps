classdef ARSRobustLocationEstimator < LocationEstimator_tdoa
    %Averaging Reference Sources to obtain a Robust Location Estimate
properties
    param_rho(1,1)       % upper bound on the NLOS bias
    b_verbose = 0
end

methods
    function [value, v_r] = computeValue(obj, v_x, v_d_in)
        % Given a point in space, range difference measurements, and
        % source positions, computes the objective value
        m_s = obj.Xenb;
        N = size(m_s, 2);
        rho = obj.param_rho;
        assert(size(m_s, 1)==2);
        assert(length(v_d_in)==N-1);
        v_d1 = [0; v_d_in(:)];
        m_D = v_d1' - v_d1;
        v_r = norms(v_x - m_s);
        
        v_valueForRefSource = nan(N, 1);
        for i_s = 1:N % index pointing at the reference source
            value_0 = 0;
            value_rho = 0;
            for k = 1:N % index pointing at the source whose TDOA we add
                if k~=i_s
                    expr_0 = v_r(i_s) - v_r(k) - m_D(k, i_s);
                    expr_rho = expr_0 - rho;
                    value_0   = value_0   + expr_0^2   + 2*rho*max(0, rho+expr_0);
                    value_rho = value_rho + expr_rho^2 + 2*rho*max(0, rho+expr_rho);
                end
            end
            v_valueForRefSource(i_s) = max(value_0, value_rho);
        end
        value = mean(v_valueForRefSource);
    end
    
    function m_map = valueMap(obj, m_grid_x1, m_grid_x2, v_d)
        % Given a grid of points as produced by meshgrid, returns a matrix 
        % containing the objective value for each point in the grid.
        
        % This function simply iterates over every point in the grid, calls
        % the method calculate_value and stores the output in a matrix.
        assert(isequal(size(m_grid_x1), size(m_grid_x2)), ...
            'size of grids not coherent')
        m_map = nan(size(m_grid_x1));
        m_s = obj.Xenb;

        ltc = LoopTimeControl(numel(m_map));
        for n = 1:numel(m_map)
            m_map(n) = obj.computeValue(...
                [m_grid_x1(n);m_grid_x2(n)], v_d, m_s);
            ltc.go(n);
        end
        
    end
    
    function [v_estimatedLocation, locUncertainty] = ...
            estimateOneLocation(obj, v_measurements)
        v_estimatedLocation = obj.solveRelaxed(v_measurements);
        locUncertainty = nan; % TODO! try to measure location uncertainty
    end
    
    function [v_x_hat, v_r_out, value] = solveRelaxed(obj, ...
                v_d_in)
        % Solve the convex relaxation of minimizing the cost function 
        % implemented by the method computeValue
                
        m_s = obj.Xenb;
    
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
            if obj.b_verbose, disp 'Conforming constraints...'; end
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

            if obj.b_verbose, disp 'Solving relaxed problem...'; end
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