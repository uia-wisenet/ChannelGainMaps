classdef SCALocationEstimator
    properties
        alpha = 100;
        beta  = 1;
        max_steps_sca = 6;
    end
    
    methods
        
        function [m_x_out, v_r_out, value] = solve_sca_loc(obj, ...
                v_d, m_s, rho, e_0, x_0, r_0)
            m_x_out = zeros(2, obj.max_steps_sca);
            m_r_out = zeros(length(r_0), obj.max_steps_sca);
            m_x_out(:,1)= x_0;
            m_r_out(:,1)= r_0;
            
            ltc = LoopTimeControl(obj.max_steps_sca);
            for k = 1:obj.max_steps_sca-1
                [m_x_out(:,k+1), v_r_out, value] = solve_convexified_loc(obj, ...
                    v_d, m_s, rho, e_0, m_x_out(:,k), m_r_out(:,k));
                ltc.go(k);
            end
        end
        
        function [x_out, v_r_out, value] = solve_convexified_loc(obj, ...
                v_d, m_s, rho, e_0, x_0, r_0)
            N = length(v_d) +1;
            
            cvx_begin quiet
                variable x(2)
                variable v_r(N)
                objective = 1/obj.alpha*sum_square(x-x_0) ...
                    + 1/obj.beta*sum_square(v_r-r_0);
                for n = 2:N
                    my_exp = v_d(n-1) - v_r(n) + v_r(1) + e_0;
                    objective = objective + ...
                        my_exp^2 + rho*max(0, rho-2*my_exp);                   
                end
                %constraints:
                for n = 1:N
                    norm(x-m_s(:,n)) <= v_r(n);
                    (x_0-m_s(:,n)).'*(x-m_s(:,n)) >= norm(x_0-m_s(:,n))*v_r(n); 
                    %linearized the concave constraint
                end
                %
                minimize(objective)
            cvx_end
            x_out = x;
            v_r_out = v_r;
            value = cvx_optval;
        end
        
    end
end