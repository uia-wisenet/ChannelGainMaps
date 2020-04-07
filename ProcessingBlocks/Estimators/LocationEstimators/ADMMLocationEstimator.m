classdef myLocationEstimator
    properties
        rho_admm = 0.1;
        lambda = 1e-6;
        max_steps_admm = 1000;
        max_steps_grad = 1000;
        tol = 0.01;
        alpha = 0.01;
        
        b_useHandles = 1;
    end
    
    methods
        
        function [c_handles] = create_handles (obj, ...
                d, e_0, rho_wc, s_i, s_0)
            x = sym('x', [2 1]);
            y = sym('y', [2 1]);
            x_bar = sym('x_bar', [2 1]);
            my_exp(x) = d - norm(x-s_i) + norm(x-s_0) + e_0;
            grad0([x; y; x_bar]) = gradient(y.'*(x-x_bar)...
                +obj.rho_admm*((x(1)-x_bar(1))^2 + (x(2)-x_bar(2))^2), x);
            grad1(x)     = gradient(my_exp^2, x);
            grad2(x)     = gradient(rho_wc*(rho_wc-2*my_exp), x);
            c_handles{1} = matlabFunction(grad0);
            c_handles{2} = matlabFunction(grad1);
            c_handles{3} = matlabFunction(grad2);
            c_handles{4} = matlabFunction(my_exp);
        end
        
        function x_out = solve_nonconvex_onesource_loc_handle(obj, ...
                c_handles, rho_wc, v_y, x_bar)
            x_now = x_bar;
            for k = 1:obj.max_steps_grad
                c_args1 = num2cell([x_now; v_y; x_bar]);
                c_x_now = num2cell(x_now);
                my_grad = feval(c_handles{1}, c_args1{:}) ...
                    + feval(c_handles{2}, c_x_now{:});
                if rho_wc > 2*feval(c_handles{4}, c_x_now{:})
                    my_grad = my_grad + feval(c_handles{3}, c_x_now{:});
                end
                x_now = x_now - obj.alpha*my_grad;
            end
            x_out = x_now;
        end
        
        function x_out = solve_nonconvex_onesource_loc_subs(obj, ...
                d, e_0, rho_wc,v_y, x_bar, s_i, s_0)
            x = sym('x', [2 1]);
            my_exp = d - norm(x-s_i) + norm(x-s_0) + e_0;
            grad0 = gradient(v_y'*(x-x_bar)...
                +obj.rho_admm*((x(1)-x_bar(1))^2 + (x(2)-x_bar(2))^2), x);
            grad1 = gradient(my_exp^2, x);
            grad2 = gradient(rho_wc*(rho_wc-2*my_exp), x);
            x_now = x_bar;
            for k = 1:obj.max_steps_grad
                my_grad = subs(grad0,  x, x_now);
                my_grad = my_grad + subs(grad1, x, x_now);
                if rho_wc > 2*subs(my_exp, x, x_now)
                    my_grad = my_grad + subs(grad2, x, x_now);
                end
                x_now = x_now - obj.alpha*eval(my_grad);
            end
            x_out = x_now;
        end
        
        
        function [x_out, b_istight] = solve_relaxed_onesource_loc(obj, ...
                d, e_0, rho_wc, v_y, x_bar, s_i, s_0)
            cvx_begin quiet
                variable r_i
                variable r_0
                variable x(2)
                my_exp = d - r_i + r_0 + e_0;
                minimize( my_exp^2 + rho_wc*(max(0, rho_wc-2*my_exp)) ...
                    + v_y'*(x-x_bar) + obj.rho_admm*sum_square(x-x_bar)...
                    + obj.lambda*norm([r_i;r_0]));
                norm(x-s_0) <= r_0 %#ok<NOPRT>
                norm(x-s_i) <= r_i %#ok<NOPRT>
            cvx_end
            x_out = x;
            b_istight = abs(norm(x-s_0)- r_0) + abs(norm(x-s_i) - r_i) < obj.tol;
        end
        
        function [x_out] = solve_onesource_loc(obj, ...
                d, e_0, rho_wc, v_y, x_bar, s_i, s_0)
            [x_out, b_istight] = solve_relaxed_onesource_loc(obj, ...
                d, e_0, rho_wc, v_y, x_bar, s_i, s_0);
            if b_istight
                return
            else
                x_out = solve_nonconvex_onesource_loc(obj, ...
                d, e_0, rho_wc,v_y, x_bar, s_i, s_0);
            end
        end
        
        function x_out = solve_admm_loc(obj, ...
                v_d, m_s, rho_wc, e_0, x_0)
            
            N = size(m_s,2);
            m_x_hats = zeros(2, N-1);
            m_y      = zeros(2, N-1);
            x_bar    = zeros(2, obj.max_steps_admm);
            x_bar(:, 1) = x_0;
            
            if obj.b_useHandles
                disp 'preparing gradients...'
                if exist('saved_handles.mat', 'file')
                    load saved_handles.mat
                    assert(exist('c_handles', 'var') && ...
                        all(size(c_handles) == [4 N-1]));
                else
                    for i = 1:N-1
                        c_handles(:, i) = obj.create_handles( ...
                            v_d(i), e_0, rho_wc, m_s(:, i+1), m_s(:,1));
                    end
                    save saved_handles.mat c_handles
                end
            end
            
            ltc = LoopTimeControl(obj.max_steps_admm);
            for k = 1:obj.max_steps_admm
                %fprintf('Iteration %d', k);
                if obj.b_useHandles
                    for i = 1:N-1 %parfor
                        m_x_hats(:,i) = obj.solve_nonconvex_onesource_loc_handle(...
                            c_handles(:,i), rho_wc, m_y(:,i), x_bar(:,k));
                    end
                else
                    for i = 1:N-1 %parfor
                        m_x_hats(:,i) = obj.solve_nonconvex_onesource_loc(...
                            v_d(i), e_0, rho_wc, m_y(:,i), x_bar(:,k), ...
                            m_s(:,i+1), m_s(:,1)); %#ok<PFBNS>
                    end
                end
                x_bar(:,k+1) = mean(m_x_hats, 2);
                for i = 1:N-1
                    m_y(:,i) = m_y(:,i) + obj.rho_admm*(...
                        m_x_hats(:, i)-x_bar(:,k+1));
                end
                ltc.go(k);
            end
            x_out = x_bar;
        end
            
    end
end