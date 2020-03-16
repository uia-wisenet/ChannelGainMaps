classdef WangLocationEstimator
    properties
        
    end
    methods
        function[x_out, v_r_out] = solve_robust_loc(obj, ...
                v_d, m_s, rho)
            % Solves the SDP Relaxation of the Robust Location Estimation
            % by Wang et al.
            % Inputs:
            % v_d: (N-1)-vector (N=number of sources) containing the
            %      measured TDOAs (assumed corrupted with small Gaussian 
            %      noise and big NLOS bias)
            % m_s: 2xN-matrix containing the source positions. The
            %      reference source is in column 1.
            % rho: upper bound on the NLOS bias
            
            
            N = length(v_d)+1;
            b = zeros(N-1, 1);
            a_= zeros(N+1, N-1);
            for i = 1:N-1
                b(i) = -v_d(i)^2-norm(m_s(:,i+1))^2+norm(m_s(:,1))^2;
                a_(:,i)=[ 2*(m_s(:,1)-m_s(:,i+1));
                          zeros(i-1,1);
                          -2*v_d(i);
                          zeros(N-i-1,1)   ];
            end
            y = sdp_ce(m_s, v_d, b, N, rho, a_);
            x_out = y(1:2);
            v_r_out = y(3:end);
        end
    end
end