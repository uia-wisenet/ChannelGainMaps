classdef HybridEstimator < Estimator
    % HybridEstimator class
    %   This class allows to train and query a LF and LB kernel
    %   machine desgined to estimate the channel gain from features of both types 
    %   (location-free features annd location features) received at sensors pairs.
    
    properties
        h_kernelLF % function handle defining the LF kernel
        h_kernelLB % function handle defining the LF kernel
        regularizationParameterLF % lambda for LocFree machine  (KRR)
        regularizationParameterLB % lambda for LocBased machine (KRR)
        max_itrs_alternating =20;
        b_verbose = 0;
        b_debugPlots = 0;
        stop_tol = 1e-3
        
    end
    
    methods
        function [v_coefficients_LF, v_coefficients_LB, intercept, wcu] ...
                = train(obj, t_locFeatures, t_estimatedLocation, ...
                m_locErrors, v_channelForPairs)
            % Train hybrid Channel-gain estimator.
            
            % GIVEN:    locF features, estimated locations, and channel
            %           gains
            % OPTIMIZE: the coefficients of the kernel machines used by the
            %           hybrid estimator of the channel gain (as a linear 
            %           combination of the outputs of the latter machines)
            
            % Inputs: 
            % t_locFeatures: MxNx2-tensor containing locF features
            %               (where M is the number of features at one
            %               sensor location, N is the number of pairs 
            %               for training, and 2 stands for pair)
            % t_estimatedLocation: DxNx2-tensor containing the estimated
            %               locations (from a LocationEstimator) for each
            %               training pair -- D is space dim: either 2 or 3
            % t_locErrors:  Nx2-matrix containing the uncertainty (error)
            %               measure for the estimated location
            % v_channelForPairs: N-vector containing channel gain for each 
            %               training pair
            
            % Outputs: 
            % v_coefficients_LF: N-vector containing the LocFree
            %               coefs resulting from training 
            % v_coefficients_LB: N-vector containing the LocBased
            %               coefs resulting from training
            % intercept: scalar intercept (mean of training channel gains)
            % v_d_out:   N-vector containing coefficients at the diagonal 
            %             of matrix D (optimal weights for the linear 
            %             combination at the training points 
            
            %%
            % Dimensionalities and coherence check
            assert(iscolumn(v_channelForPairs));
            n_ues= length(v_channelForPairs); % number of (training) UEs
            d = size(t_estimatedLocation, 1); % dim. of location vectors
            assert(d==2 || d==3);             % 2D or 3D location
            assert(size(t_estimatedLocation, 2)==n_ues);
            assert(size(t_estimatedLocation, 3)==2);%2 stands for pair
            assert(size(t_locFeatures, 2)==n_ues);
            assert(size(t_locFeatures, 3)==2); %2 stands for pair
            assert(isequal(size(m_locErrors), [n_ues, 2]));
            
            disp('Building kernel matrices row by row')
            m_K_p = obj.buildKernelMatrixbyRows(...
                obj.h_kernelLF, t_locFeatures,       obj.b_verbose);
            m_K_l = obj.buildKernelMatrixbyRows(...
                obj.h_kernelLB, t_estimatedLocation, obj.b_verbose);
             %% section that will be commented out
%             m_Ke_p = zeros(n_ues, n_ues);
%             m_Ke_l = zeros(n_ues, n_ues);           
%             ltc = LoopTimeControl(n_ues);
%             for ind_ue = 1:n_ues
%                 m_row_inputs_to_kernels = repmat([t_locFeatures(:, ind_ue, 1);...
%                     t_locFeatures(:,ind_ue,2)], [1 n_ues]);
%                 m_row_inputs_to_kernels_2 = repmat([t_estimatedLocation(:, ind_ue, 1);...
%                     t_estimatedLocation(:,ind_ue,2)], [1 n_ues]);
%                 
%                 v_my_row = feval(obj.h_kernelLF, m_row_inputs_to_kernels, ...
%                     [t_locFeatures(:,:,1);t_locFeatures(:,:,2)]);
%                 v_my_row_2 = feval(obj.h_kernelLB, m_row_inputs_to_kernels_2, ...
%                     [t_estimatedLocation(:,:,1);t_estimatedLocation(:,:,2)]);
%                 m_Ke_p(ind_ue,:) = v_my_row;
%                 m_Ke_l(ind_ue,:) = v_my_row_2;
%                 ltc.go(ind_ue);
%             end
%             norm(m_K_p(:)-m_Ke_p(:))
%             norm(m_K_l(:)-m_Ke_l(:))
%             keyboard
             %%
            intercept = mean(v_channelForPairs);
            v_zm_channelForPairs=v_channelForPairs(:)-intercept; % zero mean
            
            %% Alternating minimization 
            % (matrix D and the kernel weights in \gamma: see document sent to Seung Jun) 
            [v_locErrors_tmp, v_indices] = sort(m_locErrors(:));
            %v_h_norm_dvect=zeros(1,obj.max_itrs_alternating);
            %v_h_norm_dvecr=zeros(1,obj.max_itrs_alternating);
            v_h_obj = nan(1, obj.max_itrs_alternating);
            v_d = linspace(0, 1, n_ues);
            m_K_large = sparse(blkdiag(m_K_p, m_K_l));
            lambda_p = obj.regularizationParameterLF;
            lambda_l = obj.regularizationParameterLB;
            m_Lambda = sparse(blkdiag(lambda_p*eye(n_ues), lambda_l*eye(n_ues))); 
                     
            ltc = LoopTimeControl(obj.max_itrs_alternating);
            for ind_iters=1:obj.max_itrs_alternating 
%                 v_h_norm_dvect(ind_iters)= norm(v_d_t);
%                 v_h_norm_dvecr(ind_iters)= norm(v_d_r); %are these
%                 needed?

                % Minimization 1: Given v_d_diag, solve for 
                % kernel machine coefficients
                m_D = diag(v_d);
                m_fat = [eye(n_ues)-m_D,  m_D];
                m_toInvert = (m_fat'*m_fat)*m_K_large + n_ues*m_Lambda;
                v_gamma = m_toInvert \ ...
                    ( m_fat'* v_zm_channelForPairs);
                v_coefficients_LF = v_gamma(1:n_ues);
                v_coefficients_LB = v_gamma(n_ues+1:end);
                
                %% commented out (making sure the minimization is correct)
%                 cvx_begin quiet
%                   variables alpha_p(n_ues) alpha_l(n_ues)
%                   minimize((1/n_ues) * sum_square( v_zm_channelForPairs ...
%                       - (eye(n_ues)-m_D)*m_K_p*alpha_p ...
%                       - m_D *m_K_l * alpha_l )...
%                       + lambda_p*alpha_p'*m_K_p*alpha_p...
%                       + lambda_l*alpha_l'*m_K_l*alpha_l)
%                 cvx_end
%                 norm([alpha_p;alpha_l]-v_gamma)/norm(v_gamma)
%                 keyboard
                
                % Minimization 2: Given v_a and v_b, solve for v_d
                v_a = v_zm_channelForPairs- m_K_p*v_coefficients_LF;
                v_b = m_K_l*v_coefficients_LB - m_K_p*v_coefficients_LF;
                [m_dt_dr, obj_value] = obj.fit_D_additive(v_a, v_b, ...
                    v_indices);
                v_d = m_dt_dr*[1;1]/2;
                
                %% section that will be commented out
%                 % d_vec = diag(b_vec)\a_vec; %?
%                 % solve using CVX
%                 cvx_begin quiet
%                     variable v_d_t2(n_ues, 1);
%                     variable v_d_r2(n_ues, 1);
%                     minimize(norm(v_a - (diag(v_b) * (v_d_t2 + v_d_r2))));
%                     subject to
%                         0 <= v_d_t2 + v_d_r2;
%                         v_d_t2 + v_d_r2 <= 1;
%                         v_d_t2(indices_t(2:end)) >= v_d_t2(indices_t(1:end-1));
%                         v_d_r2(indices_r(2:end)) >= v_d_r2(indices_r(1:end-1));        
%                 cvx_end
%                 norm(v_d_t-v_d_t2)
%                 norm(v_d_r-v_d_r2)
%                 keyboard
                %%
                
                v_h_obj(ind_iters) = obj_value.^2/n_ues ...
                    + lambda_p* v_coefficients_LF'*m_K_p*v_coefficients_LF...
                    + lambda_l* v_coefficients_LB'*m_K_l*v_coefficients_LB;
                
                if obj.b_debugPlots    % TODO: refactor as method
                    figure(990); clf
                    subplot(1, 3, 1);
                    stem(v_gamma); 
                    title 'coefficients'
                    xlabel 'sample index'
                    ylabel \gamma
                    subplot(1, 3, 2);
                    plot(v_h_obj); 
                    xlabel 'iteration index' 
                    ylabel 'objective value'
                    title 'progress'
                    xlim([0 obj.max_itrs_alternating])
                    subplot(1, 3, 3)
                    plot(m_locErrors(v_indices), m_dt_dr(v_indices))   
                    title 'weighting function'
                    xlabel 'e'
                    ylabel '$\bar{g}$ (e)' interpreter latex
                end
                
                %check stopping criterion
                if ind_iters > 1 && abs( ...
                        v_h_obj(ind_iters) - v_h_obj(ind_iters-1) ...
                        ) < obj.stop_tol
                    break
                end
                
                
                ltc.go(ind_iters);
            end
            v_weights_tmp = m_dt_dr(v_indices);
            wcu = AdditiveWeightCalculatingUnit;
            % leave only unique entries:            
            [wcu.v_locErrors, v_representers] = unique(v_locErrors_tmp);
            wcu.v_weights = v_weights_tmp(v_representers);
%%            
%             plots (should be traced outside of the processing block)
%             close all
%             figure(4)
%             plot(v_h_norm_dvect+v_h_norm_dvecr,'b-o', 'LineWidth',2);   grid on
%             xlabel('Iteration index k',  'Interpreter', 'latex')
%             ylabel('\boldmath $\vert\! \vert$ $\! d_t + d_r $ $\!\vert\!\vert$', 'Interpreter', 'latex')
            
        end
        
        function predictedMeasurements = estimateGivenEstimatedDistancesAndLoc(obj, ...
                m_coefs, t_extractedLocfreeFeaturesToConsider, t_trainingLocations, v_trainingMeasurements,...
                v_funct_weights, test_locF_features, estimatedLocationVal)
            % Estimate channel gain using Hybrid method (locF + locB
            % cartography)
            
            % Inputs: (where M = no. locF features per sensor location, 
            %                N = no. pairs for evaluation;
            %                N_train = no. pairs at training 
            %                d = dimensionality of location space)
            % m_coefs:              N_train x 2- matrix containing 
            %                       coefficients gamma=[alpha, beta]
            % t_extractedLocfreeFeaturesToConsider: 
            %                       MxN_trainx2- tensor containing the
            %                       locFree features of the TRAINING data; 
            %                       2 stands for pair
            % t_trainingLocations:  dxN_trainx2- tensor, estimated locations
            %                       of the TRAINING data
            % v_trainingMeasurements: N_train- vector, measured channel gains of
            %                       TRAINING data
            % v_funct_weights:      N-vector with weights (externally
            %                       calculated for the query/test points)
            % test_locF_features:   MxNx2- array, locF features of query/test points
            % estimatedLocationVal  dxNx2- array, estimated locations of query/test points
            
            % Outputs:
            % predictedMeasurements: N-vector with predicted channel gains
            %                       for query/test points
            
            b_hybrid_option = true; 
            ed_dim = size(test_locF_features); % M x N x 2
            n_uesPair_tr = length(v_trainingMeasurements); % N_train
            predictedMeasurements = zeros(1,ed_dim(2));
            d = size(t_trainingLocations, 1);
            for ind_ue_val = 1:ed_dim(2)
                v_dist_to_source_input_to_ker_tx = ... %(such a long name could be shortened)
                    [test_locF_features(:,ind_ue_val,1); test_locF_features(:,ind_ue_val,2)];
                v_queryLocations=[estimatedLocationVal(:,ind_ue_val,1);estimatedLocationVal(:,ind_ue_val,2)];
                m_features_LF = zeros(ed_dim(1)*2, n_uesPair_tr); %2*M x N_train
                m_features_LF(1:ed_dim(1), :)     = t_extractedLocfreeFeaturesToConsider(:,:,1);
                m_features_LF(ed_dim(1)+1:end, :) = t_extractedLocfreeFeaturesToConsider(:,:,2);
                m_features_LB = zeros(d*2, n_uesPair_tr); %2*d x N_train
                m_features_LB(1:d, :)    = t_trainingLocations(:,:,1);
                m_features_LB(d+1:end,:) = t_trainingLocations(:,:,2);
                row_kernel_LF = feval(obj.h_kernelLF, v_dist_to_source_input_to_ker_tx,  m_features_LF);
                row_kernel_LB = feval(obj.h_kernelLB, v_queryLocations,                  m_features_LB);
                if b_hybrid_option==true
                    predictedMeasurements(ind_ue_val)=((1-v_funct_weights').*row_kernel_LF)*m_coefs(:,1)+...
                        (v_funct_weights'.*row_kernel_LB)*m_coefs(:,2)+mean(v_trainingMeasurements);
                else %if b_hybrid_option is false, this is just an additive model
                    predictedMeasurements(ind_ue_val)=row_kernel_LF*m_coefs(:,1)+...
                        row_kernel_LB*m_coefs(:,2)+mean(v_trainingMeasurements);
                end
                
            end
        end
    end %methods
    
    methods(Static)
        function [m_K] = buildKernelMatrixbyRows(...
                h_kernel, t_features, b_verbose)
            % Build kernel matrix, row by row.
            % Given: 
            %   - a function handle expressing the kernel function,
            %   - an MxNx2 tensor with training feature vectors
            %   - an optional verbosity flag
            % Calculate:
            %   - an NxN kernel matrix
            if ~exist('b_verbose', 'var')
                b_verbose = 0;
            end
            if b_verbose
                disp('Building kernel matrices row by row')
            end
           
            n_ues = size(t_features, 2);
            assert(size(t_features, 3) == 2);
            m_K = zeros(n_ues, n_ues); % Kernel matrix
            
            m_allFeatures = [t_features(:,:,1); t_features(:,:,2)];
            
            ltc = LoopTimeControl(n_ues);
            for ind_ue = 1:n_ues
                m_row_inputs_to_kernels = repmat(...
                    [t_features(:, ind_ue, 1); t_features(:,ind_ue,2)], ...
                        [1 n_ues]);                
                v_my_row = feval(h_kernel, m_row_inputs_to_kernels, ...
                    m_allFeatures);               
                m_K(ind_ue,:) = v_my_row;
                if b_verbose, ltc.go(ind_ue); end
            end
        end
        
        function [m_dt_dr, obj_value] ...
                = fit_D_additive(v_a, v_b, v_indices, m_d_initial)
            % Given v_a and v_b, solve for v_d
            
            % d_vec = diag(b_vec)\a_vec; %?
           
            
            n_ues = length(v_a);
            assert(length(v_b)==n_ues);
            assert(length(v_indices)==2*n_ues);
            m_dt_dr = zeros(n_ues, 2);
            
            %solve using Quadprog
            m_I = eye(2*n_ues);
            m_C = zeros(2*n_ues-1, 2*n_ues);
            for ii = 1:2*n_ues-1
                m_C(ii,:) = m_I(v_indices(ii),:) - m_I(v_indices(ii+1),:);
            end
            m_I_n = eye(n_ues);
            m_H = [m_I_n; m_I_n]*diag(v_b.^2)*[m_I_n m_I_n]/2;
            v_f = -v_a'*diag(v_b)*[m_I_n m_I_n];
            m_A = [-m_I; m_I; -m_C];
            v_bq= [zeros(2*n_ues,1); ones(2*n_ues,1); zeros(2*n_ues-1, 1)];
            optimoptions = optimset('Display', 'off');
            if ~exist('m_d_initial', 'var'), m_d_initial = []; end
            [m_dt_dr(:), fval, flag_exit] = quadprog(...
                m_H, v_f, m_A, v_bq, [], [],[],[],m_d_initial(:), optimoptions);
            obj_value = sqrt(0.5*m_dt_dr(:)'*m_H*m_dt_dr(:)+v_f*m_dt_dr(:) + v_a'*v_a);
%           obj_value == sqrt(fval + v_a'*v_a)
%           obj_value == norm( v_a - (diag(v_b) * m_dt_dr*[1; 1]/2) )
%           should hold
            
%             % solve using CVX
%             cvx_begin quiet
%               variable m_D(n_ues, 2)
%               minimize(  norm( v_a - (diag(v_b) * m_D*[1; 1]/2) )  );
%               subject to
%               0 <= m_D; %#ok<VUNUS>
%               m_D <= 1; %#ok<VUNUS>
%               m_D(v_indices(2:end)) <= m_D(v_indices(1:end-1)); %#ok<VUNUS>
%             cvx_end
%             obj_value = cvx_optval;
%                         
%             norm(m_D(:) - m_dt_dr(:))
%             keyboard
             
        end
        
    end %methods
end



