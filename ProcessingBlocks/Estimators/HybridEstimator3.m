classdef HybridEstimator3 < HybridEstimator
    % HybridEstimator class
    %   This class allows to train and query a LF and LB kernel
    %   machine desgined to estimate the channel gain from features of both types 
    %   (location-free features annd location features) received at sensors pairs.
    
%     properties
%         h_kernelLF % function handle defining the LF kernel
%         h_kernelLB % function handle defining the LF kernel
%         regularizationParameterLF % lambda for LocFree machine  (KRR)
%         regularizationParameterLB % lambda for LocBased machine (KRR)
%         max_itrs_alternating =30;
%         b_verbose = 1;
%         b_debugPlots = 0;
%         stop_tol = 1e-4
%         method_fit_D = 'twoDimensional'
%         
%         b_tryToBalance = 0;
%         
%     end

    properties
        h_kernelG % function handle defining the kernel for the gating function
        regularizationParameterG % lambda for gating function (Kernel-based)
        %kernelPar_G = 5 % parameter nu for the gating function kernel
        pgd_stepsize = 1e-2;
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
            
            %disp('Building kernel matrices row by row')
            m_K_p = obj.buildKernelMatrixbyRows(...
                obj.h_kernelLF, t_locFeatures,       obj.b_verbose);
            m_K_l = obj.buildKernelMatrixbyRows(...
                obj.h_kernelLB, t_estimatedLocation, obj.b_verbose);
            m_K_g = obj.buildKernelMatrixbyRows(...
                obj.h_kernelG, reshape(m_locErrors, [1 n_ues 2]), obj.b_verbose);
            
            intercept = mean(v_channelForPairs);
            v_zm_channelForPairs=v_channelForPairs(:)-intercept; % zero mean
            h_proj_simplex_vector = @(y) ...
                max(y-max((cumsum(sort(y,1,'descend'),1)-1)./(1:size(y,1))'),0);

            %% Alternating minimization 
            % (matrix D and the kernel weights in \gamma: see document sent to Seung Jun) 
            m_locErrors_wcu = [max(m_locErrors); m_locErrors; min(m_locErrors)];           
            [v_locErrors_tmp, v_indicesOrderErrors] = sort(m_locErrors(:));
            [m_indicesAllOrderPairs, m_indicesDAG_Pairs] = ...
                obj.build_DAG(m_locErrors);
            [m_indicesAllOrder_wcu, m_indicesDAG_wcu] = ...
                obj.build_DAG(m_locErrors_wcu);
            
            v_h_obj = nan(1, obj.max_itrs_alternating);
            m_K_large = sparse(blkdiag(m_K_p, m_K_l));
            lambda_p = obj.regularizationParameterLF;
            lambda_l = obj.regularizationParameterLB;
            lambda_g = obj.regularizationParameterG;
            
%             %initialize coefficients to the pure LocB and the pure LocF
%             % kernel machines
            v_coefficients_LF_pure = (m_K_p+lambda_p*eye(n_ues))\v_zm_channelForPairs;
            v_coefficients_LB_pure = (m_K_l+lambda_l*eye(n_ues))\v_zm_channelForPairs;
             
            v_a = v_zm_channelForPairs- m_K_p*v_coefficients_LF_pure;
            v_b = m_K_l*v_coefficients_LB_pure - m_K_p*v_coefficients_LF_pure;
            cvx_begin
                variable v_c_G(n_ues)
                minimize (sum_square(v_a - diag(v_b)*m_K_g*v_c_G))%! + lambda_g* v_c_G'*m_K_g*v_c_G)
                subject to
                v_c_G >= 0;         %#ok<VUNUS>
                sum(v_c_G) == 1;    %#ok<EQEFF>
            cvx_end
            v_d = m_K_g*v_c_G;
            
            figure(989)
            stem3(m_locErrors(:, 1), m_locErrors(:,2), v_d)
            xlabel('e_r'), ylabel('e_t')
            keyboard
            
            v_coefficients_g = v_c_G;
            
            % initialize v_d to 1/2
            % v_d = 1/2*ones(n_ues, 1);
            
            % initialize v_coefficients_G instead
%             v_coefficients_g = ones(n_ues, 1)./n_ues;
%             v_d = m_K_g*v_coefficients_g;
            m_D = diag(v_d);
            
            ltc = LoopTimeControl(obj.max_itrs_alternating);
            for ind_iters=1:obj.max_itrs_alternating 
                % Minimization 1: Given v_d_diag, solve for 
                % kernel machine coefficients
                m_fat = sparse([eye(n_ues)-m_D,  m_D]);
                m_Lambda = sparse(blkdiag(...
                    lambda_p*eye(n_ues), lambda_l*eye(n_ues)));
                if obj.b_tryToBalance
                    m_Lambda = sparse(blkdiag(...
                        lambda_p.*mean(1-v_d)*eye(n_ues),...
                        lambda_l.*mean(  v_d)*eye(n_ues) ));
                end
                m_toInvert = m_fat'*(m_fat*m_K_large) + n_ues*m_Lambda;
                v_gamma = m_toInvert \ ...
                    ( m_fat'* v_zm_channelForPairs);
                v_coefficients_LF = v_gamma(1:n_ues);
                v_coefficients_LB = v_gamma(n_ues+1:end);

                % Minimization 2: Given v_a and v_b, solve for v_d
                v_a = v_zm_channelForPairs- m_K_p*v_coefficients_LF;
                v_b = m_K_l*v_coefficients_LB - m_K_p*v_coefficients_LF;
                assert( isequal(obj.method_fit_D, 'kernel'))
                
                if obj.b_tryToBalance
                    regDiff = lambda_l * v_coefficients_LB'*...
                        m_K_l * v_coefficients_LB...
                        -  lambda_p * v_coefficients_LF'*...
                        m_K_p*v_coefficients_LF;
                else
                    regDiff = 0;
                end
                %
                % perform a step of projected gradient descent
                v_grad = (m_K_g'*diag(v_b.^2)+eye(n_ues)*lambda_g)*m_K_g * v_coefficients_g ...
                    - m_K_g'*(v_b.*v_a - regDiff);
                v_gradIterate = v_coefficients_g - obj.pgd_stepsize*v_grad;
                v_coefficients_g = h_proj_simplex_vector(v_gradIterate);
                
                
                %
                %[v_d, obj_value] = obj.fit_D_twoDimensional(...
                %    v_a, v_b, regDiff, m_indicesDAG_Pairs, v_d);
                
                v_d = m_K_g*v_coefficients_g;
                m_D = diag(v_d);
                                
                %%
                
                v_h_obj(ind_iters) = sum_square( v_zm_channelForPairs...
                    - sparse(m_D)           *(m_K_l*v_coefficients_LB) ...
                    - sparse(eye(n_ues)-m_D)*(m_K_p*v_coefficients_LF) )/n_ues ...
                    + lambda_p* mean(1-v_d)*v_coefficients_LF'*m_K_p*v_coefficients_LF...
                    + lambda_l* mean(  v_d)*v_coefficients_LB'*m_K_l*v_coefficients_LB;
                
                if obj.b_debugPlots
                    switch obj.method_fit_D
                        case 'additive'
                            obj.plotProgress(v_gamma, v_h_obj, ...
                                m_locErrors, m_dt_dr, v_indices);
                        case 'twoDimensional'
                            obj.plotProgress(v_gamma, v_h_obj, ...
                                m_locErrors, v_d);
                        case 'kernel'
                            obj.plotProgress(v_gamma, v_h_obj, ...
                                m_locErrors, v_d);
                            % TODO: perhaps show a mesh instead of a stem3
                    end
                    keyboard
                end
                
                %check stopping criterion
                if ind_iters > 1 && abs( ...
                        v_h_obj(ind_iters) - v_h_obj(ind_iters-1) ...
                        ) < obj.stop_tol
                    break
                end
                
                if obj.b_verbose
                    ltc.go(ind_iters);
                end
            end
            switch obj.method_fit_D
                case 'additive'
                    v_weights_tmp = m_dt_dr(v_indices);
                    wcu = AdditiveWeightCalculatingUnit;
                    % leave only unique entries:
                    [wcu.v_locErrors, v_representers] = unique(v_locErrors_tmp);
                    wcu.v_weights = v_weights_tmp(v_representers);
                case 'twoDimensional'
                    wcu = TwoDimensionalWeightCalculatingUnit;
                    wcu.m_locErrors_train = m_locErrors_wcu;
                    wcu.m_indicesAllOrder = m_indicesAllOrder_wcu;
                    wcu.m_indices_DAG = m_indicesDAG_wcu;
                    wcu.v_weights = [min(v_d); v_d; max(v_d)];
                case 'kernel'
                    wcu = KernelBasedWeightCalculatingUnit;
                    wcu.m_trainLocErrors = m_locErrors;
                    wcu.v_coefficients = v_coefficients_g;
                otherwise
                    error 'unrecognized curve fitting method'
            end
        end
        
        function [] = plotProgress(obj,v_gamma, v_h_obj, m_locErrors, ...
                m_dt_dr__or__v_d, v_indices)
            figure(990); clf
            subplot(1, 3, 1);
              stem(v_gamma);
              title 'coefficients'
              xlabel 'sample index'
              ylabel \gamma
            subplot(1, 3, 2);
              plot(v_h_obj, 's-');
              y_top = max(v_h_obj);
              y_down = [y_top./1.4.^(1:8) 0];
              y_bottom = max(y_down(y_down<min(v_h_obj)));
              xlabel 'iteration index'
              ylabel 'objective value'
              title 'progress'
              xlim([0 obj.max_itrs_alternating])
              ylim([y_bottom y_top])
            subplot(1, 3, 3)
              switch obj.method_fit_D
                  case 'additive'
                      m_dt_dr = m_dt_dr__or__v_d;
                      plot(m_locErrors(v_indices), m_dt_dr(v_indices))
                      xlabel 'e'
                      ylabel '$\bar{g}$ (e)' interpreter latex
                  case 'twoDimensional'
                      v_d = m_dt_dr__or__v_d;
                      stem3(m_locErrors(:, 1), m_locErrors(:,2), ...
                          v_d)
                      xlabel('e_r')
                      ylabel('e_t')
                      zlabel ('$\bar{g}$ (e)', 'interpreter', 'latex')
                  case 'kernel' % copied from twoDimensional
                      v_d = m_dt_dr__or__v_d;
                      stem3(m_locErrors(:, 1), m_locErrors(:,2), ...
                          v_d)
                      xlabel('e_r')
                      ylabel('e_t')
                      zlabel ('$\bar{g}$ (e)', 'interpreter', 'latex')
              end
            title 'gating function'
        end
        
%         function predictedMeasurements = estimateGivenEstimatedDistancesAndLoc(obj, ...
%                 m_coefs, t_extractedLocfreeFeaturesToConsider, t_trainingLocations, v_trainingMeasurements,...
%                 v_funct_weights, test_locF_features, estimatedLocationVal)
%             % Estimate channel gain using Hybrid method (locF + locB
%             % cartography)
%             
%             % Inputs: (where M = no. locF features per sensor location, 
%             %                N = no. pairs for evaluation;
%             %                N_train = no. pairs at training 
%             %                d = dimensionality of location space)
%             % m_coefs:              N_train x 2- matrix containing 
%             %                       coefficients gamma=[alpha, beta]
%             % t_extractedLocfreeFeaturesToConsider: 
%             %                       MxN_trainx2- tensor containing the
%             %                       locFree features of the TRAINING data; 
%             %                       2 stands for pair
%             % t_trainingLocations:  dxN_trainx2- tensor, estimated locations
%             %                       of the TRAINING data
%             % v_trainingMeasurements: N_train- vector, measured channel gains of
%             %                       TRAINING data
%             % v_funct_weights:      N-vector with weights (externally
%             %                       calculated for the query/test points)
%             % test_locF_features:   MxNx2- array, locF features of query/test points
%             % estimatedLocationVal  dxNx2- array, estimated locations of query/test points
%             
%             % Outputs:
%             % predictedMeasurements: N-vector with predicted channel gains
%             %                       for query/test points
%             
%             b_hybrid_option = true; 
%             ed_dim = size(test_locF_features); % M x N x 2
%             n_uesPair_tr = length(v_trainingMeasurements); % N_train
%             predictedMeasurements = zeros(1,ed_dim(2));
%             d = size(t_trainingLocations, 1);
%             for ind_ue_val = 1:ed_dim(2)
%                 v_dist_to_source_input_to_ker_tx = ... %(such a long name could be shortened)
%                     [test_locF_features(:,ind_ue_val,1); test_locF_features(:,ind_ue_val,2)];
%                 v_queryLocations=[estimatedLocationVal(:,ind_ue_val,1);estimatedLocationVal(:,ind_ue_val,2)];
%                 m_features_LF = zeros(ed_dim(1)*2, n_uesPair_tr); %2*M x N_train
%                 m_features_LF(1:ed_dim(1), :)     = t_extractedLocfreeFeaturesToConsider(:,:,1);
%                 m_features_LF(ed_dim(1)+1:end, :) = t_extractedLocfreeFeaturesToConsider(:,:,2);
%                 m_features_LB = zeros(d*2, n_uesPair_tr); %2*d x N_train
%                 m_features_LB(1:d, :)    = t_trainingLocations(:,:,1);
%                 m_features_LB(d+1:end,:) = t_trainingLocations(:,:,2);
%                 row_kernel_LF = feval(obj.h_kernelLF, v_dist_to_source_input_to_ker_tx,  m_features_LF);
%                 row_kernel_LB = feval(obj.h_kernelLB, v_queryLocations,                  m_features_LB);
%                 if b_hybrid_option==true
%                     predictedMeasurements(ind_ue_val)=((1-v_funct_weights').*row_kernel_LF)*m_coefs(:,1)+...
%                         (v_funct_weights'.*row_kernel_LB)*m_coefs(:,2)+mean(v_trainingMeasurements);
%                 else %if b_hybrid_option is false, this is just an additive model
%                     predictedMeasurements(ind_ue_val)=row_kernel_LF*m_coefs(:,1)+...
%                         row_kernel_LB*m_coefs(:,2)+mean(v_trainingMeasurements);
%                 end
%                 
%             end
%         end
    end %methods
    
    methods(Static)
        
        function [m_indices_all, m_indices_DAG] = build_DAG(m_locErrors)
            % encodes all partial order relations in an Mx2- matrix
            % where M is the total number of order relations (redundant)
            % m_locErrors(m_indices_out(k,1), :) >= ...
            %     m_locErrors(m_indices_out(k, 2), :) 
            % holds componentwise for all k
            [n_ues, shouldBeTwo] = size(m_locErrors);
            assert(shouldBeTwo==2);
            m_indices = zeros(n_ues*(n_ues-1), 2);
            k = 1;
            for ii = 1:n_ues
                for jj = ii+1:n_ues
                    if all(m_locErrors(ii, :) >= m_locErrors(jj, :))
                        m_indices(k, :) = [ii jj];
                    elseif all(m_locErrors(ii, :) <= m_locErrors(jj, :))
                        m_indices(k, :) = [jj ii];
                    else
                        continue % do not increase k
                    end
                    k = k+1;
                end
            end
            m_indices_all = m_indices(1:k-1,:);
            
            for ik = 1:k-1
                assert( isequal([1, 1], ...
                      m_locErrors(m_indices_all(ik, 1), :) ...
                  >=  m_locErrors(m_indices_all(ik, 2), :) ) )
            end            
            G = digraph(m_indices_all(:,1), m_indices_all(:,2));
            assert(isdag(G)); % make sure the result is a DAG
            G_reduced = transreduction(G);
            m_indices_DAG = G_reduced.Edges.EndNodes;
            for ik = 1:G_reduced.numedges
                assert( isequal([1, 1], ...
                      m_locErrors(m_indices_DAG(ik, 1), :) ...
                  >=  m_locErrors(m_indices_DAG(ik, 2), :) ) )
            end
        end
            
        
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
        
        function [v_d_out, obj_value] = fit_D_twoDimensional...
                (v_a, v_b, regDiff, m_indices, v_d_initial)
            n_ues = length(v_a);
            assert(length(v_b)==n_ues);
            assert(size(m_indices, 2) == 2) ;
            
            %solve using Quadprog
            m_I = eye(n_ues);
            m_C = zeros(size(m_indices,1), n_ues);
            for ii = 1:size(m_indices,1)
                m_C(ii,:) = m_I(m_indices(ii, 2),:) - m_I(m_indices(ii, 1),:);
            end
            m_H = 2*diag(v_b.^2);
            v_f = -2*v_a'*diag(v_b) + regDiff*ones(1,n_ues);
            m_A = [-m_I; m_I; -m_C];
            v_bq= [zeros(n_ues,1); ones(n_ues,1); zeros(size(m_indices,1), 1)];
            optimoptions = optimset('Display', 'off');
            if ~exist('m_d_initial', 'var'), v_d_initial = []; end
            [v_d_qp, fval, flag_exit] = quadprog(...
                sparse(m_H), v_f, sparse(m_A), v_bq, [], ...
                [],[],[],v_d_initial(:), optimoptions);
            obj_value = sqrt(0.5*v_d_qp'*m_H*v_d_qp + v_f*v_d_qp + v_a'*v_a);
%           obj_value == sqrt(fval + v_a'*v_a)
%           obj_value == norm( v_a - diag(v_b) * v_d_qp )
%           should hold
            v_d_out = v_d_qp;
            
%             % solve using CVX
%             cvx_begin quiet
%               variable v_d(n_ues)
%               minimize(  norm( v_a - (diag(v_b) * v_d) )  );
%               subject to
%               0 <= v_d; %#ok<VUNUS>
%               v_d <= 1; %#ok<VUNUS>
%               v_d(m_indices(:, 1)) <= v_d(m_indices(:, 2)); %#ok<VUNUS>
%             cvx_end
%             v_d_cvx = v_d;
%             obj_value = cvx_optval;
%                         
%             norm(v_d_cvx - v_d_qp)
%             keyboard

            
        end
        
        function [m_dt_dr, obj_value] ...
                = fit_D_additive(v_a, v_b, v_indices, m_d_initial)
            % Given v_a and v_b, solve for v_d
         
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



