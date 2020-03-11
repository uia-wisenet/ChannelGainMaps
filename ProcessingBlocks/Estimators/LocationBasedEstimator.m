classdef LocationBasedEstimator < Estimator
    %LocationBasedEstimator Location-Based cartography class.
    %   This class allows to train and query a kernel
    %   machine desgined to estimate the channel gain from
    %   the locations of UE pairs.
    
    properties
        kernel % function handle defining the kernel
        regularizationParameter % lambda (ridge regression)
        Xenb % locations of the eNBs
        locationNoiseSTD % standard deviation of the noise
    end
    
    methods
        
        function [coefficientsOut, intercept] = train(obj,estimatedLocation,channelMeasurements)
            % given a 3D tensor including the location estimates (2-by-N-by-2 array
            % where the first 2 is the area dimensionality, 
            % N is the number of pairs used to train, and the last 2 stands for pair) and
            % channel gains, optimize the coefficients
            % of the kernel machine that estimates the channel gain.
            
            n_ues= length(channelMeasurements);
            d = size(estimatedLocation, 1);
            assert(d==2 || d==3);
            assert(size(estimatedLocation, 2)==n_ues);
%             % OLD CODE
%             Ke1=zeros(n_ues,n_ues); % Kernel matrix
%             for i = 1:n_ues
%                 for k=1:n_ues
%                     x_input_to_kern_tx=zeros(2,1); %TODO: if going to higher dimensionality %(e.g. 3D), this should be variable d
%                     y_input_to_kern_tx=zeros(2,1);
%                     x_input_to_kern_rx=zeros(2,1); %TODO: if going to higher dimensionality %(e.g. 3D), this should be variable d
%                     y_input_to_kern_rx=zeros(2,1);
%                     for ind_s_1=1:2 %TODO: if going to higher dimensionality
%                         %(e.g. 3D), this should be variable d
%                         x_input_to_kern_tx(ind_s_1)=estimatedLocation(ind_s_1,i,1);
%                         y_input_to_kern_tx(ind_s_1)=estimatedLocation(ind_s_1,k,1);
%                         x_input_to_kern_rx(ind_s_1)=estimatedLocation(ind_s_1,i,2);
%                         y_input_to_kern_rx(ind_s_1)=estimatedLocation(ind_s_1,k,2);
%                     end
%                     
%                     Ke1(i,k) = feval(obj.kernel, [x_input_to_kern_tx;x_input_to_kern_rx],[y_input_to_kern_tx;y_input_to_kern_rx]);
%                 end
%             end
%             norm(Ke1-Ke2)/norm(Ke1)
% NEWER CODE, BUT ALREADY OLD
%             try
%                 t3_all_inputs_to_kernels = zeros(d*2, n_ues, n_ues);
%                 t3_all_inputs_to_kernels(1:d, :,:)     = repmat(estimatedLocation(:,:,1), [1 1 n_ues]);
%                 t3_all_inputs_to_kernels(d+1:end, :,:) = repmat(estimatedLocation(:,:,2), [1 1 n_ues]);
%                 Ke2 = squeeze(feval(obj.kernel, t3_all_inputs_to_kernels, permute(t3_all_inputs_to_kernels, [1 3 2])));
%             catch ME
%                 if isequal(ME, 'MATLAB:array:SizeLimitExceeded')
                    disp('Building array row by row')
                    Ke1 = zeros(n_ues, n_ues);
                    ltc = LoopTimeControl(n_ues);
                    for i = 1:n_ues
                        m_row_inputs_to_kernels = repmat([estimatedLocation(:, i, 1);...
                            estimatedLocation(:,i,2)], [1 n_ues]);
                        my_row = feval(obj.kernel, m_row_inputs_to_kernels, ...
                            [estimatedLocation(:,:,1);estimatedLocation(:,:,2)]);
                        Ke1(i,:) = my_row;
                        ltc.go(i);
                     end
%                 else
%                     rethrow(ME);
%                 end
%            end
% norm(Ke1-Ke2)/norm(Ke2)
% keyboard
            intercept = mean(channelMeasurements);
            channelMeasurements=channelMeasurements- intercept;
            
            coefficientsOut=(Ke1+((n_ues*obj.regularizationParameter)*eye(n_ues)))\channelMeasurements';
        end
        
        function predictedMeasurements = estimate(obj, coefficients, trainingLocations, estimatedLocationVal,intercept)
            % Estimate channel gain using Location-Based cartography.
            % This function receives as inputs the kernel weights, locations of the training and query
            % UE pairs (of size 2-by-N-by-2 array as in the training).
            n_ues= size(trainingLocations,2);
            n_ues_val= size(estimatedLocationVal,2);
            p2=zeros(1,n_ues_val);
            d = size(trainingLocations, 1);
            for ind_ue_val = 1:n_ues_val
                queryLocations=[estimatedLocationVal(:,ind_ue_val,1);estimatedLocationVal(:,ind_ue_val,2)];
                % OLD CODE:
                % row1_kernel=zeros(1,n_ues);
                % for ind_nues=1:n_ues
                %     row1_kernel(ind_nues)=feval(obj.kernel, queryLocations, [trainingLocations(:,ind_nues,1); trainingLocations(:,ind_nues,2)]);
                %     % CRITICAL LINE (executed 8000000 times)
                % end
                
                m_features_now = zeros(d*2, n_ues);
                m_features_now(1:d, :)    = trainingLocations(:,:,1);
                m_features_now(d+1:end,:) = trainingLocations(:,:,2);
                row1_kernel = feval(obj.kernel, queryLocations, m_features_now);
                %norm(row2_kernel-row1_kernel)
                %keyboard
                p2(ind_ue_val)=row1_kernel*coefficients+intercept;
            end
            
            predictedMeasurements=p2;
        end
        
        function estimatedLocation = estimateLocationFromDistances (obj,measurements)
            % Square range-based LS. Taken from the paper  beck2008exact
            
            n_ues=size(measurements,2);
            n_sources = size( measurements,1);
            A=[-2*obj.Xenb(1,1) -2*obj.Xenb(2,1) 1; -2*obj.Xenb(1,2) -2*obj.Xenb(2,2) 1; -2*obj.Xenb(1,3) -2*obj.Xenb(2,3) 1];
            D=[1 0 0;0 1 0;0 0 0]; f=[0; 0; -0.5];
            b=zeros(n_sources, n_ues);
            for ind2_s=1:n_sources
                for i=1:n_ues
                    b(ind2_s,i)=(measurements(ind2_s,i)).^2-(obj.Xenb(1,ind2_s)).^2-(obj.Xenb(2,ind2_s)).^2;
                end
            end
            
            ue_loc_calc=zeros(2,n_ues);
            for i=1:n_ues
                e= -1/max(eig(D,A'*A)); %  bisection interval lower limit
                a2=1e2; % bisection interval upper limit
                tol=1e-4; % Tolerance (final error)
                lam_opt = obj.bisection(e,a2,A,b(:,i),D,f,tol);
                if isnan(lam_opt)
                    y_est=zeros(2,1);
                else
                    y_est=(A'*A + lam_opt*D)\(A'*b(:,i)-lam_opt*f);
                end
                ue_loc_calc(:,i)= y_est(1:2);
            end
            estimatedLocation=ue_loc_calc;
        end
        
        function estimatedLocation = estimateLocationFromDistDifferences (obj,measurements)
            %TDoA Semidefinite programming. Taken from the paper  du2014semidefinite,luo2010semidefinite,
            dim=2; % consider sources and UEs in 2D
            n_ues=size(measurements,2);
            n_sources = size( measurements,1)-1;
            n_sources_tdoa=n_sources; % The first source acts as a reference source
            
            V=eye(n_sources_tdoa)+ones(n_sources_tdoa); %
            Qr=obj.locationNoiseSTD.^2*V;
            W=Qr; % TO DO: Change the way W is computed
            %            B=2*diag(dist2to4); W=B*Qr*B';
            
            A=zeros(n_sources_tdoa,dim+1,n_ues);
            b=zeros(n_sources_tdoa, n_ues);
            P=zeros(dim+2,dim+2,n_ues);
            ue_loc_calc=zeros(dim,n_ues);
            
            for i=1:n_ues
                for ind2_s=1:n_sources
                    A(ind2_s,:,i)=2*[(obj.Xenb(:,ind2_s+1)-obj.Xenb(:,1))',-measurements(ind2_s,i)];
                    b(ind2_s,i)=measurements(ind2_s,i).^2 - (norm(obj.Xenb(:,ind2_s+1)).^2) ...
                        +(norm(obj.Xenb(:,1)).^2);
                end
                P(:,:,i)=[A(:,:,i)'*(W\A(:,:,i)),  -A(:,:,i)'*(W\b(:,i)); ...
                    -(b(:,i)'/W)*A(:,:,i),   (b(:,i)'/W)*b(:,i)];
                
                
                %                cvx_solver sedumi
                cvx_solver SDPT3
                cvx_begin quiet
                variable z(dim+1)
                variable Z(dim+1,dim+1) symmetric
                minimize (trace( [Z, z; z', 1]*P(:,:,i)));
                subject to
                Z(dim+1,dim+1)==trace(Z(1:dim,1:dim))- 2*obj.Xenb(:,1)'*z(1:dim)+norm(obj.Xenb(:,1)).^2; %#ok<EQEFF>
                %                Z(end,end)==1;
                %                ZZ(dim+2,1:dim+1)==ZZ(1:dim+1,dim+2)';
                %                ZZ(1:dim+1,1:dim+1)==ZZ(1:dim+1,dim+2)*ZZ(dim+2,1:dim+1)';
                %                ZZ==semidefinite(dim+2);
                %                z(dim+1)>0;
                [Z, z; z', 1]==semidefinite(dim+2); %#ok<EQEFF>
                cvx_end
                
                ue_loc_calc(:,i)= z(1:dim);
                %                  dist2to4=zeros(n_sources_tdoa,1);
                %                  for ind2_rec=1:n_sources
                %                 dist2to4(ind2_rec)=norm(ue_loc_calc(:,i)-obj.Xenb(:,ind2_rec+1));
                %                  end
            end
            estimatedLocation=ue_loc_calc;
        end
        
    end
    
    methods (Static)
        function p = bisection(a1,a2,A,b,D,f,tol)
            y=(A'*A + a1*D)\(A'*b-a1*f);
            f_a1=y'*D*y+2*f'*y;
            y=(A'*A + a2*D)\(A'*b-a2*f);
            f_a2=y'*D*y+2*f'*y;
            if f_a2<-1e12
                warning('Interval too large. Returning NaN')
                p = nan;
                return
            elseif f_a2>0
                %warning('Choosing a larger interval')
                a2=a2*10;
                p = LocationBasedEstimator.bisection(a1,a2,A,b,D,f,tol);
            else
                p = (a1 + a2)/2;
                y=(A'*A + p*D)\(A'*b-p*f);
                f_p=y'*D*y+2*f'*y;
                
                if abs(f_p) > tol
                    if f_a1*f_p<0
                        a2 = p;
                    else
                        a1 = p;
                    end
                    p=LocationBasedEstimator.bisection(a1,a2,A,b,D,f,tol);
                end
            end
        end
    end
end

