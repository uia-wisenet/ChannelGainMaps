classdef LocationFreeEstimator < Estimator
    %LocationFreeEstimator Location-Free cartography class.
    %   This class allows to train and query a kernel
    %   machine desgined to estimate the channel gain from location 
    %   features received at sensors pairs.
    
    properties
        kernel % function handle defining the kernel
        regularizationParameter % lambda (ridge regression)
        Xenb % locations of the sources % Luismi: this property is not used: Only its size is accessed in the train method.
        enableMissingData
        completionMethodSelect % 1 for modified EGM, and 2 for Manopt
        receiverSensitivity % missing features, can be changes   ... use it for threshold
        gamma;% parameter used in the kernel of missing features: low rank and missing data
        desiredRank
        regParEval %regularization parameter for completing missing features when evaluating the map
        evalOption % choose the option of determining missing features in the map evaluation
        
        b_verbose = 1
    end
    
    methods
        function [coefficientsOut, intercept,ues_many_misses,...
                avgNumberOfMissingFeat,combin_sources,orthBasis,...
                meanFeature,featCovarianceMat,completedmeasurements] ...
                = train(obj, LocFeatures,...
                channelForPairs)
            % given a 3D tensor including the locF features (M-by-N-by-2 array
            % where M is the number of features at one sensor location,
            % N is the number of pairs for training, and 2 stands for pair ) and
            % channel gains, optimize the coefficients
            % of the kernel machine that estimates the channel gain.
            n_sources = size( LocFeatures,1); % n_sources  symbolizes number of features from all sources
            n_ues=size(LocFeatures, 2);
            %             featRank=4;% is the rank of the feature space
   
            if  obj.enableMissingData==true  % incorporates the missing feature case.
                real_n_sources=size(obj.Xenb,2);
                combin_sources=combnk(1:real_n_sources,2);
                if real_n_sources < 6
                    combin_sources=flipud(combin_sources);
                end
                identify_misses=(channelForPairs<obj.receiverSensitivity);
                X=LocFeatures(1:end-1,:);
                %                 if obj.receiverSensitivity==0
                %                     Q=orth(X);
                %                     completedmeasurements =Q(:,1:featRank)'*X;
                %                 else
                
                o=zeros(n_sources,n_ues);
                ues_many_misses=[];
                allNumberOfMissingFeat=zeros(1,n_ues);
                for ind_n=1:n_ues
                    miss_ind=find(identify_misses(:,ind_n)==1);
                    for ind_missing=1:length(miss_ind)
                        miss_indices_on_features=union(find(combin_sources(:,1)==miss_ind(ind_missing)),...
                            find(combin_sources(:,2)==miss_ind(ind_missing)));
                        X(miss_indices_on_features,ind_n)=NaN; %NaN symbolizes a missing features
                    end
                    o(:,ind_n)=double(~isnan(X(:,ind_n)));
                    numberOfMissingFeat=sum(isnan(X(:,ind_n)));
                    allNumberOfMissingFeat(ind_n)=numberOfMissingFeat;
                    if numberOfMissingFeat> size(combin_sources,1)-obj.desiredRank % Removes the measurements having a large number of missing features
                        % there must be at least 4 features available.
                        ues_many_misses=[ues_many_misses, ind_n];
                    end
                end
                avgNumberOfMissingFeat=mean(allNumberOfMissingFeat);
                X(:,ues_many_misses)=[];
                o(:,ues_many_misses)=[];
                X_nan=X;
                
                %                  for ind_n=1:n_ues
                %                     par=randperm(n_sources,obj.receiverSensitivity);
                %                     for ind_missing=1:length(par)
                %                         X(par(ind_missing),ind_n)=NaN; %NaN symbolizes a missing attribute in the the n-th observed sample
                %                     end
                %                     o(:,ind_n)=double(~isnan(X(:,ind_n)));
                %                  end
                
                X(isnan(X))=0;
                
                if obj.completionMethodSelect==1
                    completedmeasurementsToProject =singularValProj(X,o,obj.desiredRank);
                elseif obj.completionMethodSelect==2
                    completedmeasurementsToProject=low_rank_matrix_completion(X,o,obj.desiredRank);
                else
                    opts = [];
                    [X_lma,Y_lma,~] = lmafit_mc_adp(X_nan,obj.desiredRank,opts);
                    completedmeasurementsToProject=X_lma*Y_lma;
                end
                
                %Project completed measurements
                orthBasis=orth(completedmeasurementsToProject);
                completedmeasurements =orthBasis'*(completedmeasurementsToProject); %
                featRank=size(orthBasis,2);
                if obj.evalOption==1
                    meanFeature=mean(completedmeasurementsToProject,2);
                    featCovarianceMat=cov(completedmeasurementsToProject');
                else
                    Coeff_B=orthBasis'*completedmeasurementsToProject;
                    meanFeature=mean(Coeff_B,2);
                    featCovarianceMat=cov(Coeff_B');
                end
                
                %                 end
                n_ues=size(completedmeasurements,2);
                LocFeatures(:,ues_many_misses)=[];
            else
                ues_many_misses=[];
                avgNumberOfMissingFeat=0;
                combin_sources=0;
                orthBasis=0;
                meanFeature=0;
                featCovarianceMat=0;
                completedmeasurements=0;
            end
            
            Ke1=zeros(n_ues,n_ues); % Kernel matrix
            
            % OLD CODE:
%             for i = 1:n_ues
%                 for k=1:n_ues                    
%                     if  (obj.enableMissingData==true) && (obj.receiverSensitivity~=0)
%                         x_input_to_kern=zeros(featRank,1);
%                         y_input_to_kern=zeros(featRank,1);
%                         for ind_s_1=1:featRank
%                             x_input_to_kern(ind_s_1)=completedmeasurements(ind_s_1,i);
%                             y_input_to_kern(ind_s_1)=completedmeasurements(ind_s_1,k);
%                         end
%                         Ke1(i,k) = feval(obj.kernel, x_input_to_kern, y_input_to_kern);
%                         
%                     else
%                         x_input_to_kern_tx=zeros(n_sources,1);
%                         y_input_to_kern_tx=zeros(n_sources,1);
%                         x_input_to_kern_rx=zeros(n_sources,1);
%                         y_input_to_kern_rx=zeros(n_sources,1);
%                         for ind_s_1=1:n_sources
%                             x_input_to_kern_tx(ind_s_1)=LocFeatures(ind_s_1,i,1);
%                             y_input_to_kern_tx(ind_s_1)=LocFeatures(ind_s_1,k,1);
%                             x_input_to_kern_rx(ind_s_1)=LocFeatures(ind_s_1,i,2);
%                             y_input_to_kern_rx(ind_s_1)=LocFeatures(ind_s_1,k,2);
%                         end
%                         Ke1(i,k) = feval(obj.kernel,  [x_input_to_kern_tx;x_input_to_kern_rx],[y_input_to_kern_tx;y_input_to_kern_rx]);
%                     end
%                 end
%             end           
            if  (obj.enableMissingData==true) && (obj.receiverSensitivity~=0)
                for i = 1:n_ues
                    for k=1:n_ues
                        x_input_to_kern=zeros(featRank,1);
                        y_input_to_kern=zeros(featRank,1);
                        for ind_s_1=1:featRank
                            x_input_to_kern(ind_s_1)=completedmeasurements(ind_s_1,i);
                            y_input_to_kern(ind_s_1)=completedmeasurements(ind_s_1,k);
                        end
                        Ke1(i,k) = feval(obj.kernel, x_input_to_kern, y_input_to_kern);
                    end
                end
            else
%                 try
%                     t3_all_inputs_to_kernels = zeros(n_sources*2, n_ues, n_ues);
%                     t3_all_inputs_to_kernels(1:n_sources, :,:)     = ...
%                         repmat(LocFeatures(:,:,1), [1 1 n_ues]);
%                     t3_all_inputs_to_kernels(n_sources+1:end, :,:) = ...
%                         repmat(LocFeatures(:,:,2), [1 1 n_ues]);
%                     Ke1 = squeeze(feval(obj.kernel, t3_all_inputs_to_kernels, ...
%                         permute(t3_all_inputs_to_kernels, [1 3 2])));
%                 catch ME
%                     if isequal(ME, 'MATLAB:array:SizeLimitExceeded')
                        disp('Building array row by row')
                        Ke1 = zeros(n_ues, n_ues);
                        ltc = LoopTimeControl(n_ues);
                        for i = 1:n_ues
                            m_row_inputs_to_kernels = repmat([LocFeatures(:, i, 1);...
                                LocFeatures(:,i,2)], [1 n_ues]);
                            my_row = feval(obj.kernel, m_row_inputs_to_kernels, ...
                                [LocFeatures(:,:,1);LocFeatures(:,:,2)]);
                            Ke1(i,:) = my_row;
                            if obj.b_verbose, ltc.go(i); end
                        end
%                     else
%                         rethrow(ME);
%                     end
%                 end
            end
%             norm(Ke1-Ke3)
%             keyboard
            intercept = mean(channelForPairs);
            channelForPairs=channelForPairs-intercept;
            coefficientsOut=(Ke1+((n_ues*obj.regularizationParameter)*eye(n_ues)))\channelForPairs(:);
            % NOTE: the previous 3 lines were inside the for loop before. 
            % This was unnecessary and I guess it was slowing down the code.
        end
               
        function [predictedMeasurements, meanErrOnEvalFeat] = estimateGivenEstimatedDistances(obj, coefficients, extractedLocfreeFeaturesToConsider,trainingMeasurements,completedmeasurements,...
                test_locF_features,combin_sources,orthBasis,avgPower, meanFeat,CovarMat)
            % Estimate channel gain using Location Free cartography.
            % Inputs: test_locF_features M-by-N-by-2 array
            % where M is the number of features at one sensor location
            % N is the number of pairs where the kernel machine is
            % evaluated,  and 2 stands for pair
            ed_dim = size(test_locF_features);
            featNum=size(combin_sources,1);
            n_uesPair_tr=length(trainingMeasurements);
            predictedMeasurements=zeros(1,ed_dim(2));
            errAllCompletedEvalFeat=zeros(ed_dim(2),1);
            for ind_ue_val = 1:ed_dim(2)
                if  obj.enableMissingData==true
                    n_uesPair_tr=size(completedmeasurements,2);
                    identify_missesEval=(allPointAllSourcePower(:,kx, ky)<obj.receiverSensitivity);
                    evalFeat=permute(test_locF_features(kx, ky, :),[3 2 1]);
                    miss_ind=find(identify_missesEval==1);
                    for ind_missing=1:length(miss_ind)
                        miss_indices_on_features=union(find(combin_sources(:,1)==miss_ind(ind_missing)),...
                            find(combin_sources(:,2)==miss_ind(ind_missing)));
                        evalFeat(miss_indices_on_features)=NaN; %NaN symbolizes a missing features
                    end
                    numberOfMissingFeatAtEval=sum(isnan(evalFeat));
                    if ( numberOfMissingFeatAtEval >featNum -obj.desiredRank) || (size(orthBasis,2) < obj.desiredRank)
                        predictedMeasurements(kx,ky)=avgPower;
                        continue
                    end
                    
                    indOfObserved_Feat=find(~isnan(evalFeat));
                    obervedFeatNum=length(indOfObserved_Feat);
                    rowSelectMat=zeros(obervedFeatNum,featNum);
                    evalFeat(isnan(evalFeat))=0;
                    for ind_obsFeat=1:obervedFeatNum
                        rowSelectMat(ind_obsFeat,indOfObserved_Feat(ind_obsFeat))=1;
                    end
                    if obj.evalOption==1
                        weights=(orthBasis'*rowSelectMat'*rowSelectMat*orthBasis+obj.regParEval*orthBasis'*orthBasis)\...%(CovarMat\
                            (orthBasis'*rowSelectMat'*rowSelectMat*evalFeat+obj.regParEval*orthBasis'*meanFeat);
                    else
                        weights=(orthBasis'*rowSelectMat'*rowSelectMat*orthBasis+obj.regParEval*inv(CovarMat))\...
                            (orthBasis'*rowSelectMat'*rowSelectMat*evalFeat+obj.regParEval*(CovarMat\meanFeat));
                    end
                    dist_to_source_input_to_ker=orthBasis*weights;
                    %                         dist_to_source_input_to_ker(:)=test_locF_features(kx, ky, :);
                    errAllCompletedEvalFeat(kx, ky)=norm(weights-orthBasis'*permute(test_locF_features_tx(kx, ky, :),[3 2 1]));
                else
                    errAllCompletedEvalFeat(ind_ue_val)=0;
                    
                    dist_to_source_input_to_ker_tx=[test_locF_features(:,ind_ue_val,1); test_locF_features(:,ind_ue_val,2)];
                end
                
                %OLD CODE:
%                 row_kernel=zeros(1,n_uesPair_tr);
%                 for ind_nues=1:n_uesPair_tr
%                     if  obj.enableMissingData==true
%                         row_kernel(ind_nues)=feval(obj.kernel, weights, completedmeasurements(:,ind_nues)); 
%                         % here the weights=orthBasis'*dist_to_source_input_to_ker
%                         %  =orthBasis'*orthBasis*weights=weights;
%                     else
%                         row_kernel(ind_nues)=feval(obj.kernel, dist_to_source_input_to_ker_tx,...
%                             [extractedLocfreeFeaturesToConsider(:,ind_nues,1); extractedLocfreeFeaturesToConsider(:,ind_nues,2)]);
%                         % CRITICAL LINE (executed 8000000 times)
%                     end
%                 end
                if  obj.enableMissingData
                    row_kernel=zeros(1,n_uesPair_tr);
                    for ind_nues=1:n_uesPair_tr
                        row_kernel(ind_nues)=feval(obj.kernel, weights, completedmeasurements(:,ind_nues));
                        % here the weights=orthBasis'*dist_to_source_input_to_ker
                        %  =orthBasis'*orthBasis*weights=weights;
                    end
                else
                    m_features_now = zeros(ed_dim(1)*2, n_uesPair_tr);
                    m_features_now(1:ed_dim(1), :)     = extractedLocfreeFeaturesToConsider(:,:,1);
                    m_features_now(ed_dim(1)+1:end, :) = extractedLocfreeFeaturesToConsider(:,:,2);
                    row_kernel = feval(obj.kernel, dist_to_source_input_to_ker_tx, ...
                        m_features_now);
                end
                %norm(row_kernel-row_kernel2) % should be 0
                
                predictedMeasurements(ind_ue_val)=row_kernel*coefficients+mean(trainingMeasurements);
                
            end
            meanErrOnEvalFeat=mean(mean(errAllCompletedEvalFeat));
        end
        
    end %methods
end


