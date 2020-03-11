classdef Simulator
    % Wireless Cartography Simulator.
    % Contains the generator, feature extractor sampler, and estimators corresponding
    % to a single experiment with the purpose of calculating the
    % estimation error of the location-based and location-free channel gain
    % cartography estimators.
    
    properties
        generator MapGenerator = PrecomputedMapGenerator();
        sampler SpectrumMapSampler
        locFreeEstimator LocationFreeEstimator
        locBasedEstimator LocationBasedEstimator
        locEstimator
        featureExtractor;
        
        simLocBased = true % Only if set to true will evaluate the error in the location based estimation.
        features  %(1,:) vector containing the indices of the features used to estimate the C2M.
        tdoaEnable=true;
        PCAenable
        PCAdimension
        
        randomized_pairs = 1; %set to 0 to use all pairs provided as input
    end
    
    methods
        
        function [NMSE_locFree,NMSE_locBased,UELocationsPairsTest, SqErr_locFree, SqErr_locBased,...
                meanErrOnEvalFeat, locFreeMapEstimate,UEFeatures,locBasedEstimate, UELocationsPairsTr,...
                avgNumberOfMissingFeat, channelForPairsVal] = ... % SqErr_sampleAverage
                simulate( obj, N_ue_pairs,combin_pairs, frac_train,frac_val_OR_val_pairs, ~,~,...
                evaluationGrid_x, evaluationGrid_y, source_loc,~, estimatedPilotSignals)
            
            allUEPairs=[evaluationGrid_x(combin_pairs(:,1))';evaluationGrid_y(combin_pairs(:,1))';...
                evaluationGrid_x(combin_pairs(:,2))';evaluationGrid_y(combin_pairs(:,2))'];

            pairSelection=randperm(size(combin_pairs,1), N_ue_pairs);
            N_pair_train=ceil(frac_train*N_ue_pairs);
            pairSelectionTrain=pairSelection(1: N_pair_train);
            UELocationsPairsTr =allUEPairs(:,pairSelectionTrain);

            if obj.randomized_pairs
                assert(isscalar(frac_val_OR_val_pairs));
                frac_val = frac_val_OR_val_pairs;
                N_pair_val=ceil(frac_val*N_ue_pairs);
                pairSelectionVal=pairSelection(N_pair_train+1: N_pair_train+N_pair_val);
                pairSelectionTest=pairSelection( N_pair_train+N_pair_val+1: end);
                UELocationsPairsVal =allUEPairs(:,pairSelectionVal);
            else
                if frac_train<1
                    warning 'Using less than 100% of the training points in not-randomized mode'
                    N_pair_train
                else
                    warning 'Validating in all points: This may be slow!'
                end
%                 pairSelectionVal = frac_val;
%                 N_pair_val=length(pairSelectionVal);
%                 pairSelectionTest = [];
                val_pairs = frac_val_OR_val_pairs;
                assert(size(val_pairs, 2)==2);
                N_pair_val = size(val_pairs,1);
                UELocationsPairsVal=[evaluationGrid_x(val_pairs(:,1))';evaluationGrid_y(val_pairs(:,1))';...
                    evaluationGrid_x(val_pairs(:,2))';evaluationGrid_y(val_pairs(:,2))'];
                pairSelectionTest = [];
            end   
            
            UELocationsPairsTest =allUEPairs(:,pairSelectionTest);

            % This line was commented because the methods from sampler that
            % we call from here do not use the .map property of the
            % Sampler.
            %obj.sampler.map = trueMap;

            obj.sampler.UELocationsPairs=UELocationsPairsTr;
            obj.sampler.Xenb=source_loc;
            obj.sampler.evalGrid_x = evaluationGrid_x;
            obj.sampler.evalGrid_y = evaluationGrid_y;
            
            % Location Free training: 
            [~,measurement_pilots]=obj.sampler.sampleGivenPilotSignals_CG(estimatedPilotSignals); %Luismi: when sampling the pilot signals, 
            % the SNR does not depend on the distance between source and node          
            measurement_pilots=reshape(measurement_pilots, [size(source_loc,2) size(measurement_pilots,2)...
                size(measurement_pilots,3)*2]);
            channelForPairsTr=zeros(1,N_pair_train);
            for ind_pair=1: N_pair_train
                channelForPairsTr(ind_pair) = obj.generator.calculateCGBetween([UELocationsPairsTr(1, ind_pair), UELocationsPairsTr(2, ind_pair)],...
                    [UELocationsPairsTr(3, ind_pair), UELocationsPairsTr(4, ind_pair)]);
            end
            extractedLocfreeFeatures=obj.featureExtractor.locFreeExtract(measurement_pilots);
            extractedLocfreeFeatures=reshape(extractedLocfreeFeatures, [size(extractedLocfreeFeatures,1) size(extractedLocfreeFeatures,2)/2,2]);
            extractedLocfreeFeaturesToConsider=extractedLocfreeFeatures(obj.features,:,:);
            
            % Location Free evaluation:  
            obj.locFreeEstimator.Xenb = source_loc; % Luismi: This property is set but the source locations are not used in the locFree estimation, only their number is used
            obj.sampler.UELocationsPairs=UELocationsPairsVal;
             [measurement_pilots_val,~]=obj.sampler.sampleGivenPilotSignals_CG(estimatedPilotSignals);
            measurement_pilots_val=reshape(measurement_pilots_val, [size(source_loc,2) size(measurement_pilots_val,2)...
                size(measurement_pilots_val,3)*2]);
            channelForPairsVal=zeros(1,N_pair_val);
            disp('Evaluating the true map at the validation set...')
            ltc = LoopTimeControl(N_pair_val);
            for ind_pair=1: N_pair_val 
                channelForPairsVal(ind_pair) = obj.generator.calculateCGBetween([UELocationsPairsVal(1, ind_pair), UELocationsPairsVal(2, ind_pair)],...
                    [UELocationsPairsVal(3, ind_pair), UELocationsPairsVal(4, ind_pair)]);
                ltc.go(ind_pair);
            end
            disp('Extracting features...')
            extractedLocfreeFeaturesVal=obj.featureExtractor.locFreeExtract(measurement_pilots_val);
            extractedLocfreeFeaturesVal=reshape(extractedLocfreeFeaturesVal, [size(extractedLocfreeFeaturesVal,1) size(extractedLocfreeFeaturesVal,2)/2,2]);
            extractedLocfreeFeaturesToConsiderVal=extractedLocfreeFeaturesVal(obj.features,:,:);
            
            if obj.PCAenable==true
                zeroMeanExtractedLocfreeFeaturesToConsider=extractedLocfreeFeaturesToConsider;
                PCAEestimatedCentersOfMass=obj.pcaComputation(zeroMeanExtractedLocfreeFeaturesToConsider);
                measurementsLF= [PCAEestimatedCentersOfMass;measurements_channel];
                UEFeatures=measurementsLF(1:end-1,:);
                [locFreeCoefficients,~,~,...
                    avgNumberOfMissingFeat,combin_sources,orthBasis,meanFeature,...
                    featCovMat,completedmeasurements] = obj.locFreeEstimator.train(measurementsLF,powerAtLocationAllsources);
                
                % Evaluation with PCA
                %                 estimatedCentersOfMassSelected=estimatedCentersOfMassSelected;
                eval_featPCA=obj.pcaComputation(estimatedCentersOfMassSelected);
                eval_newFeatNum=size(eval_featPCA,1);
                evalGrid_PCA=reshape(eval_featPCA',[evalGridSize,evalGridSize,eval_newFeatNum]);
                
                
                [meanErrOnEvalFeat,locFreeMapEstimate] = obj.locFreeEstimator.estimateGivenEstimatedDistances(...
                    locFreeCoefficients,[extractedLocfreeFeaturesToConsider;measurements_channel],completedmeasurements,...
                    estimatedCentersOfMassToConsider,power_out,combin_sources,orthBasis,mean(channelForPairsTr), meanFeature,featCovMat);
                
            else
%                 measurementsLF=[extractedLocfreeFeaturesToConsider;channelForPairs];
                disp ('Training LocFree...')
                UEFeatures=extractedLocfreeFeaturesToConsider;
                % old code:
%                     [locFreeCoefficients,~,ues_many_misses,avgNumberOfMissingFeat,combin_sources,orthBasis,meanFeature,...
%                         featCovMat,completedmeasurements]= obj.locFreeEstimator.train(extractedLocfreeFeaturesToConsider,channelForPairsTr);
%                     disp('Evaluating LocFree at validation set...')
%                     [meanErrOnEvalFeat,locFreeMapEstimate] = obj.locFreeEstimator.estimateGivenEstimatedDistances(... % sugg: change method name to: evaluateEstimatedMap[given features]
%                         locFreeCoefficients,extractedLocfreeFeaturesToConsider,channelForPairsTr,completedmeasurements,...
%                         extractedLocfreeFeaturesToConsiderVal,combin_sources,orthBasis,mean(channelForPairsTr), meanFeature,featCovMat);
                my_locFreeTrainer = LocFreeTrainer;
                my_locFreeTrainer.estimator = obj.locFreeEstimator;
                [locFreeFKM, avgNumberOfMissingFeat] = my_locFreeTrainer.train(extractedLocfreeFeaturesToConsider,channelForPairsTr);
                disp('Evaluating LocFree at validation set...')
                [meanErrOnEvalFeat,locFreeMapEstimate] = locFreeFKM.evaluate(extractedLocfreeFeaturesToConsiderVal);
                disp('Done');
            end
            
            % Compare the estimated with true channel gains 
            SqErr_locFree = obj.compute_SqErr(channelForPairsVal, locFreeMapEstimate);
            NMSE_locFree = SqErr_locFree/sum((channelForPairsVal-mean(channelForPairsVal)).^2);
%           SqErr_sampleAverage = obj.compute_SqErr(trueMap, ones(size(trueMap))*mean(channelForPairsTr));
            
            
            if obj.simLocBased==true
                
                % Location Based: train estimator
                extractedLocBasedFeatures=obj.featureExtractor.locBasedExtract(measurement_pilots);
%                extractedLocBasedFeatures=reshape(extractedLocBasedFeatures, [size(extractedLocBasedFeatures,1) size(extractedLocBasedFeatures,2)/2,2]);
%                measurementsLB=[extractedLocBasedFeatures;channelForPairs]; 
                obj.locBasedEstimator.Xenb = source_loc;
                obj.locEstimator.Xenb = source_loc;
                if obj.tdoaEnable==false
                    estimatedLocation=obj.locBasedEstimator.estimateLocationFromDistances (extractedLocBasedFeatures);
                    estimatedLocation=reshape(estimatedLocation, [size(estimatedLocation,1) size(estimatedLocation,2)/2,2]);
                else
                    estimatedLocation=obj.locEstimator.estimateLocationIRWSRDLS(extractedLocBasedFeatures);
                    estimatedLocation=reshape(estimatedLocation, [size(estimatedLocation,1) size(estimatedLocation,2)/2,2]);
                    
                end, disp ('Training LocBased...')
              % old code:
               % [locBasedCoefficients, locBasedIntercept] = obj.locBasedEstimator.train(estimatedLocation,channelForPairsTr);
                my_locBasedTrainer = LocBasedTrainer;
                my_locBasedTrainer.estimator = obj.locBasedEstimator;
                locBasedFKM = my_locBasedTrainer.train(estimatedLocation,channelForPairsTr);
                
                % Location Based evaluation
                disp('Evaluating LocBased at validation set...')
                extractedLocBasedFeaturesVal=obj.featureExtractor.locBasedExtract(measurement_pilots_val);
                estimatedLocationVal=obj.locEstimator.estimateLocationIRWSRDLS(extractedLocBasedFeaturesVal);
                estimatedLocationVal=reshape(estimatedLocationVal, [size(estimatedLocationVal,1) size(estimatedLocationVal,2)/2,2]);
                
                % old code:
                   % locBasedEstimate = obj.locBasedEstimator.estimate(locBasedCoefficients, estimatedLocation,...
                   %                estimatedLocationVal, locBasedIntercept); %evaluateEstimateMap
                locBasedEstimate = locBasedFKM.evaluate(estimatedLocationVal);
                
                % Compare the estimated with true channel gains 
                SqErr_locBased= obj.compute_SqErr(channelForPairsVal, locBasedEstimate);
                NMSE_locBased = SqErr_locBased/sum((channelForPairsVal-mean(channelForPairsVal)).^2);
            else
                locBasedEstimate = 0;
                SqErr_locBased = 0;
            end
            
        end
        
        function [allPointEstimatedLocation,allLocFreeFeatures,allLocBasedFeatures,...
                  estimatedPilotSignals, channel_out, source_loc, evaluationGrid_x, evaluationGrid_y...
                  ] = precomputeEstimatedLocations(obj)
            % Luismi: please write a (short) description of this method
            
            [estimatedPilotSignals, channel_out, source_loc, evaluationGrid_x, evaluationGrid_y]=obj.generator.generateCGMap(obj.generator.xt);
            
            allPointSampler=obj.sampler;
            allPointSampler.map = channel_out;
            allPointSampler.Xenb=source_loc;
            allPointSampler.evalGrid_x = evaluationGrid_x;
            allPointSampler.evalGrid_y = evaluationGrid_y;
            allPointSampler.UELocationsPairs=[obj.generator.xt'.*ones(length(obj.generator.xt),numel(evaluationGrid_x(:)'));
                                             flipud(evaluationGrid_x(:)'); 
                                             flipud(evaluationGrid_y(:)');]; %every point location
            [allnoiseFreePilots,~]=allPointSampler.sampleGivenPilotSignals_CG(estimatedPilotSignals); 
            allnoiseFreePilots=reshape(allnoiseFreePilots, [size(source_loc,2) size(allnoiseFreePilots,2)...
                                size(allnoiseFreePilots,3)*2]);
            if obj.simLocBased==true
                allLocBasedFeatures= obj.featureExtractor.locBasedExtract(allnoiseFreePilots);
                % Add noise to the locB features
                allLocBasedFeatures = allLocBasedFeatures + 0*randn(size(allLocBasedFeatures));
                allLocFreeFeatures=obj.featureExtractor.locFreeExtract(allnoiseFreePilots);
                allLocFreeFeatures=reshape(allLocFreeFeatures, [size(allLocFreeFeatures,1) size(allLocFreeFeatures,2)/2,2]);
                
                obj.locBasedEstimator.Xenb = source_loc;
                obj.locEstimator.Xenb = source_loc;
                if obj.tdoaEnable==false
                    allPointEstimatedLocation=obj.locBasedEstimator.estimateLocationFromDistances(allLocBasedFeatures);
                    allPointEstimatedLocation=reshape(allPointEstimatedLocation, [size(allPointEstimatedLocation,1) size(allPointEstimatedLocation,2)/2,2]);
                else
                    allPointEstimatedLocation=obj.locEstimator.estimateLocationIRWSRDLS(allLocBasedFeatures);
                    allPointEstimatedLocation=reshape(allPointEstimatedLocation, [size(allPointEstimatedLocation,1) size(allPointEstimatedLocation,2)/2,2]);
                    %              allPointEstimatedLocation=allPointSampler.UELocations;
                end
            else
                
                allLocFreeFeatures=obj.featureExtractor.locFreeExtract(allnoiseFreePilots);
                allLocFreeFeatures=reshape(allLocFreeFeatures, [size(allLocFreeFeatures,1) size(allLocFreeFeatures,2)/2,2]);
                allLocBasedFeatures=0;
                allPointEstimatedLocation=0;
            end
        end
        
        
        function [U1,newFeat]=pcaComputationFirst(obj,centerOfMass)
            Option=1;
            if Option==1
                meanofFeat=mean(centerOfMass,2);
                centeredFeat=centerOfMass-meanofFeat;
                [U,~,~] = svd(centeredFeat);
                U1=U(:,1:obj.PCAdimension);
                newFeat=U1'*centeredFeat;
            else
                nFeatures=size(centerOfMass,1); %featSize=size(centerOfMass,2);
                %             meanCenterOfMass=mean(centerOfMass,3);
                CovFeat=zeros(nFeatures);
                for ind_feat1=1:nFeatures
                    for ind_feat2=1:nFeatures
                        tempCov=cov(centerOfMass(ind_feat1,:),....
                            centerOfMass(ind_feat2,:));
                        CovFeat(ind_feat1,ind_feat2)=tempCov(1,2);
                    end
                end
                [Eigvector,~]=eig(CovFeat);
                flipped_Eigvector=fliplr(Eigvector);
                pcaSpace=flipped_Eigvector(:,1:obj.PCAdimension); % space corresponding to the... largest eigen value
                newFeat=pcaSpace'*(centerOfMass-mean(centerOfMass,2));
            end
        end
        
        function newFeat=pcaComputation(obj,centerOfMass)
            Option=1;
            if Option==1
                meanofFeat=mean(centerOfMass,2);
                centeredFeat=centerOfMass-meanofFeat;
                [U,~,~] = svd(centeredFeat);
                U1=U(:,1:obj.PCAdimension);
                newFeat=U1'*centeredFeat;
            else
                nFeatures=size(centerOfMass,1); %featSize=size(centerOfMass,2);
                %             meanCenterOfMass=mean(centerOfMass,3);
                CovFeat=zeros(nFeatures);
                for ind_feat1=1:nFeatures
                    for ind_feat2=1:nFeatures
                        tempCov=cov(centerOfMass(ind_feat1,:),....
                            centerOfMass(ind_feat2,:));
                        CovFeat(ind_feat1,ind_feat2)=tempCov(1,2);
                    end
                end
                [Eigvector,~]=eig(CovFeat);
                flipped_Eigvector=fliplr(Eigvector);
                pcaSpace=flipped_Eigvector(:,1:obj.PCAdimension); % space corresponding to the... largest eigen value
                newFeat=pcaSpace'*(centerOfMass-mean(centerOfMass,2));
            end
        end
        
    end %methods
    
    methods (Static)
        function SqErr = compute_SqErr(X,Y)
            SqErr=sum((X-Y).^2, 'omitnan');
        end  
    end
    
end %classdef
