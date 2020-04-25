classdef CartographySimulations
    %method CartographySimulations.run is
    %executed once the simulation parameters (properties) are defined.
    
    properties % default values = baseline testcase
        rng_seed = 6; %seed of the random number generator
        forcePrecomputeDistances = false; % if set to true, it will ignore the precomputed distances file
        inParallel = true; % compute in parallel (use PARFOR)
        
        % properties pertaining the statistical analysis
        
        gridSize 	% Size of the grid
        N_ue_pairs        % Vector containing different values for the number of UEs/sensors
        n_runs           % Number of Monte Carlo runs to be evaluated in each experiment
        
        % Number of features. Set N_feat = 3 (default) for comparison of
        % locFree vs locBased;
        % set N_feat = 1:3 for evaluating impact of M
        N_feat_used % total number of features per sensor location
        N_feat % used for PCA dimension or Rx sensitivity in missing feature
        numberOfSources
        
        % Should simulate the location based estimation.
        % set simLocBased = true (default) for comparison of
        % locFree vs locBased;
        % set simLocBased = false for evaluating impact of features M, PCA, and missing features for LocF scheme.
        simLocBased
        % Hyperparameters of the kernel machines
        computeOnlyFeatures
        
        kernelSigmaLF  % Sigma width for the gaussian kernel in LocFree
        kernelSigmaLB% Sigma width for the gaussian kernel in LocBased
        lambdaLF      % Ridge regularization parameter (locFree )
        lambdaLB      % Ridge regularization parameter (locBased)
        
        % Properties of the test scenario
        
        f  %Carrier frequency
        locationNoiseSTD  % Standard deviation of location feature noise
        powerNoiseSTD     % Standard deviation of measured power noise
        pilotNoiseSTD
        maxSamplesPerPilot
        % Bandwidth of the receiver (will determine the resolution
        % of the digital impulse response)
        receiverBandwidth
        % each element of the cell array is a vector indicating which features
        % (i.e. signals from which sources) are considered for the
        featureCell % = {[2], [1 2], [1 2 3],[1 2 3 4],[1 2 4 5 6],[1:6],[1:7],[1 2 3 4 8 9 10],[4:10],[1:10]}; %#ok<NBRAK>
        
        
        testDiffSetOfFeatures;
        PCAenabler % if true, the PCA is taken into account
        considerMissingData % accommodating the misses:  see classfication with low rank and missing data
        completionMethodSelect=2 %selects either ManOpt or other specified method (Poject Gradient Descent, EGM, etc. )
        desiredRank=4 % desired rank of the matrix to be completed
        regularParEval=5.42e-0 % %regularization parameter for completing missing features when evaluating the map
        %         featMissing=4; % number of missing features, which can be changed
        evalcoeffOption=2 % choose the option of determining missing features in the map evaluation
        % 1 for computing mean  and covariance between completed features;
        % else  or computing mean  and covariance between coefficients correponding
        %to completed features;
        parGamma %parameter used in the kernel of missing features: low rank and missing data
        randomized_pairs = 1;
        
        b_new_call = 0;
        
    end
    
    
    
    methods (Static)
        
        function featureCell = featureCombinations(featNum,startingFeat)
            featureCell={};
            for indFeat=startingFeat:featNum
                selectedComb=combnk(1:featNum,indFeat);
                for indc=1:size(selectedComb,1)
                    rowSelec=selectedComb(indc,:);
                    featureCell{end+1}=rowSelec;
                end
            end
        end
        
    end %methods (Static)
    
    methods
        function [generator, source_loc] = baselineGenerator(obj, selectedWalls, x_wall_file_name, y_wall_file_name)
            warning('Using CartographySimulations.baselineGenerator, which is deprecated.')
            warning('This call can be substituted with baselineGenerator2')
            if obj.b_new_call
                generator_template =  MultiWallChannelGainGenerator;
                generator_template.f = obj.f;
                generator_template.sampling_period = 1/obj.receiverBandwidth;
                generator_template.maxSamplesPerPilot = obj.maxSamplesPerPilot;
                generator = baselineGenerator2(selectedWalls, x_wall_file_name,...
                     y_wall_file_name, obj.gridSize, generator_template);
                 
                source_loc = [generator.xt_loc; generator.yt_loc];
                return
            end
            generator = MultiWallChannelGainGenerator;
            
            generator.boundary=[
                -10,48
                -10,30
                -0,3
                ];
            x1 = generator.boundary(1,1);
            x2 = generator.boundary(1,2);
            y1 = generator.boundary(2,1);
            y2 = generator.boundary(2,2);
            
            generator.f=obj.f;           % Frequency in Hz
            
            xrange=linspace(x1,x2,obj.gridSize(1));
            yrange=linspace(y1,y2,obj.gridSize(2));
            generator.xt=  [xrange(4), yrange(4)]; %Luismi: What is the meaning of this line? Is it a random point close to the edge?
            generator.xt_loc=  [0, 28, 32 ,-3, 42]; %
            generator.yt_loc=  [0, 10, -7, 30, 30]; %
            generator.ptx=-30*ones(1, size(generator.xt,1));   % Transmitter powers in dBW
            generator.ptx_loc=-30*ones(1, length(generator.xt_loc));   % Loc anchor nodes in dBW
            
            generator.zt=1.5; %transmitter Height
            generator.x1=x1; % First boundary  along x
            generator.x2=x2; % Second boundary  along x
            generator.path_loss_exp = 2;
            generator.n_gridpoints_x=obj.gridSize(1);
            generator.n_gridpoints_y=obj.gridSize(2);
            generator.zr=1.5; %receiverHeight;
            generator.sampling_period=1/obj.receiverBandwidth; %
            generator.maxSamplesPerPilot=obj.maxSamplesPerPilot;
            generator.y_limits = [y1 y2]; %
            if selectedWalls==0
                generator.estimateLocFreeSpace=1;
            else
                generator.estimateLocFreeSpace=0;
            end
            X_wallcoord=load(x_wall_file_name);
            Y_wallcoord=load(y_wall_file_name);
            arrangedCoord=[2*selectedWalls-1;2*selectedWalls];
            if length(selectedWalls)<=2
                if selectedWalls==0
                    coordToUse=(11:18)';
                else
                    coordToUse=[arrangedCoord(:);(11:18)']; % (11 to 18)-th coordinates accounts for small structure ...
                end                                        % allowing the multiwall method to work with few walls
            else
                coordToUse=arrangedCoord(:);
            end
            generator.X=X_wallcoord.X(coordToUse);
            generator.Y=Y_wallcoord.Y(coordToUse);
            
            source_loc = [generator.xt_loc; generator.yt_loc];
        end
        %changed streetlimits from [0 300] to [0 75].
        
        function [allPointEstimatedLocation, eig_vals_ratios,  allLocFreeFeatures, ...
                  evaluationGrid_x,          evaluationGrid_y, source_loc, ...
                  trueMap,              estimatedPilotSignals, SqErr_AllRuns_locFree, ...
                  SqErr_AllRuns_locBased,   meanErrOnEvalFeat, locFreeEstimate,...
                  locBasedEstimate,                UEFeatures, UELocationsPairs,...
                  averageMisses,                 NMSE_locFree, NMSE_locBased,...
                  UELocationsPairsTest,          trueGains...
                ] = run ...
                (obj, selectedWalls, combin_pairs, frac_train,frac_val, x_wall_file_name,y_wall_file_name, simulationNumber) %"do-not-touch"
            % RUN runs simulation using the parameters stored in obj.
            % run (obj, ...) runs simulation and saves generated data
            % in the savedResults folder, in a file ending by simulationNumber
            
            mySim = Simulator;
            mySim.generator = obj.baselineGenerator(selectedWalls, x_wall_file_name,y_wall_file_name);
            mySim.features=1:obj.N_feat_used;
            
            mySim.sampler = SpectrumMapSampler;
            mySim.sampler.powerNoiseSTD =obj.powerNoiseSTD;
            mySim.sampler.pilotNoiseSTD =obj.pilotNoiseSTD;
            mySim.sampler.maxSamplesPerPilot = obj.maxSamplesPerPilot;
            mySim.PCAenable=obj.PCAenabler;
            mySim.featureExtractor=FeatureExtractor;
            mySim.featureExtractor.sampling_period=mySim.generator.sampling_period;
            gaussianKernelLF = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(obj.kernelSigmaLF^2));
            gaussianKernelLB = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(obj.kernelSigmaLB^2));
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.kernel = gaussianKernelLF;
            mySim.locFreeEstimator.regularizationParameter = obj.lambdaLF;
            mySim.locFreeEstimator.enableMissingData = obj.considerMissingData;
            mySim.locFreeEstimator.completionMethodSelect=obj.completionMethodSelect;
            mySim.locFreeEstimator.desiredRank=obj.desiredRank;
            mySim.locFreeEstimator.regParEval=obj.regularParEval;
            mySim.locFreeEstimator.evalOption=obj.evalcoeffOption;
            
            mySim.locBasedEstimator = LocationBasedEstimator;
            mySim.locEstimator = DiverseLocationEstimator('IRWSRDLS'); 
            %TODO assign location estimator type (alg) based on a property
            %of the CartographySimulations object
            mySim.locBasedEstimator.kernel = gaussianKernelLB;
            mySim.locBasedEstimator.regularizationParameter = obj.lambdaLB;
            mySim.locBasedEstimator.locationNoiseSTD = obj.locationNoiseSTD;
            
            mySim.randomized_pairs = obj.randomized_pairs;
            
            rng(obj.rng_seed)
            points_filename = sprintf('allPointEstLoc%d', simulationNumber);
            points_filename2 = sprintf('allLocFreeFeatures%d', simulationNumber);
            
            fprintf('Generating the true map based on the specified grid size \n')
            if obj.simLocBased==true
                % old code:
                % fprintf('Precomputing %s and %s \n', points_filename, points_filename2);             
                [allPointEstimatedLocation,eig_vals_ratios,~, ~, ~,...
                    ~, ~, ~] = mySim.precomputeEstimatedLocations(); %#ok<NASGU>
           
                [estimatedPilotSignals, trueMap, source_loc, evaluationGrid_x, evaluationGrid_y] = ...
                    mySim.generator.generateCGMap(mySim.generator.xt);
                allLocFreeFeatures = []; %added this line just so I don't have to change the call to .simulate
                
                %                     save(points_filename, 'allPointEstimatedLocation')
                %                     save(points_filename2,'allLocFreeFeatures')
                fprintf('Precomputation of %s and %s\n is done \n', points_filename,points_filename2)
                %                 end
            else
                disp('This simulation does not involve LocBased estimation. Not precomputing location data')
                mySim.simLocBased = false;
                fprintf('Precomputating LocF features only...%s\n', points_filename2);
                
                %                     else
                %                         error('Fatal (logics)')
                %                     end
                [~,allLocFreeFeatures,~, estimatedPilotSignals, trueMap,...
                    source_loc, evaluationGrid_x, evaluationGrid_y] = mySim.precomputeEstimatedLocations(); %#ok<NASGU>
                allPointEstimatedLocation=0;
                
                %
                %                     save(points_filename2,'allLocFreeFeatures')
                
                fprintf('Precomputation of %s is done \n ', points_filename2)
            end
            
            
            if obj.computeOnlyFeatures ==true
                SqErr_AllRuns_locFree=0;
                SqErr_AllRuns_locBased=0;
                meanErrOnEvalFeat=0;
                locFreeEstimate=0;
                locBasedEstimate=0;
                UEFeatures=0;
                UELocationsPairsTest=0;
                UELocationsPairs=0;
                averageMisses=0;
                NMSE_locFree=0;
                NMSE_locBased=0;
                trueGains=0;
                return
            end
            
            
            fprintf('Will run a total of %d simulations with %d Montecarlo runs each.\n', length(obj.N_feat)*length(obj.N_ue_pairs), obj.n_runs)
            for  ind_feat=1:length(obj.N_feat)
                if obj.PCAenabler==true
                    mySim.PCAdimension=obj.N_feat(ind_feat);
                elseif obj.considerMissingData==true
                    mySim.locFreeEstimator.receiverSensitivity=obj.N_feat(ind_feat);
                elseif obj.testDiffSetOfFeatures==true
                    mySim.features=obj.featureCell{obj.N_feat_used(ind_feat)};
                end
                for ind_nue_pairs=1:length(obj.N_ue_pairs)
                    nue_pairs = obj.N_ue_pairs(ind_nue_pairs);
                    if obj.inParallel
                        parfor ind_runs=1:obj.n_runs
                            [NMSE_locFree{ind_nue_pairs, ind_runs, ind_feat}, NMSE_locBased{ind_nue_pairs, ind_runs, ind_feat},UELocationsPairsTest{ind_nue_pairs, ind_runs, ind_feat},SqErr_AllRuns_locFree{ind_nue_pairs, ind_runs, ind_feat}, SqErr_AllRuns_locBased{ind_nue_pairs, ind_runs, ind_feat}, ...
                                meanErrOnEvalFeat{ind_nue_pairs, ind_runs, ind_feat}, locFreeEstimate{ind_nue_pairs, ind_runs, ind_feat},UEFeatures{ind_nue_pairs, ind_runs, ind_feat},...
                                locBasedEstimate{ind_nue_pairs, ind_runs, ind_feat},UELocationsPairs{ind_nue_pairs, ind_runs, ind_feat},averageMisses{ind_nue_pairs, ind_runs, ind_feat}] = ... % SqErr_AllRuns_sampleAverage(ind_nue_pairs, ind_runs, ind_feat) = ...
                                mySim.simulate(nue_pairs, combin_pairs, frac_train,frac_val, allPointEstimatedLocation,trueMap,...
                                evaluationGrid_x, evaluationGrid_y, source_loc, allLocFreeFeatures, estimatedPilotSignals);
                        end
                    else
                        disp 'Not running in parallel.'
                        for ind_runs=1:obj.n_runs
                            [NMSE_locFree{ind_nue_pairs, ind_runs, ind_feat}, NMSE_locBased{ind_nue_pairs, ind_runs, ind_feat},UELocationsPairsTest{ind_nue_pairs, ind_runs, ind_feat},SqErr_AllRuns_locFree{ind_nue_pairs, ind_runs, ind_feat}, SqErr_AllRuns_locBased{ind_nue_pairs, ind_runs, ind_feat}, ...
                                meanErrOnEvalFeat{ind_nue_pairs, ind_runs, ind_feat}, locFreeEstimate{ind_nue_pairs, ind_runs, ind_feat},UEFeatures{ind_nue_pairs, ind_runs, ind_feat},...
                                locBasedEstimate{ind_nue_pairs, ind_runs, ind_feat},UELocationsPairs{ind_nue_pairs, ind_runs, ind_feat},averageMisses{ind_nue_pairs, ind_runs, ind_feat}, ...
                                trueGains] = ... % SqErr_AllRuns_sampleAverage(ind_nue_pairs, ind_runs, ind_feat) = ...
                                mySim.simulate(nue_pairs, combin_pairs, frac_train,frac_val, [] ,trueMap,...
                                evaluationGrid_x, evaluationGrid_y, source_loc, allLocFreeFeatures, estimatedPilotSignals);
                        end
                    end
                    fprintf('Finished simulation for N=%d, M=%d\n', nue_pairs, obj.N_feat(ind_feat));
                end
            end
            %             trueMapAverage = mean(mean(trueMap));
            %             NMSE_locFree = mean(SqErr_AllRuns_locFree ,2)./sum(sum((channelForPairsVal-trueMapAverage).^2)); %#ok<NOPRT>
            %             NMSE_locBased= (mean(SqErr_AllRuns_locBased,2)./(sum(sum((channelForPairsVal-trueMapAverage).^2)))); %#ok<NOPRT>
            %             beep
        end
    end %methods
end

