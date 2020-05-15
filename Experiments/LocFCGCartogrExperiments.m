classdef LocFCGCartogrExperiments < ExperimentFunctionSet
    
    properties
        
        b_loadData = 0 % set this property to 1 if you want the experiment
        % methods to load data from the savedResults folder.
        
        % You may define here constants that you will use in the
        % experiments of this file
        
    end
    
    methods
        
        % This experiment implements cross-validations for the choice of
        % parameters
        
        function F=experiment_101(obj,niter) 
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [20 20];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue_pairs=500; % number of sensor pairs
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =false;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                %Cross validation starts here
                frac_train=0.6; % percentage of pairs to use for training
                frac_val=0.2; % percentage of pairs to use for validation
                n_points_sig_lam=50;
                kernelSigmaLF=linspace(1e-2,1e2,n_points_sig_lam);
                lambdaLF=linspace(1e-7,3e-3,n_points_sig_lam);
                kernelSigmaLB=linspace(1e-2,1e2,n_points_sig_lam);
                lambdaLB=linspace(1e-7,3e-3,n_points_sig_lam);
                
                NMSE_LF=zeros(length(kernelSigmaLF), length(lambdaLF));
                NMSE_LB=zeros(length(kernelSigmaLB), length(lambdaLB));
                ltm =LoopTimeControl(length(kernelSigmaLF)*length(lambdaLF));
                for ind_sig=1:length(kernelSigmaLF)
                    for ind_lam=1:length(lambdaLF)
                        c.kernelSigmaLF =kernelSigmaLF(ind_sig);
                        c.kernelSigmaLB =kernelSigmaLB(ind_sig);
                        c.lambdaLF=lambdaLF(ind_lam);
                        c.lambdaLB=lambdaLB(ind_lam);
                        [~,~, ~,~, ~, ~,~,...
                            ~,~,~,~,~,~,~,...
                            ~, NMSE_locFree,NMSE_locBased,~]=c.run(selectedWalls,combin_pairs, frac_train,frac_val,...
                            x_wall_file_name,y_wall_file_name,simulationNumber);
                        
                        NMSE_LF(ind_sig,ind_lam)=NMSE_locFree;
                        NMSE_LB(ind_sig,ind_lam)=NMSE_locBased;
                        ltm.go(ind_lam+(ind_sig-1)*length(lambdaLF));
                    end
                end
                
                [min_nmse_locF, J_lf]=min(NMSE_LF(:));
                [ind_sigOpt_lf, ind_lamOpt_lf]=ind2sub(size(NMSE_LF),J_lf);
                min_nmse_locF
                kernelSigmaLF_Opt=kernelSigmaLF(ind_sigOpt_lf)
                lambdaLF_Opt=lambdaLF(ind_lamOpt_lf)
                
                [min_nmse_locB, J_lb]=min(NMSE_LB(:));
                [ind_sigOpt_lb, ind_lamOpt_lb]=ind2sub(size(NMSE_LB),J_lb);
                min_nmse_locB
                kernelSigmaLB_Opt=kernelSigmaLB(ind_sigOpt_lb)
                lambdaLB_Opt=lambdaLB(ind_lamOpt_lb)
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            F = [];
        end

        function F=experiment_1011(obj,niter) %! same as 101, but toy size
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [4 4]; % was 20, 20
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = 10; %!niter;% 10
                
                c.N_ue_pairs=50; % number of sensor pairs %! was 500
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =false;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                %Cross validation starts here
                frac_train=0.6; % percentage of pairs to use for training
                frac_val=0.2; % percentage of pairs to use for validation
                n_points_sig_lam=5; %! was 50
                kernelSigmaLF=linspace(1e-2,1e2,n_points_sig_lam);
                lambdaLF=linspace(1e-7,3e-3,n_points_sig_lam);
                kernelSigmaLB=linspace(1e-2,1e2,n_points_sig_lam);
                lambdaLB=linspace(1e-7,3e-3,n_points_sig_lam);
                
                NMSE_LF=zeros(length(kernelSigmaLF), length(lambdaLF));
                NMSE_LB=zeros(length(kernelSigmaLB), length(lambdaLB));
                ltm =LoopTimeControl(length(kernelSigmaLF)*length(lambdaLF));
                for ind_sig=1:length(kernelSigmaLF)
                    for ind_lam=1:length(lambdaLF)
                        c.kernelSigmaLF =kernelSigmaLF(ind_sig);
                        c.kernelSigmaLB =kernelSigmaLB(ind_sig);
                        c.lambdaLF=lambdaLF(ind_lam);
                        c.lambdaLB=lambdaLB(ind_lam);
                        [~,~, ~,~, ~, ~,~,...
                            ~,~,~,~,~,~,~,...
                            ~, ~, NMSE_locFree,NMSE_locBased,~]=c.run(selectedWalls,combin_pairs, frac_train,frac_val,...
                            x_wall_file_name,y_wall_file_name,simulationNumber);
                        
                        NMSE_LF(ind_sig,ind_lam)=NMSE_locFree;
                        NMSE_LB(ind_sig,ind_lam)=NMSE_locBased;
                        ltm.go(ind_lam+(ind_sig-1)*length(lambdaLF));
                    end
                end
                
                [min_nmse_locF, J_lf]=min(NMSE_LF(:));
                [ind_sigOpt_lf, ind_lamOpt_lf]=ind2sub(size(NMSE_LF),J_lf);
                min_nmse_locF
                kernelSigmaLF_Opt=kernelSigmaLF(ind_sigOpt_lf)
                lambdaLF_Opt=lambdaLF(ind_lamOpt_lf)
                
                [min_nmse_locB, J_lb]=min(NMSE_LB(:));
                [ind_sigOpt_lb, ind_lamOpt_lb]=ind2sub(size(NMSE_LB),J_lb);
                min_nmse_locB
                kernelSigmaLB_Opt=kernelSigmaLB(ind_sigOpt_lb)
                lambdaLB_Opt=lambdaLB(ind_lamOpt_lb)
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            F = [];
        end

        
        % This experiment tests the optimal parameters obtained in experiment 101
        function F=experiment_102(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [10 10];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=[1:2];
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue_pairs=1000; % number of sensor pairs
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =false;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                frac_train=0.6; % percentage of pairs to use for training
                frac_val=0.4; % percentage of pairs to use for testing here
                
                c.kernelSigmaLF =100;
                c.kernelSigmaLB =59;
                c.lambdaLF=9.1844e-4;
                c.lambdaLB=1.3e-3;
                [~,~, ~,~, ~, ~,~,...
                    ~,~,~,~,~,~,~,...
                    ~, ~, NMSE_locFree,NMSE_locBased,UELocationsPairsTest]=c.run(selectedWalls,combin_pairs, frac_train,frac_val,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            F = [];
        end
        
        % This experiment plots the LocF and LocB NMSE vs the number of sensor pairs
        % 60% of the sensor pairs are used for training and the remaining
        % 40 % are used for testing
        function F=experiment_103(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [20 20];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue_pairs=1e2:1e2:2e3; % number of sensor pairs
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =true;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                frac_train=0.6; % percentage of pairs to use for training
                frac_val=0.4; % percentage of pairs to use for testing here
                
                c.kernelSigmaLF =100;
                c.kernelSigmaLB =59;
                c.lambdaLF=9.1844e-4;
                c.lambdaLB=1.3e-3;
                [~,~, ~,~, ~, ~,~,...
                    ~,~,~,~,~,~,~,...
                    ~, NMSE_locFree,NMSE_locBased,~]=c.run(selectedWalls,combin_pairs, frac_train,frac_val,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            c_legend{1} = 'LocB';
            c_styles{1} = '>r-';
            c_legend{2} = 'LocF';
            c_styles{2} = 'sb-';
            
            m_X = c.N_ue_pairs;
            m_Y =[mean(cell2mat(NMSE_locBased),2)'; mean(cell2mat(NMSE_locFree),2)'];
            F = GFigure('m_X',m_X,'m_Y',m_Y,'ch_xlabel','Number of sensor locations pairs, N','ch_ylabel','NMSE',...
                'c_legend',c_legend, 'c_styles', c_styles);
        end
        
        % This experiment runs a single instance of the training and validation,
        % putting a grid in the validation dataset with the goal of plotting the obtained maps.
        % It uses the optimal parameters obtained in experiment 101
        function F=experiment_201(obj,~)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [15 15]; %! [20 20];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = 1
                
                % c.N_ue_pairs=1000; % number of sensor pairs %! defined
                % after creating the combin_pairs variable
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =false;
                % combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2); %!
                % changed (see section "defining a trajectory")
                frac_train=0.06; % percentage of pairs to use for training
                frac_val=0; % percentage of pairs to use for testing here
                
                c.kernelSigmaLF =100;
                c.kernelSigmaLB =59;
                c.lambdaLF=9.1844e-4;
                c.lambdaLB=1.3e-3;
                
                %% Luismi:defining a trajectory
                column_begin= 4;
                rowskip = 1;
                Ntraj = ceil(c.gridSize(2)./rowskip)-1;
                final_index = column_begin+rowskip*Ntraj*c.gridSize(2);
                v_tile_sequence = column_begin:(rowskip*c.gridSize(1)):final_index;
                
                t_combin_pairs = zeros(c.gridSize(1)*c.gridSize(2),length(v_tile_sequence),2);
                v_all_points = 1:c.gridSize(1)*c.gridSize(2);
                for index = 1:length(v_tile_sequence)
                    t_combin_pairs(:, index, 1) = v_tile_sequence(index);
                    t_combin_pairs(:, index, 2) = v_all_points;
                end
                combin_pairs = reshape(t_combin_pairs, c.gridSize(1)*c.gridSize(2)*length(v_tile_sequence),2);
                c.randomized_pairs = 0;
                c.N_ue_pairs= length(v_tile_sequence)*length(v_all_points);
                [~,~,~, ~,~,~, ~,~,~ ,~,locFreeEstimate,locBasedEstimate, ~,~,~, ...
                    NMSE_locFree,NMSE_locBased,UELocationsPairsTest,  trueGains  ...
                    ]=c.run(selectedWalls,combin_pairs, frac_train,frac_val,    ...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                M_trueMap = create_animation2([combin_pairs trueGains' ],          c.gridSize(1), c.gridSize(2), 'trueMap' );
                M_locFree = create_animation2([combin_pairs locFreeEstimate{1}' ], c.gridSize(1), c.gridSize(2), 'locFree' );
                M_locBased= create_animation2([combin_pairs locBasedEstimate{1}'], c.gridSize(1), c.gridSize(2), 'locBased');
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            F = [];
        end
        
        function F=experiment_203(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [ 20 20 ];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = 1
                
                c.N_ue_pairs=1000; % number of sensor pairs
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =false;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                
                frac_train=1; % percentage of pairs to use for training
                frac_val=0; % percentage of pairs to use for testing here
                
                c.kernelSigmaLF =69;
                c.kernelSigmaLB =43;
                c.lambdaLF=1e-5;
                c.lambdaLB=6e-3;
                
                %% Luismi:defining a trajectory
                column_begin= 8;
                rowskip = 1;
                Ntraj = ceil(c.gridSize(2)./rowskip)-1;
                final_index = column_begin+rowskip*Ntraj*c.gridSize(2);
                v_tile_sequence = column_begin:(rowskip*c.gridSize(1)):final_index;
                
                t_val_pairs = zeros(c.gridSize(1)*c.gridSize(2),length(v_tile_sequence),2);
                v_all_points = 1:c.gridSize(1)*c.gridSize(2);
                for index = 1:length(v_tile_sequence)
                    t_val_pairs(:, index, 1) = v_tile_sequence(index);
                    t_val_pairs(:, index, 2) = v_all_points;
                end
                pairs_val = reshape(t_val_pairs, c.gridSize(1)*c.gridSize(2)*length(v_tile_sequence),2);
                
                c.randomized_pairs = 0;
                [~,~,~, ~,~,~, ~, ...
                    ~,~ ,~,locFreeEstimate,locBasedEstimate, ~,~,...
                    ~, ~, NMSE_locFree,NMSE_locBased,UELocationsPairsTest,  trueGains  ...
                    ]=c.run(selectedWalls,combin_pairs, frac_train, pairs_val,   ...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                %                 M_trueMap = create_animation2([combin_pairs trueGains' ],          c.gridSize(1), c.gridSize(2), 'trueMap' );
                %                 M_locFree = create_animation2([combin_pairs locFreeEstimate{1}' ], c.gridSize(1), c.gridSize(2), 'locFree' );
                %                 M_locBased= create_animation2([combin_pairs locBasedEstimate{1}'], c.gridSize(1), c.gridSize(2), 'locBased');
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            F = [];
        end

        %same as 203
        function F=experiment_204(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [ 20 20 ];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = 1
                
                c.N_ue_pairs=1000; % number of sensor pairs
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =false;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                
                frac_train=1; % percentage of pairs to use for training
                frac_val=0; % percentage of pairs to use for testing here
                
                c.kernelSigmaLF =69;
                c.kernelSigmaLB =43;
                c.lambdaLF=1e-5;
                c.lambdaLB=6e-3;
                
                %% Luismi:defining a trajectory
                column_begin= 8;
                rowskip = 1;
                Ntraj = ceil(c.gridSize(2)./rowskip)-1;
                final_index = column_begin+rowskip*Ntraj*c.gridSize(2);
                v_tile_sequence = column_begin:(rowskip*c.gridSize(1)):final_index;
                
                t_val_pairs = zeros(c.gridSize(1)*c.gridSize(2),length(v_tile_sequence),2);
                v_all_points = 1:c.gridSize(1)*c.gridSize(2);
                for index = 1:length(v_tile_sequence)
                    t_val_pairs(:, index, 1) = v_tile_sequence(index);
                    t_val_pairs(:, index, 2) = v_all_points;
                end
                pairs_val = reshape(t_val_pairs, c.gridSize(1)*c.gridSize(2)*length(v_tile_sequence),2);
                
                c.randomized_pairs = 0;
                [~,~,~, ~,~,~, ~, ...
                 ~,~ ,~,locFreeEstimate,locBasedEstimate, ~,~,...
                 ~, ~, NMSE_locFree,NMSE_locBased,UELocationsPairsTest,  trueGains  ...
                    ]=c.run(selectedWalls,combin_pairs, frac_train, pairs_val,   ...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                %                 M_trueMap = create_animation2([combin_pairs trueGains' ],          c.gridSize(1), c.gridSize(2), 'trueMap' );
                %                 M_locFree = create_animation2([combin_pairs locFreeEstimate{1}' ], c.gridSize(1), c.gridSize(2), 'locFree' );
                %                 M_locBased= create_animation2([combin_pairs locBasedEstimate{1}'], c.gridSize(1), c.gridSize(2), 'locBased');
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            F = [];
        end

        
        function F = experiment_301(obj, niter)
            % simulation using the simplified simulator class (Simulator2)
            % Same inputs as exp 203
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=1:5;
            
            kernelSigmaLF = 69;
            kernelSigmaLB = 43;
            lambdaLF = 1e-5;
            lambdaLB = 6e-3;
            
            gridSize = [40 40];
            N_pair_train = 25000; % 1%
            
            c = CartographySimulations; % we should be able to remove
            % this object from here
            c.f = 800e6;
            c.gridSize = gridSize;
            c.receiverBandwidth = 20e6;
            
            mySim = Simulator2;
            [mySim.generator, source_loc] = c.baselineGenerator(...
                selectedWalls, x_wall_file_name, y_wall_file_name);
            mySim.generator.maxSamplesPerPilot = 10;
            
            % mySim.features = 1:size(combnk(1:numberOfSources,2),1);
            % we could move .baselineGenerator to the Simulator[2] class
            mySim.sampler = SpectrumMapSampler;
            mySim.sampler.pilotNoiseSTD = 1e-5;
            mySim.sampler.powerNoiseSTD = 0.5;
            mySim.sampler.maxSamplesPerPilot=10;
            %mySim.PCAenable=obj.PCAenabler;
            mySim.featureExtractor = FeatureExtractor;
            mySim.featureExtractor.sampling_period = mySim.generator.sampling_period;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
            mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            
            mySim.locBasedEstimator = LocationBasedEstimator;
            mySim.locEstimator = DiverseLocationEstimator('IRWSRDLS');
            mySim.locBasedEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
            mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            mySim.locBasedEstimator.locationNoiseSTD = []; %?
            mySim.locBasedEstimator.Xenb = source_loc;
            
            mySim.locEstimator.Xenb = source_loc;
            
            % The original CartographySimulations has several lines here to
            % define the number of different monte carlo runs and define the
            % random parameters or each MC run
            
            % define grid:
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(mySim.generator.x2,          mySim.generator.x1,           mySim.generator.n_gridpoints_x), ...
                linspace(mySim.generator.y_limits(1), mySim.generator.y_limits(2),  mySim.generator.n_gridpoints_y));
            
            %% train and evaluation(traceMap) pairs
            combin_pairs  = combnk(1:gridSize(1)*gridSize(2),2);
            m_train_pairs = combin_pairs(randperm(size(combin_pairs, 1), N_pair_train),:);
            %TODO: problem here = m_train_pairs should be Nx2, but it is Nx1
            
            % defining a trajectory
            column_begin= 8;
            rowskip = 2;
            Ntraj = ceil(gridSize(2)./rowskip)-1;
            final_index = column_begin+rowskip*Ntraj*gridSize(2);
            v_tile_sequence = column_begin:(rowskip*gridSize(1)):final_index;
            
            t_val_pairs = zeros(gridSize(1)*gridSize(2),length(v_tile_sequence),2);
            v_all_points = 1:gridSize(1)*gridSize(2);
            for index = 1:length(v_tile_sequence)
                t_val_pairs(:, index, 1) = v_tile_sequence(index);
                t_val_pairs(:, index, 2) = v_all_points;
            end
            m_tracemap_pairs = reshape(t_val_pairs, gridSize(1)*gridSize(2)*length(v_tile_sequence),2);
            
            %%
            disp 'Running only one experiment (no monte carlo yet)'
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                m_train_pairs, m_tracemap_pairs, m_grid_x, m_grid_y);
            locFreeMapEstimate  = str_mapEstimates.locFree;
            locBasedMapEstimate = str_mapEstimates.locBased;
            
            M_trueMap = create_animation2([m_tracemap_pairs trueGains' ],  ...
                gridSize(1), gridSize(2), 'trueMap' );
            M_locFree = create_animation2([m_tracemap_pairs locFreeMapEstimate' ],...
                gridSize(1), gridSize(2), 'locFree' );
            M_locBased = create_animation2([m_tracemap_pairs locBasedMapEstimate' ],...
                gridSize(1), gridSize(2), 'locBased' );
            
            M = create_animation3([m_tracemap_pairs trueGains' locFreeMapEstimate' ...
                locBasedMapEstimate'], gridSize(1), gridSize(2), ...
                {'True', 'LocFree', 'LocBased'}, 's40x40')
            
            F = GFigure.captureCurrentFigure;
        end
        
        % experiment using simulator2 class and hybrid simulator for
        % visualizing channel gains 
        function F = experiment_302(obj, niter)
            % using Simulator2
            % Similar inputs to exp 203
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=[2 3];
            
            kernelSigmaLF = 69;
            kernelSigmaLB = 35;
            lambdaLF      = 1e-5;
            lambdaLB      = 8e-3;
        
            gridSize = [20 20];
            N_pair_train = 100; %
            
            c = CartographySimulations; % we should be able to remove
            % this object from here
            c.f = 800e6;
            c.gridSize = gridSize;
            c.receiverBandwidth = 20e6;
            
            mySim = Simulator2;
            [mySim.generator, source_loc] = c.baselineGenerator(...
                selectedWalls, x_wall_file_name, y_wall_file_name);
            mySim.generator.maxSamplesPerPilot = 10;
            
            % mySim.features = 1:size(combnk(1:numberOfSources,2),1);
            % we could move .baselineGenerator to the Simulator[2] class
            mySim.sampler = SpectrumMapSampler;
            mySim.sampler.pilotNoiseSTD = 1e-5;
            mySim.sampler.powerNoiseSTD = 0.5;
            mySim.sampler.maxSamplesPerPilot=10;
            %mySim.PCAenable=obj.PCAenabler;
            mySim.featureExtractor = FeatureExtractor;
            mySim.featureExtractor.sampling_period = mySim.generator.sampling_period;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
            mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            
            mySim.locBasedEstimator = LocationBasedEstimator;
            mySim.locEstimator = WangLocationEstimator;
            mySim.locBasedEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
            mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            mySim.locBasedEstimator.locationNoiseSTD = []; %?
            mySim.locBasedEstimator.Xenb = source_loc;
            
            mySim.locEstimator.Xenb = source_loc;
            
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
            mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
            mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
            mySim.hybridEstimator.regularizationParameterLB =lambdaLB;
            
            % The original CartographySimulations has several lines here to
            % define the number of different monte carlo runs and define the
            % random parameters or each MC run
            
            % define grid:
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(mySim.generator.x2,          mySim.generator.x1,           mySim.generator.n_gridpoints_x), ...
                linspace(mySim.generator.y_limits(1), mySim.generator.y_limits(2),  mySim.generator.n_gridpoints_y));
            
            %% train and evaluation(traceMap) pairs
            combin_pairs  = combnk(1:gridSize(1)*gridSize(2),2);
            m_train_pairs = combin_pairs(randperm(size(combin_pairs, 1), N_pair_train),:);
            %TODO: problem here = m_train_pairs should be Nx2, but it is Nx1
            
            % defining a trajectory
            column_begin= 1;
            rowskip = 2;
            Ntraj = ceil(gridSize(2)./rowskip)-1;
            final_index = column_begin+rowskip*Ntraj*gridSize(2);
            v_tile_sequence = column_begin:(rowskip*gridSize(1)):final_index;
            
            t_val_pairs = zeros(gridSize(1)*gridSize(2),length(v_tile_sequence),2);
            v_all_points = 1:gridSize(1)*gridSize(2);
            for index = 1:length(v_tile_sequence)
                t_val_pairs(:, index, 1) = v_tile_sequence(index);
                t_val_pairs(:, index, 2) = v_all_points;
            end
            m_tracemap_pairs = reshape(t_val_pairs, gridSize(1)*gridSize(2)*length(v_tile_sequence),2);
            
            %%
            disp 'Running only one experiment (no monte carlo yet)'
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                m_train_pairs, m_tracemap_pairs, m_grid_x, m_grid_y);
            locFreeMapEstimate  = str_mapEstimates.locFree;
            locBasedMapEstimate = str_mapEstimates.locBased;
            
            M_trueMap = create_animation2([m_tracemap_pairs trueGains' ],  ...
                gridSize(1), gridSize(2), 'trueMap' );
            M_locFree = create_animation2([m_tracemap_pairs locFreeMapEstimate' ],...
                gridSize(1), gridSize(2), 'locFree' );
            M_locBased = create_animation2([m_tracemap_pairs locBasedMapEstimate' ],...
                gridSize(1), gridSize(2), 'locBased' );
%             M_hybrid = create_animation2([m_tracemap_pairs hybridMapEstimate' ],...
%                 gridSize(1), gridSize(2), 'hybrid' );             
            M = create_animation3([m_tracemap_pairs trueGains' locFreeMapEstimate' ...
                locBasedMapEstimate' ], gridSize(1), gridSize(2), ...
                {'True', 'LocFree', 'LocBased'}, 's40x40') %hybridMapEstimate'  , 'hybrid'
            NMSE_locFree
            NMSE_locBased
%             NMSE_hybrid
            F = GFigure.captureCurrentFigure;
        end
        
        % experiment for obtaining localization errors along with the rank of
        % matrix Y using the robust SDR algorithm by Wang et al.
        function F=experiment_303(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(37:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found or the flag for loading data is disabled. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.gridSize = [5 5];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = true;
                
                selectedWalls=1:2;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue_pairs=1e2:1e2:2e3; % number of sensor pairs
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                
                c.inParallel =true;
                combin_pairs=combnk(1:c.gridSize(1)*c.gridSize(2),2);
                frac_train=0.6; % percentage of pairs to use for training
                frac_val=0.4; % percentage of pairs to use for testing here
                
                c.kernelSigmaLF =100;
                c.kernelSigmaLB =59;
                c.lambdaLF=9.1844e-4;
                c.lambdaLB=1.3e-3;
                [allPointEstimatedLocation,eig_vals_ratios,~, evaluationGrid_x,evaluationGrid_y, ~, ~,~,...
                    ~,~,~,~,~,~,~,...
                    ~, ~,~,~,~]=c.run(selectedWalls,combin_pairs, frac_train,frac_val,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            allTrueLocations=[flipud(evaluationGrid_x(:)'); flipud(evaluationGrid_y(:)')];
            start_point_plot=size(allTrueLocations, 2)+1;
            allPointEstimatedLocation = allPointEstimatedLocation(:,:,2);
            errors = vecnorm(allTrueLocations-allPointEstimatedLocation);
            figure(2); 
            subplot(2,1,1)
            plot(errors, 'bd-','linewidth', 2); grid on
            xlabel('Grid point index, $k$','Interpreter','latex')
            ylabel('$\Vert \boldmath{x}_k- \hat{\boldmath{x}}_k\Vert $','Interpreter','latex')
            title('Robust SDR')
            subplot(2,1,2)
            plot(eig_vals_ratios(start_point_plot:end), 'r*-.', 'linewidth', 2); grid on
            xlabel('Grid point index, $k$','Interpreter','latex')
            ylabel('$\frac{\lambda_1}{\lambda_2}$','Interpreter','latex')
            title('Ratio of the 1st and 2nd eig. values')
            min_ration=min(eig_vals_ratios(start_point_plot:end))
            F = [];
         
        end
        
        % plot histograms of the difference between range diff and
        % LocBasedFeatures
        function F = experiment_304(obj, niter)
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=1:4;
           
            gridSize = [20 20]; % 20 20
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            maxSamplesPerPilot = 10;
            
            myBiasErrorSim = BiasErrorSimulator;
           
%             c = CartographySimulations; % we should be able to remove
%             % this object from here
%             c.f = carrier_frequency;
%             c.gridSize = gridSize;
%             c.receiverBandwidth = receiverBandwidth;
%             c.b_new_call = 1;
%             
%             [myBiasErrorSim.generator, source_loc] = c.baselineGenerator(...
%                 selectedWalls, x_wall_file_name, y_wall_file_name);
%             myBiasErrorSim.generator.maxSamplesPerPilot = 10;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = 1/receiverBandwidth;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [myBiasErrorSim.generator, m_source_loc] = baselineGenerator2(...
                selectedWalls, x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            myBiasErrorSim.generator.delay_estimation_offset = 2*rand;
            
%             assert(isequal(generator2, myBiasErrorSim.generator))
%             keyboard
            
            % myBiasErrorSim.features = 1:size(combnk(1:numberOfSources,2),1);
            % we could move .baselineGenerator to the Simulator[2] class
            myBiasErrorSim.sampler = SpectrumMapSampler;
            myBiasErrorSim.sampler.pilotNoiseSTD = 0.3e-5;
            myBiasErrorSim.sampler.powerNoiseSTD = 0.5;
            myBiasErrorSim.sampler.maxSamplesPerPilot=10;
            %myBiasErrorSim.PCAenable=obj.PCAenabler;
            
            myBiasErrorSim.featureExtractor = FeatureExtractor;
            myBiasErrorSim.featureExtractor.sampling_period = myBiasErrorSim.generator.sampling_period;
        
            myBiasErrorSim.locEstimator = WangLocationEstimator;
            myBiasErrorSim.locEstimator.param_rho = 10;
            myBiasErrorSim.locEstimator.inParallel = 1;
            
            myBiasErrorSim.locEstimator.Xenb = m_source_loc;
      
            % define grid:
            % TODO: can be moved to the BiasErrorSimulator.simulate method
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(myBiasErrorSim.generator.x2,          myBiasErrorSim.generator.x1,           myBiasErrorSim.generator.n_gridpoints_x), ...
                linspace(myBiasErrorSim.generator.y_limits(1), myBiasErrorSim.generator.y_limits(2),  myBiasErrorSim.generator.n_gridpoints_y));
           
            disp 'Running only one experiment (no monte carlo yet)'
            [featDiff, ~, m_trueTDoAs, m_allLocBasedFeatures, ...
                m_estim_loc_tm_rx, m_trueLoc_all_points] =...
                myBiasErrorSim.simulate(m_grid_x, m_grid_y);
            
            %% Plot the histogram and CDF of the feature difference
            figure(999); clf           
                nbins = 100;
                [counts, bins] = hist(vec(abs(featDiff)), nbins);
                histogram(vec(abs(featDiff)), nbins);
                xlabel '|d - r|'
                ylabel count
                hold on
            F = GFigure.captureCurrentFigure();
            figure(998); clf
                plot(bins, cumsum(counts)/sum(counts));
                xlabel '|d - r|'
                ylabel cdf
            F(2) = GFigure.captureCurrentFigure();
            
            figure(997); clf
                scatter(m_allLocBasedFeatures(:), m_trueTDoAs(:));
                xlabel 'Estimated TDoA from xcorr'
                ylabel 'True TDoA'
                
            F(3) = GFigure.captureCurrentFigure();
            
            v_err_loc = norms(m_estim_loc_tm_rx - m_trueLoc_all_points);
            rms(v_err_loc)
            %histogram(v_err_loc)

            figure(996); clf 
                for i = 1:size(m_estim_loc_tm_rx, 2)
                    plot([m_estim_loc_tm_rx(1, i) m_trueLoc_all_points(1,i)],...
                         [m_estim_loc_tm_rx(2, i) m_trueLoc_all_points(2,i)])
                     hold on
                end
            F(4) = GFigure.captureCurrentFigure();
        end
           
        % plot histograms of the difference between range diff and
        % LocBasedFeatures
        % same as 304, but using ARSRobustLocationEstimator
        function F = experiment_3050(obj, niter)
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=1:4;
           
            gridSize = [20 20]; % 20 20
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            maxSamplesPerPilot = 10;
            
            myBiasErrorSim = BiasErrorSimulator;
           
%             c = CartographySimulations; % we should be able to remove
%             % this object from here
%             c.f = carrier_frequency;
%             c.gridSize = gridSize;
%             c.receiverBandwidth = receiverBandwidth;
%             c.b_new_call = 1;
%             
%             [myBiasErrorSim.generator, source_loc] = c.baselineGenerator(...
%                 selectedWalls, x_wall_file_name, y_wall_file_name);
%             myBiasErrorSim.generator.maxSamplesPerPilot = 10;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = 1/receiverBandwidth;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [myBiasErrorSim.generator, m_source_loc] = baselineGenerator2(...
                selectedWalls, x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            myBiasErrorSim.generator.delay_estimation_offset = 2*rand;
            
%             assert(isequal(generator2, myBiasErrorSim.generator))
%             keyboard
            
            % myBiasErrorSim.features = 1:size(combnk(1:numberOfSources,2),1);
            % we could move .baselineGenerator to the Simulator[2] class
            myBiasErrorSim.sampler = SpectrumMapSampler;
            myBiasErrorSim.sampler.pilotNoiseSTD = 0.3e-5;
            myBiasErrorSim.sampler.powerNoiseSTD = 0.5;
            myBiasErrorSim.sampler.maxSamplesPerPilot=10;
            %myBiasErrorSim.PCAenable=obj.PCAenabler;
            
            myBiasErrorSim.featureExtractor = FeatureExtractor;
            myBiasErrorSim.featureExtractor.sampling_period = myBiasErrorSim.generator.sampling_period;
        
            myBiasErrorSim.locEstimator = ARSRobustLocationEstimator;
            myBiasErrorSim.locEstimator.param_rho = 10;
            myBiasErrorSim.locEstimator.inParallel = 0; %! for debugging
            
            myBiasErrorSim.locEstimator.Xenb = m_source_loc;
      
            % define grid:
            % TODO: can be moved to the BiasErrorSimulator.simulate method
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(myBiasErrorSim.generator.x2,          myBiasErrorSim.generator.x1,           myBiasErrorSim.generator.n_gridpoints_x), ...
                linspace(myBiasErrorSim.generator.y_limits(1), myBiasErrorSim.generator.y_limits(2),  myBiasErrorSim.generator.n_gridpoints_y));
           
            disp 'Running only one experiment (no monte carlo yet)'
            [featDiff, ~, m_trueTDoAs, m_allLocBasedFeatures, ...
                m_estim_loc_tm_rx, m_trueLoc_all_points] =...
                myBiasErrorSim.simulate(m_grid_x, m_grid_y);
            
            %% Plot the histogram and CDF of the feature difference
            figure(999); clf           
                nbins = 100;
                [counts, bins] = hist(vec(abs(featDiff)), nbins);
                histogram(vec(abs(featDiff)), nbins);
                xlabel '|d - r|'
                ylabel count
                hold on
            F = GFigure.captureCurrentFigure();
            
            figure(998); clf
                plot(bins, cumsum(counts)/sum(counts));
                xlabel '|d - r|'
                ylabel cdf
            F(2) = GFigure.captureCurrentFigure();
            
            figure(997); clf
                scatter(m_allLocBasedFeatures(:), m_trueTDoAs(:));
                xlabel 'Estimated TDoA from xcorr'
                ylabel 'True TDoA'               
            F(3) = GFigure.captureCurrentFigure();
            
            v_err_loc = norms(m_estim_loc_tm_rx - m_trueLoc_all_points);
            rms(v_err_loc)
            %histogram(v_err_loc)

            figure(996); clf 
                for i = 1:size(m_estim_loc_tm_rx, 2)
                    plot([m_estim_loc_tm_rx(1, i) m_trueLoc_all_points(1,i)],...
                         [m_estim_loc_tm_rx(2, i) m_trueLoc_all_points(2,i)])
                     hold on
                end
            title 'Each line joins a true location and its estimate'
            F(4) = GFigure.captureCurrentFigure();
            
            figure(995); clf
              subplot(1, 2, 1);
                scatter(m_estim_loc_tm_rx(1, :), m_trueLoc_all_points(1,:));
              subplot(1, 2, 2);  
                scatter(m_estim_loc_tm_rx(2, :), m_trueLoc_all_points(2,:));
            F(5) = GFigure.captureCurrentFigure();
        end
        
        function F = experiment_4010(obj, niter)
            % copied from exp 302.
            % Added training , crossval and evaluation of hybrid
            % method with additive curve fitting.
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=2:4; %!
            
            kernelSigmaLF = 69;
            kernelSigmaLB = 35;
            lambdaLF      = 1e-3;
            lambdaLB      = 1e-3;
        
            gridSize = [20 20]; % 20 20 
            N_pair_train = 100;
            N_pair_test  = 100;
            
            mySim = Simulator2; 
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_syntheticLocError = 1;
            mySim.b_syntheticLocEstimate = 1;
            
            %% this section was adapted from 304 
            % (removal of CartographySimulations)
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            maxSamplesPerPilot = 10;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = 1/receiverBandwidth;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [mySim.generator, m_source_loc] = baselineGenerator2(...
                selectedWalls, x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            mySim.generator.delay_estimation_offset = 2*rand;
           
            %%
            % mySim.features = 1:size(combnk(1:numberOfSources,2),1);
            % we could move .baselineGenerator to the Simulator[2] class
            mySim.sampler = SpectrumMapSampler;
            mySim.sampler.pilotNoiseSTD = 1e-5;
            mySim.sampler.powerNoiseSTD = 0.5;
            mySim.sampler.maxSamplesPerPilot=10;
            %mySim.PCAenable=obj.PCAenabler;
            mySim.featureExtractor = FeatureExtractor;
            mySim.featureExtractor.sampling_period = mySim.generator.sampling_period;
            
            mySim.b_cv = 1;
            mySim.b_cvLambdas_hybrid = 1;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
            mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            
            mySim.locBasedEstimator = LocationBasedEstimator;
            mySim.locEstimator = WangLocationEstimator;
            mySim.locBasedEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
            mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            mySim.locBasedEstimator.locationNoiseSTD = []; %?
            mySim.locBasedEstimator.Xenb = m_source_loc;
            
            mySim.locEstimator.Xenb = m_source_loc;
            
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
            mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
            mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
            mySim.hybridEstimator.regularizationParameterLB =lambdaLB;
            
            % The original CartographySimulations has several lines here to
            % define the number of different monte carlo runs and define the
            % random parameters or each MC run
            
            % define grid:
            x1 = mySim.generator.boundary(1,1);   x2 = mySim.generator.boundary(1,2);
            y1 = mySim.generator.boundary(2,1);   y2 = mySim.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x2, x1, mySim.generator.n_gridpoints_x), ...
                linspace(y1, y2, mySim.generator.n_gridpoints_y));
            
            %% train and evaluation(traceMap) pairs
            combin_pairs  = combnk(1:gridSize(1)*gridSize(2),2);
            v_indices_train_test = randperm(size(combin_pairs, 1), N_pair_train + N_pair_test);
            m_train_pairs = combin_pairs(v_indices_train_test(1:N_pair_train),:);
            m_test_pairs  = combin_pairs(v_indices_train_test(N_pair_train+(1:N_pair_test)),:);
                        
            %%
            disp 'Running only one experiment (no monte carlo yet)'
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                m_train_pairs, m_test_pairs, m_grid_x, m_grid_y);
            str_NMSE
            
%             M_trueMap = create_animation2([m_tracemap_pairs trueGains' ],  ...
%                 gridSize(1), gridSize(2), 'trueMap' );
%             M_locFree = create_animation2([m_tracemap_pairs locFreeMapEstimate' ],...
%                 gridSize(1), gridSize(2), 'locFree' );
%             M_locBased = create_animation2([m_tracemap_pairs locBasedMapEstimate' ],...
%                 gridSize(1), gridSize(2), 'locBased' );
% %             M_hybrid = create_animation2([m_tracemap_pairs hybridMapEstimate' ],...
% %                 gridSize(1), gridSize(2), 'hybrid' );             
%             M = create_animation3([m_tracemap_pairs trueGains' locFreeMapEstimate' ...
%                 locBasedMapEstimate' ], gridSize(1), gridSize(2), ...
%                 {'True', 'LocFree', 'LocBased'}, 's40x40') %hybridMapEstimate'  , 'hybrid'
%             NMSE_locFree
%             NMSE_locBased
% %             NMSE_hybrid
            F = GFigure.captureCurrentFigure;
        end

        function F = experiment_4020(obj, niter)
            % copied from exp 4010
            % evaluation of hybrid
            % method with 2-dimensional curve fitting.
            % No crossvalidation of the hybrid (only xvalidated values
            % obtained from the pure LocF and pure LocB
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=2:4; %!
            
            kernelSigmaLF = 69;
            kernelSigmaLB = 35;
            lambdaLF      = 1e-3;
            lambdaLB      = 1e-3;       
        
            gridSize = [5 5 ]; %! 20 20
            N_pair_train = 100;
            N_pair_test  = 100;
            
            mySim = Simulator2; 
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_syntheticLocError = 1;
            mySim.b_syntheticLocEstimate = 1;
            
            mySim.b_cv = 1;
            mySim.b_cvLambdas_hybrid = 1;
            % TODO: accelerate the quadprog that fits the 2-dimensional v_d
            % by removing the redundant constraints

            
            %% this section was adapted from 304 
            % (removal of CartographySimulations)
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            maxSamplesPerPilot = 10;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = 1/receiverBandwidth;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [mySim.generator, m_source_loc] = baselineGenerator2(...
                selectedWalls, x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            mySim.generator.delay_estimation_offset = 2*rand;
           
            %%
            % mySim.features = 1:size(combnk(1:numberOfSources,2),1);
            % we could move .baselineGenerator to the Simulator[2] class
            mySim.sampler = SpectrumMapSampler;
            mySim.sampler.pilotNoiseSTD = 1e-5;
            mySim.sampler.powerNoiseSTD = 0.5;
            mySim.sampler.maxSamplesPerPilot=maxSamplesPerPilot;
            %mySim.PCAenable=obj.PCAenabler;
            mySim.featureExtractor = FeatureExtractor;
            mySim.featureExtractor.sampling_period = mySim.generator.sampling_period;
                        
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
            mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            
            mySim.locBasedEstimator = LocationBasedEstimator;
            mySim.locEstimator = WangLocationEstimator;
            mySim.locBasedEstimator.kernel = @(x, y) ...
                exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
            mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            mySim.locBasedEstimator.locationNoiseSTD = []; %?
            mySim.locBasedEstimator.Xenb = m_source_loc;
            
            mySim.locEstimator.Xenb = m_source_loc;
            
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 20;
            mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
            mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
            mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
            mySim.hybridEstimator.regularizationParameterLB =lambdaLB;
            
            % The original CartographySimulations has several lines here to
            % define the number of different monte carlo runs and define the
            % random parameters or each MC run
            
            % define grid:
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(mySim.generator.x2,          mySim.generator.x1,           mySim.generator.n_gridpoints_x), ...
                linspace(mySim.generator.y_limits(1), mySim.generator.y_limits(2),  mySim.generator.n_gridpoints_y));
            
            %% train and evaluation(traceMap) pairs
            combin_pairs  = combnk(1:gridSize(1)*gridSize(2),2);
            v_indices_train_test = randperm(size(combin_pairs, 1), N_pair_train + N_pair_test);
            m_train_pairs = combin_pairs(v_indices_train_test(1:N_pair_train),:);
            m_test_pairs  = combin_pairs(v_indices_train_test(N_pair_train+(1:N_pair_test)),:);
                        
            %%
            disp 'Running only one experiment (no monte carlo yet)'
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                m_train_pairs, m_test_pairs, m_grid_x, m_grid_y);
            str_NMSE
            
%             M_trueMap = create_animation2([m_tracemap_pairs trueGains' ],  ...
%                 gridSize(1), gridSize(2), 'trueMap' );
%             M_locFree = create_animation2([m_tracemap_pairs locFreeMapEstimate' ],...
%                 gridSize(1), gridSize(2), 'locFree' );
%             M_locBased = create_animation2([m_tracemap_pairs locBasedMapEstimate' ],...
%                 gridSize(1), gridSize(2), 'locBased' );
% %             M_hybrid = create_animation2([m_tracemap_pairs hybridMapEstimate' ],...
% %                 gridSize(1), gridSize(2), 'hybrid' );             
%             M = create_animation3([m_tracemap_pairs trueGains' locFreeMapEstimate' ...
%                 locBasedMapEstimate' ], gridSize(1), gridSize(2), ...
%                 {'True', 'LocFree', 'LocBased'}, 's40x40') %hybridMapEstimate'  , 'hybrid'
%             NMSE_locFree
%             NMSE_locBased
% %             NMSE_hybrid
            F = GFigure.captureCurrentFigure;
        end

    end

end
