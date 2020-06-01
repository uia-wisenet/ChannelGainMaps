clear
gsinfo;
load savedResults/results_2510.mat
v_nAllTrains = v_nTrains;
tb_allNMSE   = tb_NMSE;
s_fields = string(tb_NMSE.Properties.VariableNames);
load savedResults/results_2520.mat
for nTrain = v_nAllTrains
    if ismember(nTrain, v_nTrains)
        index_d = find(v_nAllTrains==nTrain);
        index_s = find(v_nTrains==nTrain)
%         for i_f = 1:length(s_fields)
%             tb_allNMSE.(s_fields(i_f))(index_d)...
%                 = tb_NMSE.(s_fields(i_f))(index_s);
%         end
        tb_allNMSE.locBased_fromHybrid(index_d) = ...
            tb_NMSE.locBased_fromHybrid_corrected(index_s);
    end
end

tb_allNMSE.locBased_fromHybrid(find(v_nAllTrains==2000)) = 0.6184;

F = GFigure;
F.m_X = v_nAllTrains;
F.m_Y = tb_allNMSE.Variables';
F.ch_interpreter = 'none';
F.c_legend = {'locF (step 1)', 'locB (step 2)', 'MoE (step 8)', ...
'MoE.locF ($\mathbf{f}_p$)', 'MoE.locB ($\mathbf{f}_l$)', 'MoE.locB CORRECTED'};
F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
F.ch_xlabel = 'Number of training samples';
F.ch_ylabel = 'NMSE';
F.ch_legendPosition = 'North';
F.ch_interpreter = 'Latex';
F.figureNumber = str2num(ch_expNum);

F.plot
ff = gcf;
ff.Children(2).YScale = 'log';
ylim([0.25, 0.7]);