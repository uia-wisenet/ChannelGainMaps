classdef Animator
    properties
        v_indicesTxPosition(:, 1)
        t_grid_xy(:, :, 2)
        generator MultiWallChannelGainGenerator
        v_colorAxis(2, 1) = nan(2,1)
        v_figurePosition(1, 4) = [1 1 1280 720]
        c_plotTransmitter_options = {'^k', 'MarkerFaceColor', 'w'};
        ch_colormap = 'copper'
        c_plotWall_options = {'Color', [0, 0.3, 0.1], 'LineWidth', 4}
        c_plotSource_options = {'pb', 'MarkerFaceColor', 'b', 'markersize', 4};
    end 
    
    methods
        
        function playMovie(obj, M)
            cf = figure(gcf); % to make the figure visible;
            v_currentPosition = get(cf, 'Position');
            v_moviePosition = v_currentPosition;
            v_moviePosition(3:4) = obj.v_figurePosition(3:4) - obj.v_figurePosition(1:2);
            set(cf, 'Position', v_moviePosition)
            movie(cf, M, 1, 1.5, [1 1 0 0]);
        end
        
        function n_out = n_frames(obj)
            n_out = length(obj.v_indicesTxPosition);
        end
        
        function n_out = n_gridPoints(obj)
            my_size = size(obj.t_grid_xy);
            n_out = prod(my_size(1:2));
        end
        
        function m_out = m_pairs(obj)
            m_out = [ kron(obj.v_indicesTxPosition, ...
                ones(obj.n_gridPoints, 1)), ...
                repmat((1:obj.n_gridPoints)', [obj.n_frames 1]) ];
        end
        
        function frame = draw_frame(obj, t_gains, v_xy_tx, h_f, s_titles)
            v_grid_x = obj.t_grid_xy(:,1,1);
            v_grid_y = obj.t_grid_xy(1,:,2);
            subplot_size = [size(t_gains,3), size(t_gains,4)];
            figure(h_f); clf
            for i_row = 1:subplot_size(1)
                for i_col = 1:subplot_size(2)
                    subplot(subplot_size(1), subplot_size(2), ...
                        sub2ind(fliplr(subplot_size), i_col, i_row))
                    axis off;
                    imagesc(v_grid_x, v_grid_y, t_gains(:,:,i_row, i_col)');
                    hold on
                    xlabel x; ylabel y
                    if not(isempty(obj.generator))
                        obj.generator.plot_walls(obj.c_plotWall_options{:});
                        obj.generator.plot_sources(obj.c_plotSource_options{:});
                        obj.generator.set_lims_to_plot();
                    end
                    plot(v_xy_tx(1), v_xy_tx(2), obj.c_plotTransmitter_options{:});
                    try
                        caxis(obj.v_colorAxis);
                    catch ME
                        warning('Could not set color limits')
                        disp( getReport( ME))
                    end
                    title(s_titles(i_row, i_col), 'FontSize', 14);
                end
            end
            colormap(obj.ch_colormap);
            frame = getframe(h_f);
        end
        
        function M = storyBoard(obj, m_gains_in, s_titles, v_key_frames)
            assert(isvector(s_titles));
            n_cols = length(s_titles);
            assert(isvector(v_key_frames));
            assert(max(v_key_frames)<= obj.n_frames);
            n_rows = length(v_key_frames);
            assert(n_cols==size(m_gains_in, 2))
            figure(gcf); clf;
            v_grid_x = obj.t_grid_xy(:,1,1);
            v_grid_y = obj.t_grid_xy(1,:,2);
            t_gains = reshape(m_gains_in, ...
                [length(v_grid_x), length(v_grid_y), obj.n_frames, n_cols]);
            if isnan(obj.v_colorAxis(1))
                obj.v_colorAxis(1) = min(vec(t_gains(:,:,v_key_frames, :)));
            end
            if isnan(obj.v_colorAxis(2))
                obj.v_colorAxis(2) = max(vec(t_gains(:,:,v_key_frames, :)));
            end
            for i_col = 1:n_cols
                for i_row = 1:n_rows
                    i_frame = v_key_frames(i_row);
                    i_txPosition = obj.v_indicesTxPosition(i_frame);
                    v_xy_tx = [obj.t_grid_xy(i_txPosition), ...
                        obj.t_grid_xy(i_txPosition+obj.n_gridPoints())];
                    subplot(n_rows, n_cols, ...
                        sub2ind([n_cols n_rows], i_col, i_row))
                    imagesc(v_grid_x, v_grid_y, ...
                        t_gains(:,:,v_key_frames(i_row), i_col)');
                    axis equal; axis tight
                    xlabel x, ylabel y
                    hold on
                    if not(isempty(obj.generator))
                        obj.generator.plot_walls(obj.c_plotWall_options{:});
                        obj.generator.plot_sources(obj.c_plotSource_options{:});
                        obj.generator.set_lims_to_plot();
                    end
                    plot(v_xy_tx(1), v_xy_tx(2), obj.c_plotTransmitter_options{:});
                    try
                        caxis(obj.v_colorAxis);
                    catch ME
                        warning('Could not set color limits')
                        disp( getReport( ME))
                    end
                    if i_row ==1
                        title(s_titles(i_col), 'FontSize', 14)
                    end
                end
            end
            colormap(obj.ch_colormap);
        end
        
        function M = createGivenMatrix(obj, m_gains_in, s_titles)
            assert(ismatrix(m_gains_in));
            assert(numel(s_titles)==size(m_gains_in,2));
            M = obj.create(reshape(m_gains_in, [size(m_gains_in,1), ...
                size(s_titles,1), size(s_titles,2)]), s_titles);
        end
        
        function M = create(obj, t_gains_in, s_titles)
            size_gains = size(t_gains_in);
            if(ismatrix(t_gains_in))
                size_gains(3) = 1;
            end
            size_subplot = size_gains(2:end);
            if exist('s_titles', 'var')
                assert(isequal(size(s_titles), size_subplot), ...
                    'size of s_titles must match dimensions 2 and 3 of m_gains');
            else
                s_titles = repmat("", size_subplot);
            end
            assert(size(t_gains_in, 1) == obj.n_gridPoints*obj.n_frames);
            
            if isnan(obj.v_colorAxis(1))
                obj.v_colorAxis(1) = min(t_gains_in(:));
            end
            if isnan(obj.v_colorAxis(2))
                obj.v_colorAxis(2) = max(t_gains_in(:));
            end

            m_grid_x = obj.t_grid_xy(:,:,1);
            m_grid_y = obj.t_grid_xy(:,:,2);
            t_gains = reshape(t_gains_in, ...
                [size(m_grid_x), obj.n_frames, size_subplot]);
            cf = gcf;
            clf, subplot(size_subplot(1), size_subplot(2), 1);
            %axis off
            set(cf, 'WindowStyle', 'normal');
            set(cf, 'Position', obj.v_figurePosition);
            M(obj.n_frames) = getframe(cf); %initializing vector of frames
            for i_frame = 1:obj.n_frames
                i_txPosition = obj.v_indicesTxPosition(i_frame);
                v_xy_tx = [m_grid_x(i_txPosition), m_grid_y(i_txPosition)];
                S.type = '()';
                S.subs= repmat({':'}, 1, 2+ndims(t_gains_in));
                S.subs{3} = i_frame;
                M(i_frame) = obj.draw_frame(squeeze(subsref(t_gains, S)), ...
                    v_xy_tx, cf, s_titles);
            end
        end
    end
    
    methods (Static)       
        function saveMovie(M, s_fileName_in)
            % Creates an .avi file and saves the frames gathered in the
            % structure vector M as a video.
            % Input:
            % M: each entry of M is a struct as generated by getframe();
            % 
            s_fileName = s_fileName_in;
            labelNum = 0;
            while exist(s_fileName, 'file')
                labelNum = labelNum + 1;
                s_fileName = s_fileName_in+string(labelNum);
            end
            if labelNum > 0
                warning(s_fileName_in + " already exists. Saving as " +...
                    s_fileName);
            end
            vidObj = VideoWriter([s_fileName '.avi']);
            vidObj.open
            vidObj.writeVideo(M);
            vidObj.close
        end
    end
end