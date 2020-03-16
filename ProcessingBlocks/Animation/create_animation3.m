function [M] = create_animation3(m_data, ngrid_x, ngrid_y, s_titles, s_filename)
% Given:
% - A matrix containing:
%    - First  column: index of tx position
%    - Second column: index of rx position
%    - Third and next columnts: value of the channel gains
% Output:
% - A Matlab movie showing at each instant a channel gain map between the 
%   transmitter (located at the corresponding trajectory point) and all
%   receivers (each point in the map would correspond to a receiver)


Ngrid = ngrid_x*ngrid_y;
Ntraj = ceil(size(m_data, 1)./Ngrid);
assert(size(m_data, 1)==Ngrid*Ntraj);
nPanes = size(m_data,2)-2;

if nargin<=4
    s_filename = 'Cartography';
end
if nargin<=3   
    warning 'No title strings were provided'
    s_titles = cell(nPanes,1);
    for p = 1:nPanes
        s_titles{p} = sprintf('map%d', 1:nPanes)
    end
end

cf = figure(1); clf; axis off
set(cf, 'WindowStyle', 'normal');
set(cf, 'Position', [1 1 1280 720]);
M(Ntraj) = getframe();

map_now = zeros(ngrid_x,ngrid_y);
for k = 1:Ntraj
    rows = (k-1)*Ngrid+(1:Ngrid);
    locationIndex_now = m_data(rows(1),1);
    assert( all(m_data(rows,1)==locationIndex_now) );
    [x_now, y_now]  = ind2sub([ngrid_x,ngrid_y], locationIndex_now);
    clf;
    for p = 1:nPanes
        subplot(1, nPanes, p);
        map_now(m_data(rows,2)) = m_data(rows,2+p);
        % Draw map
        cla;
        imagesc(map_now); 
        axis square; hold on
        plot(y_now, x_now, 'xk')
        title(s_titles{p});
    end
    drawnow
    M(k) = getframe(cf);    
end
vidObj = VideoWriter([s_filename '_movie.avi']);
vidObj.open
vidObj.writeVideo(M);
vidObj.close