function [M] = create_animation2(m_data, ngrid_x, ngrid_y, s_filename)
% Given:
% - A trajectory
% - A channel gain map generator OR a channel gain map estimator
% Output:
% - A Matlab movie showing at each instant a channel gain map between the 
%   transmitter (located at the corresponding trajectory point) and all
%   receivers (each point in the map would correspond to a receiver)

% assert(size(trajectory, 2) == 2, 'trajectory must be an N-by-2 matrix');
% validClasses = {'MapGenerator', 'LocationFreeEstimator', ...
%     'LocationBasedEstimator'};
% cgBlockClass = class(cgBlock);
% assert(any(strcmp(cgBlockClass, validClasses)), ...
%     'Invalid class for channel generating block')

%N = size(trajectory, 1); %number of points in trajectory
Ngrid = ngrid_x*ngrid_y;
Ntraj = ceil(size(m_data, 1)./Ngrid);
assert(size(m_data, 1)==Ngrid*Ntraj);
assert(size(m_data,2)==3);
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
    map_now(m_data(rows,2)) = m_data(rows,3);
    % Draw map
    clf;
    imagesc(map_now); hold on
    plot(y_now, x_now, 'xk')
    drawnow
    M(k) = getframe();    
end
if nargin==3
    s_filename = 'Cartography';
end
vidObj = VideoWriter([s_filename '_movie.avi']);
vidObj.open
vidObj.writeVideo(M);
vidObj.close