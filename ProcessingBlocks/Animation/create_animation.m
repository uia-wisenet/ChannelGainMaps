function [M] = create_animation(cgBlock, trajectory)
% Given:
% - A trajectory
% - A channel gain map generator OR a channel gain map estimator
% Output:
% - A Matlab movie showing at each instant a channel gain map between the 
%   transmitter (located at the corresponding trajectory point) and all
%   receivers (each point in the map would correspond to a receiver)

assert(size(trajectory, 2) == 2, 'trajectory must be an N-by-2 matrix');
validClasses = {'MapGenerator', 'LocationFreeEstimator', ...
    'LocationBasedEstimator'};
cgBlockClass = class(cgBlock);
assert(any(strcmp(cgBlockClass, validClasses)), ...
    'Invalid class for channel generating block')

N = size(trajectory, 1); %number of points in trajectory
cf = figure(1); clf; axis off
set(cf, 'WindowStyle', 'normal');
set(cf, 'Position', [1 1 1280 720]);
M(N) = getframe();

for k = 1:N
    point_now = trajectory(k,:);
    switch cgBlockClass
        case 'MapGenerator'
            [~, map_now] = cgBlock.generateCGMap(point_now);
%         case 'LocationFreeEstimator'
%         case 'LocationBasedEstimator'
%             
        otherwise
            error ('Invalid class OR not implemented yet')
    end
    % Draw map
    
    drawnow
    M(k) = getframe();    
end
s_filename = 'Cartography';
vidObj = VideoWriter([s_filename '_movie.avi']);
vidObj.open
vidObj.writeVideo(M);
vidObj.close