function [generator, source_loc] = baselineGenerator2(selectedWalls, ...
    x_wall_file_name, y_wall_file_name, gridSize, generator_template)

assert(isa(generator_template, 'MultiWallChannelGainGenerator'));
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

generator.f=generator_template.f;           % Frequency in Hz

xrange=linspace(x1,x2,gridSize(1));
yrange=linspace(y1,y2,gridSize(2));
generator.xt=  [xrange(4), yrange(4)]; %Luismi: What is the meaning of this line? Is it a random point close to the edge?
generator.xt_loc=  [0, 28, 32 ,-3, 42]; %
generator.yt_loc=  [0, 10, -7, 30, 30]; %
generator.ptx=-30*ones(1, size(generator.xt,1));   % Transmitter powers in dBW
generator.ptx_loc=-30*ones(1, length(generator.xt_loc));   % Loc anchor nodes in dBW

generator.zt=1.5; %transmitter Height
generator.x1=x1; % First boundary  along x
generator.x2=x2; % Second boundary  along x
generator.path_loss_exp = 2;
generator.n_gridpoints_x=gridSize(1);
generator.n_gridpoints_y=gridSize(2);
generator.zr=1.5; %receiverHeight;

%generator.sampling_period=1/cartSimObj.receiverBandwidth;
generator.sampling_period = generator_template.sampling_period;
generator.maxSamplesPerPilot=generator_template.maxSamplesPerPilot;
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


