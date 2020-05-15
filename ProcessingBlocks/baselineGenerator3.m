function [generator, source_loc] = baselineGenerator3(...
    environment, selectedWalls, generator_template)

assert(isa(generator_template, 'MultiWallChannelGainGenerator'));
generator = MultiWallChannelGainGenerator;

generator.boundary= environment.m_xyz_boundary;
% [x_min x_max; y_min y_max; z_min z_max]

generator.xt_loc= environment.m_xy_locSources(1,:); % [0, 28, 32 ,-3, 42];
generator.yt_loc= environment.m_xy_locSources(2,:); % [0, 10, -7, 30, 30];

generator.ptx_loc = -30*ones(1, length(generator.xt_loc));
% Loc anchor nodes transmit power (in dBW)

generator.path_loss_exp = 2; % path loss exponent

generator.z_rx_default=1.5;  % receiverHeight;

generator.f                  = generator_template.f; % Carrier Frequency (in Hz)
generator.sampling_period    = generator_template.sampling_period;
generator.maxSamplesPerPilot = generator_template.maxSamplesPerPilot;

if isempty(selectedWalls) || not(any(selectedWalls))
    generator.estimateLocFreeSpace=1;
else
    generator.estimateLocFreeSpace=0;
end

arrangedCoord=[2*selectedWalls-1;2*selectedWalls];

generator.X = environment.m_xy_walls(1, arrangedCoord(:))';
generator.Y = environment.m_xy_walls(2, arrangedCoord(:))';

source_loc = [generator.xt_loc; generator.yt_loc];
assert (isequal(source_loc, environment.m_xy_locSources));


