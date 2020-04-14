classdef PrecomputedMapGenerator < MapGenerator
	% PrecomputedMapGenerator Load power map from file
	properties
		ch_filename	
	end
	
	methods
		function [map_out, evaluationGrid_x, evaluationGrid_y, ...
                source_loc, estimatedDistances] = generateMap(obj)
			% TODO: write specification for the file containing the map
			% information
			S = load(obj.ch_filename);
            map_out=pow2db(S.power_with_shadowing_ns(:,:,3));
            evaluationGrid_x = S.evaluationGrid_x;
            evaluationGrid_y = S.evaluationGrid_y;
%             R = S.R;
            source_loc = S.source_loc;
			
			estimatedDistances = reshape( pdist2(...
				[evaluationGrid_x(:),evaluationGrid_y(:)], source_loc'),...
                [size(evaluationGrid_x), size(source_loc,2)]);
		end
	end
end
