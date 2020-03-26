function v_out = exagg_linspace(x_0,  x_1, ...
                min_z, max_z, npoints)
            % Linspace with exaggeration factors. 
            % The values min_z and max_z define a linear transformation,
            % such that if [min_z max_z] = [0 1], this does the same as 
            % linspace. 
            v_z = linspace(min_z, max_z, npoints);
            v_out = x_0 + (x_1-x_0)*v_z;
%             v_out2 = x_0*(1-v_z) + x_1*v_z;
%             v_out3 = x_0 + (x_1-x_0)*linspace(min_z, max_z, npoints);
        end