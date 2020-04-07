classdef RobustLocationEstimator
    properties
        param_rho  % upper bound assumed on the NLOS bias
        dim = 2    % dimension of search space, can also be 3
    end
    
    methods (Abstract)
        varargout = computeValue(obj, varargin)
        %computes the value for a given point in space v_x, 
        % given the vector of measured range differences v_d
        
        m_map = valueMap(obj, v_d, m_s, m_grid_x1, m_grid_x2, m_grid_x3)
        % Given a grid of points as produced by meshgrid, returns a matrix 
        % containing the objective value for each point in the grid.
        % The grid may be 2D or 3D.
    end
end