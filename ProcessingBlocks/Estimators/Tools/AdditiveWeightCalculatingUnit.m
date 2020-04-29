classdef AdditiveWeightCalculatingUnit < WeightCalculatingUnit
    properties
        v_locErrors
        v_weights
    end
    
    methods
        function v_w_out = evalWeights(obj, m_locErrors)
            assert(size(m_locErrors,2 )==2);
            v_locErrors_t = m_locErrors(:,1);
            v_locErrors_r = m_locErrors(:,2);
            v_wt = interp1(obj.v_locErrors, obj.v_weights, ...
                v_locErrors_t, 'linear', 'extrap');
            v_wr = interp1(obj.v_locErrors, obj.v_weights, ...
                v_locErrors_r, 'linear', 'extrap');
            v_w_out = (v_wt + v_wr)/2;
        end
                        
    end
end