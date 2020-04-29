classdef ForwardHybridKernelMachine
% An object of this class just contains the results of learning from a
% single training set, and provides a method to evaluate the estimated
% function for a query point.

properties        
    fkm_lf ForwardLFKernelMachine
    fkm_lb ForwardLBKernelMachine
    wcu    WeightCalculatingUnit
end

methods
    function v_estimate_hybrid = evaluate(obj, ...
            m_queryFeatures_LF, m_queryFeatures_LB, m_locationErrors)
        
        %input check
        assert(ismatrix(m_locationErrors));
        n_ues = size(m_locationErrors,1);
        assert(size(m_locationErrors, 2) == 2);
        %assert(ndims(t_queryFeatures_LF)==3);
        %assert(ndims(t_queryFeatures_LB)==3);
        assert(size(m_queryFeatures_LB, 2) == n_ues);
        assert(size(m_queryFeatures_LF, 2) == n_ues);
        assert(size(m_queryFeatures_LB, 1) == ...
            size(obj.fkm_lb.m_trainFeatures, 1));
        assert(size(m_queryFeatures_LF, 1) == ...
            size(obj.fkm_lf.m_trainFeatures, 1));
        assert(size(m_queryFeatures_LB, 3) ==2);
        assert(size(m_queryFeatures_LF, 3) ==2);
        
        % obtain LF estimate:
        v_estimate_LF = obj.fkm_lf.evaluate(m_queryFeatures_LF);
        % obtain LB estimate:
        v_estimate_LB = obj.fkm_lb.evaluate(m_queryFeatures_LB);
        % obtain eval weights from interpolating the input errors:
        v_weights = obj.wcu.evalWeights(m_locationErrors);
        v_estimate_hybrid = v_weights  .* v_estimate_LB(:) + ...
                         (1-v_weights) .* v_estimate_LF(:);
       
    end
end
end