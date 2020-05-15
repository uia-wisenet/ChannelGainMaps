classdef TwoDimensionalWeightCalculatingUnit < WeightCalculatingUnit
      
    properties
         m_locErrors_train
         m_indicesAllOrder
         m_indices_DAG
         v_weights
    end
    methods
        function v_weights_out = evalWeights(obj, m_locErrors_query)
            % Given a matrix with pairs of uncertainty measures (for tx and
            % rx), evaluates the gating function
            assert(size(m_locErrors_query,2 )==2);
            n_ues = size(m_locErrors_query, 1);
            v_weights_out = zeros(n_ues, 1);
            
            for ii = 1:n_ues
                v_e = m_locErrors_query(ii,:);
                %figure(901); clf                
                mb_upper = (obj.m_locErrors_train >= v_e);
                if not(any(all(mb_upper,2)))
                    v_weights_out(ii) = min(obj.v_weights);
                    continue
                end
                v_ind_upper  = obj.bound_triage_DAG(mb_upper, ...
                    obj.m_indices_DAG, v_e);
%                 v_ind_upper2 = obj.bound_triage(mb_upper, ...
%                     obj.m_indicesAllOrder, v_e);
%                 if not(isempty(setxor(v_ind_upper, v_ind_upper2)))
%                     keyboard
%                 end
                
                %figure(901); clf
                mb_lower = (obj.m_locErrors_train <= v_e);
                if not(any(all(mb_lower,2)))
                    v_weights_out(ii) = max(obj.v_weights);
                    continue
                end
                v_ind_lower  = obj.bound_triage_DAG(mb_lower, ...
                    fliplr(obj.m_indices_DAG), v_e);
%                 v_ind_lower2 = obj.bound_triage(mb_lower, ...
%                     fliplr(obj.m_indicesAllOrder), v_e);
%                 if not(isempty(setxor(v_ind_lower, v_ind_lower2)))
%                     keyboard
%                 end
                
                v_dists = vecnorm(obj.m_locErrors_train([v_ind_upper; ...
                    v_ind_lower],:) - v_e, 2, 2);
                if any(v_dists==0)
                    v_dists = v_dists + 1e-5;
                end
                v_coefs = 1./v_dists;
                v_weights_out(ii) = v_coefs'*obj.v_weights([v_ind_upper; ...
                    v_ind_lower])/sum(v_coefs);
            end
                    
        end
        
        function [v_indices_out] = ...
                bound_triage_DAG(obj, mb_bound, m_indices_DAG, v_e)
            %coge el subgrafo formado por los ancestros de v_e:
%             n_edges = size(m_indices_DAG,1);
%             b_isInSubgraph = false(n_edges, 1);
            vb_boundIndices = all(mb_bound, 2);
%             for i_edge = 1:n_edges
%                 b_isInSubgraph(i_edge) = vb_boundIndices(m_indices_DAG(i_edge,1)) ...
%                     && vb_boundIndices(m_indices_DAG(i_edge,2));
%             end
            b_isInSubgraph = vb_boundIndices(m_indices_DAG(:,1)) ...
                & vb_boundIndices(m_indices_DAG(:,2));
%             b_isCut = ~vb_boundIndices(m_indices_DAG(:,2)) &...
%                     vb_boundIndices(m_indices_DAG(:,1));
            % identifica los nodos que no tienen hijos
            v_indices_out = setdiff(m_indices_DAG(b_isInSubgraph,2), ...
                m_indices_DAG(b_isInSubgraph,1));
%             v_indices_out = unique(m_indices_DAG(b_isCut,1));
            if isempty(v_indices_out)
                top = setdiff(m_indices_DAG(:,1), ...
                    m_indices_DAG(:,2));
                v_indices_out = top;
            end
        end
            
           
            
            
%             % start with the leaves
%             n_train = size(mb_bound,1);
%             vb_candidates = false(n_train,1);
%             vb_candidates(setdiff(1:n_train, m_indices_DAG(:,1))) = 1;
%             vb_candidates_prev = nan;
%             while not(isequal(vb_candidates, vb_candidates_prev))
%                 vb_candidates_prev = vb_candidates
%                 for jj = 1:n_train
%                     if vb_candidates(jj)
%                         if all(
%                         % look at the children of the node
%                         v_children = m_indices_DAG(m_indices_DAG(:,2)==jj,1)
%                         
%             end
%             
%         end
        
        function [v_indices_out] = ...
                bound_triage(obj, mb_bound, m_indices, v_e)
            % Given a set of points, identifies the minimum ones
            % (points such that there is no other points in the set that
            % are lower bounds)
            
            % the input variable v_e is only for plotting purposes
            
            % if the property m_indices_DAG is defined, then this triage 
            % becomes much easier (TODO)
            
            % This function identifies the immediate ancestors of a given
            % point
            
%              plot(obj.m_locErrors_train(:,1), ...
%                  obj.m_locErrors_train(:,2),...
%                  'b. '); hold on
%              plot(v_e(1), v_e(2), 'g+');
            vb_boundIndices = all(mb_bound, 2);
            n_train = length(vb_boundIndices);
            prev_sum = n_train;
            while sum(vb_boundIndices) < prev_sum
                prev_sum = sum(vb_boundIndices);
                for jj = 1:n_train
                    if vb_boundIndices(jj)
%                          plot(obj.m_locErrors_train(jj,1),...
%                              obj.m_locErrors_train(jj,2), 'ko');
                        % select all upper bounds of jj:
                        indices_toDiscard = m_indices( ...
                            m_indices(:,2)==jj, 1  );
%                          plot(obj.m_locErrors_train(indices_toDiscard,1), ...
%                              obj.m_locErrors_train(indices_toDiscard,2), ...
%                              'rx ')
                        vb_boundIndices(indices_toDiscard) = 0;
                    end
                end
            end
            v_indices_out = find(vb_boundIndices);
        end
    end
end

