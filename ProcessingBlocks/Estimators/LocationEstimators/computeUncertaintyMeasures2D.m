function [v_diameter, v_area] = computeUncertaintyMeasures2D(...
    m_grid_x1, m_grid_x2, m_values, levels)
%COMPUTEUNCERTAINTYMEASURES2D computes the diameter and the area measure of
%the level set of the given function

if numel(levels) == 1
    v_levels = levels*ones(1, 2);
else
    v_levels = levels;
end

m_c = contour(m_grid_x1, m_grid_x2, m_values, v_levels);
index = 1;
v_contourLevels = [];
for k = 1:size(m_c,2)
    npoints_now = m_c(2, index);
    c_m_points{k} = m_c(:,index+(1:npoints_now)); %#ok<AGROW>
    v_contourLevels(k) = m_c(1, index);
    index = index + 1 + npoints_now;
    if index == size(m_c, 2)+1
        break;
    elseif index > size(m_c,2)
        error ('failed contour matrix structure')
    end
end

if length(c_m_points) > length(levels)
    error ('sorry, some level subsets seem to be disjoint')
end
m_gridArea = conv2(diff(m_grid_x1, 1, 2), 1/2*[1 1]).* ...
    conv2(diff(m_grid_x2, 1, 1), 1/2*[1;1]);
for k = 1:length(levels)
    m_mask = m_values<v_levels(k);
    if any(vec(m_mask([1 end], :))) || any(vec(m_mask(:, [1 end])))
        warning 'level set touches the boundary of the grid'
    end

    m_pdist = pdist(c_m_points{k}');
    v_diameter(k) = max(m_pdist(:));
    
    m_contrToArea = m_gridArea.*m_mask;
    v_area(k) = sum(m_contrToArea(:));
end
