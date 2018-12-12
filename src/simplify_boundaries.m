function sim_boundaries = simplify_boundaries(boundaries)

% Simplify the boundaries of flow regions by keeping only the terminal
% points of straight-line section of one boundary, and deleting all other
% points. 

% Input: boundaries. 
% Output: the simplified boundaries.

% Mengxue Hou
% Georgia Institute of Technology, 2018/12

boundary_names = fieldnames(boundaries);

for ii = 1:numel(boundary_names)
    % Iterate over all points on a boundary. Compare the slope of current
    % transect with previous transect to determine if the current point is
    % on a straight-line transect.     
    current_boundary_name = boundary_names(ii);
    current_boundary = getfield(boundaries, current_boundary_name{1});    
    slope_p = slope(current_boundary(1,:)', current_boundary(2,:)');
    
    for jj = 2:size(current_boundary,1)-1
        slope_n = slope(current_boundary(jj,:)', current_boundary(jj+1,:)');
        if slope_n == slope_p            
           current_boundary(jj,:) = [NaN, NaN];
        else
            slope_p = slope_n;
        end
    end
    
    sim_current_boundary_idx = find(~isnan(current_boundary(:,1)));
    sim_boundaries.(current_boundary_name{1}) = current_boundary(sim_current_boundary_idx,:);
end
end

