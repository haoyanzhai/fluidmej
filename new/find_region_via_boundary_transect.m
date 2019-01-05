function region_no = find_region_via_boundary_transect(boundaries, boundary_pt_x, boundary_pt_y)

% Given a point on one boundary transect, find the index of the neighboring 
% flow regions by checking if [boundary_pt_x, boundary_pt_y] belongs to a specific boundary.

% Example: 
% region_no = find_region(regions, boundary_pt_x, boundary_pt_y);

% Input: 

% regions: structure containing flow region No., Lon/Lat, vx/vy. Output of
% flow_boundaries.m
% boundary_pt_x: Lon
% boundary_pt_y: Lat

% Output: 
% region_no: Row vector. Region index number of the neighboring grid points of the boundary.
% [n_1, n_2, ..., n_N], (n_1>n_2>...>n_N)

% Mengxue Hou, Georgia Institute of Technology
% 01/2019

field_names = fieldnames(boundaries);
boundary_name = [];
region_no = [];
for ii = 1:numel(field_names)
    current_boundary = getfield(boundaries, field_names{ii});
    for jj = 1:size(current_boundary,1)-1
        p_point = current_boundary(jj,:);
        n_point = current_boundary(jj+1,:);
        if (n_point(1) - boundary_pt_x)/(n_point(2) - boundary_pt_y) ==  (p_point(1) - boundary_pt_x)/(p_point(2) - boundary_pt_y)
            %if isempty(boundary_name)
                boundary_name = field_names(ii);
            %else
            %    error(['The point ', num2str(boundary_pt_x), ', ', num2str(boundary_pt_y), ' is on multiple boundaries.']);
            %end
            region_no = [region_no, str2num(boundary_name{1}(end)), str2num(boundary_name{1}(end-1))]; 
        end
    end
end
region_no = sort(unique(region_no), 'ascend');

