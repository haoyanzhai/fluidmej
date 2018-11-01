function region_no = find_region(regions, boundary_pt_x, boundary_pt_y)

% Given a point on one boundary transect, find the number of neighboring flow regions.

% Example: 
% region_no = find_region(regions, boundary_pt_x, boundary_pt_y);

% Input: 

% regions: structure containing flow region No., Lon/Lat, vx/vy. Output of
% flow_boundaries.m
% boundary_pt_x: Lon
% boundary_pt_y: Lat

% Output: 
% region_no: Row vector. Region number of the neighboring grid points of the boundary.
% [n_1, n_2, ..., n_N], (n_1>n_2>...>n_N)

% Mengxue Hou, Georgia Institute of Technology
% 10/2018

if ceil(boundary_pt_x) == floor(boundary_pt_x)
    x_1 = boundary_pt_x - 1;
    x_2 = boundary_pt_x + 1;
else
    x_1 = floor(boundary_pt_x);
    x_2 = ceil(boundary_pt_x);
end

if ceil(boundary_pt_y) == floor(boundary_pt_y)
    y_1 = boundary_pt_y - 1;
    y_2 = boundary_pt_y + 1;
else
    y_1 = floor(boundary_pt_y);
    y_2 = ceil(boundary_pt_y);
end

neighbor_pts = [x_1, y_1; x_1, y_2; x_2, y_1; x_2, y_2];

region_no = [];
for ii = 1:size(neighbor_pts,1)
    for jj = 1:numel(regions.lon)
        if ~isempty(intersect(find(regions.lon{jj} == neighbor_pts(ii,1)), find(regions.lat{jj} == neighbor_pts(ii,2))))
            region_no = [region_no, jj];
        end
    end
end

region_no = sort(unique(region_no), 'descend');