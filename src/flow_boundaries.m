function [Boundaries, Intersection_pts, Regions] = flow_boundaries(field, starting_pos)
% From the flow field data, first divide the flow field into regions
% according to the flow speed information at each position. Then find the
% boundaries of the regions, and intersection points of the boundaries.

% Input:
% field.lon/ field.lat: Lon/Lat grid vectors;
% field.vcx/ field.vcy: horizontal/vertical flow data at each position.
% starting_pos: coloum vector of vehicle starting position.

% Output:
% Boundaries:
% Boundaries.boudary_no_n_1 n_2: boundary between region n_1 and region n_2
% (n_1 > n_2).
% Boundaries.boudary_no_n_1 n_2(:,1): Lon.
% Boundaries.boudary_no_n_1 n_2(:,2): Lat.

% Intersection_pts:
% Intersection_pts.reg_no(n,:): region number for the boundaries of the nth
% Intersection point.
% Intersection_pts.pos(n,:): position of the nth intersection point ([Lon, Lat]).  

% Regions.no: total number of regions.
% Regions.lon/lat: Lon/Lat for grid points in each of the regions.
% Regions.vx/vy: Flow speed for grid points in each of the regions.


% Mengxue Hou
% Georgia Institute of Technology, 2018/10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% Divide the flow field into regions using K-means

nx = numel(field.lon);
ny = numel(field.lat);
vx = field.vcx;
vy = field.vcy;
lon = field.lon;
lat = field.lat;
if (range(diff(lon)) == 0) && (range(diff(lat)) == 0)
    delta_lon = diff(lon(1:2));
    delta_lat = diff(lat(1:2));
else
    error('The flow map grid is not uniform!');
end

% Flow map input
figure
[xx,yy] = meshgrid(lon, lat);
quiver(xx,yy,vx,vy);
axis([min(lon), max(lon), min(lat), max(lat)]);
print(gcf, [pwd '/Generated_Plots/Flow_field'], '-r300', '-dpng')
close gcf

% Vectorize the vcx, vcy, lon/lat data
vcx = reshape(vx, [nx*ny, 1]);
vcy = reshape(vy, [nx*ny, 1]);

[xx, yy] = meshgrid(lon,lat);
gridx = reshape(xx, [nx*ny, 1]);
gridy = reshape(yy, [nx*ny, 1]);

X = [vcx, vcy, gridx/(2*nx), gridy/(2*ny)];

k = 7;
opts = statset('Display','final');
[idx,C] = kmeans(X,k,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

for ii = 1:k
    distance_to_orig(ii) = norm([C(ii,3)*2*nx, C(ii,4)*2*ny] - starting_pos');
end

[~, c_idx] = sort(distance_to_orig);

figure
Marker = jet(k);
for n_idx = 1:k
    ii = c_idx(n_idx);
    Regions.no = n_idx;
    Regions.lon{n_idx} = round(X(idx==ii,3)*2*nx);
    Regions.lat{n_idx} = round(X(idx==ii,4)*2*ny);
    Regions.vx{n_idx} = X(idx==ii,1);
    Regions.vy{n_idx} = X(idx==ii,2);
    
    plot(Regions.lon{n_idx},Regions.lat{n_idx}, 'color', Marker(n_idx,:), 'Marker', '.', 'LineStyle', 'none');hold on;
    plot(C(n_idx,3)*2*nx, C(n_idx,4)*2*ny, 'ko');hold on;
end
axis([min(lon), max(lon), min(lat), max(lat)]);

print(gcf, [pwd '/Generated_Plots/Divided_flow_field'], '-r300', '-dpng');
close gcf
% end

%% Check whether each of the divided regions are connected
Regions = check_regions(Regions, delta_lon, delta_lat);

%% Find boundaries of the flow regions
% The boundaries are the center line between two neighboring regions
[Boundaries, Intersection_pts] = boundary_points(Regions, field, delta_lon, delta_lat);

% Include intersection_pts into boundaries
for int_idx = 1:size(Intersection_pts.reg_no,1)
    reg = Intersection_pts.reg_no(int_idx, :);
    reg = sort(reg);
    b_name = ['boundary_no_', num2str(reg(3)), num2str(reg(2))];
    Boundaries.(b_name)(end+1,:) = Intersection_pts.pos(int_idx,:);
    b_name = ['boundary_no_', num2str(reg(3)), num2str(reg(1))];
    Boundaries.(b_name)(end+1,:) = Intersection_pts.pos(int_idx,:);
    b_name = ['boundary_no_', num2str(reg(2)), num2str(reg(1))];
    Boundaries.(b_name)(end+1,:) = Intersection_pts.pos(int_idx,:);
end

Boundaries = reorder_boundary(Boundaries);

figure
Line_Colors = jet(numel(fieldnames(Boundaries)));
boundary_names = fieldnames(Boundaries);
for b_idx = 1:numel(boundary_names)
    b_name = boundary_names(b_idx);
    boundary_pos = getfield(Boundaries, b_name{1});
    plot(boundary_pos(:,1), boundary_pos(:,2), '-', 'color', Line_Colors(b_idx,:),'Marker', '^', 'MarkerSize', 2);hold on;
end
plot(Intersection_pts.pos(:,1), Intersection_pts.pos(:,2),'r^','MarkerSize',6,...
    'MarkerFaceColor','r');

axis([min(lon), max(lon), min(lat), max(lat)]);

print(gcf, [pwd '/Generated_Plots/flow_region_boundaries'], '-r300', '-dpng');
close gcf
end

function [boundaries, intersection_pts] = boundary_points(Regions, field, delta_lon, delta_lat)

% Get boundary points in each region by checking whether the [i+1, j], [i-1, j],
% [i,j+1], [i,j-1] points are in the same region as [i,j]. If yes, then
% [i,j] is a boundary point.

for ii = 1:Regions.no % for each region
    for jj = 1: numel(Regions.lon{ii}) % for each point in the region
        cpos = [Regions.lon{ii}(jj), Regions.lat{ii}(jj)];
        cpos_region_no = ii;
        % neighbor positions
        neighbors(1,:) = cpos - [delta_lon, 0];
        neighbors(2,:) = cpos + [0, delta_lat];
        neighbors(3,:) = cpos + [delta_lon, 0];
        neighbors(4,:) = cpos - [0, delta_lat];
        
        [n, n_region_no] = check_neighbors(field, Regions, neighbors);
        
        if range([cpos_region_no, n_region_no]) == 0 % if the region # of current position
            % is the same as all of its neighbors
        else
            for kk = 1:numel(n_region_no)
                if cpos_region_no ~= n_region_no(kk) && cpos_region_no > n_region_no(kk)
                    % save the boundaries
                    fieldn = ['boundary_no_', num2str(cpos_region_no), num2str(n_region_no(kk))];
                    if ~exist('boundaries', 'var') || ~isfield(boundaries, fieldn)
                        boundaries.(fieldn) = 0.5*(cpos + n(kk, :));
                    else
                        boundaries.(fieldn)(end+1, :) = 0.5*(cpos + n(kk, :));
                    end
                end
            end         
            
            % Detect the intersection points
            sub_neighbors = [];
            mod_length = numel(n_region_no);
            
            for sub_n_idx = 1:numel(n_region_no)-1
                sub_neighbors(end+1,:) = [cpos_region_no, n_region_no(sub_n_idx), ...
                    n_region_no(sub_n_idx+1)];
            end
            sub_neighbors(end+1,:) = [cpos_region_no, n_region_no(end), n_region_no(1)];
            for sub_kk = 1:size(sub_neighbors, 1)
                if numel(unique(sub_neighbors(sub_kk, :))) == size(sub_neighbors, 2)
                    if sub_kk+1 > mod_length
                        n_idx = [sub_kk, mod(sub_kk+1, mod_length)];
                    else
                        n_idx = [sub_kk, sub_kk+1];
                    end
                    if ~exist('intersection_pts', 'var')
                        intersection_pts.pos = 0.5*(n(n_idx(1),:) + n(n_idx(2),:));
                        intersection_pts.reg_no = sub_neighbors(sub_kk,:);
                    else
                        intersection_pts.pos(end+1, :) = 0.5*(n(n_idx(1),:) + n(n_idx(2),:));
                        intersection_pts.reg_no(end+1, :) = sub_neighbors(sub_kk,:);
                    end
                end
            end
        end
    end
end
end

function Regions = check_regions(Regions, delta_lon, delta_lat)
% Check whether each of the regions identified by K-means is connected or
% not. If not, devide into two connected regions.

% Not implemented yet...

end

function [n, n_region_no] = check_neighbors(field, Regions, neighbors)
% First check whether the 4 neighbors of the current position is inside the
% domain. Then check the 4 neighbors belong to which region.

% n: the neighbors of current position that belong to the domain.
% n_region_no: the region of which the neighors belong to.

n = [];
for ii = 1:size(neighbors, 1)
    if (neighbors(ii, 1) >= min(field.lon)) && (neighbors(ii, 1) <= max(field.lon)) ...
            && (neighbors(ii,2) >= min(field.lat)) && (neighbors(ii,2) <= max(field.lat))
        n(end+1, :) = neighbors(ii, :);
    end
end

for ii = 1:size(n, 1) % for each of the neighbors
    for jj = 1:Regions.no % check if the neighbor belongs to region jj
        n1 = find(Regions.lon{jj} == n(ii,1));
        n2 = find(Regions.lat{jj} == n(ii,2));
        if ~isempty(intersect(n1, n2))
            n_region_no(ii) = jj;
        end
    end
    if numel(n_region_no) < ii
        error(['Cannot find the region that ', num2str(ii), 'th neighbor belongs.']);
    end
end
end

function boundaries_n = reorder_boundary(boundaries)

field_names = fieldnames(boundaries);
for b_idx = 1:numel(field_names)
    field_name = field_names{b_idx};
    boundary_position = getfield(boundaries, field_name);
    init_pos = [min(boundary_position(:,1)), max(boundary_position(:,2))];
    boundary_position_n = init_pos;
    for ii = 1:size(boundary_position,1)
        b_distance = [];
        for jj = 1:size(boundary_position,1)
            if intersect(find(boundary_position(jj,1) == boundary_position_n(:,1)), ...
                    find(boundary_position(jj,2) == boundary_position_n(:,2)))
                boundary_position(jj,:) = [NaN, NaN];
            end       
            b_distance(jj) = norm(boundary_position(jj,:) - init_pos);
        end
        b_distance(find(b_distance == 0)) = NaN;
        [~, neighbor_idx] = min(b_distance);
        init_pos = boundary_position(neighbor_idx,:);
        
        if isnan(init_pos(1)) && isnan(init_pos(2))
            break;
        end
        
%         if ii ~= 1
        if exist('boundary_position_n', 'var')
            boundary_position_n(end+1,:) = init_pos;
        else
            boundary_position_n(1,:) = init_pos;
        end
    end
    boundaries_n.(field_name) = boundary_position_n;
end
end