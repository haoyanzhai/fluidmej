function [boundaries, intersection_pts] = flow_boundaries(field)
% from the flow field data, first divide the flow field into regions
% according to the flow speed information at each position. Then find the
% boundaries of the regions.

% field.lon/ field.lat: Lon/Lat grid vectors;
% field.vcx/ field.vcy: horizontal/vertical flow data at each position.

% Mengxue Hou
% Georgia Institute of Technology, 2018/4

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

X = [vcx, vcy, gridx/nx, gridy/ny];

% for k=2:5
k = 6;
opts = statset('Display','final');
[idx,C] = kmeans(X,k,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

figure
Marker = ['r.'; 'b.';'g.';'c.';'y.';'k.'];
for ii = 1:k
    Regions.no = ii;
    Regions.lon{ii} = round(X(idx==ii,3)*nx);
    Regions.lat{ii} = round(X(idx==ii,4)*ny);
    Regions.vx{ii} = X(idx==ii,1);
    Regions.vy{ii} = X(idx==ii,2);
    
    plot(Regions.lon{ii},Regions.lat{ii}, Marker(ii,:));hold on;
    plot(C(ii,3)*nx, C(ii,4)*ny, 'ko');hold on;
end
axis([min(lon), max(lon), min(lat), max(lat)]);

print(gcf, [pwd '/Generated_Plots/Divided_flow_field'], '-r300', '-dpng');
close gcf
% end

%% Check whether each of the divided regions are connected
Regions = check_regions(Regions, delta_lon, delta_lat);

%% Find boundaries of the flow regions
% The boundaries are the center line between two neighboring regions
[boundaries, intersection_pts] = boundary_points(Regions, field, delta_lon, delta_lat);

figure
Line_Colors = jet(numel(fieldnames(boundaries)));
boundary_names = fieldnames(boundaries);
for b_idx = 1:numel(boundary_names)
    b_name = boundary_names(b_idx);
    boundary_pos = getfield(boundaries, b_name{1});
    plot(boundary_pos(:,1), boundary_pos(:,2), '.', 'color', Line_Colors(b_idx,:));hold on;
end
plot(intersection_pts.pos(:,1), intersection_pts.pos(:,2),'r^','MarkerSize',6,...
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