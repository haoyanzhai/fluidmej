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

% p,p_new:  parameters of junctions [p1 ... pn]', where the coordinates of
%           junctions can be represented as
%               xi = (1 - pi) * xi1 + pi * xi2
%               yi = (1 - pi) * yi1 + pi * yi2
%           where [xi1 xi2 yi1 yi2] are defined in config and pi\in[0,1].
%           Combined with start and end points, the total coordinate vector
%           is [0 x1 ... xn 0], [0 y1 ... yn end_y_coordinate]
% config:   config.x:   each line is two points coordinates of the
%                       endpoints of junction located line
%                       [x11 y11 x12 y12
%                        x21 y21 x22 y22
%                        ... ... ... ...
%                        x(n-1)1 y(n-1)1 x(n-1)2 y(n-1)2]
%           config.reg: [r11 r12
%                        r21 r22
%                        ... ...]
%                        r(n-1)1 r(n-1)2]
%           config.V:   the maximum speed of vehicle
%           config.u:   flow velocity
%                       [u1x u1y
%                        u2x u2y
%                        ... ...
%                        unx uny]
%           config.y:   final y coordinate


clear
clc
close all

global Boundaries
global Intersection_pts
global Regions
global epsilon

epsilon = 1e-5;

Boundaries.boundary_no_12 = [-10 4.5; -2.5 4.5];
Boundaries.boundary_no_13 = [-2.5 4.5; 2.5 4.5];
Boundaries.boundary_no_14 = [2.5 4.5; 10 4.5];
Boundaries.boundary_no_23 = [-2.5 4.5; -2.5 14.5];
Boundaries.boundary_no_34 = [2.5 4.5; 2.5 14.5];
Boundaries.boundary_no_25 = [-10 14.5; -2.5 14.5];
Boundaries.boundary_no_35 = [-2.5 14.5; 2.5 14.5];
Boundaries.boundary_no_45 = [2.5 14.5; 10 14.5];

Intersection_pts.pos = [-2.5 4.5; 2.5 4.5; -2.5 14.5; 2.5 14.5];
Intersection_pts.reg_no = {[1 2 3], [1 3 4], [2 3 5], [3 4 5]};

tmp0 = (-10 : -3)';
tmp1 = (-2 : 2)';
tmp3 = (3 : 10)';
tmp4 = (-10 : 10)';

tmp = cell(5, 1);

% region 1
for i = 0 : 4
    tmp{1} = [tmp{1}; tmp4 ones(size(tmp4))*i];
end

% region 2
for i = 5 : 14
    tmp{2} = [tmp{2}; tmp0 ones(size(tmp0))*i];
end

% region 3
for i = 5 : 14
    tmp{3} = [tmp{3}; tmp1 ones(size(tmp1))*i];
end

% region 4
for i = 5 : 14
    tmp{4} = [tmp{4}; tmp3 ones(size(tmp3))*i];
end

% region 5
for i = 15 : 20
    tmp{5} = [tmp{5}; tmp4 ones(size(tmp4))*i];
end


Regions.no = 5;
Regions.lon = cell(Regions.no, 1);
Regions.lat = cell(Regions.no, 1);
for i = 1 : Regions.no
    tt = tmp{i};
    Regions.lon{i} = tt(:,1)';
    Regions.lat{i} = tt(:,2)';
end
Regions.vx = {0, -1.5, 0, 1.5, 0};
Regions.vy = {0, 0, -2, 0, 0};

% right
config.u = [0 0; 1.5 0; 0 0];
config.x = [10 4.5 2.5 4.5; 10 14.5 2.5 14.5];
config.reg = [1 4; 4 5];

% left
% config.u = [0 -0.5; 1.5 0; 0 -0.5];
% config.x = [-10 4.5 -2.5 4.5; -10 14.5 -2.5 14.5];
% config.reg = [1 2; 2 5];

% center
% config.u = [0 -0.5; 0 2; 0 -0.5];
% config.x = [-2.5 4.5 2.5 4.5; -2.5 14.5 2.5 14.5];
% config.reg = [1 3; 3 5];


config.V = 3;
config.y = 20;
config.C = 20;
p = [0.6399; 0.81223];

alpha = 0.001;
thres = 1e-6;
for kkk = 1 : 20

    fprintf('grad\n');
    [p1, config1, ~] = GD(p, config, thres, alpha);
    p1
    config1.x
    func(p1, config1)
    coorddd = get_coord(p1, config1.x)
    timeeee = get_time_cost(p1, config1)
    fprintf('perturb\n');
    waitforbuttonpress
    [p2, config2, iter] = SGD(p1, config1, thres, 5);
    p2
    config2.x
    waitforbuttonpress
    
    p = p2;
    config = config2;

end

[x, config_final, iter] = GD(p, config, thres);

x
final = config_final.x
coorddd = get_coord(x, config_final.x)
timeeee = get_time_cost(x, config_final)



