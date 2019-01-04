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

Boundaries.boundary_no_12 = [1 1.5];
Boundaries.boundary_no_13 = [1.5 2];
Boundaries.boundary_no_23 = [1.5 1.5];

Intersection_pts.pos = [1.5 1.5];
Intersection_pts.reg_no = {[1 2 3]};

Regions.no = 3;
Regions.lon = {[1 1], [1], [2]};
Regions.lat = {[1 2], [2], [2]};
Regions.vx = {0, 0, 0};
Regions.vy = {0, 0, 0};

a = find_region(Regions, 1.5, 1.5)

% p = [0; 1];
% p_new = [-0.1; 1.2];
% 
% config.x = [4.5 2 4.5 1.5; 5 1.5 4.5 1.5];
% config.reg = [2 3; 1 3];
% config.V = 5;
% config.u = [0 0; 0 0];
% config.y = 2;

% [np, con, f] = check(p, p_new, config);




