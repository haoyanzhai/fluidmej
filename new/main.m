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

% \  ?.
%  \ ? .
%   \?  .
%    \   .
%    ?\   .
%    ? \   .
%    ?  \   .
% --------------
%    ?     .
%    ?   .
%    ? .
%    0

%  ... initial
% ? final
%  | --: boundary

clear
clc
close all

global Boundaries
global Intersection_pts
global Regions
global epsilon

epsilon = 1e-5;

Boundaries.boundary_no_12 = [3 2.5; -2 2.5; -1 2.5; 0 2.5; 1 2.5];
Boundaries.boundary_no_13 = [1 2.5; 2 2.5; 3 2.5; 4 2.5];
Boundaries.boundary_no_23 = [-3 6.5; -2.5 6; -2 5.5; -1.5 5; -1 4.5; -0.5 4; 0 3.5; 1 2.5];

Intersection_pts.pos = [1 2.5];
Intersection_pts.reg_no = {[1 2 3]};

tmp = {[-3 1; -2 1; -1 1; 0 1; 1 1; 2 1; 3 1; 4 1;
        -3 2; -2 2; -2 2; 0 2; 1 2; 2 2; 3 2; 4 2;],
    
       [-3 3; -2 3; -1 3; 0 3;
        -3 4; -2 4; -1 4;
        -3 5; -2 5;
        -3 6],
        
       [-3 8; -2 8; -1 8; 0 8; 1 8; 2 8; 3 8; 4 8;
        -3 7; -2 7; -1 7; 0 7; 1 7; 2 7; 3 7; 4 7;
        -2 6; -1 6; 0 6; 1 6; 2 6; 3 6; 4 6;
        -1 5; 0 5; 1 5; 2 5; 3 5; 4 5;
        0 4; 1 4; 2 4; 3 4; 4 4;
        1 3; 2 3; 3 3; 4 3]};

Regions.no = 3;
Regions.lon = cell(Regions.no, 1);
Regions.lat = cell(Regions.no, 1);
for i = 1 : Regions.no
    tt = tmp{i};
    Regions.lon{i} = tt(:,1)';
    Regions.lat{i} = tt(:,2)';
end
Regions.vx = {0, 0, 0};
Regions.vy = {0, 0, 0};

config.x = [2 2.5 3 2.5];
config.reg = [1 3];
config.u = [0 0; 0 0];

config.V = 3;
config.y = 8;

p = [0.5];


thres = 1e-5;
[x, iter] = GD(p, config, thres, 1);

% [np, con, f] = check(p, p_new, config);




