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

% 
% --------------------------
% --------------------------

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

tmp = (-7 : 7)';
len = length(tmp);

Boundaries.boundary_no_12 = [tmp 1.5*ones(len,1)];
Boundaries.boundary_no_23 = [tmp 4.5*ones(len,1)];

Intersection_pts.pos = [];
Intersection_pts.reg_no = {};

tmp = {[tmp zeros(len,1); tmp ones(len,1)],
    
       [tmp 2*ones(len,1); tmp 3*ones(len,1); tmp 4*ones(len,1)],
        
       [tmp 5*ones(len,1); tmp 6*ones(len,1); tmp 7*ones(len,1); tmp 8*ones(len,1)]};

Regions0.no = 3;
Regions0.lon = cell(Regions0.no, 1);
Regions0.lat = cell(Regions0.no, 1);
for i = 1 : Regions0.no
    tt = tmp{i};
    Regions0.lon{i} = tt(:,1)';
    Regions0.lat{i} = tt(:,2)';
end
Regions0.vx = {0, 2, 0};
Regions0.vy = {0, 0, 0};

Regions = Regions0;
Regions.vx{2} = 0;

p = [0.5; 0.5];
config.V = 1;
config.y = 8;

config.x = [1 1.5 2 1.5; 3 4.5 4 4.5];
config.reg = [1 2; 2 3];
config.u = [0 0; 0 0; 0 0];

thres = 1e-5;
% [p1, config1, iter] = GD(p, config, thres, 0.5);
% config1.x
% p1
% waitforbuttonpress
while Regions.vx{2} ~= Regions0.vx{2}
    
    [p1, config1, iter] = GD(p, config, thres);
    
    config1.x
    config.u
    p1
    
    Regions.vx{2} = Regions.vx{2} + 0.2;
    config1.u(2,1) = Regions.vx{2};
    config = config1;
    p = p1;
    waitforbuttonpress
end
% [np, con, f] = check(p, p_new, config);




