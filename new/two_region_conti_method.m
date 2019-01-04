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

Intersection_pts.pos = [];
Intersection_pts.reg_no = {};

tmp = {[tmp zeros(len,1); tmp ones(len,1)],
    
       [tmp 2*ones(len,1); tmp 3*ones(len,1); tmp 4*ones(len,1)]};

Regions0.no = 2;
Regions0.lon = cell(Regions0.no, 1);
Regions0.lat = cell(Regions0.no, 1);
for i = 1 : Regions0.no
    tt = tmp{i};
    Regions0.lon{i} = tt(:,1)';
    Regions0.lat{i} = tt(:,2)';
end
Regions0.vx = {0, 5};
Regions0.vy = {0, 0};
idx = 2;

Regions = Regions0;
Regions.vx{idx} = 0;

p = [0.5];
config.V = 2;
config.y = 4;

config.x = [0 1.5 1 1.5];
config.reg = [1 2];
config.u = [0 0; 0 0];

thres = 1e-6;
% config.u(2,1) = 1.5705;
% config_tmp = config;
% config_tmp.u = [0 0; 1.7 0];
% config_tmp.x = [1.0000    1.5000    2.0000    1.5000];
% [p1, config1, iter] = GD(0.5, config_tmp, thres);
% config1.x
% p1
% waitforbuttonpress

step = 0.02;
config.u = [0 0; 0 0];
while Regions.vx{idx} < Regions0.vx{idx} + step
    
    [p1, config1, iter, flg] = GD(p, config, thres);
    
    config1.x
    config1.u(idx,1)
    p1
    if flg
        fprintf('break\n');
        break;
    end
    Regions.vx{idx} = Regions.vx{idx} + step;
    config1.u(idx,1) = Regions.vx{idx};
    config = config1;
    p = p1;
    waitforbuttonpress
end
% [np, con, f] = check(p, p_new, config);

%%
config_tmp = config;
config_tmp.u(2,1) = config.u(2,1) - 0.00787;
func(p,config_tmp)


