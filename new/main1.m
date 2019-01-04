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

%      ?.   /    8
%      ?  ./
%      ?  / .
%      ? /   .
% -3 ------------4
%      ?     .
%      ?   .
%      ? .
%      0         0

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

Boundaries.boundary_no_12 = [2.5 9.5; 10 9.5];
Boundaries.boundary_no_13 = [-10 9.5; 2.5 9.5];
Boundaries.boundary_no_23 = [2.5 9.5; 10 17];

Intersection_pts.pos = [2.5 9.5];
Intersection_pts.reg_no = {[1 2 3]};

ttt = (-10:10)';
yyy = ones(size(ttt));

tttt = (-10:0)';
yyyy = ones(size(tttt));
tmp = {[ttt 9*yyy; ttt 8*yyy; ttt 7*yyy],
    
       [3 10; 4 10; 4 11; 5 11; 5 12; 6 12; 6 13; 7 13; 7 14; 8 14; 8 15; 9 15; 9 16; 10 16; 10 17],
        
       [1 11; 1 10; 2 10; 2 11; 3 11; 3 12; 4 12; 4 13; 5 13; 5 14; 6 14; 6 15; 7 15; 7 16; 8 16;
        8 17; 9 17; 9 18; 10 18; 10 19; 10 20;
        tttt 10*yyyy; tttt 11*yyyy]};

Regions.no = 3;
Regions.lon = cell(Regions.no, 1);
Regions.lat = cell(Regions.no, 1);
for i = 1 : Regions.no
    tt = tmp{i};
    Regions.lon{i} = tt(:,1)';
    Regions.lat{i} = tt(:,2)';
end
Regions.vx = {1, 1, 1};
Regions.vy = {1, 1, 1};

config.u = [1 1; 1 1; 1 1];
config.x = [2.5 9.5 10 9.5; 2.5 9.5 10 17];
config.reg = [1 2; 2 3];
config.C = 1;

config.V = 3;
config.y = 20;

p = [0.5; 0.5];


thres = 1e-5;
% alpha = 0.01;
[x, config_final, iter] = GD(p, config, thres);


p = x
func(x, config_final)
coorddd = get_coord(p, config_final.x)
timeeee = get_time_cost(p, config_final)


