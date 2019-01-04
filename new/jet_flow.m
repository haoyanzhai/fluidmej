clear
clc
close all

global Boundaries
global Intersection_pts
global Regions
global epsilon

epsilon = 1e-5;

bd = 10;
Boundaries.boundary_no_12 = [-bd 7.5; bd 7.5];
Boundaries.boundary_no_23 = [-bd 10.5; bd 10.5];

Intersection_pts.pos = [];
Intersection_pts.reg_no = {};

tmp0 = (-bd : bd)';

tmp = cell(3, 1);

% region 1
for i = 0 : 7
    tmp{1} = [tmp{1}; tmp0 ones(size(tmp0))*i];
end

% region 2
for i = 8 : 10
    tmp{2} = [tmp{2}; tmp0 ones(size(tmp0))*i];
end

% region 3
for i = 10 : 20
    tmp{3} = [tmp{3}; tmp0 ones(size(tmp0))*i];
end


Regions.no = 3;
Regions.lon = cell(Regions.no, 1);
Regions.lat = cell(Regions.no, 1);
for i = 1 : Regions.no
    tt = tmp{i};
    Regions.lon{i} = tt(:,1)';
    Regions.lat{i} = tt(:,2)';
end
Regions.vx = {0, 2.9, 0};
Regions.vy = {0, 0, 0};


% right
config.u = [0 0; 2.9 0; 0 0];
config.x = [-bd 7.5 bd 7.5; -bd 10.5 bd 10.5];
config.reg = [1 2; 2 3];

config.V = 3;
config.y = 20;
config.C = 1;
p = [0.5; 0.45];


thres = 1e-4;
alpha = 0.0001;
%% sgd
for kkk = 1 : 20
    
    kkk
    fprintf('grad\n');
    [p1, config1, iter] = GD(p, config, thres, alpha);
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

func(x, config_final)



