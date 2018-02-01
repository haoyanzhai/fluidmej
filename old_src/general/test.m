clear
clc
close all

config.ini = 0;
config.ter = 0;

config.y = {@(x) 0 * x + 0;
%             @fexp;
            @(x) x - x + 2; 
            @(x) x - x + 3};
        
% config.y = {@(x) x - x + 0;
%             @(x) 0 * x + 1;
%             @(x) -0 * x + 2;
%             @(x) x - x + 3};
config.ux = [1.5; 0];
config.uy = [0; 0];
config.max_v = 2;

thres = 1e-8;

x0 = zeros(1, 1);
% x0 = [1; 0.5];
% x0 = [rand(1); rand(1)];

[x, iter] = GD0(x0, thres, config, @f);

x1 = [config.ini; x; config.ter];

y = get_y_coord(x1, config.y);


print_bd(config.y, -1, 1);
grid on;

scatter(x1, y);







