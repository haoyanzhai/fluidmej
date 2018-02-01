clear
clc
close all

config.ini = 0;
config.ter = 0;


thres = 1e-8;

% get initial
x0 = [1];

config.ini = 0;
config.ter = 0;

total = {@(x) 0 * x + 0;
            @fexp;
            @(x) x - x + 2; 
            @(x) x - x + 3};

config.y = {@(x) 0 * x + 0;
            @(x) x - x + 2; 
            @(x) x - x + 3};

config.ux = [-1; 0];
config.uy = [0; 0];

config.max_v = 2;

% outer loop
done = false;
x = x0;
l = 1;
tc{l} = config;
while ~done
    [x, combine, end_loop, config, iter, comb] = GD1(x, thres, config, @f);
    tc{l + 1} = config;
    tmp{l} = comb;
    l = l + 1;
    done = ~end_loop;
end

% plot result
figure(1);
x1 = [config.ini; x; config.ter];
y = get_y_coord(x1, config.y);
print_bd(total, -1, 1);
grid on;

scatter(x1, y);

%% plot animation
figure(2);

for i = 1 : length(tmp)
    c = tmp{i};
    co = tc{i};
    for j = 1 : 1 : size(c, 2)
        x1 = [co.ini; c(:, j); co.ter];
        y = get_y_coord(x1, co.y);
        hold off;
        print_bd(total, -1, 1);
        grid on;
        scatter(x1, y);
        pause(0.01);
    end
end






