function [ lambda, config ] = check( p_new, config )
% check if new junction needs to be added or old ones should be combined
% 
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

global Boundaries
global Regions

start_lambda = 0.00001;

lambda = [];
reg = [];
x = [];
u = [];

current_reg = 1;
idx = 1;
while idx <= size(p_new, 1)
    pivot = [];
    sec = [];
    if p_new(idx) < 0
        p_new(idx) = 0;
        pivot = config.x(idx, 1:2);
        sec = config.x(idx, 3:4);
    elseif p_new(idx) > 1
        p_new(idx) = 1;
        pivot = config.x(idx, 3:4);
        sec = config.x(idx, 1:2);
    end
    
    if size(pivot,1) == 0
        lambda = [lambda; p_new(idx)];
        reg = [reg; config.reg(idx,:)];
        x = [x; config.x(idx, :)];
        current_reg = config.reg(idx, 2);
        u = [u; config.u(idx,:)];
    else % pivot ~= empty
        region = find_region(Regions, pivot(1), pivot(2));
        if length(region) == 2 % no need to add new junction
            lambda = [lambda; p_new(idx)];
            reg = [reg; config.reg(idx,:)];
            current_reg = config.reg(idx, 2);
            x = [x; find_next_on_bd(pivot, sec, region, Boundaries, lambda)];
            u = [u; config.u(idx,:)];
        else % add junction
%             fprintf('add junction\n');
            comb = combnk(region, 2);
            comb = delete_same(comb, config.reg(idx,:));
            while ~isempty(comb)
                [next_reg, comb] = find_next_reg(comb, current_reg);
                reg = [reg; next_reg];
                u = [u; Regions.vx{current_reg} Regions.vy{current_reg}];
                current_reg = next_reg(2);
                if p_new(idx) == 1
                    lambda = [lambda; p_new(idx)-start_lambda];
                else
                    lambda = [lambda; start_lambda];
                end
                x = [x; get_new_bd_x(Boundaries, next_reg, pivot, lambda(end))];
            end
        end
    end
    idx = idx + 1;
end
% reg
[reg, u, x, lambda] = check_cycle(reg, u, x, lambda);
% regafter = reg
config.reg = reg;
config.u = [u; config.u(end,:)];
config.x = x;

end





function next = find_next_on_bd(pivot, sec, region, Boundaries, lambda)

bd = find_boundary(Boundaries, min(region(1),region(2)), max(region(1),region(2)));
for i = 1 : size(bd, 1)
    if norm(bd(i,:) - pivot) < 1e-5 % pivot position found
        sec_is_head = false;
        if i ~= 1
            sec_is_head = norm(bd(i-1,:) - sec) < 1e-5;
        end
        if sec_is_head
            if i < size(bd, 1)
                next = bd(i+1,:);
            else
                next = sec;
            end
            break;
        else % sec_is_tail
            if i ~= 1
                next = bd(i-1,:);
            else
                next = sec;
            end
            break;
        end
    end
end

if lambda == 0
    next = [pivot next];
else
    next = [next pivot];
end

end



function comb = delete_same(comb, reg)

for i = 1 : size(comb, 1)
    if (comb(i,1)==reg(1) && comb(i,2)==reg(2)) || (comb(i,2)==reg(1) && comb(i,1)==reg(2))
        comb(i,:) = [];
        break;
    end
end

end



function [next_reg, comb] = find_next_reg(comb, current_reg)

for i = 1 : size(comb, 1)
    if current_reg == comb(i,1)
        next_reg = [current_reg comb(i,2)];
        comb(i,:) = [];
        break;
    elseif current_reg == comb(i,2)
        next_reg = [current_reg comb(i,1)];
        comb(i,:) = [];
        break;
    end
end

end



function x = get_new_bd_x(Boundaries, next_reg, pivot, lambda)

bd = find_boundary(Boundaries, min(next_reg(1),next_reg(2)), max(next_reg(1),next_reg(2)));

if norm(pivot - bd(1,:)) < 1e-5
    sec = bd(2,:);
else
    sec = bd(end-1,:);
end

if lambda == 0
    x = [pivot sec];
else
    x = [sec pivot];
end

end



function [reg, u, x, lambda] = check_cycle(reg, u, x, lambda)

map = zeros(max(max(reg)), 1);

i = 1;
while i <= size(reg,1)
    if map(reg(i,2)) ~= 0 % cycle found
        u(map(reg(i,2)):i, :) = [];
        x(map(reg(i,2)):i, :) = [];
        lambda(map(reg(i,2)):i, :) = [];
        reg(map(reg(i,2)):i, :) = [];
        map = zeros(max(max(reg)), 1);
        i = 1;
    else % not found
        map(reg(i,1)) = i;
        i = i + 1;
    end
end

end

