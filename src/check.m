function [ x_new, config, on_cross_pt ] = check( p, p_new, config )
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

x_new = p_new;

x_new(x_new > 1) = 1;
x_new(x_new < 0) = 0;
compare = x_new ~= p;
on_cross_pt = isempty(find(compare, 1));

coord = [(1-x_new).*config.x(:,1) + x_new.*config.x(:,3), ...
         (1-x_new).*config.x(:,2) + x_new.*config.x(:,4)];

global epsilon
     
combine_with_pre = [sum(abs(coord(2:end,:)-coord(1:end-1,:)),2) > epsilon; false];

% combine close junction
i = 1;
xx = [];
cx = [];
creg = [];
cu = [];
reduced = [];
while i <= length(combine_with_pre)
    if combine_with_pre(i)
        xx = [xx; x_new(i)];
        reg1 = config.reg(i,:);
        reg2 = config.reg(i+1,:);
        reg = setxor(reg1, reg2);
        if reg(1) > reg(2)
            tmpreg = reg(1); reg(1) = reg(2); reg(2) = tmpreg;
        end
        creg = [creg; reg];
        reduced = [reduced; true];
        cu = [cu; config.u(reg(1),:)];
        % find the boundary
        boundary = find_boundary(reg(1), reg(2));
        if norm(boundary(1,:) - coord(i,:)) < epsilon
            tmpx = [boundary(1,:) boundary(2,:)];
        else
            tmpx = [boundary(end,:) boundary(end-1,:)];
        end
        if x_new(i) == 1
            tmpxx = tmpx(1:2); tmpx(1:2) = tmpx(3:4); tmpx(3:4) = tmpxx;
        end
        cx = [cx; tmpx];
        i = i + 1;
    else
        xx = [xx; x_new(i)];
        cx = [cx; config.x(i,:)];
        creg = [creg; config.reg(i,:)];
        cu = [cu; config.u(i,:)];
        reduced = [reduced; false];
    end
    i = i + 1;
end

reduced_coord = [(1-xx).*cx(:,1) + xx.*cx(:,3), (1-xx).*cx(:,2) + xx.*cx(:,4)];

xxx = [];
ccx = [];
ccreg = [];
ccu = [];
% add junction
for i = 1 : length(reduced_coord)
    if ~reduced(i) && (xx(i) == 1 || xx(i) == 0)
        region = find_region(reduced_coord(i,:));
        % keep on the same boundary
        idx = -1;
        if length(region) == 2
            boundary = find_boundary(region(1), region(2));
            for j = 1 : size(boundary, 1)
                if norm(boundary(j, :) - reduced_coord(i, :)) < epsilon
                    idx = j; break;
                end
            end
            xxx = [xxx; xx(i)];
            ccreg = [ccreg; creg(i, :)];
            ccu = [ccu; cu(i, :)];
            if xx(i) == 0
                tmp = cx(i, 3:4);
            else
                tmp = cx(i, 1:2);
            end
            if idx > 1 && norm(boundary(idx-1, :) - tmp) < epsilon
                if xx(i) == 0
                    ccx = [ccx; boundary(idx, :) boundary(idx-1, :)];
                else
                    ccx = [ccx; boundary(idx-1, :) boundary(idx, :)];
                end
            else
                if xx(i) == 0
                    ccx = [ccx; boundary(idx, :) boundary(idx+1, :)];
                else
                    ccx = [ccx; boundary(idx+1, :) boundary(idx, :)];
                end
            end
        else % new region are to be added
            for j = 1 : length(region) - 1
                xxx = [xxx; xx(i)];
                ccreg = [ccreg; region(i) region(i+1)];
                ccu = [ccu; get_velocity(region(i))];
                boundary = find_boundary(region(i), region(i+1));
                if norm(boundary(1,:) - reduced_coord(i,:)) < epsilon
                    tmpx = [boundary(1,:) boundary(2,:)];
                else
                    tmpx = [boundary(end,:) boundary(end-1,:)];
                end
                if xx(i) == 1
                    tmpxx = tmpx(1:2); tmpx(1:2) = tmpx(3:4); tmpx(3:4) = tmpxx;
                end
                ccx = [ccx; tmpx];
            end
        end
    else
        xxx = [xxx; xx(i)];
        ccx = [ccx; cx(i,:)];
        ccreg = [ccreg; creg(i,:)];
        ccu = [ccu; cu(i,:)];
    end
end

x_new = xxx;
config.x = ccx;
config.reg = ccreg;
config.u = ccu;

end

