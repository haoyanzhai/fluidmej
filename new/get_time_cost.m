function vf = get_time_cost( p, config )
% this is the function to calculate the objective function value
%
% p     :   parameters of junctions [p1 ... pn-1 t1 ... tn]', where the coordinates of
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
%           config.V:   the maximum speed of vehicle
%           config.u:   flow velocity
%                       [u1x u1y
%                        u2x u2y
%                        ... ...
%                        unx uny]
%           config.y:   final y coordinate

dcx = config.x(:,3) - config.x(:,1);
dcy = config.x(:,4) - config.x(:,2);
x = config.x(:,1) + p.*dcx;
y = config.x(:,2) + p.*dcy;
sx = [x; 0] - [0; x];
sy = [y; config.y] - [0; y];
ux = config.u(:,1);
uy = config.u(:,2);

t = set_t(sx, sy, config.C, ux, uy, config.V^2);
vf = sum(t);

end




function t = set_t(sx, sy, C, ux, uy, Vsquare)

t = zeros(size(sx));
tmp = C + 2*(ux.^2+uy.^2) - 2*sqrt(C+ux.^2+uy.^2).*(sx.*ux+sy.*uy)./sqrt(sx.^2+sy.^2);
first = sqrt(sx.^2+sy.^2) ./ sqrt(C+ux.^2+uy.^2);
part1 = sx.*ux+sy.*uy;
part2 = Vsquare-ux.^2-uy.^2;
% ux
% uy
% Vsquare
second = (part1 - sqrt( part1.^2 + (sx.^2+sy.^2).*part2 )) ./ (-part2);
third = (sx.^2+sy.^2) ./ (2 * (sx.*ux + sy.*uy));
second(part2 == 0) = third(part2 == 0);
for i = 1 : length(sx)
    if tmp(i) <= Vsquare
        t(i) = first(i);
    else
        t(i) = second(i);
    end
end

end




