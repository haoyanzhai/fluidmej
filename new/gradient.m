function grad = gradient( p, config )
% this is the function to calculate the gradient of the objective function value
%
% p     :   parameters of junctions [p1 ... pn]', where the coordinates of
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

% max speed
V = config.V;

% ax ay are the boundary pts for lambda = 0 and bx by are when lambda = 1
% length = n-1
ax = config.x(:,1);
ay = config.x(:,2);
bx = config.x(:,3);
by = config.x(:,4);

% flow velocity
% length = n
ux = config.u(:,1);
uy = config.u(:,2);

% junction positions
% length = n
x = p .* (bx - ax) + ax;
y = p .* (by - ay) + ay;
x = [x; 0] - [0; x];
y = [y; config.y] - [0; y];

% t_i = C_i / (A_i - B_i)
% length = n
A = sqrt( (ux.*x+uy.*y).^2 + (x.^2+y.^2) .* (V^2-ux.^2-uy.^2) );
B = ux.*x + uy.*y;
C = x.^2 + y.^2;

% dx_i/dlambda_i
dxidpi = bx - ax;
dyidpi = by - ay;

% dx_i+1/dlambda_i
dxi1dpi = ax - bx;
dyi1dpi = ay - by;

% for dABC_i/dlambda_i and dt_i/dlambda_i
% i = 1 : n-1
ux0 = ux(1:end-1);
uy0 = uy(1:end-1);
x0 = x(1:end-1);
y0 = y(1:end-1);
A0 = A(1:end-1);
B0 = B(1:end-1);
C0 = C(1:end-1);

dAidpi = ((ux0.*x0+uy0.*y0) .* (dxidpi+dyidpi) + (x0.*dxidpi + y0.*dyidpi) .* (V^2-ux0.^2-uy0.^2)) ./ A0;
dBidpi = ux0.*dxidpi + uy0.*dyidpi;
dCidpi = 2*x0.*dxidpi + 2*y0.*dyidpi;

dtidpi = (dCidpi .* (A0-B0) - C0.*(dAidpi-dBidpi)) ./ (A0 - B0).^2;

ux1 = ux(2:end);
uy1 = uy(2:end);
x1 = x(2:end);
y1 = y(2:end);
A1 = A(2:end);
B1 = B(2:end);
C1 = C(2:end);

dAi1dpi = ((ux1.*x1+uy1.*y1).*(dxi1dpi+dyi1dpi) + (x1.*dxi1dpi+y1.*dyi1dpi).*(V^2-ux1.^2-uy1.^2)) ./ A1;
dBi1dpi = ux1.*dxi1dpi + uy1.*dyi1dpi;
dCi1dpi = 2*x1.*dxi1dpi + 2*y1.*dyi1dpi;

dti1dpi = (dCi1dpi .* (A1-B1) - C1.*(dAi1dpi-dBi1dpi)) ./ (A1 - B1).^2;

grad = dtidpi + dti1dpi;

% A
% B
% C
% dAidpi
% dBidpi
% dCidpi
% dAi1dpi
% dBi1dpi
% dCi1dpi
% dxidpi
% dyidpi
% dxi1dpi
% dyi1dpi
% ax
% ay
% bx
% by

end

