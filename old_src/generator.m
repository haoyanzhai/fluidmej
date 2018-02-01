clear
clc
close all

xx = 1 : 100;
yy = 1 : 100;

[x, y] = meshgrid(xx, yy);

vx = zeros(size(x));
vy = vx;

vx(1 : 30, :) = 0.1;
vx(31 : 60, :) = -0.1;
vx(61 : end, :) = 0.1;

field.x = x;
field.y = y;
field.vx = vx;
field.vy = vy;

quiver(x,y,vx,vy);hold on;
axis([0,100,0,100]);
grid on;





