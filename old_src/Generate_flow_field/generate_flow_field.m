clear all;
close all;
clc;

figure(1);

%% UAV speed
vg = 1;
%% grid
x = 1:100;
y = 1:100;

[xx,yy] = meshgrid(x,y);
vcx1 = zeros(100,100);
vcy1 = zeros(100,100);
%% flow 
for ii = 1:20
    for jj = 1:100
    vcx1(ii,jj) = 0.8;
    vcy1(ii,jj) = 0;
    end
end

for ii = 30:50
    for jj = 23:55
        vcx1(ii,jj) = 1.5;
        vcy1(ii,jj) = 0;
    end
end

for ii = 60:90
    for jj = round(55-sqrt(15^2-(ii-75)^2)):round(55+sqrt(15^2-(ii-75)^2))
        vcx1(ii,jj) = -1.5;
        vcy1(ii,jj) = 1;
    end
end

for ii = 70:100
    for jj = 70+(100-ii):100
        vcx1(ii,jj) = -1;
        vcy1(ii,jj) = 1;
    end
end

%% plot
quiver(xx,yy,vcx1,vcy1);hold on;
axis([0,100,0,100]);
grid on;

%% contour of the obstacles

% the circle:
% ctr = (55,75) r = 15
theta = -pi:0.1*pi:pi;
plot(15*cos(theta)+55,15*sin(theta)+75,'k-');hold on;

% the rectangular:
% ctr = (39,40) a = 32 b = 20
x_rec = [23, 55, 55, 23, 23];
y_rec = [30, 30, 50, 50, 30];
plot(x_rec, y_rec, '-k');hold on;
% the upper right corner triangle
% vertices = [70,100],[100,100],[100,70]
tri_x = [70, 100, 100, 70];
tri_y = [100, 100, 70, 100];
plot(tri_x, tri_y, '-k');hold on;
%% second figure
figure(2)
vcx2 = zeros(100,100);
vcy2 = zeros(100,100);

for ii = 1:100
    for jj = max(ii-30,1):min(ii+30,100)
        vcx2(ii,jj) = sin(1/20*(ii+jj));
        vcy2(ii,jj) = sin(1/20*(ii+jj));
    end
end

quiver(xx,yy,vcx2,vcy2);hold on;
grid on;
axis([1,100,1,100]);


%% contour

obs = zeros(100,100);
for ii = 1:100
    for jj = 1:100
        if (vcx2(ii,jj)^2+vcy2(ii,jj)^2>vg)
            obs(ii,jj)=1;
        end
    end
end
% [C,h] = contour(obs);

con1_x = [24, 40, 70, 54, 24];
con1_y = [55, 70, 40, 24, 55];

con2_x = con1_x - ones(1,5)*31;
con2_y = con1_y - ones(1,5)*31;

con3_x = con1_x + ones(1,5)*31;
con3_y = con1_y + ones(1,5)*31;

plot(con1_x,con1_y,'-k');
plot(con2_x,con2_y,'-k');
plot(con3_x,con3_y,'-k');