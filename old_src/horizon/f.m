function vf = f( x, config )

x = [config.ini; x; config.ter];
ux = config.ux;
uy = config.uy;
fy = config.y;
a = config.max_v;

y = get_y_coord(x, fy);

dx = x(2 : end) - x(1 : end-1);
dy = y(2 : end) - y(1 : end-1);

xnorm = sqrt(dx.^2 + dy.^2);
C = (uy.*dx - ux.*dy) / a ./ xnorm;

cosalpha = dy ./ xnorm;
sinalpha = dx ./ xnorm;

vx = a * (C.*cosalpha + sqrt(1-C.^2).*sinalpha);
vy = a * (sqrt(1-C.^2).*cosalpha - C.*sinalpha);

vf = sum( xnorm ./ sqrt((vx+ux).^2 + (vy+uy).^2) );

end

