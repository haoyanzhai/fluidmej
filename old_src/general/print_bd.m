function [ ] = print_bd( fy, xlb, xub )

x = linspace(xlb, xub, 100)';

for i = 1 : length(fy)
    y = fy{i}(x);
    plot(x, y); hold on;
    axis([-1 1 0 3]);
end

end

