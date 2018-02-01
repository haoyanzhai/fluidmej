function y = get_y_coord( x, fy )

y = zeros(length(x), 1);

for i = 1 : length(x)
    y(i) = fy{i}(x(i));
end

end

