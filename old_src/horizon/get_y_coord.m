function y = get_y_coord( x, fy )

y = zeros(length(fy), 1);
for i = 1 : length(fy)
    y(i) = fy{i}(x(i));
end

end

