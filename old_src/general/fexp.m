function y = fexp( x )

y = zeros(size(x));

y(x > 0.4) = 2;
y(x <= 0.4) = 2 / 1.4^2 * (x(x <= 0.4) + 1).^2;

end

