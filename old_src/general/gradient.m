function grad = gradient( func, x, vf, config )

eps = 1e-6;
J = zeros(length(vf), length(x));
for i = 1 : length(x)
    tmp = x;
    tmp(i) = tmp(i) + eps;
    J(:, i) = (func(tmp, config) - vf) / eps;
end

grad = J';

end

