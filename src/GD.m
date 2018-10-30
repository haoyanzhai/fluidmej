function [ x, iter ] = GD( x, gradient, check, config, func, thres, alpha, get_config )

if nargin == 6
    alpha = 0.5;
end

vf = func(x, config);
iter = 0;
grad = gradient(func, vf, x, config);
gap = 1;

while gap > thres
    iter = iter + 1;
    % TODO: implement line search
    x_new = x - alpha * grad;
    [x_new, config, on_cross_pt] = check(x, x_new, config, get_config); % TODO
    if on_cross_pt
        break;
    end
    vf_new = func(x_new, config);
    grad = gradient(func, vf_new, x, config);
    gap = norm(x - x_new); % abs(vf - vf_new);
    x = x_new;
%     vf = vf_new;
end

end

