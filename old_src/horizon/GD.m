function [ x, iter ] = GD( x0, thres, config, func )

x = x0;
vf = func(x, config);
grad = gradient(func, x, vf, config);

iter = 0;
while norm(grad) > thres
    iter = iter + 1;
%     alpha = line_search(func, vf, x, config, -grad);
    alpha = 0.01;
    x = x - alpha * grad;
    vf = func(x, config);
    grad = gradient(func, x, vf, config);
    fprintf('%d iteration with step size %.8f error %.8f\n', iter, alpha, norm(grad));
end

end

