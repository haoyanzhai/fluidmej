function [ x, combine, end_loop, config, iter, comb ] = GD2( x0, thres, config, func )

combine = [];
x = x0;
vf = func(x, config);
grad = gradient(func, x, vf, config);

iter = 0;
end_loop = false;
fprintf('%d iter, vf %.8f error %.8f x size %d\n', iter, vf, norm(grad), length(x));
comb(:, 1) = x;
while norm(grad) > thres && ~end_loop
    iter = iter + 1;
%     alpha = line_search(func, vf, x, config, -grad);
    alpha = 0.01;
    x = x - alpha * grad;
    comb(:, iter + 1) = x;
    vf = func(x, config);
    grad = gradient(func, x, vf, config);
    
    [end_loop, x, config] = check2(x, config);
    
    fprintf('%d iter, vf %.8f error %.8f x size %d\n', iter, vf, norm(grad), length(x));
end

end

