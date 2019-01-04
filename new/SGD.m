function [ x, config, iter, flg ] = SGD( x, config, thres, iter_max )

sigma = rand()*5;
flg = false;

% vf = func(x, config);
iter = 0;
% grad = gradient(x, config);
grad = ngradient(x, config);
gap = 1;
alpha = 0.1;

while iter <= iter_max
    
    if ~isreal(grad)
        grad
        flg = true;
        return;
    end
    
    iter = iter + 1;
%     iter
    % TODO: implement line search
%     grad
%     alpha = line_search(x, grad, config);
    dW = randn(size(grad));
    x_new = x - alpha * grad - sqrt(alpha) * sigma * dW;
%     sDW=sigma*dW'
%     x_new
%     configx = config.x
%     configreg = config.reg
%     grad
    [x_new, config] = check(x_new, config);
%     configx = config.x
%     x_new
%     vf = func(x_new, config)
%     waitforbuttonpress
%     if on_cross_pt
%         break;
%     end
%     vf_new = func(x_new, config)
%     grad = gradient(x_new, config);
    grad = ngradient(x_new, config);
%     new_grad = grad
    gap = norm(grad); % abs(vf - vf_new);
    x = x_new;
    
%     vf = vf_new;
    
end

end

