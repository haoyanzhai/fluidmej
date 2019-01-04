function [ x, config, iter ] = GD( x, config, thres, alpha )

if nargin == 3
    alpha = 0.001;
end

iter = 0;
grad = ngradient(x, config);
gap = 1;

last_reg = [];
current_reg = [];
current_p = [];
last_p = [];

while gap > thres
    
    vf = func(x, config);
    if ~isreal(vf) || isinf(vf)
        fprintf('not feasible\n');
        return
    end
    
    iter = iter + 1;
    x_new = x - alpha * grad;
    x_new
    [x_new, config] = check(x_new, config);
    grad = ngradient(x_new, config);
    gap = norm(grad);
    gap
    x = x_new;
    
    if (size(config.reg,1) ~= size(current_reg,1)) || (norm(config.reg-current_reg) > 1e-5)
        if (size(config.reg,1) == size(last_reg,1)) && (norm(config.reg-last_reg) < 1e-5)...
                && norm(last_p-x) < thres
            break;
        end
        last_reg = current_reg;
        last_p = current_p;
        current_reg = config.reg;
        current_p = x;
    end
    
%     waitforbuttonpress
end

end

