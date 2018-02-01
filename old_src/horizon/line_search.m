function alpha = line_search( func, vf, x, config, grad )

alpha = 1.0;
c = 1e-7;
rho = 0.7;

while func(x - alpha * grad, config) > vf - c * alpha * (grad' * grad)
    alpha = alpha * rho;
    fprintf('%.8f\t%.8f\t%.8f', func(x - alpha * grad, config), vf - c * alpha * (grad' * grad), alpha);
    fprintf('\n');
end

end

