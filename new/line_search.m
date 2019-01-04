function alpha = line_search( p, grad, config )

alpha = 0.9;
vf = func(p, config);
vf
c = 0.001;

while ~isreal(func(p - alpha*grad, config)) && alpha > 1e-50
    alpha = 0.8 * alpha;
end

while func(p - alpha*grad, config) - vf + c*alpha*(grad'*grad) > 0 && alpha > 1e-50
    alpha = 0.8 * alpha;
end

end

