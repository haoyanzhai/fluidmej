function grad = ngradient( p, config )

e = 1e-5;
vf = func(p, config);
grad = zeros(size(p));
for i = 1 : length(p)
    tx = p;
    tx(i) = tx(i) + e;
    grad(i) = (func(tx,config) - vf) / e;
end

end

