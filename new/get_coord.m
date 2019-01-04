function coord = get_coord( lambda, x )

coordx = (1-lambda) .* x(:,1) + lambda .* x(:,3);
coordy = (1-lambda) .* x(:,2) + lambda .* x(:,4);
coord = [coordx coordy];

end

