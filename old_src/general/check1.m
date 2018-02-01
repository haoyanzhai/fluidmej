function [ end_loop, x, config ] = check1( x, config )

end_loop = false;

if x(1) >= 0.4 && length(x) == 2
    end_loop = true;
    x = x(1);
    config.y = {@(x) 0 * x + 0;
                @(x) x - x + 2; 
                @(x) x - x + 3};
    config.ux = [-1; 0];
    config.uy = [0; 0];
    return;
end

if x(1) < 0.4 && length(x) == 1
    end_loop = true;
    x = [x(1); 0.4];
    config.y = {@(x) 0 * x + 0;
                @fexp;
                @(x) x - x + 2; 
                @(x) x - x + 3};
    config.ux = [-1; -0.9; 0];
    config.uy = [0; 0; 0];
    return;
end

end

