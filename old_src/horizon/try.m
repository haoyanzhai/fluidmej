f1 = @(x) 2*x;
f2 = @(x) x^2;

f = {@(x) 2*x; @(x) x^2};

x = [1;2];

y = [];
for i = 1 : length(f)
    y = [y; f{i}(x(i))];
end
y

