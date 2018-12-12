function slope_coeff = slope(A, B)
% Function for slope calculation. 
% A: 2*1 column vector.  
% B: 2*1 column vector.
% slope_coeff = diff(B(2)-A(2))./diff(B(1)-A(1))
   
% Mengxue Hou. Georgia Institute of Technology. 2018/12
   
if B(1) == A(1)
    slope_coeff = inf;
else
    slope_coeff = (B(2)-A(2))./(B(1)-A(1));
end
end