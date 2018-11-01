function boundary = find_boundary(boundaries, n1, n2)

% Get boundary position between two regions.

% Example: 
% boundary = find_boundary(boundaries, num2str(2), num2str(1));

% Input: 
% boundaries: structure including all boundary information
% n_1/n_2: region no. (n_1>n_2)

% Output: 
% boundary: boundary position between region n_1 and n_2 (N*2 matrix)
% boundary(:,1): Lon 
% boundary(:,2): Lat

% Mengxue Hou, Georgia Institute of Technology
% 10/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

b_name = ['boundary_no_', num2str(n1), num2str(n2)];
boundary = getfield(boundaries, b_name);

end