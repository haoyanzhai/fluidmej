function [pert_boundary_no] = perturb_region(Intersect_pts, ...
    current_transect_no, left_end, right_end, flag)

% perturb_region: given current transect of boundary, get the boundary of
% the perturbed junction given the flag.

% Input: 

% Intersection_pts:
% Intersection_pts.reg_no(n,:): region number for the boundaries of the nth
% Intersection point. Column vector.
% Intersection_pts.pos(n,:): position of the nth intersection point ([Lon, Lat]).
% --------------------------------
% current_transect_no: transect number of the current junction point, denoted by 
% regions on the two sides of this transect, n_1 n_2 (n_1>n_2). 
% --------------------------------
% left_end: left end position of the transect of current junction. 2D column
% vector, [Lon; Lat].
% right_end: right end point of the transect of current junction. 2D column
% vector, [Lon; Lat].
% --------------------------------
% flag: if flag == 0, then choose the neighboring regions of the left end  
% intersection point. Otherwise, choose the neighboring regions of the right end 
% intersection point.

% Output: 
% pert_boundary_no{ii}:  pert_boundary_no{ii} = [n1, n2], the i^th boundary after perturbation, between 
% region n_1 and region n_2 (n_1 > n_2).


% -------------------------------- 
% Example: 

%              |                    |
%              |        (5)         |
%       (7)    *--------------------*    (8)
%              |        (6)         |
%              |                    |

% Example inputs: 
% current_transect_no = [5,6];
% left_end  = [15.5, 90.5];
% right_end = [50.5, 90.5];
% flag = 0;

% outputs:
% pert_boundary_no{1} = [5,7]; % boundary between region 5 and 7
% pert_boundary_no{2} = [6,7]; % boundary between region 6 and 7
% ---------------------------------

% Mengxue Hou
% Georgia Institute of Technology, 2018/11.


n1 = current_transect_no(1);
n2 = current_transect_no(2);
kk = 1; 
intersect_no = [];

for ii = 1:size(Intersect_pts.reg_no,1)
    if numel(intersect(Intersect_pts.reg_no(ii,:), [n1, n2])) == 2
        intersect_no = [intersect_no, ii];
        kk = kk+1;
    end
end

for kk = 1:numel(intersect_no)
    if flag == 0
        if Intersect_pts.pos(intersect_no(kk),:) == left_end
           pert_no = Intersect_pts.reg_no(intersect_no(kk),:);
        end
    else 
        if flag == 1
            if Intersect_pts.pos(intersect_no(kk),:) == right_end
               pert_no = Intersect_pts.reg_no(intersect_no(kk),:);
            end
        end
    end
end

n3 = setdiff(pert_no, [n1, n2]);

for kk = 1:numel(n3)
    pert_boundary_no{2*kk-1} = sort([n1, n3(kk)],'descend');
    pert_boundary_no{2*kk}   = sort([n2, n3(kk)],'descend');
end

end

