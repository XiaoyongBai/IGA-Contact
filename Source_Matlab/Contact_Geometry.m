%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: specify contact geometry entities
%
%input: Problem=0:test problem of two blocks
%               1:Hertz problem
%       Slave_or_Master= 0:slave segment
%                        1:master segment
%
%output: p,n,Xi,P,W=NURBS parameters for the output segment(in 2D, it's a curve)
%        n_q=number of integral points within one element segment
%        patch=the patch number that the segement attached to
%        Global_ID=global ids of the control points on the curve
%        normal_type= 1:normal is 90 anti-clockise to tangent of the curve
%                     2:normal is 90 clockwise to tangent of the curve
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p, n, Xi, P, W, n_q, patch, Global_ID, normal_type] = Contact_Geometry(Problem, Slave_or_Master, dis_1, dis_2)
    
    switch Problem
        case 0
            switch Slave_or_Master
                %in problem 0, slave segment is the upper smaller block
                %             which patch=2
                case 0
                    patch=2;
                    p=2;
                    n=5;
                    Xi=[0,0,0,1,2,3,3,3];
                    
                    [~, ~, ~, ~, ~, ~, P_global, ~, ~] = Geometry(Problem, patch, dis_1, dis_2);
                    P=P_global(1:5, 1:2);

                    W=zeros(n,1)+1;
                    Global_ID=[1,2,3,4,5]+25;
                    n_q=p+1;
                    normal_type=2;
                case 1
                    patch=1;
                    p=2;
                    n=5;
                    Xi=[0,0,0,1,2,3,3,3];
                    
                    [~, ~, ~, ~, ~, ~, P_global, ~, ~] = Geometry(Problem, patch, dis_1, dis_2);
                    P=P_global(21:25, 1:2);
                    
                    W=zeros(n,1)+1;
                    Global_ID=[21,22,23,24,25];
                    n_q=p+1;
                    normal_type=1;
                otherwise
                    error('Contact_Geometry:slave_or_master can only be chose between 0 or 1');
            end
            
        case 1%pipe contact plate problem
            switch Slave_or_Master
                %in problem 0, slave segment is the upper smaller block
                %             which patch=2
                case 0
                    patch=2;
                    [p, ~, n, ~, Xi, ~, P_global, w, n_q] = Geometry(Problem, patch, dis_1, dis_2);
                    P=P_global(1:n, 1:2);

                    W=w(1:n);
                    Global_ID=1:n;
                    Global_ID=Global_ID+25;
                    normal_type=2;
                case 1
                    patch=1;
                    p=2;
                    n=5;
                    Xi=[0,0,0,1,2,3,3,3];
                    
                    [~, ~, ~, ~, ~, ~, P_global, ~, ~] = Geometry(Problem, patch, dis_1, dis_2);
                    P=P_global(21:25, 1:2);
                    
                    W=zeros(n,1)+1;
                    Global_ID=[21,22,23,24,25];
                    n_q=p+1;
                    normal_type=1;
                otherwise
                    error('Contact_Geometry:slave_or_master can only be chose between 0 or 1');
            end
        otherwise
            error('Contact_Geometry: Unsupported Problem');
    end
          

end

