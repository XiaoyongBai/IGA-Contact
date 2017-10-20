%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function: compute closet point in the master segment for a slave point
%
%Input: Slave_P=coordinates of the slave point
%       p=polynomial degree
%       Xi=Control points
%       P=coordinates of control points
%       W=weights of control points
%        normal_type= 1:normal is 90 anti-clockise to tangent of the curve
%                     2:normal is 90 clockwise to tangent of the curve
%
%Output: Closet_Point=parameter coordinate of the closet point
%        Master_P=coordinate of the closet point
%        active=0:inactive, 1:contact pair is active
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Closet_Point, Master_P, active] = Closet_Point_Curve(Slave_P, p, Xi, P, W, normal_type)

    Closet_Point=(max(Xi)+min(Xi))/2; %initial point was set to the middle of the curve
    
    criterial=1e-10;
    iteration=0;
    
    [n, sd]=size(P);
    
    C_1st=NURBS_Curve_derivatives(Closet_Point, 1, p, Xi, P, W); %tangent
    Master_P=NURBS_Curve_Point(Closet_Point, sd, p, n, Xi, P, W ); %master point

    err=dot(Slave_P-Master_P, C_1st);
    
    while abs(err)>criterial & iteration<10
        C_2nd=NURBS_Curve_derivatives(Closet_Point, 2, p, Xi, P, W);
        Dx=-err/(dot(Slave_P-Master_P, C_2nd)-dot(C_1st, C_1st));
        
        Closet_Point=Closet_Point+Dx;
        
        if Closet_Point>max(Xi)
            Closet_Point=max(Xi)-criterial;
        end
        
        if Closet_Point<min(Xi)
            Closet_Point=min(Xi)+criterial;
        end
        
        C_1st=NURBS_Curve_derivatives(Closet_Point, 1, p, Xi, P, W); %tangent
        Master_P=NURBS_Curve_Point(Closet_Point, sd, p, n, Xi, P, W ); %master point
        
        err=dot(Slave_P-Master_P, C_1st);
        
        iteration=iteration+1;
    end
    
        C_1st=NURBS_Curve_derivatives(Closet_Point, 1, p, Xi, P, W); %tangent
        switch normal_type
            case 1
                normal=[-C_1st(2), C_1st(1)];
            case 2
                normal=[-C_1st(2), C_1st(1)];
            otherwise
                error('Closet_Point_Curve:unsupported normal type');
        end
        
        if dot(Slave_P-Master_P, normal)<=0;
            active=1;
        else 
            active=0;
        end

end

