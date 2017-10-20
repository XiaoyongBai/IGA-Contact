%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function Shape_Function
%
%Input:  p_1=polynomial degree in direction 1
%        p_2=polynomial degree in direction 2
%        xi_1=shape coordinate of the quadrature point in direction 1
%        xi_2=shape coordinate of the quadrature point in direction 2
%        pbe=coordinates of Bernstein control points
%        wbe=weights of bernstein control points
%        We=weights of NURBS control points
%        Ce=local extraction operator
%
%Output:
%        Ra=vlaues of NURBS basis functions
%        Ra_x=values of NURBS basis derivatives in direction x
%        Ra_y=values of NURBS basis derivatives in direction y
%        x=coordinates of the quadrature point
%        J_det=determinant of Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, pbe, wbe, We, Ce)

    %locate and initialize arrays and matrices
    n_basis=(p_1+1)*(p_2+1);
    Ra=zeros(1, n_basis);
    Ra_x=zeros(1, n_basis);
    Ra_y=zeros(1, n_basis);
    
    Ra_z1=zeros(1,n_basis);
    Ra_z2=zeros(1,n_basis);
    
    dimension=2; %this shape function is designed for 2d problem
    x=zeros(dimension,1);
    J_matrix=zeros(dimension);
    
    w_total=0;
    w_total_z1=0; %derivative of w_total with respect to z1
    w_total_z2=0; %derivative of w_total with respect to z2
    
   
    %compute Bernstein Basis and Derivatives
    [B, d_B_z1, d_B_z2] = BernsteinBasisAndDerives(p_1, p_2, xi_1, xi_2);

    %weighting and derivatives
    for a=1:n_basis
        w_total=w_total+wbe(a)*B(a);
        w_total_z1=w_total_z1+wbe(a)*d_B_z1(a);
        w_total_z2=w_total_z2+wbe(a)*d_B_z2(a);
    end
    
    %Basis functions and Parametric Derivatives
    for a=1:n_basis
        for b=1:n_basis
            Ra(a)=Ra(a)+We(a)*Ce(a,b)*B(b)/w_total;
            
            Ra_z1(a)=Ra_z1(a)+We(a)*Ce(a,b)*(d_B_z1(b)/w_total-w_total_z1*B(b)/w_total^2);
            Ra_z2(a)=Ra_z2(a)+We(a)*Ce(a,b)*(d_B_z2(b)/w_total-w_total_z2*B(b)/w_total^2);
        end
    end
    
    

    %Physical space quantities
    for a=1:n_basis
        x(1)=x(1)+wbe(a)*pbe(a,1)*B(a)/w_total;
        x(2)=x(2)+wbe(a)*pbe(a,2)*B(a)/w_total;
        
        J_matrix(1,1)=J_matrix(1,1)+wbe(a)*pbe(a,1)*(d_B_z1(a)/w_total-w_total_z1*B(a)/w_total^2);
        J_matrix(1,2)=J_matrix(1,2)+wbe(a)*pbe(a,1)*(d_B_z2(a)/w_total-w_total_z2*B(a)/w_total^2);
        J_matrix(2,1)=J_matrix(2,1)+wbe(a)*pbe(a,2)*(d_B_z1(a)/w_total-w_total_z1*B(a)/w_total^2);
        J_matrix(2,2)=J_matrix(2,2)+wbe(a)*pbe(a,2)*(d_B_z2(a)/w_total-w_total_z2*B(a)/w_total^2);
    end
    
    J_det = det(J_matrix); 
    
    if J_det<=0
        error('Shape_function:: J<0');
    end
    
    J_inverse=inv(J_matrix);
    
    for a=1:n_basis
        Ra_x(a)=Ra_z1(a)*J_inverse(1,1)+ Ra_z2(a)*J_inverse(2,1);
        Ra_y(a)=Ra_z1(a)*J_inverse(1,2)+ Ra_z2(a)*J_inverse(2,2);
    end      



end

