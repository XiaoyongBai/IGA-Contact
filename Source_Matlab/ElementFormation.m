%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: ElementFormation
%
%Input:  p_1=polynomial degree in direction 1
%        p_2=polynomial degree in direction 2
%        xi_q=shape coordinate of the quadrature point
%        w_q=qudrature weights
%        C_operator=Extraction matrix for this element
%        Pb=coordinates of Berinstein control points
%        Wb=weights for Bernstein control points
%        We=weights of NURBS control points
%        Newmann_e=Newmann boundary condition for current element
%        Newmann_h_loc_e=stress on newmann boundary of the element
%
%Output:
%        K_e=element stiffness matrix
%        F_e=load vector of the element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K_e, F_e] = ElementFormation(p_1, p_2, xi_q, w_q, Ce, pbe, wbe, We, Newmann_e)

    sd=2;%spatial dimension
    ncp=length(wbe); %number of control points in this element
    
    K_e=zeros(2*ncp);  %stiffness matrix
    F_e=zeros(2*ncp,1); %load vector
    
    D_matrix=D_local();
    
    n_q=length(w_q);%number of qudrature points
    for qi=1:n_q
        for qj=1:n_q
            xi_1=xi_q(qi);
            xi_2=xi_q(qj);

            [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, pbe, wbe, We, Ce);

            %stiffness matrix and load vector
            for a=1:ncp
                B_a=[Ra_x(a), 0; 0, Ra_y(a); Ra_y(a), Ra_x(a)]; 
                for b=1:ncp
                    B_b=[Ra_x(b), 0; 0, Ra_y(b); Ra_y(b), Ra_x(b)];
                    K_loc=transpose(B_a)*D_matrix*B_b*w_q(qi)*w_q(qj)*J_det;
                    for A=1:sd
                        for B=1:sd
                            p=sd*(a-1)+A;
                            q=sd*(b-1)+B;
                            K_e(p,q)=K_e(p,q)+K_loc(A,B);
                        end%for B
                    end%for A
                    
                end %for b
            end%for a

        end %for qj                    
    end%for qi

    for A=1:sd
        for side=1:4
            if Newmann_e(A, side)==1
               
                for q_i=1:length(xi_q)
                    switch side
                        case 1
                            xi_1=xi_q(q_i);
                            xi_2=0;
                            tangent=[1;0];
                        case 2
                            xi_1=1;
                            xi_2=xi_q(q_i);
                            tangent=[0;1];
                        case 3
                            xi_1=xi_q(q_i);
                            xi_2=1;
                            tangent=[1;0];
                        case 4
                            xi_1=0;
                            xi_2=xi_q(q_i);
                            tangent=[0;1];
                    end%switch

                    [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, pbe, wbe, We, Ce);
                    J_side=J_matrix*tangent;
                    J_side_det=norm(J_side);
                
                   h_loc=H_loc(x); %get local heat flux on the boundary
                    for a=1:ncp
                        p_loc=sd*(a-1)+A;
                        F_e(p_loc)=F_e(p_loc)+Ra(a)*h_loc(A)*w_q(q_i)*J_side_det;
                    end %for a (nodes)                        

                
                end%quadrature point
                
            end%if newmann
            
        end%for side
    end%for spatial dimension
                        



end

