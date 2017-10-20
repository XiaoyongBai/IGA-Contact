%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function: compute the L2-norm of the error
%
%input: Problem=id of the problem
%       d=solution of the problem
%
%Output: normal=L2-norm of the error
%        energy_diff=totoal energy of the system
%        ne_l=number of elements in the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [normal, engergy_diff, n_el] = L2_norm(Problem, p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, n_q)
    
    %spatial Dimension
    sd=2;
    pi=30e6;
    ri=0.075;
    ro=0.09;
    
    global nu;
    global Ey;
    D_matrix=D_local();
    
    %initialize norm;
    normal=[0 0];
    strain_energy=0;
    external_work=0;
    
    % Extract the basis, geometry, and temperature field
    [n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2);
    [P_b,w_b] = Extract_Geometry(sd,p_1,p_2,n_el,C_operators,IEN,P,w);
    [BC1, g1,  Newmann] = Boundary_Condition(Problem, p_1, p_2, n_1, n_2, Xi_1, Xi_2);

    %set NURBS weight for element control points
    ncp_ele=(p_1+1)*(p_2+1); %number of control points in one element
    We=zeros(ncp_ele, n_el);

    for ee=1:n_el
        for cpi=1:ncp_ele
            cp=IEN(cpi,ee);
            We(cpi,ee)=w(cp);
        end
    end

    %set qudrature point and weights
    switch n_q
        case 3
            xi_q=[1/2-sqrt(3/5)/2, 1/2, 1/2+sqrt(3/5)/2];
            w_q=[5/18, 4/9, 5/18];
        case 4
            xi_q=[(1-sqrt(3/7-2/7*sqrt(6/5)))/2,  (1+sqrt(3/7-2/7*sqrt(6/5)))/2,  (1-sqrt(3/7+2/7*sqrt(6/5)))/2,  (1+sqrt(3/7+2/7*sqrt(6/5)))/2];
            w_q=[(18+sqrt(30))/36, (18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36]/2;    
        otherwise
            error('unsupported number of quadrature');
    end

    n_el_1 = length(unique(Xi_1))-1;
    n_el_2 = length(unique(Xi_2))-1;
    
    if (Problem==32 | Problem==33)
        for e_2=1:n_el_2 %loop over element       
            for e_1=1:n_el_1
                e_id=(e_2-1)*n_el_1+e_1;

                C_e=C_operators(:,:,e_id);
                Pb_e=P_b(:,:,e_id);
                wb_e=w_b(:,e_id);
                We_e=We(:,e_id);

                %loop over qudarture points
                for xi=1:n_q
                    for yi=1:n_q
                        xi_1=xi_q(xi);
                        xi_2=xi_q(yi);

                        [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, Pb_e, wb_e, We_e, C_e);

                        dx_numerical=0;
                        dy_numerical=0;

                        %stiffness matrix and load vector
                        for a=1:ncp_ele
                            n_id=IEN(a,e_id);
                            d_node=[d(2*(n_id-1)+1); d(2*n_id)];
                            dx_numerical=dx_numerical+Ra(a)*d_node(1);
                            dy_numerical=dy_numerical+Ra(a)*d_node(2);
                         end%for a
                         
                         r=sqrt(x(1)^2+x(2)^2);
                         dr_theoretical=pi*((1+nu)/Ey)*(r/((ro/ri)^2-1))*(1-2*nu+(ro/r)^2);
                         dx_theoretical=dr_theoretical*x(1)/r;
                         dy_theoretical=dr_theoretical*x(2)/r;
                         
                         normal(1)=normal(1)+(dx_numerical-dx_theoretical)^2*J_det*w_q(xi)*w_q(yi);
                         normal(2)=normal(2)+(dy_numerical-dy_theoretical)^2*J_det*w_q(xi)*w_q(yi);                         
                    end %for yi                    
                 end%for xi

            end
         end
        
        
    elseif (Problem==2 || Problem==3)
        for e_2=1:n_el_2 %loop over element       
            for e_1=1:n_el_1
                e_id=(e_2-1)*n_el_1+e_1;

                C_e=C_operators(:,:,e_id);
                Pb_e=P_b(:,:,e_id);
                wb_e=w_b(:,e_id);
                We_e=We(:,e_id);

                %compute strain energy
                %loop over qudarture points
                for xi=1:n_q
                    for yi=1:n_q
                        xi_1=xi_q(xi);
                        xi_2=xi_q(yi);

                        [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, Pb_e, wb_e, We_e, C_e);

                        strain_loc=zeros(sd+1,1);
                        %stiffness matrix and load vector
                        for a=1:ncp_ele
                            n_id=IEN(a,e_id);
                            d_node=[d(2*(n_id-1)+1); d(2*n_id)];

                             B_a=[Ra_x(a), 0; 0, Ra_y(a); Ra_y(a), Ra_x(a)];
                             strain_loc=strain_loc+B_a*d_node;
                         end%for a

                         stress_loc=D_matrix*strain_loc;
                         strain_energy=strain_energy+0.5*stress_loc'*strain_loc*J_det*w_q(xi)*w_q(yi);
                         
                    end %for yi                    
                 end%for xi

                 %compute external work
                 Newmann_e=Newmann(:,:, e_id);
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

                                [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, Pb_e, wb_e, We_e, C_e);
                                J_side=J_matrix*tangent;
                                J_side_det=norm(J_side);
                                  
                                h_loc=H_loc(x); %get local heat flux on the boundary
                                
                                dx=[0,0];%displacement at Gaussian points
                                for a=1:ncp_ele
                                    n_id=IEN(a,e_id);
                                    d_node=[d(2*(n_id-1)+1); d(2*n_id)];
                                    dx(1)=dx(1)+Ra(a)*d_node(1);
                                    dx(2)=dx(2)+Ra(a)*d_node(2);
                                    
                                    external_work=external_work+Ra(a)*h_loc(A)*dx(A)*w_q(q_i)*J_side_det;
                                end %for a (nodes)                        


                            end%quadrature point

                        end%if newmann

                    end%for side
                end%for spatial dimension
                 
            end
         end
    else
        error('L2_norm::unsupported Problem');
    end


    normal=sqrt(normal);

    engergy_diff=strain_energy-external_work;
    
    
end