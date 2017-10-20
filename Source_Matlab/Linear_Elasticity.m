%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function: solved static linear elastic problems
%
%Input: p_1=polynomial degree in direction 1
%       p_2=polynomial degree in direction 2
%       n_1=number of basis functions in direction 1
%       n_2=number of basis functions in direction 2
%       Xi_1=knot vector in direction 1
%       Xi_2=knot vector in direction 2
%       P=coordinates of control points
%       w=weight of control points
%       n_q=number of integration point per element
%       problem=intput set
%       Patch=the Patch that are computed
%
%Output: K_patch=stiffness matrix of one patch
%        F_patch=load vector of one patch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K_patch, F_patch] = Linear_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q, problem, patch )

    %Extraction
    [n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2);

    sd =2; %sptial dimension

    [P_b,w_b] = Extract_Geometry(sd,p_1,p_2,n_el,C_operators,IEN,P,w);

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

    %read boundary conditions
    [BC1, g1,  Newmann] = Boundary_Condition(problem, patch, p_1, p_2, n_1, n_2, Xi_1, Xi_2);
   
    K_patch=zeros(2*length(P));
    F_patch=zeros(2*length(P),1);
    
    
    for e_id=1:n_el %loop over element
        C_e=C_operators(:,:,e_id);
        Pb_e=P_b(:,:,e_id);
        wb_e=w_b(:,e_id);
        We_e=We(:,e_id);
        Newmann_e=Newmann(:,:, e_id);

        [K_e, F_e] = ElementFormation(p_1,p_2,xi_q, w_q, C_e, Pb_e, wb_e, We_e,  Newmann_e);

        [K_patch, F_patch]=ElementAssembly(e_id, K_e, F_e, IEN, BC1, g1, K_patch, F_patch);
    end
    
    %Account for BCs in stiffness matrix
    for ni=1:length(P)
        for A=1:sd
            if BC1(A, ni)==1
                P_global=sd*(ni-1)+A;
                K_patch(P_global,P_global)=1;
                F_patch(P_global)=g1(P_global);
            end
        end
    end    
    

end

