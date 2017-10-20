%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: incorporate contact constraint to stiffness matrix based on
%          penaly method
%
%Input: ncp_1=number of control point in patch 1
%       ncp_2=number of control point in patch 2
%       Problem=Problem set
%       Penalty=Penalty parameter
%       dis_1=displacement of patch 1
%       dis_2=displacement of patch 2
%
%Output:
%       K_Contact=Contact Stiffness matrix (finite element + contact)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K_Contact] = Contact_Formualtion(ncp_1, ncp_2, Problem, Penalty, dis_1, dis_2)
     
    K_Contact=zeros(2*(ncp_1+ncp_2)); %initialize the Output;

    SD=2;%spatial dimension
    
    %get geometry of slave segment
    [p_s, n_s, Xi_s, P_s, W_s, nq_s, patch_s, Global_ID_s, normal_type_s]=  Contact_Geometry(Problem, 0, dis_1, dis_2);

    %get geometry of master segment
    [p_m, n_m, Xi_m, P_m, W_m, nq_m, patch_m, Global_ID_m, normal_type_m]=  Contact_Geometry(Problem, 1, dis_1, dis_2);
    
    [IEN_s] = IEN1D(n_s,p_s,Xi_s);
    [IEN_m] = IEN1D(n_m,p_m,Xi_m);
    
    [nnd_s, ne_s]=size(IEN_s); %number of nodes per element, and number of elements in slave segment
    [nnd_m, ne_m]=size(IEN_m); %number of nodes per element, and number of elements in slave segment

    %get unique knot vetor, to compute the global natural coordinate from
    %                       local natural coordinate
    Xi_s_unique=unique(Xi_s);
    Xi_m_unique=unique(Xi_m);
    
    switch nq_s
        case 3
            xi_q=[1/2-sqrt(3/5)/2, 1/2, 1/2+sqrt(3/5)/2];
            w_q=[5/18, 4/9, 5/18];
        case 4
            xi_q=[(1-sqrt(3/7-2/7*sqrt(6/5)))/2,  (1+sqrt(3/7-2/7*sqrt(6/5)))/2,  (1-sqrt(3/7+2/7*sqrt(6/5)))/2,  (1+sqrt(3/7+2/7*sqrt(6/5)))/2];
            w_q=[(18+sqrt(30))/36, (18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36]/2;    
        otherwise
            error('unsupported number of quadrature');
    end
  
    %form contant stiffness matrix 
    KC_segment=zeros(SD*(n_s+n_m));%initialize contact matrix for segment
    
    for e_i=1:ne_s %loop over element in slave segment

        for q_i=1:nq_s
            xi_slave=(1-xi_q(q_i))*Xi_s_unique(e_i)+xi_q(q_i)*Xi_s_unique(e_i+1);
            x_slave=NURBS_Curve_Point(xi_slave, SD, p_s, n_s, Xi_s, P_s, W_s); %coordinate of slave point 
            
            [xi_master, x_master, active]=Closet_Point_Curve(x_slave, p_m, Xi_m, P_m, W_m, normal_type_m);%find closest projection for slave point
            
            if active==1
                Ra_s=zeros(SD,SD*nnd_s);
                Ra_m=zeros(SD,SD*nnd_m);
                for ns_i=1:nnd_s
                    node=IEN_s(ns_i, e_i);
                    Ra_s(1, SD*(ns_i-1)+1)=NURBS_1D(xi_slave, node, p_s, n_s, Xi_s, W_s);
                    Ra_s(2, SD*(ns_i-1)+2)=Ra_s(1, SD*(ns_i-1)+1);
                end
                
                %find which element the master point belongs
                m_e=1;
                for j=1:ne_m
                    if xi_master>=Xi_m_unique(j) & xi_master<=Xi_m_unique(j+1)
                        m_e=j;
                        break;
                    end
                end
                
                %compute basis functions
                for nm_i=1:nnd_m
                    node=IEN_m(nm_i, m_e);
                    Ra_m(1, SD*(nm_i-1)+1)=NURBS_1D(xi_master, node, p_m, n_m, Xi_m, W_m);
                    Ra_m(2, SD*(nm_i-1)+2)=Ra_m(1, SD*(nm_i-1)+1);
                end

                %group contact stiffness matrix for the point
                KC_cp=Penalty*[Ra_s'*Ra_s, -Ra_s'*Ra_m; -Ra_m'*Ra_s, Ra_m'*Ra_m];
                                
                %Assemeble point matrix into segment matrix
                for d_j=1:SD
                    for ns_i=1:nnd_s
                        A=ID(d_j, ns_i);
                        P=ID(d_j, IEN_s(ns_i, e_i));
                        
                        for ns_j=1:nnd_s
                            B=ID(d_j, ns_j);
                            Q=ID(d_j, IEN_s(ns_j, e_i));
                            KC_segment(P,Q)=KC_segment(P,Q)+KC_cp(A,B);
                        end
                        for nm_j=1:nnd_m
                            B=SD*nnd_s+ID(d_j, nm_j);
                            Q=SD*n_s+ID(d_j, IEN_m(nm_j, m_e));
                            KC_segment(P,Q)=KC_segment(P,Q)+KC_cp(A,B);
                        end
                    end
                    
                    for nm_i=1:nnd_m
                        A=SD*nnd_s+ID(d_j, nm_i);
                        P=SD*n_s+ID(d_j, IEN_m(nm_i, m_e));
                        
                        for ns_j=1:nnd_s
                            B=ID(d_j, ns_j);
                            Q=ID(d_j, IEN_s(ns_j, e_i));
                            KC_segment(P,Q)=KC_segment(P,Q)+KC_cp(A,B);
                        end
                        for nm_j=1:nnd_m
                            B=SD*nnd_s+ID(d_j, nm_j);
                            Q=SD*n_s+ID(d_j, IEN_m(nm_j, m_e));
                            KC_segment(P,Q)=KC_segment(P,Q)+KC_cp(A,B);
                        end
                    end
                end%end for d_j
                
            end
        end
    end %end loop over e_i
        
    K_Contact=Assemble_contact_matrix(K_Contact, KC_segment, n_s, n_m, Global_ID_s, Global_ID_m);   

 end

