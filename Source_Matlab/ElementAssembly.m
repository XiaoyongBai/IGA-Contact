%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function: element Assemble
%
%Input: e=element id
%       K_e=element stiffness matrix
%       F_e=element load vector
%       IEN=global ids of element points
%       BC=Dirichlet boundary array
%       g=Dirichlet boundary
%       K_global=global stiffness matrix
%       F_global=global load vector
%
%Output: K_global_updated=updated global stiffness matrix
%        F_global_updated=updated global load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K_global_updated, F_global_updated ] = ElementAssembly(e, K_e, F_e, IEN, BC, g, K_global, F_global)

    %initialize k_global_updated and F_global_updated
    K_global_updated=K_global;
    F_global_updated=F_global;

    sd=2;
    n_cp=size(IEN,1);%number of control point per element
    
    for A=1:sd
        for a=1:n_cp
            p_loc=sd*(a-1)+A;
            P_global=LM(A,a,e,IEN);
            i=IEN(a,e);
            
            if BC(A,i)==0
                
                for B=1:sd
                    for b=1:n_cp
                        q_loc=sd*(b-1)+B;
                        Q_global=LM(B,b,e,IEN);
                        j=IEN(b,e);
                        
                        if BC(B,j)==0
                            K_global_updated(P_global,Q_global)=K_global_updated(P_global,Q_global)+K_e(p_loc,q_loc);
                        else
                            F_global_updated(P_global)=F_global_updated(P_global)-K_e(p_loc,q_loc)*g(B,j);
                        end%if BC
                    end%forb
                end%for B
                
                F_global_updated(P_global)=F_global_updated(P_global)+F_e(p_loc);
                
            end%if BC
            
        end%for a
    end%for A
    
    


end

