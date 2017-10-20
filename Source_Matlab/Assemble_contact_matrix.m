%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: convert segment matrix to global formation
%
%Input: K_contact= original K_contact
%       KC_segment= contact matrix for segement
%       n_s=number of control points in slave segment
%       n_m=number of control points in master segment
%       Global_ID_s= global ids of slave segment
%       Global_ID_m= global ids of master segment
%
%Output: K_new:the new contact matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ K_new ] = Assemble_contact_matrix(K_Contact, KC_segment, n_s, n_m, Global_ID_s, Global_ID_m)

    K_new=K_Contact;
    
    SD=2;
   %Assemble segment contact matrix into global matrix
    for d_j=1:SD
        for ns_i=1:n_s
            A=ID(d_j, ns_i);
            P=ID(d_j, Global_ID_s(ns_i));

            for ns_j=1:n_s
                B=ID(d_j, ns_j);
                Q=ID(d_j, Global_ID_s(ns_j));
                K_new(P,Q)=K_new(P,Q)+KC_segment(A,B);
            end
            for nm_j=1:n_m
                B=SD*n_s+ID(d_j, nm_j);
                Q=ID(d_j, Global_ID_m(nm_j));
                K_new(P,Q)=K_new(P,Q)+KC_segment(A,B);
            end
        end

        for nm_i=1:n_m
            A=SD*n_s+ID(d_j, nm_i);
            P=ID(d_j, Global_ID_m(nm_i));

            for ns_j=1:n_s
                B=ID(d_j, ns_j);
                Q=ID(d_j, Global_ID_s(ns_j));
                K_new(P,Q)=K_new(P,Q)+KC_segment(A,B);
            end
            for nm_j=1:n_m
                B=SD*n_s+ID(d_j, nm_j);
                Q=ID(d_j, Global_ID_m(nm_j));
                K_new(P,Q)=K_new(P,Q)+KC_segment(A,B);
            end
        end
    end%end for d_j


end

