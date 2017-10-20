%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Extraction2D
%
% Input:  n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%
% Output: n_el = number of 2D elements
%         C_e = array of 2D extraction operators
%
% Purpose: Extract the 2D B-spline basis
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_el, C_e] = Extraction2D(n_1,n_2,p_1,p_2,Xi_1,Xi_2)

[n_el_1, C_e_1] = Extraction1D(n_1,p_1,Xi_1);
[n_el_2, C_e_2] = Extraction1D(n_2,p_2,Xi_2);

n_el = n_el_1*n_el_2;

for e1 = 1:1:n_el_1
    for e2 = 1:1:n_el_2
        e = (e1-1)*n_el_2+e2;
        C_e(:,:,e) = kron(C_e_1(:,:,e1),C_e_2(:,:,e2));
    end
end