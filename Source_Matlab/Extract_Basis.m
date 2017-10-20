%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Extract_Basis
%
% Input:  p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%
% Output: n_el = number of elements
%         C_operators = array of extraction operators
%         IEN = IEN array
%
% Purpose: Extract the basis information
%
% Notes: Calls Extraction2D and IEN2D from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2)

    [n_el,C_operators] = Extraction2D(n_1,n_2,p_1,p_2,Xi_1,Xi_2);
    [IEN] = IEN2D(n_1,n_2,p_1,p_2,Xi_1,Xi_2);
    
end