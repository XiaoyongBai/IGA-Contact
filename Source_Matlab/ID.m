%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: return the global equation number corresponding to a dof of a node
%          Option 1 is used: group unknowns by basis function number
%
%Input: A=degree of freedom
%       j=id of the node
%
%Output: P=euqation number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = ID(A,j )

  P=2*(j-1)+A;

end

