%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: return the global equation number corresponding to a dof of a
%          node in an element
%          Option 1 is used: group unknowns by basis function number
%
%Input: A=degree of freedom
%       a=local id of the node in the element
%       e=id of the element
%       IEN=IEN array
%
%Output: P=euqation number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = LM(A, a, e, IEN)

  j=IEN(a,e);
  P=2*(j-1)+A;

end
