%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: IEN1D
%
% Input:  n = number of functions
%         p = polynomial degree
%         Xi = knot vector
%
% Output: IEN = 1D IEN array
%
% Purpose: Compute the 1D IEN array
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IEN] = IEN1D(n,p,Xi)

l = p+1;
e = 1;

while l < (n+1)
    for a = 1:1:(p+1)
        IEN(a,e) = (l+a)-(p+1);
    end
    
    l = l+1;
    
    while ((Xi(l+1) == Xi(l)) && (l < n+1))
        l = l+1;
    end
    
    if l < (n+1)
        e = e+1;
    end
end