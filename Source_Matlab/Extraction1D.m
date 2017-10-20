%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Extraction1D
%
% Input:  n = number of functions
%         p = polynomial degree
%         Xi = knot vector
%
% Output: nb = number of 1D elements
%         Ce = array of 1D extraction operators
%
% Purpose: Extract the 1D B-spline basis
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nb, Ce] = Extraction1D(n,p,Xi)

a = p + 1;
b = a + 1;
nb = 1;
Ce(:,:,1) = eye(p+1);

while b < (n+p+1)
    Ce(:,:,nb+1) = eye(p+1);
    i = b;
    
    while ((b < (n+p+1)) && (Xi(b+1) == Xi(b)))
        b = b + 1;
    end
    
    mult = b - i + 1;
    
    if mult < p
        numer = Xi(b) - Xi(a);
        
        for j = p:-1:(mult+1)
            alphas(j-mult) = numer/(Xi(a+j) - Xi(a));
        end
        
        r = p - mult;
        
        for j = 1:r
            save = r-j+1;
            s = mult+j;
            
            for k = (p+1):-1:(s+1)
                alpha = alphas(k-s);
                
                Ce(:,k,nb) = alpha*Ce(:,k,nb) + (1-alpha)*Ce(:,k-1,nb);
            end
            
            if b < (n+p+1)
                Ce(save:(j+save),save,nb + 1) = Ce((p-j+1):(p+1),p+1,nb);
            end
        end
    end
        
    if b < (n+p+1)
        a = b;
        b = a + 1;
        nb = nb + 1;
    end
end