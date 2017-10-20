%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: compute derivatives of BSpline basis functions
%
%input: xi=position where the derivative was evaluated
%       d_order=order of differentiation
%       n_i=the basis function that was differentiated
%       p=the polynomial order of B-Spline
%       Xi=Knot Vector
%
%output: Derivative=derivative of the function at a fixed position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Derivative] = BSpline_Basis_derivatives(xi, d_order, n_i, p, Xi)
    
    if xi>max(Xi) | xi<min(Xi)
        error('BSpline_Basis_derivatives:xi is out of range');
    end
    
    if(n_i<0 | n_i>length(Xi)-p-1)
        error('BSpline_Basis_derivatives:n_i is out of range');
    end

    if d_order>p
        Derivative=0;
        return;
    end
    
    new_p=p-1;

    if Xi(n_i+p)-Xi(n_i)==0
        alpha=0;
    else
        alpha=p/(Xi(n_i+p)-Xi(n_i));
    end

    if Xi(n_i+p+1)-Xi(n_i+1)==0
        beta=0;
    else
        beta=p/(Xi(n_i+p+1)-Xi(n_i+1));
    end
        
    if d_order<=p & d_order>1
        Ni=BSpline_Basis_derivatives(xi, d_order-1, n_i, new_p, Xi);
        Nip1=BSpline_Basis_derivatives(xi, d_order-1, n_i+1, new_p, Xi);       
    elseif d_order==1
        n=length(Xi)-p-1; %number of basis functions
        new_n=n+1;
        w=zeros(new_n,1)+1;
        Ni=NURBS_1D(xi, n_i, new_p, new_n, Xi, w);
        Nip1=NURBS_1D(xi, n_i+1, p-1, new_n, Xi, w);       
    else
        error('BSpline_Basis_derivatives:unsupported derivative order');        
    end

    Derivative=alpha*Ni - beta*Nip1;

end

