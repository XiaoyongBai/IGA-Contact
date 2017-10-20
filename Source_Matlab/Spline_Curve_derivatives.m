%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: compute derivtatives of Spline curve
%
%Input: xi=position where the derivative was evaluated
%       d_order=order of differentiation
%       p=the polynomial order of B-Spline
%       Xi=Knot Vector
%       P=Coordinates of control points
%
%Output: Derivative=derivative of the curve at a fixed position
function [ Derivatives ] = Spline_Curve_derivatives(xi, d_order, p, Xi, P)

    if xi>max(Xi) | xi<min(Xi)
        error('BSpline_Basis_derivatives:xi is out of range');
    end
    [ncp, sd]=size(P);
    Derivatives=zeros(1, sd); %derivative of anti-perpective B-Spline curve
    
    for n_i=1:ncp
        Derivatives(1:sd)= Derivatives(1:sd)+P(n_i, 1:sd)*BSpline_Basis_derivatives(xi, d_order, n_i, p, Xi);
    end
    
    




end

