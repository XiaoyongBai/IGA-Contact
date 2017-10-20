%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: compute derivtatives of NURBS curve
%
%Input: xi=position where the derivative was evaluated
%       d_order=order of differentiation
%       p=the polynomial order of B-Spline
%       Xi=Knot Vector
%       P=Coordinates of control points
%       W=weights of control points
%
%Output: Derivative=derivative of the curve at a fixed position
function [ Derivatives ] = NURBS_Curve_derivatives(xi, d_order, p, Xi, P, W)

    [ncp, sd]=size(P);
    
    if ncp~=length(W)
        error('NURBS_Curve_derivatives:dimension mismatch between P and W');
    end
   
    if d_order>=1
        P_w=zeros(ncp, sd);
        for ncp_i=1:ncp
            P_w(ncp_i, :)=P(ncp_i, :)*W(ncp_i); %weighted control points
        end

        A_k=Spline_Curve_derivatives(xi, d_order, p, Xi, P_w); %first term in Eqn(4.8) in nurbs book 2nd edition

        B_k=zeros(1,sd); %second term in Eqn(4.8) in NURBS Book 2nd edition
        
        for d_i=1:d_order
            alpha=factorial(d_order)/(factorial(d_i)*factorial(d_order-d_i));
            w_i=Spline_Curve_derivatives(xi, d_i, p, Xi, W);
            C_ki=NURBS_Curve_derivatives(xi, d_order-d_i, p, Xi, P, W);
            B_k=B_k+alpha*w_i*C_ki;
        end
        
        w_total=0;
        for n_j=1:ncp
                w_total=w_total+NURBS_1D(xi, n_j, p, ncp, Xi, W)*W(n_j);
        end
    
        Derivatives=(A_k-B_k)/w_total;
        
    elseif d_order==0
               
        C=zeros(1, sd);
        for n_j=1:ncp
            for dim_i=1:sd
                C(dim_i)=C(dim_i)+NURBS_1D(xi, n_j, p, ncp, Xi, W)*P(n_j, dim_i);
            end
        end
        
        Derivatives=C;
    end
        





end

