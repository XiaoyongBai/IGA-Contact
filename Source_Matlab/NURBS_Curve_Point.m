function [ Point] = NURBS_Curve_Point(xi, d, p, n, Xi, P, w )
%plot NURBS curve

    [ncp, sd]=size(P);
    
    if(d ~=2 & d~=3)
        error('NURBS_Curve:: unsupported dimension');
    end

    if(n~=ncp)
        error('dimension mismatching between n and P (basis function and control knots)');
    end

    if(n~=length(w))
        error('dimension mismathcing between n and w(basis functions and weights)');
    end

    if(length(Xi) ~= n+p+1)
        error('dimension mismatching between n+p+1 and Xi(basis functions, its degree and knots)');
    end

    Point=zeros(1, sd);
    
    for ni=1:n
        N = NURBS_1D(xi, ni, p, n, Xi, w);
        for si=1:sd
             Point(si)=Point(si) + N*P(ni,si);
        end
    end
   


end

