function [new_n, new_Xi, new_P, new_w] = NURBS_Curve_Refine(d, add_Xi, p, n, Xi, P, w)
%insert multiple knots

for ii=1:length(add_Xi)
    aXi = add_Xi(ii);
    [new_n, new_Xi, new_P, new_w] = NURBS_Curve_Refine_oneknot(d, aXi, p, n, Xi, P, w);
    n=new_n;
    Xi=new_Xi;
    P=new_P;
    w=new_w;
end


end

