function [ N ] = NURBS_2D(xi_1, xi_2, i_1, i_2, p_1, p_2, n_1, n_2, Xi_1, Xi_2, w)
%compute 2d nurbs basis function

NumKnots_1 = length(Xi_1);
NumKnots_2 = length(Xi_2);

if(NumKnots_1 ~= n_1+p_1+1)
    error('n+p+1 is not equal to the number of knots, for dimension 1');
end

if(NumKnots_2 ~= n_2+p_2+1)
    error('n+p+1 is not equal to the number of knots, for dimension 2');
end

if(length(w) ~= n_1*n_2)
    error('dimension mismatching between w and n');
end

w1 = zeros(n_1, 1)+1;
w2 = zeros(n_2, 1)+1;

%initialize the basis functions
N_temp = zeros(n_1, n_2); %All non-zero basis functions are stored in this matrix

W=0;
for i1=1:n_1
    for i2=1:n_2
        N1 = NURBS_1D(xi_1, i1, p_1, n_1, Xi_1, w1);
        N2 = NURBS_1D(xi_2, i2, p_2, n_2, Xi_2, w2);
        N_temp(i1, i2) = N1*N2;
        W = W + N_temp(i1,i2)*w(i1+n_1*(i2-1));
    end
end

if(W==0)
    N=0;
    return;
else
    N = N_temp(i_1,i_2) * w(i_1 + n_1*(i_2-1))/W;
end


end

