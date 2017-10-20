function [ N ] = NURBS_1D(xi, i, p, n, Xi, w)
%compute 1d nurbs basis functions

%check to make sure the consistency of polynomial degree and number of
%basis functions

NumKnots = length(Xi);

if(NumKnots ~= n+p+1)
    error('n+p+1 is not equal to the number of knots');
end

if(length(w) ~= n)
    error('dimension mismatching between w and n');
end

if(xi<Xi(1) | xi>Xi(NumKnots))
    N=0;
    return;
elseif(xi==Xi(NumKnots) & i==n)
    N=1;
    return;
elseif(xi==Xi(NumKnots) & i ~= n)
    N=0;
    return;
end

%initialize the basis functions
N_temp = zeros(NumKnots-1, p+1); %All non-zero basis functions are stored in this matrix

%find the L set
L=1;
for dj=1:NumKnots-1
    if(xi>= Xi(dj) & xi<Xi(dj+1) & Xi(dj)~=Xi(dj+1))
        L=dj;   
        break;        
    end
end

N_temp(L,0+1) = 1;

for k=1:p
    for ni=L-k+1:L
        if Xi(ni+k)-Xi(ni)==0
            beta_ki=0;
        else
            beta_ki =(xi-Xi(ni))/(Xi(ni+k)-Xi(ni));
        end
        
        N_temp(ni,k+1) = N_temp(ni, k+1) + beta_ki*N_temp(ni, k);
        N_temp(ni-1, k+1) = N_temp(ni-1,k+1) + (1-beta_ki)*N_temp(ni,k);

    end
end

W=0;
for ni=1:n
    W= W + w(ni)*N_temp(ni,k+1);
end

N= N_temp(i,k+1) * w(i)/W;
end

