function [new_n, new_Xi, new_P, new_w] = NURBS_Curve_Refine_oneknot(d, add_Xi, p, n, Xi, P, w)
%Insert one knot into the curve

if (d~=size(P,2))
    error('dimension mismatch between d and P');
end

if(n~=size(P,1) | n~=length(Xi)-p-1)
    error('dimension mismath between n, P and Xi');
end

if(add_Xi<Xi(1) | add_Xi>Xi(length(Xi)))
    error('add xi out of range');
end

new_n = n+1;
new_Xi=zeros(1, length(Xi)+1);
new_P = zeros(size(P,1)+1, size(P,2));
new_w = zeros(1, length(w)+1);

%find the K set for the Xi to be added
K=1;
for dj=1:length(Xi)-1
    if(add_Xi>= Xi(dj) & add_Xi<=Xi(dj+1) & Xi(dj)~=Xi(dj+1))
        K=dj;   
        break;        
    end
end

%creat new knots vector
for ii=1:K
    new_Xi(ii) = Xi(ii);
end
new_Xi(K+1) = add_Xi;
for ii=K+2:length(new_Xi)
    new_Xi(ii) = Xi(ii-1);
end

%create alpha array(transforming array)
alpha = zeros(1, length(new_P));
for jj=1:K-p
    alpha(jj)=1;
end
for jj=K-p+1:K
    alpha(jj)=(add_Xi-Xi(jj))/(Xi(jj+p)-Xi(jj));
end
for jj=K+1:length(new_P)
    alpha(jj)=0;
end

%create a set of control points with higher dimension
P_HD = zeros(size(P,1), size(P,2)+1);
for ii=1:size(P,1)
    for jj=1:size(P,2)
        P_HD(ii,jj) = P(ii,jj)*w(ii);
    end
    P_HD(ii,size(P_HD,2)) = w(ii);
end

%generate the new set of control points with higher dimension
P_HD_new = zeros(size(P_HD));
P_HD_new(1, 1:size(P_HD,2)) = P_HD(1, 1:size(P_HD,2));

for ii=2:new_n-1
    for jj=1:size(P_HD,2)
        P_HD_new(ii,jj) = alpha(ii)*P_HD(ii,jj) + (1-alpha(ii))*P_HD(ii-1,jj);
    end
end
P_HD_new(new_n, 1:size(P_HD,2)) = P_HD(n, 1:size(P_HD,2));

%project to get the actual control points
for ii=1:size(P_HD_new, 1)
    new_w (ii) = P_HD_new(ii, size(P_HD_new, 2));
    for jj=1:d
        new_P(ii,jj)= P_HD_new(ii,jj)/new_w(ii);
    end
end
    

















end

