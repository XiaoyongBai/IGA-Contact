function [new_n_1,new_n_2,new_Xi_1,new_Xi_2,new_P,new_w] =...
    NURBS_Surface_Refine(d,add_Xi_1,add_Xi_2,p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w)

% Determine what new_Xi should look like:
new_Xi_1 = sort([Xi_1, add_Xi_1]);
new_Xi_2 = sort([Xi_2, add_Xi_2]);

% Do all the calculation in projective space: That is multiply the control
% points by their weights:
for i = 1:d
    P(:,:,i) = P(:,:,i).*w;
end
P(:,:,d+1) = w;

% Knot Interstion:
% We know what we want our final knot vector to look like, so we'll do knot
% insertion until our KV matches our KVF.
for i = 1:length(new_Xi_1)
    % Check to see if the current knot in KV matches the current knot in
    % KVF. If it does, great, move on to the next knot, if not, calculate
    % new control points and insert KVF(i) into KV(i).
    if new_Xi_1(i) == Xi_1(i)
        continue
    else
        
        % Correct for the algorithm indexing by zero
        ki = i-1;
        
        % Calculate the new control points
        % Loop through indexes from k-p+1 to k and calculate the new control
        % points
        for j = ki-p_1+1:ki
            a(j) = (new_Xi_1(i)-Xi_1(j))/(Xi_1(j+p_1)-Xi_1(j));
            Q(j,:,:) = (1-a(j))*P(j-1,:,:) + a(j)*P(j,:,:);
        end
        
        % Move control points below ki-1 down by 1 space to
        % make room for the p new points;
        P(ki+1:size(P,1)+1,:,:) = P(ki:size(P,1),:,:);
        % Insert the new control points
        P(ki-p_1+1:ki,:,:) = Q(ki-p_1+1:ki,:,:);
        
        % move the knots to the right of the knot to be inserted over one
        % index
        Xi_1(i+1:end+1) = Xi_1(i:end);
        
        % Insert the knot to the knot vector.
        Xi_1(i) = new_Xi_1(i);
    end
end

for i = 1:length(new_Xi_2)
    % Check to see if the current knot in KV matches the current knot in
    % KVF. If it does, great, move on to the next knot, if not, calculate
    % new control points and insert KVF(i) into KV(i).
    if new_Xi_2(i) == Xi_2(i)
        continue
    else
        
        % Correct for the algorithm indexing by zero
        ki = i-1;
        
        % Calculate the new control points
        % Loop through indexes from k-p+1 to k and calculate the new control
        % points
        for j = ki-p_2+1:ki
            a_2(j) = (new_Xi_2(i)-Xi_2(j))/(Xi_2(j+p_2)-Xi_2(j));
            Q_2(:,j,:) = (1-a_2(j))*P(:,j-1,:) + a_2(j)*P(:,j,:);
        end
        
        % Move control points below ki-1 down by 1 space to
        % make room for the p new points;
        P(:,ki+1:size(P,2)+1,:) = P(:,ki:size(P,2),:);
        % Insert the new control points
        P(:,ki-p_2+1:ki,:) = Q_2(:,ki-p_2+1:ki,:);
        
        % move the knots to the right of the knot to be inserted over one
        % index
        Xi_2(i+1:end+1) = Xi_2(i:end);
        
        % Insert the knot to the knot vector.
        Xi_2(i) = new_Xi_2(i);
    end
end
new_w = P(:,:,d+1);
for i = 1:d
    new_P(:,:,i) = P(:,:,i)./new_w;
end

new_n_1 = size(new_P,1);
new_n_2 = size(new_P,2);
return 

