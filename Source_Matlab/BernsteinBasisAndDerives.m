%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function: evaluate bernstein basis functions and its derivatives
%Input: p_1=polynomial degree in direction 1
%       p_2=polynomial degree in direction 2
%       z1=shape coordinate in direction 1
%       z2=shape coordinate in direction 2
%Output: B=value of Bernstein basis function at point (z1, z2)
%        d_B_z1=derivative of B in direction 1
%        d_B_z2=derivative of B in direction 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B, d_B_z1, d_B_z2 ] = BernsteinBasisAndDerives(p_1, p_2, xi_1, xi_2)
    
    if p_1==2 & p_2==2
        B_z1=[(1-xi_1)^2, 2*xi_1*(1-xi_1), xi_1^2]; %Bernstein basis in direction 1;
        B_z2=[(1-xi_2)^2, 2*xi_2*(1-xi_2), xi_2^2]; %Bernstein basis in direction 2;
        
        d_B1_z1=[-2*(1-xi_1),  2-4*xi_1,  2*xi_1];
        d_B2_z2=[-2*(1-xi_2),  2-4*xi_2,  2*xi_2];
                
        %Bernstein basis for 9 control points
        B =[B_z1(1)*B_z2(1), B_z1(2)*B_z2(1),B_z1(3)*B_z2(1),B_z1(1)*B_z2(2), B_z1(2)*B_z2(2),B_z1(3)*B_z2(2),B_z1(1)*B_z2(3), B_z1(2)*B_z2(3),B_z1(3)*B_z2(3)];
        d_B_z1=[d_B1_z1(1)*B_z2(1),d_B1_z1(2)*B_z2(1),d_B1_z1(3)*B_z2(1),d_B1_z1(1)*B_z2(2),d_B1_z1(2)*B_z2(2),d_B1_z1(3)*B_z2(2),d_B1_z1(1)*B_z2(3),d_B1_z1(2)*B_z2(3),d_B1_z1(3)*B_z2(3)];
        d_B_z2=[B_z1(1)*d_B2_z2(1),B_z1(2)*d_B2_z2(1),B_z1(3)*d_B2_z2(1),B_z1(1)*d_B2_z2(2),B_z1(2)*d_B2_z2(2),B_z1(3)*d_B2_z2(2),B_z1(1)*d_B2_z2(3),B_z1(2)*d_B2_z2(3),B_z1(3)*d_B2_z2(3)];
    elseif p_1==3 & p_2==3
        B_z1=[(1-xi_1)^3, 3*xi_1*(1-xi_1)^2, 3*xi_1^2*(1-xi_1), xi_1^3]; %Bernstein basis in direction 1;
        B_z2=[(1-xi_2)^3, 3*xi_2*(1-xi_2)^2, 3*xi_2^2*(1-xi_2), xi_2^3]; %Bernstein basis in direction 2;
        
        d_B1_z1=[-3*(1-xi_1)^2, 3-12*xi_1+9*xi_1^2, 6*xi_1-9*xi_1^2, 3*xi_1^2];
        d_B2_z2=[-3*(1-xi_2)^2, 3-12*xi_2+9*xi_2^2, 6*xi_2-9*xi_2^2, 3*xi_2^2];
                
        %Bernstein basis for 9 control points
        for j=1:p_2+1
            for i=1:p_1+1
                B(i+(j-1)*(p_1+1)) = B_z1(i)*B_z2(j);
                d_B_z1(i+(j-1)*(p_1+1)) = d_B1_z1(i)*B_z2(j);
                d_B_z2(i+(j-1)*(p_1+1)) = B_z1(i)*d_B2_z2(j);
            end
        end
                
    else
        error('BernsteinBasisAndDerives::unsupported polynomial degree');
    end

end

