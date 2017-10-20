%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: return the constitutive matrix
%          Voigt notation is used
%
%Input(global variable): 
%       EY=Young's modulus
%       nu=Possion's ratio
%       state= 1:Plane strain
%              2:Plane stress
%              3:3D
%
%Output: D_Matrix=Constitutive matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function [ D_Matrix ] = D_local( )

    global Ey;
    global nu;
    global state;
    
    %Lame constants
    lambda=Ey*nu/((1+nu)*(1-2*nu));
    mu=Ey/(2*(1+nu));
    
    switch state
        case 1  %plane strain         
            D_Matrix=[lambda+2*mu, lambda, 0;
                      lambda, lambda+2*mu, 0;
                      0,0,mu];
        otherwise
            error('D_local:unsupported type of problem');                      
    end


end

