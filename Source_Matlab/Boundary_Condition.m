%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: set boundary condition for simple cases
%
%Input: Problem = 0: test problem
%                 1: Hertz problem

%       p_1, p_2, n_1, n_2, Xi_1, Xi_2=NURBS control points
%
%Output: BC=Dirichlet boundary mark
%        g=Dirichlet boundary value (fixed displacement)
%        Newmann=Newmann boundary mark
%        Newmann_h_loc=Newmann boundary value (traction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ BC, g, Newmann] = Boundary_Condition(Problem, Patch, p_1, p_2, n_1, n_2, Xi_1, Xi_2)
       
    switch Problem
        case 0           
            [ne_1, Ce] = Extraction1D(n_1,p_1,Xi_1);
            [ne_2, Ce] = Extraction1D(n_2,p_2,Xi_2);
            
            ncp=n_1*n_2;
            
            BC=zeros(2,ncp);
            g=zeros(2,ncp);
            Newmann=zeros(2,4, ne_1*ne_2);
            
            switch Patch
                case 1
                    BC(1:2, 1:n_1) =1;
                    g(1:2, 1:n_1) =0;
                case 2
                    BC(1, ncp-n_1+1:ncp)=1;                   
                    g(1, ncp-n_1+1:ncp)=0;
                    
%                     BC(2, ncp-n_1+1:ncp)=1;                   
%                     g(2, ncp-n_1+1:ncp)=-0.01;
                                         
                     Newmann(2, 3, 7:9)=1;
            end
                                 
        case 1
            [ne_1, Ce] = Extraction1D(n_1,p_1,Xi_1);
            [ne_2, Ce] = Extraction1D(n_2,p_2,Xi_2);
            
            ncp=n_1*n_2;
            
            BC=zeros(2,ncp);
            g=zeros(2,ncp);
            Newmann=zeros(2,4, ne_1*ne_2);
            
            switch Patch
                case 1
                    BC(1:2, 1:n_1) =1;
                    g(1:2, 1:n_1) =0;
                case 2
                    BC(1, ncp-n_1+1:ncp)=1;                   
                    g(1, ncp-n_1+1:ncp)=0;
                    
                    BC(2, ncp-n_1+1:ncp)=1;                   
                    g(2, ncp-n_1+1:ncp)=-0.001;
                                         
%                      Newmann(2, 3, ne_1*(ne_2-1)+1:ne_1*ne_2)=1;
            end   
        otherwise 
            error('unsupported problem id');
    end


end

