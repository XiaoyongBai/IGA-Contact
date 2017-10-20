%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: generate geometry information for analysis,
%          Serving as pre-processing
%
%input: Problem = 0: test problem of two blocks
%                 1: pipe contact plate problem
%       Patch=Patch number for multi-patches problems
%       dis_1=displacement of patch 1
%       dis_2=displacement of patch 2
%
%Output: p_1=polynomial degree at direction 1
%        p_2=polynomial degree at direction 2
%        n_1=number of basis functions at direction 1
%        n_2=number of basis functions at direction 2
%        Xi_1=knots vector at direction 1
%        Xi_2=knots vector at direction 2
%        P=coordinated of control points
%        w=weights of basis functions
%        n_q=number of integration points per element
function [p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q] = Geometry(Problem, Patch, dis_1, dis_2)

    global time;
    
    gap=-0.0001;

    switch Problem %options of Problem 
        case 0 %two block problem
            switch Patch %options of Patch
                case 1 %block at the bottom
                    p_1=2;
                    p_2=2;
                    n_1=5;
                    n_2=5;
                    Xi_1=[0,0,0,1,2,3,3,3];
                    Xi_2=[0,0,0,1,2,3,3,3];                    
                    P=[0,0; 0.25,0; 0.5,0;  0.75,0; 1,0;  
                       0,0.125; 0.25,0.125; 0.5,0.125;  0.75,0.125; 1,0.125;
                       0,0.25; 0.25,0.25; 0.5,0.25;  0.75,0.25; 1,0.25;
                       0,0.375; 0.25,0.375; 0.5,0.375;  0.75,0.375; 1,0.375;
                       0,0.5; 0.25,0.5; 0.5,.50;  0.75,0.5; 1,0.5
                       ];
                   %update the coordinates using displacement
                    P=P+ dis_1;
                   w=zeros(n_1*n_2, 1)+1;
                   n_q= p_1+1;
                case 2
                    p_1=2;
                    p_2=2;
                    n_1=5;
                    n_2=5;
                    Xi_1=[0,0,0,1,2,3,3,3];
                    Xi_2=[0,0,0,1,2,3,3,3];                    
                    P=[0.2,0.5; 0.3,0.5; 0.4,0.5;  0.5,0.5; 0.6,0.5;  
                       0.2,0.55; 0.3,0.55; 0.4,0.55;  0.5,0.55; 0.6,0.55;
                       0.2,0.6; 0.3,0.6; 0.4,0.6;  0.5,0.6; 0.6,0.6;
                       0.2,0.65; 0.3,0.65; 0.4,0.65;  0.5,0.65; 0.6,0.65;
                       0.2,0.7; 0.3,0.7; 0.4,0.7;  0.5,0.7; 0.6,0.7
                       ];
                    P(:,2)=P(:,2)+gap; %create contact gap between two patches
                    P(:,1)=P(:,1);
                    
                    %update the coordinates using displacement         
                    P=P+dis_2;
                   w=zeros(n_1*n_2, 1)+1;
                   n_q= p_1+1;
                                      
                otherwise
                    error('Geometry:unsupported patch number');                 
            end
           
        case 1 %pipe contact block
            switch Patch %options of Patch
                case 1 %block at the bottom
                    p_1=2;
                    p_2=2;
                    n_1=5;
                    n_2=5;
                    Xi_1=[0,0,0,1,2,3,3,3];
                    Xi_2=[0,0,0,1,2,3,3,3];                    
                    P=[0,0; 0.25,0; 0.5,0;  0.75,0; 1,0;  
                       0,0.125; 0.25,0.125; 0.5,0.125;  0.75,0.125; 1,0.125;
                       0,0.25; 0.25,0.25; 0.5,0.25;  0.75,0.25; 1,0.25;
                       0,0.375; 0.25,0.375; 0.5,0.375;  0.75,0.375; 1,0.375;
                       0,0.5; 0.25,0.5; 0.5,.50;  0.75,0.5; 1,0.5
                       ];
                   %update the coordinates using displacement
                   P=P+ dis_1;
                   w=zeros(n_1*n_2, 1)+1;
                   n_q= p_1+1;
                case 2
                    p_1=2;
                    p_2=2;
                    n_1=3;
                    n_2=3;
                    Xi_1=[0,0,0,1,1,1];
                    Xi_2=[0,0,0,1,1,1];
                    Radius=0.3;
                    Thickness=0.1;
                    
                    P(1:3, 1:2)=[0.5-Radius/sqrt(2), 0.5+Radius/sqrt(2); 0.5, 0.5; 0.5+Radius/sqrt(2), 0.5+Radius/sqrt(2) ];
                    P(1:3, 2)=P(1:3, 2)-Radius*(sqrt(2)-1);
                    P(4:6, 1:2)=P(1:3,1:2);
                    P(4:6, 2)=P(4:6,2)+Thickness/2;
                    P(7:9, 1:2)=P(1:3,1:2);
                    P(7:9, 2)=P(7:9,2)+Thickness;
                   
                    
                    %update the coordinates using displacement         
                    w=[1;sqrt(2)/2;1;1;sqrt(2)/2;1;1;sqrt(2)/2;1];
                    n_q= p_1+1;
                   
                   refine_times=2;
                   [n_1, n_2, Xi_1, Xi_2, P, w] = Refine_Geometry(refine_times, p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w);
                   
                   P(:,2)=P(:,2)+gap; %create contact gap between two patches
                   
                   P=P+dis_2;
                   
                otherwise
                    error('Geometry:unsupported patch number');                 
            end
       otherwise
            error('Geometry::unsupported problem id');

    end




end

