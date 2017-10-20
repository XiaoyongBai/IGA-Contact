% Problem = 0: test problem
%           1: Hertzh problem

global Problem;
Problem = 1;

%number of patch in the problem
Patch_num=2;

global state; %1=Plane strain, 2= Plane stress
state=1;
global Ey; %Young's Modulus
Ey=2e6;
global nu; %Possion's Ration
nu=0;

global time; %load control
time=1;

SD=2;


%initialize the dispalcement
Patch=1;
[~, ~, ~, ~, ~, ~, P, ~, ~] = Geometry(Problem, Patch, 0, 0);
dis_1=zeros(size(P));

Patch=2;
[~, ~, ~, ~, ~, ~, P, ~, ~] = Geometry(Problem, Patch, 0, 0);
dis_2=zeros(size(P));


while(time <100)
    %Form finite element matrix for Patch 1
    Patch=1;
    [p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q] = Geometry(Problem, Patch, dis_1, dis_2);

    [K_Patch_1, F_Patch_1] = Linear_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q, Problem, Patch);
    ncp_1=n_1*n_2; %number of control points of Patch 1

    %Form Finite element matrix for Patch 2
    Patch=2;
    [p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q] = Geometry(Problem, Patch, dis_1, dis_2);
    
    [K_Patch_2, F_Patch_2] = Linear_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q, Problem, Patch);
    ncp_2=n_1*n_2; %number of control points of Patch 2

    %Assemble Finite elment matrix
    K_Global=[K_Patch_1,zeros(size(K_Patch_1, 1),size(K_Patch_2,2)) ;zeros(size(K_Patch_2, 1),size(K_Patch_1,2)) ,K_Patch_2];
    F_Global=[F_Patch_1; F_Patch_2];

                   
    %incoporate contact constraint
    Penalty=1e8;%penalty parameter
    K_Contact = Contact_Formualtion(ncp_1, ncp_2, Problem, Penalty, dis_1, dis_2);

    K_Global=K_Global+K_Contact;

    d=linsolve(K_Global, F_Global);
    d1=d(1:SD*ncp_1);
    d2=d(SD*ncp_1+1:SD*(ncp_1+ncp_2));

    dis_1=dis_1+transpose(reshape(d1, size(dis_1,2), size(dis_1,1)));
    dis_2=dis_2+transpose(reshape(d2, size(dis_2,2), size(dis_2,1)));

    time=time+1;    
end
    
    
%plot the results
figure
axis equal
hold on;

amp=1;
field=2;

%get the geoemetry of Patch 1  
Patch=1;
[p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q] = Geometry(Problem, Patch, dis_1, dis_2);

% d=zeros(2*length(P),1)+1;
[von_Mises, max_location]=Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d1, field, amp);

%get the geoemetry of Patch 1  
Patch=2;
[p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, n_q] = Geometry(Problem, Patch, dis_1, dis_2);

% d=zeros(2*length(P),1)+1;
[von_Mises, max_location]=Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d2, field, amp);

hold off;

% field=3;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);
% 
% field=4;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);
 
% field=5;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);

% field=6;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);

% field=7;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);
% 
% field=8;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);
% 
% field=9;
% Plot_Elasticity(p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, field, amp);
% 

% [L2, energy_diff, nel] = L2_norm(Problem, p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w, d, n_q);