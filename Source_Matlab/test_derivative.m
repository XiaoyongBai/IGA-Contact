p_1=2;
p_2=2;
n_1=5;
n_2=5;
Xi_1=[0,0,0,1,2,3,3,3];
Xi_2=[0,0,0,1,2,3,3,3];                    
% P=[0,0; 0.25,0; 0.5,0;  0.75,0; 1,0;  
%    0,0.125; 0.25,0.125; 0.5,0.125;  0.75,0.125; 1,0.125;
%    0,0.25; 0.25,0.25; 0.5,0.25;  0.75,0.25; 1,0.25;
%    0,0.375; 0.25,0.375; 0.5,0.375;  0.75,0.375; 1,0.375;
%    0,0.5; 0.25,0.5; 0.5,.50;  0.75,0.5; 1,0.5
%    ];

P=[0,0; 0.1,0;0.2,0;0.3,0;0.4,0];
w=zeros(length(P), 1)+1;

Slave_P=[0,0.1];
[Closet_Point, Master_P] = Closet_Point_Curve(Slave_P, p_1, Xi_1, P, w )

% 
% xi=1.5;
% d_order=1;
% 
% tan1=Spline_Curve_derivatives(xi, d_order, p_1, Xi_1, P)
% 
% tan2=NURBS_Curve_derivatives(xi, d_order, p_1, Xi_1, P, w)
% 
d_circle =2;
p_circle = 2;
n_circle =3;
knots_circle = [0,0,0,1,1,1];
P=[-1,0; -1,1; 0,1];
w=[1;1/sqrt(2);1];
xi=0.6;
d_order=2;

% tan2=NURBS_Curve_derivatives(xi, d_order, p_circle, knots_circle, P, w)
Slave_P=[0,-10];
[Closet_Point, Master_P] = Closet_Point_Curve(Slave_P, p_circle, knots_circle, P, w )


% xi=0:0.01:1;
% x_grid=zeros(length(xi),2);
% 
% for j=1:length(xi)
%     x_grid(j,:)=NURBS_Curve_Point(xi(j), d_circle, p_circle, n_circle, knots_circle, P, w );
% end
% figure;
% plot(x_grid(:,1), x_grid(:,2));
a=2;