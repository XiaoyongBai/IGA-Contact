%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Plot_Temperature
%
% Input:  p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w=NURBS control parameters
%         d=displacement result of the control points
%         field=the component that needed to be plot, choices are listed
%              within the code
%         amp=displacement amplification factor
%
%Output:
%       stress_von_max=maximum von mises stress
%       stress_von_max_position=position of maximum von mises stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stress_von_max, stress_von_max_position]=Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,field, amp)
    
    %Possion' ratio
    global nu;
    %elatic matrix
    D_matrix=D_local();
    
    %spatial Dimension
    sd = 2;

    % Extract the basis, geometry, and temperature field
    [n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2);
    [P_b,w_b] = Extract_Geometry(sd,p_1,p_2,n_el,C_operators,IEN,P,w);
    
    %set NURBS weight for element control points
    ncp_ele=(p_1+1)*(p_2+1); %number of control points in one element
    We=zeros(ncp_ele, n_el);
    for ee=1:n_el
        for cpi=1:ncp_ele
            cp=IEN(cpi,ee);
            We(cpi,ee)=w(cp);
        end
    end            
    
    %sample point in parent element
    x_plot=0:0.1:1;
    y_plot=0:0.1:1;
    
    
    x_mat=[];
    y_mat=[];
    
    %field variables that going to be computed and plotted
    dx_mat=[];
    dy_mat=[];
    dr_mat=[];   
    stress_xx_mat=[];
    stress_yy_mat=[];
    stress_xy_mat=[];
    stress_zz_mat=[];
    stress_p_max_mat=[];
    stress_p_min_mat=[];
    stress_von_mat=[];
      
    n_el_1 = length(unique(Xi_1))-1;
    n_el_2 = length(unique(Xi_2))-1;
    
    for e_2=1:n_el_2 %loop over element
        x_row=[];
        y_row=[];
        
        dx_row=[];
        dy_row=[];
        dr_row=[];
        stress_xx_row=[];
        stress_yy_row=[];
        stress_xy_row=[];
        stress_zz_row=[];
        stress_p_max_row=[];
        stress_p_min_row=[];
        stress_von_row=[];
        
        for e_1=1:n_el_1
            e_id=(e_2-1)*n_el_1+e_1;
            
            C_e=C_operators(:,:,e_id);
            Pb_e=P_b(:,:,e_id);
            wb_e=w_b(:,e_id);
            We_e=We(:,e_id);

            x_ele=zeros(length(x_plot), length(y_plot));
            y_ele=zeros(length(x_plot), length(y_plot));
            
            dx_ele=zeros(length(x_plot), length(y_plot));
            dy_ele=zeros(length(x_plot), length(y_plot));
            dr_ele=zeros(length(x_plot), length(y_plot));
            stress_xx_ele=zeros(length(x_plot), length(y_plot));
            stress_yy_ele=zeros(length(x_plot), length(y_plot));
            stress_xy_ele=zeros(length(x_plot), length(y_plot));
            stress_zz_ele=zeros(length(x_plot), length(y_plot));
            stress_p_max_ele=zeros(length(x_plot), length(y_plot));
            stress_p_min_ele=zeros(length(x_plot), length(y_plot));
            stress_von_ele=zeros(length(x_plot), length(y_plot));

            
            for xi=1:length(x_plot)
                for yi=1:length(y_plot)
                    xi_1=x_plot(xi);
                    xi_2=y_plot(yi);

                    [Ra, Ra_x, Ra_y, x, J_matrix, J_det] = Shape_Function( p_1, p_2, xi_1, xi_2, Pb_e, wb_e, We_e, C_e);
                    
                    stress_loc=zeros(sd+1,1);
                    %stiffness matrix and load vector
                    for a=1:ncp_ele
                        n_id=IEN(a,e_id);
                        d_node=[d(2*(n_id-1)+1); d(2*n_id)];
                        dx_ele(yi,xi)=dx_ele(yi,xi)+Ra(a)*d_node(1);
                        dy_ele(yi,xi)=dy_ele(yi,xi)+Ra(a)*d_node(2);
                        
                        B_a=[Ra_x(a), 0; 0, Ra_y(a); Ra_y(a), Ra_x(a)]; 
                        stress_loc=stress_loc+D_matrix*B_a*d_node;
                     end%for a
                     
                     x_ele(yi,xi)=x(1)+dx_ele(yi,xi);
                     y_ele(yi,xi)=x(2)+amp*dy_ele(yi,xi);
                                          
                     dr_ele(yi,xi)=sqrt(dx_ele(yi, xi)^2+dy_ele(yi,xi)^2);
                     stress_xx_ele(yi,xi)=stress_loc(1);
                     stress_yy_ele(yi,xi)=stress_loc(2);
                     stress_xy_ele(yi,xi)=stress_loc(3);
                     stress_zz_ele(yi,xi)=nu*(stress_xx_ele(yi,xi)+stress_yy_ele(yi,xi));
                     stress_p_max_ele(yi,xi)=max(abs(stress_xx_ele(yi,xi)), abs(stress_yy_ele(yi,xi)));
                     stress_p_min_ele(yi,xi)=min(abs(stress_xx_ele(yi,xi)), abs(stress_yy_ele(yi,xi)));
                     stress_von_ele(yi,xi)=sqrt( 0.5*((stress_xx_ele(yi,xi)-stress_yy_ele(yi,xi))^2+(stress_xx_ele(yi,xi)-stress_zz_ele(yi,xi))^2+(stress_yy_ele(yi,xi)-stress_zz_ele(yi,xi))^2 ) + 6*stress_xy_ele(yi,xi)^2);

                     
                end %for yi                    
             end%for xi
             
             x_row=[x_row x_ele];
             y_row=[y_row y_ele];
             dx_row=[dx_row dx_ele];
             dy_row=[dy_row dy_ele];
             dr_row=[dr_row dr_ele];
             stress_xx_row=[stress_xx_row stress_xx_ele];
             stress_yy_row=[stress_yy_row stress_yy_ele];
             stress_xy_row=[stress_xy_row stress_xy_ele];
             stress_zz_row=[stress_zz_row stress_zz_ele];
             stress_p_max_row=[stress_p_max_row stress_p_max_ele];
             stress_p_min_row=[stress_p_min_row stress_p_min_ele];
             stress_von_row=[stress_von_row stress_von_ele];
             
        end
        
        x_mat=[x_mat; x_row];
        y_mat=[y_mat; y_row];
        dx_mat=[dx_mat; dx_row];
        dy_mat=[dy_mat; dy_row];
        dr_mat=[dr_mat; dr_row];
        stress_xx_mat=[stress_xx_mat; stress_xx_row];
        stress_yy_mat=[stress_yy_mat; stress_yy_row];
        stress_xy_mat=[stress_xx_mat; stress_xy_row];
        stress_zz_mat=[stress_zz_mat; stress_zz_row];
        stress_p_max_mat=[stress_p_max_mat; stress_p_max_row];
        stress_p_min_mat=[stress_p_min_mat; stress_p_min_row];
        stress_von_mat=[stress_von_mat; stress_von_row];
        

    end
    
    index_1 = 1:length(x_plot);

    for e1 = 2:n_el_1
        first = (e1-1)*length(x_plot)+2;
        index_1 = [index_1 first:first+length(x_plot)-2];
    end

    index_2 = 1:length(y_plot);

    for e2 = 2:n_el_2
        first = (e2-1)*length(y_plot)+2;
        index_2 = [index_2 first:first+length(y_plot)-2];
    end

    x_mat = x_mat(index_2,index_1);
    y_mat = y_mat(index_2,index_1);
    dx_mat = dx_mat(index_2, index_1);
    dy_mat = dy_mat(index_2,index_1);
    dr_mat = dr_mat(index_2, index_1);
    stress_xx_mat=stress_xx_mat(index_2, index_1);
    stress_yy_mat=stress_yy_mat(index_2, index_1);
    stress_xy_mat=stress_xy_mat(index_2, index_1);
    stress_zz_mat=stress_zz_mat(index_2, index_1);
    stress_p_max_mat=stress_p_max_mat(index_2, index_1);
    stress_p_min_mat=stress_p_min_mat(index_2, index_1);
    stress_von_mat=stress_von_mat(index_2, index_1);


    % Contour plot of the chosen field with colorbar
    
    switch field
        case 1 %plot the x displacement
            contourf(x_mat,y_mat, dx_mat,20,'EdgeColor','none');
            h_lengend=legend('x Displacement');
        case 2 %plot the y displacement
            contourf(x_mat,y_mat, dy_mat,20,'EdgeColor','none');
             h_lengend=legend('y displacement');
        case 3 %plot the magnitude of displacement
            contourf(x_mat,y_mat, dr_mat,20,'EdgeColor','none');
            h_lengend=legend('Magnitude of Displacement');
        case 4 %plot the stress_xx
            contourf(x_mat,y_mat, stress_xx_mat,20,'EdgeColor','none');
            h_lengend=legend('Stress xx');
        case 5 %plot the stress_xy
            contourf(x_mat,y_mat, stress_xy_mat,20,'EdgeColor','none');
            h_lengend=legend('Stress xy');
        case 6 %plot the stress_yy
            contourf(x_mat,y_mat, stress_yy_mat,20,'EdgeColor','none');
            h_lengend=legend('Stress yy');
        case 7 %plot the stress_zz
            contourf(x_mat,y_mat, stress_zz_mat,20,'EdgeColor','none');
            h_lengend=legend('Stress zz');
        case 8 %plot the maximum in-plane principal stress
             contourf(x_mat,y_mat, stress_p_max_mat,20,'EdgeColor','none');
             h_lengend=legend('Maximum in plane Principal Stress');
        case 9 %plot the minimum in-plane principal stress
            contourf(x_mat,y_mat, stress_p_min_mat,20,'EdgeColor','none');
            h_lengend=legend('Minimum in plane Principal Stress');
        case 10 %plot the von Mises stress field
             contourf(x_mat,y_mat, stress_von_mat,20,'EdgeColor','none')
             h_lengend=legend('Von Mises Stress');
        otherwise
            error('Plot_Elasticity:unsupported field');
    end
            
    xlabel('x', 'fontsize', 15);
    ylabel('y');
    
     set(h_lengend,'fontsize',15);
    
%     colorbar;
    
    [stress_von_max, p_v]=max(stress_von_mat);
    [stress_von_max, p_s]=max(stress_von_max);
    max_position=[p_v(p_s), p_s];
    stress_von_max_position=[x_mat(max_position(1), max_position(2)), y_mat(max_position(1), max_position(2))];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%plot element boundaries
    %plot horizontal lines
    for line_y_i=1:n_el_2+1
        y_number_i=(line_y_i-1)*10+1;
        
        x_number=size(x_mat,2);
        
        for x_number_i=1:x_number-1
             x_1=x_mat(y_number_i, x_number_i);
             y_1=y_mat(y_number_i, x_number_i);
             x_2=x_mat(y_number_i, x_number_i+1);
             y_2=y_mat(y_number_i, x_number_i+1);
             line([x_1,x_2],[y_1,y_2],'color','k');
        end
    end
    
    %plot vertical lines
    for line_x_i=1:n_el_1+1
        x_number_i=(line_x_i-1)*10+1;
        
        y_number=size(x_mat,1);
        
        for y_number_i=1:y_number-1
             x_1=x_mat(y_number_i, x_number_i);
             y_1=y_mat(y_number_i, x_number_i);
             x_2=x_mat(y_number_i+1, x_number_i);
             y_2=y_mat(y_number_i+1, x_number_i);
             line([x_1,x_2],[y_1,y_2],'color','k');
        end
    end

    
 
    
    
end