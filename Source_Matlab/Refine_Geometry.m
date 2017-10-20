function [new_n_1, new_n_2, new_Xi_1, new_Xi_2, new_P, new_w ] = Refine_Geometry(refine_number, p_1, p_2, n_1, n_2, Xi_1, Xi_2, P, w)
   
    if refine_number >0   
        %add knots between two adjacent knots that are different
        add_Xi_1=zeros(1,0);
        for ni=1:length(Xi_1)-1
            if(Xi_1(ni) ~= Xi_1(ni+1))
                interval=(Xi_1(ni+1)-Xi_1(ni))/(refine_number+1);
                for ri=1:refine_number
                    add_Xi_1=[add_Xi_1 Xi_1(ni)+interval*ri]; 
                end
            end
        end
        
        add_Xi_2=zeros(1,0);
        for ni=1:length(Xi_2)-1
            if(Xi_2(ni) ~= Xi_2(ni+1))
                interval=(Xi_2(ni+1)-Xi_2(ni))/(refine_number+1);
                for ri=1:refine_number
                    add_Xi_2=[add_Xi_2 Xi_2(ni)+interval*ri]; 
                end
            end
        end
        
        
        %convert format of P and w to be used by NURBS_Surface_Refine
        P_temp=zeros(n_1, n_2, 2);
        w_temp=zeros(n_1, n_2);

        for ni_1=1:n_1
            for ni_2=1:n_2
                for di=1:2
                    P_temp(ni_1, ni_2, di) = P((ni_2-1)*n_1+ni_1, di);
                end
                w_temp(ni_1, ni_2) = w((ni_2-1)*n_1+ni_1);
            end
        end

        [new_n_1,new_n_2,new_Xi_1,new_Xi_2,P_temp,w_temp] = NURBS_Surface_Refine(2,add_Xi_1,add_Xi_2,p_1,p_2,n_1,n_2,Xi_1,Xi_2,P_temp,w_temp);

        new_P=zeros(new_n_1*new_n_2,2);
        new_w=zeros(new_n_1*new_n_2,1);
        %get P and w by conversion
        for ni_1=1:new_n_1
            for ni_2=1:new_n_2
                for di=1:2
                     new_P((ni_2-1)*new_n_1+ni_1, di) = P_temp(ni_1, ni_2, di);
                end
                new_w((ni_2-1)*new_n_1+ni_1)= w_temp(ni_1, ni_2) ;
            end
        end
    else
        new_n_1=n_1;
        new_n_2=n_2;
        new_Xi_1=Xi_1;
        new_Xi_2=Xi_2;
        new_P=P;
        new_w=w; 
        
    end


end

