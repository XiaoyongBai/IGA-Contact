function [ h_loc ] = H_loc(x)

    h_loc=zeros(2,1);
    global Problem;
    
    switch Problem 
        case 0
            h_loc(1)=0;
            h_loc(2)=-1e4;
        case 1
            h_loc(1)=0;
            h_loc(2)=-1e5;
        otherwise
            error('H_loc:unsupported problem type');
    end
    


end

