%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: IEN2D
%
% Input:  n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%
% Output: IEN2D = 2D IEN array
%
% Purpose: Compute the 2D IEN array
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IEN] = IEN2D(n_1,n_2,p_1,p_2,Xi_1,Xi_2)

    IEN1 = IEN1D(n_1,p_1,Xi_1);
    IEN2 = IEN1D(n_2,p_2,Xi_2);

    n_el_1 = size(IEN1,2);
    n_el_2 = size(IEN2,2);

    for e2 = 1:1:n_el_2
        for a2 = 1:1:(p_2+1)
            i2 = IEN2(a2,e2);

            for e1 = 1:1:n_el_1
                for a1 = 1:1:(p_1+1)
                    i1 = IEN1(a1,e1);

                    e = (e2-1)*n_el_1+e1;
                    a = (a2-1)*(p_1+1)+a1;
                    i = (i2-1)*n_1+i1;

                    IEN(a,e) = i;
                end
            end

        end
    end
%     
%     for e1 = 1:1:n_el_1
%         for a1 = 1:1:(p_1+1)
%             i1 = IEN1(a1,e1);
% 
%             for e2 = 1:1:n_el_2
%                 for a2 = 1:1:(p_2+1)
%                     i2 = IEN2(a2,e2);
% 
%                     e = (e1-1)*n_el_2+e2;
%                     a = (a1-1)*(p_2+1)+a2;
%                     i = (i1-1)*n_1+i2;
% 
%                     IEN(a,e) = i;
%                 end
%             end
% 
%         end
%     end

end