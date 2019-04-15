function [x_next,output_k] = state_update(build,bm,sim,x,v,u,k)
    x_next = [];
    output_k = [];
    x_k = x;
    
    for i = k:1:k+sim.Tolc-1    % apply in open-loop fashion
        Bv_help = [];
        Bx_help = [];
        for j=1:size(bm.Bvu,3)  % sum over all inputs
            Bv_help = [Bv_help,bm.Bvu(:,:,j)*v(i,:)'];
            Bx_help = [Bx_help,bm.Bxu(:,:,j)*x_k];
        end                    

        x_next = [x_next, bm.A*x_k + (bm.Bu + Bv_help + Bx_help)*u((i-k+1),:)' + bm.Bv*v(i,:)'];
        Dv_help = [];
        for j=1:size(bm.Bvu,3)
            Dv_help = [Dv_help,bm.Dvu(:,:,j)*v(i,:)'];
        end            
        output_k = [output_k, bm.C*x_k + (bm.Du+Dv_help)*u((i-k+1),:)' + bm.Dv*v(i,:)'];
        x_k = x_next(:,end);
    end
end 