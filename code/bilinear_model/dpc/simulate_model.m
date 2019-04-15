function [x_next,output_k] = simulate_model(model, x, v, u)

% A = model.A;
% Bu = model.Bu;
% Bv = model.Bv;
% Bvu = model.Bvu;
% Bxu = model.Bxu;
% C = model.C;
% Du = model.Du;
% Dv = model.Dv;
% Dvu = model.Dvu;

x_next = [];
output_k = [];
x_k = x;

Bv_help = [];
Bx_help = [];
for j=1:size(model.Bvu,3)  % sum over all inputs
    Bv_help = [Bv_help,model.Bvu(:,:,j)*v]; %#ok<*AGROW>
    Bx_help = [Bx_help,model.Bxu(:,:,j)*x_k];
end

x_next = [x_next, model.A*x_k + (model.Bu + Bv_help + Bx_help)*u + model.Bv*v];
Dv_help = [];
for j=1:size(model.Bvu,3)
    Dv_help = [Dv_help,model.Dvu(:,:,j)*v];
end
output_k = [output_k, model.C*x_k + (model.Du+Dv_help)*u + model.Dv*v];

end