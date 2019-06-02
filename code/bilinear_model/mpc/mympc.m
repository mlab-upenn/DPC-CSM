function [input, fval] = mympc(x0, xref, idk, bm, N, vPred, ymin_all, ymax_all, cost_all, umin_all, umax_all, solver, Q, R, epsc)

% collect constraints
umin = umin_all(idk:idk+N-1,:);
umax = umax_all(idk:idk+N-1,:);
ymin = ymin_all(idk+1:idk+N,:);
ymax = ymax_all(idk+1:idk+N,:);
ucost = cost_all(idk:idk+N-1,:);
    
yalmip('clear')

% Number of states and inputs
nx = size(bm.A,1);
ny = size(bm.C,1);
nu = size(bm.Bu,2);

% Variables definition
u = sdpvar(repmat(nu,1,N),ones(1,N));

% Slack variables for soft constraints
eps = sdpvar(ny*ones(1,N),ones(1,N));

% Yalmip initialization problem

constraints = [];
objective   = 0;

for jj = 1 : N  % for inputs u_0,...,u_(N_mpc-1)
    Bvu_aux{jj} = nan(nx,nu);
    Bxu_aux{jj} = nan(nx,nu);
    Dvu_aux{jj} = nan(ny,nu);
    
    for ii = 1 : nu
        Bvu_aux{jj}(:,ii) = bm.Bvu(:,:,ii) * vPred(jj,:)'; % linearized
        Bxu_aux{jj}(:,ii) = bm.Bxu(:,:,ii) * x0;   % linearized
        Dvu_aux{jj}(:,ii) = bm.Dvu(:,:,ii) * vPred(jj,:)';  % linearized
    end
end

xref = xref*ones(12,1);
x = x0;

for k = 1:N
    
    Bu_tmp = bm.Bu + Bvu_aux{k} + Bxu_aux{k};
    x = bm.A*x+ Bu_tmp*u{k}+bm.Bv*vPred(k,:)';
    
    objective = objective + (x-xref)'*Q*(x-xref) + u{k}'*R*u{k} + ucost(k,:)*u{k} + epsc*[1, 0, 0]*eps{k};
    constraints = [constraints, [ymin(k,:)' - eps{k} <= bm.C*x + (Dvu_aux{k}+bm.Du)*u{k} + bm.Dv*vPred(k,:)' <= ymax(k,:)' + eps{k};...
        umin(k,:)'    <=                     u{k}                                <= umax(k,:)';...
        eps{k}                             >= 0]];
end


switch solver
    
    case 'ipopt'
        ops = sdpsettings('solver','ipopt','verbose',1);
        diagnostics=optimize(constraints,objective,ops);
        if diagnostics.problem ~= 0
            fprintf('Yalmip error: %s\n', yalmiperror(diagnostics.problem));
        end
        
    case 'cplex'
        ops = sdpsettings('solver','cplex','verbose',0);
        diagnostics=optimize(constraints,objective,ops);
        if diagnostics.problem ~= 0
            fprintf('Yalmip error: %s\n', yalmiperror(diagnostics.problem));
        end
        
    case 'gurobi'
        ops = sdpsettings('solver','gurobi','verbose',0);
        diagnostics=optimize(constraints,objective,ops);
        if diagnostics.problem ~= 0
            fprintf('Yalmip error: %s\n', yalmiperror(diagnostics.problem));
        end
        
    case 'linprog'
        ops = sdpsettings('solver','linprog','verbose',1);
        diagnostics=optimize(constraints,objective,ops);
        if diagnostics.problem ~= 0
            fprintf('Yalmip error: %s\n', yalmiperror(diagnostics.problem));
        end
        
end

input = double(u{1});
fval = double(objective);

% instataneous correction for insufficient light
light = bm.Du(2,2)*input(2) + vPred(1,4)*input(1);
if light + 1e-10 < ymin(1,2)    % not enough light during this time
    disp(['old light input: ' num2str(input(2))]);
    input(2) = (ymin(1,2) - vPred(1,4)*input(1))/bm.Du(2,2);
    disp(['new light input: ' num2str(input(2))]);
end

end