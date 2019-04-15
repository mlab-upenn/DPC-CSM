function [u,fval,iter,time] = mydpc(models, leafmodels, k, xref, x, d, proxy, obj, N, umin, umax, xmin, xmax, solver, Q, R, epsc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solves dpc-problem in condensed form
%
% Model: x(k+1) = A*x(k) + B*u(k) + Bd*d(k)
% DPC cost = sum_0^N (x(k)-xref)'*Q*(x(k)-xref) + sum_0^(N-1) u(k)'*R*u(k)
%
% Inputs:
%  model     - regression tree model
%  Q,R,P,N   - MPC tuning parameters
%  umin,umax - box constraints for inputs
%  xmin,xmax - box constraints for states (specify -inf,inf to run without state constraints)
%  k         - current time step
%  x         - current state of the plant
%  d         - current disturbance affecting the plant
%  solver    - chosen solver to solve the dpc problem
%  z0        - starting point for some solvers
%
% Outputs:
%  u      - first input to be applied to the plant
%  fval   - value of objective function at optimal control sequence
%  iter   - number of needed iterations
%  time   - tic-toc-time for calling solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate leaf coeffs

% prepare features_d to find an appropriate leaf
dcur = d(:,k:k+N-1);
features_d = [dcur(:)', x, proxy];

coeff = zeros(4*N+1,N);
for idm = 1:N
    coeff(1:4*idm+1,idm) = find_linearmodel_in_leaves(models{idm}, leafmodels{idm}, features_d);
end

% constraints
lb = umin(k:k+N-1,:)';
ub = umax(k:k+N-1,:)';
lb = [lb(:);zeros(N,1)]; % add slack>0
ub = [ub(:);inf(N,1)];  % add slack<inf
Aineq1 = coeff(2:end,:)';
bineq1 = xmax(k+1:k+N,1)-coeff(1,:)';
Aineq2 = -coeff(2:end,:)';
bineq2 = coeff(1,:)'-xmin(k+1:k+N,1);
Aineq = [Aineq1, -eye(N);Aineq2, -eye(N)];  % add slack using eye(N)
bineq = [bineq1;bineq2];

% input rate change constraint
% UB = uprev+udotmax;
% LB = uprev-udotmax;
% umin = max(umin,LB);
% umax = min(UB,umax);

%% solve optimization problem

switch solver
    
    case 'cplex'
        
        QQ = diag(Q(1,1)*ones(N,1));
        RR = kron(eye(N),R);
        xref = xref*ones(N,1);
        C1 = coeff(1,:);
        C2 = coeff(2:end,:);
        H = 2*( RR +  C2*QQ*C2' );
        H = (H+H')/2;
        elecCost = obj(k:k+N-1,:)';
        f = 2*C2*QQ*C1'-2*C2*QQ*xref + elecCost(:);
        c = C1*QQ*C1'-C1*QQ*xref-xref'*QQ*C1'+xref'*QQ*xref;
        H = blkdiag(H, 0*eye(N)); f = [f; epsc*ones(N,1)];  % add slack
        tic;
        [U, fval] = cplexqp(H,f',Aineq,bineq,[],[],lb,ub);
        fval = fval+c;
        time = toc;
        iter = 0;
        
    case 'linprog'
        
        f = obj(k:k+N-1,:)';
        f = [f(:); 10^3*ones(N,1)]; % add slack
        options=optimset;
        options.Display='off';
        tic;
        [U, fval] = linprog(f,Aineq,bineq,[],[],lb(:),ub(:),[],options);
        time=toc;
        iter=nan;
        
end
    
u = U(1:4);

% instataneous correction for insufficient light
bldgmodel = load('data/bm.mat');
bm = bldgmodel.bm;
light = bm.Du(2,2)*u(2) + d(4,1)*u(1);
if light + 1e-10 < xmin(k+1,2)    % not enough light during this time
    disp(['old light input: ' num2str(u(2))]);
    u(2) = (xmin(k+1,2) - d(4,1)*u(1))/bm.Du(2,2);
    disp(['new light input: ' num2str(u(2))]);
end
