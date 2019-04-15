function input = mympcbackup(x0, xref, bm, N, vPred, ymin, ymax, ucost, umin, umax, solver)

    yalmip('clear')

    % Number of states and inputs

    nx = size(bm.A,1);
    ny = size(bm.C,1);
    nu = size(bm.Bu,2);
    nv = size(bm.Bv,2);
    nslack = 2*ny*N;      % size of slackMax, slackMin

    % Variables definition

    x = sdpvar(repmat(nx,1,N+1),ones(1,N+1));
    u = sdpvar(repmat(nu,1,N),ones(1,N));

    % Slack variables for soft constraints

    eps_m = sdpvar(ny*ones(1,N+1),ones(1,N+1));
    eps_M = sdpvar(ny*ones(1,N+1),ones(1,N+1));

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
    Q = sqrt(diag(1e2*ones(12,1)));
    R = sqrt(diag(1e-3*ones(4,1)));
    for k = 1:N
        
%         objective = objective + ucost(k,:)*u{k} + (10^3)*(eps_m{k} + eps_M{k});
%         objective = objective + (norm(Q*(x{k+1}-xref),2)) + norm(R*u{k},2) + ucost(k,:)*u{k} + (10^3)*(eps_m{k} + eps_M{k});
%         objective = objective + (x{k+1}-xref)'*Q*(x{k+1}-xref) + ucost(k,:)*u{k} + (10^3)*(eps_m{k} + eps_M{k});
        objective = objective + (x{k+1}-xref)'*(x{k+1}-xref) + ucost(k,:)*u{k} + (10^3)*(eps_m{k} + eps_M{k});
        
        Bu_tmp = bm.Bu + Bvu_aux{k} + Bxu_aux{k};
        x{k+1} = bm.A*x{k} + Bu_tmp*u{k}+bm.Bv*vPred(k,:)';
        
        constraints = [constraints, [ymin(k,:)' - eps_m{k} <= bm.C*x{k} + (Dvu_aux{k}+bm.Du)*u{k} + bm.Dv*vPred(k,:)' <= ymax(k,:)' + eps_M{k};...
                                             umin(k,:)'    <=                     u{k}                                <= umax(k,:)';...
                                                                                 eps_m{k}                             >= 0;...
                                                                                 eps_M{k}                             >= 0]];
    end
    constraints = [constraints , x{1} == x0];

    switch solver
        
        case 'ipopt'
            ops = sdpsettings('solver','ipopt','verbose',1);
            solution=solvesdp(constraints,objective,ops);
            if ~strcmp(solution.info,'Successfully solved (IPOPT)')
                display('Problem is: ')
                solution.info
            end

        case 'cplex'
            ops = sdpsettings('solver','quadprog','verbose',0);
            solution=solvesdp(constraints,objective,ops);
            if ~strcmp(solution.info,'Successfully solved (CPLEX-IBM)')
                display('Problem is: ')
                solution.info
                return
            end
            
        case 'linprog'
            ops = sdpsettings('solver','linprog','verbose',1);
            solution=solvesdp(constraints,objective,ops);
            if ~strcmp(solution.info,'Successfully solved (LINPROG)')
                display('Problem is: ')
                solution.info
            end
            
    end

    % Convert input variables
    input=zeros(N,nu);
    for ii=1:N
        input(ii,:)=double(u{ii});
    end


%     light = bm.Du(2,2)*input(1,2) + v_now(4)*input(1,1);
%     if light + 1e-10 < ymin(1,2)    % not enough light during this time
%         disp(['old light input: ' num2str(input(1,2))]);
%         input(1,2) = (ymin(1,2) - v_now(4)*input(1,1))/bm.Du(2,2);
%         disp(['new light input: ' num2str(input(1,2))]);
%     end

input = input(1,:)';
end