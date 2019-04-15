function input = FakeBilinearMPC(x0, bm, sim, v_now, vPred, ymin, ymax, ucost, umin, umax)

    yalmip('clear')

    % Number of states and inputs

    nx = size(bm.A,1);
    ny = size(bm.C,1);
    nu = size(bm.Bu,2);
    nv = size(bm.Bv,2);
    N  = sim.predhor;
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
    
    xref = 22*ones(12,1);
    Q = sqrt(diag(1e2*ones(12,1)));
    R = sqrt(diag(1e-3*ones(4,1)));
    for k = 1:N
        
        objective = objective + ucost(k,:)*u{k} + (10^3)*(eps_m{k} + eps_M{k});
%         objective = objective + (norm(Q*(x{k+1}-xref),2))^2 + norm(R*u{k},2) + ucost(k,:)*u{k} + (10^3)*(eps_m{k} + eps_M{k});
        Bu_tmp = bm.Bu + Bvu_aux{k} + Bxu_aux{k};
        x{k+1} = bm.A*x{k} + Bu_tmp*u{k}+bm.Bv*vPred(k,:)';
        
        constraints = [constraints, [ymin(k,:)' - eps_m{k} <= bm.C*x{k} + (Dvu_aux{k}+bm.Du)*u{k} + bm.Dv*vPred(k,:)' <= ymax(k,:)' + eps_M{k};...
                                             umin(k,:)'    <=                     u{k}                                <= umax(k,:)';...
                                                                                 eps_m{k}                             >= 0;...
                                                                                 eps_M{k}                             >= 0]];
    end
    constraints = [constraints , x{1} == x0];


    % profile on
%     ops = sdpsettings('solver','ipopt','verbose',1,'ipopt.max_iter',1000000,'ipopt.max_cpu_time',100);  
%     solution=solvesdp(constraints,objective,ops)
%     if ~strcmp(solution.info,'Successfully solved (IPOPT)')
%         display('Problem is: ')
%         solution.info
%     end

%     ops = sdpsettings('solver','cplex','verbose',0);
%     solution=solvesdp(constraints,objective,ops);
%     if ~strcmp(solution.info,'Successfully solved (CPLEX-IBM)')
%         display('Problem is: ')
%         solution.info
%         return
%     end

    ops = sdpsettings('solver','linprog','verbose',1);  
    solution=solvesdp(constraints,objective,ops);
    if ~strcmp(solution.info,'Successfully solved (LINPROG)')
        display('Problem is: ')
        solution.info
    end
    % Convert input variables


    input=zeros(N,nu);
    for ii=1:N
        input(ii,:)=double(u{ii});
    end

    % Eliminate NaN variable

%     for ii=1:N
%         nan_inputs=isnan(input{ii});
%         elements=find(nan_inputs);
%         for kk=1:length(elements)
%             input{ii}(elements(kk))=0;
%         end
%     end
%     light = bm.Du(2,2)*input(1,2) + v_now(4)*input(1,1);
%     if light + 1e-10 < ymin(1,2)    % not enough light during this time
%         disp(['old light input: ' num2str(input(1,2))]);
%         input(1,2) = (ymin(1,2) - v_now(4)*input(1,1))/bm.Du(2,2);
%         disp(['new light input: ' num2str(input(1,2))]);
%     end

end