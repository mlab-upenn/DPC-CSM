function u = get_control_input(x0, bm, sim, v_now, vPred, ymin, ymax, ucost, umin, umax)
    % Set parameters
    % the first light input parameter is adjusted automatically to have zero violations
   
    nx = size(bm.A,1);
    ny = size(bm.C,1);
    nu = size(bm.Bu,2);
    nv = size(bm.Bv,2);
    N_mpc  = sim.predhor;
    nslack = 2*ny*N_mpc;      % size of slackMax, slackMin
    
    % reformulate QP in the following, simple shape
    % X = [U slackMax slackMin]
    %       min {ucost_lin' * X} + X'*H*X
    % s.t.  [I ; -I]*U <= [umax ; umin], i.e. A_ineq1*U <= b_ineq1
    %       ymin <= y <= ymax, i.e. A_ineq2*U <= b_ineq2
    %       where, Y = tmp1*x0 + tmp2*V + A_ineq2*U
    
    ucost_lin = [reshape(ucost',N_mpc*nu,1) ; 10^3*ones(nslack,1)];   % linear part
   
    % build inequality 1 out of 4
    A_ineq1 = [eye(N_mpc*nu) ; -eye(N_mpc*nu)];
    A_ineq1 = [A_ineq1 zeros(2*nu*N_mpc,nslack)];
    b_ineq1 = [reshape(umax',N_mpc*nu,1) ; -reshape(umin',N_mpc*nu,1)];
    
    % The next cycle creates the following structure for tmp1
    
    % tmp1=y(1:N)=[C;C;C;CA;CA;CA;...;CA^N-1;CA^N-1;CA^N-1]
    
    tmp1 = nan(N_mpc*ny,nx);
    for ii = 1 : N_mpc
        tmp1(1+(ii-1)*ny:ii*ny,:) = bm.C*bm.A^(ii-1);
    end
    
    % The next cycle creates the following structure for tmp2
    
    % tmp2 = [ D_v       0         0         ...   0
    %          CBv       Dv        0         ...   0
    %          CABv      CBv       Dv        ...   0
    %           .         .         .    .
    %           .         .         .        .
    %           .         .         .            .
    %          CA^58Bv   CA^57Bv   CA^56Bv   ...   Dv ]
    
    tmp2 = zeros(N_mpc*ny,nv*N_mpc);
    for ii = 0 : N_mpc-1
        for jj = 0 : ii-1
            tmp2(1+ii*ny:(ii+1)*ny, 1+jj*nv:(jj+1)*nv) = bm.C*bm.A^(ii-1-jj)*bm.Bv;
        end
        tmp2(1+ii*ny: (ii+1)*ny, 1+ii*nv:(ii+1)*nv) = bm.Dv;
    end

    % precompute matrices   B_{u,0},...,B_{u,N-1}
    %                       D_{u,0},...,D_{u,N-1}
    
    Bvu_aux = cell(N_mpc,1);
    Bxu_aux = cell(N_mpc,1);
    Dvu_aux = cell(N_mpc,1);
    
    for jj = 1 : N_mpc  % for inputs u_0,...,u_(N_mpc-1)
        Bvu_aux{jj} = nan(nx,nu);
        Bxu_aux{jj} = nan(nx,nu);
        Dvu_aux{jj} = nan(ny,nu);
            
        for ii = 1 : nu
            Bvu_aux{jj}(:,ii) = bm.Bvu(:,:,ii) * vPred(jj,:)'; % linearized
            Bxu_aux{jj}(:,ii) = bm.Bxu(:,:,ii) * x0;   % linearized
            Dvu_aux{jj}(:,ii) = bm.Dvu(:,:,ii) * vPred(jj,:)';  % linearized
        end 
    end
    
    % The next cycle creates the following structure for A_ineq2
    
    % tmp2 = [ Dvu_aux{1} + Du    0                  0                  ...     0
    %          CBu_tmp            Dvu_aux{2} + Du    0                  ...     0
    %          CABu_tmp           CBu_tmp            Dvu_aux{3} + Du    ...     0
    %             .                  .                  .             .
    %             .                  .                  .                 .
    %             .                  .                  .                     .
    %          CA^58Bu_tmp        CA^57Bu_tmp        CA^56Bu_tmp        ...     Dvu_aux{59} + Du ]
    
    A_ineq2_tmp = zeros(N_mpc*ny,nu*N_mpc);
    for ii = 0 : N_mpc-1
        for jj = 0 : ii-1
            Bu_tmp = bm.Bu + Bvu_aux{jj+1} + Bxu_aux{jj+1};
            A_ineq2_tmp(1+ii*ny:(ii+1)*ny, 1+jj*nu:(jj+1)*nu) = bm.C*bm.A^(ii-1-jj)*Bu_tmp;
        end
        A_ineq2_tmp(1+ii*ny:(ii+1)*ny, 1+ii*nu:(ii+1)*nu) = Dvu_aux{ii+1} + bm.Du;
    end
    
    V = reshape(vPred',N_mpc*nv,1);
    
    % build inequality 2 of 4: soft constraint of upper bound
    % It basically implement the following:
    
    % line 1: (Dvu_aux{1} + Du)u{0} - slackMax <= [yM(1);yM(2);yM(3)] - [C;C;C]x0 - Dv vPred(0) ==> y(0) <= yM + slackMax
    % line 2: CBu_tmp u{0} + (Dvu_aux{1} + Du)u{1} - slackMax <= 
    %                   [yM(1);yM(2);yM(3)] - [CA;CA;CA]x0 - CBv vPred(0) - Dv vPred(1) ==> y(1) <= yM + slackMax
    % and so on
    
    A_ineq2 = [A_ineq2_tmp -eye(size(A_ineq2_tmp,1)) zeros(size(A_ineq2_tmp,1),nslack-size(A_ineq2_tmp,1))];
    b_ineq2 = reshape(ymax',N_mpc*ny,1) - tmp1*x0 - tmp2*V;
    
    %build inequality 3 of 4: soft constraint of lower bound
    % This is basically the same as inequality 2 but for the lower bound
    
    A_ineq3 = [-A_ineq2_tmp zeros(size(A_ineq2_tmp,1),size(A_ineq2_tmp,1)) -eye(size(A_ineq2_tmp,1)) zeros(size(A_ineq2_tmp,1),nslack-2*size(A_ineq2_tmp,1))];
    b_ineq3 = -reshape(ymin',N_mpc*ny,1) + tmp1*x0 + tmp2*V;
    
    % build inequality 4 of 4: slack variable must be positive
    A_ineq4 = [zeros(nslack,N_mpc*nu) -eye(nslack)];
    b_ineq4 = zeros(nslack,1);
    
    A_ineq = [A_ineq1; A_ineq2 ; A_ineq3 ; A_ineq4 ];
    b_ineq = [b_ineq1; b_ineq2 ; b_ineq3 ; b_ineq4 ];
            
    OPTIONS.verbose = 1; % verbose output
    OPTIONS.logfile = 0; % log into file
    OPTIONS.probtype = 0; % Quadratic Problem is 1
    OPTIONS.lic_rel=9000;
    PARAM.int=[1063,2]; % Tipp from Martin
        
    [ u_opt, fmin, status, details ] = cplexlp(ucost_lin, ...
        A_ineq, b_ineq);
    u = reshape(u_opt(1:nu*N_mpc), nu, N_mpc)';
    
    % simulate instantaneous light correction in first step
%     light = bm.Du(2,2)*u(1,2) + v_now(4)*u(1,1);
%     if light + 1e-10 < ymin(1,2)    % not enough light during this time
%         disp(['old light input: ' num2str(u(1,2))]);
%         u(1,2) = (ymin(1,2) - v_now(4)*u(1,1))/bm.Du(2,2);
%         disp(['new light input: ' num2str(u(1,2))]);
%     end
end