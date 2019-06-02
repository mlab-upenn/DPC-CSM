function varargout = mydpcfast(cmd, varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solves dpc-problem in condensed form
    %
    % DPC cost = sum_0^N (x(k)-xref)'*Q*(x(k)-xref) + sum_0^(N-1) u(k)'*R*u(k)
    %
    % When cmd = 'init', inputs are: (model, N, yref, umin, umax, Q, R, epsc)
    %  model     - the GP model object
    %  N         - the control horizon
    %  yref      - the constant y reference (in current setup, it's a constant scalar)
    %  umin, umax- constant bounds for control inputs (4)
    %  ymin, ymax- bounds for output (vectors of horizon length)
    %  Q, R      - MPC cost parameters
    %  epsc      - weight for slack variable in the cost function
    % and outputs are: (time, gpmpc)
    %  time      - time for constructing the GP-MPC problem (sec)
    %  gpmpc     - the GP-MPC structure
    %
    % When cmd = 'step', inputs are: (gpmpc, kstep, d, elecCost, ymin,
    % ymax, lightmin, bm)
    %  gpmpc     - GP-MPC structure
    %  kstep     - current step (starting from 1, not counting the lags)
    %  d         - disturbances in the horizon
    %  elecCost  - linear cost for control inputs
    %  ymin, ymax- bounds for output
    %  lightmin  - lower bound for lighting
    %  bm        - building model structure to adjust lighting control
    % and outputs are: (u,fval,iter,time,gpmpc)
    %  u      - first input to be applied to the plant
    %  fval   - value of objective function at optimal control sequence
    %  iter   - number of needed iterations
    %  time   - tic-toc-time for calling solver
    %  gpmpc  - the (updated) GPMPC structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch cmd
        case 'init'
            varargout = cell(1, nargout);
            [varargout{:}] = initialize_gpmpc(varargin{:});
        case 'step'
            varargout = cell(1, nargout);
            [varargout{:}] = solve_step(varargin{:});
        otherwise
            error('Unknown command: %s', cmd);
    end
    
end  % Main function

function [t, gpmpc] = initialize_gpmpc(mygpdata, N, yref, umin, umax, Q, R, epsc)
    % model - struct: input_model, hyp, cov, lik, X, Y
    
    t = tic;
    gpmpc = struct;
    
    gpmpc.input_model = mygpdata.input_model;
    
    % Get all input signals used by this GP
    allinputs = what_is_the_input(gpmpc.input_model, 1:size(mygpdata.X, 2));
    allinputs = unique(allinputs(:,1));
    
    % Construct the GP object
    gpmpc.gp = nextgp.GP.createGP(mygpdata.hyp, [], mygpdata.cov, mygpdata.lik, mygpdata.X, mygpdata.Y);
    % sn2 = exp(2*gpmpc.gp.m_hyplik);  % Variance of the default Gaussian likelihood
    
    % Set up signals
    SV = SignalsValues();
    allsignals = {
        'd1'
        'd2'
        'd3'
        'd4'
        'd5'
        'd6'
        'd7'
        'tod'
        'dow'
        'blind'
        'light'
        'heating'
        'cooling'
        'temp'};
    for k = 1:length(allsignals)
        signalname = allsignals{k};
        if isfield(mygpdata.normalization, signalname)
            normdata = mygpdata.normalization.(signalname);
            SV.addSignal(signalname, [], [normdata.min, normdata.max]);
        else
            SV.addSignal(signalname);
        end
    end
    
    SV.addSignal('success');
    SV.addSignal('t_lingp');
    SV.addSignal('t_lingp_sub');
    SV.addSignal('lingp_code');
    SV.addSignal('lingp_iters');
    SV.addSignal('pred_y');
    SV.addSignal('pred_s2');
    
    SM = SignalsModel(SV);
    SM.setInputs(gpmpc.input_model);
    SM.setOutput('temp');
    
    gpmpc.SV = SV;
    gpmpc.SM = SM;
    gpmpc.normalization = mygpdata.normalization;
    gpmpc.maxlag = max(cellfun(@(s) SM.getMaxLag(s), allsignals));
    assert(gpmpc.maxlag >= 0);
    
    % Prepare the linGPMPC
    lingp_parallel = false;  % if true, must start parpool in advance
    lingp_variance = 'full'; % We don't use the variance, so 'full' is probably faster
    lingp_adaptive = true;
    mpc = linGPMPC(SM, N, gpmpc.gp, false, lingp_variance, lingp_adaptive, lingp_parallel);
    
    % Set the signal types
    allcontrols = {'blind', 'light', 'heating', 'cooling'};
    % Only set those control inputs used by the GP
    gpmpc.includedcontrols = ismember(allcontrols, allinputs);
    mpc.setControlSignals(allcontrols(gpmpc.includedcontrols));
    mpc.setOutputSignal('temp');
    
    % User optimization variables
    % Here we use slack variables as in the DPC-tree formulation for a fair
    % comparison, not the GP variance (which is actually a better option).
    gpmpc.sdp_eps = sdpvar(N, 1);
    mpc.setUserVariables({gpmpc.sdp_eps});
    
    % User parameters
    % Electricity cost, to be multiplied with the control inputs
    gpmpc.sdp_elecCost = sdpvar(N,length(gpmpc.includedcontrols));
    % output bounds for the horizon
    gpmpc.sdp_ymin = sdpvar(N, 1);
    gpmpc.sdp_ymax = sdpvar(N, 1);
    mpc.setUserParameters({gpmpc.sdp_elecCost, gpmpc.sdp_ymin, gpmpc.sdp_ymax});
    
    % Construct the MPC problem
    subproblem_solver = 'gurobi+';
    
    Q = Q(1,1);  % only the room temperature
    
    % Select the objective function
    lingp_objfun = @(ymean, yvar, u) sdp_objfun(ymean, yvar, u, gpmpc.sdp_eps, yref, gpmpc.includedcontrols, Q, R, gpmpc.sdp_elecCost, epsc, N);
    
    % Constraints (see linGPMPC manual)
    lingp_cfun = @(ymean, yvar, u) sdp_cfun(ymean, gpmpc.sdp_eps, gpmpc.sdp_ymin, gpmpc.sdp_ymax);
    
    % User constraints (see linGPMPC manual)
    lingp_userfun = @(u) sdp_userfun(u, umin, umax, gpmpc.sdp_eps);
    
    softconstraints_weight = 100;
    
    mpc.constructMPCProblem(...
        lingp_objfun, ...
        lingp_cfun,...
        lingp_userfun,...
        softconstraints_weight, ...
        sdpsettings('solver', subproblem_solver, 'cachesolvers', 1));
    
    gpmpc.mpc = mpc;
    
    % Initialize previous solution for control
    gpmpc.prev_u = struct;
    for k = 1:numel(allcontrols)
        if isfield(mygpdata.normalization, allcontrols{k})
            normdata = mygpdata.normalization.(allcontrols{k});
            gpmpc.prev_u.(allcontrols{k}) = (normdata.max + normdata.min)/2 * ones(N, 1);
        else
            gpmpc.prev_u.(allcontrols{k}) = zeros(N, 1);
        end
    end
    
    t = toc(t);
end  % initialize_gpmpc

function [uopt, fopt, iters, t, gpmpc] = solve_step(gpmpc, kstep, d, elecCost, ymin, ymax, lightmin, bm)
    % Solve one step of the GPMPC
    if gpmpc.maxlag > 0
        kstep = kstep + gpmpc.maxlag;  % the true step, considering lag
    end
    
    lingp_options = struct('trust_size', 0.05, 'beta_fail', 0.5, 'beta_success', 2,...
        'trust_size_max', 1, 'trust_size_min', 0.0001, ...
        'tol', 1e-3, 'projection', true, 'max_iter', 200);
    
    % Adjust previous solution
    allfields = fieldnames(gpmpc.prev_u);
    for k = 1:numel(allfields)
        gpmpc.prev_u.(allfields{k}) = gpmpc.prev_u.(allfields{k})([2:end, end]);
    end
    
    t = tic;
    [U, USER, status, stats, trust_size, rho, J_exact] = ...
        gpmpc.mpc.solveMPC(kstep, gpmpc.prev_u, ...
        {elecCost, ymin(:), ymax(:)},...
        lingp_options);
    t = toc(t);
    
    % disp(trust_size(end));
    
    gpmpc.prev_u = U;
    
    uopt = zeros(4,1);
    if isfield(U, 'blind'), uopt(1) = U.blind(1); end
    if isfield(U, 'light'), uopt(2) = U.light(1); end
    if isfield(U, 'heating'), uopt(3) = U.heating(1); end
    if isfield(U, 'cooling'), uopt(4) = U.cooling(1); end
    
    % instataneous correction for insufficient light
    light = bm.Du(2,2)*uopt(2) + d(4,1)*uopt(1);
    if light + 1e-10 < lightmin    % not enough light during this time
        uopt(2) = (lightmin - d(4,1)*uopt(1))/bm.Du(2,2);
    end
    
    % Record controls
    gpmpc.SV.setSignal('blind', kstep, uopt(1));
    gpmpc.SV.setSignal('light', kstep, uopt(2));
    gpmpc.SV.setSignal('heating', kstep, uopt(3));
    gpmpc.SV.setSignal('cooling', kstep, uopt(4));
    
    % Prediction
    Xk = gpmpc.SM.getIOVectors(kstep, 'removeNaNs', true);
    [~, ~, ypred, s2pred] = gpmpc.gp.predict(Xk);
    gpmpc.SV.setSignal('pred_y', kstep, SignalsValues.postNorm(ypred, gpmpc.normalization.temp.min, gpmpc.normalization.temp.max));
    gpmpc.SV.setSignal('pred_s2', kstep, SignalsValues.postNormVar(s2pred, gpmpc.normalization.temp.min, gpmpc.normalization.temp.max));
    
    gpmpc.SV.setSignal('t_lingp', kstep, t);
    gpmpc.SV.setSignal('lingp_code', kstep, status.code);
    gpmpc.SV.setSignal('success', kstep, status.code == 0);
    gpmpc.SV.setSignal('t_lingp_sub', kstep, stats.sub_time);
    gpmpc.SV.setSignal('lingp_iters', kstep, stats.iters);
    
    
    fopt = nan;
    iters = stats.iters;
end
