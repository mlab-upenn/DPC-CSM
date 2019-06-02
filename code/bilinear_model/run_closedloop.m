

%% load data

bldgmodel = load('data/bm.mat');
bldgmodel = bldgmodel.bm;


%% close-loop control

Ts = 3600; % sampling time
t = 1*Ts:Ts:1488*Ts;
nT = length(t);
nu = 4; % no of inputs
nd = 8; % no of disturbances
nx = 12; % no of states
ny = 3;

% select disturbance
dvec = d;

% pre-allocate
xvec = nan(nx, nT);
yvec = nan(ny, nT);
uvec = nan(nu, nT);
fvalvec = nan(1, nT);
costvec = nan(1, nT);
timings = nan(1, nT-1);
iterations = nan(1, nT-1);

% initialize
x0 = [ 22 22 22 22 22 22 22 22 22 20 15 10 ]';
xvec(:,tstart-orderAR+1:tstart) = repmat(x0,[1,orderAR]);

if strcmpi(sim_type, 'gpfast')
    % Initialize the GPMPC
    [tinit, GPMPC] = mydpcfast('init', GPMPCDATA, ctrlHzn, xref, umin_all(1,:), umax_all(1,:), Q, R, epsc);
    fprintf('Initializing fast GP-MPC took %g secs.\n', tinit);
    
    % Place all non-control and non-output signals
    dpast = dvec(:, tstart-GPMPC.maxlag:tstop+ctrlHzn+1);
    dpast(3,:) = [];  % 3rd element of disturbance is always 0
    idx = 1:size(dpast,2);
    for k = 1:7
        GPMPC.SV.setSignal(sprintf('d%d', k), idx, dpast(k,:));
    end
    % TOD and DOW need to be applied by cosine function
    GPMPC.SV.setSignal('tod', idx, cos(proxy(1, tstart-GPMPC.maxlag:tstop+ctrlHzn+1)/24*2*pi));
    GPMPC.SV.setSignal('dow', idx, cos(proxy(2, tstart-GPMPC.maxlag:tstop+ctrlHzn+1)/7*2*pi));
    
    % Initialize the signals with autoregressive values
    if GPMPC.maxlag > 0
        idx = 1:GPMPC.maxlag;
        GPMPC.SV.setSignal('blind', idx, u(1, tstart-GPMPC.maxlag:tstart-1));
        GPMPC.SV.setSignal('light', idx, u(2, tstart-GPMPC.maxlag:tstart-1));
        GPMPC.SV.setSignal('heating', idx, u(3, tstart-GPMPC.maxlag:tstart-1));
        GPMPC.SV.setSignal('cooling', idx, u(4, tstart-GPMPC.maxlag:tstart-1));
        GPMPC.SV.setSignal('temp', idx, xvec(1, tstart-GPMPC.maxlag:tstart-1));
    end
end

for idk = tstart:tstop  %orderAR:nT-ctrlHzn

    dcur = dvec(:,idk:idk+ctrlHzn-1);
    xcur = fliplr(xvec(1,idk-orderAR+1:idk));
    
    switch sim_type
        case {'tree', 'forest'}
            [uopt, fvalvec(idk), ~, ] = mydpc(models,leafmodels,idk,xref,xcur,dvec,proxy(:,idk)',cost_all,ctrlHzn,umin_all,umax_all,ymin_all,ymax_all,dpcsolver,Q,R,epsc);
    
        case 'mpc'
            [uopt,fvalvec(idk)] = mympc(xvec(:,idk),xref,idk,bldgmodel,ctrlHzn,dcur',ymin_all,ymax_all,cost_all,umin_all,umax_all,mpcsolver,Q,R,epsc);
            
        case 'gpfast'
            [uopt,fvalvec(idk)] = mydpcfast('step', GPMPC, idk-tstart+1, dcur,...
                cost_all(idk:idk+ctrlHzn-1,:),...
                ymin_all(idk+1:idk+ctrlHzn,1), ymax_all(idk+1:idk+ctrlHzn,1),...
                ymin_all(idk+1,2), bldgmodel);  
    end
    
    uvec(:,idk) = uopt;
    [xvec(:,idk+1), yvec(:,idk+1)] = simulate_model(bldgmodel, xvec(:,idk), dcur(:,1), uopt);
    
    % Calculate the true cost at this step
    costvec(idk) = (xvec(1,idk)-xref)'*Q(1,1)*(xvec(1,idk)-xref) + uopt'*R*uopt + cost_all(idk,:)*uopt;
    
    if strcmpi(sim_type, 'gpfast')
        % Record output
        GPMPC.SV.setSignal('temp', idk-tstart+1+GPMPC.maxlag, yvec(1,idk+1));
    end
    
    disp(['iteration ' num2str(idk-tstart+1) ' of ' num2str(tstop-tstart+1)]);
    
end

if strcmpi(sim_type, 'gpfast')
    idx = 1:tstop-tstart+1+GPMPC.maxlag;
    lingp_success = GPMPC.SV.getSignal('success', idx);
    lingp_t_total = GPMPC.SV.getSignal('t_lingp', idx);
    lingp_t_sub = GPMPC.SV.getSignal('t_lingp_sub', idx);
    lingp_code = GPMPC.SV.getSignal('lingp_code', idx);
    lingp_iters = GPMPC.SV.getSignal('lingp_iters', idx);
    pred_y = GPMPC.SV.getSignal('pred_y', idx);
    pred_s2 = GPMPC.SV.getSignal('pred_s2', idx);
end
