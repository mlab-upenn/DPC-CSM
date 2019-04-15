

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
timings = nan(1, nT-1);
iterations = nan(1, nT-1);

% initialize
x0 = [ 22 22 22 22 22 22 22 22 22 20 15 10 ]';
xvec(:,tstart-orderAR+1:tstart) = repmat(x0,[1,orderAR]);
    
for idk = tstart:tstop  %orderAR:nT-ctrlHzn

    dcur = dvec(:,idk:idk+ctrlHzn-1);
    xcur = fliplr(xvec(1,idk-orderAR+1:idk));
    
    switch sim_type
        case {'tree', 'forest'}
            [uvec(:,idk),fvalvec(idk)] = mydpc(models,leafmodels,idk,xref,xcur,dvec,proxy(:,idk)',cost_all,ctrlHzn,umin_all,umax_all,ymin_all,ymax_all,dpcsolver,Q,R,epsc);
    
        case 'mpc'
            [uvec(:,idk),fvalvec(idk)] = mympc(xvec(:,idk),xref,idk,bldgmodel,ctrlHzn,dcur',ymin_all,ymax_all,cost_all,umin_all,umax_all,mpcsolver,Q,R,epsc);
    end
    
    [xvec(:,idk+1), yvec(:,idk+1)] = simulate_model(bldgmodel, xvec(:,idk), dcur(:,1), uvec(:,idk));
    
    disp(['iteration ' num2str(idk-tstart+1) ' of ' num2str(tstop-tstart+1)]);
    
end
