% clear;
clc;

addpath('plots/');

%% load data and run dpc

data_type = 'est';          % est/orig
sim_type = 'gpfast';           % mpc/tree/forest/gpfast

% simulation settings
% tstart = 744+19*24+1;
% tstop = 744+26*24;           % should be less than 8592
tstart = 744-3*24+1;
tstop = 744+3*24;           % should be less than 8592
% tstart = 10;
% tstop = 200;              % should be less than 8592

% mpc settings
orderAR = 6;
ctrlHzn = 6;
Q = diag(zeros(12,1)); Q(1,1) = 1e2;
R = diag(1e-3*ones(4,1));
epsc = (10^3);
xref = 22;
dpcsolver = 'quadprog'; % quadprog solves qp / linprog solves lp
mpcsolver = 'gurobi'; % cplex / ipopt / gurobi both solve same qp

% dynamic price of electricity
spotprice = ones(1488,4);
DRbegin = 64;
DRend = 68;
ecost = 1;
spotprice(tstart+DRbegin-1:tstart+DRend-1,:) = ecost; % change price for selected period

% load test data
load(['data/test-' data_type '-filtered.mat']);
cost_all = spotprice.*cost_all;
ymin_all(:,1) = ymin_all(:,1)-2;

% save results yes/no
saveResults = 1;

switch sim_type
    
    case 'mpc'
        addpath('mpc/');
        disp('running MPC...');
        run_closedloop;
        rmpath('mpc/');
        if saveResults
            save(['results_truong/mpc-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], ...
                'xvec', 'uvec', 'dvec', 'yvec', 'fvalvec', 'costvec', 'costelec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop',...
                'orderAR', 'ctrlHzn', 'Q', 'R', 'epsc', 'xref', 'DRbegin', 'DRend', 'ecost', 'spotprice');
        end
        
    case 'tree'
        addpath(genpath('dpc/regression_tree/'));
        addpath('dpc/');
        prepare_models;
        disp('running DPC with regression trees...');
        run_closedloop;
        rmpath(genpath('dpc/regression_tree/'));
        rmpath('dpc/');
        if saveResults
            save(['results_truong/dpcrt-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], ...
                'xvec', 'uvec', 'dvec', 'yvec', 'fvalvec', 'costvec', 'costelec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop',...
                'orderAR', 'ctrlHzn', 'Q', 'R', 'epsc', 'xref', 'DRbegin', 'DRend', 'ecost', 'spotprice');
        end
        
    case 'forest'
        addpath(genpath('dpc/random_forest/'));
        addpath('dpc/');
        prepare_models;
        disp('running DPC with random forests...');
        run_closedloop;
        rmpath(genpath('dpc/random_forest/'));
        rmpath(genpath('dpc/'));
        if saveResults
            save(['results_truong/dpcrf-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], ...
                'xvec', 'uvec', 'dvec', 'yvec', 'fvalvec', 'costvec', 'costelec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop',...
                'orderAR', 'ctrlHzn', 'Q', 'R', 'epsc', 'xref', 'DRbegin', 'DRend', 'ecost', 'spotprice');
        end
        
    case 'gpfast'
        addpath(genpath('gp/'));
        selected_gpmodel = 'gp-est-5000-se-meanzero-lag4';
        prepare_models;
        disp('running Fast GP-MPC ...');
        run_closedloop;
        rmpath(genpath('gp/'));
        if saveResults
            save(['results_truong/gpfast-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '-lag4-5000.mat'], ...
                'xvec', 'uvec', 'dvec', 'yvec', 'fvalvec', 'costvec', 'costelec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop',...
                'orderAR', 'ctrlHzn', 'Q', 'R', 'epsc', 'xref', 'DRbegin', 'DRend', 'ecost', 'spotprice', ...
                'selected_gpmodel', 'lingp_success', 'lingp_t_total', 'lingp_t_sub', 'lingp_code', 'lingp_iters', 'pred_y', 'pred_s2');
        end
end


%% plot results

plot_results;

