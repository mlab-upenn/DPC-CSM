% clear;
clc;

%% load data and run dpc

data_type = 'est';  % est/orig
dpc_type = 'tree';  % tree/forest
orderAR = 6;
ctrlHzn = 6;

tstart = 744-3*24+1;
tstop = 744+0*24; % should be less than 8592 
% tstart = 6;
% tstop = 200; % should be less than 8592

Q = diag(zeros(12,1)); Q(1,1) = 1e2;
R = diag(1e-3*ones(4,1));
epsc = (10^3);
xref = 22;
dpcsolver = 'cplex'; % cplex solves qp / linprog solves lp

% dynamic price of electricity
spotprice = ones(1488,4);
DRbegin = 21;
DRend = 44;
ecost = 1;
spotprice(tstart+DRbegin-1:tstart+DRend-1,:) = ecost; % change price for selected period

saveResults = 0;

switch dpc_type
    case 'tree'
        addpath(genpath('regression_tree/'));
        prepare_models;
        disp('running DPC with regression trees...');
        run_closedloop_dpc;
        rmpath(genpath('regression_tree/'));
        if saveResults;
            save(['../results/dpcrt-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], ...
            'xvec', 'uvec', 'dvec', 'yvec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop');
        end
    case 'forest'
        addpath(genpath('random_forest/'));
        prepare_models;
        disp('running DPC with random forests...');
        run_closedloop_dpc;
        rmpath(genpath('random_forest/'));
        if saveResults
            save(['../results/dpcrf-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], ...
            'xvec', 'uvec', 'dvec', 'yvec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop');
        end
end


%% plot results

plot_results;

