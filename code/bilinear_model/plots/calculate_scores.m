data_type = 'est';      % est/orig
dpc_type = 'tree';    % tree/forest
orderAR = 6;
ctrlHzn = 6;

tstart = 744-3*24+1;
tstop = 744+3*24;  
idx = tstart:tstop;
nsteps = length(idx);
ecost = 1;

load(['../data/test-' data_type '-filtered.mat']);
inputcost = cost_all(idx,:);

load(['../results/mpc-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat']);
energy_mpc = sum(sum(inputcost'.*uvec(:,idx)))/length(idx);
meandev_mpc = sum(abs(xvec(1,idx)-22))/length(idx);
objval_mpc = sum(fvalvec(idx))/length(idx);
u_mpc = uvec;

load(['../results/dpcrf-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat']);
energy_dpcrf = sum(sum(inputcost'.*uvec(:,idx)))/length(idx);
meandev_dpcrf = sum(abs(xvec(1,idx)-22))/length(idx);
objval_dpcrf = sum(fvalvec(idx))/length(idx);
u_dpcrf = uvec;

load(['../results/dpcrt-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat']);
energy_dpcrt = sum(sum(inputcost'.*uvec(:,idx)))/length(idx);
meandev_dpcrt = sum(abs(xvec(1,idx)-22))/length(idx);
objval_dpcrt = sum(fvalvec(idx))/length(idx);
u_dpcrt = uvec;

SS_tot = sum((u_mpc(3,idx)-mean(u_mpc(3,idx))).^2);
SS_reg_dpcrt = sum((u_dpcrt(3,idx)-u_mpc(3,idx)).^2);
SS_reg_dpcrf = sum((u_dpcrf(3,idx)-u_mpc(3,idx)).^2);
R2_dpcrt = 1-SS_reg_dpcrt/SS_tot;
R2_dpcrf = 1-SS_reg_dpcrf/SS_tot;

SS_tot_heat = sum((u_mpc(4,idx)-mean(u_mpc(4,idx))).^2);
SS_reg_dpcrt_heat = sum((u_dpcrt(4,idx)-u_mpc(4,idx)).^2);
SS_reg_dpcrf_heat = sum((u_dpcrf(4,idx)-u_mpc(4,idx)).^2);
R2_dpcrt_heat = 1-SS_reg_dpcrt_heat/SS_tot_heat;
R2_dpcrf_heat = 1-SS_reg_dpcrf_heat/SS_tot_heat;
