steps = 1:size(state_coll,1);
nsteps = length(steps);

%% room temperature
% load('../results/mpc2007.mat')
figure; hold on
title('output')
plot(output_coll(steps,1))
plot(ymax_coll(steps,1))
plot(ymin_coll(steps,1))
xlabel('time step')
ylabel('room temp.')

%% inputs
figure;
title('inputs')
idu = 1;
subplot(4,1,idu); hold on;
plot(input_coll(idu,steps))
plot(umax(idu)*ones(nsteps,1))
plot(umin(idu)*ones(nsteps,1))
ylabel('bPos')

idu = 2;
subplot(4,1,idu); hold on;
plot(input_coll(idu,steps))
plot(umax(idu)*ones(nsteps,1))
plot(umin(idu)*ones(nsteps,1))
ylabel('eLighting')

idu = 3;
subplot(4,1,idu); hold on;
plot(input_coll(idu,steps))
plot(umax(idu)*ones(nsteps ,1))
plot(umin(idu)*ones(nsteps ,1))
ylabel('fcUsgFact')

idu = 4;
subplot(4,1,idu); hold on; 
plot(input_coll(idu,steps))
plot(umax(idu)*ones(nsteps ,1))
plot(umin(idu)*ones(nsteps ,1))
ylabel('hPowRad')

%% disturbances
load('../data/dist.mat')
disturbance = distEst;
figure;
title('disturbances(1:4)')
idv = 1;
subplot(4,1,idv); hold on;
plot(disturbance(steps,1))
ylabel('solG')

idv = 2;
subplot(4,1,idv); hold on;
plot(disturbance(steps,2))
ylabel('dsolG')

idv = 3;
subplot(4,1,idv); hold on;
plot(disturbance(steps,3))
ylabel('illum')

idv = 4;
subplot(4,1,idv); hold on; 
plot(disturbance(steps,4))
ylabel('dIllum')

figure;
title('disturbances(5:8)')
idv = 1;
subplot(4,1,idv); hold on;
plot(disturbance(steps,5))
ylabel('persG')

idv = 2;
subplot(4,1,idv); hold on;
plot(disturbance(steps,6))
ylabel('equipG')

idv = 3;
subplot(4,1,idv); hold on;
plot(disturbance(steps,7))
ylabel('Tair')

idv = 4;
subplot(4,1,idv); hold on; 
plot(disturbance(steps,8))
ylabel('TfreeCool')