steps = tstart:tstop;
nsteps = length(steps);

%% room temperature

figure; hold on; grid on;
title('output')
plot(ymax_all(steps,1))
plot(ymin_all(steps,1))
plot(yvec(1,steps))
plot(DRbegin*ones(100,1),linspace(5,40,100), 'k')
plot(DRend*ones(100,1),linspace(5,40,100), 'k')
xlabel('time step')
ylabel('room temp.')

%% inputs
figure;
title('inputs')
idu = 1;
subplot(4,1,idu); hold on; grid on;
plot(umax_all(steps,idu))
plot(umin_all(steps,idu))
plot(DRbegin*ones(100,1),linspace(0,1,100), 'k')
plot(DRend*ones(100,1),linspace(0,1,100), 'k')
plot(uvec(idu,steps))
ylabel('bPos')

idu = 2;
subplot(4,1,idu); hold on; grid on;
plot(umax_all(steps,idu))
plot(umin_all(steps,idu))
plot(uvec(idu,steps))
ylabel('eLighting')

idu = 3;
subplot(4,1,idu); hold on; grid on;
plot(umax_all(steps,idu))
plot(umin_all(steps,idu))
plot(DRbegin*ones(100,1),linspace(0,1,100), 'k')
plot(DRend*ones(100,1),linspace(0,1,100), 'k')
plot(uvec(idu,steps))
ylabel('fcUsgFact')

idu = 4;
subplot(4,1,idu); hold on; grid on; 
plot(umax_all(steps,idu))
plot(umin_all(steps,idu))
plot(uvec(idu,steps))
ylabel('hPowRad')

%% disturbances

disturbance = dvec';
figure;
title('disturbances(1:4)')
idv = 1;
subplot(4,1,idv); hold on; grid on;
plot(disturbance(steps,1))
ylabel('solG')

idv = 2;
subplot(4,1,idv); hold on; grid on;
plot(disturbance(steps,2))
ylabel('dsolG')

idv = 3;
subplot(4,1,idv); hold on; grid on;
plot(disturbance(steps,3))
ylabel('illum')

idv = 4;
subplot(4,1,idv); hold on; grid on; 
plot(disturbance(steps,4))
ylabel('dIllum')

figure;
title('disturbances(5:8)')
idv = 1;
subplot(4,1,idv); hold on; grid on;
plot(disturbance(steps,5))
ylabel('persG')

idv = 2;
subplot(4,1,idv); hold on; grid on;
plot(disturbance(steps,6))
ylabel('equipG')

idv = 3;
subplot(4,1,idv); hold on; grid on;
plot(disturbance(steps,7))
ylabel('Tair')

idv = 4;
subplot(4,1,idv); hold on; grid on; 
plot(disturbance(steps,8))
ylabel('TfreeCool')