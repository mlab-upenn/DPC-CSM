data_type = 'est';      % est/orig
dpc_type = 'forest';    % tree/forest
orderAR = 6;
ctrlHzn = 6;

tstart = 744-3*24+1;
tstop = 744+0*24;       % should be less than 8592 
steps = tstart:tstop;
nsteps = length(steps);
DRbegin = 21;
DRend = 44;

figure(1); hold on;
figure(2); hold on;
figure(3); hold on;

for cost = [1,5,10,20,50]
    
    load(['../results/dpcrf-' data_type '-cost' num2str(cost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], ...
        'xvec', 'uvec', 'dvec', 'yvec', 'umin_all','umax_all','ymin_all','ymax_all','tstart','tstop');
    
    figure(1);
    plot(yvec(1,steps), 'LineWidth', 2);
    
    figure(2);
    plot(uvec(3,steps), 'LineWidth', 2);
    
    figure(3);
    plot(uvec(4,steps), 'LineWidth', 2);
end

figure(1);
ymax_all(ymax_all(:,1)>35,:) = 27;
ymin_all(ymin_all(:,1)<10,:) = 19;
plot(ymax_all(steps,1), 'k')
plot(ymin_all(steps,1), 'k')
plot(DRbegin*ones(100,1),linspace(min(ymin_all(:,1)),max(ymax_all(:,1)),100), 'r', 'LineWidth', 2)
plot(DRend*ones(100,1),linspace(min(ymin_all(:,1)),max(ymax_all(:,1)),100), 'r', 'LineWidth', 2)
legend(' 1x', ' 5x', '10x', '20x', '50x')

figure(2);
plot(umax_all(steps,3), 'k')
plot(umin_all(steps,3), 'k')
plot(DRbegin*ones(100,1),linspace(0,1,100), 'r', 'LineWidth', 2)
plot(DRend*ones(100,1),linspace(0,1,100), 'r', 'LineWidth', 2)
legend(' 1x', ' 5x', '10x', '20x', '50x')

figure(3);
plot(umax_all(steps,4), 'k')
plot(umin_all(steps,4), 'k')
plot(DRbegin*ones(100,1),linspace(0,umax_all(1,4),100), 'r', 'LineWidth', 2)
plot(DRend*ones(100,1),linspace(0,umax_all(1,4),100), 'r', 'LineWidth', 2)
legend(' 1x', ' 5x', '10x', '20x', '50x')