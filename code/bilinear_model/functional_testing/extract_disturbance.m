cd('../model/')

% building parameters
list = {'MSM2007',1, 'pa', 'h', 'wh', 'ih' , 60 , 0};
build.weather = list{1};        % 'MSM2007', 'WHW2007'
build.KF = 'KF_AR4_';           % '', 'KF_', 'AR2_'
build.setback = list{2};        % 1: use ymin_sb; 0: ymin_sb = ymin
build.sys = 'e01';              % {e01}
build.type = list{3};           % {Swiss Average, PAssive house}
build.env = list{4};            % {Heavy, Light}
build.win = list{5};            % {Window Low, Window High}
build.int = list{6};            % {Internal gain High, Internal gain Low}
build.dir = 'S';                % {North, South, South East, South West}
build.vac = 1;                  % code for vacancy
build.occ = 5;

% first period
build.seas = 'win';
build_stand = build;
build_stand.vac = NaN;      % use simple occupancy control
sim.wVarNeeded =  [ 1 1 0 0 1 0 ];  % must be consistent with "build.dir"
sim.igVarNeeded = [ 1 1 ];  % Needed internal gains [persG equipG]
sim.simhor = 2880;      % simulate Jan 01 - April 30
sim.predhor = 60;       % 7 days for PB = 168
[v, vPred_all] = get_disturbance(build_stand, sim, 1);  % get disturbance for Jan-April
dist = v(1:2880,:);
distEst1 = zeros(2880,8);
for ii = 1:length(distEst1)
    distEst1(ii,:) = vPred_all{ii}(1,:);
end

% second period
sim.simhor = 7303-2888+1;
build.seas = 'sum';
build_stand = build;
build_stand.vac = NaN;      % use simple occupancy control
[v, vPred_all] = get_disturbance(build_stand, sim, 2888-7);
dist = [dist;v(1:4416,:)];
distEst2 = zeros(4416,8);
for ii = 1:length(distEst2)
    distEst2(ii,:) = vPred_all{ii}(1,:);
end

% third period
build.seas = 'win';                 % winter season
build_stand = build;
build_stand.vac = NaN;      % use simple occupancy control
sim.simhor = 8767-7304-168+1;      % simulate May 01 - Okt 31
[v, vPred_all] = get_disturbance(build_stand, sim, 7304-7);  % get disturbance for May-Oct
dist = [dist ; v(1:1296,:)];
distEst3 = zeros(1296,8);
for ii = 1:length(distEst3)
    distEst3(ii,:) = vPred_all{ii}(1,:);
end
distEst = [distEst1;distEst2;distEst3];

cd('../functional_testing/')

save('../data/dist.mat','dist','distEst')