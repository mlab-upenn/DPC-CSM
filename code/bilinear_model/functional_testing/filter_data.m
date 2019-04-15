% removes data that will be used for testing

dist_type = 'orig';
matlab_data = 0;
python_data = 1;

% months different in train and test
janIdx = 1:31*24;
mayIdx = 2881:2880+31*24;
allIdx = [janIdx, mayIdx];

% training data: remove data from Jan and May
load(['../data/data-' dist_type '.mat']);
load('../data/costsAndConstraints.mat');
d(:,allIdx) = [];
u(:,allIdx) = [];
x(:,allIdx) = [];
y(:,allIdx) = [];
proxy(:,allIdx) = [];
cost_all(allIdx,:) = [];
ymin_all(allIdx,:) = [];
ymax_all(allIdx,:) = [];
umin_all(allIdx,:) = [];
umax_all(allIdx,:) = [];

if matlab_data
save(['../data/train-' dist_type '-filtered.mat'], 'x', 'u', 'd', 'y', 'proxy',...
    'cost_all','umin_all','umax_all','ymin_all','ymax_all');
end
if python_data
    csvwrite(['../../python/data/train_' dist_type '_x.csv'],x');
    csvwrite(['../../python/data/train_' dist_type '_u.csv'],u');
    csvwrite(['../../python/data/train_' dist_type '_d.csv'],d');
    csvwrite(['../../python/data/train_' dist_type '_y.csv'],y');
end

% test data: keep only Jan and May
load(['../data/data-' dist_type '.mat']);
load('../data/costsAndConstraints.mat');
d = d(:,allIdx);
u = u(:,allIdx);
x = x(:,allIdx);
y = y(:,allIdx);
proxy = proxy(:,allIdx);
cost_all = cost_all(allIdx,:);
ymin_all = ymin_all(allIdx,:);
ymax_all = ymax_all(allIdx,:);
umin_all = umin_all(allIdx,:);
umax_all = umax_all(allIdx,:);

if matlab_data
save(['../data/test-' dist_type '-filtered.mat'], 'x', 'u', 'd', 'y', 'proxy',...
    'cost_all','umin_all','umax_all','ymin_all','ymax_all');
end
if python_data
    csvwrite(['../../python/data/test_' dist_type '_x.csv'],x');
    csvwrite(['../../python/data/test_' dist_type '_u.csv'],u');
    csvwrite(['../../python/data/test_' dist_type '_d.csv'],d');
    csvwrite(['../../python/data/test_' dist_type '_y.csv'],y');
end
