function [cost_all, ymin_all, ymax_all, umin_all, umax_all] = get_constraints()

sumIdx = 2881:1:4416+2880;
winIdx = setdiff(1:8592,sumIdx);
cost_all = zeros(8592,4);
ymin_all = zeros(8592,3);
ymax_all = zeros(8592,3);

% check whether occupied or not
load('../mpc/InData/occupancy/ig_act_constr_1we.mat');    % load empirical diurnal profile of internal gain (people)

ig_Pers = repmat(ig_act_constr_1we,[51,1]); % all 0,1
ig_Pers = [ig_Pers; ig_Pers(1:24)];
ig_Pers = 0.3 * ig_Pers;    % expected 0.4
ig_Equip = ig_Pers + 0.1;   % same as in old BacLab
ig = [ig_Pers ig_Equip];

% stack output constraints together
for idx = sumIdx
    [cost_all(idx,:), ymin_all(idx,:), ymax_all(idx,:)] = getOutputConstraints('sum', ig(idx,:));
end
for idx = winIdx
    [cost_all(idx,:), ymin_all(idx,:), ymax_all(idx,:)] = getOutputConstraints('win', ig(idx,:));
end

% stack input constraints together
umin = [0      0          0           0       ];
umax = [1      7.2        1           2*11.488  ];
umin_all = repmat(umin, [8592,1]);
umax_all = repmat(umax, [8592,1]);

save('../data/costsAndConstraints.mat', 'cost_all','ymin_all','ymax_all','umin_all','umax_all');

end

function [ucost, ymin, ymax] = getOutputConstraints(season, ig)

if strcmp(season, 'sum')
    ucost = [ 0.01  3.32  7.47  1.107 ];  % 0.01
    % y = [ 1Troom  2Eroom  3TsrfCeil ];
    ymin = [ 22  500  18 ];
    ymax = [ 26  1.0E+6  1000000 ]; %26
    if ig(1) ==0
        ymin = [ 5  0  5 ];
        ymax = [40  1.0E+6  1000000 ];
    end
    
else
    ucost = [ 0  3.32  7.47  1.107 ];
    % y = [ 1Troom  2Eroom  3TsrfCeil ];
    ymin = [ 21  500  18 ];
    ymax = [ 25  1.0E+6  1000000 ];   %25
    if ig(1) ==0
        ymin = [ 5  0  5 ];
        ymax = [40  1.0E+6  1000000 ];
    end
end
end
