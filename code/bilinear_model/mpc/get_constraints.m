function [ymin_all, ymax_all] = get_constraints(sim, v, ymin, ymax, ymin_sb, ymax_sb, RGS_begin)
% cannot use v b/c we don't know v in advance
    cd('InData/occupancy')
    load('ig_act_constr_1we.mat');    % load empirical diurnal profile of internal gain (people)
    cd ..; cd ..;
    ig_Pers = repmat(ig_act_constr_1we,[52,1]); % all 0,1
    ig_Pers = 0.3 * ig_Pers;    % expected 0.4
    ig_Equip = ig_Pers + 0.1;   % same as in old BacLab
    ig = [ig_Pers ig_Equip];
    if (RGS_begin == 2881)  % start from Tuesday
        ig = ig(1+24:end,:);
    elseif RGS_begin == 7297
        ig = ig(1+3*24:end,:);
    end
    ig = ig(1:sim.simhor+sim.predhor,:);

    for i = 1:sim.simhor+sim.predhor
        if ig(i,1) == 0  % v(5) = internal gains (persons)
            ymin_var(i,:) = ymin_sb;
            ymax_var(i,:) = ymax_sb;
        else
            ymin_var(i,:) = ymin;   % ymin is if people are there
            ymax_var(i,:) = ymax;
        end
    end
    ymin_all = ymin_var(1:sim.simhor+sim.predhor,:);
    ymax_all = ymax_var(1:sim.simhor+sim.predhor,:);
end