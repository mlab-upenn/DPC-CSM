function [ig] = get_internalgains(build, sim, RGS_begin, pred_real)
    if strcmp(pred_real,'real')     % realization of occupancy
        cd('InData/occupancy')
        load('ig_act_real.mat')
        cd ..; cd ..;
        
        ig_Pers = ig_act_real(RGS_begin : RGS_begin+sim.simhor+sim.predhor-1);
        ig_Equip = ig_Pers + 0.1;
        ig = [ig_Pers ig_Equip];
        
    elseif strcmp(pred_real,'pred')
        cd('InData/occupancy')
        load('ig_act_pred_1we.mat');
        ig_Pers = repmat(ig_act_pred_1we,[52,1]); % all 0,1
        cd ..; cd ..;
        
        ig_Equip = ig_Pers + 0.1;   % same as in old BacLab
        ig = [ig_Pers ig_Equip];
        
        ig = ig(RGS_begin : RGS_begin+sim.simhor+sim.predhor-1,:); 
    else
        error('err in internal gain')
    end
    
end