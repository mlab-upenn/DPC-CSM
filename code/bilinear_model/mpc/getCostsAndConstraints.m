function [ucost, ymin,ymin_sb, ymax,ymax_sb, umin, umax] = getCostsAndConstraints(build)          
    
    % u = [1bPos 2eLighting 3fcUsgFact 4hPowRad]
    if strcmp(build.seas, 'sum')
        ucost = [ 0.01  3.32  7.47  1.107 ];  % 0.01
        % y = [ 1Troom  2Eroom  3TsrfCeil ];
        ymin = [ 22  500  18 ];
        ymin_sb = [ 5  0  5 ];
        ymax = [ 26  1.0E+6  1000000 ]; %26
        ymax_sb = [40  1.0E+6  1000000 ];
    else
        ucost = [ 0  3.32  7.47  1.107 ];
        % y = [ 1Troom  2Eroom  3TsrfCeil ];
        ymin = [ 21  500  18 ];
        ymin_sb = [ 5  0  5 ];
        ymax = [ 25  1.0E+6  1000000 ];   %25
        ymax_sb = [40  1.0E+6  1000000 ];
    end
        
    if (~build.setback)
        ymin_sb = ymin;
        ymax_sb = ymax;
    end

    umin = [0      0          0           0       ];
    umax = [1      7.2        1           2*11.488  ]; 

end