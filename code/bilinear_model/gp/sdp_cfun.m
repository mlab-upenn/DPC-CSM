function [eq, ieq] = sdp_cfun(ymean, epsvars, ymin, ymax)
    % Constraints
    ieq = [ymin - ymean - epsvars; ymean - ymax - epsvars];
    eq = [];
end